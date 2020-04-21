import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('supermongo')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

file1 = "test_opacity_file.txt"
file2 = "sph.opacities_ferguson_yes_grains_and_molecules"

rho1 = np.loadtxt(file1,max_rows=1)
data1 = np.loadtxt(file1,skiprows=1)
T1 = data1[:,0]

grid = np.loadtxt(file2,skiprows=2,max_rows=1)
grid = grid.astype(int)
d = np.loadtxt(file2,skiprows=3,max_rows=grid[0]*grid[1])

rho2 = np.unique(d[:,0])
T2 = np.unique(d[:,1])

data2 = np.zeros((grid[1],grid[0]))
for i in range(0,grid[0]):
    for j in range(0,grid[1]):
        drow = d[i*grid[1] + j]
        data2[j][i] = drow[2]
        
data2 = np.log10(data2)

vmin = min(np.amin(data1[np.isfinite(data1)]), np.amin(data2[np.isfinite(data2)]))
vmax = max(np.amax(data1[np.isfinite(data1)]), np.amax(data2[np.isfinite(data2)]))

#vmin = -5
#vmax = 10

Tmin = min(np.amin(T1),np.amin(T2))
Tmax = max(np.amax(T1),np.amax(T2))
rhomin = min(np.amin(rho1),np.amin(rho2))
rhomax = max(np.amax(rho1),np.amax(rho2))

extent1 = (np.amin(rho1),np.amax(rho1),np.amin(T1),np.amax(T1))
extent2 = (np.amin(rho2),np.amax(rho2),np.amin(T2),np.amax(T2))

my_cm = plt.cm.get_cmap('nipy_spectral')

my_cm.set_bad('magenta')

#wspace = 0.
#left = 0.05

nrows = 1
ncols = 2

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,dpi=100,figsize=(16,8))
plt.subplots_adjust(left=0.05,right=0.9)
#plt.subplots_adjust(wspace=wspace,left=left,right=0.8,top=0.95,bottom=0.15)
im = []
im.append(ax[0].imshow(data1[:,1:],extent=extent1,vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))
im.append(ax[1].imshow(data2,extent=extent2,vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))

ax[0].set_xlim(rhomin,rhomax)
ax[0].set_ylim(Tmin,Tmax)

#ax[0].yaxis.set_tick_params(which='both',direction='in')

pos = ax[-1].get_position()
cax = fig.add_axes([pos.x1+0.01,pos.y0,0.025,(pos.y1-pos.y0)])
cb = fig.colorbar(im[-1], cax=cax, orientation='vertical')

cax.set_ylabel("log opacity")

ax[0].set_ylabel("log T")
ax[0].set_xlabel("log $\\rho$")
ax[0].set_title(file1)
ax[1].set_title(file2)

"""
#logTs = [ [2.69897,3.041393],
#          [2.995635,4.0],
#          [3.75,8.7] ]
#logRs = [ [-12.999999960000000,-0.25],
#          [-9.0969100130080562,-1.0],
#          [-8,1.] ]

logTs = [0.]*8
logrhos = [0.]*8

logTs[0] = 2.69897
logTs[1] = logTs[0]
logTs[2] = 3.041393
logTs[3] = logTs[2]

logTs[4] = 2.995635
logTs[5] = logTs[4]
logTs[6] = 4.0
logTs[7] = logTs[6]

logrhos[0] = -19.+1.e-8
logrhos[1] = -7.+1.e-8
logrhos[2] = -19.+1.e-8
logrhos[3] = -7.+1.e-8
logrhos[4] = -19.+1.e-8
logrhos[5] = -7.+1.e-8
logrhos[6] = -19.+1.e-8
logrhos[7] = -7.+1.e-8

color = 'r'

#ax.plot([logrhos[0],logrhos[1]],[logTs[0],logTs[1]],color=color)
#ax.plot([logrhos[0],logrhos[2]],[logTs[0],logTs[2]],color=color)
#ax.plot([logrhos[1],logrhos[3]],[logTs[1],logTs[2]],color=color)
#ax.plot([logrhos[2],logrhos[3]],[logTs[2],logTs[3]],color=color)
#ax.plot([logrhos[3],logrhos[5]],[logTs[3],logTs[5]],color=color)
#ax.plot([logrhos[4],logrhos[5]],[logTs[4],logTs[5]],color=color)
#ax.plot([logrhos[5],logrhos[7]],[logTs[5],logTs[7]],color=color)
#ax.plot([logrhos[6],logrhos[7]],[logTs[6],logTs[7]],color=color)
#ax.plot([logrhos[6],logrhos[4]],[logTs[6],logTs[4]],color=color)

#ax.plot([-14.75,0.1],[3.75,8.7],color=color)
#ax.plot([-5.75,-14.75],[3.75,3.75],color=color)
#ax.plot([-5.75,9.1],[3.75,8.7],color=color)
#ax.plot([0.1,9.1],[8.7,8.7],color=color)
"""


plt.savefig("check_sanity_rho.png",dpi=300)

plt.show()
