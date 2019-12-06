import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('supermongo_mnras')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

rho1 = np.loadtxt("test_opacity_file.txt",max_rows=1)
data1 = np.loadtxt("test_opacity_file.txt",skiprows=1)
T1 = data1[:,0]

#rho2 = np.loadtxt("test_opacity_file1.txt",max_rows=1)
#data2 = np.loadtxt("test_opacity_file1.txt",skiprows=1)
#T2 = data2[:,0]

#vmin = min(np.amin(data1[np.isfinite(data1)]), np.amin(data2[np.isfinite(data2)]))
#vmax = max(np.amax(data1[np.isfinite(data1)]), np.amax(data2[np.isfinite(data2)]))

vmin = -5
vmax = 10

Tmin = np.amin(T1)
Tmax = np.amax(T1)
rhomin = np.amin(rho1)
rhomax = np.amax(rho1)

#dv = vmax - vmin
#buff = 1.
#vmin += dv*buff
#vmax -= dv*buff

my_cm = plt.cm.get_cmap('nipy_spectral')

my_cm.set_bad('white')

wspace = 0.3
left = 0.15

nrows = 1
ncols = 1

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,dpi=300)
plt.subplots_adjust(wspace=wspace,left=left,right=0.8,top=0.95,bottom=0.15)
im = []
im.append(ax.imshow(data1[:,1:],extent=(min(rho1),max(rho1),min(T1),max(T1)),vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))



ax.yaxis.set_tick_params(which='both',direction='in')

pos = ax.get_position()
cax = fig.add_axes([pos.x1+0.01,pos.y0,0.05,(pos.y1-pos.y0)])
cb = fig.colorbar(im[-1], cax=cax, orientation='vertical')

cax.set_ylabel("log opacity")

ax.set_ylabel("log T")
ax.set_xlabel("log $\\rho$")


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

ax.plot([logrhos[0],logrhos[1]],[logTs[0],logTs[1]],color=color)
ax.plot([logrhos[0],logrhos[2]],[logTs[0],logTs[2]],color=color)
ax.plot([logrhos[1],logrhos[3]],[logTs[1],logTs[2]],color=color)
ax.plot([logrhos[2],logrhos[3]],[logTs[2],logTs[3]],color=color)
ax.plot([logrhos[3],logrhos[5]],[logTs[3],logTs[5]],color=color)
ax.plot([logrhos[4],logrhos[5]],[logTs[4],logTs[5]],color=color)
ax.plot([logrhos[5],logrhos[7]],[logTs[5],logTs[7]],color=color)
ax.plot([logrhos[6],logrhos[7]],[logTs[6],logTs[7]],color=color)
ax.plot([logrhos[6],logrhos[4]],[logTs[6],logTs[4]],color=color)

ax.plot([-14.75,0.1],[3.75,8.7],color=color)
ax.plot([-5.75,-14.75],[3.75,3.75],color=color)
ax.plot([-5.75,9.1],[3.75,8.7],color=color)
ax.plot([0.1,9.1],[8.7,8.7],color=color)

plt.savefig("check_sanity_rho.png",dpi=300)

plt.show()
