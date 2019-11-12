import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.style.use('supermongo_mnras')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

R1 = np.loadtxt("test_opacity_file.txt",max_rows=1)
data1 = np.loadtxt("test_opacity_file.txt",skiprows=1)
T1 = data1[:,0]

#R2 = np.loadtxt("gn93_z0.02_x0.98.data",skiprows=5,max_rows=1)
#data2 = np.loadtxt("gn93_z0.02_x0.98.data",skiprows=6)
R2 = np.loadtxt("test_opacity_file1.txt",max_rows=1)
data2 = np.loadtxt("test_opacity_file1.txt",skiprows=1)
T2 = data2[:,0]

#vmin = min(np.amin(data1[np.isfinite(data1)]), np.amin(data2[np.isfinite(data2)]))
#vmax = max(np.amax(data1[np.isfinite(data1)]), np.amax(data2[np.isfinite(data2)]))

vmin = -5
vmax = 10

Tmin = max(np.amin(T1),np.amin(T2))
Tmax = min(np.amax(T1),np.amax(T2))
Rmin = max(np.amin(R1),np.amin(R2))
Rmax = min(np.amax(R1),np.amax(R2))

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

"""

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=False,sharey=False,dpi=300)
plt.subplots_adjust(bottom=0.125,top=0.9,left=left,wspace=wspace)
Rcomp_pos = int(len(R1)*0.1)
Tcomp_pos = int(len(T1)*0.5)
Rcompare_idx = find_nearest(R2,R1[Rcomp_pos])
Tcompare_idx = find_nearest(T2,T1[Tcomp_pos])

ax[0].plot(T1,data1[:,Rcomp_pos],color='k',label="My file",zorder=1.e30)
ax[0].plot(T2,data2[:,Rcompare_idx],color='r',label="Orig file",lw=1.5)
ax[0].legend()

ax[0].set_xlabel("logT")
ax[0].set_ylabel("log opacity")
ax[0].set_title("log R = %7.2E" % (R1[Rcomp_pos]))

ax[1].plot(R1,data1[Tcomp_pos][1:],color='k',label="My file",zorder=1.e30)
ax[1].plot(R2,data2[Tcompare_idx][1:],color='r',label="Orig file",lw=1.5)
ax[1].legend()

ax[1].set_xlabel("logR")
ax[1].set_ylabel("log opacity")
ax[1].set_title("log T = %7.2E" % (T1[Tcomp_pos]))

"""

fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,dpi=300)
plt.subplots_adjust(wspace=wspace,left=left,right=0.8,top=0.95,bottom=0.15)
im = []
im.append(ax.imshow(data1[:,1:],extent=(min(R1),max(R1),min(T1),max(T1)),vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))
#im.append(ax[1].imshow(data2[:,1:],extent=(min(R2),max(R2),min(T2),max(T2)),vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))


ax.yaxis.set_tick_params(which='both',direction='in')

pos = ax.get_position()
cax = fig.add_axes([pos.x1+0.01,pos.y0,0.05,(pos.y1-pos.y0)])
cb = fig.colorbar(im[-1], cax=cax, orientation='vertical')

cax.set_ylabel("log opacity")

ax.set_ylabel("log T")
ax.set_xlabel("log R")
#ax[1].set_xlabel("log R")

#ax[0].set_title("My file")
#ax[1].set_title("Original file")





#logTs = [ [2.69897,3.041393],
#          [2.995635,4.0],
#          [3.75,8.7] ]
#logRs = [ [-12.999999960000000,-0.25],
#          [-9.0969100130080562,-1.0],
#          [-8,1.] ]

logTs = [0.]*8
logRs = [0.]*8

logTs[0] = 2.69897
logTs[1] = logTs[0]
logTs[2] = 3.041393
logTs[3] = logTs[2]

logTs[4] = 2.995635
logTs[5] = logTs[4]
logTs[6] = 4.0
logTs[7] = logTs[6]

logRs[0] = (-19.+1.e-8) - 3.*logTs[0] + 18.
logRs[1] = (-7.+1.e-8) - 3.*logTs[1] + 18.
logRs[2] = (-19.+1.e-8) - 3.*logTs[2] + 18.
logRs[3] = (-7.+1.e-8) - 3.*logTs[3] + 18.
logRs[4] = (-19.+1.e-8) - 3.*logTs[4] + 18.
logRs[5] = (-7.+1.e-8) - 3.*logTs[5] + 18.
logRs[6] = (-19.+1.e-8) - 3.*logTs[6] + 18.
logRs[7] = (-7.+1.e-8) - 3.*logTs[7] + 18.

color = 'r'

ax.plot([logRs[0],logRs[1]],[logTs[0],logTs[1]],color=color)
ax.plot([logRs[0],logRs[2]],[logTs[0],logTs[2]],color=color)
ax.plot([logRs[1],logRs[3]],[logTs[1],logTs[2]],color=color)
ax.plot([logRs[2],logRs[3]],[logTs[2],logTs[3]],color=color)
ax.plot([logRs[3],logRs[5]],[logTs[3],logTs[5]],color=color)
ax.plot([logRs[4],logRs[5]],[logTs[4],logTs[5]],color=color)
ax.plot([logRs[5],logRs[7]],[logTs[5],logTs[7]],color=color)
ax.plot([logRs[6],logRs[7]],[logTs[6],logTs[7]],color=color)
ax.plot([logRs[6],logRs[4]],[logTs[6],logTs[4]],color=color)

ax.plot([-8.,1.],[3.75,3.75],color=color)
ax.plot([-8.,-8.],[3.75,8.7],color=color)
ax.plot([1.,1.],[3.75,8.7],color=color)
ax.plot([-8.,1.],[8.7,8.7],color=color)

#for i in range(1,len(logRs)):
#    ax.plot([logRs[i-1],logRs[i]],[logTs[i-1],logTs[i]],color=color)
#    ax.plot([logRs[i][0],logRs[i][0]],[logTs[i][0],logTs[i][1]],color=color)
#    ax.plot([logRs[i][0],logRs[i][1]],[logTs[i][1],logTs[i][1]],color=color)
#    ax.plot([logRs[i][1],logRs[i][1]],[logTs[i][1],logTs[i][0]],color=color)
#    ax.plot([logRs[i][1],logRs[i][0]],[logTs[i][0],logTs[i][0]],color=color)

plt.savefig("check_sanity.png",dpi=300)

plt.show()
