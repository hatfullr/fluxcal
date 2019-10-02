import numpy as np
import matplotlib.pyplot as plt

plt.style.use('supermongo_mnras')

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

R = np.loadtxt("gn93_z0.02_x0.98.data",skiprows=5,max_rows=1)
data = np.loadtxt("gn93_z0.02_x0.98.data",skiprows=6)
T = data[:,0]
data = data[:,1:]

rho = np.zeros((len(R),len(T)))
for i in range(0,len(rho)):
    for j in range(0,len(rho[i])):
        rho[i][j] = R[i] + 3.*T[j] - 18.

print(rho[0])
print(rho[1])

nrows = 1
ncols = 1
fig, ax = plt.subplots(nrows=nrows,ncols=ncols,sharex=False,sharey=False,dpi=300)

im = []
im.append(ax[0].imshow(data[:,1:],extent=(min(R),max(R),min(T),max(T)),vmin=vmin,vmax=vmax,cmap=my_cm,aspect='auto',origin='lower'))
