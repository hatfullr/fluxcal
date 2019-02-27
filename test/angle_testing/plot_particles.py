import matplotlib.pyplot as plt
import numpy as np

plt.style.use('supermongo')

data = np.loadtxt("angle_testing.dat")
x = data[:,0]
y = data[:,1]
z = data[:,2]
myx = data[:,3]
myy = data[:,4]
myz = data[:,5]


fig, ax = plt.subplots(nrows=1,ncols=2)
ax[0].scatter(x,y,color='k',s=1)
ax[1].scatter(myx,myy,color='k',s=1)

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
plt.subplots_adjust(wspace=0.2)
plt.show()
