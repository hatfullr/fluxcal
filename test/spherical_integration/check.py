import numpy as np
import matplotlib.pyplot as plt

plt.style.use('supermongo')

data = np.loadtxt("data.dat")

anglexdeg = data[:,0]
angleydeg = data[:,1]
anglezdeg = data[:,2]
x = data[:,3]
y = data[:,4]
z = data[:,5]

fig, ax = plt.subplots(1,3,figsize=(12,4))

n=8

ax[0].scatter(x[:n],y[:n],s=100,color='k')
ax[1].scatter(x[:n],z[:n],s=100,color='k')
ax[2].scatter(y[:n],z[:n],s=100,color='k')

#ax[0].scatter(x[3:7],y[3:7],s=10,color='r')
#ax[1].scatter(x[3:7],z[3:7],s=10,color='r')

ax[0].set_xlim(-1.2,1.2)
ax[0].set_ylim(-1.2,1.2)
ax[1].set_xlim(-1.2,1.2)
ax[1].set_ylim(-1.2,1.2)
ax[2].set_xlim(-1.2,1.2)
ax[2].set_ylim(-1.2,1.2)
ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

plt.show()
