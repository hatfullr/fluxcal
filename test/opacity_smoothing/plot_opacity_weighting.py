import numpy as np
import matplotlib.pyplot as plt

plt.style.use('supermongo')

data = np.loadtxt("testing.dat")
logT = data[:,0]
weight1 = data[:,1]
weight2 = data[:,2]


plt.scatter(logT,weight1,color='k')
plt.scatter(logT,weight2,color='r')
#plt.scatter(logT,opac_low,color='r')
#plt.scatter(logT,opac_high,color='b')

plt.show()
