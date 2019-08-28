import numpy as np
import matplotlib.pyplot as plt

plt.style.use('supermongo')

data = np.loadtxt("testing.dat")
logR = data[:,0]
logT = data[:,1]
opac_low = np.log10(data[:,2])
opac_high = np.log10(data[:,3])
opacity = np.log10(data[:,4])


plt.scatter(logT,opacity,color='k')
#plt.scatter(logT,opac_low,color='r')
#plt.scatter(logT,opac_high,color='b')

plt.show()
