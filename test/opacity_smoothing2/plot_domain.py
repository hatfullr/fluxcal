import numpy as np
import matplotlib.pyplot as plt

plt.style.use('supermongo_mnras')


data = np.loadtxt('test_domain.txt')
origTs = data[:,0]
origRs = data[:,1]
newTs = data[:,2]
newRs = data[:,3]

origTs = origTs[np.where(origTs > -1.e30)[0]]
origRs = origRs[np.where(origRs > -1.e30)[0]]
newTs = newTs[np.where(newTs > -1.e30)[0]]
newRs = newRs[np.where(newRs > -1.e30)[0]]


fig, ax = plt.subplots(nrows=1,ncols=2,sharey=True)
plt.subplots_adjust(top=0.8,bottom=0.2,right=0.9)

s1=10
s2=1
marker1 = '.'
marker2 = '*'
color1 = 'k'
color2 = 'r'
ax[0].scatter(np.arange(len(newTs)),newTs,color=color2,s=s2,marker=marker2)
ax[0].scatter(np.arange(len(origTs)),origTs,color=color1,s=s1,marker=marker1)
ax[1].scatter(np.arange(len(origRs)),origRs,color=color1,s=s1,marker=marker1)
ax[1].scatter(np.arange(len(newRs)),newRs,color=color2,s=s2,marker=marker2)




plt.show()
