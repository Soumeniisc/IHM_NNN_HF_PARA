import matplotlib.pylab as plt
from scipy import *
import numpy as np


delta = 0.5
t = 0.5
#list_U = [0.0,0.5,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.5]
list_U = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,2.1]
list_symbol=['-s','-o']
for i,U in enumerate(list_U):
	f1 = loadtxt("dosdelta0.5U%st20.1.txt"%U)
	#f2 = loadtxt("delta%st20.1.txt"%name)
	plt.plot(2*np.array(f1[:,0]),np.array(f1[:,1])+np.array(f1[:,2]),'-o')
	#plt.plot(2*np.array(f1[:,0]),np.array(f1[:,3])+np.array(f1[:,4]),'-o',label=r'$\downarrow$')
	#plt.plot(2*np.array(f2[:,0]),f2[:,6])

	plt.text(-2.0,0.02,r'$t2=%st,\Delta = %st,U=%st$'%(0.0/t,delta/t,U/t),size=20)
	plt.xlabel("U/t")
	plt.ylabel("A(w)")
	plt.legend()
	plt.savefig("U_%st1.eps"%(U/t))
	plt.show()
	#plt.close()
