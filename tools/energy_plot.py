import matplotlib.pylab as plt
from scipy import *
import numpy as np


t2 = 0.1

file1 = "brouden_AFM_energy.dat"
file2 = "brouden_energy.dat"
f1 = loadtxt(file1)
#f3 = loadtxt("energy_delta0.5t2%s_reverse2.txt"%t2)
f2 = loadtxt(file2)
plt.plot(2*np.array(f1[:,0]),f1[:,1],'-o',label="AFM")
plt.plot(2*np.array(f2[:,0]),f2[:,1],"-o",label="para")
#plt.plot(2*np.array(f3[:,0]),f3[:,1],"-o",label="reverse")

plt.text(0.5,0.3,r'$\Delta=%st$'%(1.0),size=20)
plt.xlabel("U/t")
plt.ylabel("energy")
#plt.legend(loc=2)
plt.legend()
plt.savefig("energy_for_para_AFM.eps")
#plt.savefig("ms_continuaous_Delta_%st.eps"%(list_f[i]/t))
plt.show()
	
