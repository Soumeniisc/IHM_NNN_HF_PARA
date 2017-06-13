import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

t2    =  0.0
t     =  0.5
delta =  0.5

lCTQMC=   loadtxt("CTQMC_t20.0.dat")
lHF   =  loadtxt("broyden_para_t20.0.dat")
t2_list = [0.0,0.1,0.2]
for t2 in t2_list:
	l     =  loadtxt("IPT_t2%s.dat"%t2)
	#plt.plot(2*np.array(lHF[:,0]),np.array(lHF[:,5]),'-o',label="HF")
	plt.plot(2*np.array(l[:,0]),np.array(l[:,7]),'-o',label="t2=%s"%t2)
	#plt.plot(np.array(lCTQMC[:,0]),0.5*(np.array(lCTQMC[:,3])-np.array(lCTQMC[:,2])),'-^',label="CTQMC")
#plt.text(2.0,0.5,r'$t2=%st,\Delta = %st$'%(0.0/t,delta/t),size=20)
plt.xlabel(r'$U/t$',size = 20)
plt.ylabel(r'$\delta n$',size=20)
plt.legend()
plt.savefig("dn_vs_UIPT.eps")

plt.show()
