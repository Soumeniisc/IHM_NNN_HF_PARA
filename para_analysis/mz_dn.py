import matplotlib.pyplot as plt
from scipy import *
import numpy as np	

t2    =  0.0
t     =  0.5
delta =  0.5
l     =  loadtxt("IPT_t20.0.dat")
lHF   =  loadtxt("broyden_para_t20.0.dat")
lCTQMC=   loadtxt("CTQMC_t20.0.dat")
lCTQMC_triqs = loadtxt("IHM_Uall_beta100_delta0.5_t20.0_mu0.0_triqs1.0.dat")
plt.plot(2*np.array(lHF[:,0]),np.array(lHF[:,5]),'-o',label="HF")
plt.plot(2*np.array(l[:,0]),np.array(l[:,7]),'-o',label="IPT")
plt.plot(np.array(lCTQMC[:,0]),0.5*(np.array(lCTQMC[:,3])-np.array(lCTQMC[:,2])),'-^',label="CTQMC")
plt.plot(np.array(lCTQMC_triqs[:,0]),np.array(lCTQMC_triqs[:,4]),'-*',label="CTQMC_triqs")
plt.text(2.0,0.5,r'$t2=%st,\Delta = %st$'%(0.0/t,delta/t),size=20)
plt.xlabel(r'$U/t$',size=20)
plt.ylabel(r'$\delta n$',size=20)
plt.legend()
plt.savefig("dn_vs_U_HFvsIPT_t2%st.eps"%(t2*2))

plt.show()


