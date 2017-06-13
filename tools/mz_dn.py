import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


lHF   = loadtxt("broyden_para_.dat")
f1 =  loadtxt("delta0.5t20_para.txt")


plt.plot(2*np.array(f1[:,0]),(-np.array(f1[:,2])-np.array(f1[:,3])+np.array(f1[:,4])+np.array(f1[:,5]))/2.0,'-*',label="Self_consistenc U inc")
plt.plot(2*np.array(lHF[:,0]),np.abs(np.array(lHF[:,5])),'-o',label="Broyden U dec")
#plt.plot([2.8,2.8],[0,0.6],"--k",label="broy_para_AFM en cross")
plt.xlabel("U/t")
plt.ylabel("dn")
plt.text(1.2,0.5,r"$\Delta = 0.5t, t2=0.0t$",size = 20)
plt.legend(loc=2)
plt.savefig("dn_t20.pdf")
plt.show()


lHF   = loadtxt("broyden_para_.dat")
f1 =  loadtxt("delta0.5t20_para.txt")


plt.plot(2*np.array(f1[:,0]),(-np.array(f1[:,2])-np.array(f1[:,3])+np.array(f1[:,4])+np.array(f1[:,5]))/2.0,'-*',label="Self_consistenc U inc")
plt.plot(2*np.array(lHF[:,0]),np.abs(np.array(lHF[:,5])),'-o',label="Broyden U dec")
#plt.plot([2.8,2.8],[0,0.6],"--k",label="broy_para_AFM en cross")
plt.xlabel("U/t")
plt.ylabel("dn")
plt.text(1.2,0.5,r"$\Delta = 0.5t, t2=0.0t$",size = 20)
plt.legend(loc=2)
plt.savefig("dn_t20.pdf")
plt.show()



