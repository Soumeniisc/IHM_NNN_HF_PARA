import matplotlib.pyplot as plt
from scipy import *
import numpy as np	



t2 = 0.1
t = 0.5
lHF   = loadtxt("broyden_para2.dat")
f1 =  loadtxt("delta0.5t20.1_para.txt")


plt.plot(2*np.array(f1[:,0]),(-np.array(f1[:,2])-np.array(f1[:,3])+np.array(f1[:,4])+np.array(f1[:,5]))/2.0,'-*',label="Self_consistenc U inc")
plt.plot(2*np.array(lHF[:,0]),np.abs(np.array(lHF[:,5])),'-o',label="Broyden U inc")
#plt.plot([2.8,2.8],[0,0.6],"--k",label="broy_para_AFM en cross")
plt.xlabel("U/t")
plt.ylabel("dn")
plt.text(1.2,0.5,r"$\Delta = 1.0t, t2=%st$"%(t2/t),size = 20)
plt.legend()
plt.savefig("dn_t20.pdf")
plt.show()

print "U",np.array(f1[:,0])
print "dn",(-np.array(f1[:,2])-np.array(f1[:,3])+np.array(f1[:,4])+np.array(f1[:,5]))/2.0
print "mu",np.array(f1[:,1])






