import matplotlib.pyplot as plt
from scipy import *
import numpy as np	


t = 0.5
'''

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
'''

U=0.1
f1 = loadtxt("energy_vs_dn_U%s.dat"%U)
#plt.plot(np.array(f1[:,0]),f1[:,1],'-o',label="U=%st"%(U/t))
U=0.1
f1 = loadtxt("energy_vs_dn_U%s.dat"%U)
plt.plot(np.array(f1[:,0]),np.array(f1[:,1])-2*np.array(f1[:,2]),'-o',label="kinetic_U=%st"%(U/t))
plt.plot(np.array(f1[:,0]),np.array(f1[:,2]),'-*',label="potetnial_U=%st"%(U/t))
plt.plot(np.array(f1[:,0]),np.array(f1[:,3]),'-^',label="total_U=%st"%(U/t))
U=2
#f1 = loadtxt("energy_vs_dn_U%s_.dat"%U)
#plt.plot(np.array(f1[:,0]),f1[:,1],'-o',label="U=%st"%(U/t))
plt.xlabel("dn")
plt.ylabel("energy")
plt.text(-0.7,0.5,r"$\Delta = 1.0t, t2=0.2t$",size = 20)
plt.legend(loc=2)
plt.savefig("energy_vs_dn_t202t_U%s.pdf"%U)
plt.show()



