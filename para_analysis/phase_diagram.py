import matplotlib.pylab as plt
from scipy import *
import numpy as np

t=0.5
delta = 0.5
U_metal_t2_list    = [0.03,0.04, 0.05, 0.075, 0.1, 0.15, 0.2,  0.25,  0.3]
U_metal            = [4.3, 3.5, 2.9, 2.1,   1.5, 0.8,  0.4,   0.0,   0.0] #there is no metalic phase for t20.025 for U upto 5.5 

U_metal_t2_listIPT    = [ 0.1,  0.2]
U_metalIPT            = [ 1.2, 0.4]
U_mottIPT             = [ 1.2,   1.5, 0.8,  0.4,  0.0,   0.0,  0.0,   0.0,  0.0]

U_metal_t2_listCTQMC    = [ 0.0, 0.1,  0.2]
U_metalCTQMC            = [ 1.2, 0.4]
U_mottCTQMC             = [ 5.0,   5.2, 5.3]

plt.plot([0.0,0.0],[2.5,10.0],'-',color='r',linewidth=6.0)
plt.plot([0.5,0.6],[0.0,0.0],'-',color='r',linewidth=6.0)
plt.plot(np.array(U_metal_t2_list)/t,np.array(U_metal)/t,'-o',markersize=8,label=r'$HF$',)
plt.text(0.3,5.0,r'$Metal$',size=20)
plt.text(0.1,1.5,r'$BI$',size=20)
#plt.plot(np.array(U_c_t2_list)/t,np.array(U_c)/t,'-*',label=r'$U_c$')
plt.plot(np.array(U_metal_t2_listIPT)/t,np.array(U_metalIPT)/t,'-^',color='y',label=r'$IPT$')
plt.plot(np.array(U_metal_t2_listCTQMC)/t,np.array(U_mottCTQMC)/t,'-s',color='m',label=r'$CTQMC$')
plt.text(0.3,11.5,r'$MI$',size=20)

#plt.plot(np.array(U_c_t2_list)/t,np.array(U_c)/t,'-s')
plt.text(0.3,6.5,r'$\Delta = %st$'%(delta/t),size=20)
plt.title("square lattice IHM para-magnetic phase digram with NNN hopping")
plt.xlabel(r'$t_2/t$',size = 20)
plt.ylabel(r'$U/t$',size = 20)
plt.ylim(0.0,13)
plt.legend()
plt.savefig("phase_diagram_delta_%st_.eps"%(delta/t))
plt.show()
