import matplotlib.pylab as plt
from scipy import *
import numpy as np

t= 0.5
delta = 0.5
list_v = [0.0,0.1]#0.2,0.3,0.4,0.5]
#list_v = [0.3,0.35,0.4,0.425,0.44]
list_f = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
#list_v = [0,0.1,0.2,0.25,0.3]
list_f = [0.5,0.5,0.5,0.5,0.5,0.5,0.5]
list_symbol=['-s','-^','-^','-*','-8','-o']
color_list = ['b','r','k']
for i,name in enumerate(list_v):
	f1 = loadtxt("broyden_para_t2%s.dat"%name)
	#f2 = loadtxt("delta%st20.1.txt"%name)
	#plt.plot(2*np.array(f1[:,0]),f1[:,6],list_symbol[i],label=r'$t2=%st$'%(list_v[i]/t))
	#plt.plot(2*np.array(f2[:,0]),f2[:,6])
	plt.plot(2*np.array(f1[:,0]),np.array(f1[:,5]),list_symbol[i],color=color_list[i],label=r'$t2=%st$'%(list_v[i]/t))
	#dn = (-np.array(f1[:,2])-np.array(f1[:,3])+np.array(f1[:,4])+np.array(f1[:,5]))/2.0
	#plt.plot(2*np.array(f1[:,0]),delta-np.array(f1[:,0])*dn/2.0,list_symbol[i],label=r'$t2=%st$'%(list_v[i]/t))

plt.text(0.5,0.3,r'$\Delta=%st$'%(list_f[i]/t),size=20)
plt.xlabel("U/t")
plt.ylabel("dn")
#plt.legend(loc=2)
plt.legend()
plt.savefig("dn_last_Delta_%stsmall.eps"%(list_f[i]/t))
#plt.savefig("ms_continuaous_Delta_%st.eps"%(list_f[i]/t))
plt.show()
	