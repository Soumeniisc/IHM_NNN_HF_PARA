import matplotlib.pyplot as plt
from scipy import *
import numpy as np


ax = [plt.subplot(2,1,i+1) for i in range(2)]
t = 0.5
t2 = 0.2
delta = 0.5
U = 0.1
f1 = loadtxt("dosdelta0.5U%st2%s.txt"%(U,t2))
ax[0].plot(2*np.array(f1[:,0]),np.array(f1[:,1])+np.array(f1[:,2]),'-',label=r'$U=%st$'%(U/t))
ax[0].set_xticklabels([])
ax[0].set_yticklabels([])
ax[0].legend(loc=2)
#plt.xlim(-8,3)
ax[0].set_xlim([-4, 8])
U = 1
f1 = loadtxt("dosdelta0.5U%st2%s.txt"%(U,t2))
ax[1].plot(2*np.array(f1[:,0]),np.array(f1[:,1])+np.array(f1[:,2]),'-',label=r'$U=%st$'%(U/t))
#ax[1].set_xticklabels([])
ax[1].set_yticklabels([])
plt.ylabel("A(w)")
ax[1].legend(loc=2)
ax[1].set_xlim([-4, 8])
U = 3.3
'''
f1 = loadtxt("dosdelta0.5U%st2%sh.txt"%(U,t2))
ax[2].plot(2*np.array(f1[:,0]),np.array(f1[:,1])+np.array(f1[:,2]),'-o',label=r'$U=%st$'%(U/t))
plt.subplots_adjust(wspace=0, hspace=0)
ax[2].set_yticklabels([])

ax[2].legend()
'''
plt.subplots_adjust(wspace=0, hspace=0)
plt.xlabel("w/t")
plt.suptitle("dos for "+r'$\Delta = %st, t2=%st$'%(delta/t,t2/t)+" in HF theory",size =20)
#plt.xlim(-8,3)
plt.savefig("dos_for_delta%st_t2%st.eps"%(delta/t,t2/t))
plt.show()
'''
	plt.text(-2.0,0.02,r'$t2=%st,\Delta = %st,U=%st$'%(0.0/t,delta/t,U/t),size=20)
	plt.xlabel("U/t")
	plt.ylabel("A(w)")
	plt.legend()
	plt.savefig("U_%st.eps"%(U/t))
	plt.show()
'''
