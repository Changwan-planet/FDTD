import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#Path
input_path1="/home/changwan/FDTD/FDTD_1D_E.txt"
input_path2="/home/changwan/FDTD/FDTD_1D_H.txt"
input_path3="/home/changwan/FDTD/Source_1D.txt"
input_path4="/home/changwan/FDTD/FDTD_1D_E1.txt"
input_path5="/home/changwan/FDTD/FDTD_1D_H1.txt"



#Laoding

#EM=np.genfromtxt(input_path,dtype=np.float)

E1=np.loadtxt(input_path1, usecols=0)
#E2=np.loadtxt(input_path1, usecols=1)
#E3=np.loadtxt(input_path1, usecols=2)

E11=np.loadtxt(input_path4, usecols=0)



M1=np.loadtxt(input_path2, usecols=0)
#M2=np.loadtxt(input_path2, usecols=1)
#M3=np.loadtxt(input_path2, usecols=2)

M11=np.loadtxt(input_path5, usecols=0)

"""
#=======================================
fig = plt.figure()
ax = fig.gca(projection='3d')

z1=E1
#z2=E2
#z3=E3
z11=E11

y1=M1
#y2=M2
#y3=M3
y11=M11

x=np.arange(len(y1))
w=np.zeros(len(y1))


#ax.plot(x,w,z1,'r',label='E_z_field')
#ax.plot(x[0:10000],w[0:10000],z1[0:10000],'r',label='E_z_field')
#ax.plot(x[18000:20000],w[18000:20000],z1[18000:20000],'r',label='E_z_field')
#ax.plot(x[18000:20000],y1[18000:20000],w[18000:20000],'b',label='H_y_field')
#ax.plot(x,w,z11,'g')
#ax.plot(x,y11,w,'g')
ax.grid()
ax.minorticks_on()
ax.tick_params(which="both", width=2)
ax.tick_params(which="major",length=7)
ax.tick_params(which="minor",length=4,color='r')
ax.legend()


#ax.set_xlabel("x [time]",fontsize=20)
#ax.set_zlabel("z [E_z]",fontsize=20)
#ax.set_ylabel("y [H_y]",fontsize=20)
#=====================================
"""

#"""
#====================================



S1=np.loadtxt(input_path3,usecols=0)

#ax.plot(x,w,S1,'g',label='Source')

plt.subplot(2,2,1)
plt.plot(S1,'b',label='Source')
#plt.title("Source")
plt.grid()
plt.minorticks_on()


plt.subplot(2,2,2)
plt.plot(E1,'r',label='E_z_field')
plt.title("E_z_field")
plt.grid()
plt.minorticks_on()

plt.subplot(2,2,3)
plt.plot(E1[0:10000],'r',label='E_z_field')
plt.title("E_z_field [0:10000]")
plt.grid()
plt.minorticks_on()

plt.subplot(2,2,4)
plt.plot(E1[18000:20000],'r',label='E_z_field')
plt.title("E_z_field [18000:20000]")

plt.grid()
plt.minorticks_on()
#=====================================
#"""

#plt.grid()
#plt.plot(z,M,'b')

plt.show()
