import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#Path
input_path1="/home/changwan/FDTD/Source_3D.txt"
input_path2="/home/changwan/FDTD/FDTD_3D_E.txt"
input_path3="/home/changwan/FDTD/FDTD_3D_E1.txt"
#input_path4="/home/changwan/FDTD/FDTD_3D_H.txt"
#input_path5="/home/changwan/FDTD/FDTD_3D_H1.txt"


#Laoding

# 1=x, 2=y, 3=z

#Source
S1=np.loadtxt(input_path1,usecols=0)
S2=np.loadtxt(input_path1,usecols=1)
S3=np.loadtxt(input_path1,usecols=2)

#"""
#E = Final output
E1=np.loadtxt(input_path2, usecols=0)
E2=np.loadtxt(input_path2, usecols=1)
E3=np.loadtxt(input_path2, usecols=2)


#E1 = middle output
E11=np.loadtxt(input_path3, usecols=0)
E12=np.loadtxt(input_path3, usecols=1)
E13=np.loadtxt(input_path3, usecols=2)
#"""

"""
#H = Final output
H1=np.loadtxt(input_path4, usecols=0)
H2=np.loadtxt(input_path4, usecols=1)
H3=np.loadtxt(input_path4, usecols=2)


#H1 = middle output
H11=np.loadtxt(input_path5, usecols=0)
H12=np.loadtxt(input_path5, usecols=1)
H13=np.loadtxt(input_path5, usecols=2)
"""


#plot the graphs
plt.subplot(2,2,1)
plt.plot(S1,'ro',label='Source_x')
#plt.plot(S2,'b-',label='Source_y')
#plt.plot(S3,'g+',label='Source_z')
#plt.title("Source")
plt.grid()
plt.minorticks_on()

#"""
plt.subplot(2,2,2)
plt.plot(E1,'ro',label='E_x_field')
#plt.plot(E2,'b',label='E_y_field')
#plt.plot(E3,'g',label='E_z_field')
plt.title("First time step")
plt.grid()
plt.minorticks_on()

plt.subplot(2,2,3)
plt.plot(E11,'r',label='E_x_field')
#plt.plot(E12,'b-',label='E_y_field')
#plt.plot(E13,'g+',label='E_z_field')
plt.title("Final time step")
plt.grid()
plt.minorticks_on()
#"""

"""
plt.subplot(2,2,2)
plt.plot(H1,'r',label='H_x_field')
plt.plot(H2,'b',label='H_y_field')
plt.plot(H3,'g',label='H_z_field')
plt.grid()
plt.minorticks_on()

plt.subplot(2,2,3)
plt.plot(H11,'r',label='H_x_field')
plt.plot(H12,'b',label='H_y_field')
plt.plot(H13,'g',label='H_z_field')
plt.title("Middle time step")
plt.grid()
plt.minorticks_on()


#plt.subplot(2,2,4)
#plt.plot(E1[18000:20000],'r',label='E_z_field')
#plt.title("E_z_field [18000:20000]")
#plt.grid()
#plt.minorticks_on()
#=====================================
"""

#plt.grid()
#plt.plot(z,M,'b')

plt.show()
