import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.gca(projection = '3d')


#FD-TD results
#x = [3,2,2,1,2,2]
#y = [4,5,4,4,3,4]
#z = [4,4,5,4,4,3]


#E_x
x1 = [3,3,3,3]
y1 = [5,3,4,4]
z1 = [4,4,3,5]
ax.scatter(x1,y1,z1)
#ax.scatter(3,4,4)


#E_yax.scatter(3,4,4)
x2 = [2,2,1,3]
y2 = [5,5,5,5]
z2 = [5,3,4,4]
ax.scatter(x2,y2,z2)
#ax.scatter(2,5,4)


#E_z
x3 = [3,1,2,2]
y3 = [4,4,3,5]
z3 = [5,5,5,5]
ax.scatter(x3,y3,z3)
#ax.scatter(2,4,5)



#H_x
x4 = [1,1,1,1]
y4 = [4,4,3,5]
z4 = [5,3,4,4]
ax.scatter(x4,y4,z4)
#ax.scatter(1,4,4)

#H_y
x5 = [3,1,2,2]
y5 = [3,3,3,3]
z5 = [4,4,3,5]
ax.scatter(x5,y5,z5)
#ax.scatter(2,3,4)

#H_z
x6 = [2,2,1,3]
y6 = [5,3,4,4]
z6 = [3,3,3,3]
ax.scatter(x6,y6,z6)
#ax.scatter(2,4,3)


ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')



plt.show()
