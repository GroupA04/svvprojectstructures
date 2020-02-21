import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numerical_functions import *

#variables
from svvprojectstructures.numerical_functions import interpolate

N_z = 81 #sections z axis (rows)
N_x = 41 #sections x axis (columns)
C_a = 0.505 #chord length aileron in [m]
l_a = 1.611 #span aileron in [m]

#start code ----------------------------------------------------------------------------------------------------------------
#import data
aerodata = np.genfromtxt('aerodynamicloadf100.dat', delimiter=',')  #aerodynamic loading data given

def q_disc(z): #discrete aerodynamic loading function over the span with input the discrete aerodynamic loading over the chord

    #trapezoidal integration chord wise
    q_x = np.array([])
    for j in range(N_x):
        q_z = np.array([])
        for i in range(N_z-1):
            s_i = (z[i+1] - z[i]) * (aerodata[i+1,j] - aerodata[i,j])/2
            q_z = np.append(q_z,s_i)

        q_x = np.append(q_x,sum(q_z))

    return q_x


#calculate the x and z coordinate mesh

zcoordinates = np.zeros([N_z])
xcoordinates = np.zeros([N_x])

for i in range(1, N_z + 1 + 1):
    theta_zi = (i-1)/N_z *np.pi
    theta_zi1 = (i+1-1)/N_z *np.pi
    z_i = -1/2 * (C_a / 2*(1-np.cos(theta_zi)) + C_a/2 *(1-np.cos(theta_zi1)))
    if i <N_z +1:
        zcoordinates[i-1] = z_i


for i in range(N_x + 1 + 1):
    theta_xi = (i-1)/N_x *np.pi
    theta_xi1 = (i - 1 + 1) / N_x * np.pi
    x_i = 1/2 * (l_a / 2*(1-np.cos(theta_xi)) + l_a/2 *(1-np.cos(theta_xi1)))
    if i <N_x +1:
        xcoordinates[i-1] = x_i

# #plotting for graphical representation

# xcoordinates,zcoordinates = np.meshgrid(xcoordinates,zcoordinates)
#
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# ax.plot_surface(xcoordinates, zcoordinates, aerodata)
# ax.set_xlabel('Spanwise axis [m]')
# ax.set_ylabel('Chordwise axis [m]')
# ax.set_zlabel('Aerodynamic loading [kN/m2]')
# plt.show()



#split up aerodynamic loading into n sections

x_list = []
q_list = []
n_sec = 100
for i in range(n_sec+1):
    x = l_a / n_sec * i
    #x = (l_a - 2 * (l_a / 100)) / n_sec * i
    x_list = np.append(x_list, x)
    q_list = np.append(q_list, interpolate(xcoordinates,q_disc(zcoordinates),x))

plt.plot(x_list,q_list)
plt.xlabel('Span [m]')
plt.ylabel('Aerodynamic load [kN/m]')
plt.show()