import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numerical_functions import *

#variables


N_z = 81 #sections z axis (rows)
N_x = 41 #sections x axis (columns)
C_a = 0.505 #chord length aileron in [m]
l_a = 1.611 #span aileron in [m]

#start code =============================================================================================================
#import data
aerodata = np.genfromtxt('aerodynamicloadf100.dat', delimiter=',')  #aerodynamic loading data given

def q_disc(z): #discrete aerodynamic loading function over the span with input the discrete aerodynamic loading over the chord

    #trapezoidal integration chord wise
    q_x = np.array([])
    z_cp = np.array([])
    for j in range(N_x):
        q_z = np.array([])
        z_cp_j = np.array([])

        for i in range(N_z-1):
            s_i = (z[i] - z[i+1]) * ((aerodata[i+1,j] - aerodata[i,j])/2 + aerodata[i,j])
            q_z = np.append(q_z,s_i)
            z_cp_i = (q_z[i] * z[i])
            z_cp_j = np.append(z_cp_j,z_cp_i)

        q_x = np.append(q_x,sum(q_z))
        z_cp = np.append(z_cp, sum(z_cp_j)/sum(q_z))
# outputs q_x, the distribution of aerodynamic load along span in [kN/m], z_cp, location of center of pressure per cross section
    return q_x, z_cp

#test
#calculate the x and z coordinate mesh ===========================================================================

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

# #plotting for graphical representation ============================================================================

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


#interpolation and integration of aerodynamic loading =========================================================

q_x, z_cp = q_disc(zcoordinates)

#split up aerodynamic loading into n sections
x_list = []
n_sec = 150

#list of q values for n sections with plot
q_list = [] #value of q for every section

for i in range(n_sec+1):
    x = l_a / n_sec * i
    x_list = np.append(x_list, x)
    q_list = np.append(q_list, interpolate(xcoordinates,q_x,x))

plt.subplot(131)
plt.plot(x_list,q_list)
plt.title('q(x)')
plt.xlabel('Span [m]')
plt.ylabel('Aerodynamic load [kN/m]')
# plt.show()

#list of center of pressure locations for n sections with plot
z_cp_list = [] #value of q for every section

for i in range(1,n_sec+2):
    x = l_a / n_sec * i
    z_cp_list = np.append(z_cp_list, interpolate(xcoordinates,z_cp,x))

plt.subplot(133)
plt.plot(x_list,z_cp_list)
plt.title('z_cp(x)')
plt.xlabel('Span [m]')
plt.ylabel('Center of Pressure location on chord [m]')
plt.show()
