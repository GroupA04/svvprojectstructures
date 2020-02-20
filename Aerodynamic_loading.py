import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate

#variables

N_z = 81 #sections z axis (rows)
N_x = 41 #sections x axis (columns)
C_a = 0.505 #chord length aileron in [m]
l_a = 1.611 #span aileron in [m]

#start code ----------------------------------------------------------------------------------------------------------------
#import data
aerodata = np.genfromtxt('aerodynamicloadf100.dat', delimiter=',')  #aerodynamic loading data given

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
#
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



def q(x):

    #variables

    N_z = 81 #sections z axis (rows)
    N_x = 41 #sections x axis (columns)
    C_a = 0.505 #chord length aileron in [m]
    l_a = 1.611 #span aileron in [m]


    #import data
    aerodata = np.genfromtxt('aerodynamicloadf100.dat', delimiter=',')  #aerodynamic loading data given

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


    #trapezoidal integration chord wise


    q_x = np.array([])
    for j in range(N_x):
        q_z = np.array([])
        for i in range(N_z-1):
            s_i = (zcoordinates[i+1] - zcoordinates[i]) * (aerodata[i+1,j] - aerodata[i,j])/2
            q_z = np.append(q_z,s_i)

        q_x = np.append(q_x,sum(q_z))

    #linear spline for span wise interpolation

    for i in range(N_x):
        if ((x >= xcoordinates[i-1]).any() and (x < xcoordinates[i]).any()):
            q = (xcoordinates[i] - x)/(xcoordinates[i] - xcoordinates[i-1])*q_x[i-1] + (x - xcoordinates[i-1])/(xcoordinates[i] - xcoordinates[i-1])*q_x[i]

    return q



x = np.arange(0,1.5,0.05)
print(q(x))
plt.plot(x,q(x))
plt.show()
