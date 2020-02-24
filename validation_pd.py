import numpy as np
import pandas as pd
from numerical_functions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#=========================================Import node and element data==========================================================================
#get node locations and elements from data files

nodes = pd.read_csv('nodes.txt', delimiter=',',
                    names = ['node_id', 'x', 'y', 'z'],
                    index_col = 0)
elements = pd.read_csv('Elements', delimiter=',',
                       names = ['number', 'node1', 'node2', 'node3', 'node4'],
                       index_col = 0)

# Check dataframes
# print(nodes.head())
# print(elements.head())

# Couple node and element data
node_defs = ['node1','node2','node3','node4']
a = [np.array(nodes.loc[elements.loc[nr, node_defs]]) for nr in elements.index]

# Verify if shape is (6634,4,3)
if np.array(a).shape == (6634,4,3):
    print('Shape is correct')
else:
    print('Shape is incorrect')
    print('Shape is', np.array(a).shape, 'Shape should be (6634,4,3)')


# Access element number 5
# print(a[4])

#take average posistion of nodes to get elements as point location
elements_ave = []
for i in range(len(a)):
    elements_ave_i = np.mean(a[i], axis = 0)
    elements_ave.append(elements_ave_i)
elements_ave = np.array(elements_ave)

#verify if averages are correct
# print(a[1])
# print(elements_ave[1,:])

#===============================================Import FEM data files from B737.rpt=================================================================

#==============================================Von Mises stress and shear stress calculations=======================================================
bending_section1 = np.array(np.genfromtxt('bending_section1.txt'))
bending_section2 = np.array(np.genfromtxt('bending_section2.txt'))

jambendingskin = np.array(np.genfromtxt('jambendingskin'))
jamstraightskin = np.array(np.genfromtxt('jamstraightskin'))

#average von Mises stress and shear stress

bending_section1ave = np.transpose(average(bending_section1))
bending_section2ave = np.transpose(average(bending_section2))

#sorted average stresses by element number
jambendingskin_ave = np.transpose(average(jambendingskin))[np.transpose(average(jambendingskin))[:,0].argsort()]
jamstraightskin_ave = np.transpose(average(jamstraightskin))[np.transpose(average(jamstraightskin))[:,0].argsort()]

bending_ave =np.concatenate((bending_section1ave, bending_section2ave))[np.concatenate((bending_section1ave, bending_section2ave))[:,0].argsort()] #combined sections into 1 array

#==================================================Deflection data===================================================================================
deflection_case1 = np.array(np.genfromtxt('Deflection1')) #bending case
deflection_case2 = np.array(np.genfromtxt('Deflection2')) #jammed bending case
deflection_case3 = np.array(np.genfromtxt('Deflection3')) #jammed straight case



#=======================================================Plotting==================================================================================
fig = plt.figure()
ax = fig.gca(projection='3d')

vonMises = ax.scatter(elements_ave[:,0], elements_ave[:,1], elements_ave[:,2], c = jambendingskin_ave[:,1], cmap = 'coolwarm')
shear = ax.scatter(elements_ave[:,0], elements_ave[:,1], elements_ave[:,2], c = jambendingskin_ave[:,2], cmap = 'coolwarm')
#deflection = ax.scatter(elements_ave[:,0] + deflection_case1[:,2], elements_ave[:,1], elements_ave[:,2], c = deflection_case1[:,1], cmap = 'coolwarm')

fig.colorbar(shear)

ax.set_xlim3d(0,2500)
ax.set_ylim3d(-1250,1250)
ax.set_zlim3d(-1250,1250)

ax.set_xlabel('Spanwise axis [mm]')
ax.set_ylabel('Height axis [mm]')
ax.set_zlabel('Chordwise axis [mm]')
plt.show()

