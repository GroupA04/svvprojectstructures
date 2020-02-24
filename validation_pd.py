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
jambendingskin_ave = np.sort(np.transpose(average(jambendingskin)), axis=0)
jamstraightskin_ave = np.sort(np.transpose(average(jamstraightskin)), axis=0)

bending_ave =np.sort(np.concatenate((bending_section1ave, bending_section2ave)), axis= 0) #combined sections into 1 array


#=======================================================Plotting==================================================================================
fig = plt.figure()
ax = fig.gca(projection='3d')

scatter = ax.scatter(elements_ave[:,0], elements_ave[:,1], elements_ave[:,2], c = jambendingskin_ave[:,1], cmap = 'hsv')

fig.colorbar(scatter)

ax.set_xlabel('Spanwise axis [mm]')
ax.set_ylabel('Height axis [mm]')
ax.set_zlabel('Chordwise axis [mm]')
plt.show()

