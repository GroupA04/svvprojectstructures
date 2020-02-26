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

nodes = np.array(nodes)

#nodes with added deflection values (10x exagerated for graphical purpouse)
x1 = nodes[:,0] + 10*deflection_case1[:,2]
y1 = nodes[:,1] + 10*deflection_case1[:,3]
z1 = nodes[:,2] + 10*deflection_case1[:,4]

x2 = nodes[:,0] + 10*deflection_case2[:,2]
y2 = nodes[:,1] + 10*deflection_case2[:,3]
z2 = nodes[:,2] + 10*deflection_case2[:,4]

x3 = nodes[:,0] + 10*deflection_case3[:,2]
y3 = nodes[:,1] + 10*deflection_case3[:,3]
z3 = nodes[:,2] + 10*deflection_case3[:,4]

#======================================Find the deflection and shear stress values==============================================================
#---------Maximum shear stress
bending_shear_max = max(bending_ave[:,2])
jambending_shear_max = max(jambendingskin_ave[:,2])
jamstraightskin_shear_max = max(jamstraightskin_ave[:,2])

#---------Deflection at hinge line
#hinge line on y = 0, z = 0 along the span/x-axis
hinge_nodes_index = np.where((nodes[:,1] == 0) & (nodes[:,2] == 0))

index = hinge_nodes_index[0]
hinge_line = np.array([nodes[index, 0], nodes[index, 1], nodes[index, 2]])
for i in hinge_nodes_index[1:]:
    hinge_line= np.vstack((hinge_line, np.array([nodes[i,0], nodes[i,1], nodes[i,2]])))
hinge_line = hinge_line.transpose()

#find deflection at hinge line
hinge_def1 = np.array([])
hinge_def2 = np.array([])
hinge_def3 = np.array([])

for i in hinge_nodes_index:
    hinge_def1 = np.vstack([nodes[i,0] + deflection_case1[i,2], nodes[i,1] + deflection_case1[i,3], nodes[i,2] + deflection_case1[i,4]])
    hinge_def2 = np.vstack([nodes[i, 0] + deflection_case2[i, 2], nodes[i, 1] + deflection_case2[i, 3], nodes[i, 2] + deflection_case2[i, 4]])
    hinge_def3 = np.vstack([nodes[i, 0] + deflection_case3[i, 2], nodes[i, 1] + deflection_case3[i, 3], nodes[i, 2] + deflection_case3[i, 4]])

#Deflection over hinge line for all load cases
hinge_def1 = hinge_def1.transpose()
hinge_def2 = hinge_def2.transpose()
hinge_def3 = hinge_def3.transpose()

#Maximum deflection over hinge line for all load cases
hinge_def1_max = [hinge_def1[np.argmax(abs(hinge_def1[:,0])),0], hinge_def1[np.argmax(abs(hinge_def1[:,1])),1], hinge_def1[np.argmax(abs(hinge_def1[:,2])),2]]
hinge_def2_max = [hinge_def2[np.argmax(abs(hinge_def2[:,0])),0], hinge_def2[np.argmax(abs(hinge_def2[:,1])),1], hinge_def2[np.argmax(abs(hinge_def2[:,2])),2]]
hinge_def3_max = [hinge_def3[np.argmax(abs(hinge_def3[:,0])),0], hinge_def3[np.argmax(abs(hinge_def3[:,1])),1], hinge_def3[np.argmax(abs(hinge_def3[:,2])),2]]


#=========================================================Plotting==================================================================================
fig = plt.figure()
# ax = fig.gca(projection='3d')

#------------------------------------Select which data to plot-------------------------------------------------------

#vonMises = ax.scatter(elements_ave[:,0], elements_ave[:,1], elements_ave[:,2], c = jambendingskin_ave[:,1], cmap = 'coolwarm')
#shear = ax.scatter(elements_ave[:,0], elements_ave[:,1], elements_ave[:,2], c = jambendingskin_ave[:,2], cmap = 'coolwarm')
# deflection = ax.scatter(x2, y2, z2, c = deflection_case2[:,1], cmap = 'coolwarm')

#-------------------------------------------------------------------------------------------------------------------
#3D plotting of deflection
# hingeline = ax.scatter(hinge_line[:,0], hinge_line[:,1], hinge_line[:,2])
# hingedef1 = ax.scatter(hinge_def1[:,0], hinge_def1[:,1], hinge_def1[:,2])
# hingedef2 = ax.scatter(hinge_def2[:,0], hinge_def2[:,1], hinge_def2[:,2])
# hingedef3 = ax.scatter(hinge_def3[:,0], hinge_def3[:,1], hinge_def3[:,2])

#2D plotting of deflection
hingedef1_2d = plt.scatter(hinge_def1[:,0], hinge_def1[:,2])
hingedef2_2d = plt.scatter(hinge_def2[:,0], hinge_def2[:,2])
hingedef3_2d = plt.scatter(hinge_def3[:,0], hinge_def3[:,2])

#------------------------------------------------------------------------------------------------------------------

# fig.colorbar(deflection)
#
# ax.set_xlim3d(0,2500)
# ax.set_ylim3d(-1250,1250)
# ax.set_zlim3d(-1250,1250)
#
# ax.set_xlabel('Spanwise axis [mm]')
# ax.set_ylabel('Height axis [mm]')
# ax.set_zlabel('Chordwise axis [mm]')

plt.show()