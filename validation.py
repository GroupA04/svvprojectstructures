import numpy as np
from numerical_functions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#get node locations and elements from data files

nodes = np.genfromtxt('nodes.txt', delimiter = ',')
elements = np.genfromtxt('elements', delimiter=',')


node_num = nodes[:,0]
x_loc = nodes[:,1]
y_loc = nodes[:,2]
z_loc = nodes[:,3]

element_nodes = np.delete(elements, 0, 1)

for i in range(len(node_num)+1):
    element_nodes = np.where(element_nodes == i, [i], element_nodes )
    element_nodes = np.where(element_nodes == i, [x_loc[i], y_loc[i], z_loc[i]], element_nodes)

print(element_nodes)


# fig = plt.figure()
# ax = fig.gca(projection='3d')
#
# ax.scatter(x_loc, y_loc, z_loc)
# ax.set_xlabel('Spanwise axis')
# ax.set_ylabel('height axis')
# ax.set_zlabel('chordwise axis')
# plt.show()

#import data files from B737.rpt
bending_section1 = np.array(np.genfromtxt('bending_section1.txt'))
bending_section2 = np.array(np.genfromtxt('bending_section2.txt'))

#average von Mises stress and shear stress

bending_section1ave = np.transpose(average(bending_section1))
bending_section2ave = np.transpose(average(bending_section2))

# index_lst = np.array([])
# for i in range(len(elements)):
#     index = np.where(bending_section1ave[:,0] == i)
#     print(index[0])
#     if index[0] > 0:
#         index_lst = np.vstack([i, index])
# print(index_lst)

