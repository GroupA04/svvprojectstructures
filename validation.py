import numpy as np
from numerical_functions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#get node locations from data file

nodes = np.genfromtxt('nodes.txt', delimiter = ',')


x_loc = nodes[:,1]
y_loc = nodes[:,2]
z_loc = nodes[:,3]



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
