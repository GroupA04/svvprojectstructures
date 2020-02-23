import numpy as np
import pandas as pd
from numerical_functions import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#get node locations and elements from data files



nodes = pd.read_csv('nodes.txt', delimiter=',',
                    names = ['node_id', 'x', 'y', 'z'],
                    index_col = 0)
elements = pd.read_csv('Elements', delimiter=',',
                       names = ['number', 'node1', 'node2', 'node3', 'node4'],
                       index_col = 0)

# Check dataframes
print(nodes.head())
print(elements.head())

# Couple
node_defs = ['node1','node2','node3','node4']
a = [np.array(nodes.loc[elements.loc[nr, node_defs]]) for nr in elements.index]

# Check shape
print(np.array(a).shape)

# Access element number 5
print(a[4])


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

