import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#get node locations from data file

nodes = np.genfromtxt('nodes.txt', delimiter = ',')
print(len(nodes))

x_loc = nodes[:,1]
y_loc = nodes[:,2]
z_loc = nodes[:,3]

# X,Y = np.meshgrid(x_loc,y_loc)


fig = plt.figure()
ax = fig.gca(projection='3d')

ax.scatter(x_loc, y_loc, z_loc)
ax.set_xlabel('Spanwise axis')
ax.set_ylabel('height axis')
ax.set_zlabel('chordwise axis')
plt.show()