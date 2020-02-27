#Reaction forces
from SVV import int1,int2,int3,int4,intb3
from math import *
import numpy as np
from scipy.sparse import linalg
import matplotlib.pyplot as plt

def macaulay(x, step):
    if (x - step) >= 0:
        return (x - step)
    else:
        return 0
#paramters
C_a = 0.505 # [m] 
l_a = 1.611 # [m]
x_1 = 0.125 # [m]
x_2 = 0.498 # [m]
x_3 = 1.494 # [m]
x_a = 0.245 # [m]
h = 0.161   # [m] 
t_sk = 0.0011 # [m] 
t_sp = 0.0024 # [m] 
t_st = 0.0012 # [m] 
h_st = 0.012  # [m] 
w_st = 0.017  # [m] 
n_st = 11
d_1 = 0.00389 # [m] 
d_3 = 0.01245 # [m] 
theta = 30 # [degrees]
P = 49.2  # [KN]
I_zz = 4.721698865*10**-6
I_yy = 4.175337875*10**-5 
E = 73.1*(10**6)

integral_1_la = 6.19
integral_2_la = 5.01 
integral_3_la = 2.538
integral_4_la = 0.96
integral_4_x1 = 1.4189e-05 
integral_4_x2 = 0.0059978 
integral_4_x3 = 0.7018
integral_B1_la = 24.0618
integral_B3_x1 = 0.00156
integral_B3_x2 = 0.19 
integral_B3_x3 = 7.57 
G = 28*(10**6)
T = 1
dtheta_dz = 1 #change
J = 0.000175317 #T/(G*dtheta_dz)
sc = 0.12
z_sc = -0.12 #change
theta_1 = theta*pi/180
x_act = x_2 - x_a/2


#Unknowns
a = ['R_1y', 'R_1z', 'R_2y', 'R_2z', 'R_3y', 'R_3z', 'R_act', 'C1','C2','C3','C4','C5']

#M_x(la), M_y(la), M_z(la), S_z(la), S_y(la), v(x1), v(x2), v(x3), w(x1), w(x2), w(x3)

M_x = [-(sc-(h/2)),0,-(sc-(h/2)),0,-(sc-(h/2)),0,((sc)*sin(pi/6)-(h/2)*cos(pi/6)),0,0,0,0,0]

M_y = [0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(pi/6)*(l_a-(x_2-x_a/2)),0,0,0,0,0]

M_z = [-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),0,-sin(pi/6)*(l_a-(x_2-x_a/2)),0,0,0,0,0]

S_z = [0,-1,0,-1,0,-1,-cos(pi/6),0,0,0,0,0]

S_y = [-1,0,-1,0,-1,0,-sin(pi/6),0,0,0,0,0]

v_x1 = [0,0,0,0,0,0,0,x_1,1,0,0,(sc-(h/2))]

v_x2 = [((x_2-x_1)**3)/(6*E*I_zz)-(1/(G*J))*((sc-(h/2))**2)*(x_2-x_1),0,0,0,0,0,(sin(theta*pi/180))*((x_a/2)**3)*(1/(6*E*I_zz))+(sc*sin(theta*pi/180))*(1/(G*J))*(x_a/2)*(sc-(h/2))-(h*cos(pi/6))*(1/(2*G*J))*(x_a/2)*(sc-(h/2)),x_2,1,0,0,(sc-(h/2))]

v_x3 = [((x_3-x_1)**3)/(6*E*I_zz)-(1/(G*J))*((sc-(h/2))**2)*(x_3-x_1),0,((x_3-x_2)**3)/(6*E*I_zz)-(1/(G*J))*((sc-(h/2))**2)*(x_3-x_2),0,0,0,((sin(theta*pi/180))*((x_3-x_2+x_a/2)**3)*(1/(6*E*I_zz))-(h*cos(theta*pi/180))/(2*G*J)*((x_3-x_2+x_a/2)**3)-(sc*sin(theta*pi/180))/(G*J)*(x_3-x_2+x_a/2)*(sc-(h/2))),x_3,1,0,0,(sc-(h/2))]

w_x1 = [0,0,0,0,0,0,0,0,0,x_1,1,0]

w_x2 = [0,((x_2-x_1)**3)/(6*E*I_yy),0,0,0,0,(cos(theta*pi/180))*((x_a/2)**3)/(6*E*I_yy),0,0,x_2,1,0]

w_x3 = [0,(x_3-x_1)**3/(6*E*I_yy),0,((x_3-x_2)**3)/(6*E*I_yy),0,0,(cos(theta*pi/180))*((x_3-x_2+x_a/2)**3)/(6*E*I_yy),0,0,x_3,1,0]

w_x4 = [0,((x_act-x_1)**3)/(6*E*I_yy),0,0,0,0,0,0,0,x_act,1,0]

#solve matrix

right_matrix = np.array([-P*((h/2)*cos(theta*pi/180)-(sc)*sin(theta*pi/180))                  -integral_B1_la,
                         -P*cos(theta*pi/180)*(l_a-(x_2+x_a/2)),
                         -P*sin(theta*pi/180)*(l_a-(x_2+x_a/2))                               -integral_2_la,
                         -P*cos(theta*pi/180),
                         -P*sin(theta*pi/180)                                                 -integral_1_la,
                         d_1*cos(theta*pi/180)                                                +(integral_4_x1)*(1/(E*I_zz))                                -(1/(G*J))*(sc-(h/2))*(integral_B3_x1),
                         (integral_4_x2)/(E*I_zz)                                             -(1/(G*J))*(sc-(h/2))*(integral_B3_x2),
                         d_3*cos(theta*pi/180)                                                -P*((h*cos(pi/6))/(2*G*J)*(x_3-x_2-x_a/2)*(sc-(h/2))-sin(theta*pi/180)*((x_3-x_2-x_a/2)**3)/(6*E*I_zz)-(sc*sin(pi/6))*(1/(G*J))*(x_3-x_2-x_a/2)*(sc-(h/2)))       +(integral_4_x3)/(E*I_zz)                     -(sc-(h/2))*(1/(G*J))*(integral_B3_x3),
                         d_1*sin(theta*pi/180),
                         0,
                         P*(cos(theta*pi/180))*((x_3-x_2-x_a/2)**3)/(6*E*I_yy)                    +d_3*sin(theta*pi/180),
                         0])

left_matrix = np.array([M_x, M_y,M_z,S_z,S_y,v_x1,v_x2,v_x3,w_x1,w_x2,w_x3,w_x4])


Reaction_forces = linalg.gmres(left_matrix,right_matrix,tol=1e-5)

#print(Reaction_forces)

# #Solutions Matrix
R_1y = Reaction_forces[0][0]
R_1z = Reaction_forces[0][1]
R_2y = Reaction_forces[0][2]
R_2z = Reaction_forces[0][3]
R_3y = Reaction_forces[0][4]
R_3z = Reaction_forces[0][5]
R_act = Reaction_forces[0][6]
C1 = Reaction_forces[0][7]
C2 = Reaction_forces[0][8]
C3 = Reaction_forces[0][9]
C4 = Reaction_forces[0][10]
C5 = Reaction_forces[0][11]

print('R_1y',R_1y)
print('R_1z', R_1z)
print('R_2y',R_2y)
print('R_2z',R_2z)
print('R_3y',R_3y)
print('R_3z',R_3z)
print('C1',C1)
print('C2',C2)
print('C3',C3)
print('C4',C4)
print('C5',C5)




#plot deflections

xlist = []
v_xlist = []
w_xlist = []
angle_of_twist = []



x=0

#function deflection
for i in range(41):
    v_x = -1/(E*I_zz)*(R_1y*macaulay(x,x_1)**3 /6-int4[i]+R_2y*macaulay(x,x_2)**3/6+R_3y*macaulay(x,x_3)**3/6-R_act*sin(theta_1)*macaulay(x,x_2-x_a/2)**3/6-P*sin(theta_1)/6*macaulay(x,x_2+x_a/2)**3)*+C1*x+C2
    w_x = (-1/(E*I_zz)*(-R_1z*macaulay(x,x_1)**3/6-R_2z*macaulay(x,x_2)**3/6-R_3z*macaulay(x,x_3)**3/6+R_act*cos(theta_1)*macaulay(x,x_2-x_a/2)**3+1/6*P*cos(theta_1)*macaulay(x,x_2+x_a/2)**3)*-1+C3*x+C4)/1000
    deflection_angle_of_twist = (1/(G*J))*((cos(theta*pi/180)-sin(theta*pi/180))*R_act*macaulay(x,x_2-x_a/2)**2+P*(sin(theta*pi/180)-cos(theta*pi/180))*macaulay(x,x_2+x_a/2)**2 + intb3[i])+C5
    v_xlist.append(v_x)
    w_xlist.append(w_x)
    angle_of_twist.append(deflection_angle_of_twist)
    xlist.append(x)
    x = x + l_a/41




deflection_wlist = []
deflection_vlist = []
for j in range(41):
    deflection_v= v_xlist[j]+(angle_of_twist[j])*(sc-(-h/2))
    deflection_vlist.append(deflection_v)
    deflection_w = w_xlist[j]
    deflection_wlist.append(deflection_w)



#Plotting the Graphs

plt.plot(xlist,deflection_vlist)
plt.title("Deflection in y-direction")
plt.xlabel("x [m]")
plt.ylabel("Deflection [m]")
plt.show()

plt.plot(xlist,deflection_wlist)
plt.title("Deflection in z-direction")
plt.xlabel("x [m]")
plt.ylabel("Deflection [m]")
plt.show()

plt.plot(xlist,angle_of_twist)
plt.title("Angle of twist")
plt.xlabel("x [m]")
plt.ylabel("Twist [-]")
plt.show()



## Attempt to Get the Correct V (attempt failed)
##
##
##v_x_tlist = []
##x_tlist = []
##x=0
##
##for i in range(41):
##    v_x_t = R_1y*   ((1/(6*E*I_zz))*(macaulay(x,x_1)**3) -  (1/(G*J))*((sc-(h/2))**2)*(macaulay(x,x_1))) + \
##            R_act*  (sin(pi/6)*(macaulay(x,(x_2-x_a/2))**3)/(6*E*I_zz) + (sc*sin(pi/6))*(macaulay(x,x_2-x_a/2))*(sc-(h/2))/(G*J) - (h*cos(pi/6)*macaulay(x,x_2-x_a/2)*(sc-(h/2))/(2*G*J)))  + \
##            R_2y*   ((1/(6*E*I_zz))*(macaulay(x,x_2))**3  - (1/(G*J))*((sc-(h/2))**2)*(macaulay(x,x_2)) )  + \
##            R_3y*   ((1/(6*E*I_zz))*(macaulay(x,x_3)**3)  - (1/(G*J))*((sc-(h/2))**2)*(macaulay(x,x_3)) )  + \
##            P*      ((h*cos(pi/6))*(macaulay(x,(x_2+x_a/2)))*(sc-(h/2))/(2*G*J)  -  (sin(pi/6))*(macaulay(x,x_2+x_a/2)**3)/(6*E*I_zz)  - (sc*sin(pi/6))*(macaulay(x,x_2+x_a/2))*(sc-(h/2))/(G*J) ) + \
##            intb3[i]*(1/(G*J))*(sc-(h/2)) - \
##            int4[i]*(1/(E*I_zz)) + \
##            C1*x + C2 + C5*(sc-(h/2))
##    v_x_tlist.append(v_x_t)
##    x_tlist.append(x)
##    x = x + l_a/41
##
##plt.plot(x_tlist,v_x_tlist)
##plt.title("Deflection in y direction")
##plt.xlabel("x [m]")
##plt.ylabel("Deflection [m]")
##plt.show()    
