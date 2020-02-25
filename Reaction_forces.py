#Reaction forces
from SVV import *

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

integral_1_la = fv_1
integral_2_la = fv_2
integral_3_la = fv_3
integral_4_la = fv_4
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
sc = -0.12 #change
theta_1 = theta*pi/180
x_act = x_2 - x_a/2


#Unknowns
a = ['R_1y', 'R_1z', 'R_2y', 'R_2z', 'R_3y', 'R_3z', 'R_act', 'C1','C2','C3','C4','C5']

#M_x(la), M_y(la), M_z(la), S_z(la), S_y(la), v(x1), v(x2), v(x3), w(x1), w(x2), w(x3)

M_x = [0,0,0,0,0,0,(l_a-x_2+x_a/2)*(cos(pi/6)-sin(pi/6)),0,0,0,0,0]

M_y = [0,-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),-cos(pi/6)*(l_a-x_2+x_a/2),0,0,0,0,0]

M_z = [-(l_a-x_1),0,-(l_a-x_2),0,-(l_a-x_3),0,-sin(pi/6)*(l_a-x_2+x_a/2),0,0,0,0,0]

S_z = [0,-1,0,-1,0,-1,-cos(pi/6),0,0,0,0,0]

S_y = [-1,0,-1,0,-1,0,-sin(pi/6),0,0,0,0,0]

v_x1 = [0,0,0,0,0,0,(cos(theta*pi/180)-sin(theta*pi/180))/(G*J)*(x_1**2/2-x_2*x_1+x_a*x_1/2)*sc,x_1,1,0,0,sc]

v_x2 = [(x_2-x_1)**3/(6*E*I_zz),0,0,0,0,0,sin(theta*pi/180)*(x_a/2)**3/(6*E*I_zz)+(cos(theta*pi/180)-sin(theta*pi/180))/(G*J)*(x_2**2/2-x_2*x_2+x_a*x_2/2)*sc,x_2,1,0,0,sc]

v_x3 = [(x_3-x_1)**3/(6*E*I_zz) , 0, (x_3-x_2)**3/(6*E*I_zz) , 0,0,0, sin(theta*pi/180)*(x_3-x_2+x_a/2)**3/(6*E*I_zz)+(cos(theta*pi/180)-sin(theta*pi/180))/(G*J)*(x_3**2/2-x_2*x_3+x_a*x_3/2)*sc ,x_3,1,0,0,sc]

w_x1 = [0,0,0,0,0,0,0,0,0,x_1,1,0]

w_x2 = [0,(x_2-x_1)**3/(6*E*I_yy),0,0,0,0,cos(theta*pi/180)*(x_a/2)**3/(6*E*I_yy),0,0,x_2,1,0]

w_x3 = [0,(x_3-x_1)**3/(6*E*I_yy),0,(x_3-x_2)**3/(6*E*I_yy),0,0,cos(theta*pi/180)*(x_3-x_2+x_a/2)**3/(6*E*I_yy),0,0,x_3,1,0]

w_x4 = [0,(x_act-x_1)**3/(6*E*I_yy),0,0,0,0,0,0,0,x_act,1,0]

#solve matrix

right_matrix = np.array([P*(l_a-x_2-x_a/2)*(cos(theta*pi/180)-sin(theta*pi/180))-integral_B1_la,
                         -P*cos(theta*pi/180)*(l_a-x_2-x_a/2),
                         -P*sin(theta*pi/180)*(l_a-x_2-x_a/2)-integral_2_la,
                         -P*cos(theta*pi/180),
                         -P*sin(theta*pi/180)-integral_1_la,
                         d_1*cos(theta*pi/180)+integral_4_x1/(E*I_zz)-P*(sin(theta*pi/180)-cos(theta*pi/180))/(G*J)*(x_1**2/2-x_2*x_1+x_a*x_1/2)*sc-1/(G*J)*integral_B3_x1*sc,
                         integral_4_x2/(E*I_zz)-P*(sin(theta*pi/180)-cos(theta*pi/180))/(G*J)*(x_2**2/2-x_2*x_2+x_a*x_2/2)*sc-1/(G*J)*integral_B3_x2*sc,
                         d_3*cos(theta*pi/180)+P*sin(theta*pi/180)*(x_3-x_2-x_a/2)**3/(6*E*I_zz)+integral_4_x3/(E*I_zz)-P*(sin(theta*pi/180)-cos(theta*pi/180))/(G*J)*(x_3**2/2-x_3*x_1+x_a*x_3/2)*sc-1/(G*J)*integral_B3_x3*sc,
                         d_1*sin(theta*pi/180),
                         0,
                         P*cos(theta*pi/180)*(x_3-x_2-x_a/2)**3/(6*E*I_yy)+d_3*sin(theta*pi/180),
                         0])

left_matrix = np.array([M_x, M_y,M_z,S_z,S_y,v_x1,v_x2,v_x3,w_x1,w_x2,w_x3,w_x4])


Reaction_forces = linalg.gmres(left_matrix,right_matrix,tol=1e-5)

#print(Reaction_forces)
#
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

#print(R_1y)
print(R_1z)
# print(R_2y)
print(R_2z)
# print(R_3y)
print(R_3z)




# x = 0
# v_x = -1/(E*I_zz)*(R_1y*macaulay(x,x_1)**3 /6-R_2y*macaulay(x,x_2)**3/6-R_3y*macaulay(x,x_3)**3/6-R_act*sin(theta_1)*macaulay(x,x_2-x_a/2)**3/6+P*sin(theta_1)/6*macaulay(x,x_2+x_a/2)**3)+C1*x+C2
# print(v_x)
#plot deflections
x = 0
xlist = []
v_xlist = []
w_xlist = []
angle_of_twist = []


#function deflection
for i in range(99):
    v_x = -1/(E*I_zz)*(R_1y*macaulay(x,x_1)**3 /6+int4-R_2y*macaulay(x,x_2)**3/6-R_3y*macaulay(x,x_3)**3/6-R_act*sin(theta_1)*macaulay(x,x_2-x_a/2)**3/6+P*sin(theta_1)/6*macaulay(x,x_2+x_a/2)**3)+C1*x+C2
    w_x = -1/(E*I_zz)*(R_1z*macaulay(x,x_1)**3/6-R_2z*macaulay(x,x_2)**3/6-R_3z*macaulay(x,x_3)**3/6)+C3*x+C4
    deflection_angle_of_twist = ((cos(theta*pi/180)-sin(theta*pi/180))*R_act+P*(sin(theta*pi/180)-cos(theta*pi/180)))*sc/(G*J)*(x**2/2-x_2*x+x_a/2*x)+sc/(G*J)*intb3+ C5*sc
    v_xlist.append(v_x)
    w_xlist.append(w_x)
    angle_of_twist.append(deflection_angle_of_twist)
    xlist.append(x)
    x = x + l_a/100

# print(len(xlist))
# print(len(angle_of_twist))
# print(len(v_xlist))

deflection_vlist = []
for j in range(99):
    deflection_v= v_xlist[j]+angle_of_twist[j]
    deflection_vlist.append(deflection_v)



# print(deflection_vlist[0])

#plt.plot(xlist,deflection_vlist)

# plt.plot(xlist,v_xlist)
# plt.plot(xlist,w_xlist)
# plt.show()

#
# #
#
#
#



