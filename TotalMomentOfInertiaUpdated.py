#All variables from Fokker 100 Data sheet
#Everything defined as meter, Newton, or deg

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
P = 49200  # [Newton]

'''
C_a = 0.605  # m
l_a = 2.661  # m
x_1 = 0.172  # m
x_2 = 1.211  # m
x_3 = 2.591  # m
x_a = 0.35   # m
h = 0.205  # m
t_sk = 1.1/1000  # m
t_sp = 2.8/1000  # m
t_st = 1.2/1000  # m
h_st = 16./1000   # m
w_st = 19./1000   # m
n_st = 15  # -
d_1 = 0.01154  # m
d_3 = 0.01840  # m
#theta = radians(28)  # rad
P = 97.4*1000  # N
'''
#Moment of Inertia calculation

import numpy as np

#Centroid calculation - centroid lies on the z axis because it is symmetrical around z 
#The position due on the z can be calculated by multypling sum of areas with distances and divide by sum of areas 
#This is done in unrotated coordinate system

#perimeter of the cross section and stiffener spacing
diagonallength = np.sqrt((h/2)**2 + (C_a - h/2)**2)
lengthskin = np.pi * h/2 + 2 * diagonallength
ST_spacing = lengthskin / n_st

#stringer area
A_st = t_st *( h_st + w_st)

#intermediate parameters for stringer coordinates (a = stringer_1 to spar and b= spar to stringer_2 ) - stringer notation starts from stringer_0
# phi is the angle of the upper corner of the triangular part of the cross-section
phi= np.arctan((C_a - h/2)/(h/2))

a = np.pi* h / 4 - ST_spacing 
b= ST_spacing - a

#stringer coordinates -  array of 11 stringers by 2 coordinates - starting from LE going around the upper skin to TE to lower skin
#coordinate system is at the TE with y positive upwards and z positive towards the LE

stringercoordinates = np.zeros([n_st ,2])
stringercoordinates[0,0] = C_a 
stringercoordinates[0,1] = 0.0
stringercoordinates[1,0] = C_a - (h/2 - h/2 * np.cos(2*ST_spacing/ h))
stringercoordinates[1,1] = h/2 * np.sin (2* ST_spacing / h)
stringercoordinates[2,0] = C_a - (h/2 + b* np.sin(phi))
stringercoordinates[2,1] = h/2 -b* np.cos(phi)
stringercoordinates[3,0] = C_a - (h/2 + (b+ ST_spacing)* np.sin(phi))
stringercoordinates[3,1] = h/2 - (b+ ST_spacing) * np.cos(phi)
stringercoordinates[4,0] = C_a - (h/2 + (b+ 2* ST_spacing) * np.sin(phi))
stringercoordinates[4,1] = h/2 - (b+ 2* ST_spacing) * np.cos(phi)
stringercoordinates[5,0] = C_a - (h/2 + (b+ 3*ST_spacing) * np.sin(phi))
stringercoordinates[5,1] = h/2 - (b+ 3 * ST_spacing) * np.cos(phi)
stringercoordinates[6,0] = C_a - (h/2 + (b+ 3*ST_spacing) * np.sin(phi))
stringercoordinates[6,1] = -(h/2 - (b+ 3 * ST_spacing) * np.cos(phi))
stringercoordinates[7,0] = C_a - (h/2 + (b+ 2* ST_spacing) * np.sin(phi))
stringercoordinates[7,1] = -(h/2-(b+ 2* ST_spacing) * np.cos(phi))
stringercoordinates[8,0] = C_a - (h/2 + (b+ ST_spacing)* np.sin(phi))
stringercoordinates[8,1] = -(h/2 - (b+ ST_spacing) * np.cos(phi)) 
stringercoordinates[9,0] = C_a - (h/2 + b* np.sin(phi))
stringercoordinates[9,1] = -(h/2 -b* np.cos(phi))
stringercoordinates[10,0] = C_a - (h/2 - h/2 * np.cos(2*ST_spacing/ h))
stringercoordinates[10,1] = -(h/2 * np.sin (2* ST_spacing / h))

#centroid calculation
y_bar= 0.0
z_bar_top = sum(stringercoordinates[:,0]*A_st) + t_sp * h * (C_a - h/2) + (C_a - h/2) / 2 * t_sk * diagonallength * 2 + (h/np.pi + C_a - h/2) * np.pi * h/2 * t_sk
z_bar_bot = A_st * n_st + t_sp * h + 2 * t_sk * diagonallength + np.pi *h/2 * t_sk
z_bar = z_bar_top / z_bar_bot
centroid = [z_bar,y_bar]

I_zz_st = A_st * sum((stringercoordinates[:,1])**2)      # MoI around z-axis due to the stringers
I_yy_st = A_st * sum(((stringercoordinates[:,0])-z_bar)**2)      # MoI around y-axis due to the stringers


I_zz_sk_diagonalpart = (t_sk * (2*diagonallength)**3 * (h/2/diagonallength)**2) / 12        # MoI around z-axis due to the skin of the diagonal part of cross-section
I_yy_sk_diagonalpart = 2* ((t_sk * (diagonallength)**3 * ((C_a-h/2)/diagonallength)**2) / 12 + t_sk * diagonallength * (z_bar - (C_a-h/2) / 2 )**2) # MoI around y-axis due to the skin of the diagonal part of cross-section


I_zz_sk_arc = np.pi * (h/2)**3 * t_sk / 2                                           # MoI around z-axis due to the skin of the arc in cross-section 
I_yy_sk_arc = I_zz_sk_arc + np.pi * h/2 * t_sk * ((C_a - h/2 - z_bar)**2 - (2*h/2/np.pi)**2)   #(4*h/2/3/np.pi)**2)     # MoI around y-axis due to the skin of the arc in cross-section

I_zz_sk = I_zz_sk_diagonalpart + I_zz_sk_arc        # Moment of Inertia around z-axis due to the skin
I_yy_sk = I_yy_sk_diagonalpart + I_yy_sk_arc        # Moment of Inertia around y-axis due to the skin

I_zz_sp = h**3 * t_sp / 12 
I_yy_sp = h * t_sp * (C_a - h/2 - z_bar)**2

I_zz = I_zz_st + I_zz_sk + I_zz_sp   # Total Moment of Inertia around z-axis
I_yy = I_yy_st + I_yy_sk + I_yy_sp   # Total Moment of Inertia around y-axis