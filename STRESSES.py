# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 21:36:50 2020

@author: Mihnea
"""

#STRESSES MAIN CODE

#IMPORTED FUNCTIONS
import matplotlib.pyplot as plt
import numpy as np
from aerodynamicloading import *
from numeric_functions import *
from Reaction_forces import *

#GEOMETRY, MATEIRLA AND GIVEN LOADS - standard units
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
A_I = np.pi/2 * (h/2)**2
A_II = h/2 * (C_a-h/2)
G = 27.1 * 10**9

#Moment of Inertia calculation

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
z_bar= ((np.pi * h/4 * t_sk) * (C_a - h/4) + h/2 * t_sp * (C_a -h/2) + diagonallength * t_sk* (C_a -h/2)/2 + A_st * sum(stringercoordinates[1:6,0])+ A_st/2 * stringercoordinates[0,0])  / (5.5* A_st + h/2 * t_sp + lengthskin/2 * t_sk)                        
centroid = [z_bar,y_bar]

I_zz_st = A_st * sum((stringercoordinates[:,1])**2)      # MoI around z-axis due to the stringers
I_yy_st = A_st * sum(((stringercoordinates[:,0])-z_bar)**2)      # MoI around y-axis due to the stringers


I_zz_sk_diagonalpart = (t_st * (2*diagonallength)**3 * (h/2/diagonallength)**2) / 12        # MoI around z-axis due to the skin of the diagonal part of cross-section
I_yy_sk_diagonalpart = 2* ((t_st * (diagonallength)**3 * ((C_a-h/2)/diagonallength)**2) / 12 + t_sk * diagonallength * (z_bar - (C_a-h/2) / 2 )**2)#+ 2 *diagonallength * t_sk * z_bar**2 # MoI around y-axis due to the skin of the diagonal part of cross-section


I_zz_sk_arc = np.pi * (h/2)**3 * t_sk / 2                                           # MoI around z-axis due to the skin of the arc in cross-section 
I_yy_sk_arc = I_zz_sk_arc + np.pi * h/2 * t_st * ((C_a - h/2 - z_bar)**2 - (4*h/2/3/np.pi)**2)     # MoI around y-axis due to the skin of the arc in cross-section

I_zz_sk = I_zz_sk_diagonalpart + I_zz_sk_arc        # Moment of Inertia around z-axis due to the skin
I_yy_sk = I_yy_sk_diagonalpart + I_yy_sk_arc        # Moment of Inertia around y-axis due to the skin

I_zz_sp = h**3 * t_sp / 12 
I_yy_sp = h * t_sp * (C_a - h/2 - z_bar)**2

I_zz = I_zz_st + I_zz_sk + I_zz_sp   # Total Moment of Inertia around z-axis
I_yy = I_yy_st + I_yy_sk + I_yy_sp   # Total Moment of Inertia around y-axis


#SHEAR FORCES ACROSS THE SPAN
SY=[]
SZ=[]
M_z=[]
M_y=[]
sum_moment_z=0
sum_moment_y=0
    #AERODYNAMIC LOAD
Q=[]
for w in range(len(q_list)-1):
    
    Q.append((q_list[w]+q_list[w+1])/2* 10**3 * l_a/(len(q_list)))

x=np.linspace(0,l_a,100)
for i in range(len(x)):
    
    if i*l_a/len(x) < x_1:
        SZ.append(0)
        SY.append(Q[i])
        
    elif x_1<i*l_a/len(x)<(x_2-x_a/2):
        SZ.append(R_1z)
        SY.append(Q[i]+R_1y)
        
    elif (x_2-x_a/2)<i*l_a/len(x)< x_2:
        SZ.append(R_1z-R_act*cos(np.pi/6))
        SY.append(Q[i]+R_1y+R_act*sin(np.pi/6))
        
    elif x_2<i*l_a/len(x)<(x_2+x_a/2):
        SZ.append(R_1z-R_act*cos(np.pi/6)+R_2z)
        SY.append(Q[i]+R_1y+R_act*sin(np.pi/6)+R_2y)
        
    elif (x_2+x_a/2) < i*l_a/len(x)<x_3:
        SZ.append(R_1z-R_act*np.cos(np.pi/6)+R_2z+ P*np.cos(np.pi/6))
        SY.append(Q[i]+R_1y+R_act*sin(np.pi/6)+R_2y+P*np.sin(np.pi/6))
        
    else:
        SZ.append(R_1z-R_act*np.cos(np.pi/6)+R_2z+ P*np.cos(np.pi/6)+R_3z)
        SY.append(Q[i]+R_1y+R_act*sin(np.pi/6)+R_2y+P*np.sin(np.pi/6)+ R_3y)
    
for i in range(len(x)-1):
    M_z.append(sum_moment_z + (SZ[i]+SZ[i+1])/2 *l_a/len(x))
    M_y.append(sum_moment_y + (SY[i]+SY[i+1])/2 *l_a/len(x))
    sum_moment_z=np.sum(M_z)
    sum_moment_y=np.sum(M_y)
    

#FUNCTIONS INSIDE INTEGRALS AND INTEGRATORS
        
#functions inside integrals under the form coordinates* thickness

def f1(s): #upper diagonal flange y value 
    return t_sk*s* np.cos(phi)
def f2(s): # upper diagonal flange z value
    return t_sk *s * np.sin(phi)

def f3(s): # semicircle y value
    return t_sk * h/2 * np.cos(2*s/h)
def f4(s): #semicircle z value
    return ((C_a- h/2) + h/2 * np.sin(2*s/h)) * t_sk

def f5(s): #lower diagonal flange y value
    return t_sk * (-h/2 + np.cos(phi) * s)
def f6(s): #lower diagonal flange z value
    return t_sk * ((C_a -h/2) - s* np.sin(phi))

def f7(s): #spar y value up down
    return (h/2 - s )*t_sp
def f8(s): #spar z value 
    return (C_a - h/2) *t_sp    
   
def f9(s): #spar y value down up
    return (-h/2 + s)*t_sp
def f10(s): #spar z value 
    return (C_a - h/2)*t_sp 

def skinintegral(a,b,f):
    value = 0 
    result = 0 
    n1=100
    for i in range(1, n1+1):
        value += f(a+((i-(1/2))*((b-a)/n1)))
    result = ((b-a)/n1)*value 
    return result

#this is adjusted to give 5 outputs along a piece of the flane with are the result of the integration from the starting point up until that point
def skinintegral_F(n2,f,a,b):
    F=0
    baseshear_array = []
    for i in range(n2):
        F = skinintegral(i/n2*(b-a)+ a, (i+1)/n2*(b-a)+ a, f) + skinintegral(0,i/n2*(b-a)+a,f)
        baseshear_array.append(F)
    return baseshear_array
#CONSTRUCTING POSITION ARRAYS
zarray = np.linspace(0, C_a-h/2, 10)
yarray = np.linspace(0, h, 10)
y2array=np.linspace(0,h/2,10)
zinversearray= np.linspace(C_a-h/2, 0 ,10)
zprimearray=np.linspace(C_a-h/2, C_a,5)
azarray=np.linspace(C_a,C_a-h/2,5)
q_semicircle_along_span= []
q_upperdiag_along_span= []
q_lowerdiag_along_span= []
q_spar_I_along_span= []
q_spar_II_along_span= []
q_spar_along_span=[]
#FINDING THE BASE SHEAR FLOW AT EACH SECTION ALONG THE SPAN 
for i in range(1):
    i=10
    SZ[10]=150000
    SY[10]=-50000
    #SZ[10]= 200000
    #SY[10]=-250000
    baseshear_upperdiag=[]
    baseshear_lowerdiag=[]
    baseshear_spar_II=[]
    baseshear_spar_I=[]
    baseshear_semicircle=[]
    previous_b=0
    previous_d=0
    previous_f=0
    previous_h=0
    previous_j=0
    previous_l=0
#CELL II
    for j in range(len(zarray)):
        if j < stringercoordinates[5,0]/C_a*len(zarray):
            B_z= 0* np.ones(5)
            B_y= 0* np.ones(5)
        if stringercoordinates[5,0]/C_a * len(zarray) < j < stringercoordinates[4,0]/C_a * len(zarray):
            B_z= A_st * sum(stringercoordinates[5:6,0]) * np.ones(5)
            B_y= A_st * sum(stringercoordinates[5:6,1]) *np.ones(5)
        if stringercoordinates[4,0]/C_a * len(zarray) < j < stringercoordinates[3,0]/C_a * len(zarray):
            B_z= A_st * sum(stringercoordinates[4:6,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[4:6,1])*np.ones(5)
        if stringercoordinates[3,0]/C_a * len(zarray) < j < stringercoordinates[2,0]/C_a * len(zarray):
            B_z= A_st * sum(stringercoordinates[3:6,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[3:6,1])*np.ones(5)
        if j> stringercoordinates[2,0]/C_a * len(zarray):
            B_z= A_st * sum(stringercoordinates[2:6,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[2:6,1])*np.ones(5)
        baseshear_upperdiag.extend(np.add([-SY[i]/I_zz *x for x in  np.add(skinintegral_F(5,f1, previous_b, (j+1)/len(zarray)*diagonallength) , B_y)] , [- SZ[i]/I_yy *x for x in np.add(skinintegral_F(5,f2, previous_b, (j+1)/len(zarray)*diagonallength), B_z)]))
        previous_b = (j+1)/len(zarray)*diagonallength
    q_b_total_upperdiag_array = [baseshear_upperdiag[-1] *x for x in np.ones(5)]
    for k in range (len(yarray)):
        baseshear_spar_II.extend(np.add(np.add(q_b_total_upperdiag_array, [- SY[i]/I_zz * x for x in  (skinintegral_F(5,f7, previous_d,  (k+1)/len(yarray)*h))]) , [- SZ[i]/I_yy * x for x in (skinintegral_F(5,f8,previous_d, (k+1)/len(yarray)*h))]))
        previous_d = (k+1) *h/len(yarray)
    q_b_total_spar_array = [baseshear_spar_II[-1] *x for x in np.ones(5)]
    for l in range(len(zinversearray)):
        if  l > stringercoordinates[9,0]/C_a*len(zinversearray) :
            B_z= 0* np.ones(5)
            B_y= 0* np.ones(5)
        if stringercoordinates[9,0]/C_a * len(zinversearray) > l > stringercoordinates[8,0]/C_a * len(zinversearray):
            B_z= A_st * sum(stringercoordinates[9:10,0]) * np.ones(5)
            B_y= A_st * sum(stringercoordinates[9:10:6,1]) *np.ones(5)
        if stringercoordinates[8,0]/C_a * len(zinversearray) > l > stringercoordinates[7,0]/C_a * len(zinversearray):
            B_z= A_st * sum(stringercoordinates[8:10,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[8:10,1])*np.ones(5)
        if stringercoordinates[7,0]/C_a * len(zinversearray) > l > stringercoordinates[6,0]/C_a * len(zinversearray):
            B_z= A_st * sum(stringercoordinates[7:10,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[7:10,1])*np.ones(5)
        if l< stringercoordinates[6,0]/C_a * len(zinversearray):
            B_z= A_st * sum(stringercoordinates[6:10,0]) *np.ones(5)
            B_y= A_st * sum(stringercoordinates[6:10,1])*np.ones(5)
        baseshear_lowerdiag.extend(np.add(np.add(q_b_total_spar_array, [-SY[i]/I_zz * x for x in  np.add(skinintegral_F(5,f5, previous_f,  (l+1)/len(zinversearray)*diagonallength), B_y)]) , [- SZ[i]/I_yy * x for x in np.add(skinintegral_F(5,f6,previous_f, (l+1)/len(zinversearray)*diagonallength),B_z)]))
        previous_f=(l+1)/len(zinversearray)*diagonallength
    q_b_total_lowerdiag_array = [baseshear_lowerdiag[-1] *x for x in np.ones(5)]
#CELL I
    for m in range(len(zprimearray)):
        if m < stringercoordinates[1,0] /C_a *len(zprimearray) :
            B_z = 0*np.ones(5)
            B_y=0*np.ones(5)
        if stringercoordinates[1,0]/C_a *len(zprimearray) < m:
            B_z=stringercoordinates[1,0]*A_st*np.ones(5)
            B_y= stringercoordinates[1,1]*A_st*np.ones(5)
        baseshear_semicircle.extend(np.add([-SY[i]/I_zz*x for x in np.add(skinintegral_F(5,f3, previous_h, (m+1)/len(zprimearray) *h/4*np.pi ), B_y)], [-SZ[i]/I_yy*x for x in np.add(skinintegral_F(5,f4,previous_h,(m+1)/len(zprimearray)*h/4*np.pi),B_z)]))
        previous_h=(m+1)/len(zprimearray)*h/2*np.pi
    
    q_b_total_quartercircle_array= [baseshear_semicircle[-1] *x for x in np.ones(5)]
    for n in range(len(azarray)):
        if n > stringercoordinates[10,0]/C_a*len(azarray):
            B_z=stringercoordinates[0,0]*A_st
            B_y=stringercoordinates[0,1]*A_st
        if n< stringercoordinates[10,0]/C_a*len(azarray):
            B_z=np.ones(5)*A_st*(stringercoordinates[0,0]+stringercoordinates[10,0])
            B_y=np.ones(5)*A_st*(stringercoordinates[0,1]+stringercoordinates[10,1])
        baseshear_semicircle.extend(np.add(np.add([-SY[i]/I_zz *x for x in np.add(skinintegral_F(5,f3,previous_j, (n+1)/len(azarray)*h/4*np.pi + h/4*np.pi) , B_y)],[-SZ[i]/I_yy *x for x in np.add(skinintegral_F(5,f4, previous_j, (n+1)/len(azarray)*h/4*np.pi + np.pi *h/4) ,B_z)]),q_b_total_quartercircle_array ))
        previous_j=(m+1)/len(zprimearray)*h/4*np.pi+h/4*np.pi
    
    q_b_total_semicircle_array =[baseshear_semicircle[-1] *x for x in np.ones(5)]
    for o in range(len(yarray)):
        baseshear_spar_I.extend(np.add(np.add(q_b_total_semicircle_array,[-SY[i]/I_zz*x for x in (skinintegral_F(5,f9,previous_l,(o+1)/len(yarray)*h))]),[-SZ[i]/I_yy*x for x in skinintegral_F(5,f10,previous_l,(o+1)/len(yarray)*h)]))
        previous_l=(o+1)/len(yarray)*h
    q_b_total_spar_cell_I_array=[baseshear_spar_I[-1]*x for x in np.ones(5)]
#now that the base shear flows are found we use 3 more equations with 3 unknowns to find the redundant shear flow
#some constants that will be used in the 3x3 redundant shear flow matrix 
    A= np.pi * (h/2)**2
    B= 2*h*(C_a-h/2)
    C=0
    for c in range(len(baseshear_semicircle)-1):
        C=C+ h/2 * (baseshear_semicircle[c] + baseshear_semicircle[c+1])/2 * (np.pi *h/2)/len(baseshear_semicircle) + h/2 * np.sin(phi) * (baseshear_upperdiag[c] + baseshear_upperdiag[c+1])/2 * diagonallength/len(baseshear_upperdiag) + h/2* np.sin(phi) * (baseshear_lowerdiag[c]+baseshear_lowerdiag[c+1])/2 * diagonallength/len(baseshear_lowerdiag)
    D=1/(2*A_I*G)
    E=0
    for e in range(len(baseshear_semicircle)-1):
        E=E+ (baseshear_semicircle[e] + baseshear_semicircle[e+1])/2/t_sk *h/2 /len(baseshear_semicircle)*np.pi
    F=0
    for f in range(len(baseshear_spar_II)-1):
        F=F+ ((baseshear_spar_I[f] + baseshear_spar_I[f+1])/2 - (baseshear_spar_II[f]+baseshear_spar_II[f+1])/2 )/t_sp /len(baseshear_spar_II) * h
    G=1/(2*A_II*G)
    H=0
    for h1 in range(len(baseshear_upperdiag)-1):
        H=H+ (baseshear_upperdiag[h1]+baseshear_upperdiag[h1+1])/2/t_sk/len(baseshear_upperdiag)* diagonallength
    I=0
    for i1 in range(len(baseshear_spar_I)-1):
        I=I+ ((baseshear_spar_II[i1] + baseshear_spar_II[i1+1])/2 - (baseshear_spar_I[i1]+baseshear_spar_I[i1+1])/2) /t_sp /len(baseshear_spar_I) * h
    J=0
    for j in range(len(baseshear_lowerdiag)-1):
        J=J+ (baseshear_lowerdiag[j] + baseshear_lowerdiag[j+1])/2 /t_sk/len(baseshear_lowerdiag) *diagonallength
#solving the matrix
    M1=np.zeros([3,3])
    M2=np.zeros([3,1])
    M3=np.zeros([3,1])
    M1[0,0] = A
    M1[0,1] = B
    M1[0,2] = 0
    M1[1,0] = D*np.pi*h/2 / t_sk +D*h/t_sp
    M1[1,1] = -D*h
    M1[1,2] = -1
    M1[2,0] = -G*h/t_sp
    M1[2,1] = G*h/t_sp + 2* G*diagonallength/t_sk 
    M1[2,2] = -1
    M3[0,0] = -C
    M3[1,0] = -D*E-D*F
    M3[2,0] = -G*(H+I+J)
    M2=np.linalg.solve(M1,M3)
    q_s_0_I = M2[0,0]
    q_s_0_II = M2[1,0]
    dtheta_dz = M2[2,0]

        
#final shear flows with base shear flow and redundant shear flow lists
    q_semicircle=[]
    q_upperdiag=[]
    q_lowerdiag=[]
    q_spar_I=[]
    q_spar_II=[]
    q_spar=[]
    
    q_semicircle= np.add(baseshear_semicircle, np.ones(50) * q_s_0_I)
    q_upperdiag= np.add(baseshear_upperdiag, np.ones(50)* q_s_0_II)
    q_lowerdiag=np.add(baseshear_lowerdiag, np.ones(50)* q_s_0_II)
    q_spar_I= np.add(baseshear_spar_I, np.ones(50)* q_s_0_I)
    q_spar_II=np.add(baseshear_spar_II, np.ones(50)* q_s_0_II)
    q_spar= np.add(q_spar_I,[-1*x for x in q_spar_II])
    
    
#shear centre from trailing edge
    
    
    SCM=0
    for c in range(len(q_semicircle)-1):
        SCM=SCM+ h/2 * (q_semicircle[c] + q_semicircle[c+1])/2 * (np.pi *h/2)/len(q_semicircle) + h/2 * np.sin(phi) * (q_upperdiag[c] + q_upperdiag[c+1])/2 * diagonallength/len(q_upperdiag) + h/2* np.sin(phi) * (q_lowerdiag[c]+q_lowerdiag[c+1])/2 * diagonallength/len(q_lowerdiag)
    z_sc=SCM/SY[10]

Sigma_x_upperdiag=[]
Sigma_x_lowerdiag=[]
Sigma_x_spar=[]
Sigma_x_semicircle=[]
zprime=[]
yprime=[]
VM=[]
q_full=[]
Sigma_x=[]
diag=np.linspace(0,diagonallength, 50)
spar=np.linspace(0,h,50)
quartercircle=np.linspace(np.pi/2,0,25)
z1=np.linspace(0,C_a-h/2,50)
y1=np.linspace(0,h/2,50)
z2=np.ones(100)*C_a-h/2
y2=np.linspace(-h/2,h/2,50)
y3=np.linspace(-h/2,0,50)
z3=np.linspace(C_a-h/2,0,50)

zn0= np.linspace(C_a,h/2,50)
yn0=np.linspace(0,h/2,50)
yn1=np.linspace(h/2,-h/2,50)
zn1=np.ones(50)*h/2
zn2=np.linspace(h/2, C_a,50)
yn2=np.linspace(-h/2,0,50)
zn3=[]
for z in range(25):
    zn3.append( h/2- np.sin(z/25 * np.pi/2 )*h/2)
yn3=[]
for w in range(25):
    yn3.append(np.cos(w/25*np.pi/2)*h/2)
zn4=[]
for u in range(25):
    zn4.append(h/2-np.cos(u/25*np.pi/2)*h/2)
yn4=[]
for v in range(25):
    yn4.append(-np.sin(v/25*np.pi/2)*h/2)
#DIRECT STRESS
for i in range(1):
    for j in range(len(diag)):
        Sigma_x_upperdiag.append((M_z[i]*I_yy* y1[j] + M_y[i]*I_zz*z1[j])/(I_yy*I_zz))
    for k in range(len(diag)):    
        Sigma_x_spar.append((M_z[i]*I_yy* y2[k] + M_y[i]*I_zz*z2[k])/(I_yy*I_zz))
    for l in range(len(diag)):
        Sigma_x_lowerdiag.append((M_z[i]*I_yy* y3[l] + M_y[i]*I_zz*z3[l])/(I_yy*I_zz))
    for m in range(len(quartercircle)):
        Sigma_x_semicircle.append(M_z[i]*I_yy* np.sin((len(quartercircle)-m)/len(quartercircle)* np.pi/2)*h/2 + M_y[i]*I_zz * (C_a + np.cos((len(quartercircle)-m)/len(quartercircle)*np.pi/2)*h/2))
    for n in range(len(quartercircle)):
        Sigma_x_semicircle.append(M_z[i]*I_yy* (-np.sin(n/len(quartercircle)*np.pi/2)*h/2 ) + M_y[i] *I_zz* (C_a + np.cos(n/len(quartercircle)*np.pi/2 *h/2 )))
#VM STRESS
    for j in range(len(diag)):
        VM.append(np.sqrt((Sigma_x_upperdiag[j]**2)/2 +3*q_upperdiag[j]**2))
        q_full.append(q_upperdiag[j])
        Sigma_x.append(Sigma_x_upperdiag[j])
        zprime.append(zn0[j])
        yprime.append(yn0[j])
    for k in range(len(spar)):    
        VM.append(np.sqrt((Sigma_x_spar[k]**2)/2 +3*q_spar[k]**2))
        zprime.append(zn1[k])
        yprime.append(yn1[k])
        q_full.append(q_spar[k])
        Sigma_x.append(Sigma_x_spar[k])
    for l in range(len(diag)):
        VM.append(np.sqrt((Sigma_x_lowerdiag[l]**2)/2 +3*q_lowerdiag[l]**2))
        q_full.append(q_lowerdiag[l])
        Sigma_x.append(Sigma_x_lowerdiag[l])
        zprime.append(zn2[l])
        yprime.append(yn2[l])

    for m in range(50):
        VM.append(np.sqrt((Sigma_x_semicircle[m]**2)/2 +3*q_semicircle[m]**2))
        q_full.append(q_semicircle[m])
        Sigma_x.append(Sigma_x_semicircle[m])
    for b in range(25):
        zprime.append(zn3[b])
        yprime.append(yn3[b])
    for s in range(25):
        zprime.append(zn4[s])
        yprime.append(yn4[s])
        
    


'''
m = plt.cm.ScalarMappable(cmap=plt.get_cmap("jet"))
m.set_array(np.array([np.min(VM), np.max(VM)]))
plt.colorbar(m)
plt.scatter(zprime,yprime,c=VM, vmin=np.min(VM), vmax=np.max(VM), s=5, cmap="jet")
'''
m0 = plt.cm.ScalarMappable(cmap=plt.get_cmap("jet"))
m0.set_array(np.array([np.min(Sigma_x), np.max(Sigma_x)]))
plt.colorbar(m0)
plt.scatter(zprime,yprime,c=Sigma_x, vmin=np.min(Sigma_x), vmax=np.max(Sigma_x), s=5, cmap="jet")
'''
m1 = plt.cm.ScalarMappable(cmap=plt.get_cmap("jet"))
m1.set_array(np.array([np.min(q_full), np.max(q_full)]))
plt.colorbar(m1)
plt.scatter(zprime,yprime,c=q_full, vmin=np.min(q_full), vmax=np.max(q_full), s=5, cmap="jet")
'''
plt.show()

print(q_s_0_I)
print(q_s_0_II)    
print(z_sc)