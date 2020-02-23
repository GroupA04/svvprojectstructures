from math import *
import numpy as np
from Aerodynamic_loading import xcoordinates, zcoordinates, q_disc
from numerical_functions import interpolate, integration

N_x = 41                                         # Number of points in span-wise direction (x-direction) --> should be 100?

q_x, z_cp = q_disc(zcoordinates)                       # This contains the values of q at each point along the span
# s_acum, s_out = integration(xcoordinates, 0, 1.611, q_x, 100)




def cumintegration(N_x , q, x):
    integral1 = []  # This will have the accumulated area at each point in x
    xcoordinates = x  # This will contain all the xcoordinates of our "processed" / "non-raw" data
    individ_integrals = []  # This will have the areas of all the little individual "bars" (will have one less value than N_x)
    for i in range(N_x - 1):                        # Minus 1 because we calculate area between i and i+1, so final point is not used
        print(i)
        s_i = (xcoordinates[i+1] - xcoordinates[i]) * ((q[i+1] - q[i])/2 + q[i])      # Calculates individual area
        print(s_i)                                                                           
        individ_integrals = np.append(individ_integrals,s_i)                                
        print(individ_integrals)

    print(individ_integrals)

    print("Then...")

    for j in range(N_x - 1):                                                                # Again, last point not used, thus minus 1
        i = j
        print("i = ",i)
        print("j = ",j)
        integral1_j = individ_integrals[i]                                                  # Start with area of individual bar in question
        while i >= 1:                                                                       # Add all areas "behind" it
            integral1_j = integral1_j + individ_integrals[i-1]                              # Keep adding cumulatively
            i = i - 1                                                                       # Reduce i to "select" previous individual bar
        print("integral1_j = ",integral1_j)
        
        integral1 = np.append(integral1 , integral1_j)                                      # Store the resultant cumulative area for every j (point in x direction)

    integral1 = np.append([0],integral1)
    
    print("integral1 = ", integral1)
    return integral1
# We will have to explain in the report that each time that we integrate we "lose" a point and this introduces a small error
# Taking more points in the x direction will reduce this error since the percentage of points "lost" goes down
cumintegration(N_x, q_x, xcoordinates)