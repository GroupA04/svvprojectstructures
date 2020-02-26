import numpy as np

#numerical function definitions


def spline(x, f): #sets up a linear spline model with inputs x coordinates and f(x)
    n = len(x)
    sp_st = []
    sp_sl = []

    for i in range(n - 1):
        sp_st = np.append(sp_st, f[i])
        slope = (f[i + 1] - f[i]) / (x[i + 1] - x[i])
        sp_sl = np.append(sp_sl, slope)
    #outputs are the starting points and slope of every spline
    return sp_st, sp_sl


def interpolate(x, f, x_t): #linear spline interpolation scheme with x_t as the target value
    n = len(x)
    sp_st, sp_sl = spline(x, f)
    left_p = n - 2
    y = 0
    for i in range(n):
        if x[left_p] > x_t:
            left_p = left_p - 1
        elif x_t < x[0]:
            y = sp_st[x[0]] + sp_sl[x[0]]*(x_t - x[0])
        else:
            y = sp_st[left_p] + sp_sl[left_p] * (x_t - x[left_p])

    #output is the value of f at any x
    return y


def integration(x, x_start,x_end,f,n): #integration using trapezoidial rule using interpolated model of x and f(x)
    s = []
    x_int = []
    f_int = []
    for i in range(1,n+1):
        x_i = (x_end-x_start)/n * i
        x_int = np.append(x_int,x_i)
        f_int = np.append(f_int, interpolate(x, f, x_i))
    for i in range(n-1):
        s_i = (x_int[i+1] - x_int[i]) * ((f_int[i+1] - f_int[i])/2 + f_int[i])
        s = np.append(s, s_i)

    s_acum = [s[0]]
    for i in range(1, len(s)):
        s_acum_i = s[i] + s_acum[i-1]
        s_acum = np.append(s_acum, s_acum_i)

    s_out = sum(s)
    #outputs are the cumulitative area distribution and the total sum of the area under the curve
    return s, s_acum, s_out

def average(case):
    case_ave = np.array([])
    for i in range(len(case)):
        ave1 = (case[:,2] + case[:,3])/2
        ave2 = (case[:,4] + case[:,5])/2
        case_ave = np.vstack([case[:,0],ave1,ave2])
    return case_ave


#=================================================Verification===================================================================================
#Test for Spline(x,f) function:
#Case1: input is a single line from (0,0) to (1,1)
# output should be starting point 0 and slope 1
x1 = [0,1]
f1 = [0,1]
sp_st1, sp_sl1 = spline(x1,f1)

if sp_st1 == 0 and sp_sl1 == 1:
    print('Spline case 1 Pass')
else:
    print('Spline case 1 Fail')

#Case 2: inputs are 2 lines from (1,2) to (2,4) to (4,10)
#output should give starting points [2,4] and slopes [2,3]
x2 = [1,2,4]
f2 = [2,4,10]
sp_st2, sp_sl2 = spline(x2, f2)

if (sp_st2[0] == 2) and (sp_st2[1] == 4) and (sp_sl2[0] == 2) and (sp_sl2[1] == 3):
    print('Spline case 2 Pass')
else:
    print('Spline case 2 Fail')

#Case 3: invalid 0 input
x3 = [0,0,0]
f3 = [1,1,1]

