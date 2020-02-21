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
            y = sp_st[x[0]] - sp_sl[x[0]]*(x[0] - x_t)
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

    s_cum = [s[0]]
    for i in range(1, len(s)):
        print(i)
        s_cum_i = s[i] + s_cum[i-1]
        s_cum = np.append(s_cum, s_cum_i)
    #outputs are the cumulitative area distribution and the total sum of the area under the curve
    return s_cum, sum(s)
