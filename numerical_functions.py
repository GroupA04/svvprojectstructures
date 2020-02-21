import numpy as np

#numerical function definitions


def spline(x, f):
    n = len(x)
    sp_st = []
    sp_sl = []

    for i in range(n - 1):
        sp_st = np.append(sp_st, f[i])
        slope = (f[i + 1] - f[i]) / (x[i + 1] - x[i])
        sp_sl = np.append(sp_sl, slope)

    return sp_st, sp_sl


def interpolate(x, f, x_t):
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

    return y


