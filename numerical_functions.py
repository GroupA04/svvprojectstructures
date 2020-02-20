import numpy as np

#numerical function definitions


def spline(x, f):
    n = len(x)
    sp_start = []
    sp_slope = []

    for i in range(n - 1):
        sp_start = np.append(sp_start, f[i])
        slope = (f[i + 1] - f[i]) / (x[i + 1] - x[i])
        sp_slope = np.append(sp_slope, slope)

    return sp_start, sp_slope


def interpolate(x, f, x_target):
    n = len(x)
    sp_start, sp_slope = spline(x, f)
    left_p = n - 2
    for i in range(n):
        if x[left_p] > x_target:
            left_p = left_p - 1

        else:
            y = sp_start[left_p] + sp_slope[left_p] * (x_target - x[left_p])

    return y


