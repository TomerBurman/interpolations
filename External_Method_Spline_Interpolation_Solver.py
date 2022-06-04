import math

from scipy import interpolate

def f(x):
    x_points = [ 0, math.pi/6, math.pi/4, math.pi/2]
    y_points = [0,0.5,0.7072,1]

    tck = interpolate.splrep(x_points, y_points)
    return interpolate.splev(x, tck)

print(f(math.pi/5)) #print(spline_interpolation([(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)], math.pi / 3 ))



