def lagrange_interpolation(table,x):
    '''
    lagrange_interpolation - builds L_i polynomial for every point in the table
    sums L_i * f(x_i) to with desired x in P(x).
    :param table: sequence of tuples where x corresponds to f(x).
    :param x: desired x
    :return: P(x) value. (estimation for f(x))
    '''

    for point in table:
        if x == point[0]:
            return point[1]
    p = 0
    for point in table:
        p += build_L(table, x, point[0]) * point[1]
    return p


def build_L(table, x, x_i):
    '''
    build_L - builds Lagrange polynomial - if we have n points builds L_polynomial of degree n-1.
    for a given point iterating through
    all points in the table.
    :param table: sequence of tuples where x corresponds to f(x).
    :param x:
    :param x_i:
    :return: L(x) value.
    '''
    lagrange = 1
    for point in table:
        if point[0] != x_i:
            lagrange *= (x-point[0]) / (x_i - point[0])
    return lagrange

print(lagrange_interpolation([(1, 1), (2, 0), (4, 1.5)],3))