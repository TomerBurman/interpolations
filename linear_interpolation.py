def linear_interpolation(table,x):
    '''
    linear_interpolation - builds a line between 2 consequent points and gives an assumption of f(x)
    :param table: iterative sequence of tuples where x corresponds to f(x).
    :param x: desired x
    :return: f(x) assumption if x is in bounds, else returns None.
    '''
    sorted_table = sorted(table) #
    print(sorted_table)
    point = sorted_table[0]
    if x < sorted_table[0][0] or x > sorted_table[-1][0]:
        print("You're asking for extra polation.")
        return None
    for point2 in sorted_table:

        if x > point2[0]:
            point = point2
            continue
        if x == point2[0]:
            return point2[1]
        else:
            print(point,point2)
            f_x = (1 / (point[0] - point2[0])) * ((x * (point[1] - point2[1]) + (point2[1] * point[0] - point[1] * point2[0])))
            return f_x




print(linear_interpolation([(0,0), (1,0.8415), (2,0.9093), (4,-0.7568), (6,-0.2794), (3,0.1411), (5,-0.9589)],-0.5))


