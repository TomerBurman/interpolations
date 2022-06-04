def nevils_interpolation(table, x):
    '''
    nevils_interpolation - numeric method to find assumption to f(x) by P(x). builds series of polynomials until we reach
    we reach 3rd degree polynomial
    :param table: array of tuples (x,f(x)) where x corresponds to f(x).
    :param x: x desired
    :return:
    '''
    if len(table) <= 3:
        print("Neville's interpolation recieves minimum of 3 points.")
        return
    sorted_table = sorted(table)
    for point in table: # if desired point is in our table
        if x == point[0]:
            return point[1]
    if x > sorted_table[-1][0] or x < sorted_table[0][0]: # extrapolation.
        print("You are asking for extrapolation.")

    i = 0
    while i < len(sorted_table) and x > sorted_table[i][0]:# discover x's location on our table.
        i += 1
    if i == len(table)-1: # between last and before last element
        p_x = build_P(x, sorted_table[i-3:i+1])[2][0]
    elif i == 1: # between first and 2nd element
        p_x = build_P(x, sorted_table[i-1:i+3])[2][0]
    elif i == len(sorted_table): # if i is greater than the last point
        p_x = build_P(x, sorted_table[i-4:i])[2][0]
    elif i == 0: # if i is smaller than the first element
        p_x = build_P(x, sorted_table[i:i+4])[2][0]
    else: # in between.
        p_x = build_P(x, sorted_table[i-2:i+2])[2][0]
    return p_x


def build_P(x,table):
    '''
    build_P - recieves table with 4 points, builds 3 linear equations from that builds 2 2nd degree polynomial equations
    and from that builds 1 3rd degree polynomial equation
    :param x: x desired
    :param table: array of tuples (x,f(x)) where x corresponds to f(x)
    :return: array of linear and polynomial equations.
    '''
    p = []
    p_1 = [((x-table[i][0]) * table[i+1][1] - (x - table[i+1][0]) * table[i][1]) / (table[i+1][0] - table[i][0]) for i in range(len(table)-1)] # linear equations
    p.append(p_1)
    p_2 = [((x - table[i][0]) * p_1[i+1] - (x - table[i+2][0]) * p_1[i]) / (table[i+2][0] - table[i][0]) for i in range(len(p_1)-1)] # 2nd degree polynomial equations.
    p.append(p_2)
    p_3 = [((x - table[i][0]) * p_2[i+1] - (x - table[i+3][0]) * p_2[i]) / (table[i+3][0] - table[i][0]) for i in range(len(p_2)-1)] # 3rd degree polynomial equations.
    p.append(p_3)
    return p



