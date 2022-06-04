from Gaussian_matrix_solver import *

def polynomial_interpolation(table,x):
    '''
    polynomial_interpolation - builds a degree n-1 polynomial from n points and gives an estimation on f(x) by P(x)
    :param table: iterative sequence of tuples where x corresponds to f(x)
    :param x:  desired x
    :return: P(x)
    '''

    for point in table:
        if x == point[0]:
            return point[1]
    sorted_table = sorted(table) # sorts the table
    matrix = []
    for point in sorted_table:
        temp = [point[0] ** i for i in range(len(sorted_table)) ] # building matrix coeffiecients.
        temp = temp + [point[1]] # adding y_i as solution.
        matrix.append(temp)
    coefficient_vector = gauss_seidel_solver(matrix,0.00001) # extracting a_0 ... a_n
    value = 0 # initiating value
    for power, elem in enumerate(coefficient_vector, 0):
        value += elem * ((x) ** power) #summing a_0 * (x)^0 + a_1 * (x)^1 ... + a_n * (x)^n-1
    return value
