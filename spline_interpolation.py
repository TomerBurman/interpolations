import math
import numpy as np
from Gaussian_matrix_solver import *

isNatural = True # indicating which method to use. either natural or full.


def set_is_natural(boolean):
    global isNatural
    isNatural = boolean


def spline_interpolation(table, x, f_tagZero=0, f_tagN=0):
    '''
    spline_interpolation - uses isNatural to decide which method to use. for natural spline isNatural is True, else it's false.
    builds a matrix with size NxN of the form Meu_i * M_i-1 + 2 * M_i + Lambda_i * M_i+1 = d_i
    solves the equations, finds M's and then builds s_x
    :param table: list of tuples where x corresponds to f(x). -> (x,f(x))
    :param x: required s(x)
    :param f_tagZero: f_tag(0)
    :param f_tagN: f_tag(N)
    :return: s(x)
    '''
    for point in table:  # if x is a point in the table.
        if x == point[0]:
            return point[1]
    sorted_table = sorted(table)
    matrix, h_list = generate_matrix(table)
    d_list = calc_d(h_list, table, f_tagZero, f_tagN)
    #Using numpy gives about the same results.
    # A =np.array(matrix)
    # B = np.array(d_list)
    # m_list = np.linalg.solve(A,B)
    m_list = gauss_seidel_solver(mergeMetrix(matrix, list(map(lambda x: [x], d_list))))
    i = 0
    if x < sorted_table[0][0]:  # extrapolation when input is smaller than first point.
        s_x = ((((sorted_table[1][0] - x) ** 3) * m_list[0]) + (((x - sorted_table[0][0]) ** 3) * m_list[1])) / (6 * h_list[0]) \
              + (((sorted_table[1][0] - x) * sorted_table[0][1] + (x - sorted_table[0][0]) * sorted_table[1][1]) / h_list[0]) - \
              ((((sorted_table[1][0] - x) * m_list[0]) + ((x - sorted_table[0][0]) * m_list[1])) / 6) * h_list[0]

    elif x > sorted_table[-1][0]:  # extrapolation when input is larger than last point
        s_x = ((((sorted_table[-1][0] - x) ** 3) * m_list[-2]) + (((x - sorted_table[-2][0]) ** 3) * m_list[-1])) / (6 * h_list[-1]) \
              + (((sorted_table[-1][0] - x) * sorted_table[-2][1] + (x - sorted_table[-2][0]) * sorted_table[-1][1]) / h_list[-1]) - \
              ((((sorted_table[-1][0] - x) * m_list[-2]) + ((x - sorted_table[-2][0]) * m_list[-1])) / 6) * h_list[-1]
    else:  # x is between 2 points.
        while sorted_table[i][0] < x:
            i += 1
        i -= 1
        s_x = ((((sorted_table[i + 1][0] - x) ** 3) * m_list[i]) + (((x - sorted_table[i][0]) ** 3) * m_list[i + 1])) / (6 * h_list[i]) +\
              (((sorted_table[i + 1][0] - x) * sorted_table[i][1] + (x - sorted_table[i][0]) * sorted_table[i+1][1]) / h_list[i]) - \
              ((((sorted_table[i + 1][0] - x) * m_list[i]) + ((x - sorted_table[i][0]) * m_list[i + 1])) / 6) * h_list[i]
    return s_x


def calc_h(table):
    '''
    calc_h - builds h_list which are h_i = x_i+1 - x_i
    :param table: list of tuples where x corresponds to f(x). -> (x,f(x))
    :return: h_list
    '''
    h_list = list()
    for i in range(len(table) - 1):
        h_list.append(table[i + 1][0] - table[i][0])
    return h_list


def calc_lambda(h_list):
    '''
    calc_lambda - builds lambda list where in Natural spline lambda_0 = 0 and in Full spline lambda_0 = 1
    lambda_i = h_list[i+1] / (h_list[i+1] + h_list[i])
    :param h_list: distances list
    :return: lambda list
    '''
    lambda_list = list()
    if isNatural:
        lambda_list.append(0)
    else:
        lambda_list.append(1)
    for i in range(len(h_list) - 1):
        lambda_list.append(h_list[i + 1] / (h_list[i + 1] + h_list[i]))
    return lambda_list


def calc_meu(lambda_list):
    '''
    calc_meu - builds meu_list where meu_i = 1-lambda_list[i] and in Natural spline meu_n = 0, in Full spline meu_n = 1
    :param lambda_list: lambda list
    :return: meu_list
    '''
    meu_list = list()
    for i in range(len(lambda_list)):
        meu_list.append(1 - lambda_list[i])
    if isNatural:
        meu_list.append(0)
    else:
        meu_list.append(1)
    return meu_list


def calc_d(h_list, table, f_tagZero, f_tagN):
    '''
    calc_d - calculates d_i for spline_interpolation matrix
    :param h_list: distances table
    :param table: table of tuples where x corresponds to f(x)
    :param f_tagZero:
    :param f_tagN:
    :return:
    '''
    d_list = list()
    #adding first element
    if isNatural:
        d_list.append(0)
    else:
        d_list.append((6 / h_list[0]) * ((((table[1][1] - table[0][1])) / h_list[0]) - f_tagZero))
    for i in range(1, len(h_list)):  # table is larger than h_list, we dont need -1
        d_list.append((6 / (h_list[i] + h_list[i - 1])) * (
                ((table[i + 1][1] - table[i][1]) / h_list[i]) - ((table[i][1] - table[i - 1][1]) / h_list[i - 1])))
    #adding last element
    if isNatural:
        d_list.append(0)
    else:
        d_list.append((6 / h_list[-1]) * (f_tagN - ((table[-1][1] - table[-2][1]) / h_list[0])))
    return d_list


def generate_matrix(table):
    '''
    builds NxN matrix in the form of :
            -------N------
    | (  2   lambda_0   0)
    N (  meu_1   2   lambda_0)   = >  Meu_i * M_i-1 + 2 * M_i + Lambda_i * M_i+1 = d_i
    | (  0   meu_2   2)
    :param table: table of tuples where x corresponds to f(x)
    :return: NXN matrix and h_list
    '''
    size = len(table)
    matrix = i_matrix_gen(size, size + 1)
    h_list = calc_h(table)
    lambda_list = calc_lambda(h_list)
    meu_list = calc_meu(lambda_list)
    matrix[0][0] = 2
    for i in range(1, size):
        matrix[i][i] = 2
        matrix[i][i - 1] = meu_list[i]
        matrix[i - 1][i] = lambda_list[i - 1]
    return matrix, h_list
#
# if __name__ == "__main__":
#     # Natural
#     print(
#         f'The partial solution is: {spline_interpolation([(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)], math.pi / 3)}')
#     set_is_natural(False)
#     # Unnatural
#     print(
#         f'The full solution is: {spline_interpolation([(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)], math.pi / 3, 1, 0)}')
#     set_is_natural(True)


