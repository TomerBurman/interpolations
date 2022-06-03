import math

from Gaussian_matrix_solver import *

isNormal = True


def natural_spline(table, x):
    global isNormal
    for point in table: # if x is a point in the table.
        if x == point[0]:
            return point[1]
    isNormal = True
    sorted_table = sorted(table)
    matrix, h_list = generate_matrix(table)
    d_list = calc_d(h_list, table)
    m_list = gauss_seidel_solver(mergeMetrix(matrix, list(map(lambda x: [x], d_list))))
    print(matrix)
    print(f'Heres d{d_list}')
    print(m_list)
    i = 0
    if x < sorted_table[0][0]:  # extrapolation when input is smaller than first point.
        s_x = ((((table[1][0] - x) ** 3) * m_list[0]) + (((x - table[0][0]) ** 3) * m_list[1])) / (6 * h_list[0]) \
              + ((((table[1][0] - x) * table[0][1]) + (x - table[0][0] * table[1][1])) / h_list[0]) - \
              ((((table[1][0] - x) * m_list[0]) + ((x - table[0][0]) * m_list[1])) / 6) * h_list[0]

    elif x > sorted_table[-1][0]: # extrapolation when input is larger than last point
        s_x = ((((table[-1][0] - x) ** 3) * m_list[-2]) + (((x - table[-2][0]) ** 3) * m_list[-1])) / (6 * h_list[-1]) \
              + ((((table[-1][0] - x) * table[-2][1]) + (x - table[-2][0] * table[-1][1])) / h_list[-1]) - \
              ((((table[-1][0] - x) * m_list[-2]) + ((x - table[-2][0]) * m_list[-1])) / 6) * h_list[-1]
    else: # x is between 2 points.
        while sorted_table[i][0] < x:
            i += 1
        i -= 1
        s_x = ((((table[i+1][0] - x) ** 3) * m_list[i]) + (((x - table[i][0]) ** 3) * m_list[i+1])) / (6 * h_list[i]) +\
              ((((table[i+1][0] - x) * table[i][1]) + (x - table[i][0] * table[i+1][1])) / h_list[i]) -\
              ((((table[i+1][0] - x) * m_list[i]) + ((x - table[i][0]) * m_list[i+1])) / 6) * h_list[i]
    return s_x

def calc_h(table):
    h_list = list()
    for i in range(len(table) - 1):
        h_list.append(table[i + 1][0] - table[i][0])
    return h_list


def calc_lambda(h_list):
    lambda_list = list()
    if isNormal:
        lambda_list.append(0)
    else:
        lambda_list.append(1)
    for i in range(len(h_list) - 1):
        lambda_list.append(h_list[i + 1] / (h_list[i + 1] + h_list[i]))
    return lambda_list


def calc_meu(lambda_list):
    meu_list = list()
    for i in range(len(lambda_list)):
        meu_list.append(1 - lambda_list[i])
    if isNormal:
        meu_list.append(0)
    else:
        meu_list.append(1)
    return meu_list


def calc_d(h_list, table, f_tagZero=0, f_tagN=0):
    '''
    calc_d - calculates d_i for spline_interpolation matrix
    :param h_list: distances table
    :param table: table of tuples where x corresponds to f(x)
    :param f_tagZero:
    :param f_tagN:
    :return:
    '''
    d_list = list()
    if isNormal:
        d_list.append(0)
    else:
        d_list.append((6 / h_list[0]) * (((table[1][1] - table[0][1])) / h_list[0]) - f_tagZero)
    for i in range(1, len(h_list)):  # table is larger than h_list, we dont need -1
        d_list.append((6 / (h_list[i] + h_list[i - 1])) * (
                    ((table[i + 1][1] - table[i][1]) / h_list[i]) - ((table[i][1] - table[i - 1][1]) / h_list[i - 1])))

    if isNormal:
        d_list.append(0)
    else:
        d_list.append((6 / h_list[-2]) * (f_tagN - ((table[-1][1] - table[-2][1]) / h_list[0])))
    return d_list


def generate_matrix(table):
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


print(natural_spline([(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)], math.pi / 3))
