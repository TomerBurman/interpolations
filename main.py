from lagrange_interpolations import *
from spline_interpolation import *
from linear_interpolation import *
from nevilles_interpolation import *
from polynomial_interpolation import *

if __name__ == "__main__":
    table_point = [(0, 0), (math.pi / 6, 0.5), (math.pi / 4, 0.7072), (math.pi / 2, 1)]
    point = math.pi / 3

    print(f'for Point {point} with the table {table_point} interpolations are : ')
    print(f'1.Linear interpolation sol : {linear_interpolation(table_point, point)}\n')
    print(f'2.Polynomial interpolation sol : {polynomial_interpolation(table_point, point)}\n')
    print(f'3.lagrange interpolation sol : {lagrange_interpolation(table_point, point)}\n')
    print(f"4.neville's interpolation sol : {nevils_interpolation(table_point, point)}\n")
    print(f'5.natural spline interpolation sol : {spline_interpolation(table_point, point)}\n')
    set_is_natural(False)  # changing to full spline
    print(f'6.full spline interpolation sol : {spline_interpolation(table_point, point)}\n') #to use full spline add f_tag(0) and f_tag(N) to parameters.






