import sys

is_main = __name__ == "__main__"

lin_sys_path = '.\\linear_systems\\' if is_main else '..\\linear_systems\\'
function_path = '.\\zeros_of_functions\\' if is_main else '..\\zeros_of_functions\\'

sys.path.append(lin_sys_path)
sys.path.append(function_path)

from Vectors import Vector
from Matrix import Matrix
from Linear_systems import Linear_System
from Functions import Func
import matplotlib.pyplot as plt
import sympy as sym
from copy import deepcopy
import math

x = sym.symbols('x') # x icognito to build functions

def build_guess_expression (parcels_functions: list, displacement: float = 0,coefficients: list= None) -> 'Func':
    """
    Builds your guessing expression as Σαi.gi(x)
    Here you send your guess functions as [g1(x), ..., gn(x)] to make sure the
    expression with α matches your order.
    if coefficients is given, it returns the function as Func object.
    """

    expression_parcels = [
        sym.symbols(f'α{i+1}') * parcel.func
        if not coefficients else
        coefficients[i] * parcel.func
        for i, parcel in enumerate(parcels_functions)
    ]

    full_expression = sym.Add(*expression_parcels, displacement, evaluate=False)

    return Func(full_expression, parcels_functions[0].x)

def sep_tuple_values (tuples_list: list, index: int) -> list:
    # Function to separate the tuples coeficients by its indexes
    # like separating all the x_coordinates from a list of points

    return [tuple_vals[index] for tuple_vals in tuples_list]

class Linear_adjust:

    def __init__(self, points: list, parcels_expressions: list = None, precision=5, custom_icognito: str= 'α') -> None:
        self.points = points
        self.num_of_points = len(self.points)
        self.prec = precision

        self.custom_icognito = custom_icognito

        self.coefficients: dict = {}

        self.parcels_funcs = parcels_expressions # The expressions from sympy will become Functions
        self.x_values = sep_tuple_values(points, 0)
        self.y_values = sep_tuple_values(points, 1)

    @property
    def parcels_funcs (self):
        return self.__parcels_funcs
    
    @parcels_funcs.setter
    def parcels_funcs (self, value):
        if not value:
            self.__parcels_funcs = value
            return
        
        x = sym.symbols('x')

        parcels_expr_as_functions = [Func(parcel, x) for parcel in value]

        self.__parcels_funcs = parcels_expr_as_functions

    @property
    def y_vector (self) -> 'Vector':
        return Vector(*self.y_values)
    
    def V_matrix (self) -> 'Matrix' or None:
        """
        Builds the matrix V where:
        
        V = |g0(x0) g1(x0) g2(x0) ... gn(x0)|
            |g0(x1) g1(x1) g2(x1) ... gn(x1)|
            |...                            |
            |g0(xm) ...               gn(xm)|

        """
        if not self.parcels_funcs:
            print('You can not build the matrix without a guess expression for the curve.')
            print('Try visualizing the points to make your guess.')
            return
        
        num_of_parcels = len(self.parcels_funcs)

        term_function = lambda i, j: round(
            self.parcels_funcs[j](self.x_values[i]),
            self.prec
        )

        list_matrix = Matrix.map_matrix (
            term_function,
            self.num_of_points,
            num_of_parcels,
        )

        matrix_V = Matrix(list_matrix)

        return matrix_V
    
    def calculate_coefficients (self) -> 'Vector':
        """
        
        """
        
        V_matrix = self.V_matrix()

        VT_matrix = V_matrix.transpose()

        if not V_matrix:
            print('You must give the guess expression functions to calculate the adjust.')
            return
        
        transf_matrix = VT_matrix * V_matrix

        transformed_y_vector = Vector(
            *((VT_matrix * self.y_vector.column_matrix).to_line_matrix())
        )

        adjust_linear_system = Linear_System(transf_matrix.matrix, transformed_y_vector, custom_icognito=self.custom_icognito)

        alphas_values = adjust_linear_system.gauss_elimination_method()

        return alphas_values
    
    def build_linear_fit_function (self, displacement: float = 0) -> 'Func':
        """
        
        """
        coefficients = self.calculate_coefficients()

        linear_fit_function = build_guess_expression(
            self.parcels_funcs,
            displacement,
            coefficients.entrys
        )

        return linear_fit_function



    def visualize_points (self, fig, axe) -> None:
        """
        Plots a figure with the points ploted in the cartesian plain to allow
        you visualize their pattern
        """

        axe.scatter(self.x_values, self.y_values, marker='.', color='r', s=50)

        fig.show()


if __name__ == '__main__':
    points_1 = [
        (0, -153),
        (0.25, 64),
        (0.5, 242),
        (0.75, 284),
        (1, 175),
    ]

    lin_ad = Linear_adjust(points_1)

    lin_ad.visualize_points()