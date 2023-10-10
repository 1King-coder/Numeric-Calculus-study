import sys

sys.path.append('.\\linear_systems\\')

from Vectors import Vector
from Matrix import Matrix
from Linear_systems import Linear_System
import sympy as sym
from copy import deepcopy

global x

x = sym.symbols('x')

memory = {}
def memoization (key: str = '', value = None, clear: bool = False):
    if clear:
        memory.clear()
        return print('Memory successfuly cleared.')
    
    if key in memory:
        return memory[key]
    
    memory[key] = value

    return value    

class Interpolate:
    """
    Class with different methods to make interpolation functions of a given group
    of points in R²
    """

    def __init__(self, points: list) -> None:
        self.points = points
        self.num_of_points = len(points)

        self.lagrange_iterations = list()
        self.vandermond_iterations = dict()
        self.newton_iterations = dict()

    @property
    def x_values (self) -> float:
        # list of x coordinates values.
        values = [point[0] for point in self.points]
        return values
    
    @property
    def y_values (self) -> float:
        # list of y coordinates values.
        values = [point[1] for point in self.points]
        return values
    
    @property
    def vandermond_x_matrix (self) -> 'Matrix':
        # Build vandermond matrix with the x coordinates
        return Matrix.map_matrix(
            lambda i, j: self.x_values[i]**j,
            self.num_of_points,
            self.num_of_points,
        )

    def build_interpolation_function (self, coefficients: list) -> 'sym.core.Add':
        # This function builds the interpolation function witgh the given coefficients
        func = 0

        for degree, coeff in enumerate(coefficients):
            func += coeff * x **degree
        
        return func
    
    def vandermond_method (self) -> 'sym.core.Add':
        """
        This method is based in solving a linear system Vc = y where V
        is the vandermond matrix where each line stands for each x coordinate
        and y is the vector with y coordinates.
        """

        vandermond_sys = Linear_System(
            self.vandermond_x_matrix,
            Vector(*self.y_values)
        )

        self.vandermond_iterations['lin_sys'] = vandermond_sys
        
        
        # Uses gaussian elimination method to solve the system and obtain the
        # coefficients
        coefficients = vandermond_sys.gauss_elimination_method()

        # Builds the function expression
        inter_function = self.build_interpolation_function(coefficients.entrys)

        return inter_function
    
    def lagrangian_coeff (self, index: int):
        """
        Function to calculate the Li(xi) lagrandian coefficient from
        the lagrange method where:
        Li(xi) = Π (x - xj) / (xi - xj) from j = 0 and j ≠ i(ndex) to n
        where n is the number of points given.
        """
        polynomial = 1

        for j in range(self.num_of_points):
            if index != j:
                # Here I have separeted the fraction by a multiplication to make
                # sure it is not evaluated and store the (x - xj) / (xi - xj) form. 
                polynomial *= sym.Mul(
                    (x - self.x_values[j]), # (x - xj)
                    sym.Rational(1 , (self.x_values[index] - self.x_values[j])), # 1 / (xi - xj)
                    evaluate=False # ensure the form and do not create any fraction inside (x - xj)
                )

        self.lagrange_iterations[index][f'L{index}(x{index})'] = polynomial

        return sym.expand(polynomial)

    def lagrange_method (self):
        """
        Lagrange method to build a interpolation polynomial in the form:
        Pn(x) = Σ(yi*Li(xi)) from i to n, where n is the number of points given.
        """
        inter_polynomial = 0
        for i in range(self.num_of_points):
            
            self.lagrange_iterations.append({'iter': i})

            inter_polynomial += self.y_values[i] * self.lagrangian_coeff(i)
        
        return inter_polynomial

    def newton_operator (self, x_coordinates: list):
        """
        Newton's operator f[xi;...;xn] for divisions between subtractions that appear
        in his method of interpolating points.
        f[xi;...;xn] = (f[xi+1;...;xn] - f[xi; xn-1]) / (xn - xi)
        As this function is recursive and has repeated terms while calculating, it was
        implemented the memoization method by a external function called "memoization".
        This memoization method helps to speedup the calculations.
        """

        num_of_points = len(x_coordinates)
        list_of_points_str = f"f{str(x_coordinates)}"

        if num_of_points == 1:
            return memoization(
                list_of_points_str,
                self.y_values[self.x_values.index(x_coordinates[0])]
            )
        
        if list_of_points_str in memory:
            return memory[list_of_points_str]
                
        return memoization(
            list_of_points_str,
            (self.newton_operator(x_coordinates[1:]) - # f[xi+1;...;[xn]]
             self.newton_operator(x_coordinates[:num_of_points-1])) / # f[xi;...;xn-1]
            (x_coordinates[num_of_points - 1] - x_coordinates[0]) # (xn - xi)
        )
        
    def newton_method (self):
        """
        Newton's method of building an interpolation polynominal function Pn(x).
        Pn(x) = f[x0] + f[x0;x1]*(x - x0) + f[x0;x1;x2]*(x - x0)*(x - x1) ... f[x0;...;xn]*...*(x - xn-1)
        where f[xi;...;xn] is the Newton's operator.
        """

        inter_polynomial = 0

        multiplications = 1

        for i in range(1, self.num_of_points + 1):
            # i-th term = f[xi-1;...;xi]
            newton_operator = self.newton_operator(self.x_values[:i])

            newton_operator *= multiplications
            
            # Multiplications = (x - xi-1)
            multiplications *= (x - self.x_values[i - 1])

            inter_polynomial += sym.expand(newton_operator)
        
        self.newton_iterations = deepcopy(memory)

        # Clear memory after usage
        memoization(clear=True)
        
        return inter_polynomial

    

if __name__ == '__main__':
    inter = Interpolate([
        (0, 1),
        (1, 6),
        (2, 5),
        (3, -8),
    ])

    # print (inter.lagrange_method(), sep="\n")
    # print (*inter.lagrange_iterations, sep="\n")

    # print (inter.vandermond_method(), sep="\n")
    # print (*inter.vandermond_iterations, sep="\n")
    dic = { str([(1, 2), (3, 4)]): 123}

    points = [0, 1, 2, 3]

    



    print (inter.lagrange_method())
    print (inter.vandermond_method())
    print (inter.newton_method())

