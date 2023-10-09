import sys

sys.path.append('.\\linear_systems\\')

from Vectors import Vector
from Matrix import Matrix
from Linear_systems import Linear_System
import sympy as sym

global x

x = sym.symbols('x')




class Interpolate:
    """
    Class with different methods to make interpolation functions of a given group
    of points in R²
    """

    def __init__(self, points: list) -> None:
        self.points = points
        self.num_of_points = len(points)

        self.lagrangian_iterations = []
        self.vandermond_iterations = []
        self.newton_iterations = []

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
                polynomial *= (x - self.x_values[j]) / (self.x_values[index] - self.x_values[j])

        return sym.expand(polynomial)
        ...


    def lagrange_method (self):
        """
        Lagrange method to build a interpolation polynomial in the form:
        Pn(x) = Σ(yi*Li(xi)) from i to n, where n is the number of points given.
        """
        inter_polynomial = 0
        
        for i in range(self.num_of_points):
            inter_polynomial += self.y_values[i] * self.lagrangian_coeff(i)
        
        return inter_polynomial

    

if __name__ == '__main__':
    inter = Interpolate([
        (0, 1),
        (1, 6),
        (2, 5),
        (3, -8),
    ])

    print (inter.lagrange_method(), sep="\n")

