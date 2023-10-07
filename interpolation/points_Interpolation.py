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
    

    def build_interpolation_function (self, coefficients: list):
        # This function builds the interpolation function witgh the given coefficients
        # calculated by any method
        func = 0

        for degree, coeff in enumerate(coefficients):
            func += coeff * x **degree
        
        return func
    
    def vandermond_method (self):
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
    

if __name__ == '__main__':
    inter = Interpolate([
        (0, 1),
        (1, 6),
        (2, 5),
        (3, -8),
    ])

    print (inter.vandermond_method(), sep="\n")

