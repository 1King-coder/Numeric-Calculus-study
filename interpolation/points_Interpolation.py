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
import sympy as sym
from copy import deepcopy
import math



x = sym.symbols('x') # x icognito to build expressions

def memoization (memory: dict = {}, key: str = '', value = None, clear: bool = False):
    """
    Memoization function that uses external global variable to store useful data
    temporaly.
    """
    if clear:
        memory.clear()
        return
    
    if key in memory:
        return memory[key]
    
    memory[key] = value

    return value

def sep_tuple_values (tuples_list: list, index: int) -> list:
    # Function to separate the tuples coeficients by its indexes
    # like separating all the x_coordinates from a list of points

    return [tuple_vals[index] for tuple_vals in tuples_list]

class Interpolation_Methods:
    """
    Class with different methods to make interpolation functions of a given group
    of points in R²
    """

    def __init__(self, points: list, original_function: 'Func' = None, precision: int = 5) -> None:
        self.points = points
        self.original_func = original_function
        self.prec = precision
        self.num_of_points = len(points)

        self.__memo = {} # memory for using memoization safetly

        self.lagrange_iterations = list()
        self.vandermond_iterations = dict()
        self.newton_iterations = dict()

    @property
    def x_values (self) -> float:
        # list of x coordinates values.
        return sep_tuple_values(self.points, 0)
    
    @property
    def y_values (self) -> float:
        # list of y coordinates values.
        return sep_tuple_values(self.points, 1)
    
    @property
    def vandermond_x_matrix (self) -> 'Matrix':
        # Build vandermond matrix with the x coordinates
        return Matrix.map_matrix(
            lambda i, j: float(self.x_values[i])**j,
            self.num_of_points,
            self.num_of_points
        )
    
    @property
    def newton_iterations (self):
        return self.__newton_iterations
    
    @newton_iterations.setter
    def newton_iterations (self, value):
        """
        Factorate the newton iterations dictionary for a better view of the
        iterations by separating the operations by order.
        """
        if not value:
            self.__newton_iterations = value
            return
        

        # Transforms the key 'f[x1,...,xn]' into the array ['x1',...,'xn']
        # to make possible the use of len and define the order of the operation
        retrieve_list_from_key = lambda key: key[2:len(key) - 1].strip().split(',')
        
        new_newton_iterations_dict = {
            f'Order {i}': []
            for i in range(self.num_of_points)
        }

        for key, val in value.items():

            factor_order = len(retrieve_list_from_key(key)) - 1

            new_newton_iterations_dict[f'Order {factor_order}'].append(
                (key, val)
            )

        self.__newton_iterations = new_newton_iterations_dict
    
    def pprint_newton_iterations (self) -> None:
        # pretty prints the newton iterations

        if not self.newton_iterations:
            print("You need to execute the Newton's interpolation method first")
            return
        
        for key, val in self.newton_iterations.items():
            string = f'{key}:\n'
            for operation, restult in val:
                string += f'\t{operation} = {restult}\n'
            print(string)

    def build_interpolation_function (self, coefficients: list) -> 'Func':
        # This function builds the interpolation function witgh the given coefficients
        func = 0

        for degree, coeff in enumerate(coefficients):
            coeff = round(coeff, self.prec)
            func += coeff * x **degree
        
        return Func(func, x)
    
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
                    sym.Rational(1 , (float(self.x_values[index]) - float(self.x_values[j]))), # 1 / (xi - xj)
                    evaluate=False # ensure the form and do not create any fraction inside (x - xj)
                )

        self.lagrange_iterations[index][f'L{index}(x{index})'] = polynomial

        return sym.expand(polynomial)

    def lagrange_method (self) -> Func:
        """
        Lagrange method to build a interpolation polynomial in the form:
        Pn(x) = Σ(yi*Li(xi)) from i to n, where n is the number of points given.
        """
        inter_polynomial = 0
        for i in range(self.num_of_points):
            
            self.lagrange_iterations.append({'iter': i})

            inter_polynomial += self.y_values[i] * self.lagrangian_coeff(i)
        
        return Func(inter_polynomial, x)

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
                self.__memo,
                list_of_points_str,
                self.y_values[self.x_values.index(x_coordinates[0])]
            )
        
        if list_of_points_str in self.__memo:
            return self.__memo[list_of_points_str]
        
        result = memoization(
            self.__memo,
            list_of_points_str,
            round((self.newton_operator(x_coordinates[1:]) - # f[xi+1;...;[xn]]
            self.newton_operator(x_coordinates[:num_of_points-1])) / # f[xi;...;xn-1]
            (x_coordinates[num_of_points - 1] - x_coordinates[0]), self.prec) # (xn - xi)
        )

        self.newton_iterations = deepcopy(self.__memo) # saves the iterations

                
        return result
        
    def newton_method (self) -> 'Func':
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
        
        # Clear memory after usage
        memoization(clear=True)
        
        return Func(inter_polynomial, x)
        
    
class Interpolate (Interpolation_Methods):

    def __init__(self, points: list, original_function: 'Func' = None, precision: int = 5) -> None:
        super().__init__(points=points, precision=precision)

        self.original_func = original_function
        self.__memo = {} # memory for using memoization safetly

        self.corollary_1_iterations: dict = dict()
        self.corollary_2_iterations: dict = dict()
        self.corollary_3_iterations: dict = dict()

        self.derivatives: dict = dict() # Allows us to retrieve the derivatives
        # calculated to use some method

    def M_factor (self, n: int) -> float:
        """
        Calculates the M[n] = max(|f^(n+1)(x)|), where f^(n + 1) corresponds to (n+1)º derivative
        of the function.
        """
        # Builds the (n + 1)º derivative to use in the corollary 
        for derivative in range(1, n + 1):
            self.derivatives['f' + ("'"*derivative) + '(x)'] = self.original_func.dfunc(derivative)

        derivative_expression = self.original_func.dfunc(n)
        derivative_func = Func(derivative_expression, x)

        # Calculates the derivatives value for every point given
        derivative_values = [
            abs(derivative_func(x_coord))
            for x_coord in self.x_values
        ]

        return max(derivative_values)
    
    def A_factor (self, order: int) -> float:     
        """
        Calculates the A = max(|f[X0;...;Xorder+1]|).
        A is the highest value of the newton's operator values from the order n + 1
        where n is the interpolation polynomial degree.
        Futhermore, x0 -> xorder+1 will be the x coordinates where the interval between them
        is equal. 
        """

        # if the number of points is only 1 higher than the polynom degree,
        # Return already the order n + 1 value 
        if self.num_of_points == order:
            return abs(self.newton_operator(self.x_values[:order]))
        
        # Create the newton_iterations dicts
        self.newton_operator(self.x_values)

        # Takes all the modules of the Newton's operator values from the order
        order_operators_values = [abs(value) for _, value in self.newton_iterations[f'Order {order}']]

        memoization(self.__memo, clear=True) # Clear memory after usage

        return max(order_operators_values)     
    
    def abs_productory (self, limit: int = 0, x_coordinates: list = []) -> 'sym.Mul':
        # Builds the productory (Π |x - xi| from i = 0 to n) expression

        if not limit:
            limit = self.num_of_points - 1

        if not x_coordinates:
            x_coordinates = self.x_values

        productory = 1
        for x_coord in x_coordinates[:limit + 1]:
            # Create the expression of the productory
            productory *= sym.Abs(x - x_coord)

        return productory

    def corollary_1 (self, polynomial_degree: int, x_coordinate: float = None):
        """
        Corollary 1 for calculating the discrepancy between the interpolation polynom
        and the original function. 
        Err(x) <= (Π |x - xi| from i = 0 to n) * M[n + 1] / (n + 1)!
        M[n + 1] = max(|f^(n+1)(x)|)
        Where f^(n + 1) corresponds to (n+1)º derivative
        of the function.

        Ps.: 1 - If the specific x_coordinate where you want to calculate the error
             is not given, it will return the error function for every x for this
             corollary.

             2 - Note that this corollary is an estimative for the error
             in the x interval where the given points is contained. 

             3 - This corollary is most used when you only have the expression of
             the original function.
        """
        if not self.original_func:
            print('Needs the original function to use this corollary.')
            return
        


        # Gets the max(f^(n+1)(x))
        M_value = self.M_factor(polynomial_degree + 1)

        # Builds the Corollary expression for any error
        corollary_1_error_expression = self.abs_productory() * (M_value / math.factorial(polynomial_degree + 1))
        self.corollary_1_iterations = {
            f'M({polynomial_degree + 1})': M_value,
            'Err(x)': corollary_1_error_expression
        }

        # Turns the expression into a function to make calculating it's values easier.
        corollary_1_error_function = Func(corollary_1_error_expression, x)

        if x_coordinate:
            err = round(corollary_1_error_function(x_coordinate), self.prec)
            self.corollary_1_iterations[f'Err({x_coordinate})'] = err
            return err
        
        return corollary_1_error_function
    
    def corollary_2 (self, polynomial_degree: int):
        """
        Corollary 2 for calculating the discrepancy between the interpolation polynom
        and the original function.
        Err(x) <=  (h^(n+1) . M[n+1]) / 4 * (n+1)
        where h = x[i] - x[i-1], ∀i (interval between every x coordinates given)

        Ps.: 1 - This corollary only aplies when the intervals between every
             x coordinate given is equal and we have the original function expression.

             2 - Note that this corollary is an estimative for the error
             in the x interval where the given points is contained.
        """
        if not self.original_func:
            print('Needs the original function to use this corollary.')
            return
        
        n = polynomial_degree # just to simplify and not polute the formulas

        M_value = self.M_factor(n + 1)

        interval_between_points = abs(self.x_values[1] - self.x_values[0])

        corollary_2_error = round(
            ((interval_between_points**(n+1)) * M_value)/ (4*(n+1)),
            self.prec
        )

        self.corollary_2_iterations = {
            f'M({n + 1})': M_value,
            'Err(x)': corollary_2_error
        }
        
        return corollary_2_error
    
    def corollary_3 (self, polynomial_degree: int, x_coordinate: float = None):
        """
        Corollary 3 for calculating the discrepancy between the interpolation polynom
        and the original function.
        Err(x) ≈  (Π |x - xi| from i = 0 to n) * A
        Where A = max(|f[X0;...;Xorder+1]|) and A is the highest value of the Newton's
        operator values from the order n + 1
        where n is the interpolation polynomial degree.

        Ps.: 1 - This corollary only works if the number of points given is higher than the
             interpolation polynom degree.

             2 - Note that this corollary is an estimative for the error
             in the x interval where the given points is contained, but this corollary.
             points to and aproximation, not an estimative.

             3 - This corollary is most used when we do not have the original function expresison. 
        """
        
        if self.num_of_points < polynomial_degree:
            print("You need to give futher points to use this corollary.")
            return
        
        A_value = self.A_factor(polynomial_degree + 1)

        self.corollary_3_iterations['Operators'] = self.newton_iterations
        self.corollary_3_iterations['A_factor'] = A_value
        
        # Takes the points where the interval between them is equal
        best_points = set() # set type to avoid repeating any coordinate
        
        for i in range(self.num_of_points - 2):
            if abs(self.x_values[i+1] - self.x_values[i]) == abs(self.x_values[i+2] - self.x_values[i+1]):
                best_points = best_points.union(set(self.x_values[i:i+3]))
        
        best_points = list(best_points) # get back to list type

        self.corollary_3_iterations['Best_points'] = best_points
        
        corollary_3_expresison = self.abs_productory(polynomial_degree + 1, best_points) * A_value
        
        self.corollary_3_iterations['Err(expr)'] = corollary_3_expresison

        corollary_3_function = Func(corollary_3_expresison, x)

        if x_coordinate:
            err = round(corollary_3_function(x_coordinate), self.prec)
            self.corollary_3_iterations[f'Err({x_coordinate})'] = err
            return err
        
        return corollary_3_function


    

if is_main:
    points_1 = [
        (0, 1),
        (1, 6),
        (2, 5),
        (3, -8),
    ]

    points_2 = [
        (-1, 1.6988),
        (-0.4, 1.6916),
        (-0.2, 1.8214),
        (0, 2),
        (0.3, 2.3408),
        (1, 3.3817),
    ]
    points_3 = [
        (-1, 0.1988),
        (-0.4, 0.6174),
        (0.1, 1.0996),
        (1, 1.4687),
    ]

    points_4 = [
        (-1, 0.1988),
        (-0.5, 0.5323),
        (0, 1),
        (0.5, 1.4469),
        (1, 1.4687),
    ]

    inter = Interpolate(points_1)

    # print (inter.lagrange_method(), sep="\n")
    # print (*inter.lagrange_iterations, sep="\n")

    # print (inter.vandermond_method(), sep="\n")
    # print (*inter.vandermond_iterations, sep="\n")

    # print (inter.lagrange_method())
    # print (inter.vandermond_method())
    or_func = Func(sym.exp(x)*sym.cos(x), x)

    err = Interpolate(points_3, or_func, precision=10)
    
    print (err.newton_method())