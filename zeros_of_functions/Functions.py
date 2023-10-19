from math import ceil, log, cos, sin, exp
from decimal import getcontext, Decimal
from random import randrange
import sympy as sym
import numpy as np


prec = 6
getcontext().prec = prec

def decimal (num: float) -> 'Decimal':
    # float to decimal conversor
    return Decimal(str(num))

class Func:
    """
    Class made to analise mathematical functions with sympy in a intuitive way
    and make it possible to visualize them in a graphic way
    """
    def __init__ (self, func: sym.core.add.Add, x: sym.core.symbol.Symbol) -> None:
        self.func = func
        self.x = x
    
    def f (self, value, function=None):
        # Makes getting the function value from a given X easy
        function_to_use = function if function else self.func

        if isinstance(value, np.ndarray):
            # Checks if the value is a ndarray and converts the results to
            # an array, making the use of it with matplotlib possible and easier
            return [float(function_to_use.subs({self.x: i}).evalf()) for i in value]

        return float(function_to_use.subs({self.x: value}).evalf())
    
    @property
    def degree (self) -> int:
        # Returns the polynomials degree. if it is no polynomial
        # it returns the highest exponent.
        
        degrees = []
        for element in self.func.args:
            degrees.append(element.as_base_exp()[1])

        return max(degrees)

    def __call__ (self, value, function=None):
        return self.f(value, function)
    
    def __str__ (self) -> str:
        return str(self.func)

    def __repr__ (self) -> str:
        return self.__str__()

    def dfunc (self, n: int = 1):
        # nº derivative of the given function
        return self.func.diff((self.x, n))
    
    def secant_factor (self, x0: 'Decimal', x1: 'Decimal') -> float:

        return (self.f(x1) - self.f(x0))/ (x1 - x0)
        
    def verify_derivate (self, a) -> bool:
        # Verify if the derivate of the function in x=a is 0
        if self.dfunc() == 0:
            print('Derivative = 0')
            return False
        
        return True

    def verify_zero (self, a: 'Decimal', b: 'Decimal') -> bool:
        # Verify if there is a 0 in a given interval

        if self.f(a)*self.f(b) > 0:
            print("Não há 0 da função no intervalo ou intervalo inválido para análise!")
            return False
        
        return True

    def bissection_method (self, a: float, b: float, TOL: float) -> dict:
        """
        Uses the bissection method to find the solution of the function
        with TOL's precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """
        a, b, TOL = decimal(a), decimal(b), decimal(TOL)

        Max_iterations = ceil(
            (decimal(log(b - a)) - decimal(log(TOL)))/decimal(log(2))
        )

        print(f'Num. Maximo de iterações: {Max_iterations}')

        if not self.verify_zero(a, b):
            return False
        
        iteration = 0
        iterations_list = []
        while iteration < Max_iterations + 1:

            x = decimal((a + b)/2)
            
            half_interval = decimal((b-a)/2)
            iterations_list.append(
                {
                    'iter': iteration,
                    f'a{iteration}': a,
                    f'b{iteration}': b,
                    f'x{iteration}': x,
                    f'f(x{iteration})': self.f(x),
                    'interval': b - a
                }
            )

            if self.f(x)*self.f(a) < 0:
                b = x
            else:
                a = x

            iteration += 1

            if self.f(x) == 0 or half_interval < TOL:                
                if iteration < Max_iterations:
                    iteration = Max_iterations - 1

        # Add x and y elements for use when building graphics with matplot  
        iterations_x, iterations_y = [], []                  
        for i, element in enumerate(iterations_list):
            iterations_x.append(element[f'x{i}'])
            iterations_y.append(element[f'f(x{i})'])

        return {
            'result': x,
            'iter_x': iterations_x,
            'iter_y': iterations_y,
            'iterations': iterations_list
        }
    
    def newton_method (self, a: float, TOL: float):
        """
        Uses the Newton's method to find the solution of the function 
        with TOL precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """

        a, TOL = decimal(a), decimal(TOL)

        if not self.verify_derivate(a):
            return False


        # Initial set
        iteration = 0
        iterations_list = []
        x = a
        f_x = self.f(x)
        
        # Max iterations to avoid infinity iterations when the method can not converge
        while abs(f_x) > TOL and iteration < 100:
            
            f_x = self.f(x)
            df_x = self.f(x, self.dfunc(1))

            x -= (f_x / df_x)

            iterations_list.append(
                {
                    'iter': iteration,
                    f'x{iteration}': x,
                    f'f(x{iteration})': f_x,
                    f"f'(x{iteration})":df_x
                }
            )

            if not self.verify_derivate(x):
                return False

            iteration += 1

            if self.f(x) == 0:
                break
            
        # x and y values for ploting graphics
        iter_x, iter_y = [], []
        
        for i, element in enumerate(iterations_list):
            iter_x.append(element[f'x{i}'])
            iter_y.append(element[f'f(x{i})'])

        return {
            'result': x,
            'iter_x': iter_x,
            'iter_y': iter_y,
            'iterations': iterations_list
        }
    
    def secant_method (self, a: float, b: float, TOL: float):
        """
        Uses the Newton's method to find the solution of the function 
        with TOL precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """

        x0, x1, TOL = decimal(a), decimal(b), decimal(TOL)

        sec_factor = self.secant_factor(x0, x1)

        if not sec_factor:
            return False
        

        iteration = 1

        f_x1 = self.f(x1)
        x2 = x1 - (f_x1 / sec_factor)

        iterations_list = []

        while abs(self.f(x2)) > TOL and iteration < 100:

            if not self.secant_factor(x0, x1):
                return False

            f_x1 = self.f(x1)
            sec_factor = self.secant_factor(x0, x1)

            x2 = x1 - (f_x1 / sec_factor)

            iterations_list.append(
                {
                    'iter': iteration,
                    f'x{iteration}': x2,
                    f'f(x{iteration})': self.f(x2),
                    f"secant factor ({iteration})": self.secant_factor(x0, x1)
                }
            )
            
            x0, x1 = x1, x2

            iteration += 1
            
            if self.f(x2) == 0:
                break
        
        # x and y values for ploting graphics
        iter_x, iter_y = [], []
        
        for i, element in enumerate(iterations_list):
            iter_x.append(element[f'x{i+1}'])
            iter_y.append(element[f'f(x{i+1})'])

        return {
            'result': x2,
            'iter_x': iter_x,
            'iter_y': iter_y,
            'iterations': iterations_list
        }
        
    def fake_position_method (self, a: float, b: float, TOL: float) -> dict:
        """
        Uses the fake position method to find the solution of the function
        with TOL's precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """
        a, b, TOL = decimal(a), decimal(b), decimal(TOL)

        if not self.verify_zero(a, b):
            return False

        # Expression used to obtain x_k in this method
        fake_pos_expr = lambda a_k, b_k: decimal(
            (self.f(b_k)*a_k - self.f(a_k)*b_k) /
            (self.f(b_k) - self.f(a_k))
        )

        iteration = 0
        iterations_list = []
        x = fake_pos_expr(a, b)
        
        while abs(self.f(x)) > TOL and iteration < 100:
            
            x = fake_pos_expr(a, b)
            
            iterations_list.append(
                {
                    'iter': iteration,
                    f'a{iteration}': a,
                    f'b{iteration}': b,
                    f'x{iteration}': x,
                    f'f(x{iteration})': self.f(x),
                    'interval': b - a
                }
            )

            if self.f(x)*self.f(a) < 0:
                b = x
            else:
                a = x

            iteration += 1

            if self.f(x) == 0:
                break

        # Add x and y elements for use when building graphics with matplot  
        iterations_x, iterations_y = [], []                  
        for i, element in enumerate(iterations_list):
            iterations_x.append(element[f'x{i}'])
            iterations_y.append(element[f'f(x{i})'])

        return {
            'result': x,
            'iter_x': iterations_x,
            'iter_y': iterations_y,
            'iterations': iterations_list
        }

    @property
    def critical_x_coordinates (self):
        """
        Method to obtain the max value of the function.
        """

        # Solves the expression f'(x) = 0 to find a critical point
        x_coordinate = sym.solve(self.dfunc(1), self.x)

        return x_coordinate # return and array with the critical points found
    
    def is_max (self, x_coordinate: float):
        # Verify if a given x coordinate is a critical max coordinate.
        return self.f(x_coordinate, self.dfunc(2)) < 0
    
    def global_max_min (self, option: str):
        """
        Return the global max or min function value, must send the option.
        if you want the global max, option='max', else option='min'
        """

        option = option.lower().strip()

        func_global_values = [
            float(self.f(x))
            for x in self.critical_x_coordinates
        ]

        global_max = max(func_global_values)
        global_min = min(func_global_values)

        return global_max if option == 'max' else global_min

    def find_x_for (self, f_value: float) -> float:
        """
        Method to find the relative x coordinate for a wanted f(x) value.
        """
        
        # Expression f(x) = value -> f(x) - value
        expression_to_solve = self.func - f_value

        result_x = sym.solve(expression_to_solve, self.x)[0]

        return result_x

if __name__ == '__main__':
    # tests
    x = sym.symbols('x')

    Vc = x**2 + 2*x - 1
    
    func = Func(
        Vc, x
    )
    
    print(
        func.dfunc(),sep='\n'
    )

    

    ...




