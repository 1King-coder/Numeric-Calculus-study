from math import ceil, log, cos, sin
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
        if isinstance(value, np.ndarray):
            # Checks if the value is a ndarray and converts the results to
            # an array, making the use of it with matplotlib possible and easier
            return [(function if function else self.func).subs({self.x: i}) for i in value]

        return (function if function else self.func).subs({self.x: value})
    
    @property
    def dfunc (self):
        # 1º derivative of the given function
        return self.func.diff()
    
    @property
    def ddfunc (self):
        # 2º derivative of the given function
        return self.dfunc.diff()
    
    @staticmethod
    def _d (function: sym.core.add.Add):
        # Derivates a given function
        return function.diff()
    
    def secant_factor (self, x0: 'Decimal', x1: 'Decimal'):

        return (self.f(x1) - self.f(x0))/ (x1 - x0)
        

    def verify_derivate (self, a):
        # Verify if the derivate of the function in x=a is 0
        if self.dfunc == 0:
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
            (decimal(log(b - a)) - decimal(log(TOL)))/decimal(log(2)))

        print(f'Num. Maximo de iterações: {Max_iterations}')

        if not self.verify_zero(a, b):
            return False
        
        iteration = 0
        x = decimal((a + b)/2)
        iterations_x = []
        iterations_y = []
        iterations_list = [
            {
                'iter': iteration,
                f'a{iteration}': a,
                f'b{iteration}': b,
                f'x{iteration}': x,
                'f(x)': self.f(x),
                'interval': b-a 
            }
        ]
        

        while iteration < Max_iterations:

            x = decimal((a + b)/2)

            if self.f(x) == 0 or decimal((b-a)/2) < TOL:
                if self.f(x)*self.f(a) < 0:
                    b = x
                else:
                    a = x

                
                iteration += 1
                iterations_x.append(x)
                iterations_y.append(self.f(x))
                iterations_list.append(
                    {
                        'iter': iteration,
                        f'a{iteration}': a,
                        f'b{iteration}': b,
                        f'x{iteration}': x,
                        'f(x)': self.f(x),
                        'interval': b - a
                    }
                )
                break
            
            if self.f(x)*self.f(a) < 0:
                b = x
            else:
                a = x

            
            iteration += 1
            iterations_x.append(x)
            iterations_y.append(self.f(x))
            iterations_list.append(
                {
                    'iter': iteration,
                    f'a{iteration}': a,
                    f'b{iteration}': b,
                    f'x{iteration}': x,
                    'f(x)': self.f(x),
                    'interval': b - a
                }
            )

            

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
        
        iteration = 0
        x = a - self.f(a)/self.f(a, self.dfunc)
        
        
        iter_x = []
        iter_y = []
        iterations = [
            {
                'iter': iteration,
                f'x{iteration}': x,
                'f(x)': self.f(x),
                #f'(x)": self.f(a, self.dfunc)
            }
        ]
        
        while abs(self.f(x)) > TOL and iteration < 30:
            iteration += 1

            if not self.verify_derivate(x):
                print(x)
                return False

            if self.f(x) == 0:
                break

            x -= self.f(x)/self.f(x, self.dfunc)

            iter_x.append(x)
            iter_y.append(self.f(x))
            iterations.append(
                {
                    'iter': iteration,
                    f'x{iteration}': x,
                    'f(x)': self.f(x),
                    "f'(x)": self.f(x, self.dfunc)
                }
            )

            

        return {
            'result': x,
            'iter_x': iter_x,
            'iter_y': iter_y,
            'iterations': iterations
        }
    
    def secant_method (self, interval:list, TOL: float):
        """
        Uses the Newton's method to find the solution of the function 
        with TOL precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """

        x0, x1, TOL = decimal(interval[0]), decimal(interval[1]), decimal(TOL)

        if not self.secant_factor(x0, x1):
            return False
        
        iteration = 1
        x2 = x1 - self.f(x1)/self.secant_factor(x0, x1)
        iter_x = []
        iter_y = []
        iterations = [
            {
                'iter': iteration,
                f'x{iteration}': x2,
                'f(x)': self.f(x2),
                "f'(x)": self.secant_factor(x0, x1),
            },
        ]
        print(abs(self.f(x2)))
        
        while abs(self.f(x2)) > TOL and iteration < 15:
            iteration += 1
            if not self.secant_factor(x0, x1):
                return False

            if self.f(x2) == 0:
                break

            x0 = x1

            x1 = x2

            x2 = x1 - self.f(x1)/self.secant_factor(x0, x1)

            iter_x.append(x2)
            iter_y.append(self.f(x2))
            iterations.append(
                {
                    'iter': iteration,
                    f'x{iteration}': x2,
                    'f(x)': self.f(x2),
                    "f'(x)": self.secant_factor(x0, x1)
                }
            )

            

        return {
            'result': x2,
            'iter_x': iter_x,
            'iter_y': iter_y,
            'iterations': iterations
        }
        
    def fake_position_method (self, a: float, b: float, TOL: float) -> dict:
        """
        Uses the bissection method to find the solution of the function
        with TOL's precision.
        Returns a dictionary with the result and lists with the x and y
        coordinates used in the method.
        """
        a, b, TOL = decimal(a), decimal(b), decimal(TOL)

        if not self.verify_zero(a, b):
            return False
        
        iteration = 0
        x = decimal((self.f(b)*a - self.f(a)*b)/(self.f(b) - self.f(a)))
        iterations_x = []
        iterations_y = []
        iterations_list = [
            {
                'iter': iteration,
                f'a{iteration}': a,
                f'b{iteration}': b,
                f'x{iteration}': x,
                'f(x)': self.f(x),
                'interval': b-a 
            }
        ]
        

        while abs(self.f(x)) > TOL and iteration < 15:


            if self.f(x) == 0 or abs(self.f(x)) < TOL:
                if self.f(x)*self.f(a) < 0:
                    b = x
                else:
                    a = x

                
                iteration += 1
                iterations_x.append(x)
                iterations_y.append(self.f(x))
                iterations_list.append(
                    {
                        'iter': iteration,
                        f'a{iteration}': a,
                        f'b{iteration}': b,
                        f'x{iteration}': x,
                        'f(x)': self.f(x),
                        'interval': b - a
                    }
                )
                break
            
            if self.f(x)*self.f(a) < 0:
                b = x
            else:
                a = x

            
            iteration += 1
            x = decimal((self.f(b)*a - self.f(a)*b)/(self.f(b) - self.f(a)))
            iterations_x.append(x)
            iterations_y.append(self.f(x))
            iterations_list.append(
                {
                    'iter': iteration,
                    f'a{iteration}': a,
                    f'b{iteration}': b,
                    f'x{iteration}': x,
                    'f(x)': self.f(x),
                    'interval': b - a
                }
            )

            

        return {
            'result': x,
            'iter_x': iterations_x,
            'iter_y': iterations_y,
            'iterations': iterations_list
        }
        




if __name__ == '__main__':
    # tests
    x = sym.symbols('x')

    f = x**3 + 4*x**2 - 2

    func = Func(
        f, x
    )
    print(
        *func.secant_method([-2, 0], 10**(-1))['iterations'], sep='\n'
    )
    ...




