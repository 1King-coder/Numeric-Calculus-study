import sympy as sym
import numpy as np

class SolvingODEs:

    def __init__(self, equation: 'sym.Add', points: list, precision: int = 5) -> None:
        self.equation: 'sym.Add' = equation
        self.points = points
        self.ode_type = "IVP" if len(points) == 1 else "BVP"
        self.prec = precision
        self.ode_order = equation.args[-1].count("'")

    def euler_equation (self, f: 'sym.Add', last_point: tuple, step: float) -> tuple:
        """
        This method calculates the next point of the solution using Euler's method.
        
        Parameters
        ----------
        f : sym.Sum
            The function of the equation.
        last_point : tuple
            The last point of the solution.
        step : float
            The step of the solution.

        Returns
        -------
        tuple
            The next point of the solution.
        """
        
        return (last_point[1], last_point[1] + step * round(f.subs({
            sym.Symbol("x"): last_point[0],
            sym.Symbol("y"): last_point[1]
        }).evalf(), self.prec))[1]

    def euler_method(self, interval: list, step: float) -> list:
        """
        This method implements Euler's method for numerical approximation of 
        solutions to initial value problems of ordinary differential equations (ODEs).

        base equation:

            y = y_i-1 + h * f(x_i-1, y_i-1)

            where y' = f(x, y)

        Parameters
        ----------
        interval : list
            The interval of the solution.
        step : float
            The step of the solution.

        Returns
        -------
        list
            A list with the solution of the ODE.
        """

        if self.ode_type != "IVP":
            raise TypeError("This method only works for initial value problems.")
        
        if self.points[0][0] != interval[0]:
            raise ValueError("The initial point is not in the interval.")
         
        f = sym.solve(self.equation, sym.Symbol("y'"))[0]
        
        
        solutions = [
            (
                round(self.points[0][0] + step, self.prec),
                self.euler_equation(f, self.points[0], step)
            )
        ]

        x_values_in_interval = np.arange(
            interval[0] + 2 * step,
            interval[1] + step,
            step
        )

        for i, x_val in enumerate(x_values_in_interval):
            solutions.append(
                (
                    round(x_val, self.prec),
                    self.euler_equation(f, solutions[i], step)
                )
            )

        return solutions
    
    
    def taylor_series_method (self, interval: list, step: float) -> list:
        if self.ode_type != "IVP":
            raise TypeError("This method only works for initial value problems.")
        
        if self.points[0][0] != interval[0]:
            raise ValueError("The initial point is not in the interval.")
        
        derivatives_funcs = {}

        derivatives_funcs[sym.Symbol("y'")] = sym.solve(
            self.equation,
            sym.Symbol("y'")
        )[0]
        



if __name__ == "__main__":
    x, y, y_1, y_2 = sym.symbols("x y y' y''")

    eq: 'sym.Add' = x + y_2 * 9 - y + 2 - 5*y_1

    print(str(eq.args[-1]).count("'"))

