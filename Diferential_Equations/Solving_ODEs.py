import sympy as sym
import numpy as np

class SolvingODEs:

    def __init__(self, equation: 'sym.Add', points: list, precision: int = 5) -> None:
        self.equation: 'sym.Add' = equation
        self.points = points
        self.ode_type = "IVP" if len(points) == 1 else "BVP"
        self.prec = precision
        self.ode_order = str(equation.args[-1]).count("'")
        

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
        
        return last_point[1] + step * round(f.subs({
            sym.Symbol("x"): last_point[0],
            sym.Function("y")(sym.Symbol("x")): last_point[1]
        }).evalf(), self.prec)

    def euler_method(self, target: float, step: float) -> list:
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
        
        if self.points[0][0] >= target:
            raise ValueError("We can only approximate the solution for points beyond the initial point.")
         
        f = sym.solve(self.equation, sym.Function("y")(sym.Symbol("x")).diff(sym.Symbol("x")))[0]

        interval = [self.points[0][0], target]
        
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
    
    def taylor_series_formulas (self, derivative: 'sym.Add', taylor_order: int, step: float) -> list:

        x = sym.Symbol("x")
        taylor_expansion = sym.Function("y")(x)
        parcels = [
            derivative
        ]
        
        for i in range(1, taylor_order + 1):
            taylor_expansion += parcels[-1] * (step**i) / sym.factorial(i)

            parcels.append(
                parcels[-1].diff(x).subs(
                    {
                        sym.Function("y")(x).diff(x): derivative
                    }
                )
            )
            
        return taylor_expansion

        ...
    
    def taylor_series_method (self, target, taylor_order: int, step: float) -> list:
        if self.ode_type != "IVP":
            raise TypeError("This method only works for initial value problems.")
        
        if self.points[0][0] >= target:
            raise ValueError("We can only approximate the solution for points beyond the initial point.")
        
        if self.ode_order > 1:
            raise NotImplementedError("This method only works for first order ODEs.")
    

        derivative: 'sym.Add' = sym.solve(
            self.equation,
            sym.Function('y')(sym.Symbol("x")).diff(sym.Symbol("x"))
        )[0]

        

        taylor_expression = self.taylor_series_formulas(derivative, taylor_order, step)
        print(taylor_expression)
        solutions = [
            (
                round(self.points[0][0] + step, self.prec),
                round(taylor_expression.subs(
                    {
                        sym.Symbol("x"): self.points[0][0],
                        sym.Function("y")(sym.Symbol("x")): self.points[0][1]
                    }
                ), self.prec)
            )
        ]
        interval = [self.points[0][0], target]

        x_values_in_interval = np.arange(
            interval[0] + 2 * step,
            interval[1] + step,
            step
        )

        for i, x_val in enumerate(x_values_in_interval):
            solutions.append(
                (
                    x_val,
                    round(taylor_expression.subs(
                    {
                        sym.Symbol("x"): x_val + step,
                        sym.Function("y")(sym.Symbol("x")): solutions[i][1]
                    }
                ), self.prec)
                )
            )

        return solutions
        



if __name__ == "__main__":
    x = sym.symbols("x")
    y = sym.Function('y')(x)

    eq: 'sym.Add' = y + y.diff(x) - 2 * x**3



    ode = SolvingODEs(eq, [(1, 0)])

    print(ode.euler_method(1.2, 0.1))

