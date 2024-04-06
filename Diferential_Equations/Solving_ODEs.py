import sympy as sym
import numpy as np

class SolvingODEs:

    def __init__(self, equation: 'sym.Add', ode_order: int, points: dict, precision: int = 5) -> None:
        x = sym.Symbol("x")
        y = sym.Function("y")(x)

        self.equation: 'sym.Add' = equation
        self.points = points
        self.ode_type = None
        self.ode_order = ode_order
        self.prec = precision

    @property
    def ode_type (self) -> str:
        return self._ode_type
        
    @ode_type.setter
    def ode_type (self, ode_type: str):
        ode_type = "IVP"

        for val in self.points.values():
            if len(val) > 1:
                ode_type = "BVP"
                return
            
        self._ode_type = ode_type

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
            x: last_point[0],
            y: last_point[1]
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
        
        x = sym.Symbol("x")
        y = sym.Function("y")(x)



        f = sym.solve(self.equation, y.diff(x))[0]

        interval = [self.points[0][0], target]
        
        solutions = [
            (
                round(self.points[y][0][0] + step, self.prec),
                self.euler_equation(f, self.points[y][0], step)
            )
        ]

        x_values_in_interval = np.arange(
            self.points[y][0][0] + 2 * step,
            self.points[y][0][1] + step,
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
        y = sym.Function("y")(x)
        taylor_expansion = 0
        parcels = [
            y
        ]
        
        for i in range(0, taylor_order + 1):
            taylor_expansion += parcels[-1] * (step**i) / sym.factorial(i)

            parcels.append(
                parcels[-1].diff(x).subs(
                    {
                        y.diff(x): derivative
                    }
                )
            )
            
        return taylor_expansion
    
    def taylor_series_method (self, target, taylor_order: int, step: float) -> list:
        
        x = sym.Symbol("x")
        y = sym.Function("y")(x)

        if self.ode_type != "IVP":
            raise TypeError("This method only works for initial value problems.")
        
        if self.points[y][0][0] >= target:
            raise ValueError("We can only approximate the solution for points beyond the initial point.")
        
        if self.ode_order > 1:
            raise NotImplementedError("This method only works for first order ODEs.")
    

        derivative: 'sym.Add' = sym.solve(
            self.equation,
            y.diff(x)
        )[0]

        taylor_expression = self.taylor_series_formulas(derivative, taylor_order, step)

        solutions = [
            (
                round(self.points[y][0][0] + step, self.prec),
                round(taylor_expression.subs(
                    {
                        x: self.points[y][0][0],
                        y: self.points[y][0][1]
                    }
                ), self.prec)
            )
        ]

        x_values_in_interval = np.arange(
            self.points[y][0][0] + 2 * step,
            target + step,
            step
        )

        for i, x_val in enumerate(x_values_in_interval):
            print (solutions)
            solutions.append(
                (
                    x_val,
                    round(taylor_expression.subs(
                    {
                        x: solutions[i][0],
                        y: solutions[i-1][1]
                    }
                ), self.prec)
                )
            )

        return solutions
    
    def runge_kutta_increment_func (self, f: 'sym.Add', rk_order: int, point: list, constants: dict, step: float) -> float:
        x = sym.Symbol("x")
        y = sym.Function("y")(x)
        x_ = point[0]
        y_ = point[1]
        m_vals = []
        full_expression = y_

        for i in range(rk_order):
            m = f.subs({ 
                x: x_,
                y: y_
            }).evalf()

            m_vals.append(m)

            x_ = x_ + step * constants["b"][i]
            y_ = y_ + sum([step * constants["a"][j + 1] * m_vals[j] for j in range(i)])
            parcels = [constants["a"][k] * m_vals[k] * step for k in range(rk_order - 1)]
            print ( parcels)
        full_expression += parcels


        return full_expression


        
        ...
        
    def runge_kutta_method (self, target: float, rk_order: int, step: float, constants: dict) -> list:
        if self.ode_type != "IVP":
            raise TypeError("This method only works for initial value problems.")
        
        x = sym.Symbol("x")
        y = sym.Function("y")(x)

        bigger_derivative = list(self.points.keys())[-1]

        f = sym.solve(self.equation, bigger_derivative)[0]

        solutions = [
           (round(self.points[y][0][0]) + step,
            self.runge_kutta_increment_func(f, rk_order, self.points[y][0], constants, step).subs(
                {
                    x: self.points[y][0][0],
                    y: self.points[y][0][1]
                }
            ))
        ]

        for i, x_val in enumerate(np.arange(self.points[y][0][0] + 2*step, target + step, step)):
            solutions.append(
                (
                    x_val,
                    self.runge_kutta_increment_func(f, rk_order, solutions[i], constants, step)
                )
            )

        return solutions    



if __name__ == "__main__":
    x = sym.symbols("x")
    y = sym.Function('y')(x)

    eq: 'sym.Add' = y + y.diff(x) - 2 * x**3



    ode = SolvingODEs(eq, 1, {y: [(1, 0)]}, 7)

    constants = {
        "a": [0.25, 0.75],
        "b": [2/3, 2/3]
    }

    print(ode.runge_kutta_method(1.2, 2, 0.1, constants))

