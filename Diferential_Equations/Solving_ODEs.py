import sympy as sym
import numpy as np

class SolvingODEs:

    def __init__(self, equation: sym.Sum, points: list) -> None:
        self.equation: sym.Sum = equation
        self.points = points
        self.ode_type = "IVP" if len(points) == 1 else "BVP"

    def euler_method(self, interval: list, step: float) -> list:
        f = sym.solve(self.equation, sym.Symbol("y'"))[0]

        solutions = [
            self.points[0]
        ]
        print(solutions[0][0], solutions[0][1])
        for i, x_val in enumerate(np.arange(interval[0], interval[1], step)):
            solutions.append(
                solutions[i-1][1] + step * f.subs({"x": solutions[i-1][0], "y": solutions[i-1][1]}) 
            )

        return solutions


if __name__ == "__main__":
    x, y, y_1 = sym.symbols("x y y'")

    eq = x - y + 2 - y_1

    ode = SolvingODEs(eq, [(0, 2)])

    print(ode.euler_method([0, 1], 0.1))

