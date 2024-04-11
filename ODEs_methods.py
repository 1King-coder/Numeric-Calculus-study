from numeric_calculus_mashup_package.Functions import Func
from numeric_calculus_mashup_package.Matrix import Matrix
from numeric_calculus_mashup_package.Vectors import Vector
from numeric_calculus_mashup_package.Linear_systems import Linear_System
from numeric_calculus_mashup_package.Linear_adjusts import Linear_adjust

import sympy as sym
class Solving_ODEs:

    def __init__(self, equation: sym.Sum, points: list) -> None:
        self.equation = equation
        self.points = points
        self.ode_type = "IVP" if len(points) == 1 else "BVP"

    
    def euler_method (self, interval: list, wanted_x: float):
        
        pass






