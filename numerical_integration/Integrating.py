import sys
is_main = __name__ == "__main__"

function_path = '.\\zeros_of_functions\\' if is_main else '..\\zeros_of_functions\\'
lin_sys_path = '.\\linear_systems\\' if is_main else '..\\linear_systems\\'
interpolation_path = '.\\interpolation\\' if is_main else '..\\interpolation\\'


sys.path.append(interpolation_path)
sys.path.append(lin_sys_path)
sys.path.append(function_path)

from points_Interpolation import Interpolate
from Functions import Func
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from copy import deepcopy
import math


class Integrate (Interpolate):

  def __init__ (self, func: 'Func', interval: list, num_of_subintervals: int, step: float = None, precision: int = 5) -> None:
    self.func = func
    self.interval = interval
    self.num_of_subintervals = num_of_subintervals
    self.step = step

    self.points = [
      (x_coord, y_coord)
      for x_coord, y_coord in zip(self.x_values, self.y_values)
    ]

    super().__init__(self.points, self.func, precision)

  @property
  def step (self) -> float:
    return self.__step

  @step.setter
  def step (self, value) -> None:
    """
    If step is not directly given, uses the number of subintervals wanted
    to calculate the step.
    """
    if value:
      self.__step = value
      return
    
    value = (self.interval[-1] - self.interval[0]) / self.num_of_subintervals
    self.__step = value
      


  @property
  def x_values (self) -> np.ndarray:
    """
    Creates the x_coordinates spaced by the step wanted.
    stop = last interval + step because np.arange always stops 1 step before.
    """
    return np.arange(self.interval[0], self.interval[1] + self.step, step=self.step)

  @property
  def y_values (self) -> list:
    # Create the y_coordinates.
    return [self.func(x_coord) for x_coord in self.x_values]

  
  def paralelogram_method (self) -> float:
    """
    This method is used to calculate a defined integral by replacing the original
    function by a interpolation polynomial with 1 degree (calculated by lagrange method).
    The final formule is: (Err(x) is the error from corollary 2 for interpolation polynoms)
    Integral(f(x), from a to b) = Σ(h/2)*(f(xi) - f(xi-1)) + Integral(Err(x), same interval)
                                = (h/2) * (2*(Σfxi from x1 to xi-1) + f(x0) + f(xi)) + Integral(Err(x), same interval)
    """
    integral =  (
      (self.step / 2) * 
      # Here I manipulated the equation to have 2*(f(x0) +f(xi)) - f(x0) - f(x0) = f(x0) + f(xi)
      (2 * sum(self.y_values) - (self.y_values[0] + self.y_values[-1]))
    )

    return round(integral, self.prec)

  def paralelogram_method_error (self) -> float:
    """
    The error estimated is the integral of corollary 2 formula between the interval.
    This integral results in:
    Err(x) <=  (b - a) * (h^(n+1) . M[n+1]) / 12
    So: Err(x) <= [(h^(n+1) . M[n+1]) / 4 * (n+1)] * (2/3) * (b - a) 
    """

    return round(
      self.corollary_2(1) * 2 * (self.interval[1] - self.interval[0]) / 3,
      self.prec
    )

  def one_third_simpson_method (self) -> float:
    """
    This method is used to calculate a defined integral by replacing the original
    function by a interpolation polynomial with 2 degree (calculated by lagrange method).
    The final formule is: (Err(x) is the error from corollary 2 for interpolation polynoms but for 3 degree polynomial)
    Integral(f(x), from a to b) = [Σ(h/3)*(f(x2i-2) + 4*f(x2i-1) + f(x2i))] + Integral(Err(x), same interval)
    """

    # In this method, f(x2i-1) will be multiplied by 4
    # and f(x2i) by 2 (except by x0 and xn), so even i's is doubled and odds is 
    # quadrupled

    images_multiplied = [
      2 * self.y_values[i] # double if i is even
      if i % 2 == 0 else
      4 * self.y_values[i] # quadruples if i is odd
      for i in range(1, len(self.y_values) - 1)
    ]
    
    integral = (
      (self.step / 3) * 
      (self.y_values[0] + self.y_values[-1] + 
       sum(images_multiplied))
    )

    return round(integral, self.prec)

  def one_third_simpson_method_error (self) -> float:
    """
    The error estimated is the integral of corollary 2 formula between the interval.
    This integral results in:
    Err(x) <=  (b - a) * (h^(n+1) . M[n+1]) / 12
    So: Err(x) <= [(h^(n+1) . M[n+1]) / 4 * (n+1)] * (2/3) * (b - a) 
    """

    return round(
      self.corollary_2(3) * 4 * (self.interval[1] - self.interval[0]) / 45,
      self.prec
    )


if is_main:
  x = sym.symbols('x')
  intergrate = Integrate(Func(sym.sin(x),x), [1, 3], num_of_subintervals=10, precision=7)

  print(intergrate.one_third_simpson_method())
  print(intergrate.one_third_simpson_method_error())
  print(intergrate.derivatives)
    