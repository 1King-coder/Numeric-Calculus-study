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