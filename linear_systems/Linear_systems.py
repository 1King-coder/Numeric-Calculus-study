from Vectors import Vector
from copy import deepcopy
from Matrix import Matrix
import sympy as sym

class Linear_System (Matrix):

    def __init__ (self, linear_sys: list = []) -> None:
        super().__init__ (linear_sys)

        
        self.num_of_changed_lines = 0
        self.icognitos = [sym.symbols(f'x{i+1}') for i in range(self.rows_num)]
        self.system = linear_sys

    @property
    def system (self) -> list:
        return self.__system
    
    @system.setter
    def system (self, value) -> None:
        if isinstance(value[0], sym.core.add.Add):
            ...
        self.__system = self.with_icognitos(value)

    def with_icognitos (self, system):
        sys = deepcopy(system)

        for j in range(self.rows_num):
            for index, icognito in enumerate(self.icognitos):
                sys[j][index] *= icognito

            sys[j] = sum(sys[j])

        return sys

    def __repr__ (self) -> str:
        lines_str = [f"l{index + 1} | {self.__system[index]} = 0" for index in range(self.rows_num)]
        string = "\n".join(lines_str)

        return string + "\n"
    
    def show_LU_decompotion (self) -> None:
        """
        Prints the LU decompotion of the system.
        """
        for key in self.__dict__.keys():
            if 'factor' in key:

                lines_str = [
                    f"\tl{index + 1} | {tuple(self.__dict__[key][index])}"
                    for index in range(self.rows_num)
                ]

                string = "\n".join(lines_str)
                print(f"{key[9:]}:\n{string}")
    
    def solve_by_gauss_elimination (self) -> Vector:
        """
        Function that utilizes the gaussian elimination method
        to triangulate the system and solve it
        """

        # Takes the partialy escalonated by gaussian elimination 
        # (U factor of the system) and reverses it to solve it starting from the last
        # variable.
        U_factor_with_icognitos = self.with_icognitos(self.U_factor)[::-1]

        # Reversed list of the icognitos
        icognitos_to_solve = deepcopy(self.icognitos[::-1])

        # dictionary to make it possible and easier to use .subs from sympy
        # and store the results.
        solved_icognitos = {
            icognito: icognito for icognito in icognitos_to_solve
        }

        for index, line in enumerate(U_factor_with_icognitos):

            # substitute the already solved variables values in the
            # line equation, making it only 1 icognito.
            line = line.subs(solved_icognitos)

            # takes the result by using sym.solve, where it takes the
            # line equation and solve it for equation = 0
            res = sym.solve(line, icognitos_to_solve[index], rational=False)[0]

            # Round and store the results
            print (res)
            solved_icognitos[icognitos_to_solve[index]] = res

        # Sorts the solutions
        solved_icognitos = {
            icognito: solution for icognito, solution in list(solved_icognitos.items())[::-1]
        }

        return Vector(*solved_icognitos.values())
    
    def fi (self):
        def fi_element (row: int, col: int):
            self.matrix[row][col]/self.matrix[row][row]
        
        fi_matrix = self.map_matrix(

        )
        ...

    def fi_line (self):
        ...

    
        
if __name__ == '__main__':

    sys_1 = [
        [2, 3, 1, -4, 0],
        [-3, 1, 1, 2, -4],
        [-1, 0, 1, -1, -2],
        [1, 1, -3, 0, 2],
    ]

    sys_2 = [
        [2, 1, 1, 0],
        [4, 3, 3, 1],
        [8, 7, 9, 5],
        [6, 7, 9, 8],
    ]

    sys_3 = [
        [2, -1, 1, 0, -1],
        [0, 3, 1, 2, -1],
        [-1, 0, 3, -1, 5],
        [1, 0, -3, 2, -6]
    ]


    linear_system = Linear_System(sys_3)
    linear_system.fi()
    



