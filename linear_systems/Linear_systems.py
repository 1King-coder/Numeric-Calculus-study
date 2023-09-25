from Vectors import Vector
from copy import deepcopy
from Matrix import Matrix
import sympy as sym

class Linear_System (Matrix):

    def __init__ (self, linear_sys: list = []) -> None:
        super().__init__ (linear_sys)

        
        self.num_of_changed_lines = 0

        self.icognitos = [sym.symbols(f'x{i+1}') for i in range(self.rows_num)]
        self.icognitos.append(1) # constant

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

            #sys[j] = sum(sys[j])

        return sys

    def __repr__ (self) -> str:
        lines_str = [f"l{index + 1} | {sum(self.__system[index])} = 0" for index in range(self.rows_num)]
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
    
    def fi_function (self, system_matrix: 'Matrix', icognitos: list):
        """
        Fi function matrix for Gauss-Jacobi solving linear system method.

        ɸ(x) = [
            x1 (x2, ..., xn),
            x2 (x1,x3, ..., xn),
            x3 (x1, x2, x4, ..., xn),
            ...
        ]
        """
        for i in range(system_matrix.rows_num):
            if not system_matrix.matrix[i][i]:
                raise InterruptedError("Can not create fi for matrices with 0 in its diagonal")

        def fi_matrix_gen (i, j):
            if not i == j:
                return -sym.Rational(system_matrix.matrix[i][j], system_matrix.matrix[i][i]) * icognitos[j]
            
        fi_matrix = self.map_matrix (
            fi_matrix_gen,
            system_matrix.rows_num,
            system_matrix.cols_num
        )

        fi_matrix = [sum(line) for line in fi_matrix]

        return fi_matrix

    def fi_diff (self, fi_func: list, icognitos: list):
        """
        Jacobian matrix of ɸ(x).
        ɸ'(x) = [
            ɸ(x)[0].diff(x1), ɸ(x)[0].diff(x2), ɸ(x)[0].diff(x3), ...
            ɸ(x)[1].diff(x1), ...
            ...
        ]
        """
        fi_diff = self.map_matrix(
            lambda i, j: fi_func[i].diff(icognitos[j]),
            self.rows_num,
            self.cols_num - 1
        )

        return fi_diff
    
    def guarantee_fi_converge (self) -> bool:
        """
        Verify if the fi function obtained has the guarantee to converge to a solution
        when using gauss-jacobi method.
        """
        for i in range(len(self.fi_diff)):
            alpha = 0
            for j in range(len(self.fi_diff[0])):
                aplha += abs(self.fi_diff[i][j])

                if alpha < 1:
                    return False
                
        return True
    
    def gauss_scibel_method (self):
        
        return




    
        
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
        [2, 1, -1, -4],
        [-1, 1, 1, 0],
        [1, 3, -1, -8],
    ]


    linear_system = Linear_System(sys_3)
    fi_func = linear_system.fi_function(Matrix(linear_system.matrix), linear_system.icognitos)
    print(linear_system.fi_diff(fi_func, linear_system.icognitos))
    print(fi_func)


    



