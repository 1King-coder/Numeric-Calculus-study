from Vectors import Vector
from copy import deepcopy
from Matrix import Matrix
import sympy as sym

class Linear_System (Matrix):

    def __init__ (self, linear_sys: list = []) -> None:
        super().__init__ (linear_sys)

        
        self.num_of_changed_lines = 0

        self.icognitos = [sym.symbols(f'x{i+1}') for i in range(self.rows_num)]
        if not self.isQuadratic:
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
        U_factor_with_icognitos = [sum(self.with_icognitos(self.U_factor)[i]) for i in range(self.rows_num)][::-1]

        # Reversed list of the icognitos
        icognitos_to_solve = deepcopy(self.icognitos[len(self.icognitos) - 2::-1])

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
            solved_icognitos[icognitos_to_solve[index]] = res

        # Sorts the solutions
        solved_icognitos = {
            icognito: solution for icognito, solution in list(solved_icognitos.items())[::-1]
        }

        return Vector(*solved_icognitos.values())
    
    @property
    def fi_function (self):
        """
        Fi function matrix for Gauss-Jacobi solving linear system method.

        ɸ(x) = [
            x1 (x2, ..., xn),
            x2 (x1,x3, ..., xn),
            x3 (x1, x2, x4, ..., xn),
            ...
        ]
        """
        for i in range(self.rows_num):
            if not self.matrix[i][i]:
                raise InterruptedError("Can not create fi for matrices with 0 in its diagonal")

        def fi_matrix_gen (i, j):
            if not i == j:
                return -self.matrix[i][j] / self.matrix[i][i] * self.icognitos[j]
            
        fi_matrix = self.map_matrix (
            fi_matrix_gen,
            self.rows_num,
            self.cols_num
        )

        fi_matrix = [sum(line) for line in fi_matrix]

        return fi_matrix

    @property
    def fi_diff (self):
        """
        Jacobian matrix of ɸ(x).
        ɸ'(x) = [
            ɸ(x)[0].diff(x1), ɸ(x)[0].diff(x2), ɸ(x)[0].diff(x3), ...
            ɸ(x)[1].diff(x1), ...
            ...
        ]
        """
        fi_diff = self.map_matrix(
            lambda i, j: self.fi_function[i].diff(self.icognitos[j]),
            self.rows_num,
            self.cols_num - 1
        )

        return fi_diff
    
    @property
    def converge_crit_alpha (self) -> dict:
        """
        Verify if the fi function obtained has the guarantee to converge to a solution
        when using gauss-jacobi method.
        """
        alphas_dict = {}

        for i in range(len(self.fi_diff)):
            alpha = 0
            for j in range(len(self.fi_diff[0])):
                alpha += abs(self.fi_diff[i][j])

            alphas_dict[f'Alpha({i+1})'] = alpha
                
        return alphas_dict
    
    @property
    def Sassenfeld_crit (self) -> dict:
        """
        Convergence criterea for gauss-Sciedel linear system solving method.
        if Beta < 1 the results will converge to a solution.
        """
        betas_dict = {}

        beta = 1

        for i in range(self.rows_num):
            beta_res = 0
            for j in range(self.cols_num - 1):
                if i == j:
                    continue

                beta = betas_dict.get(f'Beta({j + 1})')

                if not beta:
                    beta = 1


                beta_res += abs(self.matrix[i][j]) * beta

            betas_dict[f'Beta({i + 1})'] = sym.Rational(beta_res, abs(self.matrix[i][i]))

        return betas_dict
    
    def gauss_scibel_method (self, TOL) -> list:
        iterations = []

        has_no_coeff = False

        for i in range(self.rows_num):
            for j in range(self.cols_num - 1):
                if not self.matrix[i][j]:
                    has_no_coeff = True
        
        icognitos_values_dict = {
            self.icognitos[i]: 0 for i in range(self.rows_num)
        }
        
        icognitos_values_dict_0 = deepcopy(icognitos_values_dict)

        sol_vector_0 = Vector(*[0 for _ in range(self.rows_num)])

        for i in range(self.rows_num):
            icognitos_values = icognitos_values_dict
            
            if has_no_coeff:
                icognitos_values = icognitos_values_dict_0

            icognitos_values_dict[self.icognitos[i]] = round(self.fi_function[i].subs(icognitos_values), 5)


        sol_vector_1 = Vector(*icognitos_values_dict.values())

        iterations = [{'iter': 0, 'x(0)': sol_vector_0}, {'iter': 1, 'x(1)': sol_vector_1}]
        
        iteration = 1

        dif_module = (sol_vector_1 - sol_vector_0).module

        while not dif_module < TOL:
            
            iteration += 1
            sol_vector_0 = sol_vector_1

            icognitos_values_dict_0 = deepcopy(icognitos_values_dict)

            for i in range(self.rows_num):
                icognitos_values = icognitos_values_dict
                if has_no_coeff:
                    icognitos_values = icognitos_values_dict_0

                icognitos_values_dict[self.icognitos[i]] = round(self.fi_function[i].subs(icognitos_values), 5)

            sol_vector_1 = Vector(*icognitos_values_dict.values())

            dif_module = (sol_vector_1 - sol_vector_0).module

            iterations.append({
                'iter': iteration,
                f'x({iteration})': sol_vector_1,
                f'||x({iteration}) - x({iteration-1})||':  round(dif_module, 5)
            })

        return iterations




    
        
if __name__ == '__main__':

    sys_1 = [
        [2, 3, 1, -4, 0],
        [-3, 1, 1, 2, -4],
        [-1, 0, 1, -1, -2],
        [1, 1, -3, 0, 2],
    ]

    sys_2 = [
        [2, -1, 1, 0, -1],
        [0, 3, 1, 2, -1],
        [-1, 0, 3, -1, 5],
        [1, 0, -3, 2, -6],
    ]

    sys_3 = [
        [4, -1, 1, -1, -3],
        [2, 4, -1, 1, -3],
        [-1, 2, 3, -1, 0],
        [1, -1, -2, -3, 2],
    ]

    
    linear_system = Linear_System(sys_2)
    """
    fi_func = linear_system.fi_function
    fi_diff = linear_system.fi_diff
    alpha_dict = linear_system.converge_crit_alpha

    print(*fi_func, sep='\n')
    print(Matrix(fi_diff))
    print(linear_system.Sassenfeld_crit)
    """
    print(*linear_system.gauss_scibel_method(10**(-1)), sep='\n')

    linear_system.show_LU_decompotion()


    



