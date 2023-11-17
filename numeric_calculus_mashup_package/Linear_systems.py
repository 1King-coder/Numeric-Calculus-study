from .Matrix import Matrix
from .Vectors import Vector    
from copy import deepcopy
from decimal import Decimal
import sympy as sym

class Linear_System (Matrix):

    def __init__ (self, transform_matrix: list, b_vector: 'Vector', precision: int = 5, custom_icognito: str = 'x') -> None:
        
        self.rows_num = len(transform_matrix)
        self.cols_num = len(transform_matrix[0])
        self.custom_icognito = custom_icognito
        self.icognitos = [sym.symbols(f'{custom_icognito}{i+1}') for i in range(self.rows_num)]
        
        if self.rows_num != self.cols_num:
            transform_matrix = self.expression_to_matrix(transform_matrix)

        super().__init__ (transform_matrix)

        self.b_vector = b_vector
        self.precision = precision
        self.ampli_matrix = self.amplified_matrix()
        self.num_of_changed_lines = 0

        self.transform_matrix_with_icognitos = transform_matrix

    @property
    def transform_matrix_with_icognitos (self) -> list:
        return self.__transform_matrix_with_icognitos
    
    @transform_matrix_with_icognitos.setter
    def transform_matrix_with_icognitos (self, value) -> None:
        # if isinstance(value[0], sym.core.add.Add):
            # ...
        self.__transform_matrix_with_icognitos = self.with_icognitos(value)

    def expression_to_matrix (self, transform_matrix):
        
        return self.map_matrix(
            lambda i, j: transform_matrix[i][0].coeff(self.icognitos[j]),
            self.rows_num,
            self.rows_num,
        )

    def with_icognitos (self, transform_matrix):
        sys = deepcopy(transform_matrix)

        for j in range(self.rows_num):
            for index, icognito in enumerate(self.icognitos):
                sys[j][index] *= icognito

        return sys

    def __repr__ (self) -> str:
        lines_str = [f"l{index + 1} | {sum(self.transform_matrix_with_icognitos[index])} = {self.b_vector[index]}" for index in range(self.rows_num)]
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

    def amplified_matrix (self) -> list:
        return self.map_matrix (
            lambda i, j: self.matrix[i][j] if j != self.cols_num else -1 * self.b_vector[i],
            self.rows_num,
            self.cols_num + 1
        )
    
    def gauss_elimination_method (self, is_using_LU: bool = False) -> Vector:
        """
        Function that utilizes the gaussian elimination method
        to triangulate the system and solve it
        """

        upper_triangular_matrix = self.gaussian_elimination(self.ampli_matrix)
        
        # utm stands for upper triangular matrix
        utm_factor_with_icognitos = [
            sum(self.with_icognitos(upper_triangular_matrix)[i])
            for i in range(self.rows_num)
        ][::-1]
        
        # Reversed list of the icognitos
        icognitos_to_solve = deepcopy(self.icognitos)[::-1]

        if is_using_LU:
            utm_factor_with_icognitos.reverse()
            icognitos_to_solve.reverse()

        # dictionary to make it possible and easier to use .subs from sympy
        # and store the results.
        solved_icognitos = {
            icognito: icognito for icognito in icognitos_to_solve
        }

        for index, line in enumerate(utm_factor_with_icognitos):
            # substitute the already solved variables values in the
            # line equation, making it only 1 icognito.
            line = line.subs(solved_icognitos)

            # takes the result by using sym.solve, where it takes the
            # line equation and solve it for equation = 0
            res = sym.solve(line, icognitos_to_solve[index])[0]

            # Round and store the results
            solved_icognitos[icognitos_to_solve[index]] = round(res,  self.precision)

        # Sorts the solutions
        solved_icog_list = list(solved_icognitos.items())

        if not is_using_LU: 
            solved_icog_list.reverse()

        solved_icognitos = {
            icognito: solution 
            for icognito, solution in 
            solved_icog_list
        }

        return Vector(*solved_icognitos.values(), icognito=self.custom_icognito)
    
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
                term = -(round(self.matrix[i][j] / self.matrix[i][i], self.precision)) * self.icognitos[j]
                return term
            
        fi_matrix = self.map_matrix (
            fi_matrix_gen,
            self.rows_num,
            self.cols_num
        )

        b_vec_value = lambda i: round(self.b_vector[i] / self.matrix[i][i], 5)

        fi_matrix = [sum(line) - b_vec_value(i) for i, line in enumerate(fi_matrix)]

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
            lambda i, j: round(self.fi_function[i].diff(self.icognitos[j]), self.precision),
            self.rows_num,
            self.cols_num
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
    def sassenfeld_crit (self) -> dict:
        """
        Convergence criterea for gauss-Sciedel linear system solving method.
        if Beta < 1 for all i, the results will converge to a solution.
        βi = [Sum(|aij|*|βj|, from(j=1), to (i-1)) + sum(|aij|, from(i+1), to n)]/|aii|
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

                

                beta_res += abs(self.matrix[i][j] * beta)

            

            betas_dict[f'Beta({i + 1})'] = beta_res/abs(self.matrix[i][i])

        return betas_dict
    
    def verify_if_have_coeff (self) -> bool:
        """
        Verify if there is an as coeff of any icognito.
        """
        has_no_coeff = False

        for line in self.matrix:

            if 0 in line:
                has_no_coeff = True
                break

        return has_no_coeff
    
    def gauss_sciebel_will_converge (self) -> bool:
        """
        Verify if all betas is < 1, if one of them is >= 1, it is not guaranteed 
        that the method will converge.
        """
        sassenfeld_crits = [beta < 1 for beta in self.sassenfeld_crit.values()]

        will_converge = False in sassenfeld_crits

        return not will_converge

    def gauss_scibel_method (self, initial_vector: 'Vector', TOL) -> list:
        iterations = []

        there_is_zero = self.verify_if_have_coeff()

        actual_icognitos_values = {
            self.icognitos[i]: initial_vector[i] for i in range(self.rows_num)
        }
        # saves values as 0 if there is an 0 in the transform matrix
        last_icognitos_values = deepcopy(actual_icognitos_values)

        

        def calculate_icognitos_values (actual_icognitos_values: dict,
                                        last_icognitos_values: dict,
                                        there_is_zero: bool = there_is_zero) -> dict:

            """
            Function to update the icognitos values.
            """

            act_icog_vals = deepcopy(actual_icognitos_values)
            last_icog_vals = deepcopy(last_icognitos_values)


            # We can not use the last calculated value for a icognito if there is 
            # one that has 0 as coeff
            to_use_values = act_icog_vals if not there_is_zero else last_icog_vals
            
            solution_values = {}
            for i in range(self.rows_num):
                result = round(self.fi_function[i].subs(to_use_values), self.precision)
                
                # real-time updating icognitos values if there is no zero in the coefficients
                act_icog_vals[self.icognitos[i]] = result
                solution_values[self.icognitos[i]] = result

            return solution_values
        
        max_iterations = 9999999999

        if not self.gauss_sciebel_will_converge():
            print ("This method is not guaranteed to converge.")
            max_iterations = 100
        """
        very first iterations ( setting iteration )
        """
        actual_icognitos_values = calculate_icognitos_values(actual_icognitos_values, last_icognitos_values)

        sol_vector = Vector(*actual_icognitos_values.values(), icognito=self.custom_icognito)

        last_vector = initial_vector

        dif_module = (sol_vector - last_vector).module

        
        iteration = 1
        iterations = [
            {'iter': 0, f'{self.custom_icognito}(0)': last_vector},
            {'iter': 1, f'{self.custom_icognito}(1)': sol_vector, 
            f'||{self.custom_icognito}({iteration}) - {self.custom_icognito}({iteration-1})||':  round(dif_module, self.precision)}
            ]


        while dif_module > TOL and iteration < max_iterations:
            
            iteration += 1
            
            # Update the system
            last_vector = deepcopy(sol_vector)

            # Update the icognitos values
            last_icognitos_values = deepcopy(actual_icognitos_values)
            actual_icognitos_values = calculate_icognitos_values(actual_icognitos_values, last_icognitos_values)
            
            # Update solution values.
            sol_vector = Vector(*actual_icognitos_values.values(), icognito=self.custom_icognito)

            dif_module = (sol_vector - last_vector).module

            iterations.append({
                'iter': iteration,
                f'{self.custom_icognito}({iteration})': sol_vector,
                f'||{self.custom_icognito}({iteration}) - {self.custom_icognito}({iteration-1})||':  round(dif_module, self.precision)
            })
    

        return iterations

if __name__ == '__main__':

    sys_1 = [
        [2, 1, -1],
        [-1, 1, 1],
        [1, 3,- 1],
    ]

    b_vec_1 = Vector(4, 0, 8)

    sys_2 = [
        [2, -1, 1, 0],
        [0, 3, 1, 2],
        [-1, 0, 3, -1],
        [1, 0, -3, 2],
    ]

    b_vec_2 = Vector(1, 1, -5, 6)

    sys_3 = [
        [2, 1, 1, 0],
        [4, 3, 3, 1],
        [8, 7, 9, 5],
        [6, 7, 9, 8],
    ]

    b_vec_3 = Vector(1, 2, 4, 5)

    sys_4 = [
        [ 2,  3, -1,  0],
        [ 3,  1,  1, -1],
        [ 4,  2,  1,  2],
        [-2,  0,  1,  1],
    ]
    b_vec_4 = Vector(-1, 1, 4, -1)




    lin_sys = Linear_System(sys_4, b_vec_4)
    
    print(lin_sys.gauss_elimination_method(), sep='\n')
