from Matrix import Matrix
import math
from typing import Any
import sympy as sym
from copy import deepcopy

i, j, k = sym.symbols('i j k')

class Vector (list):

    def __init__ (self, *args, **kwargs):

        icognito_symbol = kwargs.get('icognito')

        if not icognito_symbol:
            icognito_symbol = 'x'
    
        for index in range(len(args)):
            self.append(args[index])
            self.__dict__[f'{icognito_symbol}{index + 1}'] = args[index]


    @property
    def dimension (self):
        return len(self.entrys)
    
    @property
    def entrys (self):
        entrys = [entry for entry in self]

        return entrys

    def __repr__(self) -> str:
        string = ' '.join([f"{key}={entry}" for key, entry in self.__dict__.items()])
        return f'Vector({string})'
    
    @property
    def module (self) -> float:
        """
        Vector's module.
        """
       
        return math.sqrt(sum([math.pow(entry, 2) for entry in self.entrys]))
    
    
    @property
    def direction (self) -> 'Vector':
        """
        Unitary version of the vector.
        """

        return Vector(*[entry/self.module for entry in self.entrys])
    
    @property
    def with_unitary_vectors (self) -> 'sym.Add':
        """
        Representation of the vectors with i, j and k.
        """
        res = sum([self.entrys[i]*sym.symbols(f'x{i}') for i in range(self.dimension)])

        return res
    
    @property
    def column_matrix (self) -> 'Matrix':
        """
        Column matrix vector representation.
        """
        column_vector =[
            [entry] for entry in self.entrys
        ] 
        return Matrix(column_vector)
    
    @property
    def row_matrix (self) -> 'Matrix':
        """
        Row matrix vector representation.
        """
        row_vector = [
            [*self.entrys]
        ]
        return Matrix(row_vector)

    def scalar (self, other_vector: 'Vector') -> float:
        """
        Scalar product between the vectors.
        """

        if other_vector.dimension != self.dimension:
            raise ArithmeticError("Can not make scalar product between vectors with different dimensions")
        
        return sum([self.entrys[i]*other_vector.entrys[i] for i in range(self.dimension)])

    def vectorial (self, other_vector: 'Vector') -> 'Vector':
        """
        Vectorial product between the vectors.
        """
        from Matrix import Determinant
        
        if self.dimension > 3 or other_vector.dimension != self.dimension or self.dimension < 2 and other_vector.dimension < 2:
            raise ArithmeticError("Can not make vetorial product in dimension greater than 3")
        

        if self.dimension == 2:
            z1 = 0
            z2 = 0
        else:
            z1 = self[2]
            z2 = other_vector[2]

        res = Determinant.det([
            [i, j, k],
            [self[0], self[1], z1],
            [other_vector[0], other_vector[1], z2],
        ])

        x, y, z = res.coeff(i), res.coeff(j), res.coeff(k)
        return Vector(x, y, z)

    def angule_between (self, other_vector: 'Vector') -> float:
        """
        Cos(0) = (u.v)/(||u||.||v||)
        """
        teta = math.acos((self * other_vector) / (self.module * other_vector.module))
        return round(math.degrees(teta), 2)

    def __mul__ (self, other_vector: 'Vector' or float or 'sym.Symbol') -> float or 'Vector':
        """
        Interpretate multiplication as a scalar product.
        """
        if isinstance(other_vector, Vector):
            return self.scalar(other_vector)
        
        return Vector(*[entry * other_vector for entry in self.entrys])
    
    
    
    def __truediv__ (self, num: float) -> float or 'Vector':
        """
        Interpretate multiplication as a scalar product.
        """
        if isinstance(num, Vector):
            raise ArithmeticError('Can not divide one vector by other')
        
        return Vector(*[entry / num for entry in deepcopy(self.entrys)])
    
    def __lt__(self, other_vector: 'Vector') -> bool:
        """
        Angle between the vectors.
        """
        return self.angule_between(other_vector)

    def __gt__(self, other_vector: 'Vector') -> bool:
        """
        Complementary angle between the vectors
        """
        return 180 - self.angule_between(other_vector)
    
    def __add__ (self, other_vector: 'Vector') -> 'Vector':
        """
        Addition of vectors.
        """
        if self.dimension != other_vector.dimension:
            raise ArithmeticError('Can not sum different dimensions vectors!')
        
        return Vector(*[self.entrys[i] + other_vector.entrys[i] for i in range(len(self.entrys))])
    
    def __sub__ (self, other_vector: 'Vector') -> 'Vector':
        """
        Subtraction of vectors.
        """
        if self.dimension != other_vector.dimension:
            raise ArithmeticError('Can not sum different dimensions vectors!')
        
        return Vector(*[self.entrys[i] - other_vector.entrys[i] for i in range(len(self.entrys))])
    
    
    def __round__ (self, decimals: int):
        return Vector(*[round(entry, decimals) for entry in self.entrys])
    
    def reverse (self) -> 'Vector':
        from copy import deepcopy
        reversed_entrys = deepcopy(self.entrys)
        reversed_entrys.reverse()
        return Vector(*reversed_entrys)

    def invert (self) -> 'Vector':
        return Vector(*(self*-1).entrys)
    
    
    
    
if __name__ == '__main__':
    
    v = Vector(*[1,2])
    u = Vector(2, 3)
    print(v.vectorial(u))
    print(v.row_matrix * u.column_matrix)

    a = sym.symbols('a')

    eq = a + 2

    eq.args