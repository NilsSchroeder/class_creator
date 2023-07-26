import numpy as np
import sage
import utils.other as other


class FF2:
    def __init__(self, mat):
        if isinstance(mat, FF2):
            self.n = mat.n
            self.m = mat.m
            self.mat = mat.mat
        else:
            if type(mat[0]) == int:
                mat = [mat]
            self.n = len(mat)
            self.m = len(mat[0])
            b = []
            for i in range(self.n):
                if len(mat[i]) != self.m:
                    print("ERROR. All rows must be of same length.")
                    return None
                else:
                    row = []
                    for j in range(self.m):
                        row.append(mat[i][j] % 2)
                    b.append(row)
            self.mat = b

    def __add__(self, y):
        z = []
        if 0 in [self.n, self.m]:
            print("ERROR")
            return None
        elif self.n != y.n:
            print("ERROR. Matrices must have the same number of rows!")
            return None
        else:
            if self.m != y.m:
                print("ERROR. Matrices must have the same number of columns!")
                return None
            else:
                for i in range(self.n):
                    row = []
                    for j in range(self.m):
                        row.append((self.mat[i][j] + y.mat[i][j]) % 2)
                    z.append(row)
                return FF2(z)
        return None

    def __iadd__(self, y):
        z = []
        if 0 in [self.n, self.m]:
            print("ERROR")
            return None
        elif self.n != y.n:
            print("ERROR. Matrices must have the same number of rows!")
            return None
        else:
            if self.m != y.m:
                print("ERROR. Matrices must have the same number of columns!")
                return None
            else:
                for i in range(self.n):
                    row = []
                    for j in range(self.m):
                        row.append((self.mat[i][j] + y.mat[i][j]) % 2)
                    z.append(row)
                return FF2(z)

    def __mul__(self, y):
        n = self.n
        m = y.m
        z = []
        if self.m != y.n:
            print("ERROR. The matrices have the wrong dimensions and can't be multiplied")
            return None
        else:
            for i in range(n):
                row = []
                for j in range(m):
                    elem = 0
                    for k in range(self.m):
                        elem += self.mat[i][k] * y.mat[k][j]
                    row.append(elem)
                z.append(row)
            return FF2(z)
        return None

    def __pow__(self, y):
        if y - int(y) != 0:
            raise ValueError("ERROR. Exponent must be an integer.")
            return None
        elif self.n != self.m:
            raise ValueError("ERROR. Matrix must be square.")
        else:
            if y == 0:
                return FF2(identity_matrix(self.n))
            elif y > 0:
                temp = self
                for i in range(1, y):
                    temp = self * temp
                return temp
            else:
                inv = self.invert()
                y *= -1
                temp = inv
                for i in range(1, y):
                    temp = inv * temp
                return temp

    def __eq__(self, y):
        return (self.mat == y.mat)

    def __ne__(self, y):
        return (self.mat != y.mat)

    def __gt__(self, y):
        if self.n != y.n:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        if self.m != y.m:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        n = self.n
        m = self.m
        for i in range(n):
            for j in range(m):
                if self.mat[i][j] > y.mat[i][j]:
                    return True
                if self.mat[i][j] < y.mat[i][j]:
                    return False
        return False

    def __ge__(self, y):
        if self.n != y.n:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        if self.m != y.m:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        n = self.n
        m = self.m
        for i in range(n):
            for j in range(m):
                if self.mat[i][j] > y.mat[i][j]:
                    return True
                if self.mat[i][j] < y.mat[i][j]:
                    return False
        return True

    def __lt__(self, y):
        if self.n != y.n:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        if self.m != y.m:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        n = self.n
        m = self.m
        for i in range(n):
            for j in range(m):
                if self.mat[i][j] < y.mat[i][j]:
                    return True
                if self.mat[i][j] > y.mat[i][j]:
                    return False
        return False

    def __le__(self, y):
        if self.n != y.n:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        if self.m != y.m:
            raise ValueError("ERROR. Matrices must be of same dimension.")
        n = self.n
        m = self.m
        for i in range(n):
            for j in range(m):
                if self.mat[i][j] < y.mat[i][j]:
                    return True
                if self.mat[i][j] > y.mat[i][j]:
                    return False
        return True


    def transpose(self):
        mat = self.mat
        n = self.n
        m = self.m
        z = []
        for j in range(m):
            row = []
            for i in range(n):
                row.append(mat[i][j])
            z.append(row)
        return FF2(z)

    def invert(self):
        if self.compute_det():
            cof = self.compute_cof()
            inv = cof.transpose()
            return inv
        else:
            print("ERROR. Matrix is not invertible")
            return None

    # computes the cofactormatrix
    def compute_cof(self):
        mat = []
        n = self.n
        for i in range(n):
            row = []
            for j in range(n):
                sub_mat = self.remove_row_and_column(i, j)
                elem = (-1) ** (i + j + 2) * sub_mat.compute_det()
                row.append(elem)
            mat.append(row)
        return FF2(mat)

    def elem(self, i, j=None):
        if self.n == 1:
            return self.mat[0][i]
        else:
            if j is None:
                return self.mat[i]
            else:
                return self.mat[i][j]

    def get_all_elems(self):
        elems = []
        for i in range(self.n):
            for j in range(self.m):
                elems.append(self.elem(i, j))

    # computes the determinant of a matrix using the Laplace expansion
    def compute_det(self):
        n = self.n
        if self.n != self.m:
            print("determinant is not defined for a non-square matrix.")
            return None
        else:
            if n == 1:
                mat = self.mat
                return mat[0][0]
            elif n == 2:
                mat = self.mat
                det = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
                return det % 2
            elif n == 3:
                mat = self.mat
                det = 0
                for i in range(self.n):
                    det += mat[i][0] * mat[(1 + i) % 3][1] * mat[(2 + i) % 3][2]
                det -= mat[0][0] * mat[2][1] * mat[1][2]
                det -= mat[2][0] * mat[1][1] * mat[0][2]
                det -= mat[1][0] * mat[0][1] * mat[2][2]
                return det % 2
            else:
                det = 0
                for j in range(self.n):
                    mat = self.remove_row_and_column(0, j)
                    det += self.mat[0][j] * (-1) ** (j + 2) * mat.compute_det()
                return det % 2

    # removes the i-th row and j-th column of a matrix
    def remove_row_and_column(self, row_num=None, col_num=None):
        mat = []
        for i in range(self.n):
            if i != row_num:
                row = []
                for j in range(self.m):
                    if j != col_num:
                        row.append(self.mat[i][j])
                mat.append(row)
        return FF2(mat)

    #return the i-th row of a matrix
    def get_row(self, i = 0):
        if (self.n-1 < i):
            print("ERROR. Row-index out of bound.")
            return None
        else:
            return self.mat[i]

    # return the j-th column of a matrix
    def get_column(self, j):
        if (self.m-1 < j):
            print("ERROR. Column-index out of bound.")
            return None
        else:
            column = []
            for i in range(self.n):
                column.append(self.mat[i][j])
            return column

    def print(self):
        for i in range(self.n):
            s = "["
            for j in range(self.m):
                s += str(self.mat[i][j])
                if j < (self.m - 1):
                    s += " "
            s += "]"
            print(s)
        print("")
        return None



    def get_dec_value(self):
        if len(self.mat) == 1:
            return other.binary_to_decimal(self.mat[0])
        else:
            raise ValueError("ERROR. Can only give decimal representation of vectors.")

def identity_matrix(n):
    mat = []
    for i in range(n):
        row = []
        for j in range(n):
            if i == j:
                row.append(1)
            else:
                row.append(0)
        mat.append(row)
    return FF2(mat)