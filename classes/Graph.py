import numpy as np
import sage
import utils.compute_stuff as cos
import utils.create_stuff as crs
import utils.other as other
from classes.FF2 import FF2
from torch.utils.data.datapipes._typing import Integer
from examples import functions as func

class Graph(FF2):
    # takes a function f:FF2^n->FF2^m and creates its graph

    def __init__(self, f, n = None, case = 0):
        # case == 0 means f is a properly defined function
        if case == 0:
            if n is None:
                raise ValueError("ERROR. Give the dimension n of the codomain.")
                return None
            V = cos.create_vectorspace(n)
            graph = []
            for v in V:
                val = FF2(f(v))
                elem = None
                if isinstance(val, FF2):
                    elem = v.mat[0] + val.mat[0]
                elif type(val) == list and type(val[0]) == Integer:
                    elem = v.mat[0] + val
                elif type(val) == list and type(val[0]) == list:
                    elem = v.mat[0] + val[0]
                else:
                    raise ValueError('Error. Check if your function has the correct return values.')
                    return None
                m = len(elem)
                graph.append(elem)

            # this step is necessary to transpose self.mat
            new_graph = []
            for j in range(m):
                row = []
                for i in range(2 ** n):
                    row.append(graph[i][j])
                new_graph.append(row)

            self.mat = new_graph
            self.f = f
            self.n = m
            self.m = 2 ** n
            self.dom_dim = n  # dimension of the domain (FF_2^(dom_dim))
            self.co_dim = m - n  # dimension of the codomain

        #case==1 means f is a graph that has not been sorted yet and is only an FF2 element.
        elif case == 1:
            if isinstance(f, FF2):

                if n is None:
                    n = 0
                    temp_n = f.m
                    while temp_n > 1:
                        temp_n = temp_n/2
                        n += 1

                f = other.sort_for_graph(f, n)

                self.mat = f.mat
                self.n = f.n
                self.m = f.m

                mat = FF2(f.mat)
                mat = mat.transpose()
                new_mat = []
                for i in range(mat.n):
                    new_mat.append(mat.mat[i][n:])
                def func(x):
                    dec = x.get_dec_value()
                    return FF2(new_mat[dec])

                self.f = func
                self.co_dim = self.n - n
                self.dom_dim = n

            else:
                raise ValueError('f is not a matrix. Check if you gave the right case.')
        # case==2 means f is a matrix f:FFn->FFm
        elif case == 2:
            return

    def __lt__(self,other):
        V = cos.create_vectorspace(self.dom_dim)
        for v in V:
            if self.f(v) < other.f(v):
                return True
            if self.f(v) > other.f(v):
                return False
        return False

    def __le__(self, other):
        V = cos.create_vectorspace(self.dom_dim)
        for v in V:
            if self.f(v) < other.f(v):
                return True
            if self.f(v) > other.f(v):
                return False
        return  True

    def get_func(self):
        temp_fun = FF2(self.mat)
        temp_fun = temp_fun.transpose()
        str = f""
        for v in temp_fun.mat:
            str = f"{str} {v[0:self.dom_dim]} --> {v[self.dom_dim:]} \n"
        return f"{str} \n"

    def print_func(self):
        str = self.get_func()
        print(str)


    #return the differential uniformity of a the function graph.f
    def get_du(self):
        f = self.f
        return cos.compute_du(f,self.dom_dim,self.co_dim)







