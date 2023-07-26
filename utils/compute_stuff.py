from utils import other
from classes.FF2 import FF2


#Creates a list containing all vectors in F2^n
def create_vectorspace(n):
    V = []
    for x in range(2**n):
        v = other.decimal_to_binary(x, n)
        V.append(FF2(v))
    return V

def compute_DDT(f ,n ,m):
    DDT = []
    V = create_vectorspace(n)
    W = create_vectorspace(m)
    for a in V:
        row = []
        for b in W:
            ddt = 0
            for x in V:
                if (f(x) +f(x+a)== b):
                    ddt += 1
            row.append(ddt)
        DDT.append(row)
    return DDT

#takes a class of graphs and computes its du
def compute_class_du(L):
    graph = L[0]
    return compute_du(graph.f, graph.dom_dim, graph.co_dim)

def compute_LAT(f, n, m):
    LAT = []
    V = create_vectorspace(n)
    W = create_vectorspace(m)
    for a in V:
        row = []
        for b in W:
            lat = 0
            for x in V:
                exp = (a * x.transpose() + b * f(x).transpose()).elem(0)
                lat += (-1) ** (exp)
            row.append(lat)
        LAT.append(row)

    return LAT


#computes the differential uniformity of a function f:FF2^n->FF2^m
def compute_du(f,n,m):
    DDT = compute_DDT(f,n,m)
    max = 0
    for i in range(1,2**n):
        for j in range(2**m):
            if DDT[i][j]>max:
                max = DDT[i][j]
    return max



# computes the linearity of a function f:FF2^n->FF2^m
def compute_lin(f, n, m):
    LAT = compute_LAT(f, n, m)
    max = 0
    for i in range(2 ** n):
        for j in range(1, 2 ** m):
            if LAT[i][j] > max:
                max = LAT[i][j]
    return max






