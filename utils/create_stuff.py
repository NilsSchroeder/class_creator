from classes.FF2 import FF2
from classes.Graph import Graph
from utils import other
from utils import compute_stuff as cos
from utils import write_files as wf

import random


#creates a list containing all function g, that are EA-equivalent to f


#creates all the EL_mappings for functions f:FFn->FFm
def create_all_EL_mappings(n, m):
    n_list = create_all_linear_permutations(n)
    m_list = create_all_linear_permutations(m)
    l_list = create_all_linear_functions(n, m)
    el_list = []
    zero = zero_matrix(n, m)
    for n_mat in n_list:
        for m_mat in m_list:
            for l_mat in l_list:
                el_list.append(create_block_matrix(n_mat, zero, l_mat, m_mat))
    return el_list


#creates all linear permutations on FFn
def create_all_linear_permutations(n):
    permutation_list = []
    for i in range(2**(n**2)):
        number = other.decimal_to_binary(i, n**2)
        mat = []
        for j in range(n):
            mat.append(number[j*n:(j*n)+n])
        mat = FF2(mat)
        if mat.compute_det() != 0:
            permutation_list.append(mat)
    return permutation_list

#creates a list containing all linear functions FFn\rightarrow\FFm
def create_all_linear_functions(n,m):
    linear_list = []
    for i in range(2**(n*m)):
        number = other.decimal_to_binary(i, n * m)
        mat = []
        for j in range(m):
            mat.append(number[(j*n):(j*n)+n])
        mat = FF2(mat)
        linear_list.append(mat)
    return linear_list

#creates a list containing the vectors of the canonical basis {e_0,e_1,...,e_(n-1)}
def create_canonical_basis(n):
    L = []
    for i in range(n):
        vec = []
        for j in range(n):
            if j==i:
                vec.append(1)
            else:
                vec.append(0)
        L.append(FF2(vec))
    return L

#creates a matrix with n rows and m columns, containing only zeroes.
def zero_matrix(n,m):
    mat = []
    for i in range(n):
        row = []
        for j in range(m):
            row.append(0)
        mat.append(row)
    return FF2(mat)


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

#creates a block matrix containing the matrices [a & b \\ c & d]
def create_block_matrix(a,b,c,d):
    if (a.n != b.n or a.m != c.m or b.m != d.m or c.n != d.n):
        raise ValueError("ERROR. Dimensions of matrices don't match!")
        return None
    else:
        mat = []
        for i in range(a.n):
            row = a.get_row(i) + b.get_row(i)
            mat.append(row)
        for i in range (c.n):
            row = c.get_row(i)  + d.get_row(i)
            mat.append(row)
        mat = FF2(mat)
        return mat

#removes redundant list (here class) elements
def remove_redundancies(sorted_list):
    new_list =  []
    i = 0
    length = len(sorted_list)
    while i < length:
        while i+1<length and sorted_list[i] == sorted_list[i+1]:
            i += 1
        new_list.append(sorted_list[i])
        i += 1
    if sorted_list[length-1] not in new_list:
        new_list.append(sorted_list[length-1])
    return new_list


#creates the EL-class for a given graph
def create_EL_class(graph):
    el_list = create_all_EL_mappings(graph.dom_dim, graph.co_dim)
    el_class = [graph]
    for el in el_list:
        new_graph = el * graph
        new_graph = Graph(new_graph, case=1)
        el_class.append(new_graph)
    el_class.sort()
    el_class = remove_redundancies(el_class)
    return el_class


#creates the EA-class for a given graph
def create_EA_class(f_graph, el_class = None):
    if el_class is None:
        el_class = create_EL_class(f_graph)
    mat_list = create_vector_matrices(f_graph.dom_dim,f_graph.co_dim)
    ea_class = []
    for mat in mat_list:
        for graph in el_class:
            new_graph = Graph(graph+mat, case=1)
            ea_class.append(new_graph)
    ea_class.sort()
    ea_class = remove_redundancies(ea_class)
    return ea_class

#creates a list containing matrices. THese can be added to a graph of f to get the function f+v.
def create_vector_matrices(n,m):
    V = cos.create_vectorspace(n+m)
    mat_list = []
    for v in V:
        mat = []
        row = v.mat[0]
        for i in range(2**n):
            mat.append(row)
        mat_list.append(FF2(mat).transpose())
    return mat_list



def create_random_function(n,m):
    func_list = []
    for i in range(2**n):
        func_list.append(create_random_vector(m))
    def f(x):
        dec = x.get_dec_value()
        return func_list[dec]
    return f


def create_random_vector(m):
    vec = []
    for i in range(m):
        vec.append(random.randint(0, 1))
    return FF2(vec)


#creates a list, containing all EL-classes on F(n,m)
def create_all_EL_classes(n, m):
    class_list = []
    functions = 0
    max_functions = (2**m)**(2**n)
    while functions < max_functions:
        f = create_random_function(n, m)
        f_graph = Graph(f, n)
        exists = False
        for class_el in class_list:
            if f_graph in class_el:
                exists = True
                break
        if not exists:
            el_class = create_EL_class(f_graph)
            el_class.sort()
            functions = functions + len(el_class)
            class_list.append(el_class)
    return class_list


#creates a list, containing all EA-classes on F(n,m)
def create_all_EA_classes(n, m, el_list = None):
    class_list = []
    functions = 0
    max_functions = (2**m)**(2**n)
    if el_list is None:
        el_list = create_all_EL_classes(n,m)
    while functions < max_functions:
        f = create_random_function(n, m)
        f_graph = Graph(f, n)
        exists = False
        for class_ea in class_list:
            if f_graph in class_ea:
                exists = True
                break
        if not exists:
            for el_class in el_list:
                if f_graph in el_class:
                    ea_class = el_class
                    break
            ea_class = create_EA_class(f_graph, ea_class)
            functions = functions + len(ea_class)
            class_list.append(ea_class)
    return class_list

def create_all_EA_classes_case_2(n,m):
    class_list = []
    functions = 0
    max_functions = (2 ** m) ** (2 ** n)
    while functions < max_functions:
        f = create_random_function(n, m)
        f_graph = Graph(f, n)
        exists = False
        for class_ea in class_list:
            if f_graph in class_ea:
                exists = True
                break
        if not exists:
            ea_class = create_EA_class(f_graph)
            functions = functions + len(ea_class)
            class_list.append(ea_class)
    return class_list


#creates the set v containing vectors (x,0)
def create_v(n,m):
    V = cos.create_vectorspace(n)
    v = []
    zero = zero_matrix(1,m)
    for x in V:
        z = x.mat[0]+zero.mat[0]
        v.append(FF2(z).transpose())
    return v

#creates the set vt containing vectors (0,y)
def create_vt(n,m):
    V = cos.create_vectorspace(m)
    vt = []
    zero = zero_matrix(1, n)
    for x in V:
        z =  zero.mat[0] + x.mat[0]
        vt.append(FF2(z).transpose())
    return vt

#creates the Walsh Zeroes of a given function f:FFn->FFm
def create_IDS(f,n,m):
    IDS = [zero_matrix(1,n+m).transpose()]
    V = cos.create_vectorspace(n)
    W = cos.create_vectorspace(m)
    for a in V:
        row = []
        for b in W:
            ddt = 0
            for x in V:
                if (f(x) + f(x + a) == b):
                    ddt += 1
            if ddt == 0:
                x = FF2(a.mat[0]+b.mat[0])
                IDS.append(x.transpose())
    return IDS

#creates the linear projection rho_(i,j):FFn->FFm.
def create_lin_proj(i, j, n, m):
    rho = []
    for k in range(m):
        row = []
        for l in range(n):
            if k == l and k >= i and k<j:
                row.append(1)
            else:
                row.append(0)
        rho.append(row)
    return(FF2(rho))

# creates the t-twist matrix M_t
def create_Mt(t, n, m):
    rho_1 = create_lin_proj(t, n, n, n)
    rho_2 = create_lin_proj(0, t, m, n)
    rho_3 = create_lin_proj(0, t, n, m)
    rho_4 = create_lin_proj(t, m, m, m)
    M_t = create_block_matrix(rho_1, rho_2, rho_3, rho_4)
    return M_t

# creates all CCZ classes on FFn->FFm. First creates all EA-classes, then uses Proposition 5 and Theorem 3.
def create_all_CCZ_classes(n,m, class_list = None):
    if class_list == None:
        class_list = create_all_EA_classes(n, m)
    du_list = []
    max_t = min(n, m)
    vt = create_vt(n, m)
    for ea_class in class_list:
        du_list.append(cos.compute_class_du(ea_class))
    fused_classes = [] # fused_classes counts all classes that have been fused to an CCZ-class already. Thus we don't need to perform t-twist on a class that has already been fused.
    ccz_classes = []
    t_twist_pairs = []
    for i, ea_class in enumerate(class_list):
        class_du = cos.compute_class_du(ea_class)
        if i not in fused_classes and class_du in du_list[i:]:
            ccz_class = ea_class
            keep_going = True
            count = 0
            while len(fused_classes) < len(class_list) and count < len(ea_class):
                f_graph = ea_class[count]
                IDS = create_IDS(f_graph.f, n, m)
                for t in range(1, max_t + 1):

                    #perm states if T_y is a permutation, thus if M_t is admissible
                    M_t = create_Mt(t, n, m)
                    perm = True
                    for z in vt:
                        y = M_t * z
                        if y not in IDS:
                            perm = False
                            break
                    if perm:
                        g_graph = Graph(M_t * f_graph, case=1)
                        if g_graph not in ccz_class:
                            t_twist_pairs.append([f_graph, g_graph])
                            for j, class_ea in enumerate(class_list):
                                if j > i:
                                    if g_graph in class_ea:
                                        ccz_class = ccz_class + class_ea
                                        fused_classes.append(j)
                keep_going = False
                count = count + 1
                for j, class_ea in enumerate(class_list[i:]):
                    if j not in fused_classes and cos.compute_class_du(class_ea) == class_du:
                        keep_going = True
                        break
                fused_classes.append(i)
            ccz_class.sort()
            ccz_classes.append(ccz_class)
    return ccz_classes


# creates a list of all t-twist pairs on a functions list.
def create_t_twist_pairs(functions):
    t_twist_pairs = []
    n = functions[0].dom_dim
    m = functions[0].co_dim
    vt = create_vt(n, m)
    max_t = min(n, m)
    t_matrices = []
    for t in range(1, max_t + 1):
        # perm states if T_y is a permutation, thus if M_t is admissible
        M_t = create_Mt(t, n, m)
        t_matrices.append(M_t)

    for f_graph in functions:
        IDS = create_IDS(f_graph.f, n, m)
        for t, M_t in enumerate(t_matrices):
            perm = True
            for z in vt:
                y = M_t * z
                if y not in IDS:
                    perm = False
                    break
            if perm:
                g_graph = Graph(M_t * f_graph, case=1)
                new = True
                for pair in t_twist_pairs:
                    if (g_graph == pair[0] and f_graph == pair[1]) or (f_graph == pair[0] and g_graph == pair[1]):
                        new = False
                        break
                if new:
                    t_twist_pairs.append([f_graph, g_graph, t])

    return t_twist_pairs

# Gives us all the results we want for given n and m
def create_everything(n,m, case = 0):
    if case == 0:
        el_list = create_all_EL_classes(n, m)
        for i, el_class in enumerate(el_list):
            wf.graphs_to_file(el_class, f"EL_class_{i}")
            t_twist_pairs = create_t_twist_pairs(el_class)
            wf.t_twist_pairs_to_file(t_twist_pairs, f"EL_t_twist_pairs_{i}", n=min(n,m))

        ea_list = create_all_EA_classes(n,m)
        for i, ea_class in enumerate(ea_list):
            wf.graphs_to_file(ea_class, f"EA_class_{i}")
            t_twist_pairs = create_t_twist_pairs(ea_class)
            wf.t_twist_pairs_to_file(t_twist_pairs, f"EA_t_twist_pairs_{i}", n=min(n,m))

        ccz_list = create_all_CCZ_classes(n, m, ea_list)
        for i, ccz_class in enumerate(ccz_list):
            wf.graphs_to_file(ccz_class, f"CCZ_class_{i}")
            t_twist_pairs = create_t_twist_pairs(ccz_class)
            wf.t_twist_pairs_to_file(t_twist_pairs, f"CCZ_t_twist_pairs_{i}", n=min(n,m))

        functions = []
        for ccz_class in ccz_list:
            functions = functions + ccz_class

        t_twist_pairs= create_t_twist_pairs(functions)
        wf.t_twist_pairs_to_file(t_twist_pairs, "all_t_twist_pairs", n=min(n,m))

    elif case == 1:
        ea_list = create_all_EA_classes_case_2(n, m)
        for i, ea_class in enumerate(ea_list):
            wf.graphs_to_file(ea_class, f"EA_class_{i}")
            t_twist_pairs = create_t_twist_pairs(ea_class)
            wf.t_twist_pairs_to_file(t_twist_pairs, f"EA_t_twist_pairs_{i}", n=min(n, m))

        ccz_list = create_all_CCZ_classes(n, m, ea_list)
        for i, ccz_class in enumerate(ccz_list):
            wf.graphs_to_file(ccz_class, f"CCZ_class_{i}")
            t_twist_pairs = create_t_twist_pairs(ccz_class)
            wf.t_twist_pairs_to_file(t_twist_pairs, f"CCZ_t_twist_pairs_{i}", n=min(n, m))

        functions = []
        for ccz_class in ccz_list:
            functions = functions + ccz_class

        t_twist_pairs = create_t_twist_pairs(functions)
        wf.t_twist_pairs_to_file(t_twist_pairs, "all_t_twist_pairs", n=min(n, m))
