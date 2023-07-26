import utils.compute_stuff as cos

import utils.create_stuff as crs
import utils.other as other
from classes.FF2 import FF2
from classes.Graph import Graph
from torch.utils.data.datapipes._typing import Integer
from examples import functions as func
import os
from pathlib import Path
from io import StringIO

def graphs_to_file(some_class, file_name='Some_class', file_type='txt'):
    n = some_class[0].dom_dim
    m = some_class[0].co_dim
    if file_type == 'txt':
        with open(f"{file_name}.txt", "w") as f:
            f.write(f"The class contains {len(some_class)} functions. \n"
                    f"It has a differential uniformity of {cos.compute_du(some_class[0].f, some_class[0].dom_dim, some_class[0].co_dim)}.\n\n")
            for graph in some_class:
                f.write(graph.get_func())
    else:
        print('ERROR. Give correct file type (txt)')


def t_twist_pairs_to_file(some_class, file_name='Some_class', file_type='txt', n = None):
    equals = 0
    if len(some_class)>0:
        dom_dim = some_class[0][0].dom_dim
        m = some_class[0][0].co_dim
        if n is not None:
            du_list = crs.zero_matrix(1, 2**(n-1)).mat[0]
            t_list = crs.zero_matrix(1, n).mat[0]
        if file_type == 'txt':
            with open(f"{file_name}.txt", "w") as f:

                count = 0
                for i, pair in enumerate(some_class):
                    if pair[0] == pair[1]:
                        equals += 1
                    t_list[pair[2]-1] += 1
                    if cos.compute_du(pair[0].f, pair[0].dom_dim, pair[0].co_dim) != cos.compute_du(pair[1].f, pair[0].dom_dim, pair[0].co_dim):
                        count = count + 1
                    f.write(f"Pair number {i} has a differential uniformity of {cos.compute_du(pair[0].f,pair[0].dom_dim, pair[0].co_dim)}.\n")
                    f.write(pair[0].get_func())
                    f.write(pair[1].get_func())
                    f.write("\n\n")
                    if n is not None:
                        for k in range(2**(n-1)):
                            if cos.compute_du(pair[0].f, pair[0].dom_dim, pair[0].co_dim) == 2*(k+1):
                                du_list[k] += 1
                if n is not None:
                    f.write(f"The list contains {len(some_class)} pairs. \n")
                    for k in range(2**(n-1)):
                        f.write(f"There are {du_list[k]} pairs with a differential uniformity of {2*(k+1)}.\n")
                    for t in range(n):
                        f.write(f"There are {t_list[t]} {t+1}-twists.\n")
                    f.write(f"{equals} functions are t-twist equivalent to themself.")
        else:
            print('ERROR. Give correct file type (txt)')



def FF2_to_file(some_list, file_name='Some_class', file_type='txt'):
    new_list = []
    comp_list = []
    for mat in some_list:
        graph = Graph(mat)
        if graph.get_func() not in comp_list:
            new_list.append(graph)
            comp_list.append(graph.get_func())
    graphs_to_file(new_list, file_name, file_type)



def old_FF2_to_file(some_list, file_name='Some_class', file_type='txt'):
    new_list =[]
    if file_type == 'txt':
        with open(f"{file_name}.txt", "w") as f:
            f.write(f"The list contains {len(some_list)} matrices. \n\n")
            for mat in some_list:
                graph = Graph(mat)
                f.write(graph.get_func())
    else:
        print('ERROR. Give correct file type (xml, tei or txt)')

