

def decimal_to_binary(x,n):
    vector = []
    for i in range(0,n):
        if (x % (2**(i+1)) == 0):
            vector.append(0)
        else:
            vector.append(1)
            x -= 2**i
    vector.reverse()
    return vector

def binary_to_decimal(x):
    dec = 0
    for i in range(len(x)):
        dec = dec + x[-(i+1)]*(2**i)
    return dec


def print_list(V):
    s = "["
    for i,v in enumerate(V):
        s += str(v.mat)
        if i < (len(V)-1):
            s += ", "
    s += "]\n"
    print(s)



    #sorts the graph which is only a FF2 element, such that it is in the right order
def sort_for_graph(graph, n):
    new_graph = graph.transpose()
    new_graph.mat.sort(key=lambda x: x[0:n])
    return new_graph.transpose()
