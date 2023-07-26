from classes.FF2 import FF2


def f1(x):
    z = [[x.elem(0)*x.elem(1), x.elem(2)]]
    return FF2(z)


#Example Function twisting f:FF⁵\rightarrowFF⁴
def f2(x):
    z = [[x.elem(0)+x.elem(2), x.elem(1)+x.elem(3), x.elem(0)*(x.elem(2)+x.elem(3)),
          x.elem(0)*x.elem(2)+x.elem(1)*x.elem(3)]]
    return FF2(z)



def g2(x):
    z = [x.elem(0)+x.elem(2), x.elem(1)+x.elem(3), (x.elem(0)+x.elem(2))*(x.elem(2)+x.elem(3)), x.elem(2)*(x.elem(0) +
         x.elem(2)) + x.elem(3)*(x.elem(1)+x.elem(3))]
    return FF2(z)

def f3(x):
    z = [x.elem(0)+x.elem(1),x.elem(0)*x.elem(1),1]
    return FF2(z)

def f4(x):
    z = [x.elem(0)*x.elem(1) + 1, x.elem(0)+x.elem(1)]
    return FF2(z)

def f5(x):
    z = [0,x.elem(0)+x.elem(1)]
    return FF2(z)
