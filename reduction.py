from polynomial import *
import numpy as np
import math
import random

def reduce(f, g):
    """
    Compute the one step reduction of polynomial f with respect to g.
    @param f: Polynomial f.
    @param g: Polynomial g.
    @return: The one step reduction represented as a polynomial object if the leading term of g divides the leading term of f,
    return False if not.
    """

    ratio = (f.lt()).divide(g.lt())
    if ratio:
        return f.subtract(ratio.multiply(g))
    else:
        return False


def reduce_lst(f, G):
    """
    Compute the complete reduction of polynomial f with respect to a set of polynomials G. At each step, randomly choose a polynomial
    g in G such that the leading term of g divides the leading term of r.
    @param f: Polynomial f.
    @param G: list of polynomials G.
    @return: The result of complete reduction represented as a polynomial object.
    """

    r = f
    lst = []
    for g in G:
        if (r.lt()).divide(g.lt()):
            lst.append(g)
    while len(lst) > 0:
        random_g = random.choice(lst)
        r = reduce(r, random_g)
        lst = []
        for g in G:
            if (r.lt()).divide(g.lt()):
                lst.append(g)
    return r
