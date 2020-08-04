import polynomial as poly
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

    assert isinstance(f, poly.Polynomial), 'The input must be a polynomial.'
    assert isinstance(g, poly.Polynomial), 'The input must be a polynomial.'

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

    assert isinstance(f, poly.Polynomial), 'The input must be a polynomial.'
    assert all([isinstance(g, poly.Polynomial) for g in G]), 'The input must be a list of polynomials.'

    num_add = 0
    r = f
    lst = [g for g in G if (r.lt()).divide(g.lt())]
    while len(lst) > 0:
        random_g = random.choice(lst)
        r = reduce(r, random_g)
        num_add += 1
        print(r)
        if not r.is_zero():
            lst = [g for g in G if (r.lt()).divide(g.lt())]
        else:
            lst = []
    return r, num_add


def S(f, g):
    """
    Compute the S polynomial of polynomials f and g.
    @param f: Polynomial f.
    @param g: Polynomial g.
    @return: The S polynomial represented as a polynomial object.
    """

    assert isinstance(f, poly.Polynomial), 'The input must be a polynomial.'
    assert isinstance(g, poly.Polynomial), 'The input must be a polynomial.'

    lt_f = f.lt()
    lt_g = g.lt()
    lcm_fg = lt_f.lcm(lt_g)

    return (f.multiply(lcm_fg.divide(lt_f))).subtract(g.multiply(lcm_fg.divide(lt_g)))
