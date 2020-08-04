import polynomial as poly
import reduction as rd
import numpy as np
from itertools import chain
import math
import random

# 还是要用finite field
# 而且，随机选择的monomial可能会相同，导致不会生成binomial

def buchberger_random(F):
    """
    The classic buchberger algorithm using random selection.
    @param F: a list of polynomials.
    @return: The Gröbner basis of the ideal generated by F represented as a list of polynomials.
    """

    assert all([isinstance(f, poly.Polynomial) for f in F]), 'The input must be a list of polynomials.'

    G = F
    P = set([frozenset([F[i], F[j]]) for i in range(len(F)) for j in range(i + 1, len(F))])
    counter = 1
    num_add = 0
    print('------------------')
    while len(P) > 0:
        fg = random.sample(P, 1)[0]
        f, g = tuple(fg)
        print('Iteration {}:'.format(counter))
        print('The choice of pair is {} and {}'.format(f, g))
        P.remove(fg)
        r, new_add = rd.reduce_lst(rd.S(f, g), G)
        num_add += new_add + 1
        print('r is {}'.format(r))
        if not r.is_zero():
            P = P.union([frozenset([f, r]) for f in G])
            G = list(set(G).union([r]))
        counter += 1
        print('Total number of additions is {}'.format(num_add))
        print('------------------')

    return G, num_add


def remove_duplicate(lst):
    """
    Remove duplicates in a list and preserve order.
    @param lst: The list from which to remove duplicates.
    @return: A list with duplicates removed and order preserved.
    """

    seen = set()
    seen_add = seen.add
    return [x for x in lst if not (x in seen or seen_add(x))]


def buchberger_first(F):
    """
    The classic buchberger algorithm using first selection.
    @param F: a list of polynomials.
    @return: The Gröbner basis of the ideal generated by F represented as a list of polynomials.
    """

    assert all([isinstance(f, poly.Polynomial) for f in F]), 'The input must be a list of polynomials.'

    G = F
    P = [frozenset([F[i], F[j]]) for i in range(len(F)) for j in range(i + 1, len(F))]
    P = remove_duplicate(P)
    counter = 1
    num_add = 0
    print('------------------')
    while len(P) > 0:
        fg = P.pop(0)
        f, g = tuple(fg)
        print('Iteration {}:'.format(counter))
        print('The choice of pair is {} and {}'.format(f, g))
        r, new_add = rd.reduce_lst(rd.S(f, g), G)
        num_add += new_add + 1
        print('r is {}'.format(r))
        if not r.is_zero():
            for f in G:
                P.append(frozenset([f, r]))
            P = remove_duplicate(P)
            G = list(set(G).union([r]))
        counter += 1
        print('Total number of additions is {}'.format(num_add))
        print('------------------')

    return G, num_add


def buchberger_degree(F):
    """
    The classic buchberger algorithm using degree selection.
    @param F: a list of polynomials.
    @return: The Gröbner basis of the ideal generated by F represented as a list of polynomials.
    """

    assert all([isinstance(f, poly.Polynomial) for f in F]), 'The input must be a list of polynomials.'

    G = F
    P = [frozenset([F[i], F[j]]) for i in range(len(F)) for j in range(i + 1, len(F))]
    P = remove_duplicate(P)
    total_degree = [(list(x)[0].lt().lcm(list(x)[1].lt())).degree for x in P]
    counter = 1
    num_add = 0
    print('------------------')
    while len(P) > 0:
        min_degree = min(total_degree)
        min_index = total_degree.index(min_degree)
        fg = P[min_index]
        P.remove(fg)
        total_degree.remove(min_degree)
        f, g = tuple(fg)
        print('Iteration {}:'.format(counter))
        print('The choice of pair is {} and {}'.format(f, g))
        r, new_add = rd.reduce_lst(rd.S(f, g), G)
        num_add += new_add + 1
        print('r is {}'.format(r))
        if not r.is_zero():
            for f in G:
                P.append(frozenset([f, r]))
            P = remove_duplicate(P)
            G = list(set(G).union([r]))
            total_degree = [(list(x)[0].lt().lcm(list(x)[1].lt())).degree for x in P]
        counter += 1
        print('Total number of additions is {}'.format(num_add))
        print('------------------')

    return G, num_add


def buchberger_benchmark(n, d, s, N, mode):
    """
    Compute the complexity (number of additions) for the three Buchberger algorithm variants above using N binomial ideals (at most 2 terms)
    in n variables, s generators with maximal degree d.
    @param n: The number of variables.
    @param d: The maximal degree of a generator.
    @param s: The number of generators.
    @param N: The number of polynomial ideals.
    @param mode: Generate ideals using "uniform" or "weighted" sampling.
    """

    assert all([isinstance(x, int) for x in [n, d, s, N]]), 'n, d, s, N should all be integers.'

    coef_max = 20
    ideals = []
    for _ in range(N):
        ideal = []
        for _ in range(s):
            if mode == 'weighted':
                new_ideal = weighted_selection(n, d, coef_max)
            elif mode == 'uniform':
                new_ideal = uniform_selection(n, d, coef_max)
            ideal.append(new_ideal)
        ideals.append(ideal)

    buch_random, buch_first, buch_degree = [], [], []
    counter = 1

    for ideal in ideals:
        buch_random.append(buchberger_random(ideal)[1])
        print('Buchberger random for ideal {} completed'.format(counter))
        buch_first.append(buchberger_first(ideal)[1])
        print('Buchberger first for ideal {} completed'.format(counter))
        buch_degree.append(buchberger_degree(ideal)[1])
        print('Buchberger degree for ideal {} completed'.format(counter))
        counter += 1

    return buch_random, buch_first, buch_degree

def weighted_selection(n, d, coef_max):
    """
    Generate a binomial using weighted selection.
    @param n: The number of variables in the binomial.
    @param d: The maximum degree of the binomial.
    @param coef_max: The maximum coefficient for both monomials in the binomial.
    @return: Generated binomial represented as a polynomial object.
    """

    polynomial = []
    for _ in range(2):
        degree = random.randint(1, d)
        coef = random.randint(1, coef_max)
        monomial = np.array(random.choice(list(all_monomials(n, degree))))
        monomial = np.concatenate([np.array([coef]), monomial])
        polynomial.append(monomial)

    return poly.Polynomial(np.array(polynomial))


def uniform_selection(n, d, coef_max):
    """
    Generate a binomial using uniform selection.
    @param n: The number of variables in the binomial.
    @param d: The maximum degree of the binomial.
    @param coef_max: The maximum coefficient for both monomials in the binomial.
    @return: Generated binomial represented as a polynomial object.
    """

    chains = list(all_monomials_up_to(n, d))
    polynomial = []
    for _ in range(2):
        coef = random.randint(1, coef_max)
        monomial = np.array(random.choice(chains))
        monomial = np.concatenate([np.array([coef]), monomial])
        polynomial.append(monomial)

    return poly.Polynomial(np.array(polynomial))


def all_monomials(n, d):
    """
    List all monomials (not including the coefficient) in n variables of degree d.
    @param n: The number of variables.
    @param d: The degree of the monomial.
    @return: An iterator object representing each monomial as a tuple of exponents.
    """

    if n == 1:
        yield (d, )
    else:
        for value in range(d + 1):
            for diff in all_monomials(n - 1, d - value):
                yield (value, ) + diff


def all_monomials_up_to(n, d):
    """
    List all monomials (not including the coefficient) in n variables of degree up to d.
    @param n: The number of variables.
    @param d: The maximum degree of the monomial.
    @return: An iterator object representing each monomial as a tuple of exponents.
    """

    if d == 1:
        return all_monomials(n, d)
    else:
        chains = all_monomials(n, 1)
        for i in range(2, d + 1):
            new_chain = all_monomials(n, i)
            chains = chain(chains, new_chain)
        return chains
