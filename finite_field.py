import numpy as np


# The choice of p
p = 23

def ff_add(a, b):
    """
    Compute the sum of two numbers in the finite field Z/pZ.
    @param a: One number.
    @param b: Another number.
    @return: The sum.
    """

    return (a + b) % p

def ff_mul(a, b):
    """
    Compute the product of two numbers in the finite field Z/pZ.
    @param a: One number.
    @param b: Another number.
    @return: The product.
    """

    return (a * b) % p

def ff_div(a, b):
    """
    Compute the ratio of two numbers in the finite field Z/pZ.
    @param a: One number.
    @param b: Another number.
    @return: The ratio.
    """

    b_inv = ff_inv(b)

    return ff_mul(a, b_inv)

def ff_inv(a):
    """
    Compute the inverse of a number in the finite field Z/pZ.
    @param a: A number.
    @return: The inverse.
    """

    t, r, new_t, new_r = 0, p, 1, a % p
    while new_r != 0:
        q = r // new_r
        t, new_t = new_t, t - q * new_t
        r, new_r = new_r, r - q * new_r

    if t < 0:
        t = t % p

    return t
