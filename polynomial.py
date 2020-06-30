import numpy as np
import math

class Polynomial:
    """
    The polynomial class defines the polynomials to be used in Buchberger's algorithm.
    """

    def __init__(self, monomials):
        """
        The constructor..
        @param monomials: Array of monomials that make up the polynomial.
        each monomial is represented as an array of exponents on variables including the coefficient.
        """

        assert len(monomials) > 0, 'There should be at least 1 monomials in the polynomial'
        assert len(set([len(x) for x in monomials])) == 1, 'All monomials should have the same number of variables, include 0 if needed'

        self.monomials = monomials
        self.nvar = len(monomials[0]) - 1
        self.nterm = len(monomials)
        self.degree = monomials[:, 1:].sum(axis = 1).max()

    def __repr__(self):
        """
        Print the polynomial in usual form in variables x1, x2, ...
        """

        result = []
        for i in range(self.nterm):
            coef = self.monomials[i][0]
            exponents = self.monomials[i][1:]
            monomial = []
            for j, k in zip(np.arange(1, self.nvar + 1), exponents):
                if k != 0:
                    monomial.append('x{}^{}'.format(j, k))
            result.append(str(coef) + '*' + '*'.join(monomial))
        return ' + '.join(result)

    def lt(self):
        """
        Get the leading term in the polynomial using the grevlex order.
        @return: leading term represented as an array of exponents.
        """

        arr = self.monomials.copy()
        mergesort(arr)
        return arr[::-1][0]

def compare(monomial1, monomial2):
    """
    Compare two monomials using the grevlex order.
    @param monomial1: The 1st monomial for comparison, represented as an array of exponents.
    @param monomial2: The 2nd monomial for comparison, represented as an array of exponents.
    @return: 1 if monomial1 is larger, -1 if monomial2 is larger, 0 if they are the same.
    """

    assert len(monomial1) == len(monomial2), 'The two monomials should have the same number variables'

    # Compare total degree of monomials first
    monomial1_degree = sum(monomial1[1:])
    monomial2_degree = sum(monomial2[1:])
    if monomial1_degree > monomial2_degree:
        return 1
    if monomial1_degree < monomial2_degree:
        return -1

    # Break the tie using reverse lexicographic order
    for i in range(len(monomial1))[1:][::-1]:
        if monomial1[i] > monomial2[i]:
            return -1
        if monomial1[i] < monomial2[i]:
            return 1

    return 0

def mergesort(arr):
    """
    A modified mergesort algorithm to sort the monomials of a polynomial in descending grevlex order.
    @param arr: Array of monomials that make up the polynomial.
    @return: sorted array of monomials.
    """
    if len(arr) > 1:
        middle = math.floor(len(arr) / 2)
        left = arr[:middle]
        right = arr[middle:]

        mergesort(left)
        mergesort(right)
        merge(left, right, arr)

def merge(left, right, arr):
    """
    The merge step in the modified mergesort algorithm.
    @param left: The left array in the mergesort algorithm.
    @param right: The right array in the mergesort algorithm.
    @param arr: The original array to be modified during the merge step.
    @return: merged array.
    """

    # Append the 'infinity' array to both left and right to assist in the merge step
    inf_arr_left = np.array([float('inf')] * left.shape[1])
    inf_arr_right = np.array([float('inf')] * right.shape[1])
    left = np.vstack([left, inf_arr_left])
    right = np.vstack([right, inf_arr_right])

    left_ix = 0
    right_ix = 0
    for ix in range(len(arr)):
        if (compare(left[left_ix], right[right_ix]) == -1):
            arr[ix] = left[left_ix]
            left_ix += 1
        else:
            arr[ix] = right[right_ix]
            right_ix += 1
