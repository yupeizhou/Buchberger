import numpy as np
import finite_field as ff
import math


class Polynomial:
    """
    The polynomial class defines the polynomials to be used in Buchberger's algorithm.
    """

    def __init__(self, monomials):
        """
        The constructor. The polynomials will be grouped by terms and sorted in grevlex order.
        @param monomials: Array of monomials that make up the polynomial.
        Each monomial is represented as an array of exponents on variables including the coefficient.
        """

        assert len(monomials) > 0, 'There should be at least 1 monomials in the polynomial'
        assert len(set([len(x) for x in monomials])) == 1, 'All monomials should have the same number of variables, include 0 if needed'

        # Sort the polynomial and group the terms
        arr = monomials.copy()
        mergesort(arr)
        arr = arr[::-1]
        arr_ix = 0
        grouped_monomials = []
        monomial = np.array([])

        while arr_ix < len(arr) - 1:
            if compare(arr[arr_ix], arr[arr_ix + 1]) != 0:
                if len(monomial) == 0:
                    grouped_monomials.append(arr[arr_ix])
                    arr_ix += 1
                else:
                    grouped_monomials.append(monomial)
                    arr_ix += 1
                    monomial = np.array([])
            else:
                if len(monomial) == 0:
                    monomial = np.concatenate([np.array([ff.ff_add(arr[arr_ix][0], arr[arr_ix + 1][0])]), arr[arr_ix][1:]])
                    arr_ix += 1
                else:
                    monomial = np.concatenate([np.array([ff.ff_add(monomial[0] + arr[arr_ix + 1][0])]), arr[arr_ix][1:]])
                    arr_ix += 1
        if len(monomial) == 0:
            grouped_monomials.append(arr[arr_ix])
        else:
            grouped_monomials.append(monomial)

        # Only keep nonzero monomials
        grouped_monomials = np.array(grouped_monomials)
        filtered_monomials = []
        for monomial in grouped_monomials:
            if (monomial[0] % ff.p != 0):
                filtered_monomials.append(np.concatenate([np.array([monomial[0] % ff.p]), monomial[1:]]))
        grouped_monomials = np.array(filtered_monomials)

        # Check if the resulting polynomial is a zero polynomial
        if len(grouped_monomials) > 0:
            self.monomials = grouped_monomials
            self.nvar = len(grouped_monomials[0]) - 1
            self.nterm = len(grouped_monomials)
            self.degree = grouped_monomials[:, 1:].sum(axis = 1).max()
        else:
            self.nvar = len(monomials[0]) - 1
            self.nterm = 0
            self.degree = 0
            self.monomials = np.array([np.zeros(self.nvar + 1)])


    def is_zero(self):
        """
        Check if the polynomial is a zero polynomial.
        @return: True if the polynomial is a zero polynomial, False if not.
        """
        if self.nterm == 0:
            return True
        else:
            return False


    def __repr__(self):
        """
        Print the polynomial in usual form in variables x1, x2, ...
        @return: Polynomial printed in string.
        """

        if not self.is_zero():
            result = []
            for i in range(self.nterm):
                coef = self.monomials[i][0]
                exponents = self.monomials[i][1:]
                monomial = [str(coef)]
                for j, k in zip(np.arange(1, self.nvar + 1), exponents):
                    if k != 0:
                        monomial.append('x{}^{}'.format(j, k))
                result.append('*'.join(monomial))
            return ' + '.join(result)
        # If the polynomial is a zero polynomial, print 0
        else:
            return '0'


    def lt(self):
        """
        Get the leading term in the polynomial under the grevlex order.
        @return: Leading term represented as a polynomial object.
        """

        leading_term = self.monomials[0]
        return Polynomial(leading_term.reshape(-1, len(leading_term)))


    def add(self, poly):
        """
        Add one polynomial to another polynomial.
        @param poly: Another polynomial to add to.
        @return: Sum of the two polynomials represented as a polynomial object.
        """

        assert isinstance(poly, Polynomial), 'Can only add to a polynomial'

        poly_sum = np.concatenate([self.monomials, poly.monomials])
        return Polynomial(poly_sum)


    def subtract(self, poly):
        """
        Add one polynomial from another polynomial.
        @param poly: Another polynomial subtract.
        @return: Difference of the two polynomials represented as a polynomial object.
        """

        assert isinstance(poly, Polynomial), 'Can only subtract a polynomial'

        return self.add(poly.scalar_multiply(-1))


    def multiply(self, poly):
        """
        Multiply one polynomial with another polynomial.
        @param poly: Another polynomial to multiply with.
        @return: Product of the two polynomials represented as a polynomial object.
        """

        assert isinstance(poly, Polynomial), 'Can only multiply with polynomial'

        poly_product = []
        for monomial1 in self.monomials:
            for monomial2 in poly.monomials:
                product = np.concatenate([np.array([ff.ff_mul(monomial1[0], monomial2[0])]), monomial1[1:] + monomial2[1:]])
                poly_product.append(product)
        poly_product = np.array(poly_product)

        return Polynomial(poly_product)


    def scalar_multiply(self, scalar):
        """
        Multiply one polynomial with a scalar.
        @param scalar: The scalar to multiply with the polynomial. Can be integer or float.
        @return: Product of the scalar with the polynomial represented as a polynomial object.
        """

        assert isinstance(scalar, (int, float)), 'Scalar multiplication only.'

        poly_scalar_prod = []
        for monomial in self.monomials:
            product = np.concatenate([np.array([monomial[0] * scalar]), monomial[1:]])
            poly_scalar_prod.append(product)
        poly_scalar_prod = np.array(poly_scalar_prod)

        return Polynomial(poly_scalar_prod)


    def divide(self, monomial):
        """
        Return the ratio of the monomial to another monomial.
        @param monomial: Another monomial being the divisor.
        @return: The ratio of the two monomial represented as a polynomial object if one divides the other,
        return false if one does not divide the other.
        """

        assert (self.nterm == 1) or (self.nterm == 0), 'This only works for monomials.'
        assert monomial.nterm == 1, 'This only works for monomials.'

        if self.nterm == 0:
            return Polynomial(np.array([np.zeros(self.nvar + 1)]))

        if (self.monomials[0][1:] >= monomial.monomials[0][1:]).all():
            ratio = np.concatenate([np.array([ff.ff_div(self.monomials[0][0], monomial.monomials[0][0])]),
            self.monomials[0][1:] - monomial.monomials[0][1:]])
            return Polynomial(np.array([ratio]))
        else:
            return False


    def lcm(self, monomial):
        """
        Return the least common multiple of two monomials.
        @param monomial: Another monomial to take lcm with.
        @return: The least common multiple represented as a polynomial object.
        """

        assert self.nterm == 1, 'This only works for monomials.'
        assert monomial.nterm == 1, 'This only works for monomials.'

        return Polynomial(np.maximum(self.monomials, monomial.monomials))


def compare(monomial1, monomial2):
    """
    Compare two monomials under the grevlex order, disregarding the coefficient.
    @param monomial1: The 1st monomial for comparison, represented as an array of exponents.
    @param monomial2: The 2nd monomial for comparison, represented as an array of exponents.
    @return: 1 if monomial1 is larger, -1 if monomial2 is larger, 0 if they are the same.
    """

    assert len(monomial1) == len(monomial2), 'The two monomials should have the same number of variables'

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
    @return: Sorted array of monomials.
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
    @return: Merged array of monomials.
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
