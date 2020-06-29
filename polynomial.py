import numpy as np

class Polynomial:

    def __init__(self, coef, monomials):
        assert len(coef) > 0, 'There should be at least 1 term in the polynomial'
        assert len(coef) == len(monomials), 'Number of coefficients should be the same as number of monomials'
        assert len(set([len(x) for x in monomials])) == 1, 'All monomials should have the same number of variables, include 0 if needed'

        self.coef = coef
        self.monomials = monomials
        self.nvar = len(monomials[0])
        self.nterm = len(monomials)
        self.degree = np.array(monomials).sum(axis = 1).max()

    def __repr__(self):
        result = []
        for i in range(self.nterm):
            exponents = self.monomials[i]
            monomial = []
            for j, k in zip(np.arange(1, self.nvar + 1), exponents):
                if k != 0:
                    monomial.append('x{}^{}'.format(j, k))
            result.append('*'.join(monomial))
        return ' + '.join(result)

    def lt(self):
