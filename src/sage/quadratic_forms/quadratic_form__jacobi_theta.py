"""
Jacobi Theta Series of Quadratic Forms

AUTHOR:

Andrew Fiori
"""
from sage.modular.jacobi.jacobi_theta import __jacobi_theta


def bilinear_form(self, x, y):
    """
    Return the associated bilinear form of a quadratic form.

    EXAMPLES:
    """
    return self(x + y) - self(x) - self(y)


def jacobi_theta(self, vec, max=10):
    """
    Compute the theta series as a power series in the variable given
    in ``var_str`` (which defaults to `q`), up to the specified precision
    `O(q^max)`.

    EXAMPLES:
    """
    return __jacobi_theta(self, vec, max)
