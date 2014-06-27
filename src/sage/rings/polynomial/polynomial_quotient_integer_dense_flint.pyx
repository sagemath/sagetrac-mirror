r"""
Dense univariate polynomials in `\\ZZ[x]/<f>` where `f` is irreducible, implemented using FLINT.

EXAMPLES::

    sage: P.<x> = ZZ['x'].quotient(x^8 + 1)
    sage: f = P.random_element()
    sage: g = P.random_element()
    sage: R = CyclotomicField(16).ring_of_integers()
    sage: R(f) * R(g) == R(f*g)
    True

    sage: n(vector(f).norm() * vector(g).norm() / vector(f*g).norm())
    1.00...

    sage: P(vector(range(8)))
    7*x^7 + 6*x^6 + 5*x^5 + 4*x^4 + 3*x^3 + 2*x^2 + x

    sage: g.matrix()
    [-95  -1  -2 -12   0   0   1  -1]
    [  1 -95  -1  -2 -12   0   0   1]
    [ -1   1 -95  -1  -2 -12   0   0]
    [  0  -1   1 -95  -1  -2 -12   0]
    [  0   0  -1   1 -95  -1  -2 -12]
    [ 12   0   0  -1   1 -95  -1  -2]
    [  2  12   0   0  -1   1 -95  -1]
    [  1   2  12   0   0  -1   1 -95]
    
"""

# We need to define this stuff before including the templating stuff
# to make sure the function get_cparent is found since it is used in
# 'polynomial_template.pxi'.

cdef inline cparent get_cparent(parent) except? NULL:
    try:
        return (<Polynomial_integer_dense_flint>parent.modulus()).__poly
    except AttributeError:
        return NULL


include "sage/libs/flint/fmpz_poly_quo_linkage.pxi"        
include "polynomial_template.pxi"

cdef class PolynomialQuotientRingElement_integer_flint(Polynomial_template):
    def __init__(self, parent, x=None, check=True, is_gen=False, construct=False):
        """
        Create a new quotient ring element.

        EXAMPLE::

            sage: P.<x> = ZZ[]
            sage: R.<X> = P.quotient(x^8 + 1)
            sage: R(x^10)
            -X^2

            sage: R(range(8))
            7*X^7 + 6*X^6 + 5*X^5 + 4*X^4 + 3*X^3 + 2*X^2 + X

        """
        if isinstance(x, Polynomial_integer_dense_flint):
            Polynomial_template.__init__(self, parent, 0, check, is_gen, construct)
            self._set_fmpz_poly((<Polynomial_integer_dense_flint>x).__poly)
            return
            
        try:
            _ = iter(x) # trigger TypeError early if we can't iterate
            Polynomial_template.__init__(self, parent, 0, check, is_gen, construct)
            self._set_list(list(x))
        except TypeError:
            pass
            
        Polynomial_template.__init__(self, parent, x, check, is_gen, construct)

    cdef Polynomial_template _new(self):
        cdef Polynomial_template e = <Polynomial_template>PY_NEW(self.__class__)
        fmpz_poly_init(&e.x)
        e._parent = self._parent
        e._cparent = self._cparent
        return e

    cpdef Polynomial _new_constant_poly(self, x, Parent P):
        """
        Construct a new constant polynomial.

        INPUT:

        - ``x`` - a scalar.
        - ``P`` - a parent.
        
        EXAMPLE::

            sage: P.<x> = ZZ[]
            sage: Q.<x> = P.quotient(P.cyclotomic_polynomial(7))
            sage: x._new_constant_poly(2, Q)
            2

            sage: x._new_constant_poly(2, P)
            Traceback (most recent call last):
            ...
            TypeError: Supplied parent does not match self's parent.

        """
        if type(P) != type(self._parent):
          raise TypeError("Supplied parent does not match self's parent.")
        
        cdef Polynomial_template r = <Polynomial_template>PY_NEW(self.__class__)
        r._parent = P
        r._cparent = get_cparent(P)
        fmpz_poly_init(&r.x)
        cdef Integer c = Integer(x)
        fmpz_poly_set_coeff_mpz(&r.x, 0, c.value)
        return r

    cdef int _set_list(self, list x) except -1:
        cdef Integer tmp
        cdef fmpz_poly_struct *modulus = self._cparent
        
        if not len(x):
            fmpz_poly_zero(&self.x)
            return 0

        sig_on()
        fmpz_poly_realloc(&self.x, fmpz_poly_degree(modulus))
        sig_off()

        sig_on()
        for i from 0 <= i < len(x):
            tmp = Integer(x[i])
            fmpz_poly_set_coeff_mpz(&self.x, i, tmp.value)
        sig_off()
        return 0

    cdef int _set_fmpz_poly(self, fmpz_poly_t x) except -1:
        sig_on()
        celement_set(&self.x, x, self._cparent)
        sig_off()
        return 0

    def __getitem__(self, i):
        """
        EXAMPLE::

            sage: P.<x> = ZZ[]
            sage: Q.<x> = P.quotient(P.cyclotomic_polynomial(7))
            sage: f = - x^4 + x^2 - 3*x + 4
            sage: f[0]
            4
            sage: f[1]
            -3
            sage: f[5]
            0
            sage: f[6]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1 of modulus.

            sage: f[-1]
            Traceback (most recent call last):
            ...
            IndexError: index must be between 0 and degree minus 1 of modulus.

        """
        cdef Integer c = PY_NEW(Integer)        
        cdef Polynomial_template r
        if 0 <= i < celement_len(&self.x, (<Polynomial_template>self)._cparent):
            fmpz_poly_get_coeff_mpz(c.value, &self.x, i)
            return c
        elif 0 <= i < self._parent.modulus().degree():
            return c
        else:
            raise IndexError("index must be between 0 and degree minus 1 of modulus.")

    def __reduce__(self):
        """
        EXAMPLE::

            sage: P.<x> = ZZ[]
            sage: Q.<x> = P.quotient(x^4+1)
            sage: loads(dumps(x)) == x
            True
        """
        return make_element, ((<Polynomial_template>self)._parent, (self.list(),))

    def matrix(self):
        """
        Return the matrix of right multiplication by `1, x, x^2, …, x^{d-1}`
        where `d` is the degree of the modulus. Thus the *rows* of this matrix
        give the images of each of the `x^i`.

        EXAMPLE::

            sage: P.<x> = ZZ[]
            sage: Q.<x> = P.quotient(P.cyclotomic_polynomial(7))
            sage: f = -22*x^5 + 378*x^4 + 3*x^2 - x - 4
            sage: f.matrix()
            [  -4   -1    3    0  378  -22]
            [  22   18   21   25   22  400]
            [-400 -378 -382 -379 -375 -378]
            [ 378  -22    0   -4   -1    3]
            [  -3  375  -25   -3   -7   -4]
            [   4    1  379  -21    1   -3]        
        """
        from sage.rings.integer_ring import ZZ
        from sage.matrix.constructor import matrix
        
        x = self._parent.gen()
        degree = self._parent.modulus().degree()
        M = matrix(ZZ, degree, degree)
        for r in xrange(degree):
            row = x**r * self
            for c in xrange(row.degree()+1):
                M[r,c] = row[c]
        return M
            