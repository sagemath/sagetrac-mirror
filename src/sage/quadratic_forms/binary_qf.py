r"""
Binary Quadratic Forms with Integer Coefficients

This module provides a specialized class for working with a binary quadratic
form `a x^2 + b x y + c y^2`, stored as a triple of integers `(a, b, c)`.

EXAMPLES::

    sage: Q = BinaryQF([1,2,3])
    sage: Q
    x^2 + 2*x*y  + 3*y^2
    sage: Q.discriminant()
    -8
    sage: Q.reduced_form()
    x^2 + 2*y^2
    sage: Q(1, 1)
    6

TESTS::

    sage: Q == loads(dumps(Q))
    True

AUTHORS:

- Jon Hanke (2006-08-08):

  - Appended to add the methods :func:`BinaryQF_reduced_representatives`,
    :meth:`~BinaryQF.is_reduced`, and ``__add__`` on 8-3-2006 for Coding Sprint
    #2.
  - Added Documentation and :meth:`~BinaryQF.complex_point` method on 8-8-2006.

- Nick Alexander: add doctests and clean code for Doc Days 2
- William Stein (2009-08-05): composition; some ReSTification.
- William Stein (2009-09-18): make immutable.
- Jon Hanke (2011-02-05):  Add spectral decomp. and some nice pictures.
"""

#*****************************************************************************
#       Copyright (C) 2006-2009 William Stein and Jon Hanke
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.pari.all import pari
from sage.rings.all import ZZ, is_fundamental_discriminant
from sage.arith.all import divisors, gcd
from sage.structure.sage_object import SageObject
from sage.misc.cachefunc import cached_method

from sage.matrix.constructor import Matrix
from sage.modules.free_module_element import vector
from sage.functions.other import sqrt, ceil
from sage.functions.trig import arccos
from sage.plot.all import Graphics, point, ellipse, text


class BinaryQF(SageObject):
    """
    A binary quadratic form over `\ZZ`.

    INPUT:

    - `v` -- a list or tuple of 3 entries:  [a,b,c], or a quadratic homogeneous
      polynomial in two variables with integer coefficients


    OUTPUT:

    the binary quadratic form a*x^2 + b*x*y + c*y^2.

    EXAMPLES::

        sage: b = BinaryQF([1,2,3])
        sage: b.discriminant()
        -8
        sage: R.<x, y> = ZZ[]
        sage: BinaryQF(x^2 + 2*x*y + 3*y^2) == b
        True
    """
    # Initializes the form with a 3-element list
    def __init__(self, abc):
        r"""
        Creates the binary quadratic form `ax^2 + bxy + cy^2` from the
        triple [a,b,c] over `\ZZ` or from a polynomial.

        INPUT:

        - ``abc`` -- 3-tuple of integers, or a quadratic homogeneous polynomial
          in two variables with integer coefficients

        EXAMPLES::

            sage: Q = BinaryQF([1,2,3]); Q
            x^2 + 2*x*y + 3*y^2
            sage: Q = BinaryQF([1,2])
            Traceback (most recent call last):
            ...
            TypeError: Binary quadratic form must be given by a list of three coefficients

            sage: R.<x, y> = ZZ[]
            sage: f = x^2 + 2*x*y + 3*y^2
            sage: BinaryQF(f)
            x^2 + 2*x*y + 3*y^2
            sage: BinaryQF(f + x)
            Traceback (most recent call last):
            ...
            TypeError: Binary quadratic form must be given by a quadratic homogeneous bivariate integer polynomial

        TESTS::

            sage: BinaryQF(0)
            0
        """
        if isinstance(abc, (list, tuple)):
            if len(abc) != 3:
                # Check we have three coefficients
                raise TypeError("Binary quadratic form must be given by a list of three coefficients")
            self._a, self._b, self._c = [ZZ(x) for x in abc]
        else:
            f = abc
            from sage.rings.polynomial.multi_polynomial_element import is_MPolynomial
            if f.is_zero():
                self._a, self._b, self._c = [ZZ(0), ZZ(0), ZZ(0)]
            elif (is_MPolynomial(f) and f.is_homogeneous() and f.base_ring() == ZZ
                    and f.degree() == 2 and f.parent().ngens() == 2):
                x, y = f.parent().gens()
                self._a, self._b, self._c = [f.monomial_coefficient(mon) for mon in [x**2, x*y, y**2]]
            else:
                raise TypeError("Binary quadratic form must be given by a quadratic homogeneous bivariate integer polynomial")

    def _pari_init_(self):
        """
        Used to convert this quadratic form to Pari.

        EXAMPLES::

            sage: f = BinaryQF([2,3,4]); f
            2*x^2 + 3*x*y + 4*y^2
            sage: f._pari_init_()
            'Qfb(2,3,4)'
            sage: pari(f)
            Qfb(2, 3, 4)
            sage: type(pari(f))
            <type 'sage.libs.pari.gen.gen'>
            sage: gp(f)
            Qfb(2, 3, 4)
            sage: type(gp(f))
            <class 'sage.interfaces.gp.GpElement'>
        """
        return 'Qfb(%s,%s,%s)'%(self._a,self._b,self._c)

    def __mul__(self, right):
        """
        Gauss composition of binary quadratic forms.  The result is
        not reduced.

        EXAMPLES:

        We explicitly compute in the group of classes of positive
        definite binary quadratic forms of discriminant -23.

        ::

            sage: R = BinaryQF_reduced_representatives(-23); R
            [x^2 + x*y + 6*y^2, 2*x^2 - x*y + 3*y^2, 2*x^2 + x*y + 3*y^2]
            sage: R[0] * R[0]
            x^2 + x*y + 6*y^2
            sage: R[1] * R[1]
            4*x^2 + 3*x*y + 2*y^2
            sage: (R[1] * R[1]).reduced_form()
            2*x^2 + x*y + 3*y^2
            sage: (R[1] * R[1] * R[1]).reduced_form()
            x^2 + x*y + 6*y^2

        """
        if not isinstance(right, BinaryQF):
            raise TypeError("both self and right must be binary quadratic forms")
        # There could be more elegant ways, but qfbcompraw isn't
        # wrapped yet in the PARI C library.  We may as well settle
        # for the below, until somebody simply implements composition
        # from scratch in Cython.
        v = list(pari('qfbcompraw(%s,%s)'%(self._pari_init_(),
                                           right._pari_init_())))
        return BinaryQF(v)

    def __getitem__(self, n):
        """
        Return the n-th component of this quadratic form.

        If this form is `a x^2 + b x y + c y^2`, the 0-th component is `a`,
        the 1-st component is `b`, and `2`-nd component is `c`.

        Indexing is like lists -- negative indices and slices are allowed.

        EXAMPLES::


            sage: Q = BinaryQF([2,3,4])
            sage: Q[0]
            2
            sage: Q[2]
            4
            sage: Q[:2]
            (2, 3)
            sage: tuple(Q)
            (2, 3, 4)
            sage: list(Q)
            [2, 3, 4]
        """
        return (self._a, self._b, self._c)[n]

    def __call__(self, *args):
        r"""
        Evaluate this quadratic form at a point.

        INPUT:

        - args -- x and y values, as a pair x, y or a list, tuple, or
          vector

        EXAMPLES::


            sage: Q = BinaryQF([2, 3, 4])
            sage: Q(1, 2)
            24

        TESTS::

            sage: Q = BinaryQF([2, 3, 4])
            sage: Q([1, 2])
            24
            sage: Q((1, 2))
            24
            sage: Q(vector([1, 2]))
            24
        """
        if len(args) == 1:
            args = args[0]
        x, y = args
        return (self._a * x + self._b * y) * x + self._c * y**2

    def __hash__(self):
        r"""
        TESTS::

            sage: hash(BinaryQF([2,2,3]))
            802
            sage: hash(BinaryQF([2,3,2]))
            562
            sage: hash(BinaryQF([3,2,2]))
            547
        """
        return hash(self._a) ^ (hash(self._b) << 4) ^ (hash(self._c) << 8)

    def __cmp__(self, right):
        """
        Returns True if self and right are identical: the same coefficients.

        EXAMPLES::


            sage: P = BinaryQF([2,2,3])
            sage: Q = BinaryQF([2,2,3])
            sage: R = BinaryQF([1,2,3])
            sage: P == Q # indirect doctest
            True
            sage: P == R # indirect doctest
            False

        TESTS::

            sage: P == P
            True
            sage: Q == P
            True
            sage: R == P
            False
            sage: P == 2
            False
        """
        if not isinstance(right, BinaryQF):
            return cmp(type(self), type(right))
        return cmp((self._a,self._b,self._c), (right._a,right._b,right._c))

    def __add__(self, Q):
        """
        Returns the component-wise sum of two forms.

        That is, given `a_1 x^2 + b_1 x y + c_1 y^2` and `a_2 x^2 + b_2 x y +
        c_2 y^2`, returns the form
        `(a_1 + a_2) x^2 + (b_1 + b_2) x y + (c_1 + c_2) y^2.`

        EXAMPLES::


            sage: P = BinaryQF([2,2,3]); P
            2*x^2 + 2*x*y + 3*y^2
            sage: Q = BinaryQF([-1,2,2]); Q
            -x^2 + 2*x*y + 2*y^2
            sage: P + Q
            x^2 + 4*x*y + 5*y^2
            sage: P + Q == BinaryQF([1,4,5]) # indirect doctest
            True

        TESTS::

            sage: Q + P == BinaryQF([1,4,5]) # indirect doctest
            True
        """
        return BinaryQF([self._a + Q._a, self._b + Q._b, self._c + Q._c])

    def __sub__(self, Q):
        """
        Returns the component-wise difference of two forms.

        That is, given `a_1 x^2 + b_1 x y + c_1 y^2` and `a_2 x^2 + b_2 x y +
        c_2 y^2`, returns the form
        `(a_1 - a_2) x^2 + (b_1 - b_2) x y + (c_1 - c_2) y^2.`

        EXAMPLES::


            sage: P = BinaryQF([2,2,3]); P
            2*x^2 + 2*x*y + 3*y^2
            sage: Q = BinaryQF([-1,2,2]); Q
            -x^2 + 2*x*y + 2*y^2
            sage: P - Q
            3*x^2 + y^2
            sage: P - Q == BinaryQF([3,0,1]) # indirect doctest
            True

        TESTS::

            sage: Q - P == BinaryQF([3,0,1]) # indirect doctest
            False
            sage: Q - P != BinaryQF([3,0,1]) # indirect doctest
            True
        """
        return BinaryQF([self._a - Q._a, self._b - Q._b, self._c - Q._c])

    def _repr_(self):
        """
        Display the quadratic form.

        EXAMPLES::


            sage: Q = BinaryQF([1,2,3]); Q # indirect doctest
            x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([-1,2,3]); Q
            -x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([0,0,0]); Q
            0
        """
        return repr(self.polynomial())

    def _latex_(self):
        """
        Return latex representation of this binary quadratic form.

        EXAMPLES::

            sage: f = BinaryQF((778,1115,400)); f
            778*x^2 + 1115*x*y + 400*y^2
            sage: latex(f) # indirect doctest
            778 x^{2} + 1115 x y + 400 y^{2}
        """
        return self.polynomial()._latex_()

    @cached_method
    def polynomial(self):
        """
        Returns the binary quadratic form as a homogeneous 2-variable
        polynomial.

        EXAMPLES::

            sage: Q = BinaryQF([1,2,3])
            sage: Q.polynomial()
            x^2 + 2*x*y + 3*y^2

            sage: Q = BinaryQF([-1,-2,3])
            sage: Q.polynomial()
            -x^2 - 2*x*y + 3*y^2

            sage: Q = BinaryQF([0,0,0])
            sage: Q.polynomial()
            0
        """
        return self(ZZ['x, y'].gens())

    @cached_method
    def discriminant(self):
        """
        Returns the discriminant `b^2 - 4ac` of the binary
        form `ax^2 + bxy + cy^2`.

        EXAMPLES::


            sage: Q = BinaryQF([1,2,3])
            sage: Q.discriminant()
            -8
        """
        return self._b**2 - 4 * self._a * self._c

    @cached_method
    def has_fundamental_discriminant(self):
        """
        Checks if the discriminant D of this form is a fundamental
        discriminant (i.e. D is the smallest element of its
        squareclass with D = 0 or 1 mod 4).

        EXAMPLES::

            sage: Q = BinaryQF([1,0,1])
            sage: Q.discriminant()
            -4
            sage: Q.has_fundamental_discriminant()
            True

            sage: Q = BinaryQF([2,0,2])
            sage: Q.discriminant()
            -16
            sage: Q.has_fundamental_discriminant()
            False
        """
        return is_fundamental_discriminant(self.discriminant())

    @cached_method
    def is_primitive(self):
        """
        Checks if the form `ax^2 + bxy + cy^2`  satisfies
        `\gcd(a,b,c)=1`, i.e., is primitive.

        EXAMPLES::


            sage: Q = BinaryQF([6,3,9])
            sage: Q.is_primitive()
            False

            sage: Q = BinaryQF([1,1,1])
            sage: Q.is_primitive()
            True

            sage: Q = BinaryQF([2,2,2])
            sage: Q.is_primitive()
            False

            sage: rqf = BinaryQF_reduced_representatives(-23*9)
            sage: [qf.is_primitive() for qf in rqf]
            [True, True, True, False, True, True, False, False, True]
            sage: rqf
            [x^2 + x*y + 52*y^2,
            2*x^2 - x*y + 26*y^2,
            2*x^2 + x*y + 26*y^2,
            3*x^2 + 3*x*y + 18*y^2,
            4*x^2 - x*y + 13*y^2,
            4*x^2 + x*y + 13*y^2,
            6*x^2 - 3*x*y + 9*y^2,
            6*x^2 + 3*x*y + 9*y^2,
            8*x^2 + 7*x*y + 8*y^2]
            sage: [qf for qf in rqf if qf.is_primitive()]
            [x^2 + x*y + 52*y^2,
            2*x^2 - x*y + 26*y^2,
            2*x^2 + x*y + 26*y^2,
            4*x^2 - x*y + 13*y^2,
            4*x^2 + x*y + 13*y^2,
            8*x^2 + 7*x*y + 8*y^2]
        """
        return gcd([self._a, self._b, self._c])==1

    @cached_method
    def is_weakly_reduced(self):
        """
        Checks if the form `ax^2 + bxy + cy^2`  satisfies
        `|b| \leq a \leq c`, i.e., is weakly reduced.

        EXAMPLES::


            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_weakly_reduced()
            False

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_weakly_reduced()
            True

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_weakly_reduced()
            True
        """
        if self.discriminant() >= 0:
            raise NotImplementedError("only implemented for negative discriminants")
        return (abs(self._b) <= self._a) and (self._a <= self._c)

    @cached_method
    def reduced_form(self):
        """
        Return the unique reduced form equivalent to ``self``. See also
        :meth:`~is_reduced`.

        EXAMPLES::

            sage: a = BinaryQF([33,11,5])
            sage: a.is_reduced()
            False
            sage: b = a.reduced_form(); b
            5*x^2 - x*y + 27*y^2
            sage: b.is_reduced()
            True

            sage: a = BinaryQF([15,0,15])
            sage: a.is_reduced()
            True
            sage: b = a.reduced_form(); b
            15*x^2 + 15*y^2
            sage: b.is_reduced()
            True
        """
        if self.discriminant() >= 0 or self._a < 0:
            raise NotImplementedError("only implemented for positive definite forms")
        if not self.is_reduced():
            v = list(pari('Vec(qfbred(Qfb(%s,%s,%s)))'%(self._a,self._b,self._c)))
            return BinaryQF(v)
        else:
            return self

    def is_equivalent(self, right):
        """
        Return true if self and right are equivalent, i.e., have the
        same reduced form.

        INPUT:

        - ``right`` -- a binary quadratic form

        EXAMPLES::

            sage: a = BinaryQF([33,11,5])
            sage: b = a.reduced_form(); b
            5*x^2 - x*y + 27*y^2
            sage: a.is_equivalent(b)
            True
            sage: a.is_equivalent(BinaryQF((3,4,5)))
            False
        """
        if not isinstance(right, BinaryQF):
            raise TypeError("right must be a binary quadratic form")
        return self.reduced_form() == right.reduced_form()

    @cached_method
    def is_reduced(self):
        """
        Checks if the quadratic form is reduced, i.e., if the form
        `ax^2 + bxy + cy^2` satisfies `|b|\leq a \leq c`, and
        that `b\geq 0` if either `a = b` or `a = c`.

        EXAMPLES::

            sage: Q = BinaryQF([1,2,3])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([2,1,3])
            sage: Q.is_reduced()
            True

            sage: Q = BinaryQF([1,-1,1])
            sage: Q.is_reduced()
            False

            sage: Q = BinaryQF([1,1,1])
            sage: Q.is_reduced()
            True
        """
        return (-self._a < self._b <= self._a < self._c) or \
               (ZZ(0) <= self._b <= self._a == self._c)

    def complex_point(self):
        r"""
        Returns the point in the complex upper half-plane associated
        to this (positive definite) quadratic form.

        For positive definite forms with negative discriminants, this is a
        root `\tau` of `a x^2 + b x + c` with the imaginary part of `\tau`
        greater than 0.

        EXAMPLES::

            sage: Q = BinaryQF([1,0,1])
            sage: Q.complex_point()
            1.00000000000000*I
        """
        if self.discriminant() >= 0:
            raise NotImplementedError("only implemented for negative discriminant")
        R = ZZ['x']
        x = R.gen()
        Q1 = R(self.polynomial()(x,1))
        return [z  for z in Q1.complex_roots()  if z.imag() > 0][0]

    def matrix_action_left(self, M):
        r"""
        Return the binary quadratic form resulting from the left action
        of the 2-by-2 matrix ``M`` on the quadratic form ``self``.

        Here the action of the matrix `M = \begin{pmatrix} a & b \\ c & d
        \end{pmatrix}` on the form `Q(x, y)` produces the form `Q(ax+cy,
        bx+dy)`.

        EXAMPLES::

            sage: Q = BinaryQF([2, 1, 3]); Q
            2*x^2 + x*y + 3*y^2
            sage: M = matrix(ZZ, [[1, 2], [3, 5]])
            sage: Q.matrix_action_left(M)
            16*x^2 + 83*x*y + 108*y^2
        """
        v, w = M.rows()
        a1 = self(v)
        c1 = self(w)
        b1 = self(v + w) - a1 - c1
        return BinaryQF([a1, b1, c1])

    def matrix_action_right(self, M):
        r"""
        Return the binary quadratic form resulting from the right action
        of the 2-by-2 matrix ``M`` on the quadratic form ``self``.

        Here the action of the matrix `M = \begin{pmatrix} a & b \\ c & d
        \end{pmatrix}` on the form `Q(x, y)` produces the form `Q(ax+by,
        cx+dy)`.

        EXAMPLES::

            sage: Q = BinaryQF([2, 1, 3]); Q
            2*x^2 + x*y + 3*y^2
            sage: M = matrix(ZZ, [[1, 2], [3, 5]])
            sage: Q.matrix_action_right(M)
            32*x^2 + 109*x*y + 93*y^2
        """
        v, w = M.columns()
        a1 = self(v)
        c1 = self(w)
        b1 = self(v + w) - a1 - c1
        return BinaryQF([a1, b1, c1])

    def small_prime_value(self, Bmax=1000):
        r"""
        Returns a prime represented by this (primitive positive definite) binary form.

        INPUT:

        - ``Bmax`` -- a positive bound on the representing integers.

        OUTPUT:

        A prime number represented by the form.

        .. note::

            This is a very elementary implementation which just substitutes
            values until a prime is found.

        EXAMPLES::

            sage: [Q.small_prime_value() for Q in BinaryQF_reduced_representatives(-23, primitive_only=True)]
            [23, 2, 2]
            sage: [Q.small_prime_value() for Q in BinaryQF_reduced_representatives(-47, primitive_only=True)]
            [47, 2, 2, 3, 3]
        """
        from sage.sets.all import Set
        from sage.misc.all import srange
        d = self.discriminant()
        B = 10
        while True:
            llist = list(Set([self(x,y) for x in srange(-B,B) for y in srange(B)]))
            llist = sorted([l for l in llist if l.is_prime()])
            if llist:
                return llist[0]
            if B >= Bmax:
                raise ValueError("Unable to find a prime value of %s" % self)
            B += 10

    def solve_integer(self, n):
        r"""
        Solve `Q(x,y) = n` in integers `x` and `y` where `Q` is this
        quadratic form.

        INPUT:

        - ``Q`` (BinaryQF) -- a positive definite primitive integral
          binary quadratic form

        - ``n`` (int) -- a positive integer

        OUTPUT:

        A tuple (x,y) of integers satisfying `Q(x,y) = n` or ``None``
        if no such `x` and `y` exist.

        EXAMPLES::

            sage: Qs = BinaryQF_reduced_representatives(-23,primitive_only=True)
            sage: Qs
            [x^2 + x*y + 6*y^2, 2*x^2 - x*y + 3*y^2, 2*x^2 + x*y + 3*y^2]
            sage: [Q.solve_integer(3) for Q in Qs]
            [None, (0, 1), (0, 1)]
            sage: [Q.solve_integer(5) for Q in Qs]
            [None, None, None]
            sage: [Q.solve_integer(6) for Q in Qs]
            [(0, 1), (-1, 1), (1, 1)]
        """
        a, b, c = self
        d = self.discriminant()
        if d >= 0 or a <= 0:
            raise ValueError("%s is not positive definite" % self)
        ad = -d
        an4 = 4*a*n
        a2 = 2*a
        from sage.misc.all import srange
        for y in srange(0, 1+an4//ad):
            z2 = an4 + d*y**2
            for z in z2.sqrt(extend=False, all=True):
                if a2.divides(z-b*y):
                    x = (z-b*y)//a2
                    return (x,y)
        return None

    def matrix_Gram(self):
        r"""
        Return the (symmetric) Gram matrix `G` of the quadratic form.

        This expresses `Q(x,y)` as the
        matrix product `Q(x,y) = [x y] * G * [x y]^t`.

        This matrix has base_ring `\QQ`.

        EXAMPLES::
        """
        return Matrix(QQ, 2, 2, [QQ(self._a), QQ(self._b) / 2,
                                 QQ(self._b) / 2, QQ(self._c)])

    def matrix_Hessian(self):
        r"""
        Return the Hessian matrix of the quadratic form.

        This is twice the Gram matrix.

        This matrix has base_ring `\ZZ`.

        EXAMPLES::
        """
        return Matrix(ZZ, 2, 2, [2 * self._a, self._b, self._b, 2 * self._c])

    def eigenspaces(self):
        """
        Return the eigenspaces of the Hessian matrix of this quadratic form.

        .. TODO::

            Cache this, and make sure the BQF is immutable!

        EXAMPLES::
        """
        return self.matrix_Hessian().eigenspaces()

    def spectral_decomposition(self):
        """
        Return the spectral decomposition matrices `(D, X)`.

        The Hessian matrix `H` (which is twice the Gram matrix of `Q`)
        can be written as `H = X^{t} * D * X` where `D` is a diagonal
        matrix and `X` is an orthogonal matrix.

        OUTPUT:

        - `D` -- a diagonal matrix

        - `X` -- an orthogonal matrix

        EXAMPLES::

            sage: B = BinaryQF([1,1,1])
            sage: B.eigenspaces()
            [
            (3, Vector space of degree 2 and dimension 1 over Rational Field
            User basis matrix:
            [1 1]),
            (1, Vector space of degree 2 and dimension 1 over Rational Field
            User basis matrix:
            [ 1 -1])
            ]

            sage: B = BinaryQF([1,0,2])
            sage: B.eigenspaces()
            [
            (4, Vector space of degree 2 and dimension 1 over Rational Field
            User basis matrix:
            [0 1]),
            (2, Vector space of degree 2 and dimension 1 over Rational Field
            User basis matrix:
            [1 0])
            ]
        """
        E = self.eigenspaces()

        # Make the diagonal matrix D
        D_list = []
        for W in E:
            D_list += [W[0][0]] * W[0][1].dimension()
        D = Matrix(QQ, 2, 2, 0)
        for i in range(2):
            D[i, i] = D_list[i]

        # Make the orthogonal matrix X
        M_rows = []
        for W in E:
            for v in W.basis():
                v_len = sqrt(v.dot_product(v))
                M_rows.append(v / v_len)
        X = Matrix(M_rows).transpose()

        # Return the spectral decomposition
        return D, X

    def major_axis(self):
        """
        Return an eigenvector with the smallest eigenvalue.

        This assumes that the quadratic form is definite. If there is
        a single eigenspace, then we return the standard basis vector
        `e_1 = (1,0)`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: B = BinaryQF([1,1,1])
            sage: B.major_axis()
            (1, -1)
            sage: B.minor_axis()
            (1, 1)

            sage: B = BinaryQF([1,0,2])
            sage: B.major_axis()
            (1, 0)
            sage: B.minor_axis()
            (0, 1)

            sage: B = BinaryQF([1,0,1])
            sage: B.major_axis()
            (1, 0)
            sage: B.minor_axis()
            (0, 1)
        """
        E = self.eigenspaces()
        if len(E) == 1:
            return vector([1, 0])
        else:
            # Choose the *smallest* eigenvalue for the major axis
            if abs(E[0][0]) < abs(E[1][0]):
                return E[0][1].basis()[0]
            else:
                return E[1][1].basis()[0]

    def minor_axis(self):
        """
        Return an eigenvector with the largest size eigenvalue.

        This assumes that the quadratic form is definite. If there is
        a single eigenspace, then we return the standard basis vector
        `e_2 = (0,1)`.

        OUTPUT:

        - a vector

        EXAMPLES::

            sage: B = BinaryQF([1,1,1])
            sage: B.major_axis()
            (1, -1)
            sage: B.minor_axis()
            (1, 1)

            sage: B = BinaryQF([1,0,2])
            sage: B.major_axis()
            (1, 0)
            sage: B.minor_axis()
            (0, 1)

            sage: B = BinaryQF([1,0,1])
            sage: B.major_axis()
            (1, 0)
            sage: B.minor_axis()
            (0, 1)
        """
        E = self.eigenspaces()
        if len(E) == 1:
            return vector([0, 1])
        else:
            # Choose the *largest* eigenvalue for the minor axis
            if abs(E[0][0]) > abs(E[1][0]):
                return E[0][1].basis()[0]
            else:
                return E[1][1].basis()[0]

    def plot_level_set(self, m, x_max=None, y_max=None, legend_center=None):
        """
        Return a graphic object showing the level set of the quadratic form.

        This displays the level sets of this binary quadratic form
        with respect to the standard basis, showing both the
        background lattice, the level set ellipse, and the lattice
        points on this level set.

        .. TODO::

            Allow arguments for:
            x_max
            y_max
            legend_center -- ok
            ---------------
            x_range
            y_range
            point_color
            point_size
            ellipse_color
            ellipse_thickness
            intersection_point_size
            intersection_point_color
            return_plot
            aspect_ratio
            legend_center
            hide_legend

        INPUT:

        - m -- an integer

        OUTPUT:

        - a graphics object

        EXAMPLES::

            B = BinaryQF([1,1,1])
            sage: B.plot_level_set(13)
            sage: B.plot_level_set(409)
            sage: B.plot_level_set(11, legend_center=(3,3.5))
        """
        # legend_label='$x^2 + xy + y^2$'
        #  legend_label = str(Q([i,j]))

                    #G += point([(i,j)], size=15, rgbcolor=(1,0,0))
        #G += circle((0,0), sqrt(13), rgbcolor=(1,0,0))

        sqrt_m = sqrt(m)

        # Find the eigenvalues and eigenvectors
        v1 = self.major_axis()
        v2 = self.minor_axis()
        euclid_len1 = sqrt(v1.dot_product(v1))
        euclid_len2 = sqrt(v2.dot_product(v2))
        Q_len1 = sqrt(self(v1))
        Q_len2 = sqrt(self(v2))

        # Find the ellipse information
        if v1[1] < 0:
            upper_v1 = -v1
            # ensure the major axis vector is in upper half-plane, so
            # Arccos gives the correct angle!
        else:
            upper_v1 = v1
        angle = arccos(RR(upper_v1[0]) / euclid_len1)
        major_radius = sqrt_m * euclid_len1 / Q_len1
        minor_radius = sqrt_m * euclid_len2 / Q_len2

        # Determine bounds for our lattice point region
        x_max = ceil(major_radius)
        y_max = x_max

        # Draw the lattice and intersection points
        G = Graphics()
        for i in range(-x_max, x_max + 1):
            for j in range(-y_max, y_max + 1):
                if self([i, j]) != m:
                    G += point([(i, j)], size=4)
                else:
                    G += point([(i, j)], size=15)

        # Draw the shaded ellipse
        G += ellipse((0, 0), major_radius, minor_radius,
                     angle, fill=True, alpha=0.3)

        # Add a label
        if legend_center is None:
            legend_center = (x_max - 1, y_max - 1 + .5)
        G += text('$' + str(self.polynomial()).replace('*', '') + ' = ' + str(m) + '$', legend_center, fontsize=15, rgbcolor=(0, 0, 1))

        # Return the plotted graphic
        G.axes_labels(['$x$', '$y$'])
        G.set_aspect_ratio(1)
        return G.show()


def BinaryQF_reduced_representatives(D, primitive_only=False):
    r"""
    Returns a list of inequivalent reduced representatives for the
    equivalence classes of positive definite binary forms of
    discriminant D.

    INPUT:

    - `D` -- (integer) A negative discriminant.

    - ``primitive_only`` -- (bool, default False) flag controlling whether only
      primitive forms are included.

    OUTPUT:

    (list) A lexicographically-ordered list of inequivalent reduced
    representatives for the equivalence classes of positive definite binary
    forms of discriminant `D`.  If ``primitive_only`` is ``True`` then
    imprimitive forms (which only exist when `D` is not fundamental) are
    omitted; otherwise they are included.

    EXAMPLES::

        sage: BinaryQF_reduced_representatives(-4)
        [x^2 + y^2]

        sage: BinaryQF_reduced_representatives(-163)
        [x^2 + x*y + 41*y^2]

        sage: BinaryQF_reduced_representatives(-12)
        [x^2 + 3*y^2, 2*x^2 + 2*x*y + 2*y^2]

        sage: BinaryQF_reduced_representatives(-16)
        [x^2 + 4*y^2, 2*x^2 + 2*y^2]

        sage: BinaryQF_reduced_representatives(-63)
        [x^2 + x*y + 16*y^2, 2*x^2 - x*y + 8*y^2, 2*x^2 + x*y + 8*y^2, 3*x^2 + 3*x*y + 6*y^2, 4*x^2 + x*y + 4*y^2]

    The number of inequivalent reduced binary forms with a fixed negative
    fundamental discriminant D is the class number of the quadratic field
    `Q(\sqrt{D})`::

        sage: len(BinaryQF_reduced_representatives(-13*4))
        2
        sage: QuadraticField(-13*4, 'a').class_number()
        2
        sage: p=next_prime(2^20); p
        1048583
        sage: len(BinaryQF_reduced_representatives(-p))
        689
        sage: QuadraticField(-p, 'a').class_number()
        689

        sage: BinaryQF_reduced_representatives(-23*9)
        [x^2 + x*y + 52*y^2,
        2*x^2 - x*y + 26*y^2,
        2*x^2 + x*y + 26*y^2,
        3*x^2 + 3*x*y + 18*y^2,
        4*x^2 - x*y + 13*y^2,
        4*x^2 + x*y + 13*y^2,
        6*x^2 - 3*x*y + 9*y^2,
        6*x^2 + 3*x*y + 9*y^2,
        8*x^2 + 7*x*y + 8*y^2]
        sage: BinaryQF_reduced_representatives(-23*9, primitive_only=True)
        [x^2 + x*y + 52*y^2,
        2*x^2 - x*y + 26*y^2,
        2*x^2 + x*y + 26*y^2,
        4*x^2 - x*y + 13*y^2,
        4*x^2 + x*y + 13*y^2,
        8*x^2 + 7*x*y + 8*y^2]

    TESTS::

        sage: BinaryQF_reduced_representatives(5)
        Traceback (most recent call last):
        ...
        ValueError: discriminant must be negative and congruent to 0 or 1 modulo 4
    """
    D = ZZ(D)
    if not ( D < 0 and (D % 4 in [0,1])):
        raise ValueError("discriminant must be negative and congruent to 0 or 1 modulo 4")

    # For a fundamental discriminant all forms are primitive so we need not check:
    if primitive_only:
        primitive_only = not is_fundamental_discriminant(D)

    form_list = []

    from sage.misc.all import xsrange

    # Only iterate over positive a and over b of the same
    # parity as D such that 4a^2 + D <= b^2 <= a^2
    for a in xsrange(1,1+((-D)//3).isqrt()):
        a4 = 4*a
        s = D + a*a4
        w = 1+(s-1).isqrt() if s > 0 else 0
        if w%2 != D%2: w += 1
        for b in xsrange(w,a+1,2):
            t = b*b-D
            if t % a4 == 0:
                c = t // a4
                if (not primitive_only) or gcd([a,b,c])==1:
                    if b>0 and a>b and c>a:
                        form_list.append(BinaryQF([a,-b,c]))
                    form_list.append(BinaryQF([a,b,c]))

    form_list.sort()
    return form_list
