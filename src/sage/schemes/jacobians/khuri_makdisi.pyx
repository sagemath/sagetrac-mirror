r"""
Khuri-Makdisi algorithms for arithmetic in Jacobians

This module implements Khuri-Makdisi's algorithms of [Khu2004]_.

EXAMPLES::

    sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
    sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
    sage: F = C.function_field()
    sage: h = C.function(y/x).divisor_of_poles()
    sage: J = C.jacobian(model='km_large', base_div=h)
    sage: G = J.group()
    sage: b = C([0,1,0]).place()
    sage: pl1 = C([-1,2,1]).place()
    sage: pl2 = C([3,7,1]).place()
    sage: p1 = G.point(pl1 - b)
    sage: p2 = G.point(pl2 - b)

AUTHORS:

- Kwankyu Lee (2022-01): initial version

"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <kwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.matrix.constructor import matrix
from sage.matrix.matrix cimport Matrix

from sage.modules.free_module_element cimport FreeModuleElement

from sage.rings.integer import Integer

cdef inline list listcat(list l):
    flat_list = []
    for sublist in l:
        flat_list.extend(sublist)
    return flat_list


cdef class KhuriMakdisi_base(object):
    cdef Matrix wL
    cdef Matrix w0
    cdef int d0, g

    cdef Matrix mu_image(self, Matrix wd, Matrix we, Matrix mu_mat, int expected_dim=0):
        """
        Lemma 2.2.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='km_large')
            sage: G = J.group()
            sage: b = C([0,0,1]).place()
            sage: pl1 = C([1,2,1]).place()
            sage: pl2 = C([3,1,1]).place()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: w1 = p1._w
            sage: w2 = p2._w
            sage: w1
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 1 0 0 1 0 2]
            [0 0 0 0 1 0 6 0 6]
            [0 0 0 0 0 1 4 0 1]
            [0 0 0 0 0 0 0 1 0]
            sage: w2
            [0 1 0 0 0 0 0 2 2]
            [0 0 1 0 0 0 0 4 6]
            [0 0 0 1 0 0 0 5 6]
            [0 0 0 0 1 0 0 3 1]
            [0 0 0 0 0 1 0 3 5]
            [0 0 0 0 0 0 1 6 3]
            sage: (p1 + p2)._w  # indirect doctest
            [1 0 0 0 0 0 2 2 4]
            [0 1 0 0 0 0 3 5 0]
            [0 0 1 0 0 0 0 3 6]
            [0 0 0 1 0 0 1 1 3]
            [0 0 0 0 1 0 5 6 3]
            [0 0 0 0 0 1 5 5 3]
        """
        cdef Matrix mat
        cdef FreeModuleElement v
        cdef Py_ssize_t n, c, r

        n = we.ncols()
        c = mu_mat.ncols() // n
        mat = matrix(0, c)
        for v in wd:
            mat = mat.stack(we * matrix(n, v * mu_mat))
            mat.echelonize()
            r = mat.rank()
            mat = mat.matrix_from_rows(range(r))
            if expected_dim and r == expected_dim:
                break

        #if expected_dim:
        #    assert r == expected_dim

        return mat

    cdef Matrix mu_preimage(self, Matrix we, Matrix wde, Matrix mu_mat_reversed, int expected_codim=0):
        """
        Lemma 2.3 (division).

        This computes

          {s: mu(s*E) subset F} = {s: s*M*E*F_perp^ == 0}
                                = {s: s*M*v*F_perp^ == 0 for v in E}
                                = {s: F_perp*(v*M^)*s == 0 for v in E}

        for E = we, F = wde, M^ = mu_mat_reversed
        """
        cdef Matrix mat, perp, vmu
        cdef FreeModuleElement v
        cdef Py_ssize_t nd, ne, nde, r

        ne = we.ncols()
        nde = wde.ncols()
        nd = nde - ne

        perp = wde.right_kernel_matrix()
        mat = matrix(0, nd)
        for v in we:
            # efficient than
            #   vmu = matrix([v * matrix(ne, row) for row in mu_mat])
            vmu = matrix(nd, v * mu_mat_reversed)
            mat = mat.stack(perp * vmu.transpose())
            mat.echelonize()
            r = mat.rank()
            mat = mat.matrix_from_rows(range(r))
            if expected_codim and r == expected_codim:
                break

        #if expected_codim:
        #    assert r == expected_codim

        return mat.right_kernel_matrix()

    cpdef Matrix negate(self, Matrix wd):
        """
        Theorem 4.4 (negation), first method.
        """
        return self.addflip(wd, self.w0)

    cpdef Matrix add(self, Matrix wd1, Matrix wd2):
        """
        Theorem 4.5 (addition).
        """
        return self.negate(self.addflip(wd1, wd2))

    cpdef Matrix subtract(self, Matrix wd1, Matrix wd2):
        """
        Theorem 4.6 (subtraction), first method.
        """
        return self.addflip(self.negate(wd1), wd2)

    cpdef Matrix multiple(self, Matrix wd, n):
        """
        Compute multiple by additive square-and-multiply method.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl = C([-1,2,1]).place()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p = G.point(pl - b)
            sage: 10*p
            Point of Jacobian determined by
            [ 1  0  0  0  0  0  5  2  2]
            [ 0  1  0  0  0  0 13  6 11]
            [ 0  0  1  0  0  0  1 11  4]
            [ 0  0  0  1  0  0  1 13  7]
            [ 0  0  0  0  1  0 12 16  2]
            [ 0  0  0  0  0  1  6  9 10]
            sage: (-10)*p
            Point of Jacobian determined by
            [ 1  0  0  0  0 13  0 10  6]
            [ 0  1  0  0  0  5  0  4 16]
            [ 0  0  1  0  0  2  0  0  4]
            [ 0  0  0  1  0  9  0  6  9]
            [ 0  0  0  0  1  6  0  0  9]
            [ 0  0  0  0  0  0  1  9  5]
            sage: 0*p
            Point of Jacobian determined by
            [1 0 0 0 0 0 0 0 0]
            [0 1 0 0 0 0 0 0 0]
            [0 0 1 0 0 0 0 0 0]
            [0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 1 0 0 0]
            [0 0 0 0 0 0 0 1 0]
        """
        cdef Matrix w
        cdef int sign, b
        cdef list bits

        if n == 0:
            return self.w0
        if n < 0:
            bits = Integer(-n).digits(2)
        else:
            bits = Integer(n).digits(2)
        bits.pop()
        mwd = None
        w = wd
        sign = 1
        for i in range(len(bits)):
            w = self.addflip(w, w)
            sign = -sign
            b = bits.pop()
            if b > 0:
                if sign < 0:
                    if mwd is None:
                        mwd = self.addflip(wd, self.w0)
                    w = self.addflip(w, mwd)
                else:
                    w = self.addflip(w, wd)
                sign = -sign
        if sign < 0 and n > 0 or sign > 0 and n < 0:
            w = self.addflip(w, self.w0)  # negate
        return w

    cpdef Matrix zero_divisor(self):
        """
        Return the matrix representing the zero divisor.
        """
        return self.wL

cdef class KhuriMakdisi_large(KhuriMakdisi_base):
    """
    Khuri-Makdisi's *large* model

    Let `C` be a geometrically smooth projective curve of genus `g`.
    Let `D_0` be a divisor of `C` of degree at least `2g + 1`.
    Let `L_0` denote the line bundle `O_C(D_0)`.

    Then `V_n` denote the linear space of global sections of the line bundle
    `L_0^\oplus n` the `n` times tensor product of `L_0`

    INPUT:

    - ``V`` -- a function that returns the space of global section

    TESTS::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
        sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
        sage: h = C.function(y/x).divisor_of_poles()
        sage: b = C([0,1,0]).place()
        sage: pl1 = C([-1,2,1]).place()
        sage: pl2 = C([3,7,1]).place()
        sage: J = C.jacobian(model='km_large', base_div=h)
        sage: G = J.group()
        sage: p1 = G.point(pl1 - b)
        sage: p2 = G.point(pl2 - b)
        sage: p1
        Point of Jacobian determined by
        [ 1  0  0  0  0  0  0 12 15]
        [ 0  1  0  0  0  0  0  0 13]
        [ 0  0  1  0  0  0  0  0  2]
        [ 0  0  0  1  0  0  0  0 16]
        [ 0  0  0  0  0  1  0  0 15]
        [ 0  0  0  0  0  0  1  0  1]
        sage: p2
        Point of Jacobian determined by
        [ 1  0  0  0  0  0  0 12  5]
        [ 0  1  0  0  0  0  0  0  2]
        [ 0  0  1  0  0  0  0  0 13]
        [ 0  0  0  1  0  0  0  0  8]
        [ 0  0  0  0  0  1  0  0 10]
        [ 0  0  0  0  0  0  1  0 14]
        sage: p1 + p2
        Point of Jacobian determined by
        [ 1  0  0  0  0 16  0  5  3]
        [ 0  1  0  0  0  6  0  8 16]
        [ 0  0  1  0  0 15  0  3 10]
        [ 0  0  0  1  0  3  0  0  0]
        [ 0  0  0  0  1 12  0 16  8]
        [ 0  0  0  0  0  0  1  3  0]
        sage: p1 - p2
        Point of Jacobian determined by
        [ 1  0  0  0  0  0 13  9  5]
        [ 0  1  0  0  0  0  2  5  8]
        [ 0  0  1  0  0  0  6  7  5]
        [ 0  0  0  1  0  0 11  3 16]
        [ 0  0  0  0  1  0  9  7 10]
        [ 0  0  0  0  0  1  4 10  5]
        sage: p1.addflip(p2) == -(p1 + p2)
        True
    """
    cdef Matrix mu_mat33

    def __init__(self, V, mu, w0, d0, g):
        """
        Initialize.
        """
        self.wL = V(3).basis_matrix()
        self.w0 = w0
        self.d0 = d0
        self.g = g
        self.mu_mat33 = matrix(listcat([list(mu(3, 3, i, j)) for j in range(3*d0-g+1)]) for i in range(3*d0-g+1))

    def equal(self, wd, we):
        """
        Theorem 4.1, second method.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(7), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: J = C.jacobian(model='km_large')
            sage: b = C([0,1,0]).place()
            sage: J.set_base_place(b)
            sage: G = J.group()
            sage: pl1 = C([3,2,1]).place()
            sage: pl2 = C([5,5,1]).place()
            sage: p1 = G(pl1)
            sage: p2 = G(pl2)
            sage: p1 + p2 == p2 + p1
            True
            sage: p1 - p2 == -(p2 - p1)
            True
            sage: zero = G.zero()
            sage: p1 + zero == p1
            True
            sage: p1 - p1 == zero
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix wf, w1, w2

        wf = matrix(wd[0])
        w1 = self.mu_image(wf, we, self.mu_mat33, 2*d0 - g + 1)
        w2 = self.mu_preimage(wd, w1, self.mu_mat33)
        return w2.rank() > 0

    cdef Matrix _add(self, Matrix wd, Matrix we):
        """
        Theorem 3.6 (addition of divisors, first method).

        This is used only for addition in Jacobian. So we assume that both wd,
        we represents divisors of degree d0.
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2

        w1 = self.mu_image(wd, we, self.mu_mat33, 4*d0 - g + 1)
        w2 = self.mu_preimage(self.wL, w1, self.mu_mat33, 2*d0)
        return w2

    cdef Matrix _flip(self, Matrix wd):
        """
        Theorem 3.10 (flipping)
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2

        # efficient than
        #   wf = matrix(wd[0])
        #   w1 = self.mu_image(wf, self.wL, mu_mat, 3*d0 - g + 1)
        w1 = matrix(3*d0 - g + 1, wd[0] * self.mu_mat33)
        w2 = self.mu_preimage(wd, w1, self.mu_mat33, d0)
        return w2

    cpdef Matrix addflip(self, Matrix wd1, Matrix wd2):
        """
        Theorem 4.3 (addflip)

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: F = C.function_field()
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1.addflip(p2)
            Point of Jacobian determined by
            [ 1  0  0  0  0  0  7 10  9]
            [ 0  1  0  0  0  0  4 14 10]
            [ 0  0  1  0  0  0  7  0  9]
            [ 0  0  0  1  0  0 10 10  6]
            [ 0  0  0  0  1  0  6  5 15]
            [ 0  0  0  0  0  1 14  9  1]
        """
        return self._flip(self._add(wd1, wd2))

    cpdef Matrix add_divisor(self, Matrix wd1, Matrix wd2, int d1, int d2):
        """
        Theorem 3.6 (addition of divisors, first method).

        We assume that wd1, wd2 represent divisors of degree at most
        3d0 - 2g - 1.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='km_large')
            sage: G = J.group()
            sage: pts = G.get_points(G.order())  # indirect doctest
            sage: len(pts)
            11
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2

        w1 = self.mu_image(wd1, wd2, self.mu_mat33, 6*d0 - d1 - d2 - g + 1)
        w2 = self.mu_preimage(self.wL, w1, self.mu_mat33, d1 + d2)
        return w2


cdef class KhuriMakdisi_medium(KhuriMakdisi_base):
    """
    Khuri-Makdisi's *medium* model

    Let `L` be a line bundle of a curve.

    Let `V_n` denote the vector space of global vectors of a line bundle...

    INPUT:

    - ``V`` -- function that returns

    TESTS::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
        sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
        sage: h = C.function(y/x).divisor_of_poles()
        sage: b = C([0,1,0]).place()
        sage: pl1 = C([-1,2,1]).place()
        sage: pl2 = C([3,7,1]).place()
        sage: J = C.jacobian(model='km_medium', base_div=h)
        sage: G = J.group()
        sage: p1 = G.point(pl1 - b)
        sage: p2 = G.point(pl2 - b)
        sage: p1
        Point of Jacobian determined by
        [ 1  0  0  0 16 12]
        [ 0  1  0  0 15  0]
        [ 0  0  1  0  1  0]
        sage: p2
        Point of Jacobian determined by
        [ 1  0  0  0  8 12]
        [ 0  1  0  0 10  0]
        [ 0  0  1  0 14  0]
        sage: p1 + p2
        Point of Jacobian determined by
        [ 1  0  0  6  3 16]
        [ 0  1  0 15 16 10]
        [ 0  0  1  3  0  0]
        sage: p1 - p2
        Point of Jacobian determined by
        [ 1  0  0  8  0 14]
        [ 0  1  0  1 10 10]
        [ 0  0  1 15  3  6]
        sage: p1.addflip(p2) == -(p1 + p2)
        True
    """
    cdef Matrix wV1, wV2, wV3, mu_mat22, mu_mat13, mu_mat32

    def __init__(self, V, mu, w0, d0, g):
        """
        Initialize.
        """
        self.wL = V(2).basis_matrix()
        self.w0 = w0
        self.d0 = d0
        self.g = g
        self.wV1 = V(1).basis_matrix()
        self.wV2 = V(2).basis_matrix()
        self.wV3 = V(3).basis_matrix()
        self.mu_mat22 = matrix(listcat([list(mu(2, 2, i, j)) for j in range(2*d0-g+1)]) for i in range(2*d0-g+1))
        self.mu_mat13 = matrix(listcat([list(mu(1, 3, i, j)) for j in range(3*d0-g+1)]) for i in range(1*d0-g+1))
        self.mu_mat32 = matrix(listcat([list(mu(3, 2, i, j)) for j in range(2*d0-g+1)]) for i in range(3*d0-g+1))

    def equal(self, wd, we):
        """
        Theorem 4.1, second method.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: J = C.jacobian(model='km_medium', base_div=h)
            sage: G = J.group()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 + p2 == p2 + p1
            True
            sage: p1 - p2 == -(p2 - p1)
            True
            sage: zero = G.zero()
            sage: p1 + zero == p1
            True
            sage: p1 - p1 == zero
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix wf, w1, w2

        wf = matrix(wd[0])
        w1 = self.mu_image(wf, we, self.mu_mat22, d0 - g + 1)
        w2 = self.mu_preimage(wd, w1, self.mu_mat22)
        return w2.rank() > 0

    cpdef Matrix addflip(self, Matrix wd1, Matrix wd2):
        """
        Theorem 5.1 (addflip in medium model).

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: h = C.function(y/x).divisor_of_poles()
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: J = C.jacobian(model='km_medium', base_div=h)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: af = p1.addflip(p2)
            sage: af
            Point of Jacobian determined by
            [ 1  0  0  6  3 16]
            [ 0  1  0  0  7  9]
            [ 0  0  1 14  2  3]

        We check the computation in other model::

            sage: J = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: G.point(af.divisor()) == p1.addflip(p2)
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2, w3, w4

        w1 = self.mu_image(wd1, wd2, self.mu_mat22, 2*d0 - g + 1)
        w2 = self.mu_preimage(self.wV1, w1, self.mu_mat13, 2*d0)
        # efficient than
        #   wf = matrix(w2[0])
        #   w3 = self.mu_image(wf, self.wV2, self.mu_mat32, 2*d0 - g + 1)
        w3 = matrix(2*d0 - g + 1, w2[0] * self.mu_mat32)
        w4 = self.mu_preimage(w2, w3, self.mu_mat32, d0)
        return w4

    cpdef Matrix add_divisor(self, Matrix wd1, Matrix wd2, int d1, int d2):
        """
        Theorem 3.6 (addition of divisors, first method).

        We assume that wd1, wd2 represent divisors of degree at most
        4d0 - 2g - 1.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='km_medium')
            sage: G = J.group()
            sage: pts = G.get_points(G.order())  # indirect doctest
            sage: len(pts)
            11
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2

        w1 = self.mu_image(wd1, wd2, self.mu_mat22, 4*d0 - d1 - d2 - g + 1)
        w2 = self.mu_preimage(self.wL, w1, self.mu_mat22, d1 + d2)
        return w2

cdef class KhuriMakdisi_small(KhuriMakdisi_base):
    """
    Khuri-Makdisi's *small* model

    Let `L` be a line bundle of a curve.

    Let `V_n` denote the vector space of global vectors of a line bundle...

    INPUT:

    - ``V`` -- function that returns

    TESTS::

        sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
        sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
        sage: b = C([0,1,0]).place()
        sage: pl1 = C([-1,2,1]).place()
        sage: pl2 = C([3,7,1]).place()
        sage: J = C.jacobian(model='km_small', base_div=2*b)
        sage: G = J.group()
        sage: p1 = G.point(pl1 - b)
        sage: p2 = G.point(pl2 - b)
        sage: p1
        Point of Jacobian determined by
        [ 1  0  0  0  0 11]
        [ 0  1  0  0  2  0]
        [ 0  0  1  0 16 10]
        [ 0  0  0  1  0  3]
        sage: p2
        Point of Jacobian determined by
        [1 0 0 0 0 3]
        [0 1 0 0 7 0]
        [0 0 1 0 3 5]
        [0 0 0 1 0 2]
        sage: p1 + p2
        Point of Jacobian determined by
        [ 1  0  0  0 10  9]
        [ 0  1  0  0  7  5]
        [ 0  0  1  0 15  4]
        [ 0  0  0  1  3 10]
        sage: p1 - p2
        Point of Jacobian determined by
        [ 1  0  0  0 10  9]
        [ 0  1  0  0  9  8]
        [ 0  0  1  0 15  4]
        [ 0  0  0  1 15 11]
        sage: p1.addflip(p2) == -(p1 + p2)
        True
    """
    cdef Matrix wV2, wV3, wV4, mu_mat22
    cdef Matrix mu_mat23, mu_mat24, mu_mat32, mu_mat33, mu_mat34, mu_mat42, mu_mat43

    def __init__(self, V, mu, w0, d0, g):
        """
        Initialize.
        """
        self.wL = V(3).basis_matrix()
        self.w0 = w0
        self.d0 = d0
        self.g = g
        self.wV2 = V(2).basis_matrix()
        self.wV3 = V(3).basis_matrix()
        self.wV4 = V(4).basis_matrix()
        self.mu_mat23 = matrix(listcat([list(mu(2, 3, i, j)) for j in range(3*d0-g+1)]) for i in range(2*d0-g+1))
        self.mu_mat24 = matrix(listcat([list(mu(2, 4, i, j)) for j in range(4*d0-g+1)]) for i in range(2*d0-g+1))
        self.mu_mat32 = matrix(listcat([list(mu(3, 2, i, j)) for j in range(2*d0-g+1)]) for i in range(3*d0-g+1))
        self.mu_mat33 = matrix(listcat([list(mu(3, 3, i, j)) for j in range(3*d0-g+1)]) for i in range(3*d0-g+1))
        self.mu_mat34 = matrix(listcat([list(mu(3, 4, i, j)) for j in range(4*d0-g+1)]) for i in range(3*d0-g+1))
        self.mu_mat42 = matrix(listcat([list(mu(4, 2, i, j)) for j in range(2*d0-g+1)]) for i in range(4*d0-g+1))
        self.mu_mat43 = matrix(listcat([list(mu(4, 3, i, j)) for j in range(3*d0-g+1)]) for i in range(4*d0-g+1))

    def equal(self, wd, we):
        """
        Theorem 4.1, second method.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: J = C.jacobian(model='km_small', base_div=2*b)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: p1 + p2 == p2 + p1
            True
            sage: p1 - p2 == -(p2 - p1)
            True
            sage: zero = G.zero()
            sage: p1 + zero == p1
            True
            sage: p1 - p1 == zero
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix wf, w1, w2

        wf = matrix(wd[0])
        w1 = self.mu_image(wf, we, self.mu_mat33, 2*d0 - g + 1)
        w2 = self.mu_preimage(wd, w1, self.mu_mat33)
        return w2.rank() > 0

    cpdef Matrix addflip(self, Matrix wd1, Matrix wd2):
        """
        Theorem 5.3 (addflip in small model), second method.

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: J = C.jacobian(model='km_small', base_div=2*b)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: af = p1.addflip(p2)
            sage: af
            Point of Jacobian determined by
            [ 1  0  0  0 10  9]
            [ 0  1  0  0 10 12]
            [ 0  0  1  0 15  4]
            [ 0  0  0  1 14  7]

        We check the computation in other model::

            sage: h = C.function(y/x).divisor_of_poles()
            sage: Jl = C.jacobian(model='km_large', base_div=h)
            sage: G = J.group()
            sage: q1 = G.point(pl1 - b)
            sage: q2 = G.point(pl2 - b)
            sage: G.point(af.divisor()) == q1.addflip(p2)
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2, w3, w4, w5

        w1 = self.mu_image(wd1, wd2, self.mu_mat33, 4*d0 - g + 1)
        w2 = self.mu_preimage(self.wV3, w1, self.mu_mat33, 2*d0)
        w3 = self.mu_preimage(self.wV2, w1, self.mu_mat24, 2*d0)
        # efficient than
        #   wf = matrix(w2[0])
        #   w4 = self.mu_image(wf, self.wV4, self.mu_mat34, 4*d0 - g + 1)
        w4 = matrix(4*d0 - g + 1, w2[0] * self.mu_mat34)
        w5 = self.mu_preimage(w3, w4, self.mu_mat43, d0)
        return w5

    cpdef Matrix negate(self, Matrix wd):
        """
        Theorem 5.4 (negation in small model).

        TESTS::

            sage: P2.<x,y,z> = ProjectiveSpace(GF(17), 2)
            sage: C = Curve(x^3 + 5*z^3- y^2*z, P2)
            sage: b = C([0,1,0]).place()
            sage: pl1 = C([-1,2,1]).place()
            sage: pl2 = C([3,7,1]).place()
            sage: J = C.jacobian(model='km_small', base_div=2*b)
            sage: G = J.group()
            sage: p1 = G.point(pl1 - b)
            sage: p2 = G.point(pl2 - b)
            sage: -(-p1) == p1
            True
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2, w3, w4

        w1 = self.mu_image(self.wV2, wd, self.mu_mat23, 4*d0 - g + 1)
        w2 = self.mu_preimage(self.wV3, w1, self.mu_mat32, d0)
        # efficient than
        #   wf = matrix(w2[0])
        #   w3 = self.mu_image(wf, self.wV4, self.mu_mat24, 4*d0 - g + 1)
        w3 = matrix(4*d0 - g + 1, w2[0] * self.mu_mat24)
        w4 = self.mu_preimage(wd, w3, self.mu_mat33, d0)
        return w4

    cpdef Matrix add_divisor(self, Matrix wd1, Matrix wd2, int d1, int d2):
        """
        Theorem 3.6 (addition of divisors, first method).

        We assume that wd1, wd2 represent divisors of degree at most
        6d0 - 2g - 1.

        TESTS::

            sage: k = GF(7)
            sage: A.<x,y> = AffineSpace(k,2)
            sage: C = Curve(y^2 + x^3 + 2*x + 1).projective_closure()
            sage: J = C.jacobian(model='km_small')
            sage: G = J.group()
            sage: pts = G.get_points(G.order())  # indirect doctest
            sage: len(pts)
            11
        """
        cdef int d0 = self.d0
        cdef int g = self.g
        cdef Matrix w1, w2

        w1 = self.mu_image(wd1, wd2, self.mu_mat33, 6*d0 - d1 - d2 - g + 1)
        w2 = self.mu_preimage(self.wL, w1, self.mu_mat33, d1 + d2)
        return w2
