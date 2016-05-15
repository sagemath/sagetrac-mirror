r"""
A dimension formula for vector-valued modular forms, and functions
that apply it to the case of Jacobi forms.

.. TODO::

    For Jacobi forms there is a much better way to implement this
    formula. This is ongoing work by Ehlen et al.  Replace this code by
    his work as soon as possible. However, for general vector valued
    modular forms with diagonalizable representation matrix `\rho([1, 1; 0,
    1])`, it is necessary to keep this method.

AUTHOR:

- Martin Raum
"""

#===============================================================================
#
# Copyright (C) 2012-2014 Martin Raum
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from sage.functions.all import exp, sqrt
from sage.matrix.all import diagonal_matrix, identity_matrix, matrix
from sage.misc.all import sum, mrange, prod, cython_lambda
from sage.modules.all import vector
from sage.rings.all import ComplexIntervalField, ZZ, QQ, lcm
from sage.quadratic_forms.all import QuadraticForm
from sage.symbolic.all import I, pi
from copy import copy
import operator


def jacobi_dimension(k, m):
    r"""
    INPUT:

    - `k` -- An integer.

    - `m` -- A quadratic form or an even symmetric matrix (over `\ZZ`).

    TESTS::

        sage: from sage.modular.jacobi.classical import _classical_jacobi_forms_as_weak_jacobi_forms
        sage: from sage.modular.jacobi.higherrank_dimension import jacobi_dimension
        sage: all(len(_classical_jacobi_forms_as_weak_jacobi_forms(k, m)) == jacobi_dimension(k, matrix([[2 * m]])) for k in range(8, 16) for m in range(1, 10) ) # long time
        True
    """
    from sage.matrix.matrix import is_Matrix

    if is_Matrix(m):
        return vector_valued_dimension(k - ZZ(m.ncols()) / 2, QuadraticForm(-m))
    else:
        return vector_valued_dimension(k - ZZ(m.dim()) / 2, m.scale_by_factor(-1))


def nmb_isotropic_vectors(k, L):
    r"""
    Compute the number of isotropic vectors, which generically give
    rise to Eisenstein series of weight `k`.

    INPUT:

    - `k` -- A half-integer.

    - `L` -- A quadratic form over `\ZZ`.

    OUTPUT:

    - An integer.

    TESTS::

        sage: from sage.modular.jacobi.higherrank_dimension import nmb_isotropic_vectors
        sage: nmb_isotropic_vectors(1/2, QuadraticForm(matrix(1, [[2]])))
        0
        sage: nmb_isotropic_vectors(-1/2, QuadraticForm(matrix(1, [[2]])))
        1
    """
    (discriminant_form_exponents, disc_quadratic, disc_bilinear) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, discriminant_form_exponents)

    plus_basis = ZZ(L.dim() + 2*k) % 4 == 0

    return len([a for a in (singls + pairs if plus_basis else pairs)
                if disc_quadratic(*a) in ZZ])


def vector_valued_dimension(k, L):
    r"""
    Compute the dimension of the space of weight `k` vector valued
    modular forms for the Weil representation attached to the lattice
    `L`.

    See [Borcherds, Borcherds - Reflection groups of Lorentzian
    lattices] for a proof of the formula that we use here.

    INPUT:

    - `k` -- A half-integer.

    - ``L`` -- An quadratic form.

    OUTPUT:

    An integer.

    TESTS::

        sage: from sage.modular.jacobi.higherrank_dimension import vector_valued_dimension
        sage: vector_valued_dimension(3, QuadraticForm(-matrix(2, [2, 1, 1, 2])))
        1
        sage: vector_valued_dimension(3, QuadraticForm(-matrix(2, [2, 0, 0, 2])))
        1
        sage: vector_valued_dimension(3, QuadraticForm(-matrix(2, [2, 0, 0, 4])))
        1
    """
    if 2 * k not in ZZ:
        raise ValueError( "Weight must be half-integral" )
    if k <= 0:
        return 0
    if k < 2:
        raise NotImplementedError( "Weight <2 is not implemented." )

    if L.matrix().rank() != L.matrix().nrows():
        raise ValueError( "The lattice (={0}) must be non-degenerate.".format(L) )

    if L.dim() % 2 != ZZ(2 * k) % 2:
        return 0

    (discriminant_form_exponents, disc_quadratic, disc_bilinear) = _discriminant_form(L)
    (singls, pairs) = _discriminant_form_pmone(L, discriminant_form_exponents)
    plus_basis = ZZ(L.dim() + 2*k) % 4 == 0

    if plus_basis:
        subspace_dimension = len(singls + pairs)
    else:
        subspace_dimension = len(pairs)

    CC_prec = 50 + subspace_dimension * 2
    while True:
        CC = ComplexIntervalField(CC_prec)

        (Smat, Tmat) = _weil_representation(CC, L, singls, pairs, plus_basis,
                                            discriminant_form_exponents, disc_quadratic, disc_bilinear)
        STmat = Smat * Tmat


        ## This function overestimates the number of eigenvalues, if it is not correct
        def eigenvalue_multiplicity(mat, ev):
            mat = matrix(CC, mat - ev * identity_matrix(subspace_dimension))
            return len(filter( lambda row: all( e.contains_zero() for e in row), _qr(mat).rows() ))

        rti = CC(exp(2 * pi * I / 8))
        S_ev_multiplicity = [eigenvalue_multiplicity(Smat, rti**n) for n in range(8)]
        ## Together with the fact that eigenvalue_multiplicity overestimates the multiplicities
        ## this asserts that the computed multiplicities are correct
        assert sum(S_ev_multiplicity) == subspace_dimension

        rho = CC(exp(2 * pi * I / 12))
        ST_ev_multiplicity = [eigenvalue_multiplicity(STmat, rho**n) for n in range(12)]
        ## Together with the fact that eigenvalue_multiplicity overestimates the multiplicities

        if sum(ST_ev_multiplicity) == subspace_dimension:
            break
        else:
            CC_prec += 100


    normalize_QQ = lambda a: (a.numerator() % a.denominator()) / a.denominator()
    T_evs = [ normalize_QQ(disc_quadratic(*a))
              for a in (singls + pairs if plus_basis else pairs) ]


    return (subspace_dimension * (1 + QQ(k) / 12)
           - ZZ(sum( (ST_ev_multiplicity[n] * ((-2 * k - n) % 12)) for n in range(12) )) / 12
           - ZZ(sum( (S_ev_multiplicity[n] * ((2 * k + n) % 8)) for n in range(8) )) / 8
           - sum(T_evs))

def _discriminant_form(L):
    r"""
    Compute the discriminant form attached to a quadratic form `L`.

    INPUT:

    - `L` -- A quadratic form.

    OUTPUT:

    - A triple ``(elementary_divisors_inv, disc_quadratic,
      disc_bilinear)``.  The first component is a list of positive
      integers, say, length `l` that correspond to elementary divisors
      of the discriminant form.  The second and third component are
      functions on `l` and `2l` arguments, respectively.

    EXAMPLES::

        sage: L = QuadraticForm(ZZ, 2, [1,1,1])
        sage: from sage.modular.jacobi.higherrank_dimension import _discriminant_form
        sage: _discriminant_form(L)
        ([3], <function ..., <function ...)

    TESTS:

    See test_higherrank_dimension:test__discrimant_form.
    """
    if L.matrix().rank() != L.matrix().nrows():
        raise ValueError("The lattice (={0}) must be non-degenerate.".format(L))

    ## The bilinear and the quadratic form attached to L
    # quadratic = lambda x: L(x) // 2
    # bilinear = lambda x,y: L(x + y) - L(x) - L(y)


    ## A dual basis for L
    (elementary_divisors, dual_basis_pre, _) = L.matrix().smith_form()
    elementary_divisors = elementary_divisors.diagonal()
    dual_basis = map(operator.div, list(dual_basis_pre), elementary_divisors)

    L_level = ZZ(lcm([ b.denominator() for b in dual_basis ]))

    (elementary_divisors, _, discriminant_basis_pre) = (L_level * matrix(dual_basis)).change_ring(ZZ).smith_form()
    elementary_divisors = filter( lambda d: d not in ZZ, (elementary_divisors / L_level).diagonal() )
    elementary_divisors_inv = map(ZZ, [ed ** -1 for ed in elementary_divisors])
    discriminant_basis = matrix(map( operator.mul,
                                     discriminant_basis_pre.inverse().rows()[:len(elementary_divisors)],
                                     elementary_divisors )).transpose()

    ## This is a form over QQ, so that we cannot use an instance of QuadraticForm
    discriminant_form = discriminant_basis.transpose() * L.matrix() * discriminant_basis
    if prod(elementary_divisors_inv) > 100:
        disc_den = discriminant_form.denominator()
        disc_bilinear_pre = \
            cython_lambda( ', '.join(   ['int a{0}'.format(i) for i in range(discriminant_form.nrows())]
                                        + ['int b{0}'.format(i) for i in range(discriminant_form.nrows())] ),
                           ' + '.join('{0} * a{1} * b{2}'.format(disc_den * discriminant_form[i,j], i, j)
                                      for i in range(discriminant_form.nrows())
                                      for j in range(discriminant_form.nrows())) )
        disc_bilinear = lambda *a: disc_bilinear_pre(*a) / disc_den
    else:
        disc_bilinear = lambda *xy: vector(ZZ, xy[:discriminant_form.nrows()]) * discriminant_form * vector(ZZ, xy[discriminant_form.nrows():])

    disc_quadratic = lambda *a: disc_bilinear(*(2 * a)) / 2

    return (elementary_divisors_inv, disc_quadratic, disc_bilinear)

def _discriminant_form_pmone(L, discriminant_form_exponents):
    r"""
    Partition the discriminant form into elements that are fixed by
    multiplication by `-1` and those which are not.

    INPUT:

    - `L` -- A quadratic form over `\ZZ`.

    - ``discriminant_form_exponents`` -- A list of positive integers
                                         that correspond to the
                                         exponents of the Jordan
                                         components.

    OUTPUT:

    - A pair ``(singls, pairs)``, the first component of which
      contains a list of elements that are fixed by `-1` and the
      second component of which consists of repepresentatives of pairs
      of elements which are not.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank_dimension import _discriminant_form_pmone
        sage: L = QuadraticForm(ZZ, 2, [1,1,1])
        sage: _discriminant_form_pmone(L, [3])
        ([[0]], [[2]])

    TESTS:

    See test_higherrank_dimension:test__discrimant_form_pmone.
    """
    ## red gives a normal form for elements in the discriminant group
    red = lambda x: map(operator.mod, x, discriminant_form_exponents)

    ## singls and pairs are elements of the discriminant group that are, respectively,
    ## fixed and not fixed by negation.
    singls = list()
    pairs = list()
    for x in mrange(discriminant_form_exponents):
        y = red(map(operator.neg, x))
        for (e, f) in zip(x, y):
            if e < f:
                si = -1
                break
            elif e > f:
                si = 1
                break
        else:
            singls.append(x)
            continue

        if si == 1:
            pairs.append(x)

    return (singls, pairs)


def _weil_representation(CC, L, singls, pairs, plus_basis,
                         discriminant_form_exponents, disc_quadratic,
                         disc_bilinear):
    r"""
    Construct the Weil representation with values in a complex field
    (or interval field).

    INPUT:

    - ``CC`` -- A complex field, a complex interval field, or any
                field that allows for conversion of symbolic
                expressions.

    - `L` -- A quadratic form over `\ZZ`.

    - ``singls`` -- A list of tuples of integers.  See meth:`_discriminant_form_pmone`.

    - ``pairs`` -- A list of tuples of integers.  See meth:`_discriminant_form_pmone`.

    - ``plus_basis`` -- Boolean.  Compute the Weil representation on
                        the span of `e_v + e_{-v}` or on the span of
                        `e_v - e_{-v}`.

    - ``discriminant_form_exponents`` -- A list of positive integers.
                                         See `meth:_discriminant_form`.

    - ``disc_quadratic`` -- A function.  See `meth:_discriminant_form`.

    - ``disc_bilinear`` -- A function.  See `meth:_discriminant_form`.

    OUTPUT:

    - A pair of matrices `S` and `T` that corrspond to the image under
      the Weil representation of `((0,-1; 1,0), \sqrt{\tau})` and
      `((1,1; 0,1), 1)`.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank_dimension import _weil_representation, _discriminant_form, _discriminant_form_pmone
        sage: CI = ComplexIntervalField(10)
        sage: L = QuadraticForm(ZZ, 2, [1,1,1])
        sage: (eds, quad, bil) = _discriminant_form(L)
        sage: (singls, pairs) = _discriminant_form_pmone(L, eds)
        sage: _weil_representation(CI, L, singls, pairs, True, eds, quad, bil)
        (
        [0.0? + 0.6?*I 0.0? + 0.8?*I]  [           1            0]
        [0.0? + 0.8?*I   0.? - 1.?*I], [           0 -1.? + 1.?*I]
        )
        sage: _weil_representation(CI, L, singls, pairs, False, eds, quad, bil)
        ([-1.? + 0.?*I], [-1.? + 1.?*I])

    TESTS:

    See test_higherrank_dimension:test__weil_representation.
    """
    zeta_order = ZZ(lcm([8, 12] + map(lambda ex: 2 * ex,
                                      discriminant_form_exponents)))

    zeta = CC(exp(2 * pi * I / zeta_order))
    sqrt2 = CC(sqrt(2))
    drt = CC(sqrt(abs(L.det())))

    Tmat = diagonal_matrix(CC, [zeta ** (zeta_order * disc_quadratic(*a))
                                for a in (singls + pairs if plus_basis
                                          else pairs)])

    neg = lambda v: map(operator.neg, v)

    if plus_basis:
        Smat = zeta ** (zeta_order / 8 * L.dim()) / drt * matrix(CC, [[zeta**(-zeta_order * disc_bilinear(*(gamma + delta))) for delta in singls]
                              + [sqrt2 * zeta**(-zeta_order * disc_bilinear(*(gamma + delta))) for delta in pairs]
                              for gamma in singls]
                           + [  [sqrt2 * zeta**(-zeta_order * disc_bilinear(*(gamma + delta))) for delta in singls]
                              + [zeta**(-zeta_order * disc_bilinear(*(gamma + delta))) + zeta**(-zeta_order * disc_bilinear(*(gamma + neg(delta)))) for delta in pairs]
                              for gamma in pairs] )
    else:
        Smat = zeta ** (zeta_order / 8 * L.dim()) / drt * matrix( CC, [  [  zeta**(-zeta_order * disc_bilinear(*(gamma + delta)))
                                 - zeta ** (-zeta_order * disc_bilinear(*(gamma + neg(delta))))
                                 for delta in pairs]
                                for gamma in pairs ])

    return (Smat, Tmat)


def _qr(mat):
    r"""
    Compute the R matrix in QR decomposition using Housholder reflections.

    INPUT:

    - ``mat`` -- A matrix over a complex interval field.

    OUTPUT:

    - A matrix in row echelon form.  The algorithm raises an
      ArithmeticError if the precision is insufficient to compute the
      decomposition.

    .. NOTE::

        This is an adoption of the implementation in mpmath by Andreas
        Strombergson.

    .. TODO::

        Naturally, this is a stub implementation.  It belongs into the
        corresponding matrix class, but for the infrastructure of the
        matrix class seems a bit too obscrure to me.  Can an expert help?
        Ask Rob Beezer?

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank_dimension import _qr
        sage: CI = ComplexIntervalField(100)
        sage: m = matrix(CI, 2, [2,1,1,2])
        sage: _qr(m)
        [2.23606797749978969640917366873?  1.7888543819998317571273389350?]
        [                               0 -1.3416407864998738178455042012?]
        sage: m = MatrixSpace(CI, 3,3).random_element()
        sage: mqr = _qr(m)
        sage: mqr[1,0] == 0 and mqr[2,0] == 0 and mqr[2,1] == 0
        True
        sage: mqr[2,2].contains_zero() or ((m * mqr.inverse()) * (m * mqr.inverse()).transpose().conjugate()).is_one()
        True

    TESTS::

        sage: _qr(matrix(CC, [[2]]))
        Traceback (most recent call last):
        ...
        TypeError: This QR decomposition is implemented only over complex inverval fields
    """
    from sage.rings.complex_interval_field import is_ComplexIntervalField

    CC = mat.base_ring()
    if not is_ComplexIntervalField(CC):
        raise TypeError("This QR decomposition is implemented only over complex inverval fields")

    mat = copy(mat)
    m = mat.nrows()
    n = mat.ncols()

    cur_row = 0
    for j in range(0, n):
        if all(mat[i, j].contains_zero() for i in xrange(cur_row + 1, m)):
            if not mat[cur_row, j].contains_zero():
                cur_row += 1
            continue

        s = sum((abs(mat[i, j])) ** 2 for i in xrange(cur_row, m))
        if s.contains_zero():
            raise ArithmeticError("Cannot handle sums of elements that are too imprecise")

        p = sqrt(s)
        if (s - p * mat[cur_row, j]).contains_zero():
            raise ArithmeticError("Cannot handle sums of elements that are too imprecise")
        kappa = 1 / (s - p * mat[cur_row, j])

        mat[cur_row, j] -= p
        for k in range(j + 1, n):
            y = sum(mat[i, j].conjugate() * mat[i, k]
                    for i in xrange(cur_row, m)) * kappa
            for i in range(cur_row, m):
                mat[i, k] -= mat[i, j] * y

        mat[cur_row, j] = p
        for i in range(cur_row + 1, m):
            mat[i, j] = CC(0)

        cur_row += 1

    return mat
