r"""
Unitary Groups `GU(n,q)` and `SU(n,q)`

These are `n \times n` unitary matrices with entries in
`GF(q^2)`.

EXAMPLES::

    sage: G = SU(3,5)
    sage: G.order()
    378000
    sage: G
    Special Unitary Group of degree 3 over Finite Field in a of size 5^2
    sage: G.gens()
    (
    [      a       0       0]  [4*a   4   1]
    [      0 2*a + 2       0]  [  4   4   0]
    [      0       0     3*a], [  1   0   0]
    )
    sage: G.base_ring()
    Finite Field in a of size 5^2

AUTHORS:

- David Joyner (2006-03): initial version, modified from
  special_linear (by W. Stein)

- David Joyner (2006-05): minor additions (examples, _latex_, __str__,
  gens)

- William Stein (2006-12): rewrite

- Volker Braun (2013-1) port to new Parent, libGAP, extreme refactoring.

- Sebastian Oehms (2018-8) add  ``_UG``,
  :meth:`~sage.groups.matrix_gps.unitary.UnitaryMatrixGroup_generic.invariant_form`,
  option for user defined invariant bilinear form, and bug-fix in
  ``_check_matrix`` (see :trac:`26028`)
"""

# ****************************************************************************
#       Copyright (C) 2006 David Joyner and William Stein
#       Copyright (C) 2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.all import GF
from sage.rings.finite_rings.finite_field_base import is_FiniteField
from sage.misc.latex import latex
from sage.misc.cachefunc import cached_method
from sage.groups.matrix_gps.named_group import (
    normalize_args_vectorspace, normalize_args_invariant_form,
    NamedMatrixGroup_generic, NamedMatrixGroup_gap )
from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap


def finite_field_sqrt(ring):
    """
    Helper function.

    INPUT:

    A ring.

    OUTPUT:

    Integer q such that ``ring`` is the finite field with `q^2` elements.

    EXAMPLES::

        sage: from sage.groups.matrix_gps.unitary import finite_field_sqrt
        sage: finite_field_sqrt(GF(4, 'a'))
        2
    """
    if not is_FiniteField(ring):
        raise ValueError('not a finite field')
    q, rem = ring.cardinality().sqrtrem()
    if rem:
        raise ValueError('cardinality not a square')
    return q


###############################################################################
# Unitary Group: common Code for both GU and SU
###############################################################################

def _UG(n, R, special, var='a', invariant_form=None):
    r"""
    This function is commonly used by the functions :func:`GU` and :func:`SU`
    to avoid duplicated code. For documentation and examples
    see the individual functions.

    TESTS::

        sage: GU(3,25).order()  # indirect doctest
        3961191000000
    """
    prefix = 'General'
    latex_prefix ='G'
    if special:
        prefix = 'Special'
        latex_prefix ='S'

    degree, ring = normalize_args_vectorspace(n, R, var=var)
    if is_FiniteField(ring):
        q = ring.cardinality()
        ring = GF(q**2, name=var)
        if invariant_form is not None:
            raise NotImplementedError("invariant_form for finite groups is fixed by GAP")

    if invariant_form is not None:
        invariant_form = normalize_args_invariant_form(ring, degree, invariant_form)
        if not invariant_form.is_hermitian():
            raise ValueError("invariant_form must be hermitian")

        inserted_text = 'with respect to hermitian form'
        try:
            if not invariant_form.is_positive_definite():
               inserted_text = 'with respect to non positive definite hermitian form'
        except:
            pass

        name = '{0} Unitary Group of degree {1} over {2} {3}\n{4}'.format(prefix,
                                 degree, ring, inserted_text,invariant_form)
        ltx  = r'\text{{{0}U}}_{{{1}}}({2})\text{{ {3} }}{4}'.format(latex_prefix,
                     degree, latex(ring), inserted_text, latex(invariant_form))
    else:
        name = '{0} Unitary Group of degree {1} over {2}'.format(prefix, degree, ring)
        ltx  = r'\text{{{0}U}}_{{{1}}}({2})'.format(latex_prefix, degree, latex(ring))

    if is_FiniteField(ring):
        cmd = '{0}U({1}, {2})'.format(latex_prefix, degree, q)
        return UnitaryMatrixGroup_gap(degree, ring, special, name, ltx, cmd)
    else:
        return UnitaryMatrixGroup_generic(degree, ring, special, name, ltx, invariant_form=invariant_form)



###############################################################################
# General Unitary Group
###############################################################################

def GU(n, R, var='a', invariant_form=None):
    r"""
    Return the general unitary group.

    The general unitary group `GU( d, R )` consists of all `d \times d`
    matrices that preserve a nondegenerate sesquilinear form over the
    ring `R`.

    .. NOTE::

        For a finite field the matrices that preserve a sesquilinear
        form over `F_q` live over `F_{q^2}`. So ``GU(n,q)`` for
        a prime power ``q`` constructs the matrix group over the base
        ring ``GF(q^2)``.

    .. NOTE::

        This group is also available via ``groups.matrix.GU()``.

    INPUT:

    - ``n`` -- a positive integer

    - ``R`` -- ring or an integer; if an integer is specified, the
      corresponding finite field is used

    - ``var`` -- (optional, default: ``'a'``) variable used to
      represent generator of the finite field, if needed

    - ``invariant_form`` -- (optional) instances being accepted by
      the matrix-constructor which define a `n \times n` square matrix
      over R describing the hermitian form to be kept invariant
      by the unitary group; the form is checked to be
      non-degenerate and hermitian but not to be positive definite

    OUTPUT:

    Return the general unitary group.

    EXAMPLES::

        sage: G = GU(3, 7); G
        General Unitary Group of degree 3 over Finite Field in a of size 7^2
        sage: G.gens()
        (
        [  a   0   0]  [6*a   6   1]
        [  0   1   0]  [  6   6   0]
        [  0   0 5*a], [  1   0   0]
        )
        sage: GU(2,QQ)
        General Unitary Group of degree 2 over Rational Field

        sage: G = GU(3, 5, var='beta')
        sage: G.base_ring()
        Finite Field in beta of size 5^2
        sage: G.gens()
        (
        [  beta      0      0]  [4*beta      4      1]
        [     0      1      0]  [     4      4      0]
        [     0      0 3*beta], [     1      0      0]
        )

    Using the ``invariant_form`` option::

        sage: UCF = UniversalCyclotomicField(); e5=UCF.gen(5)
        sage: m=matrix(UCF, 3,3, [[1,e5,0],[e5.conjugate(),2,0],[0,0,1]])
        sage: G  = GU(3, UCF)
        sage: Gm = GU(3, UCF, invariant_form=m)
        sage: G == Gm
        False
        sage: G.invariant_form()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: Gm.invariant_form()
        [     1   E(5)      0]
        [E(5)^4      2      0]
        [     0      0      1]
        sage: pm=Permutation((1,2,3)).to_matrix()
        sage: g = G(pm); g in G; g
        True
        [0 0 1]
        [1 0 0]
        [0 1 0]
        sage: Gm(pm)
        Traceback (most recent call last):
        ...
        TypeError: matrix must be unitary with respect to the hermitian form
        [     1   E(5)      0]
        [E(5)^4      2      0]
        [     0      0      1]

        sage: GU(3,3, invariant_form=[[1,0,0],[0,2,0],[0,0,1]])
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP

        sage: GU(2,QQ, invariant_form=[[1,0],[2,0]])
        Traceback (most recent call last):
        ...
        ValueError: invariant_form must be non-degenerate

    TESTS::

        sage: TestSuite(G).run()
        sage: groups.matrix.GU(2, 3)
        General Unitary Group of degree 2 over Finite Field in a of size 3^2
    """
    return _UG(n, R, False, var=var, invariant_form=invariant_form)



###############################################################################
# Special Unitary Group
###############################################################################

def SU(n, R, var='a', invariant_form=None):
    r"""
    The special unitary group `SU( d, R )` consists of all `d \times d`
    matrices that preserve a nondegenerate sesquilinear form over the
    ring `R` and have determinant `1`.

    .. NOTE::

        For a finite field the matrices that preserve a sesquilinear
        form over `F_q` live over `F_{q^2}`. So ``SU(n,q)`` for
        a prime power ``q`` constructs the matrix group over the base
        ring ``GF(q^2)``.

    .. NOTE::

        This group is also available via ``groups.matrix.SU()``.

    INPUT:

    - ``n`` -- a positive integer

    - ``R`` -- ring or an integer; if an integer is specified, the
      corresponding finite field is used

    - ``var`` -- (optional, default: ``'a'``) variable used to
      represent generator of the finite field, if needed

    - ``invariant_form`` -- (optional) instances being accepted by
      the matrix-constructor which define a `n \times n` square matrix
      over R describing the hermitian form to be kept invariant
      by the unitary group; the form is checked to be
      non-degenerate and hermitian but not to be positive definite

    OUTPUT:

    Return the special unitary group.

    EXAMPLES::

        sage: SU(3,5)
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: SU(3, GF(5))
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: SU(3,QQ)
        Special Unitary Group of degree 3 over Rational Field

    Using the ``invariant_form`` option::

        sage: CF3 = CyclotomicField(3); e3 = CF3.gen()
        sage: m=matrix(CF3, 3,3, [[1,e3,0],[e3.conjugate(),2,0],[0,0,1]])
        sage: G  = SU(3, CF3)
        sage: Gm = SU(3, CF3, invariant_form=m)
        sage: G == Gm
        False
        sage: G.invariant_form()
        [1 0 0]
        [0 1 0]
        [0 0 1]
        sage: Gm.invariant_form()
        [         1      zeta3          0]
        [-zeta3 - 1          2          0]
        [         0          0          1]
        sage: pm=Permutation((1,2,3)).to_matrix()
        sage: G(pm)
        [0 0 1]
        [1 0 0]
        [0 1 0]
        sage: Gm(pm)
        Traceback (most recent call last):
        ...
        TypeError: matrix must be unitary with respect to the hermitian form
        [         1      zeta3          0]
        [-zeta3 - 1          2          0]
        [         0          0          1]

        sage: SU(3,5, invariant_form=[[1,0,0],[0,2,0],[0,0,3]])
        Traceback (most recent call last):
        ...
        NotImplementedError: invariant_form for finite groups is fixed by GAP

    TESTS::

        sage: TestSuite(Gm).run()
        sage: groups.matrix.SU(2, 3)
        Special Unitary Group of degree 2 over Finite Field in a of size 3^2
    """
    return _UG(n, R, True, var=var, invariant_form=invariant_form)



########################################################################
# Unitary Group class
########################################################################

class UnitaryMatrixGroup_generic(NamedMatrixGroup_generic):
    r"""
    General Unitary Group over arbitrary rings.

    EXAMPLES::

        sage: G = GU(3, GF(7)); G
        General Unitary Group of degree 3 over Finite Field in a of size 7^2
        sage: latex(G)
        \text{GU}_{3}(\Bold{F}_{7^{2}})

        sage: G = SU(3, GF(5));  G
        Special Unitary Group of degree 3 over Finite Field in a of size 5^2
        sage: latex(G)
        \text{SU}_{3}(\Bold{F}_{5^{2}})

        sage: CF3 = CyclotomicField(3); e3 = CF3.gen()
        sage: m=matrix(CF3, 3,3, [[1,e3,0],[e3.conjugate(),2,0],[0,0,1]])
        sage: G = SU(3, CF3, invariant_form=m)
        sage: latex(G)
        \text{SU}_{3}(\Bold{Q}(\zeta_{3}))\text{ with respect to hermitian form }\left(\begin{array}{rrr}
        1 & \zeta_{3} & 0 \\
        -\zeta_{3} - 1 & 2 & 0 \\
        0 & 0 & 1
        \end{array}\right)
    """

    @cached_method
    def invariant_form(self):
        """
        Return the hermitian form preserved by the unitary
        group.

        OUTPUT:

        A square matrix describing the bilinear form

        EXAMPLES::

            sage: SU4 = SU(4,QQ)
            sage: SU4.invariant_form()
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
        """
        if self._invariant_form is not None:
            return self._invariant_form

        from sage.matrix.constructor import identity_matrix
        m = identity_matrix(self.base_ring(), self.degree())
        m.set_immutable()
        return m


    def _check_matrix(self, x, *args):
        """
        Check whether the matrix ``x`` is unitary.

        See :meth:`~sage.groups.matrix_gps.matrix_group._check_matrix`
        for details.

        EXAMPLES::

            sage: G = GU(2, GF(5))
            sage: G._check_matrix(G.an_element().matrix())
            sage: G = SU(2, GF(5))
            sage: G._check_matrix(G.an_element().matrix())
        """
        if self._special and x.determinant() != 1:
            raise TypeError('matrix must have determinant one')

        H = self.invariant_form()
        if x * H * x.conjugate_transpose() != H:
            if H == self.one().matrix():
                raise TypeError('matrix must be unitary')
            else:
                raise TypeError('matrix must be unitary with respect to the hermitian form\n{}'.format(H))

class UnitaryMatrixGroup_gap(UnitaryMatrixGroup_generic, NamedMatrixGroup_gap, FinitelyGeneratedMatrixGroup_gap):
    r"""
    The general or special unitary group in GAP.

    TESTS:

    Check that :trac:`20867` is fixed::

        sage: from sage.groups.matrix_gps.finitely_generated import FinitelyGeneratedMatrixGroup_gap
        sage: G = GU(3,3)
        sage: isinstance(G, FinitelyGeneratedMatrixGroup_gap)
        True
    """

    @cached_method
    def invariant_form(self):
        """
        Return the hermitian form preserved by the unitary group.

        OUTPUT:

        A square matrix describing the bilinear form

        EXAMPLES::

            sage: G32=GU(3,2)
            sage: G32.invariant_form()
            [0 0 1]
            [0 1 0]
            [1 0 0]
        """
        d = self.degree()
        R = self.base_ring()
        # note that self.gap().InvariantSesquilinearForm()['matrix'].matrix().base_ring() != R for example for self = GU(3.2)
        # therefore we have to coerce into the right matrix space
        from sage.matrix.constructor import matrix
        m = matrix(R, d, d, self.gap().InvariantSesquilinearForm()['matrix'].matrix())
        m.set_immutable()
        return m

