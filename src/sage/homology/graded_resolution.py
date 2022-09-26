r"""
Graded free resolutions

Let `R` be a commutative ring. A graded free resolution of a graded
`R`-module `M` is a :mod:`free resolution <sage.homology.free_resolution>`
such that all maps are homogeneous module homomorphisms.

EXAMPLES::

    sage: S.<x,y,z,w> = PolynomialRing(QQ)
    sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
    sage: r = I.graded_free_resolution(algorithm='minimal')
    sage: r
    S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
    sage: I.graded_free_resolution(algorithm='shreyer')
    S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
    sage: I.graded_free_resolution(algorithm='standard')
    S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
    sage: I.graded_free_resolution(algorithm='heuristic')
    S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0

::

    sage: d = r.differential(2)
    sage: d
    Free module morphism defined as left-multiplication by the matrix
    [ y  x]
    [-z -y]
    [ w  z]
    Domain: Ambient free module of rank 2 over the integral domain Multivariate Polynomial Ring
    in x, y, z, w over Rational Field
    Codomain: Ambient free module of rank 3 over the integral domain Multivariate Polynomial Ring
    in x, y, z, w over Rational Field
    sage: d.image()
    Submodule of Ambient free module of rank 3 over the integral domain Multivariate Polynomial Ring
    in x, y, z, w over Rational Field
    Generated by the rows of the matrix:
    [ y -z  w]
    [ x -y  z]
    sage: m = d.image()
    sage: m.graded_free_resolution(shifts=(2,2,2))
    S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0

An example of multigraded resolution from Example 9.1 of [MilStu2005]_::

    sage: R.<s,t> = QQ[]
    sage: S.<a,b,c,d> = QQ[]
    sage: phi = S.hom([s, s*t, s*t^2, s*t^3])
    sage: I = phi.kernel(); I
    Ideal (c^2 - b*d, b*c - a*d, b^2 - a*c) of Multivariate Polynomial Ring in a, b, c, d over Rational Field
    sage: P3 = ProjectiveSpace(S)
    sage: C = P3.subscheme(I)  # twisted cubic curve
    sage: r = I.graded_free_resolution(degrees=[(1,0), (1,1), (1,2), (1,3)])
    sage: r
    S((0, 0)) <-- S((-2, -4))⊕S((-2, -3))⊕S((-2, -2)) <-- S((-3, -5))⊕S((-3, -4)) <-- 0
    sage: r.K_polynomial(names='s,t')
    s^3*t^5 + s^3*t^4 - s^2*t^4 - s^2*t^3 - s^2*t^2 + 1

AUTHORS:

- Kwankyu Lee (2022-05): initial version
- Travis Scrimshaw (2022-08-23): refactored for free module inputs
"""

# ****************************************************************************
#       Copyright (C) 2022 Kwankyu Lee <ekwankyu@gmail.com>
#                 (C) 2022 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.libs.singular.singular import si2sa_resolution_graded
from sage.libs.singular.function import singular_function
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.element import Matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import Module_free_ambient
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.rings.ideal import Ideal_generic

from sage.homology.free_resolution import (FiniteFreeResolution,
                                           FiniteFreeResolution_free_module,
                                           FiniteFreeResolution_singular)

class GradedFiniteFreeResolution(FiniteFreeResolution):
    r"""
    Graded finite free resolutions.

    INPUT:

    - ``module`` -- a homogeneous submodule of a free module `M` of rank `n`
      over `S` or a homogeneous ideal of a multivariate polynomial ring `S`
    - ``degrees`` -- (default: a list with all entries `1`) a list of integers
      or integer vectors giving degrees of variables of `S`
    - ``shifts`` -- a list of integers or integer vectors giving shifts of
      degrees of `n` summands of the free module `M`; this is a list of zero
      degrees of length `n` by default
    - ``name`` -- a string; name of the base ring

    .. WARNING::

        This does not check that the module is homogeneous.
    """
    def __init__(self, module, degrees=None, shifts=None, name='S', **kwds):
        r"""
        Initialize ``self``.

        TESTS::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: TestSuite(r).run(skip=['_test_pickling'])

        An overdetermined system over a PID::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: M = matrix([[x^2, 2*x^2],
            ....:             [3*x^2, 5*x^2],
            ....:             [5*x^2, 4*x^2]])
            sage: res = FreeResolution(M, graded=True)
            sage: res
            S(0)⊕S(0) <-- S(-2)⊕S(-2) <-- 0
            sage: res._res_shifts
            [[2, 2]]
        """
        super().__init__(module, name=name, **kwds)

        nvars = self._base_ring.ngens()

        if degrees is None:
            degrees = nvars * (1,)  # standard grading

        if len(degrees) != nvars:
            raise ValueError('the length of degrees does not match the number of generators')

        if degrees[0] in ZZ:
            zero_deg = 0
            multigrade = False
        else: # degrees are integer vectors
            degrees = tuple([vector(v) for v in degrees])
            zero_deg = degrees[0].parent().zero()
            multigrade = True

        if shifts is None:
            if isinstance(self._module, Ideal_generic):
                rank = 1
            elif isinstance(self._module, Module_free_ambient):
                rank = self._m().nrows()
            elif isinstance(self._module, Matrix):
                rank = self._module.ncols()

            shifts = rank * [zero_deg]

        self._shifts = shifts
        self._degrees = tuple(degrees)
        self._multigrade = multigrade
        self._zero_deg = zero_deg

    def _repr_module(self, i):
        """
        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r._repr_module(0)
            'S(0)'
            sage: r._repr_module(1)
            'S(-2)⊕S(-2)⊕S(-2)'
            sage: r._repr_module(2)
            'S(-3)⊕S(-3)'
            sage: r._repr_module(3)
            '0'

            sage: r = I.graded_free_resolution(shifts=[-1])
            sage: r._repr_module(0)
            'S(1)'
        """
        self._maps  # to set _res_shifts
        if i > len(self):
            return '0'

        if i == 0:
            shifts = self._shifts
        else:
            shifts = self._res_shifts[i - 1]

        if not shifts:
            return '0'

        return '\u2295'.join(f'{self._name}' + '({})'.format(-sh)
                             for sh in shifts)

    def shifts(self, i):
        r"""
        Return the shifts of ``self``.

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r.shifts(0)
            [0]
            sage: r.shifts(1)
            [2, 2, 2]
            sage: r.shifts(2)
            [3, 3]
            sage: r.shifts(3)
            []
        """
        if i < 0:
            raise IndexError('invalid index')
        elif i == 0:
            shifts = self._shifts
        elif i > len(self):
            shifts = []
        else:
            self._maps  # to set _res_shifts
            shifts = self._res_shifts[i - 1]

        return shifts

    def betti(self, i, a=None):
        r"""
        Return the `i`-th Betti number in degree `a`.

        INPUT:

        - ``i`` -- nonnegative integer

        - ``a`` -- a degree; if ``None``, return Betti numbers in all degrees

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r.betti(0)
            {0: 1}
            sage: r.betti(1)
            {2: 3}
            sage: r.betti(2)
            {3: 2}
            sage: r.betti(1, 0)
            0
            sage: r.betti(1, 1)
            0
            sage: r.betti(1, 2)
            3
        """
        shifts = self.shifts(i)

        if a is None:
            degrees = shifts
        else:
            degrees = [a]

        betti = {}
        for s in degrees:
            betti[s] = len([d for d in shifts if d == s])

        if a is None:
            return betti
        else:
            return betti[a] if a in betti else 0

    def K_polynomial(self, names=None):
        r"""
        Return the K-polynomial of this resolution.

        INPUT:

        - ``names`` -- (optional) a string of names of the variables
          of the K-polynomial

        EXAMPLES::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r.K_polynomial()
            2*t^3 - 3*t^2 + 1
        """
        if self._multigrade:
            n = self._degrees[0].degree()
        else:
            n = 1

        if names is not None:
            L = LaurentPolynomialRing(ZZ, names=names)
        else:
            L = LaurentPolynomialRing(ZZ, 't', n)

        kpoly = 1
        sign = -1
        self._maps  # to set _res_shifts
        for j in range(len(self)):
            for v in self._res_shifts[j]:
                if self._multigrade:
                    kpoly += sign * L.monomial(*list(v))
                else:
                    kpoly += sign * L.monomial(v)
            sign = -sign

        return kpoly


class GradedFiniteFreeResolution_free_module(GradedFiniteFreeResolution, FiniteFreeResolution_free_module):
    r"""
    Graded free resolution of free modules.

    .. WARNING::

        This does not check that the module is homogeneous.

    EXAMPLES::

        sage: from sage.homology.free_resolution import FreeResolution
        sage: R.<x> = QQ[]
        sage: M = matrix([[x^3, 3*x^3, 5*x^3],
        ....:             [0, x, 2*x]])
        sage: res = FreeResolution(M, graded=True)
        sage: res
        S(0)⊕S(0)⊕S(0) <-- S(-3)⊕S(-1) <-- 0
    """
    def __init__(self, module, degrees=None, *args, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: R.<x> = QQ[]
            sage: M = matrix([[x^3, 3*x^3, 5*x^3],
            ....:             [0, x, 2*x]])
            sage: res = FreeResolution(M, graded=True)
            sage: TestSuite(res).run(skip="_test_pickling")
        """
        super().__init__(module, degrees=degrees, *args, **kwds)

        if len(self._degrees) > 1 and any(d != 1 for d in self._degrees):
            raise NotImplementedError("only the natural grading supported "
                                      "when more than one generator")

    @lazy_attribute
    def _maps(self):
        r"""
        The maps that define ``self``.

        This also sets the attribute ``_res_shifts``.

        TESTS::

            sage: from sage.homology.free_resolution import FreeResolution
            sage: R.<x> = QQ[]
            sage: M = matrix([[x^3, 3*x^3, 5*x^3],
            ....:             [0, x, 2*x]])
            sage: res = FreeResolution(M, graded=True)
            sage: res
            S(0)⊕S(0)⊕S(0) <-- S(-3)⊕S(-1) <-- 0
            sage: res._maps
            [
            [  x^3     0]
            [3*x^3     x]
            [5*x^3   2*x]
            ]
            sage: res._res_shifts
            [[3, 1]]

            sage: I = R.ideal([x^4])
            sage: res = I.graded_free_resolution(shifts=[1], degrees=[2])
            sage: res
            S(-1) <-- S(-9) <-- 0
            sage: res._maps
            [[x^4]]
            sage: res._res_shifts
            [[9]]
        """
        def compute_degree(base, i):
            """
            Compute the degree by ``base * deg + shift``,
            where ``*`` is entry-wise multiplication, ``deg`` and
            ``shift`` are the ``i``-th index.
            """
            deg = self._degrees[0]
            shift = self._shifts[i]
            if self._multigrade:
                return vector([val * d + s for val, d, s in zip(base, deg, shift)])
            return base * deg + shift

        if isinstance(self._module, Ideal_generic):
            from sage.matrix.constructor import matrix
            val = self._module.gen(0)
            self._res_shifts = [[compute_degree(val.degree(), 0)]]
            return [matrix([[val]])]

        M = self._m()

        def find_deg(i):
            for j in range(M.nrows()):
                ret = M[j,i].degree()
                if ret != -1:
                    return ret
            raise NotImplementedError("a generator maps to 0")

        self._res_shifts = [[compute_degree(find_deg(i), i)
                             for i in range(M.ncols())]]
        return [M]


class GradedFiniteFreeResolution_singular(GradedFiniteFreeResolution, FiniteFreeResolution_singular):
    r"""
    Graded free resolutions of submodules and ideals of multivariate
    polynomial rings implemented using Singular.

    INPUT:

    - ``module`` -- a homogeneous submodule of a free module `M` of rank `n`
      over `S` or a homogeneous ideal of a multivariate polynomial ring `S`

    - ``degrees`` -- (default: a list with all entries `1`) a list of integers
      or integer vectors giving degrees of variables of `S`

    - ``shifts`` -- a list of integers or integer vectors giving shifts of
      degrees of `n` summands of the free module `M`; this is a list of zero
      degrees of length `n` by default

    - ``name`` -- a string; name of the base ring

    - ``algorithm`` -- Singular algorithm to compute a resolution of ``ideal``

    If ``module`` is an ideal of `S`, it is considered as a submodule of a
    free module of rank `1` over `S`.

    The degrees given to the variables of `S` are integers or integer vectors of
    the same length. In the latter case, `S` is said to be multigraded, and the
    resolution is a multigraded free resolution. The standard grading where all
    variables have degree `1` is used if the degrees are not specified.

    A summand of the graded free module `M` is a shifted (or twisted) module of
    rank one over `S`, denoted `S(-d)` with shift `d`.

    The computation of the resolution is done by using ``libSingular``.
    Different Singular algorithms can be chosen for best performance.

    OUTPUT: a graded minimal free resolution of ``ideal``

    The available algorithms and the corresponding Singular commands are shown
    below:

        ============= ============================
        algorithm     Singular commands
        ============= ============================
        ``minimal``   ``mres(ideal)``
        ``shreyer``   ``minres(sres(std(ideal)))``
        ``standard``  ``minres(nres(std(ideal)))``
        ``heuristic`` ``minres(res(std(ideal)))``
        ============= ============================

    .. WARNING::

        This does not check that the module is homogeneous.

    EXAMPLES::

        sage: S.<x,y,z,w> = PolynomialRing(QQ)
        sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
        sage: r = I.graded_free_resolution()
        sage: r
        S(0) <-- S(-2)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3) <-- 0
        sage: len(r)
        2

        sage: I = S.ideal([z^2 - y*w, y*z - x*w, y - x])
        sage: I.is_homogeneous()
        True
        sage: r = I.graded_free_resolution()
        sage: r
        S(0) <-- S(-1)⊕S(-2)⊕S(-2) <-- S(-3)⊕S(-3)⊕S(-4) <-- S(-5) <-- 0
    """
    def __init__(self, module, degrees=None, shifts=None, name='S', algorithm='heuristic', **kwds):
        """
        Initialize.

        TESTS::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: TestSuite(r).run(skip=['_test_pickling'])
        """
        super().__init__(module, degrees=degrees, shifts=shifts, name=name, **kwds)
        self._algorithm = algorithm

    @lazy_attribute
    def _maps(self):
        """
        The maps that define ``self``.

        This also sets the attribute ``_res_shifts``.

        TESTS::

            sage: S.<x,y,z,w> = PolynomialRing(QQ)
            sage: I = S.ideal([y*w - z^2, -x*w + y*z, x*z - y^2])
            sage: r = I.graded_free_resolution()
            sage: r._maps
            [
                                             [-y  x]
                                             [ z -y]
            [z^2 - y*w y*z - x*w y^2 - x*z], [-w  z]
            ]
            sage: r._res_shifts
            [[2, 2, 2], [3, 3]]
        """
        #cdef int i, j, k, ncols, nrows
        #cdef list res_shifts, prev_shifts, new_shifts

        # This ensures the first component of the Singular resolution to be a
        # module, like the later components. This is important when the
        # components are converted to Sage modules.
        module = singular_function("module")
        mod = module(self._m())

        if self._algorithm == 'minimal':
            mres = singular_function('mres')  # syzygy method
            r = mres(mod, 0)
        elif self._algorithm == 'shreyer':
            std = singular_function('std')
            sres = singular_function('sres')  # Shreyer method
            minres = singular_function('minres')
            r = minres(sres(std(mod), 0))
        elif self._algorithm == 'standard':
            nres = singular_function('nres')  # standard basis method
            minres = singular_function('minres')
            r = minres(nres(mod, 0))
        elif self._algorithm == 'heuristic':
            std = singular_function('std')
            res = singular_function('res')    # heuristic method
            minres = singular_function('minres')
            r = minres(res(std(mod), 0))

        res_mats, res_degs = si2sa_resolution_graded(r, self._degrees)

        # compute shifts of free modules in the resolution
        res_shifts = []
        prev_shifts = list(self._shifts)
        for k in range(len(res_degs)):
            new_shifts = []
            degs = res_degs[k]
            ncols = len(degs)
            for j in range(ncols):
                col = degs[j]
                nrows = len(col)
                # should be enough to compute the new shifts
                # from any one entry of the column vector
                for i in range(nrows):
                    d = col[i]
                    if d is not None:
                        e = prev_shifts[i]
                        new_shifts.append(d + e)
                        break
            res_shifts.append(new_shifts)
            prev_shifts = new_shifts

        self._res_shifts = res_shifts
        return res_mats

