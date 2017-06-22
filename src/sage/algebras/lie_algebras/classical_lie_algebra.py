"""
Classical Lie Algebras

These are the Lie algebras corresponding to types `A_n`, `B_n`, `C_n`,
and `D_n`. We also include support for the exceptional types
`E_{6,7,8}`, `F_4`, and `G_2` in the Chevalley basis, and we
give the matrix representation given in [HRT2000]_.

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2013-2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.misc import repr_lincomb
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import RingElement
from sage.structure.indexed_generators import IndexedGenerators
from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.algebra import Algebra
from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra, LieAlgebraFromAssociative
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.algebras.lie_algebras.structure_coefficients import LieAlgebraWithStructureCoefficients
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.cartan_matrix import CartanMatrix
from sage.combinat.root_system.dynkin_diagram import DynkinDiagram_class
from sage.combinat.free_module import CombinatorialFreeModule
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.arith import binomial
from sage.sets.family import Family

class ClassicalMatrixLieAlgebra(LieAlgebraFromAssociative):
    """
    A classical Lie algebra represented using matrices.

    INPUT:

    - ``R`` -- the base ring
    - ``ct`` -- the finite Cartan type

    EXAMPLES::

        sage: lie_algebras.ClassicalMatrix(QQ, ['A', 4])
        Special linear Lie algebra of rank 5 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, CartanType(['B',4]))
        Special orthogonal Lie algebra of rank 9 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, 'C4')
        Symplectic Lie algebra of rank 8 over Rational Field
        sage: lie_algebras.ClassicalMatrix(QQ, cartan_type=['D',4])
        Special orthogonal Lie algebra of rank 8 over Rational Field
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Return the correct parent based on input.

        EXAMPLES::

            sage: lie_algebras.ClassicalMatrix(QQ, ['A', 4])
            Special linear Lie algebra of rank 5 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, CartanType(['B',4]))
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, 'C4')
            Symplectic Lie algebra of rank 8 over Rational Field
            sage: lie_algebras.ClassicalMatrix(QQ, cartan_type=['D',4])
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        if isinstance(cartan_type, (CartanMatrix, DynkinDiagram_class)):
            cartan_type = cartan_type.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)

        if not cartan_type.is_finite():
            raise ValueError("only for finite types")

        if cartan_type.type() == 'A':
            return sl(R, cartan_type.rank() + 1)
        if cartan_type.type() == 'B':
            return so(R, 2*cartan_type.rank() + 1)
        if cartan_type.type() == 'C':
            return sp(R, 2*cartan_type.rank())
        if cartan_type.type() == 'D':
            return so(R, 2*cartan_type.rank())
        if cartan_type.type() == 'E':
            if cartan_type.rank() == 6:
                return e6(R)
            if cartan_type.rank() in [7,8]:
                raise NotImplementedError("not yet implemented")
        if cartan_type.type() == 'F' and cartan_type.rank() == 4:
            return f4(R)
        if cartan_type.type() == 'G' and cartan_type.rank() == 2:
            return g2(R)
        raise ValueError("invalid Cartan type")

    def __init__(self, R, ct, e, f, h):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: TestSuite(g).run()

        TESTS:

        Check that :trac:`23266` is fixed::

            sage: sl2 = lie_algebras.sl(QQ, 2, 'matrix')
            sage: isinstance(sl2.indices(), FiniteEnumeratedSet)
            True
        """
        n = len(e)
        names = ['e%s'%i for i in range(1, n+1)]
        names += ['f%s'%i for i in range(1, n+1)]
        names += ['h%s'%i for i in range(1, n+1)]
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(names)
        LieAlgebraFromAssociative.__init__(self, e[0].parent(),
                                           gens=tuple(e + f + h),
                                           names=tuple(names),
                                           index_set=index_set,
                                           category=category)
        self._cartan_type = ct

        gens = tuple(self.gens())
        i_set = ct.index_set()
        self._e = Family(dict( (i, gens[c]) for c,i in enumerate(i_set) ))
        self._f = Family(dict( (i, gens[n+c]) for c,i in enumerate(i_set) ))
        self._h = Family(dict( (i, gens[2*n+c]) for c,i in enumerate(i_set) ))

    def e(self, i):
        r"""
        Return the generator `e_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.e(2)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._e[i]

    def f(self, i):
        r"""
        Return the generator `f_i`.-

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.f(2)
            [0 0 0]
            [0 0 0]
            [0 1 0]
        """
        return self._f[i]

    def h(self, i):
        """
        Return the generator `h_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.h(2)
            [ 0  0  0]
            [ 0  1  0]
            [ 0  0 -1]
        """
        return self._h[i]

    @cached_method
    def index_set(self):
        """
        Return the index_set of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.index_set()
            (1, 2)
        """
        return self._cartan_type.index_set()

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def epsilon(self, i, h):
        r"""
        Return the action of the functional
        `\varepsilon_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.epsilon(1, g.h(1))
            1
            sage: g.epsilon(2, g.h(1))
            -1
            sage: g.epsilon(3, g.h(1))
            0
        """
        return h[i-1,i-1]

    # Do we want this to be optional or requried?
    # There probably is a generic implementation we can do.
    @abstract_method(optional=True)
    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.simple_root(1, g.h(1))
            2
            sage: g.simple_root(1, g.h(2))
            -1
        """

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.
        If ``pos`` is ``True``, then returns `e_{\theta}`, otherwise it
        returns `f_{\theta}`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3, representation='matrix')
            sage: g.highest_root_basis_elt()
            [0 0 1]
            [0 0 0]
            [0 0 0]
        """
        RL = self._cartan_type.root_system().root_lattice()
        coroots = RL.simple_coroots()
        theta = RL.highest_root()
        i,w = theta.to_simple_root(True)
        r = RL.simple_root(i)
        if pos:
            gens = self._e
        else:
            gens = self._f
        cur = gens[i]
        for j in reversed(w):
            for k in range(-r.scalar(coroots[j])):
                cur = self.bracket(gens[j], cur)
            r = r.reflection(coroots[j], True)
        return cur

    @cached_method
    def basis(self):
        """
        Return a basis of ``self``.

        EXAMPLES::

            sage: M = LieAlgebra(ZZ, cartan_type=['A',2], representation='matrix')
            sage: list(M.basis())
            [
            [ 1  0  0]  [0 1 0]  [0 0 1]  [0 0 0]  [ 0  0  0]  [0 0 0]  [0 0 0]
            [ 0  0  0]  [0 0 0]  [0 0 0]  [1 0 0]  [ 0  1  0]  [0 0 1]  [0 0 0]
            [ 0  0 -1], [0 0 0], [0 0 0], [0 0 0], [ 0  0 -1], [0 0 0], [1 0 0],
            <BLANKLINE>
            [0 0 0]
            [0 0 0]
            [0 1 0]
            ]
        """
        # This is a fairly generic method of constructing a basis
        from sage.matrix.constructor import matrix

        R = self.base_ring()
        basis = list(self.lie_algebra_generators())
        expanded = True
        while expanded:
            expanded = False
            mat = []
            for i,x in enumerate(basis):
                mat.append(x.value.list())
                for y in basis[i+1:]:
                    mat.append(x.bracket(y).value.list())
            mat = matrix(R, mat)
            mat.echelonize()
            if mat.rank() != len(basis):
                basis = [self.element_class( self, self._assoc(mat[i].list()) )
                         for i in range(mat.rank())]
                expanded = True
        return Family(basis)

    def affine(self, kac_moody=False):
        """
        Return the affine (Kac-Moody) Lie algebra of ``self``.

        EXAMPLES::

            sage: so5 = lie_algebras.so(QQ, 5, 'matrix')
            sage: so5
            Special orthogonal Lie algebra of rank 5 over Rational Field
            sage: so5.affine()
            Affine Special orthogonal Kac-Moody algebra of rank 5 over Rational Field
        """
        from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
        return AffineLieAlgebra(self, kac_moody)

class gl(LieAlgebraFromAssociative):
    r"""
    The matrix Lie algebra `\mathfrak{gl}_n`.

    The Lie algebra `\mathfrak{gl}_n` which consists of all `n \times n`
    matrices.

    INPUT:

    - ``R`` -- the base ring
    - ``n`` -- the size of the matrix
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: TestSuite(g).run()

        TESTS:

        Check that :trac:`23266` is fixed::

            sage: gl2 = lie_algebras.gl(QQ, 2)
            sage: isinstance(gl2.basis().keys(), FiniteEnumeratedSet)
            True
            sage: Ugl2 = gl2.pbw_basis()
            sage: prod(Ugl2.gens())
            PBW['E_0_0']*PBW['E_0_1']*PBW['E_1_0']*PBW['E_1_1']
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        names = []
        gens = []
        for i in range(n):
            for j in range(n):
                names.append('E_{0}_{1}'.format(i,j))
                mat = MS({(i,j):one})
                mat.set_immutable()
                gens.append(mat)
        self._n = n
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(names)
        LieAlgebraFromAssociative.__init__(self, MS, tuple(gens),
                                           names=tuple(names),
                                           index_set=index_set,
                                           category=category)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.gl(QQ, 4)
            General linear Lie algebra of rank 4 over Rational Field
        """
        return "General linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{gl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy) - 2 \mathrm{tr}(x)
            \mathrm{tr}(y).

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: x = g.an_element()
            sage: y = g.gens()[1]
            sage: g.killing_form(x, y)
            8
        """
        return (2 * self._n * (x.value * y.value).trace()
                - 2 * x.value.trace() * y.value.trace())

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 2)
            sage: tuple(g.basis())
            (
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            )
        """
        G = self.gens()
        return Family(self._indices, lambda i: G[self._indices.index(i)])

    def monomial(self, i):
        r"""
        Return the basis element indexed by ``i``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: gl4 = lie_algebras.gl(QQ, 4)
            sage: gl4.monomial('E_2_1')
            [0 0 0 0]
            [0 0 0 0]
            [0 1 0 0]
            [0 0 0 0]
            sage: gl4.monomial((2,1))
            [0 0 0 0]
            [0 0 0 0]
            [0 1 0 0]
            [0 0 0 0]
        """
        if isinstance(i, tuple):
            return self.basis()['E_{}_{}'.format(*i)]
        return self.basis()[i]

    class Element(LieAlgebraFromAssociative.Element):
        def monomial_coefficients(self, copy=True):
            r"""
            Return the monomial coefficients of ``self``.

            EXAMPLES::

                sage: gl4 = lie_algebras.gl(QQ, 4)
                sage: x = gl4.monomial('E_2_1') + 3*gl4.monomial('E_0_3')
                sage: x.monomial_coefficients()
                {'E_0_3': 3, 'E_2_1': 1}
            """
            d = {}
            for k in self.value.nonzero_positions():
                d['E_{}_{}'.format(*k)] = self.value[k]
            return d

class sl(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{sl}_n`.

    The Lie algebra `\mathfrak{sl}_n`, which consists of all `n \times n`
    matrices with trace 0. This is the Lie algebra of type `A_{n-1}`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        e = [MS({(i,i+1):one}) for i in range(n-1)]
        f = [MS({(i+1,i):one}) for i in range(n-1)]
        h = [MS({(i,i):one, (i+1,i+1):-one}) for i in range(n-1)]
        self._n = n
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['A', n-1]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sl(QQ, 5, representation='matrix')
            Special linear Lie algebra of rank 5 over Rational Field
        """
        return "Special linear Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sl}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            10
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        return h[i,i] - h[i+1,i+1]

class so(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{so}_n`.

    The Lie algebra `\mathfrak{so}_n`, which consists of all real
    anti-symmetric `n \times n` matrices. This is the Lie algebra of
    type `B_{(n-1)/2}` or `D_{n/2}` if `n` is odd or even respectively.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        self._n = n
        if n % 2 == 0: # Even
            m = n / 2 - 1 # -1 for indexing
            n -= 1
            e = [MS({(m-1,n):one, (m,n-1):-one})]
            f = [MS({(n,m-1):one, (n-1,m):-one})]
            h = [MS({(m-1,m-1):one, (m,m):one, (n-1,n-1):-one, (n,n):-one})]
            m += 1
            ct = CartanType(['D', m])
        else: # Odd
            m = (n-1) / 2 - 1 # -1 for indexing
            n -= 1
            e = [MS({(m,n):2, (n,n-1):-2})]
            f = [MS({(n,m):one, (n-1,n):-one})]
            h = [MS({(m,m):2, (n-1,n-1):-2})]
            m += 1
            ct = CartanType(['B', m])
        e = [MS({(i,i+1):one, (m+i+1,m+i):-one}) for i in range(m-1)] + e
        f = [MS({(i+1,i):one, (m+i,m+i+1):-one}) for i in range(m-1)] + f
        h = [MS({(i,i):one, (i+1,i+1):-one, (m+i,m+i):-one, (m+i+1,m+i+1):one}) for i in range(m-1)] + h
        ClassicalMatrixLieAlgebra.__init__(self, R, ct, e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['B', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: LieAlgebra(QQ, cartan_type=['D', 4], representation='matrix')
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        return "Special orthogonal Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{so}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (n - 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            12
            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            14
        """
        return (self._n - 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES:

        The even or type `D` case::

            sage: g = lie_algebras.so(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1 -1]
            [ 0 -1  2  0]
            [ 0 -1  0  2]

        The odd or type `B` case::

            sage: g = lie_algebras.so(QQ, 9, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -1]
            [ 0  0 -2  2]
        """
        i = self.index_set().index(i)
        if i == len(self.index_set()) - 1:
            if self._n % 2 == 0:
                return h[i-1,i-1] + h[i,i]
            # otherwise we are odd
            return h[i,i]
        return h[i,i] - h[i+1,i+1]

class sp(ClassicalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{sp}_n`.

    The Lie algebra `\mathfrak{sp}_{2k}`, which consists of all
    `2k \times 2k` matrices `X` that satisfy the equation:

    .. MATH::

        X^T M - M X = 0

    where

    .. MATH::

        M = \begin{pmatrix}
        0 & I_k \\
        -I_k & 0
        \end{pmatrix}.

    This is the Lie algebra of type `C_k`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, n, sparse=True)
        one = R.one()
        self._n = n
        n = n // 2
        e = [MS({(i,i+1):one, (n+i+1,n+i):-one}) for i in range(n-1)]
        e.append(MS({(n-1,2*n-1):one})) # -1 for indexing
        f = [MS({(i+1,i):one, (n+i,n+i+1):-one}) for i in range(n-1)]
        f.append(MS({(2*n-1,n-1):one})) # -1 for indexing
        h = [MS({(i,i):one, (i+1,i+1):-one, (n+i,n+i):-one, (n+i+1,n+i+1):one}) for i in range(n-1)]
        h.append(MS({(n-1,n-1):one, (2*n-1,2*n-1):-one})) # -1 for indexing
        ClassicalMatrixLieAlgebra.__init__(self, R, CartanType(['C', n]), e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: lie_algebras.sp(QQ, 8, representation='matrix')
            Symplectic Lie algebra of rank 8 over Rational Field
        """
        return "Symplectic Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{sp}_n` is:

        .. MATH::

            \langle x \mid y \rangle = (2n + 2) \mathrm{tr}(xy).

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: x = g.an_element()
            sage: y = g.lie_algebra_generators()['e1']
            sage: g.killing_form(x, y)
            36
        """
        return (2 * self._n + 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8, representation='matrix')
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1  0]
            [ 0 -1  2 -2]
            [ 0  0 -1  2]
        """
        i = self.index_set().index(i)
        if i == self._n / 2 - 1:
            return 2*h[i,i]
        return h[i,i] - h[i+1,i+1]

class ExceptionalMatrixLieAlgebra(ClassicalMatrixLieAlgebra):
    """
    A matrix Lie algebra of exceptional type.
    """
    def __init__(self, R, cartan_type, e, f, h=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E',6], representation='matrix')
            sage: all(g.h(i) == g.e(i).bracket(g.f(i)) for i in range(1,7))
            True
        """
        if h is None:
            h = [e[i] * f[i] - f[i] * e[i] for i in range(len(e))]
        ClassicalMatrixLieAlgebra.__init__(self, R, cartan_type, e, f, h)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['G',2], representation='matrix')
            Simple matrix Lie algebra of type ['G', 2] over Rational Field
        """
        return "Simple matrix Lie algebra of type {} over {}".format(self.cartan_type(), self.base_ring())

class e6(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{e}_6`.

    The simple Lie algebra `\mathfrak{e}_6` of type `E_6`. The matrix
    representation is given following [HRT2000]_.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['E',6], representation='matrix')
            sage: TestSuite(g).run()  # long time
        """
        MS = MatrixSpace(R, 27, sparse=True)
        one = R.one()
        coords = [[(0,1), (10,12), (13,15), (16,17), (18,19), (20,21)],
                  [(3,4), (5,6), (7,9), (18,20), (19,21), (22,23)],
                  [(1,2), (8,10), (11,13), (14,16), (19,22), (21,23)],
                  [(2,3), (6,8), (9,11), (16,18), (17,19), (23,24)],
                  [(3,5), (4,6), (11,14), (13,16), (15,17), (24,25)],
                  [(5,7), (6,9), (8,11), (10,13), (12,15), (25,26)]]
        e = [MS({c: one for c in coord}) for coord in coords]
        f = [MS({(c[1],c[0]): one for c in coord}) for coord in coords]
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['E', 6]), e, f)

class f4(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{f}_4`.

    The simple Lie algebra `\mathfrak{f}_f` of type `F_4`. The matrix
    representation is given following [HRT2000]_ but indexed in the
    reversed order (i.e., interchange 1 with 4 and 2 with 3).
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['F',4], representation='matrix')
            sage: TestSuite(g).run()  # long time
        """
        MS = MatrixSpace(R, 26, sparse=True)
        one = R.one()

        coords = [[(0,1), (5,7), (6,9), (8,11), (10,12), (10,13), (12,14),
                   (15,16), (17,18), (19,20), (24,25)],
                  [(1,2), (3,5), (4,6), (8,10), (11,12), (11,13), (13,15),
                   (14,16), (18,21), (20,22), (23,24)],
                  [(2,3), (6,8), (9,11), (15,17), (16,18), (22,23)],
                  [(3,4), (5,6), (7,9), (17,19), (18,20), (21,22)]]
        e = [MS({c: one for c in coord}) for coord in coords]
        # Double (10, 12) in e1 and (11,13) in e2
        e[0][10,12] = 2*one
        e[1][11,13] = 2*one

        coords = [[(1,0), (7,5), (9,6), (11,8), (12,10), (14,12), (14,13),
                   (16,15), (18,17), (20,19), (25,24)],
                  [(2,1), (5,3), (6,4), (10,8), (13,11), (15,12), (15,13),
                   (16,14), (21,18), (22,20), (24,23)],
                  [(3,2), (8,6), (11,9), (17,15), (18,16), (23,22)],
                  [(4,3), (6,5), (9,7), (19,17), (20,18), (22,21)]]
        f = [MS({c: one for c in coord}) for coord in coords]
        # Double (14, 12) in f1 and (15,13) in f2
        f[0][14,12] = 2*one
        f[1][15,13] = 2*one

        # Our Cartan matrix convention is dual to that of [HRT2000]_
        e.reverse()
        f.reverse()
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['F', 4]), e, f)

class g2(ExceptionalMatrixLieAlgebra):
    r"""
    The matrix Lie algebra `\mathfrak{g}_2`.

    The simple Lie algebra `\mathfrak{g}_2` of type `G_2`. The matrix
    representation is given following [HRT2000]_.
    """
    def __init__(self, R):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = LieAlgebra(QQ, cartan_type=['G',2], representation='matrix')
            sage: TestSuite(g).run()
        """
        MS = MatrixSpace(R, 7, sparse=True)
        one = R.one()
        e = [MS({(0,1): one, (2,3): 2*one, (3,4): one, (5,6): one}),
             MS({(1,2): one, (4,5): one})]
        f = [MS({(1,0): one, (3,2): one, (4,3): 2*one, (6,5): one}),
             MS({(2,1): one, (5,4): one})]
        h = [MS({(0,0): one, (1,1): -one, (2,2): 2*one, (4,4): -2*one, (5,5): one, (6,6): -one}),
             MS({(1,1): one, (2,2): -one, (4,4): one, (5,5): -one})]
        ExceptionalMatrixLieAlgebra.__init__(self, R, CartanType(['G', 2]), e, f, h)

#######################################
## Chevalley Basis

class LieAlgebraChevalleyBasis(LieAlgebraWithStructureCoefficients):
    r"""
    A simple finite dimensional Lie algebra in the Chevalley basis.

    Let `L` be a simple (complex) Lie algebra with roots `\Phi`, then the
    Chevalley basis is given by `e_{\alpha}` for all `\alpha \in \Phi` and
    `h_{\alpha_i} := h_i` where `\alpha_i` is a simple root subject. These
    generators are subject to the relations:

    .. MATH::

        \begin{aligned}
        [h_i, h_j] & = 0
        \\ [h_i, e_{\beta}] & = A_{\alpha_i, \beta} e_{\beta}
        \\ [e_{\beta}, e_{-\beta}] & = \sum_i A_{\beta, \alpha_i} h_i
        \\ [e_{\beta}, e_{\gamma}] & = \begin{cases}
        N_{\beta,\gamma} e_{\beta + \gamma} & \beta + \gamma \in \Phi \\
        0 & \text{otherwise.} \end{cases}
        \end{aligned}

    where `A_{\alpha, \beta} = \frac{2 (\alpha, \beta)}{(\alpha, \alpha)}` and
    `N_{\alpha, \beta}` is the maximum such that
    `\alpha - N_{\alpha, \beta} \beta \in \Phi`.

    For computing the signs of the coefficients, see Section 3 of [CMT2003]_.
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Normalize ``self`` to ensure a unique represntation.

        TESTS::

            sage: L1 = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L2 = LieAlgebra(QQ, cartan_type=CartanType(['A', 2]))
            sage: L3 = LieAlgebra(QQ, cartan_type=CartanMatrix(['A', 2]))
            sage: L1 is L2 and L2 is L3
            True
        """
        if isinstance(cartan_type, (CartanMatrix, DynkinDiagram_class)):
            cartan_type = cartan_type.cartan_type()
        else:
            cartan_type = CartanType(cartan_type)
        return super(LieAlgebraChevalleyBasis, cls).__classcall__(
            cls, R, cartan_type)

    def __init__(self, R, cartan_type):
        r"""
        Initialize ``self``.

        TESTS::

            sage: L = LieAlgebra(QQ, cartan_type=['A',2])
            sage: TestSuite(L).run()  # long time
        """
        self._cartan_type = cartan_type
        RL = cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        p_roots = list(RL.positive_roots_by_height())
        n_roots = [-x for x in p_roots]
        self._p_roots_index = {al: i for i,al in enumerate(p_roots)}
        alphacheck = RL.simple_coroots()
        roots = frozenset(RL.roots())
        num_sroots = len(alpha)
        one = R.one()

        # Determine the signs for the structure coefficients from the root system
        # We first create the special roots
        sp_sign = {}
        for i,a in enumerate(p_roots):
            for b in p_roots[i+1:]:
                if a + b not in p_roots:
                    continue

                # Compute the sign for the extra special pair
                x, y = (a + b).extraspecial_pair()

                if (x, y) == (a, b): # If it already is an extra special pair
                    if (x, y) not in sp_sign:
                        # This swap is so the structure coefficients match with GAP
                        if (sum(x.coefficients()) == sum(y.coefficients())
                            and str(x) > str(y)):
                            y,x = x,y
                        sp_sign[(x, y)] = -one
                        sp_sign[(y, x)] = one
                    continue

                if b - x in roots:
                    t1 = ((b-x).norm_squared() / b.norm_squared()
                          * sp_sign[(x, b-x)] * sp_sign[(a, y-a)])
                else:
                    t1 = 0
                if a - x in roots:
                    t2 = ((a-x).norm_squared() / a.norm_squared()
                          * sp_sign[(x, a-x)] * sp_sign[(b, y-b)])
                else:
                    t2 = 0

                if t1 - t2 > 0:
                    sp_sign[(a,b)] = -one
                elif t2 - t1 > 0:
                    sp_sign[(a,b)] = one
                sp_sign[(b,a)] = -sp_sign[(a,b)]

        # Function to construct the structure coefficients (up to sign)
        def e_coeff(r, s):
            p = 1
            while r - p*s in roots:
                p += 1
            return p

        # Now we can compute all necessary structure coefficients
        s_coeffs = {}
        for i,r in enumerate(p_roots):
            # [e_r, h_i] and [h_i, f_r]
            for ac in alphacheck:
                c = r.scalar(ac)
                if c == 0:
                    continue
                s_coeffs[(r, ac)] = {r: -c}
                s_coeffs[(ac, -r)] = {-r: -c}

            # [e_r, f_r]
            s_coeffs[(r, -r)] = {alphacheck[j]: c
                                 for j, c in r.associated_coroot()}

            # [e_r, e_s] and [e_r, f_s] with r != +/-s
            # We assume s is positive, as otherwise we negate
            #   both r and s and the resulting coefficient
            for j, s in enumerate(p_roots[i+1:]):
                j += i+1 # Offset
                # Since h(s) >= h(r), we have s - r > 0 when s - r is a root
                # [f_r, e_s]
                if s - r in p_roots:
                    c = e_coeff(r, -s)
                    a, b = s-r, r
                    if self._p_roots_index[a] > self._p_roots_index[b]: # Note a != b
                        c *= -sp_sign[(b, a)]
                    else:
                        c *= sp_sign[(a, b)]
                    s_coeffs[(-r, s)] = {a: -c}
                    s_coeffs[(r, -s)] = {-a: c}

                # [e_r, e_s]
                a = r + s
                if a in p_roots:
                    # (r, s) is a special pair
                    c = e_coeff(r, s) * sp_sign[(r, s)]
                    s_coeffs[(r, s)] = {a: c}
                    s_coeffs[(-r, -s)] = {-a: -c}

        # Lastly, make sure a < b for all (a, b) in the coefficients and flip if necessary
        for k in s_coeffs.keys():
            a,b = k[0], k[1]
            if self._basis_key(a) > self._basis_key(b):
                s_coeffs[(b,a)] = [(index, -v) for index,v in s_coeffs[k].items()]
                del s_coeffs[k]
            else:
                s_coeffs[k] = s_coeffs[k].items()

        names = ['e{}'.format(i) for i in range(1, num_sroots+1)]
        names += ['f{}'.format(i) for i in range(1, num_sroots+1)]
        names += ['h{}'.format(i) for i in range(1, num_sroots+1)]
        category = LieAlgebras(R).FiniteDimensional().WithBasis()
        index_set = p_roots + n_roots + list(alphacheck)
        names = tuple(names)
        from sage.sets.finite_enumerated_set import FiniteEnumeratedSet
        index_set = FiniteEnumeratedSet(index_set)
        LieAlgebraWithStructureCoefficients.__init__(self, R, s_coeffs, names, index_set,
                                                     category, prefix='E', bracket='[',
                                                     sorting_key=self._basis_key)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, cartan_type=['A', 2])
            Lie algebra of ['A', 2] in the Chevalley basis
        """
        return "Lie algebra of {} in the Chevalley basis".format(self._cartan_type)

    def _test_structure_coeffs(self, **options):
        """
        Check the structure coefficients against the GAP implementation.

        EXAMPLES::

            sage: L = LieAlgebra(ZZ, cartan_type=['G',2])
            sage: L._test_structure_coeffs()
        """
        tester = self._tester(**options)
        ct = self.cartan_type()

        # Setup the GAP objects
        from sage.libs.gap.libgap import libgap
        L = libgap.SimpleLieAlgebra(ct.letter, ct.n, libgap(self.base_ring()))
        pos_B, neg_B, h_B = libgap.ChevalleyBasis(L)
        gap_p_roots = libgap.PositiveRoots(libgap.RootSystem(L)).sage()
        #E, F, H = libgap.CanonicalGenerators(L)

        # Setup the conversion between the Sage roots and GAP roots.
        #   The GAP roots are given in terms of the weight lattice.
        p_roots = list(ct.root_system().root_lattice().positive_roots_by_height())
        WL = ct.root_system().weight_lattice()
        La = WL.fundamental_weights()
        convert = {WL(root): root for root in p_roots}
        index = {convert[sum(c*La[j+1] for j,c in enumerate(rt))]: i
                 for i, rt in enumerate(gap_p_roots)}

        # Run the check
        basis = self.basis()
        roots = frozenset(p_roots)
        for i,x in enumerate(p_roots):
            for y in p_roots[i+1:]:
                if x + y in roots:
                    c = basis[x].bracket(basis[y]).leading_coefficient()
                    a, b = (x + y).extraspecial_pair()
                    if (x, y) == (a, b): # If it already is an extra special pair
                        tester.assertEqual(pos_B[index[x]] * pos_B[index[y]],
                                           c * pos_B[index[x+y]],
                                           "extra special pair differ for [{}, {}]".format(x, y))
                    else:
                        tester.assertEqual(pos_B[index[x]] * pos_B[index[y]],
                                           c * pos_B[index[x+y]],
                                           "incorrect structure coefficient for [{}, {}]".format(x, y))
                if x - y in roots: # This must be a negative root if it is a root
                    c = basis[x].bracket(basis[-y]).leading_coefficient()
                    tester.assertEqual(pos_B[index[x]] * neg_B[index[y]],
                                       c * neg_B[index[x-y]],
                                       "incorrect structure coefficient for [{}, {}]".format(x, y))

    def _repr_generator(self, m):
        """
        Return a string representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: K = L.basis().keys()
            sage: L._repr_generator(K[0])
            'E[alpha[2]]'
            sage: L._repr_generator(K[-1])
            'h2'
        """
        if m in self._cartan_type.root_system().root_lattice().simple_coroots():
            return "h{}".format(m.support()[0])
        return IndexedGenerators._repr_generator(self, m)

    def _latex_generator(self, m):
        r"""
        Return a latex representation of the basis element indexed by ``m``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: K = L.basis().keys()
            sage: L._latex_generator(K[0])
            'E_{\\alpha_{2}}'
            sage: L._latex_generator(K[-1])
            'h_{2}'
        """
        if m in self._cartan_type.root_system().root_lattice().simple_coroots():
            return "h_{{{}}}".format(m.support()[0])
        return IndexedGenerators._latex_generator(self, m)

    def _basis_key(self, x):
        """
        Compare two basis element indices. We order the basis elements by
        positive roots, coroots, and negative roots and then according to
        height.

        OUTPUT:

        If ``x == y``, return 0. If ``x < y``, return -1. Else return 1.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['B', 2])
            sage: K = L.basis().keys()
            sage: S = sorted(K, key=L._basis_key); S
            [alpha[2],
             alpha[1],
             alpha[1] + alpha[2],
             alpha[1] + 2*alpha[2],
             alphacheck[1],
             alphacheck[2],
             -alpha[2],
             -alpha[1],
             -alpha[1] - alpha[2],
             -alpha[1] - 2*alpha[2]]
            sage: S == K
            False
        """
        if x in self._p_roots_index:
            return self._p_roots_index[x]
        if -x in self._p_roots_index:
            return (len(self._p_roots_index)
                    + self._cartan_type.rank()
                    + self._p_roots_index[-x])
        RL = self._cartan_type.root_system().root_lattice()
        alphacheck = list(RL.simple_coroots())
        return len(self._p_roots_index) + alphacheck.index(x)

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.cartan_type()
            ['A', 2]
        """
        return self._cartan_type

    def weyl_group(self):
        """
        Return the Weyl group of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.weyl_group()
            Weyl Group of type ['A', 2] (as a matrix group acting on the ambient space)
        """
        from sage.combinat.root_system.weyl_group import WeylGroup
        return WeylGroup(self._cartan_type)

    def affine(self, kac_moody=False):
        """
        Return the affine Lie algebra of ``self``.

        EXAMPLES::

            sage: sp6 = lie_algebras.sp(QQ, 6)
            sage: sp6
            Lie algebra of ['C', 3] in the Chevalley basis
            sage: sp6.affine()
            Affine Kac-Moody algebra of ['C', 3] in the Chevalley basis
        """
        from sage.algebras.lie_algebras.affine_lie_algebra import AffineLieAlgebra
        return AffineLieAlgebra(self, kac_moody)

    # Useful in creating the UEA
    @cached_method
    def indices_to_positive_roots_map(self):
        """
        Return the map from indices to positive roots.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.indices_to_positive_roots_map()
            {1: alpha[1], 2: alpha[2], 3: alpha[1] + alpha[2]}
        """
        RL = self._cartan_type.root_system().root_lattice()
        return {i+1: r for i,r in enumerate(RL.positive_roots())}

    def lie_algebra_generators(self, str_keys=False):
        r"""
        Return the Chevalley Lie algebra generators of ``self``.

        INPUT:

        - ``str_keys`` -- (default: ``False``) set to ``True`` to have the
          indices indexed by strings instead of simple (co)roots

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 1])
            sage: L.lie_algebra_generators()
            Finite family {-alpha[1]: E[-alpha[1]], alpha[1]: E[alpha[1]], alphacheck[1]: h1}
            sage: L.lie_algebra_generators(True)
            Finite family {'f1': E[-alpha[1]], 'h1': h1, 'e1': E[alpha[1]]}
        """
        index_set = self._cartan_type.index_set()
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        B = self.basis()
        ret = {}

        if str_keys:
            for i in index_set:
                al = alpha[i]
                ret['e{}'.format(i)] = B[al]
                ret['f{}'.format(i)] = B[-al]
                ret['h{}'.format(i)] = B[alphacheck[i]]
        else:
            for i in index_set:
                al = alpha[i]
                ret[al] = B[al]
                ret[-al] = B[-al]
                ret[alphacheck[i]] = B[alphacheck[i]]

        return Family(ret)

    @cached_method
    def gens(self):
        """
        Return the generators of ``self`` in the order of `e_i`, `f_i`,
        and `h_i`.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.gens()
            (E[alpha[1]], E[alpha[2]], E[-alpha[1]], E[-alpha[2]], h1, h2)
        """
        index_set = self._cartan_type.index_set()
        RL = self._cartan_type.root_system().root_lattice()
        alpha = RL.simple_roots()
        alphacheck = RL.simple_coroots()
        B = self.basis()

        ret = []
        for i in index_set:
            ret.append(B[alpha[i]])
        for i in index_set:
            ret.append(B[-alpha[i]])
        for i in index_set:
            ret.append(B[alphacheck[i]])
        return tuple(ret)

    def highest_root_basis_elt(self, pos=True):
        r"""
        Return the basis element corresponding to the highest root `\theta`.

        INPUT:

        - ``pos`` -- (default: ``True``) if ``True``, then return
          `e_{\theta}`, otherwise return `f_{\theta}`

        EXAMPLES::

            sage: L = LieAlgebra(QQ, cartan_type=['A', 2])
            sage: L.highest_root_basis_elt()
            E[alpha[1] + alpha[2]]
            sage: L.highest_root_basis_elt(False)
            E[-alpha[1] - alpha[2]]
        """
        RL = self._cartan_type.root_system().root_lattice()
        theta = RL.highest_root()
        B = self.basis()
        if pos:
            return B[theta]
        return B[-theta]

