"""
Classical Lie Algebras

These are the Lie algebras corresponding to types `A_n`, `B_n`, `C_n`, `D_n`.

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version

EXAMPLES::
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.algebras.algebra import Algebra
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra, LieAlgebraFromAssociative
from sage.algebras.lie_algebras.lie_algebra_element import LieAlgebraElement
from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.lie_algebras import LieAlgebras
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.family import Family

class ClassicalMatrixLieAlgebra(LieAlgebraFromAssociative):
    """
    A classical Lie algebra represented using matrices.
    """
    @staticmethod
    def __classcall_private__(cls, R, cartan_type):
        """
        Return the correct parent based on input.

        INPUT:

        - ``R`` -- the base ring
        - ``ct`` -- the Cartan type

        EXAMPLES::

            sage: from sage.algebras.lie_algebras.classical_lie_algebra import ClassicalMatrixLieAlgebra
            sage: ClassicalMatrixLieAlgebra(QQ, ['A', 4])
            Special linear Lie algebra of rank 5 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, CartanType(['B',4]))
            Special orthogonal Lie algebra of rank 9 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, 'C4')
            Symplectic Lie algebra of rank 8 over Rational Field
            sage: ClassicalMatrixLieAlgebra(QQ, cartan_type=['D',4])
            Special orthogonal Lie algebra of rank 8 over Rational Field
        """
        cartan_type = CartanType(cartan_type)

        if cartan_type.is_affine():
            if not cartan_type.is_untwisted_affine():
                raise NotImplementedError("Only currently implemented for untwisted affine types")
            classical = ct.classical()
            return AffineLieAlgebra(ClassicalMatrixLieAlgebra(R, classical))

        if cartan_type.type() == 'A':
            return sl(R, cartan_type.rank() + 1)
        if cartan_type.type() == 'B':
            return so(R, 2*cartan_type.rank() + 1)
        if cartan_type.type() == 'C':
            return sp(R, 2*cartan_type.rank())
        if cartan_type.type() == 'D':
            return so(R, 2*cartan_type.rank())
        raise NotImplementedError

    def __init__(self, R, ct, e, f, h):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
            sage: TestSuite(g).run()
        """
        n = len(e)
        names = ['e%s'%i for i in range(1, n+1)]
        names += ['f%s'%i for i in range(1, n+1)]
        names += ['h%s'%i for i in range(1, n+1)]
        LieAlgebraFromAssociative.__init__(self, R, e[0].parent(),
                                           tuple(e + f + h), tuple(names),
                                           category=LieAlgebras(R).FiniteDimensional().WithBasis())
        self._cartan_type = ct

        gens = self.gens()
        i_set = ct.index_set()
        self._e = Family(dict( (i, gens[c]) for c,i in enumerate(i_set) ))
        self._f = Family(dict( (i, gens[n+c]) for c,i in enumerate(i_set) ))
        self._h = Family(dict( (i, gens[2*n+c]) for c,i in enumerate(i_set) ))

    def e(self, i):
        r"""
        Return the generator `e_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
            sage: g.e(2)
            [0 0 0]
            [0 0 1]
            [0 0 0]
        """
        return self._e[i]

    def f(self, i):
        r"""
        Return the generator `f_i`.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
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

            sage: g = lie_algebras.sl(QQ, 3)
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

            sage: g = lie_algebras.sl(QQ, 3)
            sage: g.index_set()
            (1, 2)
        """
        return self._cartan_type.index_set()

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
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

            sage: g = lie_algebras.sl(QQ, 3)
            sage: g.epsilon(1, g.h(1))
            1
            sage: g.epsilon(2, g.h(1))
            -1
            sage: g.epsilon(3, g.h(1))
            0
        """
        return h[i-1,i-1]

    @abstract_method
    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``h``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
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

            sage: g = lie_algebras.sl(QQ, 3)
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
        for j in w:
            for k in range(-r.scalar(coroots[j])):
                cur = self.bracket(gens[j], cur)
        return cur

    @cached_method
    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
            sage: tuple(g.basis())
            (
            [0 0 0]  [0 0 0]  [ 0  0  0]  [ 1  0  0]  [0 1 0]  [0 0 0]
            [1 0 0]  [0 0 0]  [ 0  1  0]  [ 0 -1  0]  [0 0 0]  [0 0 1]
            [0 0 0], [0 1 0], [ 0  0 -1], [ 0  0  0], [0 0 0], [0 0 0]
            )
        """
        d = {}
        for i in self.index_set():
            d['e{}'.format(i)] = self._e[i]
            d['f{}'.format(i)] = self._f[i]
            d['h{}'.format(i)] = self._h[i]
        return Family(d)

    def affine(self, kac_moody=False):
        """
        Return the (untwisted) affine Lie algebra of ``self``.
        """
        return AffineLieAlgebra(self, kac_moody)

class gl(LieAlgebraFromAssociative):
    r"""
    The Lie algebra `\mathfrak{gl}_n` which consists of all `n \times n`
    matrices.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 4)
            sage: TestSuite(g).run()
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
        LieAlgebraFromAssociative.__init__(self, R, MS, tuple(gens), tuple(names),
                                           category=LieAlgebras(R).FiniteDimensional().WithBasis())

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
        """
        return 2 * self._n * (x.value * y.value).trace() \
            - 2 * x.value.trace() * y.value.trace()

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.gl(QQ, 2)
            sage: g.basis()
            (
            [1 0]  [0 1]  [0 0]  [0 0]
            [0 0], [0 0], [1 0], [0 1]
            )
        """
        return self.gens()

class sl(ClassicalMatrixLieAlgebra):
    r"""
    The Lie algebra `\mathfrak{sl}_n` which consists of `n \times n` matrices
    with trace 0. This is the Lie algebra of type `A_{n-1}`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5)
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

            sage: lie_algebras.sl(QQ, 5)
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

            sage: g = lie_algebras.sl(QQ, 5)
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 5)
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
    The Lie algebra `\mathfrak{so}_n` which consists of orthogonal `n \times n`
    matrices. This is the Lie algebra of type `B_{(n-1)/2}` or `D_{n/2}` if `n`
    is odd or even respectively.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.so(QQ, 8)
            sage: TestSuite(g).run()
            sage: g = lie_algebras.so(QQ, 9)
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

            sage: g = lie_algebras.so(QQ, 8)
            sage: g = lie_algebras.so(QQ, 9)
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES:

        The even or type `D` case::

            sage: g = lie_algebras.so(QQ, 8)
            sage: matrix([[g.simple_root(i, g.h(j)) for i in g.index_set()] for j in g.index_set()])
            [ 2 -1  0  0]
            [-1  2 -1 -1]
            [ 0 -1  2  0]
            [ 0 -1  0  2]

        The odd or type `B` case::

            sage: g = lie_algebras.so(QQ, 9)
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
    The Lie algebra `\mathfrak{sp}_{2n}` which consists of `2n \times 2n`
    matrices `X` that satisfy the equation:

    .. MATH::

        X^T M - M X = 0

    where

    .. MATH::

        M = \begin{pmatrix}
        0 & I_n \\
        -I_n & 0
        \end{pmatrix}.

    This is the Lie algebra of type `C_n`.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8)
            sage: TestSuite(g).run()
        """
        if n % 2 != 0:
            raise ValueError("the rank must be even")
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

            sage: lie_algebras.sp(QQ, 8)
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

            sage: g = lie_algebras.sp(QQ, 8)
        """
        return (2 * self._n + 2) * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sp(QQ, 8)
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

class su(ClassicalMatrixLieAlgebra):
    r"""
    The Lie algebra `\mathfrak{su}_n` which consists of unitary `n \times n`
    matrices.
    """
    def __init__(self, R, n):
        """
        Initialize ``self``.
        """
        raise NotImplementedError
        self._n = n
        one = R.one()
        #gens = {('E[%s,%s]'%i,j): MS({(i,j):one}) for i in range(n) for j in range(n)}
        #e = [MS({(i,i+1):one}) for i in range(n-1)]
        #f = [MS({(i+1,i):one}) for i in range(n-1)]
        h = [MS({(i,i):one, (i+1,i+1):-one}) for i in range(n-1)]
        gens = {}
        LieAlgebraFromAssociative.__init__(self, R, gens)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: g = lie_algebras.su(QQ, 4) # not tested -- todo
            Special orthogonal Lie algebra of rank 9 over Rational Field
        """
        return "Special unitary Lie algebra of rank {} over {}".format(self._n, self.base_ring())

    def killing_form(self, x, y):
        r"""
        Return the Killing form on ``x`` and ``y``.

        The Killing form on `\mathfrak{su}_n` is:

        .. MATH::

            \langle x \mid y \rangle = 2n \mathrm{tr}(xy).
        """
        return 2 * self._n * (x.value * y.value).trace()

    def simple_root(self, i, h):
        r"""
        Return the action of the simple root
        `\alpha_i \colon \mathfrak{h} \to R`, where `R` is the base
        ring of ``self``, on the element ``j``.

        EXAMPLES::

            sage: g = lie_algebras.sl(QQ, 3)
            sage: g.simple_root(1, g.h(1))
            2
            sage: g.simple_root(1, g.h(2))
            -1
        """
        i = self.index_set().index(i)
        return h[i,i] - h[i+1,i+1]

class AffineLieAlgebra(FinitelyGeneratedLieAlgebra):
    r"""
    Affine Lie algebra.

    Given a finite dimensional simple Lie algebra `\mathfrak{g}` over `R`,
    we construct an affine Lie algebra `\hat{\mathfrak{g}}` as

    .. MATH::

        \hat{\mathfrak{g}} = \left( \mathfrak{g} \otimes R[t, t^{-1}] \right)
        \oplus R c

    where `c` is the canonical central element and `R[t, t^{-1}]` is the
    Laurent polynomial ring over `R`. We define the Lie bracket as

    .. MATH::

        [a \otimes t^n + \alpha c, b \otimes t^m + \beta c] =
        [a, b] \otimes t^{n+m} + \delta_{n+m,0} \langle a \mid b \rangle n c

    where `\langle a \mid b \rangle` is the Killing form on `\mathfrak{g}`.

    We can also form the affine Kac--Moody algebra by adding the additional
    generator `d` such that `[d, x] = \delta(x)` where `\delta` is the
    Lie derivative.
    """
    def __init__(self, g, kac_moody=False):
        """
        Initalize ``self``.
        """
        self._g = g
        self._cartan_type = g.cartan_type().affine()
        R = g.base_ring()
        names = list(g.variable_names()) + ['e0', 'f0', 'c']
        if kac_moody:
            names += ['d']
        self._kac_moody = kac_moody
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, LieAlgebras(R).WithBasis())

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        base = "Affine "
        rep = repr(self._g)
        if self._kac_moody:
            old_len = len(rep)
            rep = rep.replace("Lie", "Kac-Moody")
            if len(rep) == old_len: # We did not replace anything
                base += "Kac-Moody "
        return base + rep

    def derived_subalgebra(self):
        """
        Return the derived subalgebra of ``self``.
        """
        if self._kac_moody:
            return AffineLieAlgebra(self._g)
        raise NotImplementedError

    def cartan_type(self):
        """
        Return the Cartan type of ``self``.
        """
        return self._cartan_type

    def classical(self):
        """
        Return the classical Lie algebra of ``self``.
        """
        return self._g

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.
        """
        # These are placeholders and will require something more complex
        if self._kac_moody:
            PolynomialRing(self._g.universal_enveloping_algebra(), 't,c,d')
        return PolynomialRing(self._g.universal_enveloping_algebra(), 't,c')

    def gen(self, i):
        """
        Return the `i`-th generator of ``self``.
        """
        n = self.ngens()
        if self._kac_moody:
            if i == n - 1:
                return self.element_class(self, {'d': 1})
            n -= 1
        # If it is a standard generator
        if i < n - 3:
            return self.element_class(self, {(self._g.gen(i), 0): 1})

        if i == n - 3: # e_0 = f_{\theta} t
            return self.element_class(self, {(self._g.highest_root_basis_elt(False), 1):1})
        if i == n - 2: # f_0 = e_{\theta} t^-1
            return self.element_class(self, {(self._g.highest_root_basis_elt(True), -1):1})

        if i == n - 1: # c
            return self.element_class(self, {'c': 1})
        raise ValueError("i is too large")

    # A summand is either pairs (a, i) where a is an element in the classical
    #   Lie algebra and i is in ZZ representing the power of t, or it's 'c'
    #   for the canonical central element.
    class Element(LieAlgebraElement):
        """
        Element of an affine Lie algebra.
        """
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.
            """
            s_mon = sorted(self.list())
            y_mon = sorted(y.list())
            if not self or not y or s_mon == y_mon:
                return self.parent().zero()
            d = {}
            gd = self.parent()._g.gens_dict()
            for ml,cl in sorted(s_mon): # The left monomials
                for mr,cr in sorted(y_mon): # The right monomials
                    if ml == mr or ml == 'c' or mr == 'c':
                        continue
                    if ml == 'd' and mr != 'd' and mr[1] != 0:
                        d[mr] = cr * mr[1]
                        continue
                    if mr == 'd' and ml[1] != 0:
                        d[ml] = -cl * ml[1]
                        continue
                    gl,tl = ml
                    gr,tr = mr
                    b = gl._bracket_(gr)
                    if b:
                        for m, c in b:
                            d[(gd[m], tl+tr)] = cl * cr * c
                    if tl != 0 and tr + tl == 0:
                        d['c'] = gl.killing_form(gr) * cl * cr * tl
            if len(d) == 0:
                return self.parent().zero()
            return self.__class__(self.parent(), d)

        def lie_derivative(self):
            r"""
            Return the Lie derivative `\delta` of ``self``.

            The Lie derivative `\delta` is defined as

            .. MATH::

                \delta(a \otimes t^m + \alpha c) = a \otimes m t^m.

            Another formulation is by `\delta = t \frac{d}{dt}`.
            """
            d = {}
            for m, c in self.list():
                if m != 'c' and m[1] != 0:
                    d[m] = c*m[1]
            return self.__class__(self.parent(), d)

