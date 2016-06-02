r"""
Representations of Lie Algebras

There are the following representations for Lie algebras:

- highest weight modules,
- polynomials for `\mathfrak{sl}_n`,

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
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

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.coerce_actions import GenericAction
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.categories.category_o import CategoryOInt

from sage.rings.all import ZZ, QQ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.matrix_space import MatrixSpace
from sage.modules.free_module_element import vector
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.root_system.root_system import RootSystem
from sage.algebras.lie_algebras.classical_lie_algebra import gl, ClassicalMatrixLieAlgebra

class HighestWeightRepresentation(Parent, UniqueRepresentation):
    """
    A highest weight representation of `\mathfrak{g}`.

    INPUT:

    - ``g`` -- a Lie algebra
    - ``hw`` -- a highest weight
    """
    def __init__(self, g, hw):
        """
        Initialize ``self``.
        """
        self._g = g
        self._hw = hw
        Parent.__init__( self, base=g, category=CategoryOInt(g) )

class RepresentationFromCrystal(CombinatorialFreeModule):
    """
    A representation of a Lie algebra from a crystal.

    INPUT:

    - ``crystal`` -- a crystal
    """
    def __init__(self, g, crystal, prefix='v'):
        """
        Initialize ``self``.
        """
        if g.cartan_type() != crystal.cartan_type():
            raise ValueError("mismatched Cartan type")
        self._g = g
        self._crystal = crystal
        CombinatorialFreeModule.__init__(self, base_ring, category=CategoryOInt(self._g),
                                         prefix=prefix, bracket="")

    class Element(CombinatorialFreeModuleElement):
        def e(self, i):
            """
            Return the action of `e_i` on ``self``.
            """
            return self.parent().sum_of_terms([(m.e(i), c) for m,c in self])

        def f(self, i):
            """
            Return the action of `e_i` on ``self``.
            """
            return self.parent().sum_of_terms([(m.f(i), c) for m,c in self])

        def h(self, i):
            """
            Return the action of `h_i` on ``self``.
            """
            h = self._crystal.cartan_type().weight_lattice_realization().simple_coroots()
            return self.parent().sum_of_terms([(m, h[i].scalar(m.weight())*c) for m,c in self])

class LaurentPolynomialRepresentation(Parent, UniqueRepresentation):
    r"""
    The representation of `\mathfrak{sl}_n` or `\mathfrak{gl}_n`
    using Laurent polynomials in `\{ x_1, \ldots, x_n \}`. This is also
    known as the Weyl representation.

    Consider `(a_{ij})_{i,j} = A \in \mathfrak{gl}_n`, then the Weyl
    representation `\mathcal{D}` is defined by:

    .. MATH::

        \mathcal{D}(A) = \sum_{i,j} a_{ij} x_i \frac{\partial}{\partial x_j}.

    INPUT:

    - ``g`` -- the Lie algebra `\mathfrak{sl}_n` or `\mathfrak{gl}_n`
    - ``la`` -- the highest weight
    - ``var_name`` -- (default: ``'x'``) the variable name
    """
    def __init__(self, g, la, var_name='z'):
        """
        Initialize ``self``.
        """
        if hasattr(g, "cartan_type"):
            if g.cartan_type().type() != 'A' or not g.cartan_type().is_finite():
                raise ValueError("only for finite type A")
            self._n = g.cartan_type().rank() + 1 # A_n => sl_{n+1}
        elif isinstance(g, gl):
            self._n = g._n
        else:
            raise ValueError("invalid Lie algebra; can only be sl_n or gl_n")
        self._amb = RootSystem(['A', self._n-1]).ambient_space()
        self._la = la
        self._g = g
        self._lp_ring = LaurentPolynomialRing(g.base_ring(), self._n, var_name)

        Parent.__init__( self, base=g, category=CategoryOInt(g) )

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "The Laurent polynomial representation of {}".format(self._g)

    def highest_weight_vector(self):
        """
        Return the higehst weight vector of ``self``.
        """
        return self(prod( self._lp_ring.gen(i)**c for i,c in self._amb(self._la) ))

    class Element(ElementWrapper):
        def __eq__(self, rhs):
            if rhs == 0 and not self.value:
                return True
            if isinstance(rhs, LaurentPolynomialRepresentation.Element):
                return rhs.parent() is self.parent() and rhs.value == self.value
            return False

        def normalize(self):
            return self.__class__(self.parent(), self.value / self.value.coefficients()[0])

        def derivative(self, i):
            """
            Take the derivative of ``self`` with respect to the ``i``-th
            variable.
            """
            v = self.parent()._lp_ring.gen(i)
            return self.parent()( sum(c*m.exponents()[0][i] * m * v**(-1) for c,m in self.value) )

        def e(self, i):
            """
            Apply the action of `e_i` on ``self``.
            """
            P = self.parent()
            return P( P._lp_ring.gen(i-1) * self.derivative(i).value )

        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.
            """
            P = self.parent()
            return P( P._lp_ring.gen(i) * self.derivative(i-1).value )

        def h(self, i):
            """
            Apply the action of `h_i` on ``self``.
            """
            P = self.parent()
            return P( P._lp_ring.gen(i-1) * self.derivative(i-1).value
                      - P._lp_ring.gen(i) * self.derivative(i).value )

        @cached_method
        def weight(self):
            r"""
            Return the weight of ``self``.
            """
            La = self.parent()._g.cartan_type().root_system().weight_space().fundamental_weights()
            return sum(La[i] * QQ(self.h(i).value / self.value) for i in range(1, self.parent()._n))

        def _l_action_(self, x):
            """
            Return the left action of `\mathfrak{sl}_n` on ``self``.
            """
            if x in self._g:
                x = x.lift()
            uea = self._g.universal_enveloping_algebra()
            if x not in uea:
                raise ValueError('cannot convert {} to an action'.format(x))

            # If it is a matrix
            if x in MatrixSpace(self.base_ring(), self._n):
                P = self.parent()
                return P(sum( x[i,j] * P._lp_ring.gen(j) * self.derivative(i).value
                                 for i in range(self._n) for j in range(self._n) ))

            s = P.zero()
            for cg, monomial in g:
                c = self
                for x in monomial.variables():
                    if x in uea._e:
                        c = c.e(uea._e.index(x))
                    if x in uea._f:
                        c = c.f(uea._f.index(x))
                    else:
                        c = c.h(uea._h.index(x))
                s += cg * c
            return s

class PolynomialRepresentation(Parent, UniqueRepresentation):
    """
    The polynomial representation of `\mathfrak{sl}_n` using the array of
    variables `z_{ij}` for `1 \leq i < j \leq n`.

    INPUT:

    - ``g`` -- the Lie algebra `\mathfrak{sl}_n`
    - ``la`` -- the highest weight
    - ``var_name`` -- (default: ``'z'``) the variable name
    """
    def __init__(self, g, la, var_name='z'):
        """
        Initialize ``self``.
        """
        ct = g.cartan_type()
        if ct.type() != 'A' or not ct.is_finite():
            raise ValueError("Only for finite type A")
        self._la = la
        self._g = g
        n = ct.rank() + 1
        names = [var_name + repr(i) + "c" + repr(j) for i in range(1, n+1) for j in range(i+1, n+1)]
        self._p_ring = PolynomialRing(g.base_ring(), names)
        self._vars = []
        v = self._p_ring.gens()
        cur = 0
        for i in range(1, n):
            self._vars.append([None]*i + [v[cur+j] for j in range(n-i)])
            cur += n - i
        self._n = n
        Parent.__init__( self, base=g, category=CategoryOInt(g) )

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "The polynomial representation of {}".format(self._g)

    def highest_weight_vector(self):
        """
        Return the higehst weight vector of ``self``.
        """
        return self(self._p_ring.one())

    class Element(ElementWrapper):
        def __eq__(self, rhs):
            if rhs == 0 and not self.value:
                return True
            if isinstance(rhs, PolynomialRepresentation.Element):
                return rhs.parent() is self.parent() and rhs.value == self.value
            return False

        def normalize(self):
            return self.__class__(self.parent(), self.value / self.value.coefficients()[0])

        def e(self, i):
            """
            Apply the action of `e_i` on ``self``.
            """
            P = self.parent()
            v = P._vars
            i -= 1 # For indexing
            return P( self.value.derivative(v[i][i+1])
                      + sum(v[j][i] * self.value.derivative(v[j][i+1])
                            for j in range(i)) )

        def f(self, i):
            """
            Apply the action of `f_i` on ``self``.
            """
            P = self.parent()
            n = P._n
            v = P._vars
            i -= 1 # For indexing
            return P( P._la[i+1] * v[i][i+1] * self.value
                      + sum(v[j][i+1] * self.value.derivative(v[j][i])
                            for j in range(i))
                      - sum(v[i][j] * self.value.derivative(v[i+1][j])
                            for j in range(i+2, n))
                      - v[i][i+1] * ( sum(v[i][j] * self.value.derivative(v[i][j])
                                          for j in range(i+1, n))
                                     - sum(v[i+1][j] * self.value.derivative(v[i+1][j])
                                           for j in range(i+2, n)) ) )

        def h(self, i):
            """
            Apply the action of `h_i` on ``self``.
            """
            P = self.parent()
            n = P._n
            v = P._vars
            i -= 1 # For indexing
            return P( P._la[i+1] * self.value
                      - 2 * v[i][i+1] * self.value.derivative(v[i][i+1])
                      + sum(v[j][i] * self.value.derivative(v[j][i])
                            - v[j][i+1] * self.value.derivative(v[j][i+1])
                            for j in range(i))
                      + sum(v[i+1][j] * self.value.derivative(v[i+1][j])
                            - v[i][j] * self.value.derivative(v[i][j])
                            for j in range(i+2, n)) )

        @cached_method
        def weight(self):
            r"""
            Return the weight of ``self`` as a list.

            The weight of a vector `v` in a (highest weight) representation
            `V` is defined as `\lambda in \mathfrak{h}^*` such that
            `hv = \lambda(h) v` for all `h \in \mathfrak{h}`. Here we express
            the weight as `\sum_{i=1}^n \lambda(h_i) \Lambda_i` where
            `\{ h_i \}` are the simple coroots and `\{ \Lambda_i \}` are the
            fundamental weights.
            """
            La = self.parent()._g.cartan_type().root_system().weight_space().fundamental_weights()
            return sum(La[i] * QQ(self.h(i).value / self.value) for i in range(1, self.parent()._n))

        @cached_method
        def root(self):
            """
            Return the root of ``self``.
            """
            cm = self.parent()._g.cartan_type().cartan_matrix()
            Q = cm.root_system().root_space()
            alpha = Q.simple_roots()
            # +1 for indexing
            vec = self.weight().to_vector()
            return sum( alpha[i+1] * x for i,x in enumerate(cm.solve_right(vec)) )

        def _l_action_(self, x):
            """
            Return the left action of `\mathfrak{sl}_n` on ``self``.
            """
            if x in self._g:
                x = x.lift()
            uea = self._g.universal_enveloping_algebra()
            if x not in uea:
                raise ValueError('cannot convert {} to an action'.format(x))

            s = self.zero()
            for cg, monomial in g:
                c = self
                for x in monomial.variables():
                    if x in uea._e:
                        c = c.e(uea._e.index(x))
                    if x in uea._f:
                        c = c.f(uea._f.index(x))
                    else:
                        c = c.h(uea._h.index(x))
                s += cg * c
            return s

