"""
Lie Algebras Given By Structure Coefficients

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

from copy import copy
from sage.misc.cachefunc import cached_method
from sage.misc.indexed_generators import IndexedGenerators
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper

from sage.categories.algebras import Algebras
from sage.categories.lie_algebras import LieAlgebras

from sage.algebras.free_algebra import FreeAlgebra
from sage.algebras.lie_algebras.lie_algebra_element import LieGenerator, \
    LieBracket, LieAlgebraElement
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.subalgebra import LieSubalgebra
from sage.algebras.lie_algebras.ideal import LieAlgebraIdeal
from sage.algebras.lie_algebras.quotient import QuotientLieAlgebra
from sage.rings.all import ZZ
from sage.rings.ring import Ring
from sage.rings.integer import Integer
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.infinity import infinity
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.modules.free_module import FreeModule, span
from sage.sets.family import Family, AbstractFamily

# Move to FiniteDimLieAlgWithBasis category
class LieAlgebraWithStructureCoefficients(FinitelyGeneratedLieAlgebra, IndexedGenerators):
    r"""
    A Lie algebra with a set of specified structure coefficients.

    The structure coefficients are specified as a dictionary whose keys are
    pairs of generators and values are dictionaries of generators mapped
    to coefficients.

    EXAMPLES:

    We create the Lie algebra of `\QQ^3` under the Lie bracket defined
    by `\times` (cross-product)::

        sage: L = LieAlgebra(QQ, 'x,y,z', {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
        sage: (x,y,z) = L.gens()
        sage: L.bracket(x, y)
        z
        sage: L.bracket(y, x)
        -z
    """
    @staticmethod
    def __classcall_private__(cls, R, s_coeff, names=None, index_set=None):
        """
        Normalize input to ensure a unique representation.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L2 = LieAlgebra(QQ, 'x,y', {('y','x'):{'x':-1}})
            sage: L is L2
            True
        """
        s_coeff = LieAlgebraWithStructureCoefficients._standardize_s_coeff(s_coeff)
        if len(s_coeff) == 0:
            return AbelianLieAlgebra(R, names, index_set)

        if names is None:
            if index_set is None:
                raise ValueError("either the names or the index set must be specified")
            if len(index_set) <= 1:
                return AbelianLieAlgebra(R, names, index_set)
        elif len(names) <= 1:
            return AbelianLieAlgebra(R, names, index_set)

        return super(LieAlgebraWithStructureCoefficients, cls).__classcall__(
            cls, R, s_coeff, tuple(names), index_set)

    @staticmethod
    def _standardize_s_coeff(s_coeff):
        """
        Helper function to standardize ``s_coeff`` into the appropriate tuple
        of tuples. Strips items with coefficients of 0 and duplicate entries.
        This does not check the Jacobi relation (nor antisymmetry if the
        cardinality is infinite).
        """
        if isinstance(sc, AbstractFamily) and sc.cardinality() == infinity:
            return sc
        sc = {}
        if isinstance(s_coeff, dict):
            s_coeff = s_coeff.iteritems()
        # Make sure the first gen is smaller than the second in each key
        for k,v in s_coeff:
            if isinstance(v, dict):
                v = v.items()
            if k[0] > k[1]:
                vals = tuple((g, -val) for g, val in v if val != 0)
            else:
                key = LieBracket(*k)
                vals = tuple((g, val) for g, val in v if val != 0)

            if key in sc.keys() and sorted(sc[key]) != sorted(vals):
                raise ValueError("non-equal brackets")

            if len(vals) > 0:
                sc[key] = vals
        return Family(sc)

    def __init__(self, R, s_coeff, names, index_set):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: TestSuite(L).run()
        """
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        FinitelyGeneratedLieAlgebra.__init__(self, R, names, index_set, cat)
        self.__s_coeff = s_coeff

    def basis(self):
        """
        Return the basis of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.basis()
            (x, y)
        """
        return Family({i: self.monomial(i) for i in self._indices})

    def get_order(self):
        """
        Return the order of the elements in the basis.
        """
        return sorted(map(lambda x: x.support()[0], self.basis()))

    def structure_coefficients(self):
        """
        Return the dictonary of structure coefficients of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.structure_coefficients()
            Finite family {[x, y]: ((x, 1),)}
        """
        return self.__s_coeff

    def dimension(self):
        """
        Return the dimension of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y', {('x','y'):{'x':1}})
            sage: L.dimension()
            2
        """
        return self.basis().cardinality()

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket of ``[self, y]``.

        EXAMPLES::

            sage: L.<x,y,z> = LieAlgebra(QQ, {('x','y'):{'z':1}, ('y','z'):{'x':1}, ('z','x'):{'y':1}})
            sage: L.bracket(x, y) # indirect doctest
            z
            sage: L.bracket(y, x)
            -z
            sage: L.bracket(x + y - z, x - y + z)
            -2*y - 2*z
        """
        b = LieBracket(x, y)
        if b not in self.__s_coeff:
            return self.zero()
        return self.element_class(self, dict(self.__s_coeff[b]))

    def free_module(self, sparse=True):
        """
        Return the ``self`` as a free module.
        """
        return FreeModule(self.base_ring(), self.dimension(), sparse=sparse)

    def subalgebra(self, gens, names=None):
        """
        Return the subalgebra of ``self`` generated by ``gens``.
        """
        return LieSubalgebraWithStructureCoefficients(self, gens, names)

    class Element(LieAlgebraElement):
        """
        An element of a Lie algebra given by structure coefficients.
        """
        def to_vector(self):
            """
            Return ``self`` as a vector.
            """
            V = self.parent().free_module()
            B = V.basis()
            return V.sum(B[k]*c for k,c in self)

class LieSubalgebraWithStructureCoefficients(LieSubalgebra):
    """
    A Lie subalgebra of a Lie algebra given by structure coefficients.
    """
    def subalgebra(self, gens, names=None):
        """
        Return the subalgebra of ``self`` generated by ``gens``.
        """
        return LieSubalgebraWithStructureCoefficients(self._ambient, map(lambda x: x.value, gens), names)

class LieAlgebraIdealWithStructureCoefficients(LieAlgebraIdeal,
                        LieSubalgebraWithStructureCoefficients):
    """
    A Lie algebra ideal of a Lie algebra given by structure coefficients.
    """
    def reduce(self, y):
        """
        Return ``y`` modulo ``self``.
        """
        M = self._ambient.free_module()
        I = M.subspace(map(lambda x: x.to_vector(), self._gens))
        R = I.complement()
        ld = M.linear_dependence(R.basis() + I.basis() + [y.to_vector()])[0]
        normalizer = -ld[-1]
        coeffs = map(lambda x: x / normalizer, ld[:len(R)])
        vec = M(R.linear_combination_of_basis(coeffs))
        return self._ambient(vec)

# This should either not inherit from QuotientLieAlgebra or refactor out common code
class QuotientLieAlgebraWithStructureCoefficients(QuotientLieAlgebra):
    """
    A quotient Lie algebra of a Lie algebra given by structure coefficients.
    """
    def __init__(self, lie, I, names=None, index_set=None, category=None):
        """
        Initialize ``self``.
        """
        R = lie.base_ring()
        cat = LieAlgebras(R).FiniteDimensional().WithBasis()
        QuotientLieAlgebra.__init__(self, lie, I, names, index_set, cat)

        # Construct the structure coefficients
        M = lie.free_module()
        IM = M.subspace(map(lambda x: x.to_vector(), I._gens))
        RM = IM.complement()
        B = RM.basis() + IM.basis()
        dim = len(RM)
        self.__s_coeff = {}
        for i,kl,zl in enumerate(gens):
            for kr,zr in gens[i+1:]:
                b = LieBracket(zl, zr)
                ld = M.linear_dependence(B + [y.to_vector()])[0]
                normalizer = -ld[-1]
                self.__s_coeffs[b] = {indices[j]: val / normalizer
                                      for j,val in enumerate(ld[:dim])}

    Element = LieAlgebraElement

class AbelianLieAlgebra(LieAlgebraWithStructureCoefficients):
    r"""
    An abelian Lie algebra.

    A Lie algebra `\mathfrak{g}` is abelian if `[x, y] = 0` for all
    `x, y \in \mathfrak{g}`.

    EXAMPLES::

        sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
        sage: L.bracket(x, y)
        0
    """
    def __init__(self, R, names=None, index_set=None):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: TestSuite(L).run()
        """
        LieAlgebraWithStructureCoefficients.__init__(self, R, Family({}), names, index_set)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LieAlgebra(QQ, 3, 'x', abelian=True)
            Abelian Lie algebra on 3 generators (x0, x1, x2) over Rational Field
        """
        gens = self.lie_algebra_generators()
        return "Abelian Lie algebra on {} generators {} over {}".format(
            gens.cardinality(), gens.values(), self.base_ring())

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L._construct_UEA()
            Multivariate Polynomial Ring in x0, x1, x2 over Rational Field
        """
        return PolynomialRing(self.base_ring(), self.variable_names())

    def is_abelian(self):
        """
        Return ``True`` since this is an abelian Lie algebra.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 3, 'x', abelian=True)
            sage: L.is_abelian()
            True
        """
        return True

    def bracket_on_basis(self, x, y):
        """
        Return the Lie bracket of basis elements indexed by ``x`` and ``y``.

        EXAMPLES::

            sage: L.<x,y> = LieAlgebra(QQ, abelian=True)
            sage: L.bracket_on_basis(x.leading_support(), y.leading_support())
            0
        """
        return self.zero()

    class Element(LieAlgebraElement):
        def _bracket_(self, y):
            """
            Return the Lie bracket ``[self, y]``.

            EXAMPLES::

                sage: L.<x, y> = LieAlgebra(QQ, abelian=True)
                sage: L.bracket(x, y)
                0
            """
            return self.parent().zero()

