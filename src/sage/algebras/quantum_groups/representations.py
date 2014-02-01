r"""
Representations of Quantum Groups

There are the following representations for quantum groups:

- highest weight modules,
- representations from crystals,
- matrices for finite types.

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
from sage.structure.element_wrapper import ElementWrapper
from sage.misc.cachefunc import cached_method
from sage.categories.algebras_with_basis import AlgebrasWithBasis
from sage.categories.category_o import CategoryOInt
from sage.categories.highest_weight_crystals import HighestWeightCrystals

from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
from sage.algebras.lie_algebras.classical_lie_algebra import ClassicalMatrixLieAlgebra
from sage.algebras.quantum_groups.q_numbers import q_factorial

class HighestWeightRepresentation(Parent, UniqueRepresentation):
    """
    A highest weight representation of `U_q(\mathfrak{g})`.

    INPUT:

    - ``g`` -- a Lie algebra or quantum group
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
    A representation of a quantum group from a crystal given in the (local)
    crystal basis.

    INPUT:

    - ``g`` -- a quantum group
    - ``crystal`` -- a crystal
    """
    @staticmethod
    def __classcall_private__(cls, g, crystal, prefix='v'):
        """
        Normalize input to ensure a unique representation.
        """
        if isinstance(crystal, TensorProductOfCrystals):
            return RepresentationFromTensorProductOfCrystals(g, crystal, prefix)
        return super(RepresentationFromCrystal, cls).__classcall__(cls, g, crystal, prefix)

    def __init__(self, g, crystal, prefix='v'):
        """
        Initialize ``self``.
        """
        if g.cartan_type() != crystal.cartan_type():
            raise ValueError("mismatched Cartan type")
        #if crystal not in HighestWeightCrystals():
        #    raise ValueError("only for highest weight crystals")
        self._g = g
        # This is the wrong category, should be CategoryO_q
        #   and needs checks for Int (ex. if a normal crystal, it is in Int)
        CombinatorialFreeModule.__init__(self, g.base_ring(), crystal, category=CategoryOInt(self._g),
                                         prefix=prefix, bracket=False)

    @cached_method
    def highest_weight_vector(self, b, index_set=None):
        """
        Return the highest weight vector corresponding to the highest weight
        crystal element ``b``.
        """
        if index_set is None:
            if not b.is_highest_weight():
                raise ValueError("{} is not highest weight".format(b))
            index_set = self._g.cartan_type().index_set()
        elif not b.is_highest_weight(index_set):
            raise ValueError("{} is not {} highest weight".format(b, index_set))
        cur = self.monomial(b)
        for i in index_set:
            next = cur.e(i)
            while next != 0:
                m,c = next.leading_item()
                elt = self.monomial(m).f(i)
                temp = elt - elt[b] * cur
                cur -= c / temp.e(i).leading_coefficient() * temp
                next = cur.e(i)
        return cur

    def e_on_basis(self, i, b):
        """
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        q = self._g.q(i)
        h = self.basis().keys().weight_lattice_realization().simple_coroots()
        hw_elt, path = b.to_highest_weight()
        ret = []
        for j in path:
            b = b.e(j)
            if i == j:
                p = h[i].scalar(b.weight())
                ret.append( (b, (q**p - q**-p)/(q - q**-1)) )
        return self.sum_of_terms(ret)

    def f_on_basis(self, i, b):
        """
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        next = b.f(i)
        if next is None:
            return self.zero()
        return self.monomial(next)

    class Element(CombinatorialFreeModuleElement):
        def h(self, i):
            """
            Return the action of `h_i` on ``self``.
            """
            P = self.parent()
            q = P._g.q(i)
            h = P.basis().keys().weight_lattice_realization().simple_coroots()
            return P.sum_of_terms([(m, c * q**h[i].scalar(m.weight())) for m,c in self])

class RepresentationFromTensorProductOfCrystals(RepresentationFromCrystal):
    """
    Natural representation of a tensor product of `U_q(\mathfrak{g})`-crystals.
    """
    def e_on_basis(self, i, b):
        """
        Return the action of `e_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self._g.q(i)
        for n in range(k): # We are using the anti-Kashiwara convention
            next = RepresentationFromCrystal.e_on_basis(self, i, b[n])
            for m,c in next:
                ret.append((TP(*(b[:n] + [m] + b[n+1:])), c * coeff))
            coeff *= q**-b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

    def f_on_basis(self, i, b):
        """
        Return the action of `f_i` on the basis element indexed by ``b``.
        """
        TP = self.basis().keys()
        h = TP.weight_lattice_realization().simple_coroots()

        ret = []
        k = len(TP.crystals)
        coeff = self.base_ring().one()
        q = self._g.q(i)
        for n in reversed(range(k)): # We are using the anti-Kashiwara convention
            next = b[n].f(i)
            if next is not None:
                ret.append((TP(*(b[:n] + [next] + b[n+1:])), coeff))
            coeff *= q**b[n].weight().scalar(h[i])
        return self.sum_of_terms(ret)

class GlobalCrystalBasis(CombinatorialFreeModule):
    """
    The global crystal basis of (a tensor product of) highest weight modules
    in `\mathcal{O}^q_{int}`.
    """
    def __init__(self, g, B, prefix='G'):
        """
        Initialize ``self``.
        """
        if g.cartan_type() != B.cartan_type():
            raise ValueError("mismatched Cartan type")
        if B not in HighestWeightCrystals():
            raise ValueError("only for highest weight crystals")
        self._g = g
        CombinatorialFreeModule.__init__(self, g.base_ring(), B, category=CategoryOInt(self._g),
                                         prefix=prefix, bracket=False)

        V = RepresentationFromCrystal(self._g, B)
        self.module_morphism(self.to_natural_basis, codomain=V,
                             triangular='lower').register_as_coercion()

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Global crystal basis of {}".format(self.basis().keys())

    def _latex_(self):
        r"""
        Return a `\LaTeX` representation of ``self``.
        """
        from sage.misc.latex import latex
        return "G\\left( {} \\right)".format(latex(self.basis().keys()))

    def to_natural_basis(self, m):
        """
        Return the basis element indexed by ``m`` in the natural basis
        of the tensor product.

        EXAMPLES:

        We construct the examples in [HK02]_ Section 6.2::

            from sage.algebras.quantum_groups.representations import GlobalCrystalBasis
            sage: Q = QuantumGroup(['A',1])
            sage: B = CrystalOfTableaux(['A',1],shape=[2])
            sage: T = TensorProductOfCrystals(B,B)
            sage: G = GlobalCrystalBasis(Q, T)
            sage: map(lambda x: G.to_natural_basis(x), T)
        """
        B = self.basis().keys()
        R = RepresentationFromCrystal(self._g, B)
        # Special case for a single factor
        if not isinstance(B, TensorProductOfCrystals):
            return R.monomial(m)

        # Construct the maximal vector
        hw, path = m.to_highest_weight()
        cur = R.highest_weight_vector(hw)

        # Now build down the crystal until we reach the desired node
        if len(path) == 0:
            return cur

        last = path[0]
        p = 0
        coeff = self.base_ring().one()
        for i in reversed(path):
            if i != last:
                coeff *= q_factorial(p, self._g.q(last))
                p = 1
                last = i
            else:
                p += 1
            cur = cur.f(i)
        # FIXME: Absorb in all divided powers?
        return cur / (coeff * q_factorial(p, self._g.q(last)))

class QuantumVectorRepresentation(Parent, UniqueRepresentation):
    r"""
    Representation of finite type `U_q(\mathfrak{g})` as matrices.

    INPUT:

    - ``R`` -- the base ring or a matrix Lie algebra `\mathfrak{g}`
    - ``q`` -- the deformation parameter
    - ``ct`` -- a Cartan type if ``R`` is not a matrix Lie algebra
    """
    @staticmethod
    def __classcall_private__(cls, R, q, ct=None):
        """
        Normalize input to ensure a unique reprensetation.
        """
        if isinstance(R, ClassicalMatrixLieAlgebra):
            g = R
        else:
            #R = R.fraction_field()
            g = ClassicalMatrixLieAlgebra(R, ct)
        return super(QuantumVectorRepresentation, cls).__classcall__(cls, g, q)

    def __init__(self, g, q):
        """
        Initialize ``self``.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: V = QuantumVectorRepresentation(R, q, ['A',3])
            sage: TestSuite(V).run()
        """
        self._base_repr = g
        self._n = g._n
        self._q = q
        Parent.__init__(self, base=g.base_ring(), category=AlgebrasWithBasis(g.base_ring()))

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: QuantumVectorRepresentation(R, q, ['A',3])
            Vector representation of the quantum group of
             Special linear Lie algebra of rank 4 over
             Univariate Polynomial Ring in q over Rational Field
        """
        return "Vector representation of the quantum group of {}".format(self._base_repr)

    def _an_element_(self):
        """
        Return an element of ``self``.
        """
        return self.e(0)

    def some_elements(self):
        """
        Return some elements of ``self``.
        """
        return [self.e(0), self.f(0), self.K(0)]

    def base_representation(self):
        """
        Return the base representation of ``self``.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: V = QuantumVectorRepresentation(R, q, ['A',3])
            sage: V.base_representation()
            Special linear Lie algebra of rank 4 over
             Univariate Polynomial Ring in q over Rational Field
        """
        return self._base_repr

    def e(self, i=None):
        """
        Return the generator `e_i`.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: V = QuantumVectorRepresentation(R, q, ['A',3])
            sage: V.e()
            [[0 1 0 0]
             [0 0 0 0]
             [0 0 0 0]
             [0 0 0 0],
             [0 0 0 0]
             [0 0 1 0]
             [0 0 0 0]
             [0 0 0 0],
             [0 0 0 0]
             [0 0 0 0]
             [0 0 0 1]
             [0 0 0 0]]
        """
        if i is None:
            return map(lambda x: self(x.value), self._base_repr._e)
        return self(self._base_repr._e[i].value)

    def f(self, i=None):
        """
        Return the generator `f_i`.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: V = QuantumVectorRepresentation(R, q, ['A',3])
            sage: V.f()
            [[0 0 0 0]
             [1 0 0 0]
             [0 0 0 0]
             [0 0 0 0],
             [0 0 0 0]
             [0 0 0 0]
             [0 1 0 0]
             [0 0 0 0],
             [0 0 0 0]
             [0 0 0 0]
             [0 0 0 0]
             [0 0 1 0]]
        """
        if i is None:
            return map(lambda x: self(x.value), self._base_repr._f)
        return self(self._base_repr._f[i].value)

    def K(self, i=None):
        """
        Return the generator `K_i`.

        EXAMPLES::

            from sage.algebras.quantum_groups.representations import QuantumVectorRepresentation
            sage: R.<q> = QQ[]
            sage: V = QuantumVectorRepresentation(R, q, ['A',3])
            sage: V.K()
            [[  q   0   0   0]
             [  0 1/q   0   0]
             [  0   0   1   0]
             [  0   0   0   1],
             [  1   0   0   0]
             [  0   q   0   0]
             [  0   0 1/q   0]
             [  0   0   0   1],
             [  1   0   0   0]
             [  0   1   0   0]
             [  0   0   q   0]
             [  0   0   0 1/q]]
        """
        if i is None:
            return [self.K(i) for i in range(len(self._base_repr._h))]
        q = self._q
        h = self._base_repr._h[i].value
        MS = h.parent()
        # Make sure we are in the fraction field
        MS = MS.change_ring(MS.base_ring().fraction_field())
        return self(MS( {(i,i): q**h[i,i] for i in range(self._n)} ))

    @cached_method
    def zero(self):
        """
        Return the additive identity `0`.
        """
        return self(self._base_repr.universal_enveloping_algebra().zero())

    @cached_method
    def one(self):
        """
        Return the multiplicative identity `1`.
        """
        return self(self._base_repr.universal_enveloping_algebra().one())

    class Element(ElementWrapper):
        def _add_(self, rhs):
            """
            Add ``self`` to ``rhs``.
            """
            return self.__class__(self.parent(), self.value + rhs.value)

        def __neg__(self, rhs):
            """
            Negate ``self``.
            """
            return self.__class__(self.parent(), -self.value)

        def _mul_(self, rhs):
            """
            Multiply ``self`` to ``rhs``.
            """
            return self.__class__(self.parent(), self.value * rhs.value)

        def _lmul_(self, x):
            """
            Return ``self * x``.
            """
            return self.__class__(self.parent(), self.value._lmul_(x))

        def _rmul_(self, x):
            """
            Return ``x * self``.
            """
            return self.__class__(self.parent(), self.value._rmul_(x))

