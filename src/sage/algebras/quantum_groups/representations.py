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
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.highest_weight_crystals import HighestWeightCrystals

from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
from sage.algebras.quantum_groups.q_numbers import q_factorial

class RepresentationFromCrystal(CombinatorialFreeModule):
    """
    A representation of a quantum group in the global (canonical) basis
    constructed from a (local) crystal (basis).

    INPUT:

    - ``q`` -- the variable `q`
    - ``crystal`` -- a crystal
    - ``prefix`` -- (default: ``'v'``) the prefix for the basis
    - ``base`` -- (optional) the base ring
    """
    @staticmethod
    def __classcall_private__(cls, q, crystal, prefix='v', base=None):
        """
        Normalize input to ensure a unique representation.
        """
        if base is None:
            base = q.parent()
        # Make sure it's in the fraction field
        base = base.fraction_field()
        q = base(q)

        if isinstance(crystal, TensorProductOfCrystals):
            return RepresentationFromTensorProductOfCrystals(q, crystal, prefix)
        return super(RepresentationFromCrystal, cls).__classcall__(cls, q, crystal, prefix)

    def __init__(self, q, crystal, prefix='v'):
        """
        Initialize ``self``.
        """
        #if crystal not in HighestWeightCrystals():
        #    raise ValueError("only for highest weight crystals")
        self._q = q
        CombinatorialFreeModule.__init__(self, q.parent(), crystal,
                                         prefix=prefix, bracket=False,
                                         category=ModulesWithBasis(R))

    @cached_method
    def highest_weight_vector(self, b, index_set=None):
        """
        Return the highest weight vector corresponding to the highest weight
        crystal element ``b``.
        """
        if index_set is None:
            if not b.is_highest_weight():
                raise ValueError("{} is not highest weight".format(b))
            index_set = self._basis_keys.cartan_type().index_set()
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
        q = self._q
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
        q = self._q
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
        q = self._q
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

