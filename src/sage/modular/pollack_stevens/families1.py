from sage.structure.element import ModuleElement
from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.rings.integer_ring import ZZ
from manin_map import ManinMap
import operator
from sage.misc.cachefunc import cached_method
from sage.rings.padics.factory import Qp
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.arith import next_prime
from sage.rings.infinity import infinity
from sage.misc.misc import verbose
from sage.rings.padics.precision_error import PrecisionError

from sage.categories.action import Action
from sage.modular.pollack_stevens.modsym import PSModularSymbolElement

class PSModularSymbolElement_fam(PSModularSymbolElement):

    def _show_malformed_dist(self, location_str):
        malformed = []
        gens = self.parent().source().gens()
        for j, g in enumerate(gens):
            val = self._map[g]
            if val._is_malformed():
                malformed.append((j, val))
        return location_str + ": (%s/%s malformed)%s"%(len(malformed), len(gens), ", %s -- %s"%(malformed[0][0], str(malformed[0][1])) if len(malformed) > 0 else "")

    def reduce_precision(self, M):
        r"""
        Only holds on to `M` moments of each value of self
        """
        return self.__class__(self._map.reduce_precision(M), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the number of moments of each value of self
        """
        return min([a.precision_absolute() for a in self._map])

    def specialize(self, new_base_ring=None):
        r"""
        Returns the underlying classical symbol of weight `k` -- i.e.,
        applies the canonical map `D_k --> Sym^k` to all values of
        self.
        
        EXAMPLES::

            sage: D = Distributions(0, 5, 10);  M = PSModularSymbols(Gamma0(2), coefficients=D); M
            Space of overconvergent modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Space of 5-adic distributions with k=0 action and precision cap 10
            sage: f = M(1)
            sage: f.specialize()
            Modular symbol with values in Sym^0 Z_5^2
            sage: f.specialize().values()
            [1 + O(5^10), 1 + O(5^10)]
            sage: f.values()
            [1, 1]
            sage: f.specialize().parent()
            Space of modular symbols for Congruence Subgroup Gamma0(2) with sign 0 and values in Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module()
            Sym^0 Z_5^2
            sage: f.specialize().parent().coefficient_module().is_symk()
            True

            sage: f.specialize(QQ)
            Modular symbol with values in Sym^0 Q^2
            sage: f.specialize(QQ).values()
            [1, 1]
            sage: f.specialize(QQ).parent().coefficient_module()
            Sym^0 Q^2
        """
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return self.__class__(self._map.specialize(new_base_ring),
                              self.parent()._specialize_parent_space(new_base_ring), construct=True)

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        """
        rels = self.parent()._grab_relations()
        # TODO: no clue how to do this until this object fully works again...
        raise NotImplementedError
