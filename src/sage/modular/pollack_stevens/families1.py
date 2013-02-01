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
        """
        I guess, it should check whether the Family is given correctly (in the right filtration)
        """

    def reduce_precision(self, M, w):
        r"""
        Reduce the number of moments and also reduce w
        """
        return self.__class__(self._map.reduce_precision(M, W), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the number of moments and w
        """
        return min([a.precision_absolute() for a in self._map])
        """
        We need to think here, we have two kinds of precision
        """

    def specialize(self, new_base_ring=None):
        """
        We want to specialize to k.
        """

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        """
        rels = self.parent()._grab_relations()
        # TODO: no clue how to do this until this object fully works again...
        raise NotImplementedError
