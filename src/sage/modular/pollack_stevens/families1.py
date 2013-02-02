#*****************************************************************************
#       Copyright (C) 2013 Nathan Clement <clement.nathan@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

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
    """
    A note on precision in families.  Precision can now be specified by M (as before) and 
    a new parameter D.  This means that mu(1) will be known modulo p^M + w^D, mu(z) will
    known modulo p^(M-1) + w^D, and in general, mu(z^k) will be stored modulo p^(M-k) + w^D
    """
    
    def _show_malformed_dist(self, location_str):
        r"""
        This checks if the distribution is malformed.
        """

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

    def reduce_precision(self, M, D):
        r"""
        Reduce the number of moments and also reduce w
        """
        return self.__class__(self._map.reduce_precision(M, D), self.parent(), construct=True)

    def precision_absolute(self):
        r"""
        Returns the smallest p-adic/moment precision (M) and the smallest w-adic
        precision (D).
        """
        precision_l = []
        for a in self._map:
            precision_l.append(a.precision_absolute())
        minM = min(p[0] for p in precision_l)
        minD = min(p[1] for p in precision_l)
        return (minM,minD)
        """
        We need to think here, we have two kinds of precision
        The simplest, to me, seems to be return the minimum of M and the minimum of D, taken over
        the whole map.
        """

    def specialize(self, weight):
        """
        This is different than the specialize() in PSModularSymbolElement_fam
        """
        raise NotImplementedError

    def _consistency_check(self):
        """
        Check that the map really does satisfy the Manin relations loop (for debugging).
        """
        rels = self.parent()._grab_relations()
        # TODO: no clue how to do this until this object fully works again...
        raise NotImplementedError
