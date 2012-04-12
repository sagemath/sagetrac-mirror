from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.rings.padics.factory import ZpCA
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.categories.action import PrecomposedAction
from sage.matrix.all import MatrixSpace
from sage.rings.fast_arith import prime_range
from sage.modular.overconvergent.pollack.dist import get_dist_classes

class Distributions(Module):
    def __init__(self, k, p, prec_cap, base=None, character=None):
        if base is None:
            base = ZpCA(p,prec_cap)
        else:
            assert (isinstance(base, pAdicGeneric) and base.prime() == p) or (p.is_prime() and (base is ZZ or base is QQ))
        from dist import Dist_vector, WeightKAction_vector, Dist_long, WeightKAction_long, 
        from sage.rings.padics.pow_computer import PowComputer_long
        # should eventually be the PowComputer on ZpCA once that uses longs.
        p = ZZ(p)
        if 7*p**
        self._element_constructor_ = Dist_vector
        Parent.__init__(self, base)
        self._k = k
        self._p = ZZ(p)
        self._prec_cap = prec_cap
        self._approx_modules = {} # indexed by precision
        act = WeightKAction_vector(self, character)
        self._act = act
        self._populate_coercion_lists_(action_list=[act])

    @cached_method
    def approx_module(self, M):
        assert M <= self._prec_cap
        return self.base_ring()**M

    def random_element(self, M):
        return self(self.approx_module(M).random_element())

    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()
