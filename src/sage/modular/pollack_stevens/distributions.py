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
from sage.modular.overconvergent.pollack.dist import get_dist_classes, Dist_long

class Distributions(Module):
    def __init__(self, k, p=None, prec_cap=None, base=None, character=None, tuplegen=None, act_on_left=False):
        """
        - ``character`` --
          - None (default)
          - (chi, None)
          - (None, n) (n integral)
          - (chi, n)
          - lambda (for n half-integral use this form)
        """
        if p is not None:
            p = ZZ(p)
        if base is None:
            if p is None: raise ValueError("specify p or a base")
            base = ZpCA(p,prec_cap)
        elif isinstance(base, pAdicGeneric):
            if p is None:
                p = base.prime()
            elif base.prime() != p:
                raise ValueError("p must be the same as the prime of base")
            if prec_cap is None:
                prec_cap = base.precision_cap()
            elif base.precision_cap() != prec_cap:
                raise ValueError("prec_cap must match the precision cap of base")
        elif prec_cap is not None and prec_cap > k+1: # non-classical
            if p is None or not p.is_prime(): raise ValueError("p must be prime for non-classical weight")
        from sage.rings.padics.pow_computer import PowComputer_long
        # should eventually be the PowComputer on ZpCA once that uses longs.
        Dist, WeightKAction = get_dist_classes(p, prec_cap, base)
        self.Element = Dist
        if Dist is Dist_long:
            self.prime_pow = PowComputer_long(p, prec_cap, prec_cap, prec_cap, 0)
        Parent.__init__(self, base)
        self._k = k
        self._p = p
        self._prec_cap = prec_cap
        act = WeightKAction(self, character, tuplegen, act_on_left)
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

    @cached_method
    def basis(self, M):
        V = self.approx_module(M)
        return [self(v) for v in V.basis()]

    def _an_element_(self):
        if self._prec_cap > 1:
            return self([2,1])
        else:
            return self([1])

    #def get_action(self):
    #    return self._act
