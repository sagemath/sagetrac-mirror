from sage.modules.module import Module
from sage.structure.parent import Parent
from sage.rings.padics.factory import ZpCA
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.categories.action import PrecomposedAction
from sage.structure.coerce_actions import LeftModuleAction, RightModuleAction
from sage.matrix.all import MatrixSpace
from sage.rings.fast_arith import prime_range
from sage.modular.pollack_stevens.dist import get_dist_classes, Dist_long, iScale
import operator

class Distributions(Module):
    """
    Parent object for distributions.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions
        sage: Distributions(2, 17, 100)
        Space of 17-adic distributions with k=2 action and precision cap 100
    """
    def __init__(self, k, p=None, prec_cap=None, base=None, character=None, tuplegen=None, act_on_left=False):
        """
        INPUT:

        - `k` -- integer; k is the usual modular forms weight minus 2
        - `p` -- None or prime
        - ``prec_cap`` -- None or positive integer
        - ``base`` -- None or TODO
        - ``character`` --
          - None (default)
          - (chi, None)
          - (None, n) (n integral)
          - (chi, n)
          - lambda (for n half-integral use this form)
        - ``tuplegen`` -- None or TODO
        - ``act_on_left`` -- bool (default: False)

        EXAMPLES::

            sage: D = Distributions(2, 3, 5); D
            Space of 3-adic distributions with k=2 action and precision cap 5
            sage: type(D)
            <class 'sage.modular.pollack_stevens.distributions.Distributions'>

        p must be a prime, but p=6 below, which is not prime::
        
            sage: Distributions(k=0, p=6, prec_cap=10)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
        """
        if p is not None:
            p = ZZ(p)
        if base is None:
            if p is None:
                if prec_cap is None:
                    base = QQ
                    prec_cap = k + 1
                else:
                    raise ValueError("prec cap cannot be specified if p and base are both None")
            else:
                if prec_cap is None:
                    prec_cap = 20
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
        if prec_cap is None:
            self._prec_cap = k + 1
        else:
            self._prec_cap = prec_cap
        act = WeightKAction(self, character, tuplegen, act_on_left)
        self._act = act
        self._populate_coercion_lists_(action_list=[iScale(self, act_on_left), act])

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: Distributions(0, 5, 10)._repr_()
            'Space of 5-adic distributions with k=0 action and precision cap 10'
            sage: Distributions(0, 5, 10)
            Space of 5-adic distributions with k=0 action and precision cap 10

            sage: Distributions(0)
            Sym^0 Q^2
        """
        # TODO: maybe account for character, etc.
        if self._p is None:
            if self.base_ring() is QQ:
                V = 'Q^2'
            elif self.base_ring() is ZZ:
                V = 'Z^2'
            elif isinstance(self.base_ring(), pAdicGeneric) and self.base_ring().degree() == 1:
                if self.base_ring().is_field():
                    V = 'Q_%s^2'%(self._p)
                else:
                    V = 'Z_%s^2'%(self._p)
            else:
                V = '(%s)^2'%(self.base_ring())
            return "Sym^%s %s"%(self._k, V)
        else:
            return "Space of %s-adic distributions with k=%s action and precision cap %s"%(
                self._p, self._k, self._prec_cap)

    @cached_method
    def approx_module(self, M=None):
        """
        Return the M-th approximation module, or if M is not specified,
        return the largest approximation module.

        INPUT::

        - `M` -- None or nonnegative integer that is at most the precision cap

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(0, 5, 10)
            sage: D.approx_module()
            Ambient free module of rank 10 over the principal ideal domain 5-adic Ring with capped absolute precision 10
            sage: D.approx_module(1)
            Ambient free module of rank 1 over the principal ideal domain 5-adic Ring with capped absolute precision 10
            sage: D.approx_module(0)
            Ambient free module of rank 0 over the principal ideal domain 5-adic Ring with capped absolute precision 10

        Note that M must be at most the precision cap, and must be nonnegative::
        
            sage: D.approx_module(11)
            Traceback (most recent call last):
            ...
            ValueError: M must be less than the precision cap
            sage: D.approx_module(-1)
            Traceback (most recent call last):
            ...
            ValueError: rank (=-1) must be nonnegative
        """
        if M is None:
            M = self._prec_cap
        elif M > self._prec_cap:
            raise ValueError("M must be less than the precision cap")
        return self.base_ring()**M

    def random_element(self, M=None):
        """
        Return a random element of the M-th approximation module.

        INPUT:

        - `M` -- None or a nonnegative integer

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions
            sage: D = Distributions(0, 5, 10)
            sage: D.random_element()
            (..., ..., ..., ..., ..., ..., ..., ..., ..., ...)
            sage: D.random_element(0)
            ()
            sage: D.random_element(5)
            (..., ..., ..., ..., ...)
            sage: D.random_element(-1)
            Traceback (most recent call last):
            ...
            ValueError: rank (=-1) must be nonnegative
            sage: D.random_element(11)
            Traceback (most recent call last):
            ...
            ValueError: M must be less than the precision cap
        """
        return self(self.approx_module(M).random_element())

    def clear_cache(self):
        self.approx_module.clear_cache()
        self._act.clear_cache()

    @cached_method
    def basis(self, M=None):
        V = self.approx_module(M)
        return [self(v) for v in V.basis()]

    def _an_element_(self):
        if self._prec_cap > 1:
            return self([2,1])
        else:
            return self([1])

    def zero_element(self, M=None):
        return self(self.approx_module(M)(0))

#    def _get_action_(self, S, op, self_on_left):
#        if S is self.base_ring():
#            if self_on_left:
#                return LeftModuleAction(self.base_ring(), self)
#            else:
#                return RightModuleAction(self.base_ring(), self)
#        f = self.base_ring().coerce_map_from(S)
#        if op is operator.mul and f is not None:
#            A = self.get_action(self.base_ring(), op, self_on_left)
#            if self_on_left:
#                return PrecomposedAction(A, f, None)
#            else:
#                return PrecomposedAction(A, None, f)

    #def get_action(self):
    #    return self._act
