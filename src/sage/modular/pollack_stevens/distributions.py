"""
Spaces of Distributions

"""

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
from sage.structure.factory import UniqueFactory
from sage.structure.unique_representation import UniqueRepresentation
import operator
import sage.rings.ring as ring

# Need this to be pickleable
class _default_tuplegen(UniqueRepresentation):
    """
    Callable object that turns matrices into 4-tuples.

    EXAMPLES::

        sage: sage.modular.pollack_stevens.distributions._default_tuplegen()
        <sage.modular.pollack_stevens.distributions._default_tuplegen object at 0x...>
    """
    def __call__(self, g):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: T = sage.modular.pollack_stevens.distributions._default_tuplegen()
            sage: T(matrix(ZZ,2,[1..4]))
            (1, 2, 3, 4)
        """
        return g[0,0], g[0,1], g[1,0], g[1,1]

class Distributions_factory(UniqueFactory):
    def create_key(self, k, p=None, prec_cap=None, base=None, symk=None, character=None, tuplegen=None, act_on_left=False):
        """
        INPUT:

        - `k` -- nonnegative integer
        - `p` -- prime number or None
        - ``prec_cap`` -- positive integer or None
        - ``base`` -- ring or None
        - ``symk`` -- bool or None
        - ``character`` -- a dirichlet character or None
        - ``tuplegen`` -- None or callable that turns 2x2 matrices into a 4-tuple
        - ``act_on_left`` -- bool (default: False)

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: Distributions(20, 3, 10)              # indirect doctest
            Space of 3-adic distributions with k=20 action and precision cap 10
        """
        k = ZZ(k)
        if tuplegen is None:
            tuplegen = _default_tuplegen()
        if p is None:
            try:
                p = base.prime()
            except AttributeError:
                raise ValueError("You must specify a prime")
        else:
            p = ZZ(p)
        if base is None:
            if prec_cap is None:
                base = ZpCA(p)
            else:
                base = ZpCA(p, prec_cap)
        if prec_cap is None:
            try:
                prec_cap = base.precision_cap()
            except AttributeError:
                raise ValueError("You must specify a base or precision cap")
        return (k, p, prec_cap, base, character, tuplegen, act_on_left, symk)

    def create_object(self, version, key):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: Distributions(0, 7, 5)              # indirect doctest
            Space of 7-adic distributions with k=0 action and precision cap 5
        """
        return Distributions_class(*key)

class Symk_factory(UniqueFactory):
    def create_key(self, k, base=None, character=None, tuplegen=None, act_on_left=False):
        k = ZZ(k)
        if tuplegen is None:
            tuplegen = _default_tuplegen()
        prec_cap = k+1
        if base is None:
            base = QQ
        if isinstance(base, pAdicGeneric):
            p = base.prime()
        return (k, base, character, tuplegen, act_on_left)

    def create_object(self, version, key):
        return Symk_class(*key)

Distributions = Distributions_factory('Distributions')
Symk = Symk_factory('Symk')

class Distributions_abstract(Module):
    """
    Parent object for distributions.

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions
        sage: Distributions(2, 17, 100)
        Space of 17-adic distributions with k=2 action and precision cap 100
    """
    def __init__(self, k, p=None, prec_cap=None, base=None, character=None, tuplegen=None, act_on_left=False, symk=False):
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

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(2, 3, 5); D
            Space of 3-adic distributions with k=2 action and precision cap 5
            sage: type(D)
            <class 'sage.modular.pollack_stevens.distributions.Distributions_class'>

        p must be a prime, but p=6 below, which is not prime::

            sage: Distributions(k=0, p=6, prec_cap=10)
            Traceback (most recent call last):
            ...
            ValueError: p must be prime
        """
        if not isinstance(base, ring.Ring):
            raise TypeError("base must be a ring")
        from sage.rings.padics.pow_computer import PowComputer_long
        # should eventually be the PowComputer on ZpCA once that uses longs.
        Dist, WeightKAction = get_dist_classes(p, prec_cap, base, symk)
        self.Element = Dist
        if Dist is Dist_long:
            self.prime_pow = PowComputer_long(p, prec_cap, prec_cap, prec_cap, 0)
        Parent.__init__(self, base)
        self._k = k
        self._p = p
        self._prec_cap = prec_cap
        self._character = character
        self._symk = symk
        act = WeightKAction(self, character, tuplegen, act_on_left)
        self._act = act
        self._populate_coercion_lists_(action_list=[iScale(self, act_on_left), act])

    def prime(self):
        """
        Return prime `p` such that this is a space of `p`-adic distributions.

        In case this space is Symk of a non-padic field, this makes no
        sense, and we raise a ValueError.

        OUTPUT:

        - bool

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7); D
            Space of 7-adic distributions with k=0 action and precision cap 20
            sage: D.prime()
            7
            sage: D = Symk(4, base=GF(7)); D
            Sym^4 (Finite Field of size 7)^2
            sage: D.prime()
            Traceback (most recent call last):
            ...
            ValueError: not a space of p-adic distributions

        But Symk of a `p`-adic field does work::

            sage: D = Symk(4, base=Qp(7)); D
            Sym^4 Q_7^2
            sage: D.prime()
            7
            sage: D.is_symk()
            True
        """
        if self._p is None:
            raise ValueError, "not a space of p-adic distributions"
        return self._p

    def weight(self):
        """
        Return the weight of this distribution space.  The standard
        caveat applies, namely that the weight of `Sym^k` is
        defined to be `k`, not `k+2`.

        OUTPUT:

        - nonnegative integer

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7); D
            Space of 7-adic distributions with k=0 action and precision cap 20
            sage: D.weight()
            0
            sage: Distributions(389, 7).weight()
            389
        """
        return self._k

    def precision_cap(self):
        """
        Return the precision cap on distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 10); D
            Space of 7-adic distributions with k=0 action and precision cap 10
            sage: D.precision_cap()
            10
            sage: D = Symk(389, base=QQ); D
            Sym^389 Q^2
            sage: D.precision_cap()
            390
        """
        return self._prec_cap

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
            ValueError: M must be less than or equal to the precision cap
            sage: D.approx_module(-1)
            Traceback (most recent call last):
            ...
            ValueError: rank (=-1) must be nonnegative
        """
        if M is None:
            M = self._prec_cap
        elif M > self._prec_cap:
            raise ValueError("M must be less than or equal to the precision cap")
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
            ValueError: M must be less than or equal to the precision cap
        """
        return self(self.approx_module(M).random_element())

    def clear_cache(self):
        """
        Clear some caches that are created only for speed purposes.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 10)
            sage: D.clear_cache()
        """
        self.approx_module.clear_cache()
        self._act.clear_cache()

    @cached_method
    def basis(self, M=None):
        """
        Return a basis for this space of distributions.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.basis()
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]

            sage: D = Symk(3, base=QQ); D
            Sym^3 Q^2
            sage: D.basis()
            [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
        """
        V = self.approx_module(M)
        return [self(v) for v in V.basis()]

    def _an_element_(self):
        """
        Return a typical element of self.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Symk(3, base=QQ); D
            Sym^3 Q^2
            sage: D.an_element()                  # indirect doctest
            (2, 1)
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.an_element()
            (2, 1)
        """
        if self._prec_cap > 1:
            return self([2,1])
        else:
            return self([1])

    def zero_element(self, M=None):
        """
        Return zero element in the M-th approximating module.

        INPUT:

        - `M` -- None (default), or a nonnegative integer, less than or equal to the precision cap

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.zero_element()
            (0, 0, 0, 0)
            sage: D.zero_element(0)
            ()
            sage: D.zero_element(1)
            0
            sage: D.zero_element(2)
            (0, 0)
            sage: D.zero_element(3)
            (0, 0, 0)
            sage: D.zero_element(4)
            (0, 0, 0, 0)
            sage: D.zero_element(5)
            Traceback (most recent call last):
            ...
            ValueError: M must be less than or equal to the precision cap
        """
        return self(self.approx_module(M)(0))

class Symk_class(Distributions_abstract):
    def __init__(self, k, base, character, tuplegen, act_on_left):
        if hasattr(base, 'prime'):
            p = base.prime()
        else:
            p = None
        Distributions_abstract.__init__(self, k, p, k+1, base, character, tuplegen, act_on_left)

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Symk
            sage: Distributions(0, 5, 10)._repr_()
            'Space of 5-adic distributions with k=0 action and precision cap 10'
            sage: Distributions(0, 5, 10)
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: Symk(0)
            Sym^0 Q^2
        """
        # TODO: maybe account for character, etc.
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

    def is_symk(self):
        """
        Whether or not this distributions space is Sym^k (ring).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(4, 17, 10); D
            Space of 17-adic distributions with k=4 action and precision cap 10
            sage: D.is_symk()
            False
            sage: D = Symk(4); D
            Sym^4 Q^2
            sage: D.is_symk()
            True
            sage: D = Symk(4, base=GF(7)); D
            Sym^4 (Finite Field of size 7)^2
            sage: D.is_symk()
            True
        """
        return True

    def change_ring(self, new_base_ring):
        """
        Return a Symk with the same k but a different base ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.base_ring()
            7-adic Ring with capped absolute precision 4
            sage: D2 = D.change_ring(QpCR(7)); D2
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D2.base_ring()
            7-adic Field with capped relative precision 20
        """
        return Symk(k=self._k, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

    def lift(self, p=None, M=None, new_base_ring=None):
        """
        Return distribution space that contains lifts with given p,
        precision cap M, and base ring new_base_ring.

        INPUT:

        - `p` -- prime or None
        - `M` -- nonnegative integer or None
        - ``new_base_ring`` -- ring or None

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Symk(0, Qp(7)); D
            Sym^0 Q_7^2
            sage: D.lift(M=20)
            Space of 7-adic distributions with k=0 action and precision cap 20
            sage: D.lift(p=7, M=10)
            Space of 7-adic distributions with k=0 action and precision cap 10
            sage: D.lift(p=7, M=10, new_base_ring=QpCR(7,15)).base_ring()
            7-adic Field with capped relative precision 15
        """
        if self._character is not None:
            # need to change coefficient ring for character
            raise NotImplementedError
        if M is None:
            M = self._prec_cap + 1
        if p is None:
            try:
                p = self.base_ring().prime()
            except AttributeError:
                raise ValueError("You must specify a prime")
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return Distributions(k=self._k, p=p, prec_cap=M, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

class Distributions_class(Distributions_abstract):
    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: Distributions(0, 5, 10)._repr_()
            'Space of 5-adic distributions with k=0 action and precision cap 10'
            sage: Distributions(0, 5, 10)
            Space of 5-adic distributions with k=0 action and precision cap 10
            sage: Symk(0)
            Sym^0 Q^2
        """
        # TODO: maybe account for character, etc.
        return "Space of %s-adic distributions with k=%s action and precision cap %s"%(self._p, self._k, self._prec_cap)

    def is_symk(self):
        """
        Whether or not this distributions space is Sym^k (ring).

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(4, 17, 10); D
            Space of 17-adic distributions with k=4 action and precision cap 10
            sage: D.is_symk()
            False
            sage: D = Symk(4); D
            Sym^4 Q^2
            sage: D.is_symk()
            True
            sage: D = Symk(4, base=GF(7)); D
            Sym^4 (Finite Field of size 7)^2
            sage: D.is_symk()
            True
        """
        return False

    def change_ring(self, new_base_ring):
        """
        Return space of distributions like this one, but with the base ring changed.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.base_ring()
            7-adic Ring with capped absolute precision 4
            sage: D2 = D.change_ring(QpCR(7)); D2
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D2.base_ring()
            7-adic Field with capped relative precision 20
        """
        return Distributions(k=self._k, p=self._p, prec_cap=self._prec_cap, base=new_base_ring, character=self._character, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())

    def specialize(self, new_base_ring=None):
        """
        Return distribution space got by specializing to Sym^k, over
        the new_base_ring.  If new_base_ring is not given, use current
        base_ring.

        EXAMPLES::

            sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
            sage: D = Distributions(0, 7, 4); D
            Space of 7-adic distributions with k=0 action and precision cap 4
            sage: D.is_symk()
            False
            sage: D2 = D.specialize(); D2
            Sym^0 Z_7^2
            sage: D2.is_symk()
            True
            sage: D2 = D.specialize(QQ); D2
            Sym^0 Q^2
        """
        if self._character is not None:
            raise NotImplementedError
        if new_base_ring is None:
            new_base_ring = self.base_ring()
        return Symk(k=self._k, base=new_base_ring, tuplegen=self._act._tuplegen, act_on_left=self._act.is_left())
