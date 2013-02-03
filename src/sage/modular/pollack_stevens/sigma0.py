r"""
The monoid `\Sigma_0(p)`
"""

from sage.matrix.matrix_integer_2x2 import MatrixSpace_ZZ_2x2
from sage.matrix.matrix_space import MatrixSpace
from sage.structure.factory import UniqueFactory
from sage.structure.element import MonoidElement
from sage.categories.monoids import Monoids
from sage.structure.parent import Parent
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity
from sage.structure.unique_representation import UniqueRepresentation

# Need this to be pickleable
class _default_tuplegen(UniqueRepresentation):
    """
    Callable object that turns matrices into 4-tuples.

    EXAMPLES::

        sage: A = sage.modular.pollack_stevens.distributions._default_tuplegen(); A
        <sage.modular.pollack_stevens.distributions._default_tuplegen object at 0x...>
        sage: TestSuite(A).run()
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

class Sigma0_factory(UniqueFactory):

    def create_key(self, p, base_ring=ZZ,tuplegen = None):
        p = ZZ(p)
        if not (p == 0 or p.is_prime()):
            raise ValueError("Primes should be prime")
        if tuplegen is None:
            tuplegen = _default_tuplegen()

        if base_ring not in (QQ, ZZ):
            try:
                if p!= 0 and base_ring.prime() != p:
                    raise ValueError("Base ring and prime do not match")
            except AttributeError:
                raise ValueError("Base ring must be QQ, ZZ or a p-adic field")
        return (p, base_ring,tuplegen)

    def create_object(self, version, key):
        return Sigma0_class(*key)

Sigma0 = Sigma0_factory('sage.modular.pollack_stevens.sigma0.Sigma0')

class Sigma0Element(MonoidElement):

    def __init__(self, parent, mat):
        self._mat = mat
        MonoidElement.__init__(self, parent)

    def __hash__(self):
        return hash(self.matrix())

    def _mul_(self, other):
        return self.parent()(self._mat * other._mat, check=False)

    def __cmp__(self, other):
        return cmp(self._mat, other._mat)

    def _repr_(self):
        return self.matrix().__repr__()

    def matrix(self):
        return self._mat
    
    def __getitem__(self, *args):
        return self._mat.__getitem__(*args)

    def inverse(self):
        return self.parent()(self._mat._invert_unit())

class Sigma0_class(Parent):

    Element = Sigma0Element

    def __init__(self, p, base_ring,tuplegen):
        self._prime = p
        self._base_ring = base_ring
        self._tuplegen = tuplegen
        if base_ring == ZZ:
            self._matrix_space = MatrixSpace_ZZ_2x2()
        else:
            self._matrix_space = MatrixSpace(base_ring, 2)
        Parent.__init__(self, category=Monoids())

    def _an_element_(self):
        return self([1,0,0,1])

    def __cmp__(self, other):
        return cmp(type(self), type(other)) \
            or cmp(self.prime(), other.prime())

    def prime(self):
        return self._prime

    def base_ring(self):
        return self._base_ring

    def _coerce_map_from_(self, other):
        r"""
        The *only* thing that coerces canonically into `\Sigma_0` is another
        `\Sigma_0`. It is *very bad* if integers are allowed to coerce in, as
        this leads to a noncommutative coercion diagram.
        """
        if isinstance(other, Sigma0_class) \
            and self.prime().divides(other.prime()) \
            and self.base_ring().has_coerce_map_from(other.base_ring()): 
            return True
        else:
            return False

    def _element_constructor_(self, x, check=True):
        if isinstance(x, Sigma0Element):
            x = x.matrix()
        if check:
            x = self._matrix_space(x)
            if self.prime() != 0:
                a,b,c,d = self._tuplegen(x)
                if c.valuation(self.prime()) <= 0:
                    raise ValueError("p must divide c")
                if a.valuation(self.prime()) != 0:
                    raise ValueError("a must be a p-adic unit")
            if x.det() == 0:
                raise ValueError("matrix must be nonsingular")
        x.set_immutable()
        return self.element_class(self, x)

    def _repr_(self):
        if self.prime() == 0:
            if self.base_ring() is ZZ:
                return 'Monoid of invertible 2x2 integer matrices'
            elif self.base_ring() is QQ:
                return 'Monoid of invertible 2x2 rational matrices'
        return 'Monoid of invertible 2x2 matrices over %s, upper-triangular modulo %s' % (self.base_ring(), self.prime())
