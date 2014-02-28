r"""
Arithmetic (finite index) subgroups of `{\rm SL}_2(\mathcal{O})` and
`{\rm GL}_2(\mathcal{O})`) where `\mathcal{O}` is an (maximal) order
in a number field `K`.

AUTHORS:

- Fredrik Stromberg (2013): initial version based on arithgroup_generic.py
"""

#############################################################################
#
#       Copyright (C) 2009, The Sage Group -- http://www.sagemath.org/
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
#############################################################################

from sage.rings.all import ZZ, Integer, is_Ring, QQ
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from sage.groups.matrix_gps.linear import LinearMatrixGroup_generic


class ArithmeticSubgroup_NF_class(LinearMatrixGroup_generic):
    r"""
    Base class for arithmetic subgroups of `{\rm SL}_2(K)`. Not
    intended to be used directly, but still includes quite a few
    general-purpose routines which compute data about an arithmetic subgroup.

    """
    def __init__(self, ring, group='SL', name='', ltx=''):
        r"""
        Standard init routine.

        INPUT:

        - `ring` -- ring
        - `special` -- bool (True for SL and False for GL)
        - `name` -- string.
        - `ltx` -- string

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O)
            Hilbert modular group `SL_{2}(O)`
            sage: HilbertModularGroup(O,group='GL')
            Hilbert modular group `GL_{2}(O)`
        """
        degree = 2
        if group not in ['SL', 'GL', 'GL+']:
            raise NotImplementedError("Only groups SL, GL and GL+ implemented!")
        self._group = group
        if name == '':
            name = 'Arithmetic Subgroup of the Group {0} of degree {1} over {2}'.format(group, degree, ring)
        if ltx == '':
            ltx = 'GL({0}, {1})'.format(degree, latex(ring))
        assert is_Ring(ring)
        if group == 'SL':
            special = True
        else:
            special = False
        super(ArithmeticSubgroup_NF_class, self).__init__(Integer(degree),
                                                          ring, special,
                                                          name, ltx)
        self._base_ring = ring
        self._cusps = None
        self._fundamental_domain = None
        if ring == ZZ:
            self._number_field = QQ
        elif hasattr(ring, "number_field"):
            self._number_field = ring.number_field()
        else:
            raise NotImplementedError

    def __reduce__(self):
        r"""
        Used for pickling self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H=sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)
            sage: H.__reduce__()
            Traceback (most recent call last):
            ...
            NotImplementedError: all subclasses must define a __reduce__ method
        """
        raise NotImplementedError("all subclasses must define a __reduce__ method")

    def __hash__(self):
        r"""
        Return a hash of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: hash(HilbertModularGroup(O))
            -6209004311378228744
        """
        return hash(str(self))

    def group(self):
        r"""
        Return the group type of self. Either 'SL', 'GL' or 'GL+'

        EXAMPLES:

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O).group()
            'SL'
            sage: HilbertModularGroup(O,group='GL').group()
            'GL'
        """
        return self._group

    def is_special(self):
        r"""
        Check if self is special, i.e. of type 'SL' or not.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.is_special()
            True
            sage: HilbertModularGroup(QuadraticField(41).ring_of_integers(),group='GL').is_special()
            False
        """
        return self._special

    def number_field(self):
        r"""
        Return the number field over which self is defined.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: HilbertModularGroup(QuadraticField(41).ring_of_integers()).number_field()
            Number Field in a with defining polynomial x^2 - 41
        """
        return self._number_field

    def coset_reps(self, G=None):
        r"""
        Return right coset representatives for self \\ G, where G is another
        arithmetic subgroup that contains self.  If G = None, default to G =
        SL_2(O) where O is the ring of integers of the number field
        over which self is defined.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H=sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)
            sage: H.coset_reps()
            Traceback (most recent call last):
            ...
            NotImplementedError: All subclasses should implement coset representatives
            sage: G = HilbertModularGroup(O)
            sage: G.coset_reps()
            [1 0]
            [0 1]
        """
        raise NotImplementedError("All subclasses should implement coset representatives")

    def fundamental_domain(self):
        r"""
        Return a fundamental domain of self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.fundamental_domain()
            Traceback (most recent call last):
            ...
            NotImplementedError:
        """
        if self._fundamental_domain is not None:
            return self._fundamental_domain
        raise NotImplementedError

    def nu(self, order=2):
        r"""
        Return the number of orbits of elliptic points of given order
        for this arithmetic subgroup.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.nu(6)
            Traceback (most recent call last):
            ...
            NotImplementedError:
        """
        raise NotImplementedError

    def nu2(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.nu2()
            Traceback (most recent call last):
            ...
            NotImplementedError:
        """
        return self.nu(2)

    def nu3(self):
        r"""
        Return the number of orbits of elliptic points of order 3 for this
        arithmetic subgroup.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.nu3()
            Traceback (most recent call last):
            ...
            NotImplementedError:
        """
        return self.nu(3)

    def orders_of_elliptic_elements(self):
        r"""
        Returns the possible orders of elliptic elements.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.orders_of_elliptic_elements()
            Traceback (most recent call last):
            ...
            NotImplementedError: Should be implemented in subclasses!
        """
        raise NotImplementedError("Should be implemented in subclasses!")

    def __cmp__(self, other):
        r"""
        Compare self to other.

        ..NOTE: This function must be overridden by all subclasses.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).__cmp__(G)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __contains__(self, elt):
        r"""
        Check if self contains elt

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)
            sage: 1 in G
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        if elt.parent() == self:
            return True
        raise NotImplementedError

    def is_abelian(self):
        r"""
        Return True if this arithmetic subgroup is abelian.

        Since arithmetic subgroups are always nonabelian, this always
        returns False.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.is_abelian()
            False
        """
        return False

    def is_finite(self):
        r"""
        Return True if this arithmetic subgroup is finite.

        Since arithmetic subgroups are always infinite, this always
        returns False.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.is_finite()
            False
        """
        return False

    def is_subgroup(self, right):
        r"""
        Return True if self is a subgroup of right, and False otherwise. For
        generic arithmetic subgroups this is done by the absurdly slow
        algorithm of checking all of the generators of self to see if they are
        in right.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.is_subgroup(G)
            True
        """
        if self == right:
            return True
        # ridiculously slow generic algorithm
        w = self.gens()
        for g in w:
            if not (g in right):
                return False
        return True

    def is_normal(self, G=None):
        r"""
        Return True precisely if this subgroup is a normal subgroup of G

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: G = HilbertModularGroup(O)
            sage: G.is_normal()
            True
        """
        if G is None:
            G = ArithmeticSubgroup_NF_class(self.base_ring(), self.group())
        if self.index(G) == 1:
            return True
        for x in self.gens():
            for y in G.gens():
                if y * G(x) * (~y) not in self:
                    return False
        return True

    def is_odd(self):
        r"""
        Return True precisely if this subgroup does not contain the
        matrix -1.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O).is_odd()
            False
            """
        return not self.is_even()

    def is_even(self):
        r"""
        Return True precisely if this subgroup contains the matrix -1.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O).is_even()
            True
        """
        minus_one = self([-1, 0, 0, -1])
        return not minus_one in self

    def to_even_subgroup(self):
        r"""
        Return the smallest even subgroup of `SL(2, \ZZ)` containing self.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O).to_even_subgroup()
            Hilbert modular group `SL_{2}(O)`
        """
        if self.is_even():
            return self
        raise NotImplementedError

    def order(self):
        r"""
        Return the number of elements in this arithmetic subgroup.

        Since arithmetic subgroups are always infinite, this always returns
        infinity.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: HilbertModularGroup(O).order()
            +Infinity
        """
        from sage.rings.infinity import infinity
        return infinity

    def reduce_cusp(self, c):
        r"""
        Given a cusp `c \in \mathbb{P}^1(\QQ)`, return the unique reduced cusp
        equivalent to c under the action of self, where a reduced cusp is an
        element `\tfrac{r}{s}` with r,s coprime non-negative integers, s as
        small as possible, and r as small as possible for that s.

        ..NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: c = NFCusp(K,0)
            sage: HilbertModularGroup(O).reduce_cusp(c)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def cusps(self):
        r"""
        Return a sorted list of inequivalent cusps for self, i.e. a set of
        representatives for the orbits of self on `\mathbb{P}^1(\QQ)`.

        These should be returned in a reduced form where this makes sense.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).cusps()
            Traceback (most recent call last):
            ...
            NotImplementedError
            """
        if self._cusps is None:
            self._cusps = self._find_cusps()
        return self._cusps

    def _find_cusps(self):
        r"""
        Calculate a list of inequivalent cusps.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)._find_cusps()
            Traceback (most recent call last):
            ...
            NotImplementedError

        ..NOTE: This should be implemented in subclasses.
              (or if a generic algorithm)
        """
        raise NotImplementedError

    def are_equivalent(self, x, y, trans=False):
        r"""
        Test whether or not cusps x and y are equivalent modulo self.

        If self has a reduce_cusp() method, use that; otherwise do a
        slow explicit test.

        If trans = False, returns True or False. If trans = True, then
        return either False or an element of self mapping x onto y.

        EXAMPLE::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: a = NFCusp(K,0); b = NFCusp(K,1)
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).are_equivalent(a,b)
            Traceback (most recent call last):
            ...
            NotImplementedError
            """
        if hasattr(self, 'reduce_cusp'):
            if self.reduce_cusp(x) == self.reduce_cusp(y):
                return True
            return False
        else:
            return NotImplementedError

    def cusp_data(self, c):
        r"""
        Return a triple (g, w, t) where g is an element of self
        generating the stabiliser of the given cusp, w is the width of
        the cusp, and t is 1 if the cusp is regular and -1 if not.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: c = NFCusps(K,0)
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).cusp_data(c)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def index(self, G=None):
        r"""
        Return the index of self in G (default SL(2,O))

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: H = sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O)
            sage: H.index()
            Traceback (most recent call last):
            ...
            NotImplementedError: All subclasses should implement coset representatives
        """
        return len(list(self.coset_reps()))

    # Question: Does any of these algorithms apply for number fields.
    #
    # def is_regular_cusp(self, c):
    #     r"""
    #     Return True if the orbit of the given cusp is a regular cusp for self,
    #     otherwise False. This is automatically true if -1 is in self.

    #     EXAMPLES::

    #         sage: Gamma1(4).is_regular_cusp(Cusps(1/2))
    #         False
    #         sage: Gamma1(4).is_regular_cusp(Cusps(oo))
    #         True
    #     """
    #     if self.is_even(): return True
    #     return (self.cusp_data(c)[2] == 1)
    # def cusp_width(self, c):
    #     r"""
    #     Return the width of the orbit of cusps represented by c.

    #     EXAMPLES::

    #         sage: Gamma0(11).cusp_width(Cusps(oo))
    #         1
    #         sage: Gamma0(11).cusp_width(0)
    #         11
    #         sage: [Gamma0(100).cusp_width(c) for c in Gamma0(100).cusps()]
    #         [100, 1, 4, 1, 1, 1, 4, 25, 1, 1, 4, 1, 25, 4, 1, 4, 1, 1]
    #     """
    #     return self.cusp_data(c)[1]


    # def generalised_level(self):
    #     r"""
    #     Return the generalised level of self, i.e. the least common multiple of
    #     the widths of all cusps.

    #     If self is *even*, Wohlfart's theorem tells us that this is equal to
    #     the (conventional) level of self when self is a congruence subgroup.
    #     This can fail if self is odd, but the actual level is at most twice the
    #     generalised level. See the paper by Kiming, Schuett and Verrill for
    #     more examples.

    #     EXAMPLE::

    #         sage: Gamma0(18).generalised_level()
    #         18
    #         sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).index()
    #         Traceback (most recent call last):
    #         ...
    #         NotImplementedError

    #     In the following example, the actual level is twice the generalised
    #     level. This is the group `G_2` from Example 17 of K-S-V.

    #     ::

    #         sage: G = CongruenceSubgroup(8, [ [1,1,0,1], [3,-1,4,-1] ])
    #         sage: G.level()
    #         8
    #         sage: G.generalised_level()
    #         4
    #     """
    #     return arith.lcm([self.cusp_width(c) for c in self.cusps()])

    # def projective_index(self):
    #     r"""
    #     Return the index of the image of self in `{\rm PSL}_2(\ZZ)`. This is equal
    #     to the index of self if self contains -1, and half of this otherwise.

    #     This is equal to the degree of the natural map from the modular curve
    #     of self to the `j`-line.

    #     EXAMPLE::

    #         sage: Gamma0(5).projective_index()
    #         6
    #         sage: Gamma1(5).projective_index()
    #         12
    #     """

    #     if self.is_even():
    #         return self.index()
    #     else:
    #         return self.index() // 2

    def is_congruence(self):
        r"""
        Return True if self is a congruence subgroup.

        EXAMPLE::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: ArithmeticSubgroup_NF_class(O).is_congruence()
            False
        """
        return False

    @cached_method
    def generators(self):
        r"""
        Return a list of generators for this arithmetic subgroup. The result is cached.

        EXAMPLE::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).generators()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def gens(self, *args, **kwds):
        r"""
        Return a tuple of generators for this arithmetic subgroup.

        The generators need not be minimal. For arguments, see
        :meth:`~generators`.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).gens()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return tuple(self.generators(*args, **kwds))

    def gen(self, i):
        r"""
        Return the i-th generator of self, i.e. the i-th element of
        the tuple self.gens().

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).gen(1)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return self.generators()[i]

    def ngens(self):
        r"""
        Return the size of the minimal generating set of self returned by
        :meth:`generators`.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).ngens()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return len(self.generators())

    def ncusps(self):
        r"""
        Return the number of cusps of this arithmetic subgroup.

        This is provided as a separate function since for dimension
        formulae in even weight all we need to know is the number of
        cusps, and this can be calculated very quickly, while
        enumerating all cusps is much slower.

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).ncusps()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return ZZ(len(self.cusps()))

    def dimension_modular_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k modular forms for this
        group.

        EXAMPLE::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).dimension_modular_forms(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def dimension_cusp_forms(self, k=2):
        r"""
        Return the dimension of the space of weight k cusp forms for this
        group.

        EXAMPLE::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).dimension_cusp_forms(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def dimension_eis(self, weight=(2,)):
        r"""
        Return the dimension of the space of Eisenstein series of the
        given weight for this group, which is a subspace of the space
        of modular forms complementary to the space of cusp forms.

        INPUT:

        - ``k`` - an integer (default 2).

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: sage.rings.number_field.arithgroup_nf.all.ArithmeticSubgroup_NF_class(O).dimension_eis(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def sturm_bound(self, weight=(2,)):
        r"""
        Returns the Sturm bound for modular forms of the given weight and level
        this subgroup.

        INPUT:

        -  ``weight`` - an tuple of integers `\geq 2` (default: 2)

        EXAMPLES::

            sage: from sage.rings.number_field.arithgroup_nf.all import *
            sage: K = QuadraticField(41)
            sage: O = K.ring_of_integers()
            sage: ArithmeticSubgroup_NF_class(O).sturm_bound(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError
