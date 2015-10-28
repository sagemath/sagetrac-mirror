r"""
Congruence arithmetic subgroups of `{\rm SL}_2(\ZZ)`

Sage can compute extensively with the standard congruence subgroups
`\Gamma_0(N)`, `\Gamma_1(N)`, and `\Gamma_H(N)`.

AUTHORS:

- William Stein
- David Loeffler (2009, 10) -- modifications to work with more general arithmetic subgroups
"""

################################################################################
#
#       Copyright (C) 2004, 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#
################################################################################

from sage.rings.all import QQ, ZZ, Zmod
from sage.rings.arith import gcd
from sage.sets.set import Set
from sage.groups.matrix_gps.all import MatrixGroup
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.misc_c import prod
from arithgroup_generic import ArithmeticSubgroup


def CongruenceSubgroup_constructor(*args):
    r"""
    Attempt to create a congruence subgroup from the given data.

    The allowed inputs are as follows:

    - A :class:`~sage.groups.matrix_gps.matrix_group.MatrixGroup` object. This
      must be a group of matrices over `\ZZ / N\ZZ` for some `N`, with
      determinant 1, in which case the function will return the group of
      matrices in `SL(2, \ZZ)` whose reduction mod `N` is in the given group.

    - A list of matrices over `\ZZ / N\ZZ` for some `N`. The function will then
      compute the subgroup of `SL(2, \ZZ)` generated by these matrices, and
      proceed as above.

    - An integer `N` and a list of matrices (over any ring coercible to `\ZZ /
      N\ZZ`, e.g. over `\ZZ`). The matrices will then be coerced to `\ZZ /
      N\ZZ`.

    The function checks that the input G is valid. It then tests to see if
    `G` is the preimage mod `N` of some group of matrices modulo a proper
    divisor `M` of `N`, in which case it replaces `G` with this group before
    continuing.

    EXAMPLES::

        sage: from sage.modular.arithgroup.congroup_generic import CongruenceSubgroup_constructor as CS
        sage: CS(2, [[1,1,0,1]])
        Congruence subgroup of SL(2,Z) of level 2, preimage of:
         Matrix group over Ring of integers modulo 2 with 1 generators (
        [1 1]
        [0 1]
        )
        sage: CS([matrix(Zmod(2), 2, [1,1,0,1])])
        Congruence subgroup of SL(2,Z) of level 2, preimage of:
         Matrix group over Ring of integers modulo 2 with 1 generators (
        [1 1]
        [0 1]
        )
        sage: CS(MatrixGroup([matrix(Zmod(2), 2, [1,1,0,1])]))
        Congruence subgroup of SL(2,Z) of level 2, preimage of:
         Matrix group over Ring of integers modulo 2 with 1 generators (
        [1 1]
        [0 1]
        )
        sage: CS(SL(2, 2))
        Modular Group SL(2,Z)

    Some invalid inputs::

        sage: CS(SU(2, 7))
        Traceback (most recent call last):
        ...
        TypeError: Ring of definition must be Z / NZ for some N
    """
    from sage.groups.matrix_gps.matrix_group import is_MatrixGroup
    if is_MatrixGroup(args[0]):
        G = args[0]

    elif isinstance(args[0], list):
        G = MatrixGroup(args[0])

    elif args[0] in ZZ:
        M = MatrixSpace(Zmod(args[0]), 2)
        G = MatrixGroup([M(x) for x in args[1]])

    R = G.matrix_space().base_ring()
    if not hasattr(R, "cover_ring") or R.cover_ring() != ZZ:
        raise TypeError("Ring of definition must be Z / NZ for some N")

    if not all([x.matrix().det() == 1 for x in G.gens()]):
        raise ValueError("Group must be contained in SL(2, Z / N)")
    GG = _minimize_level(G)
    if GG in ZZ:
        from all import Gamma
        return Gamma(GG)
    else:
        return CongruenceSubgroupFromGroup(GG)

def is_CongruenceSubgroup(x):
    """
    Return True if x is of type CongruenceSubgroup.

    Note that this may be False even if `x` really is a congruence subgroup --
    it tests whether `x` is "obviously" congruence, i.e.~whether it has a
    congruence subgroup datatype. To test whether or not an arithmetic subgroup
    of `SL(2, \ZZ)` is congruence, use the ``is_congruence()`` method instead.

    EXAMPLES::

        sage: from sage.modular.arithgroup.congroup_generic import is_CongruenceSubgroup
        sage: is_CongruenceSubgroup(SL2Z)
        True
        sage: is_CongruenceSubgroup(Gamma0(13))
        True
        sage: is_CongruenceSubgroup(Gamma1(6))
        True
        sage: is_CongruenceSubgroup(GammaH(11, [3]))
        True
        sage: G = ArithmeticSubgroup_Permutation(L = "(1, 2)", R = "(1, 2)"); is_CongruenceSubgroup(G)
        False
        sage: G.is_congruence()
        True
        sage: is_CongruenceSubgroup(SymmetricGroup(3))
        False
    """
    return isinstance(x, CongruenceSubgroupBase)

class CongruenceSubgroupBase(ArithmeticSubgroup):

    def __init__(self, level):
        """
        Create a congruence subgroup with given level.

        EXAMPLES::

            sage: Gamma0(500)
            Congruence Subgroup Gamma0(500)
        """
        level = ZZ(level)
        if level <= 0:
            raise ArithmeticError("Congruence groups only defined for positive levels.")
        self.__level = level
        ArithmeticSubgroup.__init__(self)

    def _an_element_(self):
        r"""
        Return an element of self (mainly for use by the test suite).

        EXAMPLE::

            sage: Gamma(3).an_element() # indirect doctest
            [-2 -3]
            [ 3  4]
        """
        N = self.level()
        return self([1-N, -N, N, 1+N])

    def is_congruence(self):
        r"""
        Return True, since this is a congruence subgroup.

        EXAMPLE::

            sage: Gamma0(7).is_congruence()
            True
        """

        return True

    def level(self):
        """
        Return the level of this congruence subgroup.

        EXAMPLES::

            sage: SL2Z.level()
            1
            sage: Gamma0(20).level()
            20
            sage: Gamma1(11).level()
            11
            sage: GammaH(14, [5]).level()
            14
        """
        return self.__level

    def __cmp__(self, other):
        r"""
        EXAMPLE::

            sage: CongruenceSubgroup(3,[ [1,1,0,1] ]) == Gamma1(3)
            True
            sage: CongruenceSubgroup(3,[ [1,1,0,1] ]) == Gamma(3)
            False
            sage: CongruenceSubgroup(3,[ [1,1,0,1] ]) == QQ
            False
        """
        # This is carefully laid out so it can be called early on in the Sage
        # startup process when we want to create the standard generators of
        # SL2Z for use in arithgroup_perm. Hence it must work in this case
        # without being able to import the arithgroup_perm module. That's why
        # the most general case is *first*, not last.
        # Note that lazy_import doesn't work here, because it doesn't play
        # nicely with isinstance().
        if not isinstance(other, ArithmeticSubgroup):
            return cmp(type(self), type(other))

        elif is_CongruenceSubgroup(other):
            t = cmp(self.level(), other.level())
            if t: return t
            if self.level() == 1: return 0 # shouldn't come up except with pickling/unpickling
            t = cmp(self.index(), other.index())
            if t: return t
            return cmp(self.image_mod_n(),other.image_mod_n())

        from sage.modular.arithgroup.arithgroup_perm import ArithmeticSubgroup_Permutation_class
        if isinstance(other, ArithmeticSubgroup_Permutation_class):
            return cmp(self.as_permutation_group(), other)

        else:
            # we shouldn't ever get here
            raise NotImplementedError

class CongruenceSubgroupFromGroup(CongruenceSubgroupBase):
    r"""
    A congruence subgroup, defined by the data of an integer `N` and a subgroup
    `G` of the finite group `SL(2, \ZZ / N\ZZ)`; the congruence subgroup
    consists of all the matrices in `SL(2, \ZZ)` whose reduction modulo `N`
    lies in `G`.

    This class should not be instantiated directly, but created using the
    factory function
    :func:`~sage.modular.arithgroup.congroup_generic.CongruenceSubgroup_constructor`,
    which accepts much more flexible input, and checks the input to make sure
    it is valid.

    TESTS::

        sage: G = CongruenceSubgroup(5, [[0,-1,1,0]]); G
        Congruence subgroup of SL(2,Z) of level 5, preimage of:
         Matrix group over Ring of integers modulo 5 with 1 generators (
        [0 4]
        [1 0]
        )
        sage: TestSuite(G).run()
    """

    def __init__(self, G):
        r"""
        Standard init function.

        TESTS::

            sage: from sage.modular.arithgroup.congroup_generic import CongruenceSubgroupFromGroup
            sage: G = MatrixGroup([matrix(Zmod(2), 2, [1,1,1,0])])
            sage: CongruenceSubgroupFromGroup(G).index() # indirect doctest
            2
        """
        N = G.base_ring().characteristic()
        self.__G = G
        CongruenceSubgroupBase.__init__(self, N)

    def __reduce__(self):
        r"""
        Data defining self (for pickling).

        EXAMPLE::

            sage: G = CongruenceSubgroup(5, [[0,-1,1,0]])
            sage: G.__reduce__()
            (<function CongruenceSubgroup_constructor at ...>,
             (Matrix group over Ring of integers modulo 5 with 1 generators (
             [0 4]
             [1 0]
             ),))
        """
        return CongruenceSubgroup_constructor, (self.image_mod_n(),)

    def _contains_sl2(self, a,b,c,d):
        r"""
        Test whether ``[a,b;c,d]`` is an element of self.

        EXAMPLE::

            sage: G = MatrixGroup([matrix(Zmod(2), 2, [1,1,1,0])])
            sage: H = sage.modular.arithgroup.congroup_generic.CongruenceSubgroupFromGroup(G)
            sage: H(1)
            [1 0]
            [0 1]
            sage: H([0,-1,1,0])
            Traceback (most recent call last):
            ...
            TypeError: matrix [ 0 -1]
            [ 1  0] is not an element of Congruence subgroup of SL(2,Z) of level 2, preimage of:
             Matrix group over Ring of integers modulo 2 with 1 generators (
            [1 1]
            [1 0]
            )
            sage: H([1,2,0,1])
            [1 2]
            [0 1]
            sage: H(SL2Z([0,-1,1,0]), check=False)
            [ 0 -1]
            [ 1  0]
            sage: H([1,2,0,1]).parent()
            Modular Group SL(2,Z)
        """
        return ([a,b,c,d] in self.image_mod_n())

    def to_even_subgroup(self):
        r"""
        Return the smallest even subgroup of `SL(2, \ZZ)` containing self.

        EXAMPLE::

            sage: G = Gamma(3)
            sage: G.to_even_subgroup()
            Congruence subgroup of SL(2,Z) of level 3, preimage of:
             Matrix group over Ring of integers modulo 3 with 1 generators (
            [2 0]
            [0 2]
            )
        """
        if self.is_even():
            return self
        else:
            G = self.image_mod_n()
            H = MatrixGroup([ g.matrix() for g in G.gens()] + [G.matrix_space()(-1)])
            return CongruenceSubgroup_constructor(H)

    def _repr_(self):
        r"""
        String representation of self.

        EXAMPLE::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroupFromGroup(MatrixGroup([matrix(Zmod(2), 2, [1,1,1,0])]))._repr_()
            'Congruence subgroup of SL(2,Z) of level 2, preimage of:\n Matrix group over Ring of integers modulo 2 with 1 generators (\n[1 1]\n[1 0]\n)'
        """
        return "Congruence subgroup of SL(2,Z) of level %s, preimage of:\n %s" % (self.level(), self.image_mod_n())

    def index(self):
        r"""
        Return the index of self in the full modular group. This is equal to
        the index in `SL(2, \ZZ / N\ZZ)` of the image of this group modulo
        `\Gamma(N)`.

        EXAMPLE::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroupFromGroup(MatrixGroup([matrix(Zmod(2), 2, [1,1,1,0])])).index()
            2
        """
        return prod([p**(3*e-2)*(p*p-1) for (p,e) in self.level().factor()]) // self.image_mod_n().order()

    def image_mod_n(self):
        r"""
        Return the subgroup of `SL(2, \ZZ / N\ZZ)` of which this is the preimage, where `N` is the level of self.

        EXAMPLE::

            sage: G = MatrixGroup([matrix(Zmod(2), 2, [1,1,1,0])])
            sage: H = sage.modular.arithgroup.congroup_generic.CongruenceSubgroupFromGroup(G); H.image_mod_n()
            Matrix group over Ring of integers modulo 2 with 1 generators (
            [1 1]
            [1 0]
            )
            sage: H.image_mod_n() == G
            True
        """
        return self.__G

class CongruenceSubgroup(CongruenceSubgroupFromGroup):
    r"""
    One of the "standard" congruence subgroups `\Gamma_0(N)`, `\Gamma_1(N)`,
    `\Gamma(N)`, or `\Gamma_H(N)` (for some `H`).

    This class is not intended to be instantiated directly. Derived subclasses
    must override ``_contains_sl2``, ``_repr_``, and ``image_mod_n``.
    """

    def image_mod_n(self):
        r"""
        Raise an error: all derived subclasses should override this function.

        EXAMPLE::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5).image_mod_n()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __init__(self,*args, **kwds):
        r"""
        Bypass the init function of the CongruenceSubgroupFromGroup class.

        EXAMPLE::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5) # indirect doctest
            Generic congruence subgroup of level 5
        """
        return CongruenceSubgroupBase.__init__(self, *args, **kwds)

    def _repr_(self):
        """
        Return the string representation of self.

        NOTE: This function should be overridden by all subclasses.

        EXAMPLES::

            sage: sage.modular.arithgroup.congroup_generic.CongruenceSubgroup(5)._repr_()
            'Generic congruence subgroup of level 5'
        """
        return "Generic congruence subgroup of level %s" % self.level()

    def modular_symbols(self, sign=0, weight=2, base_ring=QQ):
        """
        Return the space of modular symbols of the specified weight and sign
        on the congruence subgroup self.

        EXAMPLES::

            sage: G = Gamma0(23)
            sage: G.modular_symbols()
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Rational Field
            sage: G.modular_symbols(weight=4)
            Modular Symbols space of dimension 12 for Gamma_0(23) of weight 4 with sign 0 over Rational Field
            sage: G.modular_symbols(base_ring=GF(7))
            Modular Symbols space of dimension 5 for Gamma_0(23) of weight 2 with sign 0 over Finite Field of size 7
            sage: G.modular_symbols(sign=1)
            Modular Symbols space of dimension 3 for Gamma_0(23) of weight 2 with sign 1 over Rational Field
        """
        from sage.modular.modsym.modsym import ModularSymbols
        return ModularSymbols(self, sign=sign, weight=weight, base_ring=base_ring)

    def modular_abelian_variety(self):
        """
        Return the modular abelian variety corresponding to the congruence
        subgroup self.

        EXAMPLES::

            sage: Gamma0(11).modular_abelian_variety()
            Abelian variety J0(11) of dimension 1
            sage: Gamma1(11).modular_abelian_variety()
            Abelian variety J1(11) of dimension 1
            sage: GammaH(11,[3]).modular_abelian_variety()
            Abelian variety JH(11,[3]) of dimension 1
        """
        from sage.modular.abvar.abvar_ambient_jacobian import ModAbVar_ambient_jacobian
        return ModAbVar_ambient_jacobian(self)

    def _new_group_from_level(self, level):
        r"""
        Return a new group of the same type (Gamma0, Gamma1, or
        GammaH) as self of the given level. In the case that self is of type
        GammaH, we take the largest H inside `(\ZZ/ \text{level}\ZZ)^\times`
        which maps to H, namely its inverse image under the natural reduction
        map.

        EXAMPLES::

            sage: G = Gamma0(20)
            sage: G._new_group_from_level(4)
            Congruence Subgroup Gamma0(4)
            sage: G._new_group_from_level(40)
            Congruence Subgroup Gamma0(40)

            sage: G = Gamma1(10)
            sage: G._new_group_from_level(6)
            Traceback (most recent call last):
            ...
            ValueError: one level must divide the other

            sage: G = GammaH(50,[7]); G
            Congruence Subgroup Gamma_H(50) with H generated by [7]
            sage: G._new_group_from_level(25)
            Congruence Subgroup Gamma_H(25) with H generated by [7]
            sage: G._new_group_from_level(10)
            Congruence Subgroup Gamma0(10)
            sage: G._new_group_from_level(100)
            Congruence Subgroup Gamma_H(100) with H generated by [7, 57]
        """
        from congroup_gamma0 import is_Gamma0
        from congroup_gamma1 import is_Gamma1
        from congroup_gammaH import is_GammaH
        from all import Gamma0, Gamma1, GammaH
        N = self.level()
        if (level%N) and (N%level):
            raise ValueError("one level must divide the other")
        if is_Gamma0(self):
            return Gamma0(level)
        elif is_Gamma1(self):
            return Gamma1(level)
        elif is_GammaH(self):
            H = self._generators_for_H()
            if level > N:
                d = level // N
                diffs = [ N*i for i in range(d) ]
                newH = [ h + diff for h in H for diff in diffs ]
                return GammaH(level, [x for x in newH if gcd(level, x) == 1])
            else:
                return GammaH(level, [ h%level for h in H ])
        else:
            raise NotImplementedError

def _minimize_level(G):
    r"""
    Utility function. Given a matrix group `G` contained in `SL(2, \ZZ / N\ZZ)`
    for some `N`, test whether or not `G` is the preimage of a subgroup of
    smaller level, and if so, return that subgroup.

    The trivial group is handled specially: instead of returning a group, it
    returns an integer `N`, representing the trivial subgroup of `SL(2, \ZZ /
    N\ZZ)`.

    EXAMPLE::

        sage: M = MatrixSpace(Zmod(9), 2, 2)
        sage: G = MatrixGroup([M(x) for x in [[1,1,0,1],[1,3,0,1],[1,0,3,1],[4,0,0,7]]]); G
        Matrix group over Ring of integers modulo 9 with 4 generators (
        [1 1]  [1 3]  [1 0]  [4 0]
        [0 1], [0 1], [3 1], [0 7]
        )
        sage: sage.modular.arithgroup.congroup_generic._minimize_level(G)
        Matrix group over Ring of integers modulo 3 with 1 generators (
        [1 1]
        [0 1]
        )
        sage: G = MatrixGroup([M(x) for x in [[1,3,0,1],[1,0,3,1],[4,0,0,7]]]); G
        Matrix group over Ring of integers modulo 9 with 3 generators (
        [1 3]  [1 0]  [4 0]
        [0 1], [3 1], [0 7]
        )
        sage: sage.modular.arithgroup.congroup_generic._minimize_level(G)
        3
    """
    from congroup_gamma import Gamma_constructor as Gamma
    Glist = list(G)
    N = G.base_ring().characteristic()
    i = Gamma(N).index()

    for d in N.divisors()[:-1]:
        j = Gamma(d).index()
        k = len([g for g in Glist if g.matrix().change_ring(Zmod(d)) == 1])
        if k == i // j:
            if d == 1:
                return ZZ(1)
            G = MatrixGroup([g.matrix().change_ring(Zmod(d)) for g in G.gens()])
            N = d
            break

    # now sanitize the generators (remove duplicates and copies of the identity)
    new_gens = [x.matrix() for x in G.gens() if x.matrix() != 1]
    all([x.set_immutable() for x in new_gens])
    new_gens = list(Set(new_gens))
    if new_gens == []:
        return ZZ(N)
    return MatrixGroup(new_gens)
