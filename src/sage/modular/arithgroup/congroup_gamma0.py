r"""
Congruence Subgroup `\Gamma_0(N)`
"""
from __future__ import absolute_import

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from six.moves import range

from .congroup_gammaH import GammaH_class
from .congroup_gamma1 import is_Gamma1
from sage.modular.modsym.p1list import lift_to_sl2z
from .congroup_generic import CongruenceSubgroup

from sage.modular.cusps import Cusp
from sage.misc.cachefunc import cached_method
from sage.rings.all import IntegerModRing, ZZ
from sage.arith.all import kronecker_symbol
from sage.misc.all import prod
import sage.modular.modsym.p1list
import sage.arith.all as arith


def is_Gamma0(x):
    """
    Return True if x is a congruence subgroup of type Gamma0.

    EXAMPLES::

        sage: from sage.modular.arithgroup.all import is_Gamma0
        sage: is_Gamma0(SL2Z)
        True
        sage: is_Gamma0(Gamma0(13))
        True
        sage: is_Gamma0(Gamma1(6))
        False
    """
    return isinstance(x, Gamma0_class)

_gamma0_cache = {}
def Gamma0_constructor(N):
    """
    Return the congruence subgroup Gamma0(N).

    EXAMPLES::

        sage: G = Gamma0(51) ; G # indirect doctest
        Congruence Subgroup Gamma0(51)
        sage: G == Gamma0(51)
        True
        sage: G is Gamma0(51)
        True
    """
    from .all import SL2Z
    if N == 1: return SL2Z
    try:
        return _gamma0_cache[N]
    except KeyError:
        _gamma0_cache[N] = Gamma0_class(N)
        return _gamma0_cache[N]

class Gamma0_class(GammaH_class):
    r"""
    The congruence subgroup `\Gamma_0(N)`.

    TESTS::

        sage: Gamma0(11).dimension_cusp_forms(2)
        1
        sage: a = Gamma0(1).dimension_cusp_forms(2); a
        0
        sage: type(a)
        <type 'sage.rings.integer.Integer'>
        sage: Gamma0(5).dimension_cusp_forms(0)
        0
        sage: Gamma0(20).dimension_cusp_forms(1)
        0
        sage: Gamma0(20).dimension_cusp_forms(4)
        6

        sage: Gamma0(23).dimension_cusp_forms(2)
        2
        sage: Gamma0(1).dimension_cusp_forms(24)
        2
        sage: Gamma0(3).dimension_cusp_forms(3)
        0
        sage: Gamma0(11).dimension_cusp_forms(-1)
        0

        sage: Gamma0(22).dimension_new_cusp_forms()
        0
        sage: Gamma0(100).dimension_new_cusp_forms(2, 5)
        5

    Independently compute the dimension 5 above::

        sage: m = ModularSymbols(100, 2,sign=1).cuspidal_subspace()
        sage: m.new_subspace(5)
        Modular Symbols subspace of dimension 5 of Modular Symbols space of dimension 18 for Gamma_0(100) of weight 2 with sign 1 over Rational Field

    """


    def __init__(self, level):
        r"""
        The congruence subgroup `\Gamma_0(N)`.

        EXAMPLES::

            sage: G = Gamma0(11); G
            Congruence Subgroup Gamma0(11)
            sage: TestSuite(G).run()
            sage: G is loads(dumps(G))
            True

        TESTS::

            sage: g = Gamma0(5)([1,1,0,1])
            sage: g in Gamma0(7)
            True
            sage: g = Gamma0(5)([1,0,5,1])
            sage: g in Gamma0(7)
            False
            sage: g = Gamma0(2)([1,0,0,1])
            sage: g in SL2Z
            True
        """
        CongruenceSubgroup.__init__(self, level)

        # We *don't* call the GammaH init script, as this requires calculating
        # generators for the units modulo N which is time-consuming; this will
        # be done if needed by the _generators_for_H and _list_of_elements_in_H
        # methods.
        #
        #GammaH_class.__init__(self, level, [int(x) for x in IntegerModRing(level).unit_gens()])

    def _repr_(self):
        """
        Return the string representation of self.

        EXAMPLES::

            sage: Gamma0(98)._repr_()
            'Congruence Subgroup Gamma0(98)'
        """
        return "Congruence Subgroup Gamma0(%s)"%self.level()

    def __reduce__(self):
        """
        Used for pickling self.

        EXAMPLES::

            sage: Gamma0(22).__reduce__()
            (<function Gamma0_constructor at ...>, (22,))
        """
        return Gamma0_constructor, (self.level(),)

    def _latex_(self):
        r"""
        Return the \LaTeX representation of self.

        EXAMPLES::

            sage: Gamma0(20)._latex_()
            '\\Gamma_0(20)'
            sage: latex(Gamma0(20))
            \Gamma_0(20)
        """
        return "\\Gamma_0(%s)"%self.level()

    @cached_method
    def _generators_for_H(self):
        """
        Return generators for the subgroup H of the units mod
        self.level() that defines self.

        EXAMPLES::

            sage: Gamma0(15)._generators_for_H()
            [11, 7]
        """
        if self.level() in [1, 2]:
            return []
        return [ZZ(x) for x in IntegerModRing(self.level()).unit_gens()]

    @cached_method
    def _list_of_elements_in_H(self):
        """
        Returns a sorted list of Python ints that are representatives
        between 0 and N-1 of the elements of H.

        EXAMPLES::

            sage: G = Gamma0(11)
            sage: G._list_of_elements_in_H()
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

            sage: G = Gamma0(6)
            sage: G._list_of_elements_in_H()
            [1, 5]

            sage: G = Gamma0(1)
            sage: G._list_of_elements_in_H()
            [1]
        """
        N = self.level()
        if N != 1:
            gcd = arith.gcd
            H = [x for x in range(1, N) if gcd(x, N) == 1]
        else:
            H = [1]

        return H

    def divisor_subgroups(self):
        r"""
        Return the subgroups of SL2Z of the form Gamma0(M) that contain this subgroup,
        i.e. those for M a divisor of N.

        EXAMPLES::

            sage: Gamma0(24).divisor_subgroups()
            [Modular Group SL(2,Z),
            Congruence Subgroup Gamma0(2),
            Congruence Subgroup Gamma0(3),
            Congruence Subgroup Gamma0(4),
            Congruence Subgroup Gamma0(6),
            Congruence Subgroup Gamma0(8),
            Congruence Subgroup Gamma0(12),
            Congruence Subgroup Gamma0(24)]
        """
        return [Gamma0_constructor(M) for M in self.level().divisors()]

    def is_even(self):
        r"""
        Return True precisely if this subgroup contains the matrix -1.

        Since `\Gamma0(N)` always contains the matrix -1, this always
        returns True.

        EXAMPLES::

            sage: Gamma0(12).is_even()
            True
            sage: SL2Z.is_even()
            True
        """
        return True

    def is_subgroup(self, right):
        """
        Return True if self is a subgroup of right.

        EXAMPLES::

            sage: G = Gamma0(20)
            sage: G.is_subgroup(SL2Z)
            True
            sage: G.is_subgroup(Gamma0(4))
            True
            sage: G.is_subgroup(Gamma0(20))
            True
            sage: G.is_subgroup(Gamma0(7))
            False
            sage: G.is_subgroup(Gamma1(20))
            False
            sage: G.is_subgroup(GammaH(40, []))
            False
            sage: Gamma0(80).is_subgroup(GammaH(40, [31, 21, 17]))
            True
            sage: Gamma0(2).is_subgroup(Gamma1(2))
            True
        """
        if right.level() == 1:
            return True
        if is_Gamma0(right):
            return self.level() % right.level() == 0
        if is_Gamma1(right):
            if right.level() >= 3:
                return False
            elif right.level() == 2:
                return self.level() == 2
            # case level 1 dealt with above
        else:
            return GammaH_class.is_subgroup(self, right)

    def coset_reps(self):
        r"""
        Return representatives for the right cosets of this congruence
        subgroup in `{\rm SL}_2(\ZZ)` as a generator object.

        Use ``list(self.coset_reps())`` to obtain coset reps as a
        list.

        EXAMPLES::

            sage: list(Gamma0(5).coset_reps())
            [
            [1 0]  [ 0 -1]  [1 0]  [ 0 -1]  [ 0 -1]  [ 0 -1]
            [0 1], [ 1  0], [1 1], [ 1  2], [ 1  3], [ 1  4]
            ]
            sage: list(Gamma0(4).coset_reps())
            [
            [1 0]  [ 0 -1]  [1 0]  [ 0 -1]  [ 0 -1]  [1 0]
            [0 1], [ 1  0], [1 1], [ 1  2], [ 1  3], [2 1]
            ]
            sage: list(Gamma0(1).coset_reps())
            [
            [1 0]
            [0 1]
            ]
        """
        from .all import SL2Z
        N = self.level()
        if N == 1: # P1List isn't very happy working modulo 1
            yield SL2Z([1,0,0,1])
        else:
            for z in sage.modular.modsym.p1list.P1List(N):
                yield SL2Z(lift_to_sl2z(z[0], z[1], N))

    @cached_method
    def generators(self, algorithm="farey"):
        r"""
        Return generators for this congruence subgroup.

        INPUT:

        - ``algorithm`` (string): either ``farey`` (default) or
          ``todd-coxeter``.

        If ``algorithm`` is set to ``"farey"``, then the generators will be
        calculated using Farey symbols, which will always return a *minimal*
        generating set. See :mod:`~sage.modular.arithgroup.farey_symbol` for
        more information.

        If ``algorithm`` is set to ``"todd-coxeter"``, a simpler algorithm
        based on Todd-Coxeter enumeration will be used. This tends to return
        far larger sets of generators.

        EXAMPLE::

            sage: Gamma0(3).generators()
            [
            [1 1]  [-1  1]
            [0 1], [-3  2]
            ]
            sage: Gamma0(3).generators(algorithm="todd-coxeter")
            [
            [1 1]  [-1  0]  [ 1 -1]  [1 0]  [1 1]  [-1  0]  [ 1  0]
            [0 1], [ 0 -1], [ 0  1], [3 1], [0 1], [ 3 -1], [-3  1]
            ]
            sage: SL2Z.gens()
            (
            [ 0 -1]  [1 1]
            [ 1  0], [0 1]
            )
        """
        if self.level() == 1:
            # we return a fixed set of generators for SL2Z, for historical
            # reasons, which aren't the ones the Farey symbol code gives
            return [ self([0,-1,1,0]), self([1,1,0,1]) ]

        elif algorithm=="farey":
            return self.farey_symbol().generators()

        elif algorithm=="todd-coxeter":
            from sage.modular.modsym.p1list import P1List
            from .congroup import generators_helper
            level = self.level()
            if level == 1: # P1List isn't very happy working mod 1
                return [ self([0,-1,1,0]), self([1,1,0,1]) ]
            gen_list = generators_helper(P1List(level), level)
            return [self(g, check=False) for g in gen_list]

        else:
            raise ValueError("Unknown algorithm '%s' (should be either 'farey' or 'todd-coxeter')" % algorithm)

    def gamma_h_subgroups(self):
        r"""
        Return the subgroups of the form `\Gamma_H(N)` contained
        in self, where `N` is the level of self.

        EXAMPLES::

            sage: G = Gamma0(11)
            sage: G.gamma_h_subgroups()
            [Congruence Subgroup Gamma0(11), Congruence Subgroup Gamma_H(11) with H generated by [4], Congruence Subgroup Gamma_H(11) with H generated by [10], Congruence Subgroup Gamma1(11)]
            sage: G = Gamma0(12)
            sage: G.gamma_h_subgroups()
            [Congruence Subgroup Gamma0(12), Congruence Subgroup Gamma_H(12) with H generated by [7], Congruence Subgroup Gamma_H(12) with H generated by [11], Congruence Subgroup Gamma_H(12) with H generated by [5], Congruence Subgroup Gamma1(12)]
        """
        from .all import GammaH
        N = self.level()
        R = IntegerModRing(N)
        return [GammaH(N, H) for H in R.multiplicative_subgroups()]

    def _contains_sl2(self, a,b,c,d):
        r"""
        Test whether x is an element of this group.

        EXAMPLES::

            sage: G = Gamma0(12)
            sage: [1, 0, 24, 1] in G
            True
            sage: matrix(ZZ, 2, [1, 1, -12, -11]) in G
            True
            sage: SL2Z([0,-1,1,0]) in G
            False
            sage: 1 in G
            True
            sage: -1 in G
            True
            sage: 2 in G
            False

        The _element_constructor_ method calls this method::

            sage: G([1, 0, 23, 1]) # indirect doctest
            Traceback (most recent call last):
            ...
            TypeError: matrix [ 1  0]
            [23  1] is not an element of Congruence Subgroup Gamma0(12)
        """
        return (c % self.level() == 0)

    def _find_cusps(self):
        r"""
        Return an ordered list of inequivalent cusps for self, i.e. a
        set of representatives for the orbits of self on
        `\mathbb{P}^1(\QQ)`.  These are returned in a reduced
        form; see self.reduce_cusp for the definition of reduced.

        ALGORITHM:
            Uses explicit formulae specific to `\Gamma_0(N)`: a reduced cusp on
            `\Gamma_0(N)` is always of the form `a/d` where `d | N`, and `a_1/d
            \sim a_2/d` if and only if `a_1 \cong a_2 \bmod {\rm gcd}(d,
            N/d)`.

        EXAMPLES::

            sage: Gamma0(90)._find_cusps()
            [0, 1/45, 1/30, 1/18, 1/15, 1/10, 1/9, 2/15, 1/6, 1/5, 1/3, 11/30, 1/2, 2/3, 5/6, Infinity]
            sage: Gamma0(1).cusps()
            [Infinity]
            sage: Gamma0(180).cusps() == Gamma0(180).cusps(algorithm='modsym')
            True
        """
        N = self.level()
        s = []

        for d in arith.divisors(N):
            w = arith.gcd(d, N//d)
            if w == 1:
                if d == 1:
                    s.append(Cusp(1,0))
                elif d == N:
                    s.append(Cusp(0,1))
                else:
                    s.append(Cusp(1,d))
            else:
                for a in range(1, w):
                    if arith.gcd(a, w) == 1:
                        while arith.gcd(a, d//w) != 1:
                            a += w
                        s.append(Cusp(a,d))
        return sorted(s)

    def ncusps(self):
        r"""
        Return the number of cusps of this subgroup `\Gamma_0(N)`.

        EXAMPLES::

            sage: [Gamma0(n).ncusps() for n in [1..19]]
            [1, 2, 2, 3, 2, 4, 2, 4, 4, 4, 2, 6, 2, 4, 4, 6, 2, 8, 2]
            sage: [Gamma0(n).ncusps() for n in prime_range(2,100)]
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        """
        n = self.level()
        return sum([arith.euler_phi(arith.gcd(d,n//d)) for d in n.divisors()])


    def nu2(self):
        r"""
        Return the number of elliptic points of order 2 for this congruence
        subgroup `\Gamma_0(N)`. The number of these is given by a standard formula:
        0 if `N` is divisible by 4 or any prime congruent to -1 mod 4, and
        otherwise `2^d` where d is the number of odd primes dividing `N`.

        EXAMPLE::

            sage: Gamma0(2).nu2()
            1
            sage: Gamma0(4).nu2()
            0
            sage: Gamma0(21).nu2()
            0
            sage: Gamma0(1105).nu2()
            8
            sage: [Gamma0(n).nu2() for n in [1..19]]
            [1, 1, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 0, 2, 0, 0]
        """
        n = self.level()
        if n%4 == 0:
            return ZZ(0)
        return prod([ 1 + kronecker_symbol(-4, p) for p, _ in n.factor()])

    def nu3(self):
        r"""
        Return the number of elliptic points of order 3 for this congruence
        subgroup `\Gamma_0(N)`. The number of these is given by a standard formula:
        0 if `N` is divisible by 9 or any prime congruent to -1 mod 3, and
        otherwise `2^d` where d is the number of primes other than 3 dividing `N`.

        EXAMPLE::

            sage: Gamma0(2).nu3()
            0
            sage: Gamma0(3).nu3()
            1
            sage: Gamma0(9).nu3()
            0
            sage: Gamma0(7).nu3()
            2
            sage: Gamma0(21).nu3()
            2
            sage: Gamma0(1729).nu3()
            8
        """
        n = self.level()
        if (n % 9 == 0):
            return ZZ(0)
        return prod([ 1 + kronecker_symbol(-3, p) for p, _ in n.factor()])

    def index(self):
        r"""
        Return the index of self in the full modular group.

        This is given by

        .. MATH::

            N \prod_{\substack{p \mid N \\ \text{$p$ prime}}}\left(1 + \frac{1}{p}\right).

        EXAMPLES::

            sage: [Gamma0(n).index() for n in [1..19]]
            [1, 3, 4, 6, 6, 12, 8, 12, 12, 18, 12, 24, 14, 24, 24, 24, 18, 36, 20]
            sage: Gamma0(32041).index()
            32220
        """
        return prod([p**e + p**(e-1) for (p,e) in self.level().factor()])
