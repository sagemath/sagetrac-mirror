"""
Torsion subgroups of modular abelian varieties

Sage can compute information about the structure of the torsion
subgroup of a modular abelian variety. Sage computes a multiple of
the order by computing the greatest common divisor of the orders of
the torsion subgroup of the reduction of the abelian variety modulo
p for various primes p. Sage computes a divisor of the order by
computing the rational cuspidal subgroup. When these two bounds
agree (which is often the case), we determine the exact structure
of the torsion subgroup.

AUTHORS:

- William Stein (2007-03)

EXAMPLES: First we consider `J_0(50)` where everything
works out nicely::

    sage: J = J0(50)
    sage: T = J.rational_torsion_subgroup(); T
    Torsion subgroup of Abelian variety J0(50) of dimension 2
    sage: T.multiple_of_order()
    15
    sage: T.divisor_of_order()
    15
    sage: T.gens()
    [[(1/15, 3/5, 2/5, 14/15)]]
    sage: T.invariants()
    [15]
    sage: d = J.decomposition(); d
    [
    Simple abelian subvariety 50a(1,50) of dimension 1 of J0(50),
    Simple abelian subvariety 50b(1,50) of dimension 1 of J0(50)
    ]
    sage: d[0].rational_torsion_subgroup().order()
    3
    sage: d[1].rational_torsion_subgroup().order()
    5

Next we make a table of the upper and lower bounds for each new
factor.

::

    sage: for N in range(1,38):
    ...    for A in J0(N).new_subvariety().decomposition():
    ...        T = A.rational_torsion_subgroup()
    ...        print '%-5s%-5s%-5s%-5s'%(N, A.dimension(), T.divisor_of_order(), T.multiple_of_order())
    11   1    5    5
    14   1    6    6
    15   1    8    8
    17   1    4    4
    19   1    3    3
    20   1    6    6
    21   1    8    8
    23   2    11   11
    24   1    8    8
    26   1    3    3
    26   1    7    7
    27   1    3    3
    29   2    7    7
    30   1    6    12
    31   2    5    5
    32   1    4    4
    33   1    4    4
    34   1    6    6
    35   1    3    3
    35   2    16   16
    36   1    6    6
    37   1    1    1
    37   1    3    3

TESTS::

    sage: T = J0(54).rational_torsion_subgroup()
    sage: loads(dumps(T)) == T
    True
"""

#*****************************************************************************
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.modular.abvar.torsion_point import TorsionPoint
from sage.modules.module            import Module
from finite_subgroup                import FiniteSubgroup
from sage.rings.all                 import ZZ
from sage.sets.primes               import Primes
from sage.modular.arithgroup.all    import is_Gamma0
from sage.arith.all import divisors, gcd, prime_range

class RationalTorsionSubgroup(FiniteSubgroup):
    """
    The torsion subgroup of a modular abelian variety.
    """
    def __init__(self, abvar):
        """
        Create the torsion subgroup.

        INPUT:


        -  ``abvar`` - a modular abelian variety


        EXAMPLES::

            sage: T = J0(14).rational_torsion_subgroup(); T
            Torsion subgroup of Abelian variety J0(14) of dimension 1
            sage: type(T)
            <class 'sage.modular.abvar.torsion_subgroup.RationalTorsionSubgroup_with_category'>
        """
        FiniteSubgroup.__init__(self, abvar)

    def _repr_(self):
        """
        Return string representation of this torsion subgroup.

        EXAMPLES::

            sage: T = J1(13).rational_torsion_subgroup(); T
            Torsion subgroup of Abelian variety J1(13) of dimension 2
            sage: T._repr_()
            'Torsion subgroup of Abelian variety J1(13) of dimension 2'
        """
        return "Torsion subgroup of %s"%self.abelian_variety()

    def __cmp__(self, other):
        """
        Compare torsion subgroups.

        INPUT:


        -  ``other`` - an object


        If other is a torsion subgroup, the abelian varieties are compared.
        Otherwise, the generic behavior for finite abelian variety
        subgroups is used.

        EXAMPLE::

            sage: G = J0(11).rational_torsion_subgroup(); H = J0(13).rational_torsion_subgroup()
            sage: G == G
            True
            sage: G < H   # since 11 < 13
            True
            sage: G > H
            False
            sage: G < 5   # random (meaningless since it depends on memory layout)
            False
        """
        if isinstance(other, RationalTorsionSubgroup):
            return cmp(self.abelian_variety(), other.abelian_variety())
        return FiniteSubgroup.__cmp__(self, other)

    def order(self):
        """
        Return the order of the torsion subgroup of this modular abelian
        variety.

        This may fail if the multiple obtained by counting points modulo
        `p` exceeds the divisor obtained from the rational cuspidal
        subgroup.

        EXAMPLES::

            sage: a = J0(11)
            sage: a.rational_torsion_subgroup().order()
            5
            sage: a = J0(23)
            sage: a.rational_torsion_subgroup().order()
            11
            sage: t = J0(37)[1].rational_torsion_subgroup()
            sage: t.order()
            3
        """
        try:
            return self._order
        except AttributeError:
            pass
        O = self.possible_orders()
        if len(O) == 1:
            n = O[0]
            self._order = n
            return n
        raise RuntimeError("Unable to compute order of torsion subgroup (it is in %s)"%O)

    def lattice(self):
        """
        Return lattice that defines this torsion subgroup, if possible.

        .. warning::

           There is no known algorithm in general to compute the
           rational torsion subgroup. Use rational_cusp_group to
           obtain a subgroup of the rational torsion subgroup in
           general.

        EXAMPLES::

            sage: J0(11).rational_torsion_subgroup().lattice()
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [  1   0]
            [  0 1/5]

        The following fails because in fact I know of no (reasonable)
        algorithm to provably compute the torsion subgroup in general.

        ::

            sage: T = J0(33).rational_torsion_subgroup()
            sage: T.lattice()
            Traceback (most recent call last):
            ...
            NotImplementedError: unable to compute the rational torsion subgroup in this case (there is no known general algorithm yet)

        The problem is that the multiple of the order obtained by counting
        points over finite fields is twice the divisor of the order got
        from the rational cuspidal subgroup.

        ::

            sage: T.multiple_of_order(30)
            200
            sage: J0(33).rational_cusp_subgroup().order()
            100
        """
        A = self.abelian_variety()
        if A.dimension() == 0:
            return []
        R = A.rational_cusp_subgroup()
        if R.order() == self.multiple_of_order():
            return R.lattice()
        else:
            raise NotImplementedError("unable to compute the rational torsion subgroup in this case (there is no known general algorithm yet)")

    def possible_orders(self):
        """
        Return the possible orders of this torsion subgroup, computed from
        a known divisor and multiple of the order.

        EXAMPLES::

            sage: J0(11).rational_torsion_subgroup().possible_orders()
            [5]
            sage: J0(33).rational_torsion_subgroup().possible_orders()
            [100, 200]

        Note that this function has not been implemented for `J_1(N)`,
        though it should be reasonably easy to do so soon (see Conrad,
        Edixhoven, Stein)::

            sage: J1(13).rational_torsion_subgroup().possible_orders()
            Traceback (most recent call last):
            ...
            NotImplementedError: torsion multiple only implemented for Gamma0
        """
        try:
            return self._possible_orders
        except AttributeError:
            pass
        u = self.multiple_of_order()
        l = self.divisor_of_order()
        assert u % l == 0
        O = [l * d for d in divisors(u//l)]
        self._possible_orders = O
        return O

    def divisor_of_order(self):
        """
        Return a divisor of the order of this torsion subgroup of a modular
        abelian variety.

        EXAMPLES::

            sage: t = J0(37)[1].rational_torsion_subgroup()
            sage: t.divisor_of_order()
            3
        """
        A = self.abelian_variety()
        if A.dimension() == 0:
            return ZZ(1)
        R = A.rational_cusp_subgroup()
        return R.order()

    def multiple_of_order(self, maxp=None):
        """
        Return a multiple of the order of this torsion group.

        The multiple is computed using characteristic polynomials of Hecke
        operators of odd index not dividing the level.

        INPUT:


        -  ``maxp`` - (default: None) If maxp is None (the
           default), return gcd of best bound computed so far with bound
           obtained by computing GCD's of orders modulo p until this gcd
           stabilizes for 3 successive primes. If maxp is given, just use all
           primes up to and including maxp.


        EXAMPLES::

            sage: J = J0(11)
            sage: G = J.rational_torsion_subgroup()
            sage: G.multiple_of_order(11)
            5
            sage: J = J0(389)
            sage: G = J.rational_torsion_subgroup(); G
            Torsion subgroup of Abelian variety J0(389) of dimension 32
            sage: G.multiple_of_order()
            97
            sage: [G.multiple_of_order(p) for p in prime_range(3,11)]
            [92645296242160800, 7275, 291]
            sage: [G.multiple_of_order(p) for p in prime_range(3,13)]
            [92645296242160800, 7275, 291, 97]
            sage: [G.multiple_of_order(p) for p in prime_range(3,19)]
            [92645296242160800, 7275, 291, 97, 97, 97]

        ::

            sage: J = J0(33) * J0(11) ; J.rational_torsion_subgroup().order()
            Traceback (most recent call last):
            ...
            NotImplementedError: torsion multiple only implemented for Gamma0

        The next example illustrates calling this function with a larger
        input and how the result may be cached when maxp is None::

            sage: T = J0(43)[1].rational_torsion_subgroup()
            sage: T.multiple_of_order()
            14
            sage: T.multiple_of_order(50)
            7
            sage: T.multiple_of_order()
            7
        """
        if maxp is None:
            try:
                return self.__multiple_of_order
            except AttributeError:
                pass
        bnd = ZZ(0)
        A = self.abelian_variety()
        if A.dimension() == 0:
            T = ZZ(1)
            self.__multiple_of_order = T
            return T
        N = A.level()
        if not (len(A.groups()) == 1 and is_Gamma0(A.groups()[0])):
            # to generalize to this case, you'll need to
            # (1) define a charpoly_of_frob function:
            #       this is tricky because I don't know a simple
            #       way to do this for Gamma1 and GammaH.  Will
            #       probably have to compute explicit matrix for
            #       <p> operator (add to modular symbols code),
            #       then compute some charpoly involving
            #       that directly...
            # (2) use (1) -- see my MAGMA code.
            raise NotImplementedError("torsion multiple only implemented for Gamma0")
        cnt = 0
        if maxp is None:
            X = Primes()
        else:
            X = prime_range(maxp+1)
        for p in X:
            if (2*N) % p == 0:
                continue

            f = A.hecke_polynomial(p)
            b = ZZ(f(p+1))

            if bnd == 0:
                bnd = b
            else:
                bnd_last = bnd
                bnd = ZZ(gcd(bnd, b))
                if bnd == bnd_last:
                    cnt += 1
                else:
                    cnt = 0
                if maxp is None and cnt >= 2:
                    break

        # The code below caches the computed bound and
        # will be used if this function is called
        # again with maxp equal to None (the default).
        if maxp is None:
            # maxp is None but self.__multiple_of_order  is
            # not set, since otherwise we would have immediately
            # returned at the top of this function
            self.__multiple_of_order = bnd
        else:
            # maxp is given -- record new info we get as
            # a gcd...
            try:
                self.__multiple_of_order = gcd(self.__multiple_of_order, bnd)
            except AttributeError:
                # ... except in the case when self.__multiple_of_order
                # was never set.  In this case, we just set
                # it as long as the gcd stabilized for 3 in a row.
                if cnt >= 2:
                    self.__multiple_of_order = bnd
        return bnd


class QQbarTorsionSubgroup(Module):

    Element = TorsionPoint

    def __init__(self, abvar):
        """
        Group of all torsion points over the algebraic closure on an
        abelian variety.

        INPUT:


        -  ``abvar`` - an abelian variety


        EXAMPLES::

            sage: A = J0(23)
            sage: A.qbar_torsion_subgroup()
            Group of all torsion points in QQbar on Abelian variety J0(23) of dimension 2
        """
        self.__abvar = abvar
        Module.__init__(self, ZZ)

    def _repr_(self):
        """
        Print representation of QQbar points.

        OUTPUT: string

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup()._repr_()
            'Group of all torsion points in QQbar on Abelian variety J0(23) of dimension 2'
        """
        return 'Group of all torsion points in QQbar on %s'%self.__abvar

    def field_of_definition(self):
        """
        Return the field of definition of this subgroup. Since this is the
        group of all torsion it is defined over the base field of this
        abelian variety.

        OUTPUT: a field

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup().field_of_definition()
            Rational Field
        """
        return self.__abvar.base_field()

    def _element_constructor_(self, x):
        r"""
        Create an element in this torsion subgroup.

        INPUT:

        - ``x`` -- vector in `\QQ^{2d}`

        OUTPUT: torsion point

        EXAMPLES::

            sage: P = J0(23).qbar_torsion_subgroup()([1,1/2,3/4,2]); P
            [(1, 1/2, 3/4, 2)]
            sage: P.order()
            4
        """
        v = self.__abvar.vector_space()(x)
        return self.element_class(self, v)

    def abelian_variety(self):
        """
        Return the abelian variety that this is the set of all torsion
        points on.

        OUTPUT: abelian variety

        EXAMPLES::

            sage: J0(23).qbar_torsion_subgroup().abelian_variety()
            Abelian variety J0(23) of dimension 2
        """
        return self.__abvar

