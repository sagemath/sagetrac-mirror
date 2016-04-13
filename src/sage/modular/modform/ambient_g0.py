r"""
Modular Forms for `\Gamma_0(N)` over `\QQ`

TESTS::

    sage: m = ModularForms(Gamma0(389),6)
    sage: loads(dumps(m)) == m
    True
"""

#########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

import sage.rings.all as rings

import sage.modular.arithgroup.all as arithgroup

import ambient
import cuspidal_submodule
import eisenstein_submodule
from itertools import combinations

from sage.sets.set import Set
from sage.misc.misc_c import prod
from sage.functions.other import floor
from sage.arith.all import GCD, factor
from sage.rings.integer_ring import ZZ
from sage.misc.mrange import cartesian_product_iterator

from sage.rings.number_field.order import primitive_ideal_number


class ModularFormsAmbient_g0_Q(ambient.ModularFormsAmbient):
    """
    A space of modular forms for `\Gamma_0(N)` over `\QQ`.
    """
    def __init__(self, level, weight):
        r"""
        Create a space of modular symbols for `\Gamma_0(N)` of given
        weight defined over `\QQ`.

        EXAMPLES::

            sage: m = ModularForms(Gamma0(11),4); m
            Modular Forms space of dimension 4 for Congruence Subgroup Gamma0(11) of weight 4 over Rational Field
            sage: type(m)
            <class 'sage.modular.modform.ambient_g0.ModularFormsAmbient_g0_Q_with_category'>
        """
        ambient.ModularFormsAmbient.__init__(self, arithgroup.Gamma0(level), weight, rings.QQ)

    ####################################################################
    # Computation of Special Submodules
    ####################################################################
    def cuspidal_submodule(self):
        r"""
        Return the cuspidal submodule of this space of modular forms for
        `\Gamma_0(N)`.

        EXAMPLES::

            sage: m = ModularForms(Gamma0(33),4)
            sage: s = m.cuspidal_submodule(); s
            Cuspidal subspace of dimension 10 of Modular Forms space of dimension 14 for Congruence Subgroup Gamma0(33) of weight 4 over Rational Field
            sage: type(s)
            <class 'sage.modular.modform.cuspidal_submodule.CuspidalSubmodule_g0_Q_with_category'>
        """
        try:
            return self.__cuspidal_submodule
        except AttributeError:
            if self.level() == 1:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_level1_Q(self)
            else:
                self.__cuspidal_submodule = cuspidal_submodule.CuspidalSubmodule_g0_Q(self)
        return self.__cuspidal_submodule

    def eisenstein_submodule(self):
        r"""
        Return the Eisenstein submodule of this space of modular forms for
        `\Gamma_0(N)`.

        EXAMPLES::

            sage: m = ModularForms(Gamma0(389),6)
            sage: m.eisenstein_submodule()
            Eisenstein subspace of dimension 2 of Modular Forms space of dimension 163 for Congruence Subgroup Gamma0(389) of weight 6 over Rational Field
        """
        try:
            return self.__eisenstein_submodule
        except AttributeError:
            self.__eisenstein_submodule = eisenstein_submodule.EisensteinSubmodule_g0_Q(self)
        return self.__eisenstein_submodule

    def _compute_atkin_lehner_matrix(self, d):
        r"""
        Compute the matrix of the Atkin-Lehner involution W_d acting on self,
        where d is a divisor of the level.  This is only implemented in the
        (trivial) level 1 case.

        EXAMPLE::

            sage: ModularForms(1, 30).atkin_lehner_operator()
            Hecke module morphism Atkin-Lehner operator W_1 defined by the matrix
            [1 0 0]
            [0 1 0]
            [0 0 1]
            Domain: Modular Forms space of dimension 3 for Modular Group SL(2,Z) ...
            Codomain: Modular Forms space of dimension 3 for Modular Group SL(2,Z) ...
        """
        if self.level() == 1:
            from sage.matrix.matrix_space import MatrixSpace
            return MatrixSpace(self.base_ring(), self.rank())(1)
        else:
            raise NotImplementedError

    ####################################################################
    # dimension of eigenspaces of Atkin-Lehner involutions
    ####################################################################

    def atkin_lehner_ramification(self, Q):
        r"""
        Return the number of points stabilized by the Atkin-Lehner
        involution `w_Q` on `X_0(N)`.

        EXAMPLES::

            sage: ModularForms(100, 2).atkin_lehner_ramification(4)
            8
            sage: ModularForms(144169).atkin_lehner_ramification(144169)
            208

        TESTS::

            sage: ModularForms(100, 2).atkin_lehner_ramification(5)
            Traceback (most recent call last):
            ...
            ValueError: argument Q must exactly divide the level
        """
        N = self.level()
        M = floor(N / Q)
        if (N % Q) != 0 or GCD(M, Q) != 1:
            raise ValueError("argument Q must exactly divide the level")
        else:
            r = Q.squarefree_part()

        R = ZZ.zero()
        if Q == 1:
            return R
        elif Q == 2:
            # D = -4, will catch D = -8 below.
            R += primitive_ideal_number(-4, M)
        elif Q == 3:
            # D = -3 and D = -12 and return.
            if M == 1:
                return 2
            elif (M % 4) == 2:
                # D = -3 to D = -12 inclusion and vice versa.
                return 2 * primitive_ideal_number(-3, (M / 2).floor())
            else:
                return (primitive_ideal_number(-3, M)
                        + primitive_ideal_number(-12, M))
        elif Q == 4:
            # Add in the ramified cusps.
            R += (arithgroup.Gamma0((N / 4).floor()).ncusps())
        elif (Q % 4) == 3:
            if (M % 2) == 0:
                # Contribution of (-Q |--> -Q) ideals, and of ideals
                # (-4Q |--> -Q) and duals (-Q |--> -4Q).
                R += ((primitive_ideal_number(-Q, M) +
                       2 * primitive_ideal_number(-Q, floor(M / 2)) *
                       (2 - (-Q).kronecker(2))) * (-Q).class_number())
            else:
                R += primitive_ideal_number(-r, M) * ((-Q).class_number())
        return R + primitive_ideal_number(-4 * Q, M) * ((-4 * Q).class_number())

    def atkin_lehner_ramification_B(self):
        r"""
        Return the sequence of the number of points stabilized by
        the Atkin-Lehner involution `w_Q` on `X_0(N)`, for each positive
        integer `Q` exactly dividing `N`.

        EXAMPLES::

            sage: ModularForms(100, 2).atkin_lehner_ramification_B()
            [0, 0, 8, 4]
            sage: ModularForms(144169).atkin_lehner_ramification_B()
            [0, 208]
        """
        qs = [d ** e for d, e in factor(self.level())]
        ex = cartesian_product_iterator([[0, 1]] * len(qs))
        Qs = [prod([de ** e[i] for i, de in enumerate(qs)]) for e in ex]
        return [self.atkin_lehner_ramification(Q) for Q in Qs]

    def _atkin_lehner_closure(self, S):
        """
        Helper function for :meth:`atkin_lehner_quotient_genus_X0`.

        INPUT:

        - S -- a set of numbers

        EXAMPLES::

            sage: ModularForms(100, 2)._atkin_lehner_closure(set([4]))
            {1, 4}
            sage: ModularForms(120, 2)._atkin_lehner_closure(set([24, 20]))
            {1, 4, 5, 6, 20, 24, 30, 120}
            sage: ModularForms(120, 2)._atkin_lehner_closure(set([6, 20]))
            {1, 4, 5, 6, 20, 30}
            sage: ModularForms(120, 2)._atkin_lehner_closure(set([6, 7, 20]))
            {1, 4, 5, 6, 7, 20, 28, 30, 35, 42, 140, 210}
        """
        stop = False
        while not stop:
            T = set([P * Q // P.gcd(Q) ** 2 for P, Q in combinations(S, 2)])
            ClosureS = S.union(T)
            if len(ClosureS) == len(S):
                stop = True
            else:
                S = ClosureS
        ClosureS.add(ZZ.one())
        return ClosureS

    def atkin_lehner_quotient_genus_X0(self, Q):
        """
        Return the genus of the quotient of `X_0(N)` by the Atkin-Lehner
        operator `w_Q`.

        EXAMPLES::

            sage: ModularForms(100, 2).atkin_lehner_quotient_genus_X0(4)
            2
            sage: ModularForms(100, 2).atkin_lehner_quotient_genus_X0(25)
            4
        """
        N = self.level()
        if N == 1:
            return 0
        R = self.atkin_lehner_ramification(Q)
        return floor((arithgroup.Gamma0(N).genus() + 1 - floor(R / 2)) / 2)

    def atkin_lehner_quotient_many_genus_X0(self, S):
        """
        Return the genus of the quotient of `X_0(N)` by the Atkin-Lehner
        operators `w_Q` for `Q` in `S`.

        EXAMPLES::

            sage: ModularForms(100, 2).atkin_lehner_quotient_many_genus_X0([4])
            2
            sage: ModularForms(100, 2).atkin_lehner_quotient_many_genus_X0([4, 25])
            1
        """
        N = self.level()
        if N == 1:
            return 0
        S = self._atkin_lehner_closure(set(S))
        t = ZZ(len(S)).valuation(2)
        R = sum([self.atkin_lehner_ramification(Q) for Q in S])
        return floor((arithgroup.Gamma0(N).genus()
                      - 1 - floor(R / 2)) / 2 ** t) + 1

    def atkin_lehner_quotient_genus_X0_list(self):
        """
        Return the genus of the quotient of `X_0(N)` by the Atkin-Lehner
        operators.

        EXAMPLES::

            sage: ModularForms(100, 2).atkin_lehner_quotient_genus_X0_list()
            1
        """
        N = self.level()
        return self.atkin_lehner_quotient_many_genus_X0([d ** e
                                                         for d, e in factor(N)])

    def _character_kernel(self, w):
        """
        Helper function for :meth:`modular_genus_X0`.

        EXAMPLES::

            sage: ModularForms(100, 2)._character_kernel([1,2,4])
            [(0, 0, 0), (1, 0, 0)]
            sage: ModularForms(100, 2)._character_kernel([1,2,4,6])
            [(0, 0, 0, 0), (1, 0, 0, 0)]
        """
        r = len(w)
        return [d for d in cartesian_product_iterator([[0, 1]] * r)
                if 1 == prod([w[i] ** d[i] for i in range(r)])]

    def modular_genus_X0(self, w):
        r"""
        Return the dimension of the space of differentials on `X_0(N)`
        with character `w` under the Atkin-Lehner operators.

        INPUT:

        - w -- a list of elements of [-1, 1]

        EXAMPLES::

            sage: ModularForms(71,2).modular_genus_X0([1]) # largest prime for which this holds
            0
            sage: ModularForms(10**5,2).modular_genus_X0([-1, -1])
            3675

        TESTS::

            sage: ModularForms(10**5,2).modular_genus_X0([-1, -1, -1])
            Traceback (most recent call last):
            ...
            ValueError: length of w is not the number of prime divisors
        """
        PrPowDiv = [d ** ex for d, ex in factor(self.level())]
        if not len(w) == len(PrPowDiv):
            raise ValueError('length of w is not the number of prime divisors')
        gplus = self.atkin_lehner_quotient_genus_X0_list()
        if all([e == 1 for e in w]):
            return gplus
        D = self._character_kernel(w)
        W = [prod([PrPowDiv[i] ** d[i] for i in range(len(PrPowDiv))])
             for d in D]
        return self.atkin_lehner_quotient_many_genus_X0(W) - gplus

    def atkin_lehner_eigenspace_dimensions(self):
        r"""
        Return the dimensions of the Atkin-Lehner eigenspaces of
        weight 2 modular forms for `\Gamma_0(N)`.

        EXAMPLES::

            sage: ModularForms(1000003).atkin_lehner_eigenspace_dimensions()
            [41562, 41771]
            sage: ModularForms(12345).atkin_lehner_eigenspace_dimensions()
            [185, 226, 214, 198, 204, 206, 220, 192]
        """
        N = self.level()
        r = len(factor(N))
        eigs = [[w[i] for i in range(r)]
                for w in cartesian_product_iterator([[1, -1]] * r)]
        return [self.modular_genus_X0(w) for w in eigs]

    def new_subspace_dimensionX0(self):
        r"""
        Return the dimension of the new subspace of weight 2 modular
        forms on `X_0(N)`.

        EXAMPLES::

            sage: ModularForms(100000, 2).new_subspace_dimensionX0()
            2400
            sage: ModularForms(54321).new_subspace_dimensionX0()
            2855
            sage: ModularForms(12345).new_subspace_dimensionX0()
            547
        """
        N = self.level()
        Genera = []
        DivSeq = N.divisors()
        NewDims = [0] * len(DivSeq)
        for k in range(len(DivSeq)):
            M = DivSeq[k]
            Genera.append(arithgroup.Gamma0(M).genus())
            if M in range(1, 11) + [12, 13, 16, 18, 25]:
                NewDims[k] = 0
            else:
                OldDims = sum([NewDims[j] * len(floor(M / DivSeq[j]).divisors())
                               for j in range(k) if (M % DivSeq[j]) == 0])
                NewDims[k] = Genera[k] - OldDims
        return NewDims[len(NewDims) - 1]

    def atkin_lehner_new_eigenspace_dimension(self, w):
        r"""
        Return the dimension of the eigenspace of `w` in the new subspace
        of level `N`.

        EXAMPLES::

            sage: ModularForms(315, 2).atkin_lehner_new_eigenspace_dimension([-1,-1,-1])
            3
            sage: ModularForms(88, 2).atkin_lehner_new_eigenspace_dimension([-1,1])
            2
        """
        N = self.level()
        if not len(w) == len(N.prime_divisors()):
            print w, N.prime_divisors()
            raise ValueError("Length of argument w must equal the number of "
                             "primes dividing the level")
        DivSeq = N.divisors()
        return (self.modular_genus_X0(w)
                - sum([self.old_subspace_dimension(M, R, w)
                       for R in DivSeq for M in DivSeq
                       if (N % (M * R ** 2)) == 0]))

    def old_subspace_dimension(self, M, R, w):
        r"""
        Return the dimension of the (M, R)-old subspace of `S_2(\Gamma_0(N))`.

        EXAMPLES::

            sage: ModularForms(77,2).old_subspace_dimension(11,1,[-1,-1])
            1
            sage: ModularForms(98,2).old_subspace_dimension(49,1,[1,1])
            0
        """
        N = self.level()
        if not (N % (M * R ** 2)) == 0:
            raise ValueError("bad input")
        if (N == M) or (M == 1):
            return 0
        S = N // (M * R ** 2)
        SuppN = N.prime_divisors()
        if not len(w) == len(SuppN):
            raise ValueError('length of w is not the number of prime divisors')
        SuppM = [p for p in SuppN if (M % p) == 0]
        SuppS = [p for p in SuppN if (S % p) == 0]
        TogglePrimes = [p for p in SuppN if not(p in SuppS or p in SuppM)]
        if any([w[SuppN.index(p)] != 1 for p in TogglePrimes]):
            return 0
        BoundPrimes = [p for p in SuppM if not(p in SuppS)]
        cartesian_power = cartesian_product_iterator([[1, -1]] * len(SuppM))
        KernChars_seq = [[u[i] for i in range(len(SuppM))]
                         for u in cartesian_power
                         if all([1 == u[SuppM.index(p)] for p in BoundPrimes])]
        pull_back_char = [w[SuppN.index(p)] for p in SuppM]
        # here below the level is M
        modf_M = ModularFormsAmbient_g0_Q(M, 2)

        char_coset = [[e * t[i] for i, e in enumerate(pull_back_char)]
                      for t in KernChars_seq]
        return sum([modf_M.atkin_lehner_new_eigenspace_dimension(v)
                    for v in char_coset])

    def atkin_lehner_new_eigenspace_dimensions(self):
        r"""
        Return the dimensions of the new Atkin-Lehner eigenspaces
        for `\Gamma_0(N)` and weight 2.

        EXAMPLES::

            sage: ModularForms(5077, 2).atkin_lehner_new_eigenspace_dimensions()
            [206, 216]
            sage: ModularForms(315, 2).atkin_lehner_new_eigenspace_dimensions() # optional, long
            [2, 0, 2, 0, 2, 1, 0, 3]
        """
        N = self.level()
        cartesian_power = cartesian_product_iterator([[1, -1]] *
                                                     len(N.prime_divisors()))
        return [self.atkin_lehner_new_eigenspace_dimension(v)
                for v in cartesian_power]
