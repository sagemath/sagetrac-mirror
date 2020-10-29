# -*- coding: utf-8 -*-
r"""
Galois representations attached to Jacoians of hyperelliptic curves over
the rationals

More text here.

The following are the most useful functions for the class ``GaloisRepresentation``.

For the reducibility:

- ``is_reducible(p)``
- ``is_irreducible(p)``
- ``reducible_primes()``

For the image:

- ``is_surjective(p)``
- ``non_surjective()``
- ``image_type(p)``

For the classification of the representation

- ``is_semistable(p)``
- ``is_unramified(p, ell)``
- ``is_crystalline(p)``

EXAMPLES::

    sage: E = EllipticCurve('196a1')
    sage: rho = E.galois_representation()
    sage: rho.is_irreducible(7)
    True
    sage: rho.is_reducible(3)
    True
    sage: rho.is_irreducible(2)
    True
    sage: rho.is_surjective(2)
    False
    sage: rho.is_surjective(3)
    False
    sage: rho.is_surjective(5)
    True
    sage: rho.reducible_primes()
    [3]
    sage: rho.non_surjective()
    [2, 3]
    sage: rho.image_type(2)
    'The image is cyclic of order 3.'
    sage: rho.image_type(3)
    'The image is contained in a Borel subgroup as there is a 3-isogeny.'
    sage: rho.image_type(5)
    'The image is all of GL_2(F_5).'

For semi-stable curve it is known that the representation is
surjective if and only if it is irreducible::

    sage: E = EllipticCurve('11a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [5]
    sage: rho.reducible_primes()
    [5]

For cm curves it is not true that there are only finitely many primes for which the
Galois representation mod p is surjective onto `GL_2(\GF{p})`::

    sage: E = EllipticCurve('27a1')
    sage: rho = E.galois_representation()
    sage: rho.non_surjective()
    [0]
    sage: rho.reducible_primes()
    [3]
    sage: E.has_cm()
    True
    sage: rho.image_type(11)
    'The image is contained in the normalizer of a non-split Cartan group. (cm)'

REFERENCES:

- [Ser1972]_
- [Ser1987]_
- [Coj2005]_

AUTHORS:

- Barinder Singh Banwait, Armand Brumer, Hyun Jong Kim, Zev Klagsbrun,
  Jacob Mayle, Padmavathi Srinivasan, Isabel Vogt

"""

######################################################################
#                 Copyright (C) 2020 The Authors
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
######################################################################
from __future__ import print_function, absolute_import

from sage.structure.sage_object import SageObject
import sage.arith.all as arith
from sage.rings.fast_arith import prime_range
import sage.misc.all as misc
import sage.rings.all as rings
from sage.rings.all import RealField, GF

from math import sqrt
from sage.libs.pari.all import pari


class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representations
    attached to the Jacobian of a  hyperelliptic curve over the rational
    numbers.

    More text here.

    EXAMPLES::

        sage: R.<x>=QQ[]
        sage: f = x^5 + 17
        sage: C = HyperellipticCurve(f)
        sage: J = C.jacobian()
        sage: rho = J.galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Jacobian of Hyperelliptic Curve defined by y^2 + y = x^3 - x^2 - 10*x - 20 over Rational Field

    """

    def __init__(self, A):
        r"""
        see ``GaloisRepresentation`` for documentation

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: loads(rho.dumps()) == rho
            True

        """
        self.__image_type = {}
        self._A = A

    def __repr__(self):
        r"""
        string representation of the class

        EXAMPLES::

            sage: rho = EllipticCurve([0,1]).galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Elliptic Curve defined by y^2 = x^3 + 1 over Rational Field

        """
        return "Compatible family of Galois representations associated to the " + repr(self._E)

    def __eq__(self,other):
        r"""
        Compares two Galois representations.
        We define tho compatible families of representations
        attached to elliptic curves to be isomorphic if the curves are equal.

        EXAMPLES::

            sage: rho = EllipticCurve('11a1').galois_representation()
            sage: rho2 = EllipticCurve('11a2').galois_representation()
            sage: rho == rho
            True
            sage: rho == rho2
            False
            sage: rho == 34
            False

        """
        # if rho_E = rho_E' then the L-functions agree,
        # so E and E' are isogenous
        # except for p=2
        # the mod p representations will be different
        # for p dividing the degree of the isogeny
        # anyway, there should not be a _compatible_
        # isomorphism between rho and rho' unless E
        # is isomorphic to E'
        # Note that rho can not depend on the Weierstrass model
        if type(self) is not type(other):
            return False
        return self._E.is_isomorphic(other._E)

    def elliptic_curve(self):
        r"""
        The elliptic curve associated to this representation.

        EXAMPLES::

            sage: E = EllipticCurve('11a1')
            sage: rho = E.galois_representation()
            sage: rho.elliptic_curve() == E
            True

        """
        from copy import copy
        return copy(self._E)

#####################################################################
# surjectivity
#####################################################################

    def _init_exps():
        """
        Return a dictionary with keys l = 3, 5, and 7; for each l, the associated value is the list of characteristic polynomials of the matrices in the exceptional subgroup of GSp(4,l).
        """
        #char3 is the list of characteristic polynomials of matrices in the one subgroup of GSp(4,3) (up to conjugation) that isn't ruled out by surj_tests
        R.<x> = PolynomialRing(Zmod(3))
        char3 = [x^4 + 2*x^3 + x^2 + 2*x + 1,
            x^4 + 1,
            x^4 + x^3 + 2*x^2 + x + 1,
            x^4 + 2*x^3 + 2*x + 1,
            x^4 + x^3 + 2*x^2 + 2*x + 1,
            x^4 + x^3 + x^2 + x + 1,
            x^4 + x^3 + x^2 + 2*x + 1,
            x^4 + 2*x^2 + 1,
            x^4 + x^3 + x + 1,
            x^4 + 2*x^3 + 2*x^2 + x + 1,
            x^4 + 2*x^3 + 2*x^2 + 2*x + 1,
            x^4 + x^2 + 1,
            x^4 + 2*x^3 + x^2 + x + 1]
        #char5 Is the list of characteristic polynomials of matrices in the one subgroup of GSp(4,5) (up to conjugation) that isn't ruled out by surj_tests
        R.<x> = PolynomialRing(Zmod(5))
        char5 = [x^4 + x^3 + 2*x^2 + x + 1,
            x^4 + x^3 + 4*x + 1,
            x^4 + x^3 + x^2 + x + 1,
            x^4 + x^3 + x + 1,
            x^4 + 4*x^2 + 1,
            x^4 + 4*x^3 + 4*x^2 + 4*x + 1,
            x^4 + 3*x^2 + 1,
            x^4 + 4*x^3 + 3*x^2 + 4*x + 1,
            x^4 + 2*x^2 + 1,
            x^4 + 4*x^3 + 4*x^2 + x + 1,
            x^4 + 4*x^3 + 2*x^2 + 4*x + 1,
            x^4 + x^2 + 1,
            x^4 + 2*x^3 + 4*x^2 + 3*x + 1,
            x^4 + 4*x^3 + 3*x^2 + x + 1,
            x^4 + 1,
            x^4 + 4*x^3 + x^2 + 4*x + 1,
            x^4 + 2*x^3 + 3*x^2 + 3*x + 1,
            x^4 + 4*x^3 + 2*x^2 + x + 1,
            x^4 + 4*x^3 + 4*x + 1,
            x^4 + 2*x^3 + 2*x^2 + 3*x + 1,
            x^4 + 4*x^3 + x^2 + x + 1,
            x^4 + 3*x^3 + 4*x^2 + 3*x + 1,
            x^4 + 2*x^3 + x^2 + 3*x + 1,
            x^4 + 4*x^3 + x + 1,
            x^4 + 3*x^3 + 3*x^2 + 3*x + 1,
            x^4 + 2*x^3 + 3*x + 1,
            x^4 + 3*x^3 + 2*x^2 + 3*x + 1,
            x^4 + 3*x^3 + x^2 + 3*x + 1,
            x^4 + 3*x^3 + 3*x + 1,
            x^4 + 2*x^3 + 4*x^2 + 2*x + 1,
            x^4 + 2*x^3 + 3*x^2 + 2*x + 1,
            x^4 + 2*x^3 + 2*x^2 + 2*x + 1,
            x^4 + 3*x^3 + 4*x^2 + 2*x + 1,
            x^4 + 2*x^3 + x^2 + 2*x + 1,
            x^4 + x^3 + 4*x^2 + 4*x + 1,
            x^4 + 3*x^3 + 3*x^2 + 2*x + 1,
            x^4 + 2*x^3 + 2*x + 1,
            x^4 + x^3 + 3*x^2 + 4*x + 1,
            x^4 + 3*x^3 + 2*x^2 + 2*x + 1,
            x^4 + x^3 + 4*x^2 + x + 1,
            x^4 + 3*x^3 + x^2 + 2*x + 1,
            x^4 + x^3 + 2*x^2 + 4*x + 1,
            x^4 + x^3 + 3*x^2 + x + 1,
            x^4 + x^3 + x^2 + 4*x + 1,
            x^4 + 3*x^3 + 2*x + 1]
        #char7 Is the list of characteristic polynomials of matrices in the one subgroup of GSp(4,7) (up to conjugation) that isn't ruled out by surj_tests
        R.<x> = PolynomialRing(Zmod(7))
        char7 = [x^4 + 2*x^3 + 5*x^2 + 5*x + 1,
            x^4 + 5*x^3 + 5*x^2 + 2*x + 1,
            x^4 + x^3 + 6*x^2 + 2*x + 4,
            x^4 + x^3 + 3*x^2 + 4*x + 2,
            x^4 + 5*x^3 + 3*x^2 + 5*x + 1,
            x^4 + 1,
            x^4 + 2*x^2 + 1,
            x^4 + 6*x^2 + 1,
            x^4 + 4*x^3 + x + 4,
            x^4 + 4*x^3 + 2*x^2 + x + 4,
            x^4 + 6*x^3 + x^2 + 6*x + 1,
            x^4 + 4*x^3 + 5*x^2 + 2*x + 2,
            x^4 + x^3 + x + 1,
            x^4 + 3*x^3 + 5*x^2 + 5*x + 2,
            x^4 + 3*x^3 + 6*x + 4,
            x^4 + 3*x^3 + 2*x^2 + 6*x + 4,
            x^4 + 6*x^3 + 3*x^2 + 3*x + 2,
            x^4 + 2,
            x^4 + x^3 + 4*x^2 + 6*x + 1,
            x^4 + 3*x^2 + 4,
            x^4 + 6*x^2 + 2,
            x^4 + 5*x^2 + 4,
            x^4 + 3*x^3 + 6*x^2 + 3*x + 1,
            x^4 + 6*x^3 + 6*x^2 + 5*x + 4,
            x^4 + 3*x^3 + 4*x^2 + 4*x + 1,
            x^4 + 4*x^3 + 4*x^2 + 3*x + 1,
            x^4 + 4*x^3 + 6*x^2 + 4*x + 1,
            x^4 + 2*x^3 + x + 2,
            x^4 + x^3 + 2*x^2 + 3*x + 2,
            x^4 + 2*x^3 + 4*x^2 + x + 2,
            x^4 + 6*x^3 + 4*x^2 + x + 1,
            x^4 + 3*x^3 + x^2 + x + 4,
            x^4 + 2*x^3 + x^2 + 3*x + 4,
            x^4 + x^3 + 3*x^2 + 5*x + 4,
            x^4 + 5*x^2 + 1,
            x^4 + 2*x^3 + 5*x^2 + 4*x + 4,
            x^4 + 3*x^3 + 6*x^2 + 2*x + 2,
            x^4 + 6*x^3 + 6*x + 1,
            x^4 + 2*x^3 + 2*x^2 + 6*x + 2,
            x^4 + x^3 + x^2 + x + 1,
            x^4 + 5*x^3 + 2*x^2 + x + 2,
            x^4 + 5*x^3 + 5*x^2 + 3*x + 4,
            x^4 + 4*x^3 + 6*x^2 + 5*x + 2,
            x^4 + 4*x^3 + x^2 + 6*x + 4,
            x^4 + 5*x^3 + x^2 + 4*x + 4,
            x^4 + 2*x^3 + 3*x^2 + 2*x + 1,
            x^4 + 6*x^3 + 3*x^2 + 2*x + 4,
            x^4 + x^2 + 2,
            x^4 + 4,
            x^4 + 5*x^3 + 6*x + 2,
            x^4 + 3*x^2 + 2,
            x^4 + 6*x^3 + 2*x^2 + 4*x + 2,
            x^4 + 4*x^2 + 4,
            x^4 + 5*x^3 + 4*x^2 + 6*x + 2]
        return {3 : char3, 5 : char5, 7 : char7}


    def _init_wit(L):
        """
        Return a list for witnesses with all entries initially all set to zero, in the following format:
            2: [_] <-> [_is_surj_at_2 ]
            3: [_,_,_] <-> [witness for _surj_test_A, witness for _surj_test_B, witness for _surj_test_exp]
            5: [_,_,_] <-> [witness for _surj_test_A, witness for _surj_test_B, witness for _surj_test_exp]
            7: [_,_,_] <-> [witness for _surj_test_A, witness for _surj_test_B, witness for _surj_test_exp]
            3: [_,_,_] <-> [witness for _surj_test_A, witness for _surj_test_B]
        """
        witnesses = {}
        for l in L:
            if l == 2:
                witnesses[l] = [0]
            elif l in [3,5,7]:
                witnesses[l] = [0,0,0]
            else:
                witnesses[l] = [0,0]
        return witnesses


    def _is_surj_at_2(f,h):
        """
        Return True if and only if the mod 2 Galois image of the Jacobian of y^2 + h(x) y = f(x) is surjective, i.e. if and
        only if the Galois group of the polynomial 4*f+h^2 is all of S_6.
        """
        F = 4*f + h^2
        return F.is_irreducible() and F.galois_group().order() == 720


    def _surj_test_A(frob_mod):
        """
        Return True if frob_mod is irreducible.
        """
        return frob_mod.is_irreducible()


    def _surj_test_B(frob_mod):
        """
        Return True if frob_mod has nonzero trace and has a linear factor with multiplicity one.
        """
        if -frob_mod[3] != 0:
            for fact in frob_mod.factor():
                if fact[0].degree() == 1 and fact[1] == 1:
                    return True
        return False


    def _surj_test_exp(l, frob_mod, exps):
        """
        Return True if frob_mod is the characteristic polynomial of a matrix that is not in the exceptional subgroup mod l
        """
        return not frob_mod in exps[l]


    def _update_wit(l, p, frob, f, h, exps, wit):
        """
        Return an updated list of witnesses, based on surjectivity tests for frob at p.
        """
        frob_mod = frob.change_ring(Zmod(l))
        for i in range(0,len(wit)):
            if wit[i] == 0:
                if l == 2:
                    if _is_surj_at_2(f,h):
                        wit[i] = 1
                    else:
                        wit[i] = -1
                elif i == 0 and _surj_test_A(frob_mod):
                    wit[i] = p
                elif i == 1 and _surj_test_B(frob_mod):
                    wit[i] = p
                elif i == 2 and _surj_test_exp(l, frob_mod, exps):
                    wit[i] = p
        return wit


    def find_surj_from_list(self, L=list(primes(1000)), bound=1000, verbose=False):
        r"""
        Return a list of primes in L at which the residual Galois representation of the Jacobian of `H` might not be surjective.
        Outside of the returned list, all primes in L are surjective. The primes in the returned list are likely non-surjective.

        INPUT:

        - ``H`` -- hyperelliptic curve with typical Jacobian

        - ``L`` -- list of primes (default: list of primes less than 1000)

        - ``bound`` -- integer (default: `1000`)

        - ``verbose`` -- boolean (default: `False`)

        OUTPUT: a list

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([0, 1, 1]), R([1, 0, 0, 1]));
            sage: find_surj_from_list(H)
            [2, 7]

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurve(R([0, 21, -5, -9, 1, 1]), R([1, 1])); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x + 1)*y = x^5 + x^4 - 9*x^3 - 5*x^2 + 21*x
            sage: find_surj_from_list(H)
            [2, 13]

        """
        f,h = H.hyperelliptic_polynomials()
        # C = 2 * genus2reduction(h, f).conductor # An integer which agrees up with the conductor of H: y^2 + h y = f, except possibly at two. Bad primes of Jac(H) divide it.
        C = 2 * prod(genus2reduction(h,f).local_data.keys())
        witnesses = _init_wit(L)
        exps = _init_exps()
        to_check = L.copy() # to_check is the list of primes for which we still need to determine surjectivity. Initially, it equals L and we remove primes as their status is determined.
        for p in primes(3,bound):
            if C % p != 0:
                Hp = H.change_ring(GF(p))
                frob = Hp.frobenius_polynomial()
                to_remove = []
                for l in to_check:
                    if p != l and 0 in witnesses[l]:
                        witnesses[l] = _update_wit(l, p, frob, f, h, exps, witnesses[l])
                        if not 0 in witnesses[l]:
                            to_remove.append(l)
                for l in to_remove:
                    to_check.remove(l)
                if len(to_check) == 0:
                    break
        if verbose:
            return witnesses
        probably_non_surj_primes = []
        for l in L:
            if -1 in witnesses[l] or 0 in witnesses[l]:
                probably_non_surj_primes.append(l)
        return probably_non_surj_primes


    def is_surjective(self, p, bound=1000, verbose=False):
        r"""
        Return True if the mod-p representation is
        surjective onto `Aut(A[p]) = GSp_4(\GF{p})`.

        False if it is not, or None if we were unable to
        determine whether it is or not.

        INPUT:

        -  ``p`` - int (a prime number)

        -  ``A`` - int (a bound on the number of a_p to use)

        OUTPUT:

        - boolean. True if the mod-p representation is surjective
          and False if not.

        The answer is cached.

        EXAMPLES::

            sage: rho = EllipticCurve('37b').galois_representation()
            sage: rho.is_surjective(2)
            True
            sage: rho.is_surjective(3)
            False

        ::

            sage: rho = EllipticCurve('121a1').galois_representation()
            sage: rho.non_surjective()
            [11]
            sage: rho.is_surjective(5)
            True
            sage: rho.is_surjective(11)
            False

            sage: rho = EllipticCurve('121d1').galois_representation()
            sage: rho.is_surjective(5)
            False
            sage: rho.is_surjective(11)
            True

        Here is a case, in which the algorithm does not return an answer::

            sage: rho = EllipticCurve([0,0,1,2580,549326]).galois_representation()
            sage: rho.is_surjective(7)

        In these cases, one can use image_type to get more information about the image::

            sage: rho.image_type(7)
            'The image is contained in the normalizer of a split Cartan group.'


        REMARKS:

        1. If `p \geq 5` then the mod-p representation is
           surjective if and only if the p-adic representation is
           surjective. When `p = 2, 3` there are
           counterexamples. See papers of Dokchitsers and Elkies
           for more details.

        2. For the primes `p=2` and 3, this will always answer either
           True or False. For larger primes it might give None.

        """

        ans = find_surj_from_list(self, L=[p], bound=1000, verbose=False)

        if ans:
            return True
        else:
            return False

    #########################################################
    #                            #
    #                   Auxiliary functions            #
    #                            #
    #########################################################

    def p_part(p, N):
        if N != 0:
            return p^valuation(N,p)
        else:
            return 1


    def poor_mans_conductor(C):
        f,h = C.hyperelliptic_polynomials()
        red_data = genus2reduction(h,f)
        N = red_data.conductor
        if red_data.prime_to_2_conductor_only:
            N = 2*N
        return N

    def maximal_square_divisor(N):
        PP = prime_factors(N)
        n = 1
        for p in PP:
            n = n * p^(floor(valuation(N,p)/2))

        return n

    def true_conductor(C):
        print("Warning. This hasn't yet been implemented. The return value isn't actually the conductor.")
        return poor_mans_conductor(C)


    #########################################################
    #                            #
    #          Governed by a quadratic character        #
    #                            #
    #########################################################


    def maximal_quadratic_conductor(N):
        if N % 2 == 0:
            return 4*radical(N)
        else:
            return radical(N)

    def character_list(N):
        #N = poor_mans_conductor(C)
        c = maximal_quadratic_conductor(N)
        D = DirichletGroup(c,base_ring=QQ,zeta_order=2)
        return [phi for phi in D if phi.conductor() != 1]


    def set_up_quadratic_chars(N):
        return [(phi,0,0) for phi in character_list(N)]



    def rule_out_quadratic_ell_via_Frob_p(p,fp,MM):
        """Provide a summary of what this method is doing.

        Args:
            p (int): new prime
            fp (integer poly): charpoly of frobenius at p on a hyperelliptic curve
            MM (list): list of the form <phi,M,y>, where phi is a non-trivial
            quadratic character, all primes ell for which there is a quadratic
            obstruction associated with phi must divide M, y is a counter for the
        the number of nontrivial Frobenius constraints going into M

        Returns:
            (list): TODO
        """
        ap = -fp.coefficients(sparse=false)[3]
        if ap == 0:
            return MM
        else:
            MM0 = []
            for phi,M,y in MM:
                if (M == 1 or phi(p) != -1 or y>1):
                    MM0.append((phi,M,y))
                else:
                    MM0.append((phi,gcd(M,p*ap), y+1))

    #        MM0 = [(phi,M/(p_part(2,M)*p_part(3,M)), y) for phi,M,y in MM0]
            return MM0



    #########################################################
    #                            #
    #             Reducible (easy cases)            #
    #                            #
    #########################################################

    #This should probably be update with the GCD of the return value and 120
    def compute_f_from_d(d,p):
        return Integers(d)(p).multiplicative_order()

    def rule_out_one_dim_ell(p,fp,d,M):
        """Zev's direct implementation of what is in Dieulefait SS 3.1 and 3.2.
        Warning: some bugs because the conductor used here is wrong at 2

        Args:
            p ([type]): [description]
            fp ([type]): [description]
            d ([type]): [description]
            M ([type]): [description]

        Returns:
            [type]: [description]
        """
        if M != 1:
            f = compute_f_from_d(d,p)
            x = fp.variables()[0]
            M = gcd(M,p*fp.resultant(x^f-1))

        return M


    def rule_out_related_two_dim_ell_case1(p,fp,d,M):
        if M != 1:
            f = compute_f_from_d(d,p)
            ap = -fp.coefficients(sparse=false)[3]
            bp = fp.coefficients(sparse=false)[2]
            x = fp.variables()[0]

            Q = (x*bp - 1 - p^2*x^2)*(p*x + 1)^2 - ap^2*p*x^2
            M = gcd(M,p*Q.resultant(x^f-1))

        return M

    def rule_out_related_two_dim_ell_case2(p,fp,d,M):
        if M != 1:
            f = compute_f_from_d(d,p)
            ap = -fp.coefficients(sparse=false)[3]
            bp = fp.coefficients(sparse=false)[2]
            x = fp.variables()[0]

            Q = (x*bp - p - p*x^2)*(x + 1)^2 - ap^2*x^2
            M = gcd(M,p*Q.resultant(x^f-1))

        return M

    """
    Isabel's code using the fact that the exponent of the determinant
    character divides 120.

    the following three functions implement the following for n=2,3,5:
    f(x) = x^4 - t*x^3 + s*x^2 - p^alpha*t*x + p^(2*alpha) is a
    degree 4 polynomial whose roots multiply in pairs to p^alpha
    returns the tuple (p, tn, sn, alphan) of the polynomial
    f^(n)(x) = x^4 - tn*x^3 + sn*x^2 - p^(alphan)*tn*x + p^(2*alphan)
    whose roots are the nth powers of the roots of f
    """

    def power_roots2(ptsa):
        p, t, s, alpha = ptsa
        return (p, t^2 - 2*s, s^2 - 2*p^alpha*t^2 + 2*p^(2*alpha), 2*alpha)


    def power_roots3(ptsa):
        p, t, s, alpha = ptsa
        return (p, t^3 - 3*s*t + 3*p^alpha*t, s^3 - 3*p^alpha*s*t^2 + 3*p^(2*alpha)
                *t^2 + 3*p^(2*alpha)*t^2 - 3*p^(2*alpha)*s, 3*alpha)


    def power_roots5(ptsa):
        p, t, s, alpha = ptsa
        return (p, t^5 - 5*s*t^3 + 5*s^2*t + 5*p^alpha*t^3 - 5*p^alpha*s*t -
        5*p^(2*alpha)*t, s^5 - 5*p^alpha*s^3*t^2 + 5*p^(2*alpha)*s*t^4 +
        5*p^(2*alpha)*s^2*t^2 - 5*p^(3*alpha)*t^4 + 5*p^(2*alpha)*s^2*t^2 -
        5*p^(2*alpha)*s^3 - 5*p^(3*alpha)*t^4 - 5*p^(3*alpha)*s*t^2 +
        5*p^(4*alpha)*t^2 + 5*p^(4*alpha)*t^2 + 5*p^(4*alpha)*s, 5*alpha)


    #put these together to do any power dividing 120 that we actually need
    #c is the power
    def power_roots(cptsa):
        c, p, t, s, alpha = cptsa
        if 120 % c != 0:
            raise ValueError("can't raise to this power")

        ptsa = (p, t, s, alpha)

        while c % 2 == 0:
            c, ptsa = c/2, power_roots2(ptsa)

        while c % 3 == 0:
            c, ptsa = c/3, power_roots3(ptsa)

        while c % 5 == 0:
            c, ptsa = c/5, power_roots5(ptsa)

        return ptsa


    #given a quartic f whose roots multiply to p^alpha in pairs,
    #returns the quartic whose roots are the products of roots
    #of f that DO NOT multiply to p^alpha
    def roots_pairs_not_p(ptsa):
        p, t, s, alpha = ptsa
        return (p, s - 2*p, p*t^2 - 2*p*s + 2*p^2, 2*alpha)

    #t and s are the first and second elementary symmetric functions in the
    #roots of the characteristic polynomial of Frobenius at p for a curve C
    #M is an integer such that every prime ell for which J_C[ell] could be
    #1 + 3 reducible divides M
    #y is a counter for the number of nontrivial Frobenius conditions going
    #into M
    def rule_out_1_plus_3_via_Frob_p(c, p, t, s, M=0, y=0):
        p, tnew, snew, alphanew = power_roots((c, p, t, s, 1))
        Pnew(x) = x^4 - tnew*x^3 + snew*x^2 - p^alphanew*tnew*x + p^(2*alphanew)
        Pval = Pnew(1)
        if Pval != 0:
            return ZZ(gcd(M, p*Pval)), y+1
        else:
            return M, y

    #t and s are the first and second elementary symmetric functions in the
    #roots of the characteristic polynomial of Frobenius at p for a curve C
    #M is an integer such that every prime ell for which J_C[ell] could be
    #2+2 non-self-dual type reducible divides M
    #y is a counter for the number of nontrivial Frobenius conditions going
    #into M
    def rule_out_2_plus_2_nonselfdual_via_Frob_p(c, p, t, s, M=0, y=0):
        p, tnew, snew, alphanew = roots_pairs_not_p((p, t, s, 1))
        p, tnew, snew, alphanew = power_roots((c, p, tnew, snew, alphanew))
        Pnew(x) = x^4 - tnew*x^3 + snew*x^2 - p^alphanew*tnew*x + p^(2*alphanew)
        Pval = Pnew(1)*Pnew(p^c)
        if Pval != 0:
            return ZZ(gcd(M, p*Pval)), y+1
        else:
            return M, y



    #########################################################
    #                            #
    #          Reducible (modular forms case)        #
    #                            #
    #########################################################


    def special_divisors(N):
        D0 = [d for d in divisors(N) if d <= sqrt(N)]
        D0.reverse()
        D = []
        for d0 in D0:
                if all([d % d0 != 0 for d in D]):
                        D.append(d0)

        D.reverse()
        return D


    def get_cuspidal_levels(N, max_cond_exp_2=None):

        if max_cond_exp_2 is not None:
            # if we're here, then N is the even poor mans conductor
            conductor_away_two = N/2  # recall we put a 2 in the poor mans conductor
            possible_conductors = [conductor_away_two * 2^i for i in range(max_cond_exp_2 + 1)]
            return list(set([d for N in possible_conductors for d in special_divisors(N)]))  # not ordered, hopefully not a problem.
        else:
            return special_divisors(N)


    def create_polynomial_database(path_to_datafile, levels_of_interest):

        df = pd.read_csv(path_to_datafile, sep=":", header=None, names=["N", "p", "coeffs"])
        actual_levels_of_interest = [i for i,j,k in levels_of_interest]
        df_relevant = df.loc[df["N"].isin(actual_levels_of_interest)].copy()

        df_relevant["coeffs"] = df_relevant["coeffs"].apply(ast.literal_eval)

        return df_relevant


    def set_up_cuspidal_spaces(N, path_to_datafile=None, max_cond_exp_2=None):
        D = get_cuspidal_levels(N, max_cond_exp_2)
        if path_to_datafile is not None:
            levels_of_interest = [(d,0,0) for d in D]

            # There are no modular forms of level < 11
            bad_levels = [(i,0,0) for (i,0,0) in levels_of_interest if i in range(11)]

            levels_of_interest = [z for z in levels_of_interest if z not in bad_levels]

            DB = create_polynomial_database(path_to_datafile, levels_of_interest)
            return levels_of_interest, DB
        else:
            return [(CuspForms(d),0,0) for d in D], None


    def reconstruct_hecke_poly_from_trace_polynomial(cusp_form_space, p):
        """Implement Zev and Joe Wetherell's idea"""

        char_T_x = R(cusp_form_space.hecke_polynomial(p))
        S.<a,b> = QQ[]
        char_T_a_b = S(char_T_x(x=a)).homogenize(var='b')
        substitute_poly = char_T_a_b(a=1+p*b^2)

        return R(substitute_poly(a=0, b=x))


    def get_hecke_characteristic_polynomial(cusp_form_space, p, coeff_table = None):
        """This should return the left hand side of Equation 3.8.

        Args:
            cusp_form_space ([type]): either a space of weight 2 cuspforms with trivial
            Nebentypus or a level (given as an integer)
            p (int): prime number
        coeff_table: a filename with a list of characteristic polyomials for spaces of modular forms

        Returns:
            [pol]: an integer polynomial of twice the dimension of the cuspform
                space
        """

        if coeff_table is None:
            return reconstruct_hecke_poly_from_trace_polynomial(cusp_form_space, p)
        else:

            slice_of_coeff_table = coeff_table.loc[(coeff_table["N"] == cusp_form_space)
                                                & (coeff_table["p"] == p)]

            if slice_of_coeff_table.shape[0] == 1:
                print("doing pandas stuff for level {}".format(cusp_form_space))
                hecke_charpoly_coeffs = slice_of_coeff_table.iloc[int(0)]["coeffs"]
                hecke_charpoly = sum([hecke_charpoly_coeffs[i]*(R.0)^i for i in range(len(hecke_charpoly_coeffs))])
            else:  # i.e., can't find data in database
                warning_msg = ("Warning: couldn't find level {} and prime {} in DB.\n"
                            "Reconstructing on the fly...").format(cusp_form_space, p)
                print(warning_msg)
                CuspFormSpaceOnFly = CuspForms(cusp_form_space)
                hecke_charpoly = reconstruct_hecke_poly_from_trace_polynomial(CuspFormSpaceOnFly, p)

            return hecke_charpoly


    def rule_out_cuspidal_space_using_Frob_p(S,p,fp,M,y,coeff_table=None):
        if M != 1 and y<2:
            Tp = get_hecke_characteristic_polynomial(S,p, coeff_table=coeff_table)
            res = fp.resultant(Tp)
            if res != 0:
                return gcd(M,p*res), y+1
        return M, y


    def rule_out_cuspidal_spaces_using_Frob_p(p,fp,MC, coeff_table = None):
        MC0 = []
        for S,M,y in MC:
            Mm, yy = rule_out_cuspidal_space_using_Frob_p(S,p,fp,M,y,coeff_table=coeff_table)
            MC0.append((S,Mm,yy))
        return MC0


    def non_surjective(self, N=None, path_to_datafile=None, bound=1000):
        r"""
        Returns a list of primes p such that the mod-p representation
        *might* not be surjective. If `p` is not in the returned list,
        then the mod-p representation is provably surjective.

        By a theorem of Serre, there are only finitely
        many primes in this list, except when the curve has
        complex multiplication.

        If the curve has CM, we simply return the
        sequence [0] and do no further computation.

        INPUT:

        - ``A`` - an integer (default 1000). By increasing this parameter
          the resulting set might get smaller.

        OUTPUT:

        -  ``list`` - if the curve has CM, returns [0].
           Otherwise, returns a list of primes where mod-p representation is
           very likely not surjective. At any prime not in this list, the
           representation is definitely surjective.

        EXAMPLES::

            sage: E = EllipticCurve([0, 0, 1, -38, 90])  # 361A
            sage: E.galois_representation().non_surjective()   # CM curve
            [0]

        ::

            sage: E = EllipticCurve([0, -1, 1, 0, 0]) # X_1(11)
            sage: E.galois_representation().non_surjective()
            [5]

            sage: E = EllipticCurve([0, 0, 1, -1, 0]) # 37A
            sage: E.galois_representation().non_surjective()
            []

            sage: E = EllipticCurve([0,-1,1,-2,-1])   # 141C
            sage: E.galois_representation().non_surjective()
            [13]

        ::

            sage: E = EllipticCurve([1,-1,1,-9965,385220]) # 9999a1
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [2]

            sage: E = EllipticCurve('324b1')
            sage: rho = E.galois_representation()
            sage: rho.non_surjective()
            [3, 5]

        ALGORITHM:
        We first find an upper bound `B` on the possible primes. If `E`
        is semi-stable, we can take `B=11` by a result of Mazur. There is
        a bound by Serre in the case that the `j`-invariant is not integral
        in terms of the smallest prime of good reduction. Finally
        there is an unconditional bound by Cojocaru, but which depends
        on the conductor of `E`.
        For the prime below that bound we call ``is_surjective``.

        """

        C = self._A.curve()

        M1p3 = 0
        y1p3 = 0
        M2p2nsd = 0
        y2p2nsd = 0

        if N is None:
            f,h = C.hyperelliptic_polynomials()
            red_data = genus2reduction(h,f)
            N = red_data.conductor  # is this the true conductor if red_data.prime_to_2_conductor_only is False?
            max_cond_exp_2 = None
            if red_data.prime_to_2_conductor_only:
                # I think this is the case where we don't know exactly the two-part of conductor
                N = 2*N
                max_cond_exp_2 = red_data.minimal_disc.valuation(2)

        #MCusp is a list of the form <S,M,y>, where S is either a space of cusp forms or
        #a level, M is an integer such that all primes with a reducible sub isomorphic to the
        #rep of a cusp form in S divide M, y is a counter for the number of nontrivial Frobenius
        #conditions go into M
        MCusp, DB = set_up_cuspidal_spaces(N, path_to_datafile=path_to_datafile, max_cond_exp_2=None)

        #MQuad is a list of the form <phi,M,y>, where phi is a quadratic character, M is the integer
        #all nonsurjective primes governed by phi must divide, and y is counter for the number of nontrivial
        #Frobenius conditions going into M
        MQuad = set_up_quadratic_chars(N)

        d = maximal_square_divisor(N)

        #we'll test as many p as we need to get at least 2 nontrivial Frobenius conditions for every
        #possible cause of non-surjectivity
        sufficient_p = False

        p = 1

        while not sufficient_p:
                p = next_prime(p)
                if N % p != 0:
                    Cp = C.change_ring(FiniteField(p))
                    fp = Cp.frobenius_polynomial()
                    fp_rev = Cp.zeta_function().numerator()

                    #M31 = rule_out_one_dim_ell(p,fp,d,M31);
                    #M32A = rule_out_related_two_dim_ell_case1(p,fp,d,M32A)
                    #M32B = rule_out_related_two_dim_ell_case2(p,fp,d,M32B)

                    f = Integers(d)(p).multiplicative_order()
                    c = gcd(f, 120)
                    c = lcm(c, 8)  # adding in the max power of 2
                    tp = - fp.coefficients(sparse=false)[3]
                    sp = fp.coefficients(sparse=false)[2]

                    M1p3, y1p3 = rule_out_1_plus_3_via_Frob_p(c, p, tp, sp, M1p3, y1p3)
                    M2p2nsd, y2p2nsd = rule_out_2_plus_2_nonselfdual_via_Frob_p(c, p, tp, sp, M2p2nsd, y2p2nsd)
                    MCusp = rule_out_cuspidal_spaces_using_Frob_p(p,fp_rev,MCusp,coeff_table=DB)
                    MQuad = rule_out_quadratic_ell_via_Frob_p(p,fp,MQuad)

                if (M1p3 == 1) or (y1p3 > 1):
                    if (M2p2nsd == 1) or (y2p2nsd > 1):
                        if all((Mc == 1 or yc>1) for S, Mc, yc in MCusp):
                            if all((Mq == 1 or yq > 1) for phi, Mq, yq in MQuad):
                                sufficient_p = True


        #ell_red_easy = [prime_factors(M31), prime_factors(M32A), prime_factors(M32B)]

        # we will always include 2, 3, 5, 7 and the non-semistable primes.
        non_maximal_primes = {2,3,5,7}.union(set([p[0] for p in list(N.factor()) if p[1]>1]))

        ell_red_easy = [M1p3.prime_factors(), M2p2nsd.prime_factors()]
        non_maximal_primes = non_maximal_primes.union(set([p for j in ell_red_easy for p in j]))

        if path_to_datafile is not None:
            ell_red_cusp = [(S,prime_factors(M)) for S,M,y in MCusp]
        else:
            ell_red_cusp = [(S.level(),prime_factors(M)) for S,M,y in MCusp]

        non_maximal_primes = non_maximal_primes.union(set([p for a,j in ell_red_cusp for p in j]))

        ell_irred = [(phi,prime_factors(M)) for phi,M,t in MQuad]
        non_maximal_primes = non_maximal_primes.union(set([p for a,j in ell_irred for p in j]))

        return self.find_surj_from_list(L=non_maximal_primes, bound)