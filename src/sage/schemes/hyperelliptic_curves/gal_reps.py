# -*- coding: utf-8 -*-
r"""
Galois representations attached to Jacoians of hyperelliptic curves of
genus 2 over the rationals.

The Jacobian variety of a genus 2 hyperelliptic curve over `\QQ` defines
an abelian surface whose `p`-torsion Galois module defines a
`4`-dimensional Galois representation of the absolute Galois group
`G_{\QQ}` of `\QQ`. The image of this representation is a subgroup of
`\text{GSp}_4(\GF{p})`, and it is known by an analogue of Serre's Open Image
Theorem that this representation is surjective for almost all primes `p`
under the additional assumption that the Jacobian is "generic", meaning
that the ring of endomorphisms defined over `\overline{\QQ}` is `\ZZ`.

Currently sage can decide whether or not the image of this representation
associated to a generic jacobian is surjective, and moreover can
determine exactly the finitely many primes at which the representation
is not surjective. This is based on Dieulefait's algorithm described in
[Di2002]_

For the surjectivity at one prime:

- ``is_surjective(p)``

For the list of non-surjective primes:

- ``non_surjective()``

EXAMPLES::

    sage: R.<x>=QQ[]
    sage: f = x^6 + 2*x^3 + 4*x^2 + 4*x + 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: rho = A.galois_representation()
    sage: rho.non_surjective()  # long time
    [2, 7]
    sage: rho.is_surjective(7)  # long time
    False
    sage: rho.is_surjective(2)
    False
    sage: rho.is_surjective(3)
    True
    sage: rho.is_surjective(13)
    True

The curve 587.a.587.1 on the LMFDB has surjective mod-$\ell$ Galois
representation at all primes `\ell`. ::

    sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, -1, -1]), R([1, 1, 0, 1]));
    sage: A = C.jacobian()
    sage: rho = A.galois_representation()
    sage: rho.non_surjective()
    []

If the Jacobian has any non-trivial endomorphisms, we raise an error:

    sage: R.<x>=QQ[]
    sage: f = x^6 - 2*x^4 + 2*x^2 - 1
    sage: C = HyperellipticCurve(f)
    sage: A = C.jacobian()
    sage: rho = A.galois_representation()
    sage: rho.non_surjective()
    Traceback (most recent call last):
    ...
    NotImplementedError: Computation of non-surjective primes currently only works for Jacobians whose geometric endomorphism ring is Z.

REFERENCES:

- [Di2002]_

AUTHORS:

- Barinder S. Banwait, Armand Brumer, Hyun Jong Kim, Zev Klagsbrun,
  Jacob Mayle, Padmavathi Srinivasan, Isabel Vogt

"""

######################################################################
#                 Copyright (C) 2022 The Authors
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
from sage.arith.all import valuation, lcm, gcd
from sage.rings.fast_arith import prime_range
from sage.misc.lazy_import import lazy_import
from sage.modular.all import CuspForms
from sage.modular.dirichlet import DirichletGroup
from sage.misc.all import prod
from sage.rings.all import GF, ZZ, QQ, Zmod, PolynomialRing
from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve

from math import sqrt, floor


lazy_import('sage.interfaces.genus2reduction',
            ['genus2reduction', 'Genus2reduction'])


class GaloisRepresentation(SageObject):
    r"""
    The compatible family of Galois representations
    attached to the Jacobian of a  hyperelliptic curve over the rational
    numbers.

    EXAMPLES::

        sage: R.<x>=QQ[]
        sage: f = x**5 + 17
        sage: C = HyperellipticCurve(f)
        sage: J = C.jacobian()
        sage: rho = J.galois_representation()
        sage: rho
        Compatible family of Galois representations associated to the Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + 17

    """

    # Written by: TODO
    def __init__(self, A):
        r"""
        see ``GaloisRepresentation`` for documentation

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: f = x**5 + 17
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: loads(rho.dumps()) == rho
            False

        """
        self.__image_type = {}
        self._A = A
        self.non_surjective_primes = None

    def __repr__(self):
        r"""
        string representation of the class

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: f = x**5 + 17
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho
            Compatible family of Galois representations associated to the Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + 17

        """
        return ("Compatible family of Galois representations "
                f"associated to the {repr(self._A)}")

#####################################################################
# surjectivity
#####################################################################

    # Written by: TODO
    def _init_exps(self):
        r"""
        Return a dictionary of lists of characteristic polynomials of the
        matrices in the exceptional subgroup of `\operatorname{GSp}(4,\ell)`.

        OUTPUT:
        A dictionary whose keys are `\ell = 3,5, 7`; for each `\ell`,
        the associated value is the list of characteristic polynomials of the
        matrices in the exceptional subgroup of `\operatorname{GSp}(4,\ell)`.

        TESTS::

            sage: R.<x>=QQ[]
            sage: f = x**5 + 17
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: init_exps = rho._init_exps()
            sage: init_exps.keys()
            dict_keys([3, 5, 7])
            sage: assert isinstance(init_exps[3], list)
            sage: assert isinstance(init_exps[5], list)
            sage: assert isinstance(init_exps[7], list)
            sage: assert init_exps[3][0] in PolynomialRing(Zmod(3), "x")
            sage: assert init_exps[5][0] in PolynomialRing(Zmod(5), "x")
            sage: assert init_exps[7][0] in PolynomialRing(Zmod(7), "x")

        """
        # TODO consider keeping the dictionary as a class attribute rather
        # than having this as a function.
        # char3 is the list of characteristic polynomials of matrices in the
        # one subgroup of GSp(4,3) (up to conjugation) that isn't ruled out by
        # surj_tests
        R = PolynomialRing(Zmod(3), "x")
        x = R.gen()
        char3 = [
            x**4 + 2*x**3 + x**2 + 2*x + 1,
            x**4 + 1,
            x**4 + x**3 + 2*x**2 + x + 1,
            x**4 + 2*x**3 + 2*x + 1,
            x**4 + x**3 + 2*x**2 + 2*x + 1,
            x**4 + x**3 + x**2 + x + 1,
            x**4 + x**3 + x**2 + 2*x + 1,
            x**4 + 2*x**2 + 1,
            x**4 + x**3 + x + 1,
            x**4 + 2*x**3 + 2*x**2 + x + 1,
            x**4 + 2*x**3 + 2*x**2 + 2*x + 1,
            x**4 + x**2 + 1,
            x**4 + 2*x**3 + x**2 + x + 1]
        # char5 Is the list of characteristic polynomials of matrices in the
        # one subgroup of GSp(4,5) (up to conjugation) that isn't ruled out by
        # surj_tests
        R = PolynomialRing(Zmod(5), "x")
        x = R.gen()
        char5 = [
            x**4 + x**3 + 2*x**2 + x + 1,
            x**4 + x**3 + 4*x + 1,
            x**4 + x**3 + x**2 + x + 1,
            x**4 + x**3 + x + 1,
            x**4 + 4*x**2 + 1,
            x**4 + 4*x**3 + 4*x**2 + 4*x + 1,
            x**4 + 3*x**2 + 1,
            x**4 + 4*x**3 + 3*x**2 + 4*x + 1,
            x**4 + 2*x**2 + 1,
            x**4 + 4*x**3 + 4*x**2 + x + 1,
            x**4 + 4*x**3 + 2*x**2 + 4*x + 1,
            x**4 + x**2 + 1,
            x**4 + 2*x**3 + 4*x**2 + 3*x + 1,
            x**4 + 4*x**3 + 3*x**2 + x + 1,
            x**4 + 1,
            x**4 + 4*x**3 + x**2 + 4*x + 1,
            x**4 + 2*x**3 + 3*x**2 + 3*x + 1,
            x**4 + 4*x**3 + 2*x**2 + x + 1,
            x**4 + 4*x**3 + 4*x + 1,
            x**4 + 2*x**3 + 2*x**2 + 3*x + 1,
            x**4 + 4*x**3 + x**2 + x + 1,
            x**4 + 3*x**3 + 4*x**2 + 3*x + 1,
            x**4 + 2*x**3 + x**2 + 3*x + 1,
            x**4 + 4*x**3 + x + 1,
            x**4 + 3*x**3 + 3*x**2 + 3*x + 1,
            x**4 + 2*x**3 + 3*x + 1,
            x**4 + 3*x**3 + 2*x**2 + 3*x + 1,
            x**4 + 3*x**3 + x**2 + 3*x + 1,
            x**4 + 3*x**3 + 3*x + 1,
            x**4 + 2*x**3 + 4*x**2 + 2*x + 1,
            x**4 + 2*x**3 + 3*x**2 + 2*x + 1,
            x**4 + 2*x**3 + 2*x**2 + 2*x + 1,
            x**4 + 3*x**3 + 4*x**2 + 2*x + 1,
            x**4 + 2*x**3 + x**2 + 2*x + 1,
            x**4 + x**3 + 4*x**2 + 4*x + 1,
            x**4 + 3*x**3 + 3*x**2 + 2*x + 1,
            x**4 + 2*x**3 + 2*x + 1,
            x**4 + x**3 + 3*x**2 + 4*x + 1,
            x**4 + 3*x**3 + 2*x**2 + 2*x + 1,
            x**4 + x**3 + 4*x**2 + x + 1,
            x**4 + 3*x**3 + x**2 + 2*x + 1,
            x**4 + x**3 + 2*x**2 + 4*x + 1,
            x**4 + x**3 + 3*x**2 + x + 1,
            x**4 + x**3 + x**2 + 4*x + 1,
            x**4 + 3*x**3 + 2*x + 1]
        # char7 Is the list of characteristic polynomials of matrices in the
        # one subgroup of GSp(4,7) (up to conjugation) that isn't ruled out by
        # surj_tests
        R = PolynomialRing(Zmod(7), "x")
        x = R.gen()
        char7 = [
            x**4 + 2*x**3 + 5*x**2 + 5*x + 1,
            x**4 + 5*x**3 + 5*x**2 + 2*x + 1,
            x**4 + x**3 + 6*x**2 + 2*x + 4,
            x**4 + x**3 + 3*x**2 + 4*x + 2,
            x**4 + 5*x**3 + 3*x**2 + 5*x + 1,
            x**4 + 1,
            x**4 + 2*x**2 + 1,
            x**4 + 6*x**2 + 1,
            x**4 + 4*x**3 + x + 4,
            x**4 + 4*x**3 + 2*x**2 + x + 4,
            x**4 + 6*x**3 + x**2 + 6*x + 1,
            x**4 + 4*x**3 + 5*x**2 + 2*x + 2,
            x**4 + x**3 + x + 1,
            x**4 + 3*x**3 + 5*x**2 + 5*x + 2,
            x**4 + 3*x**3 + 6*x + 4,
            x**4 + 3*x**3 + 2*x**2 + 6*x + 4,
            x**4 + 6*x**3 + 3*x**2 + 3*x + 2,
            x**4 + 2,
            x**4 + x**3 + 4*x**2 + 6*x + 1,
            x**4 + 3*x**2 + 4,
            x**4 + 6*x**2 + 2,
            x**4 + 5*x**2 + 4,
            x**4 + 3*x**3 + 6*x**2 + 3*x + 1,
            x**4 + 6*x**3 + 6*x**2 + 5*x + 4,
            x**4 + 3*x**3 + 4*x**2 + 4*x + 1,
            x**4 + 4*x**3 + 4*x**2 + 3*x + 1,
            x**4 + 4*x**3 + 6*x**2 + 4*x + 1,
            x**4 + 2*x**3 + x + 2,
            x**4 + x**3 + 2*x**2 + 3*x + 2,
            x**4 + 2*x**3 + 4*x**2 + x + 2,
            x**4 + 6*x**3 + 4*x**2 + x + 1,
            x**4 + 3*x**3 + x**2 + x + 4,
            x**4 + 2*x**3 + x**2 + 3*x + 4,
            x**4 + x**3 + 3*x**2 + 5*x + 4,
            x**4 + 5*x**2 + 1,
            x**4 + 2*x**3 + 5*x**2 + 4*x + 4,
            x**4 + 3*x**3 + 6*x**2 + 2*x + 2,
            x**4 + 6*x**3 + 6*x + 1,
            x**4 + 2*x**3 + 2*x**2 + 6*x + 2,
            x**4 + x**3 + x**2 + x + 1,
            x**4 + 5*x**3 + 2*x**2 + x + 2,
            x**4 + 5*x**3 + 5*x**2 + 3*x + 4,
            x**4 + 4*x**3 + 6*x**2 + 5*x + 2,
            x**4 + 4*x**3 + x**2 + 6*x + 4,
            x**4 + 5*x**3 + x**2 + 4*x + 4,
            x**4 + 2*x**3 + 3*x**2 + 2*x + 1,
            x**4 + 6*x**3 + 3*x**2 + 2*x + 4,
            x**4 + x**2 + 2,
            x**4 + 4,
            x**4 + 5*x**3 + 6*x + 2,
            x**4 + 3*x**2 + 2,
            x**4 + 6*x**3 + 2*x**2 + 4*x + 2,
            x**4 + 4*x**2 + 4,
            x**4 + 5*x**3 + 4*x**2 + 6*x + 2]
        return {3: char3, 5: char5, 7: char7}

    # Written by: TODO
    def _init_wit(self, L):
        r"""
        Initialize a dictionary of lists for witnesses of surjectivitiy tests.

        INPUT:

        - ``L`` - list of primes;

        OUTPUT:

        A dictionary whose keys are the items `\ell` of ``L``. The
        corresponding values are lists of lengths 1, 2, or 3 whose items are
        all `0`. These lists are intended to eventually carry the data of
        results for surjectivity tests or witnessnes for surjectivity tests
        in the following key-value formats:

        - `\ell = 2`: [_is_surj_at_2]
        - `\ell = 3, 5, 7`: [witness for _surj_test_A,
                             witness for _surj_test_B,
                             witness for _surj_test_exp]
        - `\ell > 7`: [witness for surj_test_A, witness for _surj_test_B]

        TESTS::

            sage: R.<x>=QQ[]
            sage: f = x**5 + 17
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho._init_wit([2, 3, 5, 7])
            {2: [0], 3: [0, 0, 0], 5: [0, 0, 0], 7: [0, 0, 0]}
            sage: rho._init_wit([2, 3, 11])
            {2: [0], 3: [0, 0, 0], 11: [0, 0]}

        """
        witnesses = {}
        for l in L:
            if l == 2:
                witnesses[l] = [0]
            elif l in [3, 5, 7]:
                witnesses[l] = [0, 0, 0]
            else:
                witnesses[l] = [0, 0]
        return witnesses

    # Written by: TODO
    def _is_surj_at_2(self, f, h):
        """
        Return ``True`` if the mod-`2` Galois image of the Jacobian of the
        hyperelliptic curve specified by polynomials is surjective.

        The hyperelliptic curve is expected to be ``self._A.curve()``.
        The hyperelliptic curve is given by the equation `y^2 + h(x)y = f(x)`.

        INPUT:

        - ``f`` -- rational polynomial

        - ``h`` -- rational polynomial

        ALGORITHM:

        The following algorithm is adapted from

        The mod-`2` Galois representation of the Jacobian of the hyperelliptic
        curve given by `y^2 + h(x)y = f(x)` is surjective if and only if the
        Galois group of the polynomial `4f + h^2` is `S_6`.

        EXAMPLES:

            sage: R.<x>=QQ[]
            sage: f = x**5 + 17
            sage: C = HyperellipticCurve(f)
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho._is_surj_at_2(f, 0)
            False

        The following example is for the curve 1253.b.1253.1 on the LMFDB. ::

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 0, 1]), R([1, 1, 0, 1]));
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: f, h = C.hyperelliptic_polynomials()
            sage: rho._is_surj_at_2(f, h)
            True

        """
        # TODO: Add the paper to Sage's master bibliography file, and cite
        # the part that talks about the mod-2 Galois image in the ALGORITHM
        # section.
        F = 4*f + h**2
        return F.is_irreducible() and F.galois_group().order() == 720

    # Written by: TODO
    def _surj_test_A(self, frob_mod):
        r"""
        Return ``True`` if the specified mod-`\ell` characteristic polynomial
        of Frobenius is irreducible.

        INPUT:

        - ``frob_mod`` -- polynomial over `\GF{\ell}` for some prime `\ell`;
          this must be the mod-`\ell` reduction of the characteristic
          polynomial of the Frobenius element `\operatorname{Frob}_p` of
          `\operatorname{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})` acting on
          the `\ell`-adic Tate module `T_\ell A`, where `p \neq \ell` is
          a prime of good reduction for the abelian surface `A/\mathbb{Q}`.

        OUTPUT:
        ``True`` if ``frob_mod`` is irreducible over `\GF{\ell}`.
        ``False`` otherwise.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 0, 1]), R([1, 1, 0, 1]));
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho._surj_test_A(PolynomialRing(GF(5), 'x')(x^4 + x^3 + 3*x^2 + 3*x + 4))
            False
            sage: rho._surj_test_A(PolynomialRing(GF(7), 'x')(x^4 + x^3 + 3*x^2 + 3*x + 2))
            True
            sage: rho._surj_test_A(PolynomialRing(GF(3), 'x')(x^4 + x^3 + 2*x + 1))
            True
            sage: rho._surj_test_A(PolynomialRing(GF(5), 'x')(x^4 + x^3 + x^2 + x + 1))
            False
            sage: rho._surj_test_A(PolynomialRing(GF(5), 'x')(x^4 + x^3 + 2*x^2 + 3*x + 4))
            True
        """
        return frob_mod.is_irreducible()

    # Written by: TODO
    def _surj_test_B(self, frob_mod):
        r"""
        Return ``True`` if the specified mod-`\ell` characteristic polynomial
        of Frobenius has nonzero trace and has a linear factor with
        multiplicity one.

        INPUT:

        - ``frob_mod`` -- polynomial over `\GF{\ell}` for some prime `\ell`;
          this must be the mod-`\ell` reduction of the characteristic
          polynomial of the Frobenius element `\operatorname{Frob}_p` of
          `\operatorname{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})` acting on
          the `\ell`-adic Tate module `T_\ell A`, where `p \neq \ell` is
          a prime of good reduction for the abelian surface `A/\mathbb{Q}`.

        OUTPUT:
        ``True`` if ``frob_mod`` has nonzero trace over `\GF{\ell}` and
        has a linear factor with multiplicity one over `\GF{\ell}`.
        ``False`` otherwise.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 0, 1]), R([1, 1, 0, 1]));
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho._surj_test_B(PolynomialRing(GF(5), 'x')(x^4 + x^3 + 3*x^2 + 3*x + 4))
            False
            sage: rho._surj_test_B(PolynomialRing(GF(7), 'x')(x^4 + x^3 + 3*x^2 + 3*x + 2))
            False
            sage: rho._surj_test_B(PolynomialRing(GF(3), 'x')(x^4 + x^3 + 2*x + 1))
            False
            sage: rho._surj_test_B(PolynomialRing(GF(7), 'x')(x^4 + 4*x^3 + 2*x^2 + 6*x + 4))
            False
            sage: rho._surj_test_B(PolynomialRing(GF(3), 'x')(x^4 + x^3 + 2*x + 1))
            False
            sage: rho._surj_test_B(PolynomialRing(GF(5), 'x')(x^4 + x^3 + x^2 + x + 1))
            False

        The LMFDB curve 1091.a.1091.a has a non-surjective mod-`7` Galois
        representation as it has a rational torsion point of order
        `7`. ::

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, -1, -2, 0, 1]), R([1, 1, 1]));            
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()            
            sage: rho._surj_test_B(PolynomialRing(GF(7), 'x')(x^4 + 3*x^3 + 6*x^2 + 2*x + 2))
            True

        """
        if -frob_mod[3] != 0:
            for fact in frob_mod.factor():
                if fact[0].degree() == 1 and fact[1] == 1:
                    return True
        return False

    # Written by: TODO
    def _surj_test_exp(self, l, frob_mod, exps):
        r"""
        Return ``True`` if the specified mod-`\ell` characteristic polynomial
        of Frobenius is the characteristic polynomial of a matrix that is
        not in the exceptional subgroup of
        `\operatorname{GSp}_4(\mathbb{F}_\ell)`.

        INPUT:

        - ``l`` -- prime integer

        - ``frob_mod`` -- polynomial over `\GF{\ell}`;
          this must be the mod-`\ell` reduction of the characteristic
          polynomial of the Frobenius element `\operatorname{Frob}_p` of
          `\operatorname{Gal}(\overline{\mathbb{Q}}/\mathbb{Q})` acting on
          the `\ell`-adic Tate module `T_\ell A`, where `p \neq \ell` is
          a prime of good reduction for the abelian surface `A/\mathbb{Q}`.

        - ``exps`` -- dictionary; the keys are prime integers. The value
          corresponding to a key `\ell` is the list of characteristic
          polynomials of the exceptional subgroup of
          `\operatorname{GSp}_4(\mathbb{F}_\ell)`.

        OUTPUT:
        ``True`` if ``frob_mod`` is the characteristic polynomial of a matrix
        that is not in the exceptional subgroup of
        `\operatorname{GSp}_4(\mathbb{F}_\ell)`.
        ``False`` otherwise.

        EXAMPLES:

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([1, 4, 4, 2, 0, 0, 1]))
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: exps = rho._init_exps()
            sage: rho._surj_test_exp(5, PolynomialRing(GF(5), 'x')(x^4 + 3*x^3 + x^2 + 4*x + 4), exps)
            True
            sage: rho._surj_test_exp(7, PolynomialRing(GF(7), 'x')(x^4 + 3*x^3 + 6*x^2 + 2*x + 2), exps)
            False
            sage: rho._surj_test_exp(3, PolynomialRing(GF(3), 'x')(x^4 + 2*x^2 + 1), exps)
            False
            sage: rho._surj_test_exp(7, PolynomialRing(GF(7), 'x')(x^4 + 3*x^3 + 5*x^2 + x + 4), exps)
            True
            sage: rho._surj_test_exp(3, PolynomialRing(GF(3), 'x')(x^4 + 1), exps)
            False

        """
        # TODO: check that the inputs are specified correctly in the docstring
        return frob_mod not in exps[l]

    # Written by: TODO
    def _update_wit(self, l, p, frob, f, h, exps, wit):
        """
        Return an updated list of witnesses, based on surjectivity tests for
        ``frob`` at p.

        TESTS::

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([1, 4, 4, 2, 0, 0, 1]))
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: f, h = rho._A.curve().hyperelliptic_polynomials()
            sage: f, h
            (x^6 + 2*x^3 + 4*x^2 + 4*x + 1, 0)
            sage: L = prime_range(1000)
            sage: witnesses = rho._init_wit(L)
            sage: exps = rho._init_exps()
            sage: rho._update_wit(2, 5, x^4 + 2*x^2 + 25, f, h, exps, witnesses[2])
            [-1]
        """
        # print(f"_update_wit({l}, {p}, {frob}, {f}, {h}, {exps}, {wit})")
        frob_mod = frob.change_ring(Zmod(l))
        for i in range(0, len(wit)):
            if wit[i] == 0:
                if l == 2:
                    if self._is_surj_at_2(f,h):
                        wit[i] = 1
                    else:
                        wit[i] = -1
                elif i == 0 and self._surj_test_A(frob_mod):
                    wit[i] = p
                elif i == 1 and self._surj_test_B(frob_mod):
                    wit[i] = p
                elif i == 2 and self._surj_test_exp(l, frob_mod, exps):
                    wit[i] = p
        return wit

    # Written by: TODO
    def find_surj_from_list(
            self, L=prime_range(1000), bound=1000, verbose=False):
        r"""
        Return a list of primes `\ell` in ``L`` for which the mod-`\ell` Galois
        representation of the Jacobian of the hyperelliptic curve might
        not be surjective.

        INPUT:

        - ``L`` -- list of primes (default: list of primes less than 1000);
          the primes `\ell` to consider the surjectivity of the mod-`\ell`
          Galois representations for.

        - ``bound`` -- integer (default: `1000`); the exclusive upper bound
          for the primes `p \neq \ell` whose Frobenii `\operatorname{Frob}_p`
          are checked. Only the primes of good reduction which are at least `3`
          are checked.

        - ``verbose`` -- boolean (default: `False`)

        OUTPUT:

        A list of prime numbers. The mod-`\ell` Galois representations
        are provably surjective for all primes `\ell` which are in ``L``
        but not in the returned list.

        EXAMPLES:

        In the following example, the hyperelliptic curve is given by
        `y^2 + (x^3+1)y = x^2 + x`. The function returns the list ``[2, 7]``,
        which tells us that `2` and `7` are the only primes `\ell < 1000`
        whose mod-`\ell` Galois representation for the (Jacobian of the)
        hyperelliptic curve might not be surjective. ::

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([0, 1, 1]), R([1, 0, 0, 1]));
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.find_surj_from_list()  # long time
            [2, 7]

        The hyperelliptic curve is given by
        `y^2 + (x+1)y = x^5 + x^4 - 9x^3 - 5x^2 + 21x` in the following
        example. ::

            sage: R.<x> = PolynomialRing(QQ)
            sage: H = HyperellipticCurve(R([0, 21, -5, -9, 1, 1]), R([1, 1])); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x + 1)*y = x^5 + x^4 - 9*x^3 - 5*x^2 + 21*x
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.find_surj_from_list()  # long time
            [2, 13]

        """
        # TODO: Add examples with varying `L` and `bound``
        # TODO: add verbose examples in the docstring

        H = self._A.curve()
        f, h = H.hyperelliptic_polynomials()
        # C = 2 * genus2reduction(h, f).conductor
        # C is an integer which agrees up with the conductor of
        # H: y**2 + h y = f, # except possibly at 2.
        # Bad primes of Jac(H) divide it.
        C = 2 * prod(genus2reduction(h, f).local_data.keys())
        witnesses = self._init_wit(L)
        exps = self._init_exps()
        # to_check is the list of primes for which we still need to determine
        # surjectivity. Initially, it equals L and we remove primes as
        # their status is determined.
        to_check = L.copy()
        for p in prime_range(3, bound):
            if C % p != 0:
                Hp = H.change_ring(GF(p))
                frob = Hp.frobenius_polynomial()
                to_remove = []
                for l in to_check:
                    if p != l and 0 in witnesses[l]:
                        witnesses[l] = self._update_wit(
                            l, p, frob, f, h, exps, witnesses[l])
                        if 0 not in witnesses[l]:
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

    # Written by: TODO
    def is_surjective(self, l, bound=1000, verbose=False):
        r"""
        Return whether the mod-`\ell` representation is surjective onto
        `\operatorname{Aut}(A[\ell]) = \operatorname{GSp}_4(\GF{\ell})`.

        For the primes `\ell=2` and `3`, this function will always return
        either ``True`` or ``False``. For larger primes it might return
        ``None``.

        The output of this function is cached.

        INPUT:

        -  ``\ell`` -- prime integer

        - ``bound`` -- integer (default: `1000`); the exclusive upper bound
          for the primes `p \neq \ell` whose Frobenii `\operatorname{Frob}_p`
          are checked. Only the primes of good reduction which are at least `3`
          are checked.

        - ``verbose`` -- boolean (default: ``False``)

        OUTPUT: ``True`` if the mod-`\ell` representation is determined to be
        surjective, ``False`` if the representation is determined to be
        not surjection, and ``None`` if the representation is neither
        determined to be surjective nor determined to be not surjective.

        EXAMPLES::

            sage: R.<x>=QQ[]
            sage: f = x^6 + 2*x^3 + 4*x^2 + 4*x + 1
            sage: C = HyperellipticCurve(f)
            sage: A = C.jacobian()
            sage: rho = A.galois_representation()
            sage: rho.is_surjective(7)  # long time
            False
            sage: rho.is_surjective(2)
            False
            sage: rho.is_surjective(3)
            True
            sage: rho.is_surjective(13)
            True

        TESTS::

            sage: R.<x> = PolynomialRing(Rationals())
            sage: H = HyperellipticCurve(R([0, 21, -5, -9, 1, 1]), R([1, 1])); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x + 1)*y = x^5 + x^4 - 9*x^3 - 5*x^2 + 21*x
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.is_surjective(2)
            False
            sage: rho.is_surjective(3)
            True
            sage: rho.is_surjective(5)
            True
            sage: rho.is_surjective(7)
            True
            sage: rho.is_surjective(11)
            True
            sage: rho.is_surjective(13)  # long time
            False
            sage: rho.is_surjective(17)
            True

        ::

            sage: R.<x> = PolynomialRing(Rationals())
            sage: H = HyperellipticCurve(5*x^6-4*x^5+20*x^4-2*x^3+24*x^2+20*x+5)
            sage: H
            Hyperelliptic Curve over Rational Field defined by y^2 = 5*x^6 - 4*x^5 + 20*x^4 - 2*x^3 + 24*x^2 + 20*x + 5
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.is_surjective(2)
            False
            sage: rho.is_surjective(3)  # long time
            False
            sage: rho.is_surjective(5)  # long time
            False
            sage: rho.is_surjective(7)  # long time
            False
            sage: rho.is_surjective(11)  # long time
            False
            sage: rho.is_surjective(13)  # long time
            False

        ::

            sage: R.<x> = PolynomialRing(Rationals())  #249.a.249.1
            sage: H = HyperellipticCurve(R([2, 3, 1, 1, 0, -1]), R([1, 0, 0, 1])); H
            Hyperelliptic Curve over Rational Field defined by y^2 + (x^3 + 1)*y = -x^5 + x^3 + x^2 + 3*x + 2
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.is_surjective(2)
            False
            sage: rho.is_surjective(3)
            True
            sage: rho.is_surjective(5)
            True
            sage: rho.is_surjective(7)  # long time
            False

        """
        # TODO: the two tests are the first two examples in the
        # find_surj_from_list.sage file. Maybe the actual outputs
        # should be checked.
        # TODO: add hyperelliptic curve examples
        # TODO: add verbose examples
        # TODO: Add the papers of Dokchitsers and Elklies in the master
        # bibliography file; consider adding the reference to these papers
        # in the module docstring instead.
        if self.non_surjective_primes is not None:
            if not l.is_prime():
                raise ValueError("l must be prime")
            return (l not in self.non_surjective_primes)

        ans = self.find_surj_from_list(L=[l], bound=1000, verbose=False)

        if ans:
            return False
        else:
            return True

    # Written by: TODO
    def non_surjective(self, N=None, bound=1000):
        r"""
        Return a list of primes `\ell` such that the mod-`\ell` representation
        is probably not surjective. If `\ell` is not in the returned list,
        then the mod-p representation is provably surjective.

        By a theorem of Serre, there are only finitely many primes in this
        list, except when the jacobian of the curve has nontrivial geometric
        endomorphisms.

        INPUT:

        - ``N`` -- an integer (default: None); this is the conductor of the
        jacobian; providing this can speed up the computation.

        - ``bound`` -- integer (default: `1000`); the exclusive upper bound
          for the primes `p \neq \ell` whose Frobenii `\operatorname{Frob}_p`
          are checked. Only the primes of good reduction which are at least `3`
          are checked.

        OUTPUT:

        A list; if the curve has nontrivial geometric endomorphisms,
        returns ``[0]``. Otherwise, returns a list of primes `\ell` where
        mod-`\ell` representation is very likely not surjective.
        The mod-`\ell` representation is definitely surjective for any prime
        `\ell` not in this list.

        EXAMPLES::

            sage: R.<x> = PolynomialRing(QQ);
            sage: H = HyperellipticCurve(R([0, 1, 1]), R([1, 0, 0, 1]));
            sage: J = H.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.non_surjective()  # long time
            [2, 7]

        The curve 743.a.743.1 on the LMFDB has non-surjective mod-`\ell` representations
        for all primes $\ell$. ::

            sage: R.<x> = PolynomialRing(QQ); C = HyperellipticCurve(R([0, 0, 1, 0, -1]), R([1, 1, 0, 1]));
            sage: J = C.jacobian()
            sage: rho = J.galois_representation()
            sage: rho.non_surjective()
            []

        If the Jacobian has any non-trivial endomorphisms, we raise an error. ::

            sage: R.<x>=QQ[]
            sage: f = x^6 - 2*x^4 + 2*x^2 - 1
            sage: C = HyperellipticCurve(f)
            sage: A = C.jacobian()
            sage: rho = A.galois_representation()
            sage: rho.non_surjective()
            Traceback (most recent call last):
            ...
            NotImplementedError: Computation of non-surjective primes currently only works for Jacobians whose geometric endomorphism ring is Z.

        ALGORITHM:

        We first find an upper bound `B` on the possible primes. If `E`
        is semi-stable, we can take `B=11` by a result of Mazur. There is
        a bound by Serre in the case that the `j`-invariant is not integral
        in terms of the smallest prime of good reduction. Finally
        there is an unconditional bound by Cojocaru, but which depends
        on the conductor of `E`.
        For the prime below that bound we call ``is_surjective``.

        """
        # TODO Add the results of Mazur, Serre, and Cojocaru in Sage's master
        # bibliography file and cite them here.
        # TODO Add an explanation of why the primes are likely to not be
        # surjective.
        # TODO: add the attribute non_surjective_primes to the class level
        # docstring
        if self.non_surjective_primes is not None:
            return self.non_surjective_primes

        C = self._A.curve()
        f, h = C.hyperelliptic_polynomials()
        simplified_hyperelliptic_model = HyperellipticCurve(4*f + h**2)
        jacobian = simplified_hyperelliptic_model.jacobian()
        if not jacobian.geometric_endomorphism_ring_is_ZZ():
            raise NotImplementedError(
                "Computation of non-surjective primes currently only works "
                "for Jacobians whose geometric endomorphism ring is Z."
            )

        M1p3 = 0
        y1p3 = 0
        M2p2nsd = 0
        y2p2nsd = 0

        if N is None:
            f, h = C.hyperelliptic_polynomials()
            red_data = genus2reduction(h, f)
            N = red_data.conductor  # is this the true conductor if red_data.prime_to_2_conductor_only is False?
            max_cond_exp_2 = None
            if red_data.prime_to_2_conductor_only:
                # I think this is the case where we don't know exactly the
                # two-part of conductor
                N = 2*N
                max_cond_exp_2 = red_data.minimal_disc.valuation(2)

        # MCusp is a list of the form <S,M,y>, where S is either a space
        # of cusp forms or a level, M is an integer such that all primes with
        # a reducible sub isomorphic to the rep of a cusp form in S divide M,
        # y is a counter for the number of nontrivial Frobenius conditions
        # go into M
        MCusp = _set_up_cuspidal_spaces(N, max_cond_exp_2=max_cond_exp_2)

        # MQuad is a list of the form <phi,M,y>, where phi is a
        # quadratic character, M is the integer all nonsurjective primes
        # governed by phi must divide, and y is counter for the number of
        # nontrivial Frobenius conditions going into M
        MQuad = _set_up_quadratic_chars(N)

        d = _maximal_square_divisor(N)

        # we'll test as many p as we need to get at least 2 nontrivial
        # Frobenius conditions for every possible cause of non-surjectivity
        sufficient_p = False

        p = 1

        while not sufficient_p:
            p = ZZ(p).next_prime()
            if N % p != 0:
                Cp = C.change_ring(GF(p))
                fp = Cp.frobenius_polynomial()
                fp_rev = Cp.zeta_function().numerator()

                f = Zmod(d)(p).multiplicative_order()
                c = gcd(f, 120)
                c = lcm(c, 8)  # adding in the max power of 2
                tp = - fp.coefficients(sparse=False)[3]
                sp = fp.coefficients(sparse=False)[2]

                M1p3, y1p3 = _rule_out_1_plus_3_via_Frob_p(
                    c, p, tp, sp, M1p3, y1p3)
                M2p2nsd, y2p2nsd = _rule_out_2_plus_2_nonselfdual_via_Frob_p(
                    c, p, tp, sp, M2p2nsd, y2p2nsd)
                MCusp = _rule_out_cuspidal_spaces_using_Frob_p(
                    p, fp_rev, MCusp)
                MQuad = _rule_out_quadratic_ell_via_Frob_p(p, fp, MQuad)

            if (M1p3 == 1) or (y1p3 > 1):
                if (M2p2nsd == 1) or (y2p2nsd > 1):
                    if all((Mc == 1 or yc>1) for S, Mc, yc in MCusp):
                        if all((Mq == 1 or yq > 1) for phi, Mq, yq in MQuad):
                            sufficient_p = True

        # we will always include 2, 3, 5, 7 and the non-semistable primes.
        non_maximal_primes = {2, 3, 5, 7}.union(
            set([p[0] for p in list(N.factor()) if p[1] > 1]))

        ell_red_easy = [M1p3.prime_factors(), M2p2nsd.prime_factors()]
        non_maximal_primes = non_maximal_primes.union(
            set([p for j in ell_red_easy for p in j]))

        ell_red_cusp = [(S.level(), ZZ(M).prime_factors())
                        for S, M, y in MCusp]

        non_maximal_primes = non_maximal_primes.union(
            set([p for a, j in ell_red_cusp for p in j]))

        ell_irred = [(phi, ZZ(M).prime_factors()) for phi, M, t in MQuad]
        non_maximal_primes = non_maximal_primes.union(
            set([p for a, j in ell_irred for p in j]))
        self.non_surjective_primes = self.find_surj_from_list(
            L=non_maximal_primes, bound=bound)
        return self.non_surjective_primes

#########################################################
#                            #
#                   Auxiliary functions            #
#                            #
#########################################################


# Written by: TODO
def _maximal_square_divisor(N):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _maximal_square_divisor
        sage: _maximal_square_divisor(1)
        1
        sage: _maximal_square_divisor(2)
        1
        sage: _maximal_square_divisor(30)
        1
        sage: _maximal_square_divisor(54)
        3
        sage: _maximal_square_divisor(2176)
        8
        sage: _maximal_square_divisor(252)
        6
        sage: _maximal_square_divisor(0)
        Traceback (most recent call last):
        ...
        ArithmeticError: factorization of 0 is not defined

    """
    PP = ZZ(N).prime_divisors()
    n = 1
    for p in PP:
        n = n * p**(floor(valuation(N, p)/2))
    return n

#########################################################
#                            #
#          Governed by a quadratic character        #
#                            #
#########################################################


# Written by: TODO
def _maximal_quadratic_conductor(N):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _maximal_quadratic_conductor
        sage: _maximal_quadratic_conductor(1)
        1
        sage: _maximal_quadratic_conductor(2)
        8
        sage: _maximal_quadratic_conductor(3)
        3
        sage: _maximal_quadratic_conductor(15)
        15
        sage: _maximal_quadratic_conductor(30)
        120
        sage: _maximal_quadratic_conductor(72)
        24
        sage: _maximal_quadratic_conductor(343)
        7
        sage: _maximal_quadratic_conductor(495)
        165
    """
    if N % 2 == 0:
        return 4 * ZZ(N).radical()
    else:
        return ZZ(N).radical()


# Written by: TODO
def _character_list(N):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _character_list
        sage: _character_list(1)
        []
        sage: _character_list(2)
        [Dirichlet character modulo 8 of conductor 4 mapping 7 |--> -1, 5 |--> 1,
        Dirichlet character modulo 8 of conductor 8 mapping 7 |--> 1, 5 |--> -1,
        Dirichlet character modulo 8 of conductor 8 mapping 7 |--> -1, 5 |--> -1]
        sage: _character_list(3)
        [Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1]
        sage: _character_list(15)
        [Dirichlet character modulo 15 of conductor 3 mapping 11 |--> -1, 7 |--> 1,
        Dirichlet character modulo 15 of conductor 5 mapping 11 |--> 1, 7 |--> -1,
        Dirichlet character modulo 15 of conductor 15 mapping 11 |--> -1, 7 |--> -1]
        sage: len(_character_list(30))
        15

    """
    c = _maximal_quadratic_conductor(N)
    D = DirichletGroup(c, base_ring=QQ, zeta_order=2)
    return [phi for phi in D if phi.conductor() != 1]


# Written by: TODO
def _set_up_quadratic_chars(N):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _set_up_quadratic_chars
        sage: assert len(_set_up_quadratic_chars(1)) == 0
        sage: assert len(_set_up_quadratic_chars(2)) == 3
        sage: assert len(_set_up_quadratic_chars(3)) == 1
        sage: assert len(_set_up_quadratic_chars(15)) == 3
        sage: assert len(_set_up_quadratic_chars(30)) == 15
        sage: _set_up_quadratic_chars(3)[0]
        (Dirichlet character modulo 3 of conductor 3 mapping 2 |--> -1, 0, 0)

    """
    return [(phi, 0, 0) for phi in _character_list(N)]


# Written by: TODO
def _rule_out_quadratic_ell_via_Frob_p(p, fp, MM):
    r"""Provide a summary of what this method is doing.

    INPUT:

    - ``p`` -- prime integer; new prime.

    - ``fp`` -- polynomial over the integers; the characteristic polynomial
      of Frobenius at ``p`` on a hyperelliptic curve.

    - ``MM`` -- list; the items are tuples of the form ``(\phi, M, y)``,
      where ``\phi`` is a non-trivial quadratic character, all primes
      ``\ell`` for which there is a quadratic obstruction associated with
      ``\phi`` must divide ``M``, ``y`` is a counter for the number of
      nontrivial Frobenius constraints going into ``M``.

    OUTPUT: a list

    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _rule_out_quadratic_ell_via_Frob_p, _set_up_quadratic_chars
        sage: R.<x> = PolynomialRing(QQ)
        sage: MQuad = _set_up_quadratic_chars(498)
        sage: len(MQuad)
        15
        sage: assert isinstance(MQuad[0], tuple)
        sage: assert len(MQuad[0]) == 3
        sage: assert MQuad[0][1] == 0
        sage: assert MQuad[0][2] == 0
        sage: MQuad[0][0]
        Dirichlet character modulo 1992 of conductor 4 mapping 1495 |--> -1, 997 |--> 1, 665 |--> 1, 1081 |--> 1
        sage: MQuad = _rule_out_quadratic_ell_via_Frob_p(5, x^4 + 2*x^2 + 25, MQuad)
        sage: len(MQuad)
        15
        sage: MQuad = _rule_out_quadratic_ell_via_Frob_p(7, x^4 + x^3 - 2*x^2 + 7*x + 49, MQuad)
        sage: len(MQuad)
        15


    """
    # TODO: Complete the description and output description.
    # print(f"_rule_out_quadratic_ell_via_Frob_p({p}, {fp}, See above for input MM)")
    ap = -fp.coefficients(sparse=False)[3]
    if ap == 0:
        return MM
    else:
        MM0 = []
        for phi, M, y in MM:
            if (M == 1 or phi(p) != -1 or y > 1):
                MM0.append((phi, M, y))
            else:
                MM0.append((phi, gcd(M, p*ap), y+1))
        return MM0


#########################################################
#                            #
#             Reducible (easy cases)            #
#                            #
#########################################################

# Isabel's code using the fact that the exponent of the determinant
# character divides 120.

# the following three functions implement the following for n=2,3,5:
# f(x) = x**4 - t*x**3 + s*x**2 - p**alpha*t*x + p**(2*alpha) is a
# degree 4 polynomial whose roots multiply in pairs to p**alpha
# returns the tuple (p, tn, sn, alphan) of the polynomial
# f**(n)(x) = x**4 - tn*x**3 + sn*x**2 - p**(alphan)*tn*x + p**(2*alphan)
# whose roots are the nth powers of the roots of f

# Written by: TODO
def _power_roots2(ptsa):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _power_roots2
        sage: _power_roots2((5, 0, 2, 1))
        (5, -4, 54, 2)
        sage: _power_roots2((7, -1, -2, 1))
        (7, 5, 88, 2)
        sage: _power_roots2((11, 1, 2, 1))
        (11, -3, 224, 2)
        sage: _power_roots2((43, -604, 6075814, 4))
        (43, -11786812, 57797449706566, 8)

    """
    p, t, s, alpha = ptsa
    return (p, t**2 - 2*s, s**2 - 2*p**alpha*t**2 + 2*p**(2*alpha),
            2*alpha)


# Written by: TODO
def _power_roots3(ptsa):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _power_roots3
        sage: _power_roots3((5, 0, 2, 1))
        (5, 0, -142, 3)
        sage: _power_roots3((7, -1, -2, 1))
        (7, -28, 622, 3)
        sage: _power_roots3((11, 1, 2, 1))
        (11, 28, -58, 3)
        sage: _power_roots3((43, -604, 6075814, 4))
        (43, 4594158692, 14096196939501130726, 12)

    """
    p, t, s, alpha = ptsa
    return (
        p,
        t**3 - 3*s*t + 3*p**alpha*t,
        s**3 - 3*p**alpha*s*t**2 + 3*p**(2*alpha)*t**2
        + 3*p**(2*alpha)*t**2 - 3*p**(2*alpha)*s,
        3*alpha)


# Written by: TODO
def _power_roots5(ptsa):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _power_roots5
        sage: _power_roots5((5, 0, 2, 1))
        (5, 0, 5282, 5)
        sage: _power_roots5((7, -1, -2, 1))
        (7, 109, 3678, 5)
        sage: _power_roots5((11, 1, 2, 1))
        (11, -649, 267002, 5)
        sage: _power_roots5((43, -604, 6075814, 4))
        (43, -10608093658315484, -452809576478300004206054711642586, 20)

    """
    p, t, s, alpha = ptsa
    return (
        p,
        t**5 - 5*s*t**3 + 5*s**2*t + 5*p**alpha*t**3 - 5*p**alpha*s*t -
        5*p**(2*alpha)*t, s**5 - 5*p**alpha*s**3*t**2
        + 5*p**(2*alpha)*s*t**4 + 5*p**(2*alpha)*s**2*t**2
        - 5*p**(3*alpha)*t**4 + 5*p**(2*alpha)*s**2*t**2
        - 5*p**(2*alpha)*s**3 - 5*p**(3*alpha)*t**4
        - 5*p**(3*alpha)*s*t**2 + 5*p**(4*alpha)*t**2
        + 5*p**(4*alpha)*t**2 + 5*p**(4*alpha)*s,
        5*alpha)


# put these together to do any power dividing 120 that we actually need
# c is the power
# Written by: TODO
def _power_roots(cptsa):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _power_roots
        sage: _power_roots((30, 5, 0, -142, 3))
        (5,
        -40439319944794656777824091572164,
        2024421783291998908109426253166728620498898914410821047074441974,
        90)
        sage: _power_roots((11, 3, 5, -27, 11))
        Traceback (most recent call last):
        ...
        ValueError: can't raise to this power

    """
    c, p, t, s, alpha = cptsa
    if 120 % c != 0:
        raise ValueError("can't raise to this power")

    ptsa = (p, t, s, alpha)

    while c % 2 == 0:
        c, ptsa = c/2, _power_roots2(ptsa)

    while c % 3 == 0:
        c, ptsa = c/3, _power_roots3(ptsa)

    while c % 5 == 0:
        c, ptsa = c/5, _power_roots5(ptsa)

    return ptsa


# given a quartic f whose roots multiply to p**alpha in pairs,
# returns the quartic whose roots are the products of roots
# of f that DO NOT multiply to p**alpha
# Written by: TODO
def _roots_pairs_not_p(ptsa):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _roots_pairs_not_p
        sage: _roots_pairs_not_p((5, 0, 2, 1))
        (5, -8, 30, 2)
        sage: _roots_pairs_not_p((7, -1, -2, 1))
        (7, -16, 133, 2)
        sage: _roots_pairs_not_p((11, 1, 2, 1))
        (11, -20, 209, 2)
        sage: _roots_pairs_not_p((43, -604, 6075814, 4))
        (43, 6075728, -506829218, 8)

    """
    p, t, s, alpha = ptsa
    return (p, s - 2*p, p*t**2 - 2*p*s + 2*p**2, 2*alpha)


# t and s are the first and second elementary symmetric functions in the
# roots of the characteristic polynomial of Frobenius at p for a curve C
# M is an integer such that every prime ell for which J_C[ell] could be
# 1 + 3 reducible divides M
# y is a counter for the number of nontrivial Frobenius conditions going
# into M
# Written by: TODO
def _rule_out_1_plus_3_via_Frob_p(c, p, t, s, M=0, y=0):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _rule_out_1_plus_3_via_Frob_p
        sage: _rule_out_1_plus_3_via_Frob_p(8, 5, 0, 2, 0, 0)
        (759564288000, 1)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 43, 8, 70, 7168, 11)
        (7168, 12)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 41, 4, -2, 7168, 10)
        (7168, 11)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 31, 3, 30, 7168, 8)
        (7168, 9)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 29, 3, 32, 7168, 7)
        (7168, 8)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 23, 0, -26, 7168, 6)
        (7168, 7)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 19, -6, 22, 7168, 5)
        (7168, 6)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 17, -1, 28, 7168, 4)
        (7168, 5)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 13, 0, -2, 179200, 3)
        (7168, 4)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 11, 1, 2, 1254400, 2)
        (179200, 3)
        sage: _rule_out_1_plus_3_via_Frob_p(8, 7, -1, -2, 759564288000, 1)
        (1254400, 2)
    """
    # print(f"_rule_out_1_plus_3_via_Frob_p({c}, {p}, {t}, {s}, {M}, {y})")
    p, tnew, snew, alphanew = _power_roots((c, p, t, s, 1))
    x = PolynomialRing(QQ, "x").gen()
    Pnew = (x**4 - tnew*x**3 + snew*x**2 - p**alphanew*tnew*x
            + p**(2*alphanew))
    Pval = Pnew(1)
    if Pval != 0:
        return ZZ(gcd(M, p*Pval)), y+1
    else:
        return M, y


# t and s are the first and second elementary symmetric functions in the
# roots of the characteristic polynomial of Frobenius at p for a curve C
# M is an integer such that every prime ell for which J_C[ell] could be
# 2+2 non-self-dual type reducible divides M
# y is a counter for the number of nontrivial Frobenius conditions going
# into M
# Written by: TODO
def _rule_out_2_plus_2_nonselfdual_via_Frob_p(c, p, t, s, M=0, y=0):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _rule_out_2_plus_2_nonselfdual_via_Frob_p
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 5, 0, 2, 0, 0)
        (0, 0)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 43, 8, 70, 1, 8)
        (1, 9)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 37, -5, 8, 1, 6)
        (1, 7)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 29, 3, 32, 1, 4)
        (1, 5)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 19, -6, 22, 1, 3)
        (1, 4)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 17, -1, 28, 7, 2)
        (1, 3)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 13, 0, -2, 7, 2)
        (7, 2)
        sage: _rule_out_2_plus_2_nonselfdual_via_Frob_p(8, 11, 1, 2, 54877666934581886389580413794274402125175221246505984375, 1)
        (7, 2)
    """
    # print(f"_rule_out_2_plus_2_nonselfdual_via_Frob_p({c}, {p}, {t}, {s}, {M}, {y})")
    p, tnew, snew, alphanew = _roots_pairs_not_p((p, t, s, 1))
    p, tnew, snew, alphanew = _power_roots(
        (c, p, tnew, snew, alphanew))
    x = PolynomialRing(QQ, "x").gen()
    Pnew = (x**4 - tnew*x**3 + snew*x**2 - p**alphanew*tnew*x
            + p**(2*alphanew))
    Pval = Pnew(1)*Pnew(p**c)
    if Pval != 0:
        return ZZ(gcd(M, p*Pval)), y+1
    else:
        return M, y

#########################################################
#                            #
#          Reducible (modular forms case)        #
#                            #
#########################################################


# Written by: TODO
def _special_divisors(N):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _special_divisors
        sage: _special_divisors(1)
        [1]
        sage: _special_divisors(6)
        [2]
        sage: _special_divisors(15)
        [3]
        sage: _special_divisors(8)
        [2]
        sage: _special_divisors(30)
        [2, 3, 5]
        sage: _special_divisors(57)
        [3]

    """
    D0 = [d for d in ZZ(N).divisors() if d <= sqrt(N)]
    D0.reverse()
    D = []
    for d0 in D0:
        if all([d % d0 != 0 for d in D]):
            D.append(d0)
    D.reverse()
    return D


# Written by: TODO
def _get_cuspidal_levels(N, max_cond_exp_2=None):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _get_cuspidal_levels
        sage: _get_cuspidal_levels(498, 0)
        [3]

    """
    # print(f"_get_cuspidal_levels({N}, {max_cond_exp_2}")
    if max_cond_exp_2 is not None:
        # if we're here, then N is the even poor mans conductor
        # recall we put a 2 in the poor mans conductor
        conductor_away_two = N/2
        possible_conductors = [conductor_away_two * 2 ** i
                                for i in range(max_cond_exp_2 + 1)]
        # not ordered, hopefully not a problem.
        return list(set([d for N in possible_conductors
                        for d in _special_divisors(N)]))
    else:
        return _special_divisors(N)


# Written by: TODO
def _set_up_cuspidal_spaces(N, max_cond_exp_2=None):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _set_up_cuspidal_spaces
        sage: _set_up_cuspidal_spaces(498, 0)
        [(Cuspidal subspace of dimension 0 of Modular Forms space of dimension 1 for Congruence Subgroup Gamma0(3) of weight 2 over Rational Field,
          0,
          0)]

    """
    # print(f"_set_up_cuspidal_spaces({N}, {max_cond_exp_2})")
    D = _get_cuspidal_levels(N, max_cond_exp_2)
    return [(CuspForms(d), 0, 0) for d in D]


# Written by: TODO
def _reconstruct_hecke_poly_from_trace_polynomial(cusp_form_space, p):
    """
    Implement Zev and Joe Wetherell's idea

    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _reconstruct_hecke_poly_from_trace_polynomial
        sage: cusp_form_space = CuspForms(Gamma0(3), 2)
        sage: _reconstruct_hecke_poly_from_trace_polynomial(cusp_form_space, 5)
        1
        sage: _reconstruct_hecke_poly_from_trace_polynomial(cusp_form_space, 7)
        1
    """
    # print(f"_reconstruct_hecke_poly_from_trace_polynomial(See above for input cusp_form_space, {p})")
    R = PolynomialRing(QQ, "x")
    x = R.gen()
    char_T_x = R(cusp_form_space.hecke_polynomial(p))
    S = PolynomialRing(QQ, 2, "ab")
    a, b = S.gens()
    char_T_a_b = S(char_T_x(x=a)).homogenize(var='b')
    substitute_poly = char_T_a_b(a=1+p*b**2)
    return R(substitute_poly(a=0, b=x))


# Written by: TODO
def _rule_out_cuspidal_space_using_Frob_p(S, p, fp, M, y):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _rule_out_cuspidal_space_using_Frob_p
        sage: S = CuspForms(Gamma0(3), 2)
        sage: R.<x> = PolynomialRing(QQ)
        sage: _rule_out_cuspidal_space_using_Frob_p(S, 5, 25*x^4 + 2*x^2 + 1, 0, 0)
        (5, 1)
        sage: _rule_out_cuspidal_space_using_Frob_p(S, 7, 49*x^4 + 7*x^3 - 2*x^2 + x + 1, 5, 1)
        (1, 2)

    """
    # print(f"_rule_out_cuspidal_space_using_Frob_p(See above for input S, {p}, {fp}, {M}, {y})")
    if M != 1 and y < 2:
        Tp = _reconstruct_hecke_poly_from_trace_polynomial(S, p)
        res = fp.resultant(Tp)
        if res != 0:
            return gcd(M, p*res), y+1
    return M, y


# Written by: TODO
def _rule_out_cuspidal_spaces_using_Frob_p(p, fp, MC):
    """
    TESTS::

        sage: from sage.schemes.hyperelliptic_curves.gal_reps import _set_up_cuspidal_spaces, _rule_out_cuspidal_spaces_using_Frob_p
        sage: MCusp = _set_up_cuspidal_spaces(498, max_cond_exp_2=0)
        sage: MCusp
        [(Cuspidal subspace of dimension 0 of Modular Forms space of dimension 1 for Congruence Subgroup Gamma0(3) of weight 2 over Rational Field,
          0,
          0)]
        sage: R.<x> = PolynomialRing(QQ)
        sage: MCusp = _rule_out_cuspidal_spaces_using_Frob_p(5, 25*x^4 + 2*x^2 + 1, MCusp)
        sage: MCusp
        [(Cuspidal subspace of dimension 0 of Modular Forms space of dimension 1 for Congruence Subgroup Gamma0(3) of weight 2 over Rational Field,
          5,
          1)]
        sage: MCusp = _rule_out_cuspidal_spaces_using_Frob_p(7, 49*x^4 + 7*x^3 - 2*x^2 + x + 1, MCusp)
        sage: MCusp
        [(Cuspidal subspace of dimension 0 of Modular Forms space of dimension 1 for Congruence Subgroup Gamma0(3) of weight 2 over Rational Field,
          1,
          2)]

    """
    # print(f"_rule_out_cuspidal_spaces_using_Frob_p({p}, {fp}, See above for input MC)")
    MC0 = []
    for S, M, y in MC:
        Mm, yy = _rule_out_cuspidal_space_using_Frob_p(S, p, fp, M, y)
        MC0.append((S, Mm, yy))
    return MC0
