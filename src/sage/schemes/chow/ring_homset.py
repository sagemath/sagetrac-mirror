r"""
The set of morphism between ChowRings

Space of ChowRing homomorphisms.

Derives from the space of quotient ring homomorphisms implement in
`sage.ring.homset.RingHomset_quo_ring` in order to explicitly allow
"no generators", e.g. im_gens = [].


EXAMPLES::

    sage: A.<h> = ChowRing('h', 1, 'h^3')
    sage: B.<k> = ChowRing('k', 1, 'k^6')
    sage: phi = B.hom([2*h], A); phi
    Ring morphism:
      From: Quotient of Multivariate Polynomial Ring in k over Rational Field by the ideal (k^6)
      To:   Quotient of Multivariate Polynomial Ring in h over Rational Field by the ideal (h^3)
      Defn: k |--> 2*h
    sage: phi(2)
    2
    sage: phi(k)
    2*h

    sage: A = ChowRing()
    sage: f = A.hom([], A)
    sage: f(1)
    1
    sage: B.<h> = ChowRing('h', 1, 'h^3')
    sage: f = B.hom([0], A)
    sage: f(h)
    0
    sage: f(3)
    3

TESTS::

    sage: A = ChowRing()
    sage: B.<h> = ChowRing('h', 1, 'h^3')
    sage: H = B.Hom(A)
    sage: H == loads(dumps(H))
    True

    sage: TestSuite(phi).run()

AUTHORS:

- Manfred Lehn (2013)
- Christoph Sorger (2013)
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2013 Manfred Lehn <lehn@mathematik.uni-mainz.de>
#       Copyright (C) 2013 Christoph Sorger <christoph.sorger@univ-nantes.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.homset import RingHomset_generic
from sage.rings import morphism


class ChowRingHomSet(RingHomset_generic):

    def __call__(self, im_gens, check=True):
        r"""
        Build the homomorphism.

        TESTS::

            sage: A = ChowRing()
            sage: f = A.hom([], A)
            sage: f(1)
            1

        """
        if isinstance(im_gens, morphism.RingHomomorphism_from_quotient):
            return morphism.RingHomomorphism_from_quotient(self, im_gens._phi())
        try:
            pi = self.domain().cover()
            #
            # Remark: we explicitly allow "no generators", e.g. im_gens = [].
            # which is not the case in the original RingHomset_quo_ring
            # implementation.
            #
            if not im_gens:
                phi = pi.domain().hom([], self.codomain())
            else:
                if isinstance(im_gens, str):
                    im_gens = [im_gens]
                im_gens = [self.codomain()(x) for x in im_gens]
                phi = pi.domain().hom(im_gens, check=check)
            return morphism.RingHomomorphism_from_quotient(self, phi)
        except (NotImplementedError, ValueError):
            try:
                return self._coerce_impl(im_gens)
            except TypeError:
                raise TypeError("images do not define a valid homomorphism")

    def _coerce_impl(self, x):
        """
        Overrides _coerce_impl.
        """
        if not isinstance(x, morphism.RingHomomorphism_from_quotient):
            raise TypeError
        if x.parent() is self:
            return x
        if x.parent() == self:
            return morphism.RingHomomorphism_from_quotient(self, x._phi())
        raise TypeError
