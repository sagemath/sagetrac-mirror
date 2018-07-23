r"""
The OM (Okutsu-Montes/Ore-Mac Lane) tree of an monic and square free polynomial over a fixed modulus pi-adic ring.

AUTHORS:

- Brian Sinclair and Sebastian Pauli (2012-02-22): initial version (as factoring.py)
- Brian Sinclair and Sebastian Pauli (2017-07-18): class OMTree

"""
from sage.rings.polynomial.padics.omtree.frame import Frame
from sage.structure.factorization import Factorization
from sage.structure.sage_object import SageObject

class OMTree(SageObject):
    """
    The OM (Okutsu-Montes/Ore-Mac Lane) tree of an monic and square free polynomial over a fixed modulus pi-adic ring.

    REFERENCES:

    [GMN12] J. Guàrdia, J. Montes, E. Nart, Newton polygons of higher order in algebraic number theory, Transactions of the AMS, Volume 364, Number 1, 2012.

    [Sin15] B. Sinclair,  Algorithms for Enumerating Invariants and Extensions of Local Fields, PhD thesis, UNCG 2015, http://www.uncg.edu/mat/numbertheory/publications/Sinclair_2015_Algorithms-for-Enumerating-Invariants-and-Extensions-of-Local-Fields.pdf

    [Pau10] S. Pauli, Factoring polynomials over local fields II. Algorithmic Number Theory, 9th International Symposium, ANTS-IX, Nancy, France, July 2010, LNCS 6197, 301-315, Springer Verlag 2010.

    AUTHORS:

    - Brian Sinclair and Sebastian Pauli (2012-02-22): initial version
    - Brian Sinclair and Sebastian Pauli (2017-07-18): class OMTree

    """
    def __init__(self, Phi):
        r"""
        Construct a tree of OM (Okutsu-Montes/Ore-Mac Lane) representations for Phi.

        INPUT:

        - ``Phi`` -- squarefree, monic padic polynomial with fixed precision
          coefficients

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: Phi = ZpFM(2, 20, 'terse')['x'](x^32 + 16)
            sage: OMTree(Phi)
            OM Tree of (1 + O(2^20))*x^32 + (0 + O(2^20))*x^31 + (0 + O(2^20))*x^30 + (0 + O(2^20))*x^29 + (0 + O(2^20))*x^28 + (0 + O(2^20))*x^27 + (0 + O(2^20))*x^26 + (0 + O(2^20))*x^25 + (0 + O(2^20))*x^24 + (0 + O(2^20))*x^23 + (0 + O(2^20))*x^22 + (0 + O(2^20))*x^21 + (0 + O(2^20))*x^20 + (0 + O(2^20))*x^19 + (0 + O(2^20))*x^18 + (0 + O(2^20))*x^17 + (0 + O(2^20))*x^16 + (0 + O(2^20))*x^15 + (0 + O(2^20))*x^14 + (0 + O(2^20))*x^13 + (0 + O(2^20))*x^12 + (0 + O(2^20))*x^11 + (0 + O(2^20))*x^10 + (0 + O(2^20))*x^9 + (0 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (16 + O(2^20)) with 1 leaves
        """
        from sage.misc.flatten import flatten

        self._Phi=Phi

        def followsegment(next, Phi):
            """
            Returns next if it corresponds to an irreducible factor of $\Phi$ 
            and follows every branch if not.

            """
            # Handle the unlikely event that our approximation is actually a factor
            if next.is_root() == False and next.phi == next.prev_frame().phi:
                return [next]
            if next.phi_divides_Phi():
                return [next] + [[followsegment(fact.next_frame(fact.multiplicity + 1), Phi)
                                for fact in seg.factors] for seg in next.polygon[1:]]
            # Check if we have an irreducible factor
            if sum([seg.length for seg in next.polygon]) == 1:
                return next
            return [[followsegment(fact.next_frame(fact.multiplicity + 1), Phi)
                     for fact in seg.factors] for seg in next.polygon]

        # Construct and initialize the first frame (phi = x)
        next = Frame(Phi)
        next.seed(Phi.parent().gen())

        # Handle the special case wherein our initial approximation (phi = x) is a factor
        if next.phi_divides_Phi():
            tree = [next] + [[followsegment(fact.next_frame(fact.multiplicity + 1), Phi)
                              for fact in seg.factors] for seg in next.polygon[1:]]

        # Construct the next level of the tree by following every factor of the
        # residual polynomials of every Newton polygon segment in our frame
        else:
            tree = [[followsegment(fact.next_frame(fact.multiplicity + 1), Phi)
                     for fact in seg.factors] for seg in next.polygon]

        # tree contains the leaves of the tree of frames where each leaf corresponds
        # to an irreducible factor of Phi.  We flatten the list to remove unneeded nesting.
        self._leaves = flatten(tree)

    def Phi(self):
        """
        The polynomial ``Phi`` underlying the OM Tree.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: Phi = ZpFM(2, 20, 'terse')['x'](x^8 + 2)
        sage: T = OMTree(Phi)
        sage: T.Phi()
        (1 + O(2^20))*x^8 + (0 + O(2^20))*x^7 + (0 + O(2^20))*x^6 + (0 + O(2^20))*x^5 + (0 + O(2^20))*x^4 + (0 + O(2^20))*x^3 + (0 + O(2^20))*x^2 + (0 + O(2^20))*x + (2 + O(2^20))
        """
        return self._Phi

    def leaves(self):
        """
        The leaves of the OM tree.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<x> = ZpFM(3, 7)[]
        sage: f=x^4 + 234
        sage: om=OMTree(f)
        sage: om.leaves()
        [Frame with phi (1 + O(3^7))*x^2 + (O(3^7))*x + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7)),
         Frame with phi (1 + O(3^7))*x^2 + (O(3^7))*x + (2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))]
        """
        return self._leaves

    def ramification_indices(self):
        """
        The ramification indices of the extensions generated by the irreducible factors of the polynomial ``Phi`` underlying the OMTree self.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<c> = ZqFM(125, prec = 30)
        sage: Rz.<z>=R[]
        sage: g=(z^3 + 2)^5 + 5
        sage: om=OMTree(g)
        sage: om.ramification_indices()
        [5, 5]
        """
        return [fr.E for fr in self.leaves()]

    def inertia_degrees(self):
        """
        The inertia degrees of the extensions generated by the irreducible factors of the polynomial ``Phi`` underlying the OMTree self.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<c> = ZqFM(125, prec = 30)
        sage: Rz.<z>=R[]
        sage: g=(z^3 + 2)^5 + 5
        sage: om=OMTree(g)
        sage: om.inertia_degrees()
        [1, 2]
        """
        return [fr.F for fr in self.leaves()]

    def degrees_of_factors(self):
        """
        The degrees of the irreducible factors of the polynomial ``Phi`` underlying the OMTree.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<c> = ZqFM(125, prec = 30)
        sage: Rz.<z>=R[]
        sage: g=(z^6 + 2)^25 + 5
        sage: om=OMTree(g)
        sage: om.degrees_of_factors()
        [50, 50, 50]

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R=ZpFM(5, 50)
        sage: S.<x>=R[]
        sage: f=x^5 + 75*x^3-15*x^2 + 125*x-5
        sage: W.<w>=R.ext(f)
        sage: Wz.<z>=W[]
        sage: g=(z^5 + 75*z^3-15*z^2 + 125*z-5)*(z^25 + 5)*(z^5 + w) + 125
        sage: om=OMTree(g)
        sage: om.degrees_of_factors()
        [1, 2, 2, 5, 25]
        """
        return [fr.phi.degree() for fr in self.leaves()]

    def approximations(self):
        """
        Liftable/Montes approximations to the irreducible factors of the polynomial ``Phi`` underlying the OMTree.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<x> = ZpFM(3, 7)[]
        sage: f=x^4 + 234
        sage: om=OMTree(f)
        sage: om.approximations()
        [(1 + O(3^7))*x^2 + (O(3^7))*x + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7)), (1 + O(3^7))*x^2 + (O(3^7))*x + (2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))]
        """
        return [fr.phi for fr in self.leaves()]

    def roots(self):
        """
        The roots of the trees whose leaves are in the OMTree.

        As a note, the leaves for a polynomial's OM Tree may have different roots.
        Becuase optimal OM Trees remove 'improvement frames' (those for which
        neither ramificaiton or inertia is found), an OM Tree may technically
        be a forest.

        See :meth: `Frame.root()` which is used to find the root for each leaf.

        EXAMPLES::

        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
        sage: R.<c> = ZqFM(125, prec = 30, print_mode='terse')
        sage: Rz.<z>=R[]
        sage: g=(z^3 + 2)^5 + 5
        sage: om=OMTree(g)
        sage: om.roots()
        [Frame with phi (1 + O(5^30))*z + (931322574615478515623 + O(5^30)),
         Frame with phi (1 + O(5^30))*z + (0 + O(5^30))]
        """
        return [fr.root() for fr in self.leaves()]

    def depths(self):
        """
        The depths of the leaves of the OMTree

        EXAMPLES::
        sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree

        sage: k = ZpFM(2, prec = 20)
        sage: kx.<x> = k[]
        sage: f = x^32 + 16
        sage: om=OMTree(f)
        sage: om.depths()
        [3]

        sage: R.<c> = ZqFM(125, prec = 30)
        sage: Rz.<z>=R[]
        sage: g=(z^3 + 2)^5 + 5
        sage: om=OMTree(g)
        sage: om.depths()
        [1, 2]
        """
        return [fr.depth for fr in self.leaves()]

    def factorization(self):
        r"""
        Return a factorization of Phi using an OM (Okutsu-Montes/Ore-Mac Lane) algorithm.

        This is accomplished by constructing a tree of Frames of approximations
        of factors.

        The factorization of the polynomial ``Phi`` underlying the OM Tree.

        EXAMPLES:

        Factoring polynomials over Zp(2)[x]::

           sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
           sage: f = ZpFM(2, 50, 'terse')['x']( (x^32 + 16)*(x^32 + 16 + 2^16*x^2) + 2^34 )
           sage: omtree = OMTree(f); factors = omtree.factorization(); len(factors) # long time (4.5s)
           2

        See the irreducibility of x^32 + 16 in Zp(2)[x]::

           sage: OMTree(ZpFM(2)['x'](x^32 + 16)).factorization()
           (1 + O(2^20))*x^32 + (O(2^20))*x^31 + (O(2^20))*x^30 + (O(2^20))*x^29 + (O(2^20))*x^28 + (O(2^20))*x^27 + (O(2^20))*x^26 + (O(2^20))*x^25 + (O(2^20))*x^24 + (O(2^20))*x^23 + (O(2^20))*x^22 + (O(2^20))*x^21 + (O(2^20))*x^20 + (O(2^20))*x^19 + (O(2^20))*x^18 + (O(2^20))*x^17 + (O(2^20))*x^16 + (O(2^20))*x^15 + (O(2^20))*x^14 + (O(2^20))*x^13 + (O(2^20))*x^12 + (O(2^20))*x^11 + (O(2^20))*x^10 + (O(2^20))*x^9 + (O(2^20))*x^8 + (O(2^20))*x^7 + (O(2^20))*x^6 + (O(2^20))*x^5 + (O(2^20))*x^4 + (O(2^20))*x^3 + (O(2^20))*x^2 + (O(2^20))*x + (2^4 + O(2^20))

        We look at the appoximations to the factors and lift them::

            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(3, 7)[]
            sage: f=x^4 + 234
            sage: om=OMTree(f)
            sage: om.approximations()
            [(1 + O(3^7))*x^2 + (O(3^7))*x + (3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7)), (1 + O(3^7))*x^2 + (O(3^7))*x + (2*3 + 2*3^2 + 2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))]
            sage: om.factorization()
            ((1 + O(3^7))*x^2 + (O(3^7))*x + (3 + 3^4 + 3^5 + 3^6 + O(3^7))) * ((1 + O(3^7))*x^2 + (O(3^7))*x + (2*3 + 2*3^2 + 2*3^3 + 3^4 + 3^5 + 3^6 + O(3^7)))

        AUTHORS:

        - Brian Sinclair and Sebastian Pauli (2012-02-22): initial version
        - Brian Sinclair and Sebastian Pauli (2017-07-18): class OMTree
        """

        Phi = self.Phi()

        # Handle the situation that x is a factor of $\Phi(x)$
        if Phi.constant_coefficient() == 0:
            x_divides = True
            Phi = Phi >> 1
        else:
            x_divides = False

        if Phi == Phi.base_ring().one():
            return Factorization([(Phi.parent().gen(), 1)] if x_divides else [])

        # the leaves correspond to the factors of Phi
        tree = self.leaves()

        # If we only have one leaf, Phi is irreducible, so we do not lift it.
        if len(tree) == 1:
            return Factorization(([(Phi.parent().gen(), 1)] if x_divides else []) +
                                 [(Phi, 1)])
        # quo_rem is faster than single_factor_lift, so Phi = f*g is specially handled
        if len(tree) == 2:
            if tree[0].phi.degree() < tree[1].phi.degree():
                fact = tree[0].single_factor_lift()
                return Factorization(([(Phi.parent().gen(), 1)] if x_divides else []) +
                                     [(fact, 1), (Phi.quo_rem(fact)[0], 1)])
            else:
                fact = tree[1].single_factor_lift()
                return Factorization(([(Phi.parent().gen(), 1)] if x_divides else []) +
                                     [(fact, 1), (Phi.quo_rem(fact)[0], 1)])
        # Phi has more than two factors, so we lift them all
        return Factorization(([(Phi.parent().gen(), 1)] if x_divides else []) +
                             [(frame.single_factor_lift(), 1) for frame in tree])

    def __repr__(self):
        """
        Representation of the OM tree self.

        EXAMPLES:
            sage: from sage.rings.polynomial.padics.omtree.omtree import OMTree
            sage: R.<x> = ZpFM(3, 7)[];f=(x^2 - 15)*(x^2-12) + 81
            sage: om = OMTree(f)
            sage: om
            OM Tree of (1 + O(3^7))*x^4 + (O(3^7))*x^3 + (2*3^3 + 2*3^4 + 2*3^5 + 2*3^6 + O(3^7))*x^2 + (O(3^7))*x + (2*3^2 + 3^5 + O(3^7)) with 2 leaves
        """
        return "OM Tree of " + str(self.Phi()) + " with " + str(len(self.leaves())) + " leaves"
