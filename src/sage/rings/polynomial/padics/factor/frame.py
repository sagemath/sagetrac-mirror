r"""
Factoring local field polynomials using an OM algorithm

AUTHORS:

- Brian Sinclair (2012-02-22): initial version

"""
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.padics.factor.frameelt import FrameElt, FrameEltTerm
from sage.rings.polynomial.padics.factor.segment import Segment
from sage.rings.infinity import infinity
from sage.misc.functional import denominator

class Frame:
    """
    All data for one iteration of local field polynomial factorization.

    A Frame object corresponds to a single node of the tree of OM
    representations. This will, when seeded with an approximation,
    contain a newton polygon of higher order (as a list of Segment
    objects), which will in turn contain associated polynomials and
    their factors (as lists of AssociatedFactors). Each of these
    branchings is a partitioning of the factors of ``Phi`` and will thus
    correspond to a branching of the OM tree.

    INPUT:

    - ``Phi`` - The local field polynomial being factored.

    - ``prev`` - AssociatedFactor, default None; The associated polynomial
      factor of this Frame's parent in the tree of OM representations
      that created this Frame.

    - ``iteration_count`` - Integer, default 0; The number of previous
      iterations of the algorithm run to this point. This will count steps
      of all types, including those discarded from the tree.

    OUTPUT:

    - A Frame representing a node in the OM tree of ``Phi`` with parent
      Frame containing ``prev``.

    EXAMPLES::

    Creating a Frame leaves it unseeded, not having an approximation::

        sage: from sage.rings.polynomial.padics.factor.frame import *
        sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
        sage: f = Frame(Phi)
        sage: f
        Unseeded Frame regarding (1 + O(2^20))*x^32 + (16 + O(2^20))

    Each Frame needs to be seeded to compute intermediate values::

        sage: f.seed(Phi.parent().gen())
        sage: f
        Frame with phi (1 + O(2^20))*x

    A Seeded Frame has a newton polygon as a list of Segment objects::

        sage: f.polygon
        [Segment of length 32 and slope 1/8]
        sage: f.polygon[0].psi.polynomial()
        2 + O(2^20)
        sage: f.polygon[0].factors
        [AssociatedFactor of rho z + 1]

    New Frames in the tree can be created from each AssociatedFactor::

        sage: f.polygon[0].factors[0].next_frame()
        Frame with phi (1 + O(2^20))*x^8 + (1048574 + O(2^20))

    Notice that Frames created by AssociatedFactors are seeded.    

    """
    def __init__(self,Phi,prev=None,iteration_count=0):
        """
        Initializes self.

        See ``Frame`` for full documentation.

        """
        self.prev = prev
        self.Phi = Phi
        self.O = Phi.base_ring()
        self.Ox = Phi.parent()
        self.x = self.Ox.gen()
        self.R = self.O.residue_class_field()
        self.Rz = PolynomialRing(self.R,names='z')
        self.phi = None
        self.iteration = iteration_count + 1
        if self.is_first(): # that is self.prev is None
            self.E = 1
            self.F = 1
            self.depth = 0
        else:
            self.E = self.prev_frame().E * self.prev.segment.Eplus
            self.F = self.prev_frame().F * self.prev.Fplus
            self.depth = self.prev_frame().depth + 1

    def seed(self,phi,length=infinity):
        """
        Seed all of the intermediate values of the Frame based on the new
        approximation to a factor ``phi``.

        In seeding, the Frame will compute a ``phi``-expansion of ``Phi`` and
        the newton polygon of higher order. Then Segments will construct their
        assocated polynomials and factors.

        INPUT:

        - ``phi`` - an approximation to a factor of ``Phi``

        - ``length`` - Integer or infinity; default infinity. the length of
          the segment that resulted in this frame.  This provides an extra
          stopping point for finding phi-expansions if finite as we only
          need up to ``length`` terms.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
            sage: f = Frame(Phi); f
            Unseeded Frame regarding (1 + O(2^20))*x^32 + (16 + O(2^20))
            sage: f.phi is None
            True
            sage: f.seed(Phi.parent().gen()); f
            Frame with phi (1 + O(2^20))*x
            sage: f.phi
            (1 + O(2^20))*x

        """
        self.phi = phi

        # Construct the phi expansion of Phi
        self._phi_expansion_as_polynomials = []
        if self.phi.degree() > self.Phi.degree():
            self._phi_expansion_as_polynomials = [self.Phi]
        q, r = self.Phi.quo_rem(self.phi)
        self._phi_expansion_as_polynomials.append(r)
        while q != 0 and length > len(self._phi_expansion_as_polynomials):
            q, r = q.quo_rem(self.phi)
            self._phi_expansion_as_polynomials.append(r)
        self._phi_expansion_as_elts = [FrameElt(self,a) for a in self._phi_expansion_as_polynomials]

        # If phi divides Phi, we may be in a leaf that resembles its parent
        # and need to break recursion
        if self.is_first() == False and self.phi == self.prev_frame().phi and self.phi_divides_Phi():
            return
        else:
            self.polygon = self._newton_polygon([e.valuation() for e in self._phi_expansion_as_elts]) # list of segments

    def find_psi(self,val):
        """
        Find a polynomial (as a FrameElt) with given valuation

        INPUT:

        - ``val`` - Rational. The desired valuation. The denominator of ``val``
          must divide the current level of ramification (``E``).

        OUTPUT:

        - A FrameElt with respect to the current frame with valuation ``val``.

        EXAMPLES::

        First we need an appropriate Frame::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^16+16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f
            Frame with phi (1 + O(2^20))*x^4 + (1048574 + O(2^20))

        We get a valid FrameElt with integer exponents as long as the
        denominator of ``val`` divides the current ramification::

            sage: f.E
            4
            sage: f.prev.segment.slope
            1/4
            sage: f.find_psi(7/4)
            [[1*2^1]phi1^3]
            sage: f.find_psi(7/4).polynomial()
            (2 + O(2^20))*x^3

        If the denominator does not divide the ramification, then we cannot
        construct a polynomial of this valuation and an error is raised::

            sage: f.find_psi(3/8)
            Traceback (most recent call last):
            ...
            ValueError: Denominator of given valuation does not divide E

        """

        if not self.E % denominator(val) == 0:
            raise ValueError, "Denominator of given valuation does not divide E"
        psielt = FrameElt(self)
        if self.prev is None:
            psielt.terms = [FrameEltTerm(psielt,self.O(1),val)]
        else:
            vphi = self.prev.segment.slope
            d = self.prev_frame().E
            vprime = val*d
            e = vphi * d
            psimod = denominator(e)
            s = 0
            if not psimod == 1:
                s = vprime / e
                if denominator(s) == 1:
                    s = s % psimod
                else:
                    s = int(s % psimod)
            val = val - s * vphi
            psielt.terms = [FrameEltTerm(psielt,self.prev_frame().find_psi(val),s)]
        return psielt

    def _newton_polygon(self,a,xoffset=0):
        """
        Compute the newton polygon of higher order of ``Phi`` with respect to
        the valuations of the phi-expansion of ``Phi``.

        This method constructs finds the slopes of each segment as well as
        the points on the line of each segment. These values are used to
        initialize the output Segments.

        INPUT:

        - ``a`` -- A List of the valuations of the coefficients in the 
          phi-expansion of ``Phi``.

        - ``xoffset`` -- Integer, default 0; The first x value to consider
          while constructing the newton polygon. The area to the left of
          ``xoffset`` will be ignored. This is mainly used for recursion.

        OUTPUT:

        - A List of Segments comprising the newton polygon of higher order.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x*(x+1)*(x+2))
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: v = [e.valuation() for e in f._phi_expansion_as_elts];v
            [20, 1, 0, 0]
            sage: f._newton_polygon(v)
            [Segment of length 1 and slope +Infinity,
             Segment of length 1 and slope 1,
             Segment of length 1 and slope 0]

        """
        if len(a) == 0:
            raise ValueError, "Cannot compute Newton polygon from empty list"

        # Handle the case where the first segment has infinite slope
        # (This will occur iff phi divides Phi)
        if self.phi_divides_Phi() and xoffset == 0:
            for i in range(1,len(a)):
                if a[i] < self.O.precision_cap():
                    verts = [(0,infinity),(i,a[i])]
                    slope = infinity
                    length = i
                    if i == len(a)-1:
                        return [Segment(self,slope,verts,length)]
                    else:
                        return [Segment(self,slope,verts,length)]+self._newton_polygon(a[verts[len(verts)-1][0]-xoffset:],verts[len(verts)-1][0])

        slope = (a[0]-a[len(a)-1]) / (len(a)-1)
        verts = [(xoffset,a[0])]
        length = 0
        for i in range(1,len(a)):
            y = a[0] - i*slope
            if a[i] == y:
                verts.append((xoffset+i,y))
                length = i
            elif a[i] < y:
                verts = [(xoffset,a[0]),(xoffset+i,a[i])]
                slope = (a[0]-a[i]) / i
                length = i
            elif y < a[len(a)-1]:
                if len(a[verts[len(verts)-1][0]-xoffset:]) == 0:
                    return [Segment(self,slope,verts,length)]
                else:
                    return [Segment(self,slope,verts,length)]+self._newton_polygon(a[verts[len(verts)-1][0]-xoffset:],verts[len(verts)-1][0])
        return [Segment(self,slope,verts,length)]

    # Data Access Methods

    def prev_frame(self):
        """
        Returns the parent of this frame in the OM tree.

        OUTPUT:

        - The Frame object that is the parent of self in the OM tree.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f
            Frame with phi (1 + O(2^20))*x
            sage: len(f.polygon)           
            1
            sage: len(f.polygon[0].factors)
            1
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f
            Frame with phi (1 + O(2^20))*x^8 + (1048574 + O(2^20))
            sage: f.prev_frame()
            Frame with phi (1 + O(2^20))*x

        """
        if self.prev is None:
            return None
        else:
            return self.prev.segment.frame      

    def is_first(self):
        """
        Returns ``True`` if self is the root of the OM tree, otherwise
        ``False``.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
            sage: f = Frame(Phi)
            sage: f.is_first()
            True
            sage: f.seed(Phi.parent().gen())
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f.is_first()
            False

        """
        return self.prev is None

    def phi_divides_Phi(self):
        """
        Returns ``True`` if phi divides Phi, otherwise ``False``.

        This is done by checking if the constant term of the phi-expansion
        of Phi is 0.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^3+x)          
            sage: f = Frame(Phi)
            sage: f.seed(Phi.parent().gen())
            sage: f.Phi
            (1 + O(2^20))*x^3 + (1 + O(2^20))*x
            sage: f.phi
            (1 + O(2^20))*x
            sage: f.phi_divides_Phi()
            True

        """
        return self._phi_expansion_as_polynomials[0] == 0

    def __repr__(self):
        """
        Representation of self.
        
        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import *
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)        
            sage: f = Frame(Phi)
            sage: f.__repr__()
            'Unseeded Frame regarding (1 + O(2^20))*x^32 + (16 + O(2^20))'
            sage: f.seed(Phi.parent().gen())
            sage: f.__repr__()
            'Frame with phi (1 + O(2^20))*x'

        """
        if self.phi:
            return 'Frame with phi '+repr(self.phi)
        else:
            return 'Unseeded Frame regarding '+repr(self.Phi)
