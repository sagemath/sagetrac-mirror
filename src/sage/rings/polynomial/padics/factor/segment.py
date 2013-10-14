r"""
Segment of a newton polygon of higher order needed for OM computations.

AUTHORS:

- Brian Sinclair (2012-02-22): initial version

"""
from sage.rings.polynomial.padics.factor.associatedfactor import AssociatedFactor
from sage.rings.arith import gcd
from sage.rings.infinity import infinity
from sage.misc.functional import denominator

class Segment:
    r"""
    Segment of a newton polygon of higher order needed for OM computations.

    Segments are typically contained in a list with other Segments comprising
    the newton polygon of a frame. The Segment needs to be capable of refering
    back to its orignating frame for information about the polygon.  Further,
    most of the important data is computed while building the polygon, and is
    thus passed to the Segment at initialization rather than have it compute
    these for itself, namely: slope, points on the line, and horizontal length.

    Each Segment represents a partitioning of the roots of the original
    polynomial, and thus each is responsible for branching the OM tree.

    Also, Segments can find their own associated polynomials, whose factors
    represent further branching of the OM tree.

    INPUT:

    - ``frame`` -- Frame; the Frame to whose newton polygon this segment
      belongs.

    - ``slope`` -- Rational or infinity; The slope of the segment.

    - ``verts`` -- List of tuples; The list of vertices of points of the
      associated polynomial found on this segment. Most notably, this
      needs to include the endpoints of the segment.

    - ``length`` -- Integer; The horizontal length of the segment.

    EXAMPLES::

    Polygons will only have one segment if they cannot show reducibility::

        sage: from sage.rings.polynomial.padics.factor.frame import Frame
        sage: from sage.rings.polynomial.padics.factor.segment import Segment
        sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)
        sage: f = Frame(Phi); f.seed(Phi.parent().gen())
        sage: f.polygon
        [Segment of length 32 and slope 1/8]
        sage: f = f.polygon[0].factors[0].next_frame()
        sage: f.polygon                                 
        [Segment of length 4 and slope 5/4]
        sage: f = f.polygon[0].factors[0].next_frame()  
        sage: f.polygon                               
        [Segment of length 4 and slope 21/16]
        sage: f = f.polygon[0].factors[0].next_frame()
        sage: f.polygon                               
        [Segment of length 2 and slope 85/32]

    Segments note ramification increases (as ``Eplus``) and the valuation
    of the approximation (as ``slope``)::

        sage: f = Frame(Phi); f.seed(Phi.parent().gen())
        sage: f.polygon[0].Eplus
        8
        sage: f.polygon[0].slope
        1/8
        sage: f = f.polygon[0].factors[0].next_frame()  
        sage: f.polygon[0].Eplus                      
        1
        sage: f = f.polygon[0].factors[0].next_frame()
        sage: f.polygon[0].Eplus                      
        2
        sage: f = f.polygon[0].factors[0].next_frame()
        sage: f.polygon[0].Eplus                      
        2

    Associate polynomials for each segment check for possible inertia and
    the residue field is extended by their irreducible factors::

        sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
        sage: Phi = x^4+20*x^3+44*x^2+80*x+1040   
        sage: f = Frame(Phi); f.seed(Phi.parent().gen())
        sage: f.polygon
        [Segment of length 4 and slope 1]
        sage: f.polygon[0].associate_polynomial()
        z^4 + z^2 + 1
        sage: f.polygon[0].factors
        [AssociatedFactor of rho z^2 + z + 1]
        sage: f = f.polygon[0].factors[0].next_frame()
        sage: f.polygon
        [Segment of length 2 and slope 5]
        sage: f.polygon[0].associate_polynomial()
        z0^2 + a0*z0 + 1
        sage: f.polygon[0].factors
        [AssociatedFactor of rho z0^2 + a0*z0 + 1]

    """
    def __init__(self,frame,slope,verts,length):
        """
        Initialises self.

        See ``Segment`` for full documentation.

        """
        self.frame = frame
        self.verts = verts
        self.slope = slope
        self.length = length
        if slope != infinity:
            self.Eplus = (denominator(self.slope) /
                          gcd(denominator(self.slope),int(self.frame.E)))
            self.psi = self.frame.find_psi(self.slope*self.Eplus)
        else:
            self.Eplus = 1
        self._associate_polynomial = self.associate_polynomial(cached=False)
        self.factors = [AssociatedFactor(self,afact[0],afact[1]) 
                        for afact in list(self._associate_polynomial.factor())]

    def associate_polynomial(self,cached=True):
        """
        Return the associated polynomial of this segment.

        The associated polynomial is found by taking the points on the
        segment, shortening the segment by the discovered ramification
        (the increase in slope denominator) and taking their residues.

        INPUT:

        - ``cached`` - Boolean; default TRUE. If TRUE, returns the cached associated
          polynomial. If FALSE, complutes the associated polynomial anew.

        OUTPUT:

        - The associated polynomial of the segment, a polynomial over the
          residue field, which may have been extended in previous Frames.

        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: from sage.rings.polynomial.padics.factor.segment import Segment
            sage: Phi = ZpFM(2,20,'terse')['x'](x^32+16)    
            sage: f = Frame(Phi); f.seed(Phi.parent().gen())
            sage: seg = Segment(f,1/8,[(0,4),(32,0)],32); seg
            Segment of length 32 and slope 1/8
            sage: seg.associate_polynomial()
            z^4 + 1

        """
        if cached:
            return self._associate_polynomial

        if self.slope == infinity:
            if self.frame.prev == None:
                Az = self.frame.Rz.gen() ** self.length
            else:
                Az = self.frame.prev.FFz.gen() ** self.length
            return Az

        a = self.frame._phi_expansion_as_elts
        vertx = [v[0] for v in self.verts]
        chiex = [int((v-vertx[0]) // self.Eplus) for v in vertx]
        chi = [a[vertx[i]] * self.psi**chiex[i] for i in range(len(vertx))]
        psitilde = self.frame.find_psi(chi[0].valuation())
        Ahat = [(c/psitilde).reduce() for c in chi]
        if self.frame.prev == None:
            Az = sum([(Ahat[i].residue())*self.frame.Rz.gen()**chiex[i]
                       for i in range(len(Ahat))])
        else:
            Az = sum([(Ahat[i].residue())*self.frame.prev.FFz.gen()**chiex[i]
                       for i in range(len(Ahat))])
        return Az

    def __repr__(self):
        """
        Representation of self.
        
        EXAMPLES::

            sage: from sage.rings.polynomial.padics.factor.frame import Frame
            sage: k = ZpFM(2,20,'terse'); kx.<x> = k[]
            sage: Phi = x^4+20*x^3+44*x^2+80*x+1040
            sage: f = Frame(Phi); f.seed(Phi.parent().gen())
            sage: f.polygon[0].__repr__()
            'Segment of length 4 and slope 1'
            sage: f = f.polygon[0].factors[0].next_frame()
            sage: f.polygon[0].__repr__()
            'Segment of length 2 and slope 5'

        """
        return 'Segment of length '+repr(self.length)+' and slope '+repr(self.slope)
