#*****************************************************************************
#      Copyright (C) 2006 - 2017
#      Mark Bauer <bauerm@ucalgary.ca>
#      Avi Kulkarni <akulkarn@sfu.ca> <trac username: akulkarn>
#      Colin Weir <Colin.Weir@cse-cst.gc.ca>
#
#      Special thanks to Julian Ruth
#
# Code developed at Sage days 86.5
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRing_generic
from sage.rings.polynomial.polynomial_quotient_ring import PolynomialQuotientRingElement
from sage.rings.finite_rings.finite_field_base import *
from sage.categories.fields import Fields

# Element class
class RelativeFiniteFieldElement(PolynomialQuotientRingElement):
    def minpoly(self, var):
        """
        The minimal polynomial of self over the base field of the parent.

        INPUT:
        
        - 'var' -- a string. The name of the variable of the minimal polynomial.

        """
        return self.matrix().minpoly(var)

    
    def coefficients(self):
        """
        Returns the list of coefficients defining this element. Alias for list().
        """
        return self.list()

class RelativeFiniteField(PolynomialQuotientRing_generic):
    """
    Relative finite field class of some kind. What do we want this to do?
    Relative finite field class
    init:
    There are two constructors for RelativeFiniteField. The first creates an extension of Fq by an irreducible polynomial over Fq.
    The second returns an extension of a preexisting RelativeFiniteField L by a polynomial irreducible over L.

    class variables: 

    ground_field   -- Returns Fq. i.e, the first FiniteField object after succesive calls to base_ring()
    base_field     -- Returns the base_ring. The base_ring is either the ground field or the RelativeFiniteField from which self is constructed.
    base degree    -- The degree of self over its ground field
    GF             -- Constructs the FiniteField object isomorphic to the field corresponding to self. This variable is crucial to the genral operation of this class.
    GF_map         -- A map to ``self.GF()``
    GF_map_inverse -- The inverse to ``self.__GF``
    base_GF        -- The finite field isomorphic to the base field of ``self``
    """
    def __init__(self, ring, polynomial, name=None, category=None):
        """
        Initialize ``self``.

        Extensions of RelativeFiniteFields should be constructed using the ``extension()`` method.

        EXAMPLES::

        sage: k = GF(49);
        sage: R.<x> = PolynomialRing(k);
        sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
        sage: L = RelativeFiniteField(R,g,'y'); L
        Relative Finite Field with defining polynomial x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3 over Finite Field in z2 of size 7^2

        """
        assert ring.base_ring().is_field(), 'Base ring not a field'
        assert ring.base_ring().is_finite(), 'Base ring not finite'
        assert polynomial.parent() is ring, 'Parent of polynomial is not given ring'

        from sage.rings.polynomial.polynomial_element import Polynomial
        if not isinstance(polynomial, Polynomial):
            try:
                polynomial = polynomial.polynomial(self)
            except (AttributeError, TypeError):
                raise TypeError("polynomial (=%s) must be a polynomial." % repr(polynomial))
                            
        # Assign names
        if name is None:
            name = str(polynomial.parent().gen(0))
        
        if isinstance(ring.base_ring(), FiniteField):
            self._base_GF = ring.base_ring()
            self._ground_field = self._base_GF
            base_GF_map = self._base_GF.Hom(self._base_GF).identity()
        else:
            self._base_GF = ring.base_ring().GF()
            self._ground_field = ring.base_ring().ground_field()
            base_GF_map = ring.base_ring().GF_map()
        if not polynomial.map_coefficients(base_GF_map).is_irreducible():
            raise ValueError("Polynomial is reducible.")
            
        self._GF = self._base_GF.extension(polynomial.degree())
        base_GF_to_self_GF = self._GF.coerce_map_from(self._base_GF) ## This coercion map actually does exist.
        PolynomialQuotientRing_generic.__init__(self, ring, polynomial, name, category)
        self._refine_category_(Fields())
    
        #Ensure we can map elements from self to self.__GF
        from sage.categories.homset import Hom
        from sage.categories.morphism import SetMorphism

        alph = polynomial.map_coefficients(base_GF_map).change_ring(self.GF()).any_root()
        
        homspace = Hom(self, self.GF())
        self._GF_map = homspace.__make_element_class__(SetMorphism)(homspace,
                lambda f: f.lift().map_coefficients(base_GF_map)(alph))

        # We use Julian's code here with some optimizations for this class. 

        # POSSIBLE IMPROVEMENTS:
        # - replace Hom object with a vector space morphism over the prime field.
        
        basis, basis_in_GF = self._basis_over_prime_field_()
        assert(len(basis) == self.GF().degree())
        from sage.matrix.constructor import matrix
        
        A = matrix([b._vector_() for b in basis_in_GF])
        assert(A.is_square())
        x = A.solve_left(A.column_space().basis()[1])
        primitive_element = sum(c*b for c,b in zip(x.list(), basis))
        self._GF_map_inverse = self.GF().hom([primitive_element], self, check=False)

    # Element class    
    Element = RelativeFiniteFieldElement

    @cached_method
    def _basis_over_prime_field_(self ):
        """
        Compute a basis of self over the prime field as well as a corresponding basis of self.GF()
        """
        g = self.GF_map()(self.gen())
        if isinstance(self.base_field(), FiniteField):
            A = [self.base_field().gen() ** j for j in range(self.base_field().degree()) ]
            B = [self.base_field().gen() ** j for j in range(self.base_field().degree()) ]
        else:
            A,B = self.base_field()._basis_over_prime_field_()
            
        return ([self.gen() ** i * a for i in range(self.degree()) for a in A],
             [g ** i * b for i in range(self.degree()) for b in B])
    
    def drop_base(self, name=None):
        """
        Construct a new relative extension tower with the base field removed from the chain. The ground field of the tower cannot be removed.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL_over_k = LL.drop_base(); LL_over_k
            Relative Finite Field with defining polynomial 6*t^14 + 6*t^8 + (4*z2 + 2)*t^6 + (2*z2 + 1)*t^4 + (4*z2 + 1)*t^2 + 2*z2 over Finite Field in z2 of size 7^2
            sage: LL.base_field() == L
            True
            sage: LL_over_k.base_field() == k
            True
            sage: LL_over_k.relative_degree()
            14
            
            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: RRR.<z> = PolynomialRing(LL)
            sage: LLL = LLL = LL.extension(z^5 - LL.0 - k.0) 
            sage: wubwubwub = LLL.drop_base()
            sage: wubwubwub.base_field() == L
            True

        """
        if isinstance(self.base_ring(), FiniteField):
            raise TypeError("Cannot remove ground field from tower")
        current_base_field = self.base_field()
        base_of_base_field = current_base_field.base_field()

        
        base_modulus_coefficients = current_base_field.defining_polynomial().list()
        self_modulus_coefficients = self.defining_polynomial().list()
        self_modulus_coefficient_lifts = [ coeff.lift() for coeff in self_modulus_coefficients]
        
        # Construct a bivariate polynomial ring containing both defining polynomials
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        R = PolynomialRing(base_of_base_field, ['x','y'])
        x = R.gens()[0]
        y = R.gens()[1]
        
        base_modulus_R = sum([(x ** i) * base_modulus_coefficients[i] for i in range(current_base_field.relative_degree()+1)])
        
        self_modulus_R = sum([(y ** j) * sum([ (x ** i) * self_modulus_coefficient_lifts[j][i]
                    for i in range(self_modulus_coefficient_lifts[j].degree()+1)])
                    for j in range(len(self_modulus_coefficient_lifts))])
        
        new_modulus = self_modulus_R.resultant(base_modulus_R)
        ring = PolynomialRing(base_of_base_field, 't')
        t = ring.gen()
        
        if isinstance(base_of_base_field, FiniteField):
            return RelativeFiniteField(ring, new_modulus(0,t), name)
        return base_of_base_field.extension(new_modulus(0,t), name)
        
    def GF_map(self):
        """
        Returns the map from this relative finite field to the isomorphic finite field.

        The map is computed upon initialization of the RelativeFiniteField object. It is not canonical, and its construction involves the (random) selection of a root of the defining polynomial in self.GF()

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.GF_map()
            
            Generic morphism:
              From: Relative Finite Field with defining polynomial y^2 + 6*x + 6*z2 over Relative Finite Field with defining polynomial x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3 over Finite Field in z2 of size 7^2
              To:   Finite Field in z28 of size 7^28
            sage: LL.GF_map()(LL.0)
            2*z28^27 + z28^26 + 3*z28^25 + 4*z28^24 + 2*z28^23 + 6*z28^22 + 5*z28^21 + 6*z28^20 + 3*z28^19 + 4*z28^18 + 5*z28^17 + 5*z28^16 + z28^15 + z28^14 + 4*z28^12 + 6*z28^11 + 6*z28^9 + 4*z28^8 + 4*z28^7 + 2*z28^5 + z28^4 + 6*z28^3 + 5*z28 + 2
        """
        return self._GF_map

    def GF_map_inverse(self):
        """
        Return the map to ``self`` from ``self.GF()``, the isomorphic finite field. This map is the inverse to ``self.GF_map()``.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.GF_map_inverse()
            
            Ring morphism:
              From: Finite Field in z28 of size 7^28
              To:   Relative Finite Field with defining polynomial y^2 + 6*x + 6*z2 over Relative Finite Field with defining polynomial x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3 over Finite Field in z2 of size 7^2
              Defn: z28 |--> ((4*z2 + 6)*x^6 + 6*x^5 + (z2 + 6)*x^4 + 4*z2*x^3 + (4*z2 + 6)*x^2 + (6*z2 + 6)*x + 2*z2 + 4)*y + (2*z2 + 6)*x^6 + 3*x^5 + (4*z2 + 2)*x^4 + 3*z2*x^3 + (z2 + 4)*x^2 + (5*z2 + 3)*x + 2*z2 + 
            sage: (LL.GF_map_inverse() * LL.GF_map())(y)
            y
            sage: (LL.GF_map() * LL.GF_map_inverse()) ( LL.GF().0 ) 
            z28

        """
        return self._GF_map_inverse

    # Override the isomorphic ring method. We compute this data upon initialization
    def _isomorphic_ring(self):
        return (self.GF_map_inverse(), self.GF_map(), self.GF())
    
    def _repr_(self):
        try:
            return self._cached_repr
        except AttributeError:
            pass
        s = "Relative Finite Field with defining polynomial %s over %s"%(
                self.defining_polynomial(), self.base_ring())
        self._cached_repr = s
        return s
    
    def extension(self, polynomial, name=None):
        """
        Extend this relative finite field by joining the root of a polynomial with coefficients in this ring. 

        INPUT:
        
        -``polynomial`` -- a polynomial with coefficients in ``self``.

        -``name`` -- string (default None). The name for the generator of this extension.

        OUTPUTS:

            Outputs a new Relative finite field whose base ring is ``self`` and whose defining polynomial is ``polynomial``.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0); LL
            Relative Finite Field with defining polynomial y^2 + 6*x + 6*z2 over Relative Finite Field with defining polynomial x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3 over Finite Field in z2 of size 7^2

        """
        assert polynomial.parent().base_ring() == self, 'Polynomial not defined over the base field'
        return RelativeFiniteField(polynomial.parent(), polynomial, name)
        
    def degree_over_ground(self):
        """
        Returns the degree over the ground field.

        EXAMPLES::
            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.ground_field()
            Finite Field in z2 of size 7^2
            sage: LL.degree_over_ground()
            14

        """
        if self.ground_field() == self.base_ring():
            return self.degree()
        else:
            return self.base_field().degree_over_ground() * self.degree()

    def absolute_degree(self):
        """
        Return the degree of this relative finite field over its prime field.

        EXAMPLES::
            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.absolute_degree()
            28
        """
        if self.ground_field() == self.base_ring():
            return self.ground_field().degree() * self.degree()
        else:
            return self.base_field().absolute_degree() * self.degree()
    
    def relative_degree(self):
        """
        Return the degree of this relative finite field over its base field.
        Alias for degree()

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.relative_degree()
            2
        """
        return self.degree()
    
    def is_field(self, proof=True):
        """
        As a relative finite field is always a finite field, returns true.

        """
        return True

    def is_finite(self, proof=True):
        """
        As a relative finite field is always a finite field, returns true.

        """
        return True
    
    def ground_field(self):
        """
        Returns the ground field. That is, the field at the bottom of the relative finite field tower.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.ground_field()
            Finite Field in z2 of size 7^2
        """
        return self._ground_field
    
    def base_field(self):
        """
        Returns the field immediately below this field in the tower. Alias for base_ring().

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.base_field()
            Relative Finite Field with defining polynomial x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3 over Finite Field in z2 of size 7^2
            sage: LL.base_field == L
            True

        """
        return self.base_ring()
    
    def GF(self):
        """
        Returns the explicit finite field isomorphic to this relative finite field.

        OUTPUTS:

        A FiniteField object of the same cardinality as this ring.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.GF()
            Finite Field in z28 of size 7^28
        """
        return self._GF

    def defining_polynomial(self):
        """
        Returns the polynomial used to define this extension.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.defining_polynomial()
            y^2 + 6*x + 6*z2
        """
        return self.modulus()
    
    def base_GF(self):
        """
        Returns the explicit finite field isomorphic to ``self.base_field()``.

        OUTPUTS:

        A FiniteField object of the same cardinality as the base ring of this ring.

        EXAMPLES::

            sage: k = GF(49);
            sage: R.<x> = PolynomialRing(k);
            sage: g = x^7 + x^4 + 5*x^3 + 3*x^2 + 4*x + 3;
            sage: L = RelativeFiniteField(R,g)
            sage: RR.<y> = PolynomialRing(L)
            sage: LL = L.extension(y^2 - L.0 - k.0)
            sage: LL.base_GF()
            Finite Field in z14 of size 7^14
            sage: LL.base_GF() == L.GF()
            True
        """        
        return self._base_GF
           

        
