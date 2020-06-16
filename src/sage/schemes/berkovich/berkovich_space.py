"""
A framework for implementing the Berkovich construction over a scheme.
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element import Element
from sage.categories.topological_spaces import TopologicalSpaces
from sage.symbolic.expression import is_Expression
from sage.rings.real_mpfr import RR
from sage.rings.padics.padic_generic_element import pAdicGenericElement
from sage.rings.padics.factory import Qp
from sage.schemes.projective.projective_space import ProjectiveSpace
from sage.schemes.projective.projective_point import SchemeMorphism_point_projective_field
from sage.schemes.generic.morphism import is_SchemeMorphism
from sage.schemes.affine.affine_space import AffineSpace
from sage.schemes.generic.scheme import Scheme
from sage.rings.rational_field import QQ



class Berkovich_Element(Element):
    """
    The parent class for any element of a Berkovich space
    """
    pass

class Berkovich_Element_Cp(Berkovich_Element):
    """
    The parent class for any element of Berkovich space over ``Cp``. 
    This class should never be instantiated, instead Berkovich_Element_Cp_Affine
    or Berkovich_Element_Cp_Projective should be used.
    """

    def __init__(self, parent, center, radius=None, power=None, prec=20, child=None, error_check=True):
        Element.__init__(self, parent)
        from sage.rings.function_field.element import is_FunctionFieldElement
        from sage.rings.polynomial.polynomial_element import is_Polynomial
        from sage.rings.fraction_field_element import FractionFieldElement_1poly_field
        self._type = None

        #if center is a list, this is a Type 4 point
        if isinstance(center, list):
            if error_check:
                if not isinstance(radius, list):
                    raise ValueError("center was passed a list but radius was not a list")
                if len(radius) != len(center):
                    raise ValueError("The same number of centers and radii must be specified to create \
                        a Type IV point")
            self._center_lst = list(center)
            self._radius_lst = list(radius)
            self._prec = len(self._radius_lst)
            self._center_func = None
            self._radius_func = None
            self._type = 4
            if not error_check:
                return

        #is_FunctionFieldElement calls .parent
        elif hasattr(center, "parent"):
            #if center is a supported univariate function, this is a Type 4 point
            if is_FunctionFieldElement(center) or is_Polynomial(center) or\
                isinstance(center, FractionFieldElement_1poly_field) or is_Expression(center):
                #check that both center and radii are supported univariate function
                center_expr_check = False
                radius_expr_check = False
                if error_check:
                    if is_Expression(center):
                        if len(center.variables()) != 1:
                            raise ValueError("An expression with %s variables cannot define \
                                the centers approximating a Type IV point" %(len(center.variables())))
                        else:
                            #we do this since .subs is currently bugged for polynomials but not expressions
                            center_expr_check = True
                    if not (is_FunctionFieldElement(radius) or is_Polynomial(radius) or\
                        isinstance(radius, FractionFieldElement_1poly_field) or is_Expression(radius)):
                        raise ValueError("center was passed a function but radius was not a function")
                    if is_Expression(radius):
                        if len(radius.variables()) != 1:
                            raise ValueError("An expression with %s variables cannot define the radii \
                                approximating a Type IV point" %(len(radius.variables())))
                        else:
                            radius_expr_check = True
                else:
                    if is_Expression(center):
                        center_expr_check = True
                    if is_Expression(radius):
                        radius_expr_check = True
                self._type = 4
                self._prec = prec
                center_lst = []
                radius_lst = []
                self._center_func = center
                self._radius_func = radius
                if center_expr_check:
                    x = self._center_func.variables()[0]
                if radius_expr_check:
                    y = self._radius_func.variables()[0]
                for i in range(1,self._prec+1):
                    if center_expr_check:
                        #we use .subs for expressions to avoid deprecation
                        center_lst.append(self._center_func.subs({x:i}))
                    else:
                        #.subs for polynomials is currently bugged
                        center_lst.append(self._center_func(i))
                    if radius_expr_check:
                        radius_lst.append(self._radius_func.subs({y:i}))
                    else:
                        radius_lst.append(self._radius_func(i))
                self._center_lst = center_lst
                self._radius_lst = radius_lst
                if not error_check:
                    return

        if self._type == 4 and error_check:
            from sage.rings.real_mpfr import is_RealNumber
            if child == "projective":
                for i in range(len(self._center_lst)):
                    center = self._center_lst[i]
                    radius = self._radius_lst[i]
                    #make sure the center is a point of projective space and not the point at infinity
                    if not isinstance(center, SchemeMorphism_point_projective_field):
                        try:
                            center = (self._base_space)(center)
                            self._center_lst[i] = center
                        except:
                            raise ValueError("The center of a point of Projective Berkovich \
                                space must be a point of P^1(Cp(%s)), not of %s"%(self._p, center.parent()))
                    if center.scheme().ambient_space() != center.scheme():
                        raise ValueError("The center of a point of Berkovich space over P^1(Cp(%s)) \
                            cannot be a point of %s" %(self._p,center.scheme()))
                    if (center.scheme()).base_ring().prime() != self._p:
                        raise ValueError("The center of a disk approximating a Type IV point of Berkovich space\
                             over P^1(Cp(%s)) cannot be a point of %s") %(self._p, center.scheme())
                    if center == (center.scheme())((1,0)):
                        raise ValueError("The center of a disk approximating a Type IV point of Berkovich \
                            space cannot be centered at %s" %((center.scheme())((1,0))))
                    #make sure the radius coerces into the reals
                    if not is_RealNumber(radius):
                        if is_Expression(radius):
                            radius = RR(radius)
                        elif RR.has_coerce_map_from(radius.parent()):
                            radius = RR(radius)
                            self._radius_lst[i] = radius
                        else:
                            raise ValueError("The radius of a disk approximating a Type IV point \
                                must coerce into the real numbers, %s does not coerce" %(radius))
                    if i != 0:
                        #check containment for the sequence of disks
                        previous_center = self._center_lst[i-1]
                        previous_radius = self._radius_lst[i-1]
                        dist = (center[0] - previous_center[0]).abs()
                        if previous_radius < radius or dist > radius:
                            raise ValueError("Sequence of disks does not define a Type IV point as \
                                containment is not proper")
                return
            elif child == "affine":
                for i in range(len(self._center_lst)):
                    center = self._center_lst[i]
                    radius = self._radius_lst[i]
                    #make sure the center is in Cp
                    if not isinstance(center, pAdicGenericElement):
                        try:
                            center = (self._base_space)(center)
                            self._center_lst[i] = center
                        except:
                            raise ValueError("The center of a disk approximating a Type IV point must \
                                be padic or coerce into Qp, %s does not coerse into Qp" %(center))
                    if (center.parent()).prime() != self._p:
                        raise ValueError("The center of a disk approximating a Type IV point of Berkovich \
                            space over Cp(%s) cannot be a point of %s") %(self._p, center.parent())
                    #make sure the radius coerces into the reals
                    if not is_RealNumber(radius):
                        if is_Expression(radius):
                            radius = RR(radius)
                        elif RR.has_coerce_map_from(radius.parent()):
                            radius = RR(radius)
                            self._radius_lst[i] = radius
                        else:
                            raise ValueError("The radius of a disk approximating a Type IV point must \
                                coerce into the real numbers, %s does not coerce" %(radius))
                    if i != 0:
                        #check containment for the sequence of disks
                        previous_center = self._center_lst[i-1]
                        previous_radius = self._radius_lst[i-1]
                        dist = (center - previous_center).abs()
                        if previous_radius < radius or dist > radius:
                            raise ValueError("Sequence of disks does not define a Type IV point as \
                                containment is not proper")
                return
            else:
                raise ValueError("bad value %s passed to child. Do not initialize Berkovich_Element_Cp \
                    directly" %(child))
            return

        #the point must now be Type 1, 2, or 3, so we check that the center is of the appropriate type
        if child == "affine" and error_check:
            if not isinstance(center, pAdicGenericElement):
                try:
                    center = (self._base_space)(center)
                except:
                    raise ValueError("The center of a point of Affine Berkovich space over %s must convert \
                        to %s" %(self._base_space,self._base_space))
            if (center.parent()).prime() != self._p:
                raise ValueError("The center of a point of Berkovich space over Cp(%s) cannot be a point of %s")\
                    %(self._p, center.parent())
        elif child == "projective" and error_check:
            if not isinstance(center, SchemeMorphism_point_projective_field):
                try:
                    center = (self._base_space)(center)
                except:
                    raise ValueError("The center of a point of Projective Berkovich space must be a \
                        point of P^1(Cp) or coerce into P^1(Qp)") #TODO change to Qpbar
            if not(center.scheme().ambient_space() is center.scheme()):
                    raise ValueError("The center of a point of Projective Berkovich space cannot be \
                        a point of %s" %(center.scheme()))
            if (center.scheme()).base_ring().prime() != self._p:
                raise ValueError("The center of a point of Berkovich space over P^1(Cp(%s)) cannot \
                    be a point of %s") %(self._p, center.scheme())
        elif child != "affine" and child != "projective":
            raise ValueError("bad value %s passed to child. Do not initialize Berkovich_Element_Cp \
                    directly" %(child))
        self._center = center
        if (radius == None and power == None) or radius == 0:
            self._type = 1
            self._radius = 0
            self._power = None
            return
        #In order to simplify our representation, Type II and III points cannot be centered at infinity
        if child == "projective":
            try:
                center.dehomogenize(1)
            except ValueError:
                raise ValueError("Type II and III points cannot be centered at (1 : 0)")
        if power != None:
            if error_check:
                if power.parent() != QQ:
                    try:
                        power = QQ(power)
                    except:
                        raise ValueError("power must convert to rationals")
                if radius != None:
                    if radius != self._p**power:
                        raise ValueError("Conflicting inputs for power and radius")
            self._power = power
            self._radius = RR(self._p**power)
            self._type = 2
            return
        if radius != None:
            if is_Expression(radius):
                try:
                    power = QQ((radius.log(self._p)).expand_log())
                except:
                    pass
                try:
                    radius = RR(radius)
                    self._radius = radius
                except:
                    raise ValueError("Symbolic radius must be a real number")
            from sage.rings.real_mpfr import is_RealNumber
            if (not is_RealNumber(radius)) and power == None:
                if RR.has_coerce_map_from(radius.parent()):
                    self._radius = RR(radius)
                else:
                    raise ValueError("Radius must coerce into real numbers")
            else:
                self._radius = radius
            if power != None:
                self._power = power
                self._type = 2
                return
            power = radius.log(self._p)
            if power.is_integer():
                self._power = QQ(power)
                self._type = 2
            else:
                self._type = 3
                self._power = None

    def radius(self):
        """
        Radius of the corresponding disk (or sequence of disks) in ``Cp``.

        OUTPUT:

        - For Type I, II, and III points, returns a list of a single non-negative 
          real number. For Type IV points, returns a list of non-negative real numbers.
        """
        if self._type == 4:
            return self._radius_lst[:]
        return self._radius

    def diameter(self):
        """
        Diameter function on Berkovich space.

        For Type I, II, and III points, returns the radius.
        For Type IV points returns either the last radius
        in the finite approximation, or if a generating function
        was given for the radii, the diameter is computed
        as the limit of the function as it's variable tends
        to infinity.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(3)
            sage: Q1.diameter()
            0

            sage: Q2 = B(1/2,9)
            sage: Q2.diameter()
            9.00000000000000

        """
        if self._type == 4:
            if self._radius_func == None:
                return self._radius_lst[len(self._radius_lst)-1]
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(QQ,names="x")
            x = R.gens()[0]
            if is_Expression(self._radius_func):
                radius_func_variable = self._radius_func.variables()[0]
                radius_expr = self._radius_func.subs({radius_func_variable:x})
            else:
                radius_expr = self._radius_func(x)
                from sage.symbolic.ring import SymbolicRing as SR
                radius_expr = SR(RR)(radius_expr)
            return radius_expr.limit(x="oo")
        return self._radius

    def hyperbolic_metric(self, other):
        """
        Returns the hypebolic distance metric distance between ``self`` and ``other``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(1/4,4)
            sage: Q2 = B(1/4,6)
            sage: Q1.hyperbolic_metric(Q2)
            0.369070246428542

        """
        if not isinstance(other,Berkovich_Element_Cp):
            raise ValueError("Hyperbolic metric not defined between point of Berkovich space and %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("Input to hyperbolic metric was an element of a different Berkovich space")
        if self.type_of_point() == 1 or other.type_of_point() == 1:
            if self == other:
                return 0
            else:
                return "oo"
        return 2*((self.join(other)).diameter().log(self.prime())) - self.diameter().log(self.prime()) - other.diameter().log(other.prime())

    def small_metric(self, other):
        """
        Returns the small metric distance between ``self`` and ``other``.

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(1/4,4)
            sage: Q2 = B(1/4,6)
            sage: Q1.small_metric(Q2)
            2.00000000000000

        """
        if not isinstance(other,Berkovich_Element_Cp):
            raise ValueError("Hyperbolic metric not defined between point of Berkovich space and %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("Input to hyperbolic metric was an element of a different Berkovich space")
        return 2*((self.join(other)).diameter()) - self.diameter() - other.diameter()

    def Hsia_kernel_infinity(self, other):
        """
        Return the Hsia kernel at infinity evaluated at the point (``self``,``other``).

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(1/4,4)
            sage: Q2 = B(1/4,6)
            sage: Q1.Hsia_kernel_infinity(Q2)
            6.00000000000000

        """
        return self.join(other).diameter()

    def center(self):
        """
        Returns the center of the corresponding disk (or sequence of disks) in ``Cp``

        OUTPUT:

        - For Type I, II, and III points, returns an element of ``Qp`` or a finite extension. For
          Type IV points, returns a list of points of ``Qp`` or a finite extension
        """
        if self._type == 4:
            return self._center_lst[:]
        return self._center

    def type_of_point(self):
        """
        Returns the Type of this point of Berkovich space over ``Cp``
        """
        return self._type

    def prime(self):
        """
        Shorthand for residue characteristic of the parent
        """
        return self._p
    
    def __ne__(self, other):
        """
        Non-equality operator.
        """
        return not (self == other)

    def _repr_(self):
        if self._type == 1:
            return "Type I point centered at " + format(self._center)
        elif self._type == 2:
            return "Type II point centered at " + format(self._center) + " of radius %s^%s" %(self._p,self._power)
        elif self._type == 3:
            return "Type III point centered at " + format(self._center) + " of radius " + format(self._radius)
        else:
            if self._center_func != None and self._radius_func != None:
                return "Type IV point of precision %s with centers given by %s and radii given by %s" %(self._prec, self._center_func,self._radius_func)
            else:
                return "Type IV point of precision %s, approximated by disks centered at %s ... with radii %s ..." %(self._prec, self._center_lst[:2],self._radius_lst[:2])

    def _latex_(self):
        from sage.misc.latex import latex
        if self._type == 1:
            text = r"the point %s of } \Bold{C}_%s" %(self._center, self._p)
        elif self._type in [2,3]:
            text = r"the disk centered at %s of radius %s in } \Bold{C}_%s" %(self._center, self._radius, self._p)
        else:
            text = r"the sequence of disks with centers %s } \ldots \text{ and radii %s } \ldots" %(self._center_lst[:2],self._radius_lst[:2])
        return r"\text{Type %s Point of }" %(self._type) + latex(self.parent()) + r"\text{equivalent to " + text

class Berkovich_Element_Cp_Affine(Berkovich_Element_Cp):
    r"""
    Element class of the Berkovich affine line over `\CC_p`

    Elements are categorized into four Types, represented by specific data:

     - Type I points are represented by a center in `\QQ_p` or a finite extension

     - Type II points are represented by a center in `\QQ_p` and a rational power of `p`

     - Type III points are represented by a center in `\QQ_p` and a radius in `[0,\infty)`

     - Type IV points are represented by a finite list of centers in `\QQ_p` and 
       a finite list of radii in `[0,\infty)`. 

    INPUT:

    - ``center`` -- For Type I, II, and III points, the center of the 
     corresponding disk in `\CC_p`. Must be an element of `\QQ_p`, a finite extension
     of `\QQ_p`, or coerce into `\QQ_p`. For Type IV points, can be a list of centers
     used to approximate the point or a univariate function that computes the centers 
     (computation starts at 1).

    - ``radius`` -- (optional) For Type I, II, and III points, the radius of the 
     corresponding disk in ``Cp``. Must coerce into the real numbers. For Type IV points,
     can be a list of radii used to approximate the point or a univariate function that
     computes the radii (computation starts at 1). 

    - ``power`` -- (optional) Rational number. Used for constructing Type II points; specifies
     the power of ``p`` such that p^power = radius

    - ``prec`` -- (default: 20) The number of disks to be used to approximate a Type IV point

    - ``error_check`` -- (default: True) If error checking should be run on input. If
     input is correctly formatted, can be set to ``False`` for better performance. 
     WARNING: with error check set to ``False``, any error in the input will lead to 
     incorrect results.

    EXAMPLES:

    Type I points can be created by specifying the corresponding point of ``Cp``::

        sage: B = Berkovich_Cp_Affine(Qp(3))
        sage: a = B(4)
        sage: a
        Type I point centered at 1 + 3 + O(3^20)

    The center of a point can be an element of a finite extension of ``Qp``::

        sage: A.<t> = Qq(27)
        sage: a = B(1+t)
        sage: a
        Type I point centered at (t + 1) + O(3^20)

    Type II and III points can be created by specifying a center and a radius::

        sage: b = B(2, 3**(1/2)); b
        Type II point centered at 2 + O(3^20) of radius 3^1/2
        sage: c = B(2,1.6); c
        Type III point centered at 2 + O(3^20) of radius 1.60000000000000
    
    Some Type II points may be mistaken for Type III points::

        sage: b = B(3, 3**0.5); b
        Type III point centered at 3 + O(3^21) of radius 1.73205080756888

    To avoid these errors, specify the power instead of the radius::

        sage: b = B(3, power=RR(1/100000)); b
        Type II point centered at 3 + O(3^21) of radius 3^1/100000
    
    Type IV points can be constructed in a number of ways, the first being
    from a list of centers and radii used to approximate the point::

        sage: d = B([Qp(3)(2), Qp(3)(2), Qp(3)(2)], [1.761, 1.123, 1.112]); d
        Type IV point of precision 3, approximated by disks centered at 
        [2 + O(3^20), 2 + O(3^20)] ... with radii [1.76100000000000, 1.12300000000000] ...

    Type IV points can be constructed from univariate functions, with arbitrary precision::

        sage: A.<t> = Qq(27)
        sage: R.<x> = PolynomialRing(A)
        sage: f = (1+t)^2*x
        sage: S.<y> = PolynomialRing(RR)
        sage: S = FractionField(S)
        sage: g = (y+1)/y
        sage: d = B(f,g,prec=100); d
        Type IV point of precision 100 with centers given by 
        ((t^2 + 2*t + 1) + O(3^20))*x and radii given by (y + 1.00000000000000)/y

    For increased performance, error_check can be set to ``False``. WARNING: with error check set
    to ``False``, any error in the input will lead to incorrect results::

        sage: d = B(f,g,prec=100,error_check=False); d
        Type IV point of precision 100 with centers given by 
        ((t^2 + 2*t + 1) + O(3^20))*x and radii given by (y + 1.00000000000000)/y
        
    TESTS::

        sage: A= Berkovich_Cp_Affine(Qp(3))
        sage: Q1=A(3,1); Q1
        Type II point centered at 3 + O(3^21) of radius 3^0

        sage: Q2=A(2.5,1); Q2
        Type II point centered at 1 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 3^7 + 
        3^8 + 3^9 + 3^10 + 3^11 + 3^12 + 3^13 + 3^14 + 3^15 + 3^16 + 3^17 + 3^18 + 3^19 + O(3^20) of radius 3^0

        sage: Q5=A(3,0); Q5
        Type I point centered at 3 + O(3^21)

        sage: Q1 == Q2
        True

        sage: Q1 == Q5
        False

        sage: Q3 = A(Qp(3)(3),power=0,error_check=False); Q3
        Type II point centered at 3 + O(3^21) of radius 3^0

        sage: Q4 = A(3,3**0); Q4
        Type II point centered at 3 + O(3^21) of radius 3^0

        sage: Q5 = A(3, power = 1/2); Q5
        Type II point centered at 3 + O(3^21) of radius 3^1/2

        sage: Q6 = A(3, RR(3**(1/2))); Q6
        Type III point centered at 3 + O(3^21) of radius 1.73205080756888

        sage: Q5 == Q6
        True

        sage: k = Qp(5)
        sage: R.<x> = k[]
        sage: l.<w> = k.extension(x^2-5)
        sage: B=Berkovich_Cp_Affine(Qp(5))
        sage: B(w,power=1)
        Type II point centered at w + O(w^41) of radius 5^1

    """
    def __init__(self, parent, center, radius=None, power=None, prec=20, error_check=True):
        #we call Berkovich_Element_Cp constructor which is shared with projective Berkovich space
        #unless we are passed a point of projective Berkovich space
        self._p = parent.prime()
        self._base_space = parent.base()

        #if this is a point of projective berkovich space, we do the conversion
        if isinstance(center, Berkovich_Element_Cp_Projective):
            if (center.prime() == self._p):
                Element.__init__(self, parent)
                try:
                    center.center().dehomogenize(1)
                except:
                    raise ValueError("The point at infinity of Berkovich Projective space does not \
                        convert to Berkovich Affine space")
                self._center = center.center()[0]
                self._radius = center._radius
                self._power = center._power
                self._type = center._type
                if self._type == 4:
                    self._prec = center._prec
                    self._center_func = center._center_func
                    self._radius_func = center._center_func
                    self._radius_lst = center._radius_lst
                    center_lst = []
                    for i in center._center_lst:
                        center_lst.append(i[0])
                    self._center_lst = center_lst
                return
            else:
                raise ValueError("Tried to convert from a point of Berkovich space over Cp(%s) to a \
                    point of Berkovich space over Cp(%s)" %(center._base_space.base_ring().prime(), \
                    parent.base().prime()))

        #otherwise we call the Berkovich_Element_Cp constructor with child=affine
        Berkovich_Element_Cp.__init__(self,parent=parent,center=center,radius=radius,power=power,\
            prec=prec,child="affine",error_check=error_check)

    def __eq__(self, other):
        """
        Equality operator.
        """
        if other is self:
            return True
        if not isinstance(other, Berkovich_Element_Cp_Affine):
            return False
        if other.parent() != self.parent():
            return False
        stype = self.type_of_point()
        otype = other.type_of_point()
        if stype == otype and stype == 1:
            return self.center() == other.center()
        elif stype == otype and stype == 4:
            raise NotImplementedError("Equality for Type IV points not yet implemented")
        elif stype in [2,3] and otype in [2,3]:
            if self.radius() != other.radius():
                return False
            center_dist = (self.center() - other.center()).abs()
            return center_dist <= self.radius()
        else:
            return False

    def partial_order(self,other):
        """
        The standard partial order on Berkovich space

        INPUT:
         - ``other`` - another point of the same Berkovich space

        OUTPUT:
         - ``True`` - If self > other in the standard partial order
         - ``False`` - If self < other in the standard partial order
         - ``None`` - If the two points are not comparable

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(2,4)
            sage: Q2 = B(2,6)
            sage: Q1.partial_order(Q2)
            False

            sage: Q3 = B(1/2)
            sage: Q1.partial_order(Q3)
            True

            sage: Q4 = B(1/81,1)
            sage: print(Q4.partial_order(Q1))
            None

        """
        if not isinstance(other,Berkovich_Element_Cp_Affine):
            raise ValueError("partial order takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("partial order takes another point of the same Berkovich Affine line")
        if self.type_of_point() in [1,4]:
            if other.type_of_point() in [1,4]:
                if self == other:
                    return True
                else:
                    return None
            else:
                return False
        else:
            if other.type_of_point() in [1,4]:
                return True
            else:
                dist = (self.center()-other.center()).abs()
                if(dist <= self.radius()):
                    if(other.radius() <= self.radius()):
                        return True
                    else:
                        return False
                else:
                    if(dist <= other.radius()):
                        return False
                    else:
                        return None

    def join(self, other, basepoint="infty"):
        """
        Computes the join of ``self`` and ``other`` with respect to ``basepoint``.

        The join is first point that lies on the interesection
        of the path from self to basepoint and the path from other to
        basepoint.

        INPUT:

         - ``other`` -- the second point with which to take the join
         - ``basepoint`` -- (default: Infinity) the basepoint with which
         to take the join with respect to

        EXAMPLES::

            sage: B = Berkovich_Cp_Affine(Qp(3))
            sage: Q1 = B(2,1)
            sage: Q2 = B(2,2)
            sage: Q1.join(Q2)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000

            sage: Q3 = B(5)
            sage: Q3.join(Q1)
            Type II point centered at 2 + 3 + O(3^20) of radius 3^0

            sage: Q3.join(Q1, basepoint=Q2)
            Type II point centered at 2 + O(3^20) of radius 3^0

        TESTS::

            sage: Q4 = B(1/3**8+2,1)
            sage: Q2.join(Q4,basepoint = Q1)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000

            sage: Q5 = B(2,1/9)
            sage: Q6 = B(1,1/27)
            sage: Q4.join(Q5,basepoint=Q6)
            Type II point centered at 1 + O(3^20) of radius 3^0

            sage: Q7 = B(1/27,1/27)
            sage: Q1.join(Q7,Q2)
            Type III point centered at 2 + O(3^20) of radius 2.00000000000000
        """
        if not isinstance(other,Berkovich_Element_Cp_Affine):
            raise ValueError("join of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("join takes two points in the same Berkovich Affine line")
        if self.type_of_point() == 4 or other.type_of_point() == 4:
            raise NotImplementedError("join with Type IV points not implemented")
        if basepoint == "infty":
            dist = (self.center() - other.center()).abs()
            return self.parent()(self.center(),max(dist,self.radius(),other.radius()))
        else:
            if not isinstance(basepoint, Berkovich_Element_Cp_Affine):
                raise ValueError("basepoint for join must be a point of the Berkovich Affine line over Cp")
            if basepoint.parent() != self.parent():
                raise ValueError("basepoint must be a point of the same Berkovich Affine line")
            if basepoint.type_of_point() == 4:
                raise NotImplementedError("join not implemented for Type IV basepoint")
            b_greater_than_s = basepoint.partial_order(self)
            b_greater_than_o = basepoint.partial_order(other)
            if b_greater_than_s == None and b_greater_than_o == None:
                dist_b_s = (self.center() - basepoint.center()).abs()
                dist_b_o = (other.center() - basepoint.center()).abs()
                return self.parent()(basepoint.center(),min(max(dist_b_o,other.radius(),basepoint.radius()),max(dist_b_s,self.radius(),basepoint.radius())))
            elif b_greater_than_s != None and b_greater_than_o == None:
                dist_b_o = (other.center() - basepoint.center()).abs()
                if dist_b_o > basepoint.radius() and dist_b_o < self.radius(): #TODO remove the check on basepoint.radius()
                    return (self.parent())(basepoint.center(),dist_b_o)
                elif dist_b_o > basepoint.radius() and dist_b_o >= self.radius():
                    if b_greater_than_s:
                        return basepoint
                    else:
                        return self
                else:
                    raise NotImplementedError("impossible case in join") #TODO remove this
            elif b_greater_than_s == None and b_greater_than_o != None:
                return other.join(self,basepoint)
            else:
                if b_greater_than_s:
                    if b_greater_than_o:
                        if self.partial_order(other):
                            return self
                        else:
                            return other
                    else:
                        return basepoint
                else:
                    if b_greater_than_o:
                        return basepoint
                    else:
                        if self.partial_order(other):
                            return other
                        else:
                            return self

    def potential_kernel(self, other, basepoint):
        """
        The potential kernel of ``self`` with ``other``, with basepoint ``basepoint``.
        """
        if not isinstance(other,Berkovich_Element_Cp_Affine):
            raise ValueError("potential kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("potential kernel takes two points in the same Berkovich Affine line")
        if not isinstance(basepoint, Berkovich_Element_Cp_Affine):
            raise ValueError("basepoint must be a point of the Berkovich Affine line over Cp")
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich Affine line")
        return basepoint.hyperbolic_metric((self.join(other,basepoint)))

    def spherical_kernel(self,other):
        """
        The spherical kernel of ``self`` with ``other``.

        Equivalently, the Hsia kernel with basepoint the Gauss
        point.
        """
        if not isinstance(other,Berkovich_Element_Cp_Affine):
            raise ValueError("spherical kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("spherical kernel takes two points in the same Berkovich Affine line")
        gauss_point = (self.parent())(0,1)
        w = self.join(other,gauss_point)
        dist = gauss_point.hyperbolic_metric(w)
        if dist == "oo":
            return 0
        return (self.prime())**(-1*dist)

    def Hsia_kernel(self, other, basepoint):
        """
        The Hsia kernel of ``self`` and ``other``, with basepoint ``basepoint``.
        """
        if not isinstance(other,Berkovich_Element_Cp_Affine):
            raise ValueError("Hsia kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("Hsia kernel takes two points in the same Berkovich Affine line")
        if not isinstance(basepoint, Berkovich_Element_Cp_Affine):
            raise ValueError("basepoint must be a point of the Berkovich Affine line over Cp")
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich Affine line")
        if basepoint.type_of_point() == 1:
            if self == basepoint or other == basepoint:
                return "oo"
        return (self.spherical_kernel(other))/(self.spherical_kernel(basepoint)*other.spherical_kernel(basepoint))
    
    def diameter(self, basepoint="oo"):
        """
        Generalized diameter function.

        If the basepoint is infinity, the diameter is equal to 
        the limit of the radii of the corresponding disks in ``Cp``.

        If the basepoint is not infinity, the diameter
        is the Hsia kernel of ``self`` with itself at
        basepoint ``basepoint``.
        """
        if basepoint == "oo":
            return super().diameter()
        else:
            return self.Hsia_kernel(self,basepoint)

class Berkovich_Element_Cp_Projective(Berkovich_Element_Cp):
    r"""
    Element class of the Berkovich Projective line over `\CC_p`
    
    Elements are categorized into four Types, which are represented as follows:

     - Type I points are represented by a center in Projective Space of dimension 1 over
       `\QQ_p` or over a finite extension of `\QQ_p`

     - Type II points are represented by a center in Projective Space of dimension 1 over
       `\QQ_p` or over a finite extension of `\QQ_p` and a rational power of `p`.
       Type II points cannot be centered at the point at infinity.

     - Type III points are represented by a center in Projective Space of dimension 1 over
       `\QQ_p` or over a finite extension of `\QQ_p` and a radius in [0,`\infty`).
       Type III points are cannot be centered at the point at infinity.

     - Type IV points are represented by finite list of centers in Projective Space of dimension 1 over
       `\QQ_p` or over a finite extension of `\QQ_p` and by a finite list of radii in [0,`\infty`).
       None of the centers can be the point at infinity.

    INPUT:

    - ``center`` -- For Type I, II, and III points, the center of the 
     corresponding disk in `P^1(\CC_p)`. Must be an element of `\QQ_p`, a finite extension
     of `\QQ_p`, or coerce into `\QQ_p`. For Type IV points, can be a list of centers
     used to approximate the point or a univariate function that computes the centers 
     (computation starts at 1).

    - ``radius`` -- (optional) For Type I, II, and III points, the radius of the 
     corresponding disk in ``Cp``. Must coerce into the real numbers. For Type IV points,
     can be a list of radii used to approximate the point or a univariate function that
     computes the radii (computation starts at 1). 

    - ``power`` -- (optional) Rational number. Used for constructing Type II points; specifies
     the power of ``p`` such that p^power = radius

    - ``prec`` -- (default: 20) The number of disks to be used to approximate a Type IV point

    - ``error_check`` -- (default: True) If error checking should be run on input. If
     input is correctly formatted, can be set to ``False`` for better performance.
     WARNING: Setting error_check to ``False`` can lead to incorrect results.

    EXAMPLES::

    Type I points can be created by specifying the corresponding point of `P^1(\CC_p)`::

        sage: S = ProjectiveSpace(Qp(5),1)
        sage: P = Berkovich_Cp_Projective(S); P
        Projective Berkovich line over Cp(5) of precision 20

        sage: a = S((0,1))
        sage: Q1 = P(a); Q1
        Type I point centered at (0 : 1 + O(5^20))

        sage: Q2 = P((1,0)); Q2
        Type I point centered at (1 + O(5^20) : 0)

    Type II and III points can be created by specifying a center and a radius::

        sage: Q3 = P((0,5),5**(3/2)); Q3
        Type II point centered at (0 : 1 + O(5^20)) of radius 5^3/2

        sage: Q4 = P(0,3**(3/2)); Q4
        Type III point centered at (0 : 1 + O(5^20)) of radius 5.19615242270663

    Type IV points can be created from lists of centers and radii::

        sage: b = S((3,2)) #create centers
        sage: c = S((4,3))
        sage: d = S((2,3))
        sage: L = [b,c,d]
        sage: R = [1.761, 1.123, 1.112]
        sage: Q5 = P(L,R); Q5
        Type IV point of precision 3, approximated by disks centered at 
        [(4 + 2*5 + 2*5^2 + 2*5^3 + 2*5^4 + 2*5^5 + 2*5^6 + 2*5^7 + 2*5^8 + 2*5^9 + 2*5^10 +
         2*5^11 + 2*5^12 + 2*5^13 + 2*5^14 + 2*5^15 + 2*5^16 + 2*5^17 + 2*5^18 + 2*5^19 + O(5^20) : 
         1 + O(5^20)), (3 + 3*5 + 5^2 + 3*5^3 + 5^4 + 3*5^5 + 5^6 + 3*5^7 + 5^8 + 3*5^9 + 
         5^10 + 3*5^11 + 5^12 + 3*5^13 + 5^14 + 3*5^15 + 5^16 + 3*5^17 + 5^18 + 3*5^19 + O(5^20) : 
         1 + O(5^20))] ... with radii [1.76100000000000, 1.12300000000000] ...

    Type IV points can also be created from univariate functions. Since the centers of 
    the sequence of disks can not be the point at infinity in `P^1(\CC_p)`, only functions
    into `\CC_p` are supported::

        sage: L.<t> = PolynomialRing(Qp(5))
        sage: T = FractionField(L)
        sage: f = T(1/t)
        sage: R.<x> = RR[]
        sage: Y = FractionField(R)
        sage: g = (40*pi)/x
        sage: Q6 = P(f,g); Q6
        Type IV point of precision 20 with centers given by (1 + O(5^20))/((1 + O(5^20))*t)
         and radii given by 40.0000000000000*pi/x

    TEST::

        sage: P((1,0),3)
        Traceback (most recent call last):
        ...
        ValueError: Type II and III points cannot be centered at (1 : 0)

    """
    def __init__(self, parent, center, radius=None, power=None, prec=20, error_check=True):
        #if we are given a point of Affine Berkovich Space, we do the conversion
        #otherwise we call the Berkovich_Element_Cp constructor with child="projective"

        self._p = parent.prime()
        self._base_space = parent.base()

        #conversion from Affine points is handled in this constructor
        if isinstance(center, Berkovich_Element_Cp_Affine):
            if (center.prime() == self._p):
                Element.__init__(self, parent)
                self._center = (self._base_space)(center._center)
                self._radius = center._radius
                self._power = center._power
                self._type = center._type
                if self._type == 4:
                    self._prec = center._prec
                    self._center_func = center._center_func
                    self._radius_func = center._center_func
                    self._radius_lst = center._radius_lst
                    center_lst = []
                    for i in center._center_lst:
                        center_lst.append((self._base_space)(i))
                    self._center_lst = center_lst
                return
            else:
                raise ValueError("Tried to convert from a point of Berkovich space over Cp(%s) to a point \
                    of Berkovich space over Cp(%s)" %(center._base_space.base_ring().prime(),parent.base().prime()))

        Berkovich_Element_Cp.__init__(self,parent=parent,center=center,radius=radius,power=power,\
                prec=prec,child="projective",error_check=error_check)

    def _field_center(self):
        """
        Returns the center of ``self`` as an element 
        of a field, allowing for arithmetic.
        Used to make handling Type II and III points
        easier, as the centers of these points are never (1 : 0).
        """
        return self.center()[0]

    def __eq__(self, other):
        """
        Equality operator.
        """
        if other is self:
            return True
        if not isinstance(other, Berkovich_Element_Cp_Projective):
            return False
        if other.parent() != self.parent():
            return False
        stype = self.type_of_point()
        otype = other.type_of_point()
        if stype == otype and stype == 1:
            return self.center() == other.center()
        elif stype == otype and stype == 4:
            raise NotImplementedError("Equality for Type IV points not yet implemented")
        elif stype in [2,3] and otype in [2,3]:
            if self.radius() != other.radius():
                return False
            scent = self._field_center()
            ocent = other._field_center()
            center_dist = (scent - ocent).abs()
            return  center_dist <= self.radius()
        else:
            return False

    def partial_order(self,other):
        """
        The standard partial order on Berkovich space

        INPUT:
         - ``other`` - another point of the same Berkovich space

        OUTPUT:
         - ``True`` - If self > other in the standard partial order
         - ``False`` - If other > self in the standard partial order
         - ``None`` - If the two points are not comparable

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(ProjectiveSpace(Qp(3),1))
            sage: Q1 = B(2,4)
            sage: Q2 = B(2,6)
            sage: Q1.partial_order(Q2)
            False

            sage: Q3 = B(1/2)
            sage: Q1.partial_order(Q3)
            True

            sage: Q4 = B(1/81,1)
            sage: print(Q4.partial_order(Q1))
            None

        """
        if not isinstance(other,Berkovich_Element_Cp_Projective):
            raise ValueError("partial order takes another point of the Berkovich Projective line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("partial order takes another point of the same Berkovich Projective line")
        if self.type_of_point() in [1,4]:
            if other.type_of_point() in [1,4]:
                if self == other:
                    return True
                else:
                    return None
            else:
                return False
        else:
            if other.type_of_point() in [1,4]:
                return True
            else:
                dist = (self._field_center()-other._field_center()).abs()
                if(dist <= self.radius()):
                    if(other.radius() <= self.radius()):
                        return True
                    else:
                        return False
                else:
                    if(dist <= other.radius()):
                        return False
                    else:
                        return None

    def join(self, other, basepoint="infty"):
        """
        Computes the join of ``self`` and ``other``, with respect to ``basepoint``. 

        The join is first point that lies on the interesection
        of the path from self to basepoint and the path from other to
        basepoint.

        INPUT:
         
         - ``other`` -- the second point with which to take the join
         - ``basepoint`` -- (default: Infinity) the basepoint with which
         to take the join with respect to

        EXAMPLES::

            sage: B = Berkovich_Cp_Projective(ProjectiveSpace(Qp(3),1))
            sage: Q1 = B(2,1)
            sage: Q2 = B(2,2)
            sage: Q1.join(Q2)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000

            sage: Q3 = B(5)
            sage: Q3.join(Q1)
            Type II point centered at (2 + 3 + O(3^20) : 1 + O(3^20)) of radius 3^0

            sage: Q3.join(Q1, basepoint=Q2)
            Type II point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 3^0

        TESTS::

            sage: Q4 = B(1/3**8+2,1)
            sage: Q2.join(Q4,basepoint = Q1)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000

            sage: Q5 = B(2,1/9)
            sage: Q6 = B(1,1/27)
            sage: Q4.join(Q5,basepoint=Q6)
            Type II point centered at (1 + O(3^20) : 1 + O(3^20)) of radius 3^0

            sage: Q7 = B(1/27,1/27)
            sage: Q1.join(Q7,Q2)
            Type III point centered at (2 + O(3^20) : 1 + O(3^20)) of radius 2.00000000000000
        """
        if not isinstance(other,Berkovich_Element_Cp_Projective):
            raise ValueError("join of a point of the Berkovich Projective line takes another point of the Berkovich Projective line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("join takes two points in the same Berkovich Projective line")
        if self.type_of_point() == 4 or other.type_of_point() == 4:
            raise NotImplementedError("join with Type IV points not implemented")
        infty = (self.parent())((1,0))
        if basepoint == "infty" or basepoint == infty:
            if self == infty or other == infty:
                return infty
            else:
                dist = (self._field_center() - other._field_center()).abs()
                return self.parent()(self.center(),max(dist,self.radius(),other.radius()))
        else:
            if not isinstance(basepoint, Berkovich_Element_Cp_Projective):
                raise ValueError("basepoint for join must be a point of the Berkovich Affine line over Cp")
            if basepoint.parent() != self.parent():
                raise ValueError("basepoint must be a point of the same Berkovich Affine line")
            if basepoint.type_of_point() == 4:
                raise NotImplementedError("join not implemented for Type IV basepoint")
            if self == infty or other == infty:
                if self == other:
                    return infty
                else:
                    if self == infty:
                        return other.join(basepoint)
                    else:
                        return self.join(basepoint)
            else:
                base = self._base_space.base_ring()
                B = Berkovich_Cp_Affine(base)
                affine_self = B(self)
                affine_other = B(other)
                affine_basepoint = B(basepoint)
                affine_join = affine_self.join(affine_other,affine_basepoint)
                return (self.parent())(affine_join)

    def potential_kernel(self, other, basepoint):
        """
        The potential kernel of ``self`` with ``other``,
        with basepoint ``basepoint``.
        """
        if not isinstance(other,Berkovich_Element_Cp_Projective):
            raise ValueError("potential kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("potential kernel takes two points in the same Berkovich Affine line")
        if not isinstance(basepoint, Berkovich_Element_Cp_Projective):
            raise ValueError("basepoint must be a point of the Berkovich Affine line over Cp")
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich Affine line")
        return basepoint.hyperbolic_metric((self.join(other,basepoint)))

    def spherical_kernel(self,other):
        """
        The spherical kernel of ``self`` with ``other``.

        Equivalently, the Hsia kernel with basepoint the Gauss
        point.
        """
        if not isinstance(other,Berkovich_Element_Cp_Projective):
            raise ValueError("spherical kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("spherical kernel takes two points in the same Berkovich Affine line")
        gauss_point = (self.parent())(0,1)
        w = self.join(other,gauss_point)
        dist = gauss_point.hyperbolic_metric(w)
        if dist == "oo":
            return 0
        return (self.prime())**(-1*dist)

    def Hsia_kernel(self, other, basepoint):
        """
        The Hsia kernel of ``self`` and ``other``,
        with basepoint ``basepoint``
        """
        if not isinstance(other,Berkovich_Element_Cp_Projective):
            raise ValueError("Hsia kernel of a point of the Berkovich Affine line takes another point of the Berkovich Affine line, not %s" %(other))
        if self.parent() != other.parent():
            raise ValueError("Hsia kernel takes two points in the same Berkovich Affine line")
        if not isinstance(basepoint, Berkovich_Element_Cp_Projective):
            raise ValueError("basepoint must be a point of the Berkovich Affine line over Cp")
        if basepoint.parent() != self.parent():
            raise ValueError("basepoint must be a point of the same Berkovich Affine line")
        if basepoint.type_of_point() == 1:
            if self == basepoint or other == basepoint:
                return "oo"
        return (self.spherical_kernel(other))/(self.spherical_kernel(basepoint)*other.spherical_kernel(basepoint))
    
    def diameter(self, basepoint="oo"):
        """
        Generalized diameter function.

        If the basepoint is infinity, the diameter is equal to 
        the limit of the radii of the corresponding disks in ``Cp``.

        If the basepoint is not infinity, the diameter
        is the Hsia kernel of ``self`` with itself at
        basepoint ``basepoint``.
        """
        if basepoint == "oo":
            return super().diameter()
        else:
            return self.Hsia_kernel(self,basepoint)

class Berkovich(UniqueRepresentation,Parent):
    """
    The parent class for any Berkovich space
    """
    pass

class Berkovich_Cp_Affine(Berkovich):
    r"""
    The Berkovich Affine line over `\CC_p`.
    
    The Berkovich Affine line can be thought of as the set of seminorms on `\CC_p[x]`,
    with the weakest topology that makes the map `| \cdot | \to |f|` continuous
    for all `f \in \CC_p[x]`.

    INPUT:

     - ``base`` - The prime ``p``. Alternative, can be `\QQ_p` or a finite extension.
       This allows for more control over automated conversion of centers of points.

    EXAMPLES::

        sage: B = Berkovich_Cp_Affine(3); B
        Affine Berkovich line over Cp(3) of precision 20

    Initializing by passing in Qp space looks the same::

        sage: B = Berkovich_Cp_Affine(Qp(3)); B
        Affine Berkovich line over Cp(3) of precision 20

    However, this method allows for more control over behind-the-scenes conversion::

        sage: B = Berkovich_Cp_Affine(Qp(3,1)); B
        Affine Berkovich line over Cp(3) of precision 1

        sage: Q1 = B(1/2); Q1
        Type I point centered at 2 + O(3)

    Note that this point has very low precision, as B was initialized
    with a padic field of capped-relative precision one. For high precision,
    pass in a high precision padic field::

        sage: B = Berkovich_Cp_Affine(Qp(3,1000)); B
        Affine Berkovich line over Cp(3) of precision 1000
    """
    def __init__(self,base):
        from sage.rings.integer_ring import ZZ
        if base in ZZ:
            if base.is_prime():
                base = Qp(base) #TODO chance to Qpbar
            else:
                raise ValueError("Non-prime pased into Berkovich space")
        from sage.rings.padics.generic_nodes import is_pAdicField
        if not is_pAdicField(base): #TODO change base to Qpbar(prime)
            raise ValueError("Base of Berkovich Space must be a padic field")
        self._p = base.prime()
        Parent.__init__(self, base = base, category=TopologicalSpaces()) 

    def residue_characteristic(self):
        return self._p

    def prime(self):
        return self._p

    def _coerce_map_from_(self,S):
        if isinstance(S, Berkovich_Cp_Affine):
            if S.prime() == self.prime():
                return True
        return False

    def _repr_(self):
        return "Affine Berkovich line over Cp(%s) of precision %s" %(self.prime(),self.base().precision_cap())
    
    def _latex_(self):
        return r"\text{Affine Berkovich line over } \Bold{C}_{%s}" %(self.prime())

    Element = Berkovich_Element_Cp_Affine

class Berkovich_Cp_Projective(Berkovich):
    r"""
    The Berkovich Projective line over `\CC_p`.

    The Berkovich Projective line can be thought of as the one-point compactification
    of the Berkovich Affine line.

    INPUT:

     - base - The prime number `p`. Alternatively, can be a Projective Space over
       `\QQ_p` or a finite extension of `\QQ_p`. This allows for more control of conversion of centers.

    EXAMPLES::

        sage: B = Berkovich_Cp_Projective(3); B
        Projective Berkovich line over Cp(3) of precision 20

    Initializing by passing in a projective space looks the same::

        sage: S = ProjectiveSpace(Qp(3),1)
        sage: B = Berkovich_Cp_Projective(S); B
        Projective Berkovich line over Cp(3) of precision 20

    However, this method allows for more control over behind-the-scenes conversion::

        sage: S = ProjectiveSpace(Qp(3,1),1)
        sage: B = Berkovich_Cp_Projective(S); B
        Projective Berkovich line over Cp(3) of precision 1

        sage: Q1 = B(1/2); Q1
        Type I point centered at (2 + O(3) : 1 + O(3))

    Note that this point has very low precision, as S is a scheme over
    Qp(3) of capped-relative precision one.

    """
    def __init__(self,base):
        from sage.rings.integer_ring import ZZ
        if base in ZZ:
            if base.is_prime():
                base = ProjectiveSpace(Qp(base),1)
            else:
                raise ValueError("Non-prime pased into Berkovich space")
        else:
            from sage.schemes.projective.projective_space import is_ProjectiveSpace
            if not is_ProjectiveSpace(base):
                raise ValueError("Base of Projective Berkovich Space must be Projective Space")
            from sage.rings.padics.generic_nodes import is_pAdicField
            if not (is_pAdicField(base.base_ring())):
                raise ValueError("Base of Projective Berkovich Space must be Projective Space over Qp")
            if base.ambient_space() != base:
                raise ValueError("Base of Projective Berkovich Space must be Projective Space over Qp")
            if base.dimension() != 1:
                raise ValueError("Base of Projective Berkovich Space must be Projective Space of \
                    dimension 1 over Qp")
        self._p = base.base_ring().prime()
        Parent.__init__(self, base = base, category=TopologicalSpaces())

    def residue_characteristic(self):
        return self._p

    def prime(self):
        return self._p
        
    def _coerce_map_from_(self,S):
        if isinstance(S, Berkovich_Cp_Affine) or isinstance(S,Berkovich_Cp_Projective):
            if S.prime() == self.prime():
                return True
        return False
    
    def _repr_(self):
        return "Projective Berkovich line over Cp(%s) of precision %s" %(self.prime(),\
            self.base().base_ring().precision_cap())
    
    def _latex_(self):
        return r"\text{Projective Berkovich line over } \Bold{C}_{%s}" %(self.prime())

    Element = Berkovich_Element_Cp_Projective

    
