"""
Base class for Jacobians of curves
"""

#*******************************************************************************
#  Copyright (C) 2005 William Stein
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*******************************************************************************

from sage.categories.fields import Fields
_Fields = Fields()
from sage.schemes.generic.scheme import Scheme, is_Scheme

def is_Jacobian(J):
    """
    Return True if `J` is of type Jacobian_generic.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian, is_Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C)
        sage: is_Jacobian(J)
        True

    ::

        sage: E = EllipticCurve('37a1')
        sage: is_Jacobian(E)
        False
    """
    return isinstance(J, Jacobian_generic)

def Jacobian(C):
    """
    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: Jacobian(C)
        Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    try:
        return C.jacobian()
    except AttributeError:
        return Jacobian_generic(C)

class Jacobian_generic(Scheme):
    """
    Base class for Jacobians of projective curves.

    The input must be a projective curve over a field.

    EXAMPLES::

        sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
        sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
        sage: C = Curve(x^3 + y^3 + z^3)
        sage: J = Jacobian(C); J
        Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
    """
    def __init__(self, C):
        """
        TESTS::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian_generic
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: C = Curve(x^3 + y^3 + z^3)
            sage: J = Jacobian_generic(C); J
            Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: type(J)
            <class 'sage.schemes.jacobians.abstract_jacobian.Jacobian_generic_with_category'>

        Note: this is an abstract parent, so we skip element tests::

            sage: TestSuite(J).run(skip =["_test_an_element",\
                                          "_test_elements",\
                                          "_test_elements_eq_reflexive",\
                                          "_test_elements_eq_symmetric",\
                                          "_test_elements_eq_transitive",\
                                          "_test_elements_neq",\
                                          "_test_some_elements"])

        ::

            sage: Jacobian_generic(ZZ)
            Traceback (most recent call last):
            ...
            TypeError: Argument (=Integer Ring) must be a scheme.
            sage: Jacobian_generic(P2)
            Traceback (most recent call last):
            ...
            ValueError: C (=Projective Space of dimension 2 over Rational Field) must have dimension 1.
            sage: P2.<x, y, z> = ProjectiveSpace(Zmod(6), 2)
            sage: C = Curve(x + y + z)
            sage: Jacobian_generic(C)
            Traceback (most recent call last):
            ...
            TypeError: C (=Projective Plane Curve over Ring of integers modulo 6 defined by x + y + z) must be defined over a field.
        """
        if not is_Scheme(C):
            raise TypeError("Argument (=%s) must be a scheme."%C)
        if C.base_ring() not in _Fields:
            raise TypeError("C (=%s) must be defined over a field."%C)
        if C.dimension() != 1:
            raise ValueError("C (=%s) must have dimension 1."%C)
        self.__curve = C
        Scheme.__init__(self, C.base_scheme())

    def __cmp__(self, J):
        """
        Compare the Jacobian self to `J`.  If `J` is a Jacobian, then
        self and `J` are equal if and only if their curves are equal.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J1 = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J1 == J1
            True
            sage: J1 == P2
            False
            sage: J1 != P2
            True
            sage: J2 = Jacobian(Curve(x + y + z))
            sage: J1 == J2
            False
            sage: J1 != J2
            True
        """
        if not is_Jacobian(J):
            return cmp(type(self), type(J))
        return cmp(self.curve(), J.curve())

    def _repr_(self):
        """
        Return a string representation of this Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3)); J
            Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
            sage: J._repr_()
            'Jacobian of Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3'
        """
        return "Jacobian of %s"%self.__curve

    def _point(self):
        """
        Return the Hom-set from some affine scheme to ``self``.

        OUTPUT:

        This method always raises a ``NotImplementedError``; it is
        only abstract.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J._point()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def curve(self):
        """
        Return the curve of which self is the Jacobian.

        EXAMPLES::

            sage: from sage.schemes.jacobians.abstract_jacobian import Jacobian
            sage: P2.<x, y, z> = ProjectiveSpace(QQ, 2)
            sage: J = Jacobian(Curve(x^3 + y^3 + z^3))
            sage: J.curve()
            Projective Plane Curve over Rational Field defined by x^3 + y^3 + z^3
        """
        return self.__curve

    def change_ring(self, R):
        r"""
        Return the Jacobian over the ring `R`.

        INPUT:

        - ``R`` -- a field. The new base ring.

        OUTPUT:

        The Jacobian over the ring `R`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: Jac = H.jacobian();   Jac
            Jacobian of Hyperelliptic Curve over Rational
            Field defined by y^2 = x^3 - 10*x + 9
            sage: Jac.change_ring(RDF)
            Jacobian of Hyperelliptic Curve over Real Double
            Field defined by y^2 = x^3 - 10.0*x + 9.0
        """
        return self.curve().change_ring(R).jacobian()

    def base_extend(self, R):
        r"""
        Return the natural extension of ``self`` over `R`

        INPUT:

        - ``R`` -- a field. The new base field.

        OUTPUT:

        The Jacobian over the ring `R`.

        EXAMPLES::

            sage: R.<x> = QQ['x']
            sage: H = HyperellipticCurve(x^3-10*x+9)
            sage: Jac = H.jacobian();   Jac
            Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^3 - 10*x + 9
            sage: F.<a> = QQ.extension(x^2+1)
            sage: Jac.base_extend(F)
            Jacobian of Hyperelliptic Curve over Number Field in a with defining
            polynomial x^2 + 1 defined by y^2 = x^3 - 10*x + 9
        """
        if R not in _Fields:
            raise ValueError('Not a field: '+str(R))
        if self.base_ring() is R:
            return self
        if not R.has_coerce_map_from(self.base_ring()):
            raise ValueError('no natural map from the base ring (=%s) to R (=%s)!'
                             % (self.base_ring(), R))
        return self.change_ring(R)
