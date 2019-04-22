# -*- coding: utf-8 -*-
"""
Function Fields implemented with an expect interface to the optional
package kash.

This implementation is incomplete.  Many standard methods have not
yet been implemented.

A major feature of this implemention is its support for QQbar as
a base field.  Since kash only supports number fields, this is
implemented by keeping a working_constant_field, which is a number
field, and extending it as needed to include new algebraic numbers.
Most dependent objects keep a local variable kash_constant_field
and check it against the base function field's kash_constant_field
any time an operation is requested, recomputing all of the kash
objects if the kash_constant_field has changed.

EXAMPLES::

    # silences warnings about toy Buchberger implementation due to our use
    # of variety() and dimension() on ideals over QQbar
    sage: set_verbose(-1)                                     # optional - kash

    sage: F.<x> = FunctionField(QQ, implementation='kash')    # optional - kash
    sage: R.<Y> = F[]                                         # optional - kash
    sage: L.<y> = F.extension(Y^2 - x^8 - 1)                  # optional - kash
    sage: O = L.maximal_order()                               # optional - kash
    sage: I = O.ideal(x, y-1)                                 # optional - kash
    sage: P = I.place()                                       # optional - kash
    sage: D = P.divisor()                                     # optional - kash
    sage: D.basis_function_space()                            # optional - kash
    [1]
    sage: (2*D).basis_function_space()                        # optional - kash
    [1]
    sage: (3*D).basis_function_space()                        # optional - kash
    [1]
    sage: (4*D).basis_function_space()                        # optional - kash
    [1, 1/x^4*y + 1/x^4]


    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K) # optional - kash
    sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)            # optional - kash
    sage: O = F.maximal_order()                               # optional - kash
    sage: I = O.ideal(y)                                      # optional - kash
    sage: I.divisor()                                         # optional - kash
    2*Place (x, (1/(x^3 + x^2 + x))*y^2)
     + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)

    sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
    sage: L.<y> = K.extension(Y^2+Y+x+1/x)                    # optional - kash
    sage: O = L.maximal_order()                               # optional - kash
    sage: I = O.ideal(y)                                      # optional - kash
    sage: I.divisor()                                         # optional - kash
    - Place (x, x*y)
     + Place (x^2 + 1, x*y)

    sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
    sage: L.<y> = R[]                                         # optional - kash
    sage: F.<y> = R.extension(y^2 - (x^2+1))                  # optional - kash
    sage: D = (x/y).divisor()                                 # optional - kash
    sage: D                                                   # optional - kash
    - Place (x - I, y)
     + Place (x, y + x - 1)
     + Place (x, y + x + 1)
     - Place (x + I, y)
    sage: D2 = (x/y*x.differential()).divisor()               # optional - kash
    sage: pl1 = D.support()[0]                                # optional - kash
    sage: m = F.completion(pl1)                               # optional - kash
    sage: m(x/y)                                              # optional - kash
    1/2*s^-1 + 1/2*s + O(s^20)
    sage: m(x/y*x.differential())                             # optional - kash
    [2*I + 6*I*s^2 + 10*I*s^4 + 14*I*s^6 + 18*I*s^8 + 22*I*s^10 + 26*I*s^12 + 30*I*s^14 + 34*I*s^16 + 38*I*s^18 + O(s^20)] ds
    sage: pl2 = D.support()[1]                                # optional - kash
    sage: m2 = F.completion(pl2)                              # optional - kash
    sage: m2(x)                                               # optional - kash
    s + 1/4*s^3 + 1/16*s^5 + 1/64*s^7 + 1/256*s^9 + 1/1024*s^11 + 1/4096*s^13 + 1/16384*s^15 + 1/65536*s^17 + 1/262144*s^19 + O(s^20)
    sage: m2(x/y)                                             # optional - kash
    s - 1/4*s^3 + 1/16*s^5 - 1/64*s^7 + 1/256*s^9 - 1/1024*s^11 + 1/4096*s^13 - 1/16384*s^15 + 1/65536*s^17 - 1/262144*s^19 + O(s^20)
    sage: m2(x.differential())                                # optional - kash
    [1 + 3/4*s^2 + 5/16*s^4 + 7/64*s^6 + 9/256*s^8 + 11/1024*s^10 + 13/4096*s^12 + 15/16384*s^14 + 17/65536*s^16 + 19/262144*s^18 + O(s^20)] ds
    sage: m2(x * x.differential())                            # optional - kash
    [s + s^3 + 9/16*s^5 + 1/4*s^7 + 25/256*s^9 + 9/256*s^11 + 49/4096*s^13 + 1/256*s^15 + 81/65536*s^17 + 25/65536*s^19 + O(s^20)] ds
"""

import itertools
import operator

from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.misc.latex import latex

from sage.interfaces.kash import kash, KashElement

from sage.structure.richcmp import richcmp
from sage.structure.factorization import Factorization

from sage.rings.rational_field import QQ
from sage.rings.fraction_field import FractionField
from sage.rings.finite_rings.finite_field_base import FiniteField
from sage.rings.number_field.number_field_base import NumberField
from sage.rings.qqbar import QQbar, number_field_elements_from_algebraics
from sage.rings.ideal import Ideal

from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

from sage.symbolic.ring import SR

from .function_field import RationalFunctionField
from .function_field import FunctionFieldElement, FunctionFieldElement_rational
from .function_field import FunctionField_polymod
from .function_field import FunctionFieldElement_polymod
from .differential import DifferentialsSpace, FunctionFieldDifferential_global
from .order import FunctionFieldMaximalOrder
from .order import FunctionFieldMaximalOrderInfinite
from .ideal import FunctionFieldIdeal
from .place import FunctionFieldPlace
from .divisor import FunctionFieldDivisor
from .constructor import FunctionField
from .maps import FunctionFieldCompletion

class RationalFunctionField_kash(RationalFunctionField):
    """
    Rational function field `K(t)` in one variable, over the rationals,
    implemented using kash.

    EXAMPLES::

        sage: K.<t> = FunctionField(QQ, implementation='kash'); K # optional - kash
        Rational function field in t over Rational Field
        sage: K.gen()                                         # optional - kash
        t
        sage: 1/t + t^3 + 5                                   # optional - kash
        (t^4 + 5*t + 1)/t

    There are various ways to get at the underlying fields and rings
    associated to a rational function field::

        sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
        sage: K.base_field()                                  # optional - kash
        Rational function field in t over Rational Field
        sage: K.field()                                       # optional - kash
        Fraction Field of Univariate Polynomial Ring in t over Rational Field
        sage: K.constant_field()                              # optional - kash
        Rational Field
        sage: K.maximal_order()                               # optional - kash
        Maximal order of Rational function field in t over Rational Field

    We define a morphism::

        sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
        sage: L = FunctionField(QQ, 'tbar', implementation='kash') # give variable name as second input, optional - kash
        sage: K.hom(L.gen())                                  # optional - kash
        Function Field morphism:
          From: Rational function field in t over Rational Field
          To:   Rational function field in tbar over Rational Field
          Defn: t |--> tbar

    Here's a calculation over a number field::

        sage: R.<x> = FunctionField(QQ, implementation='kash') # optional - kash
        sage: L.<y> = R[]                                     # optional - kash
        sage: F.<y> = R.extension(y^2 - (x^2+1))              # optional - kash
        sage: (y/x).divisor()                                 # optional - kash
        - Place (x, y + x - 1)
         - Place (x, y + x + 1)
         + Place (x^2 + 1, y)

        sage: A.<z> = QQ[]                                    # optional - kash
        sage: NF.<i> = NumberField(z^2+1)                     # optional - kash
        sage: R.<x> = FunctionField(NF, implementation='kash') # optional - kash
        sage: L.<y> = R[]                                     # optional - kash
        sage: F.<y> = R.extension(y^2 - (x^2+1))              # optional - kash

        sage: (x/y*x.differential()).divisor()                # optional - kash
        -2*Place (1/x, 1/x*y + (-x + 1)/x)
         - 2*Place (1/x, 1/x*y + (x + 1)/x)
         + Place (x, y + x - 1)
         + Place (x, y + x + 1)

        sage: (x/y).divisor()                                 # optional - kash
        - Place (x - i, y)
         + Place (x, y + x - 1)
         + Place (x, y + x + 1)
         - Place (x + i, y)

    We try the same calculation over QQbar::

        sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
        sage: L.<y> = R[]                                     # optional - kash
        sage: F.<y> = R.extension(y^2 - (x^2+1))              # optional - kash
        sage: (y/x).divisor()                                 # optional - kash
        Place (x - I, y)
         - Place (x, y + x - 1)
         - Place (x, y + x + 1)
         + Place (x + I, y)

    """

    def __init__(self, constant_field, names, category=None):
        """
        Create a rational function field in one variable over the rational
        field.

        INPUT:

        - ``constant_field`` -- QQ

        - ``names`` -- string or tuple of length 1

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K # optional - kash
            Rational function field in t over Rational Field
            sage: K.category()                                # optional - kash
            Category of function fields
            sage: FunctionField(QQ[I], 'alpha', implementation='kash')   # not implemented, optional - kash
            Rational function field in alpha over Number Field in I with defining polynomial x^2 + 1

        Must be over a field::

            sage: FunctionField(ZZ, 't', implementation='kash') # optional - kash
            Traceback (most recent call last):
            ...
            TypeError: constant_field must be a field
        """

        RationalFunctionField.__init__(self, constant_field, names, category)

        # the field, without the function field structure (referenced by divisor.py)
        self._field = constant_field[names[0]].fraction_field()

        if constant_field is QQbar:
            self._set_working_constant_field(QQ)
        else:
            self._set_working_constant_field(constant_field)

    def _set_working_constant_field(self, constant_field):

        self.working_constant_field = constant_field

        if constant_field is QQ:
            self.kash_constant_field = kash.RationalField()
            self.reverse_map = {}
        elif isinstance(constant_field, FiniteField):
            assert constant_field.order().is_prime()
            self.kash_constant_field = kash.GaloisField(constant_field.order())
            self.reverse_map = {}
        elif isinstance(constant_field, NumberField):
            kZa = kash.IntegerRing().PolynomialAlgebra()
            ka = kZa.Element(list(constant_field.defining_polynomial()))
            self.kash_constant_field = ka.NumberField()
            # If constant_field is equipped with an embedding (almost surely to QQbar),
            # we map the kash generator to that target.  Otherwise, we map the kash
            # generator to the number field generator.
            if constant_field.gen_embedding():
                self.reverse_map = {self.kash_constant_field.gen(1) : constant_field.gen_embedding()}
            else:
                self.reverse_map = {self.kash_constant_field.gen(1) : constant_field.gen(0)}
        else:
            raise ValueError("The constant field must be either QQ, QQbar, a number field, or a finite field.")

        # we seem to need this to avoid getting variable names like '$.1' Kash's output
        self.kash_constant_field.PolynomialAlgebra().AssignNames_(['"x"'])

        self._kash_ = self.kash_constant_field.RationalFunctionField()

    def _extend_constant_field(self, arg):
        """
        Extend the constant field to express `arg`.
        (only for function fields over QQbar)

        Currently, `arg` can only be a list of polynomials or a FunctionFieldDivisor.

        TESTS::

        basis_function_space() creates an ideal(1), which broke this code before it
        checked if g was a FunctionFieldElement before calling g.parent()

            sage: F.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: R.<Y> = F[]                                 # optional - kash
            sage: L.<y> = F.extension(Y^2 - x^8 - 1)          # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(x, y-1)                         # optional - kash
            sage: P = I.place()                               # optional - kash
            sage: D = P.divisor()                             # optional - kash
            sage: D.basis_function_space()                    # optional - kash
            [1]

        """

        algebraics = []
        if self.working_constant_field != QQ:
            algebraics.append(self.working_constant_field.gen())

        if isinstance(arg, FunctionFieldDivisor):
            ideals = [pls.prime_ideal().gens() for pls in arg.support()]
        else:
            ideals = [arg]

        # We've now got a list of ideals, each specified as a list of
        # generators.  For zero-dimensional ideals, we find their
        # solution points and extend working_constant_field enough to
        # express their coordinates.  For higher dimensional ideals,
        # we currently do nothing, though I'm not sure that's right.

        # Function field ideals are not polynomial ideals, so we can't
        # just form the generators into an ideal; they have to be
        # converted into a polynomial ring, using our numerator()
        # function.

        # _to_bivariate_polynomial() doesn't work on polynomials over
        # the underlying function field, only on polynomials in the
        # extension.  Also, ideal(1) needs to be converted to
        # ideal(Integer(1)), or it will fail.

        from sage.rings.integer import Integer

        def numerator(x):
            if isinstance(x, int):
                return Integer(x)
            elif not isinstance(x, FunctionFieldElement):
                return x
            elif x.parent() is self:
                return x.numerator()
            else:
                return self._to_bivariate_polynomial(x)[0]

        from sage.rings.polynomial.multi_polynomial_ring_base import MPolynomialRing_base
        from sage.rings.polynomial.polynomial_ring import PolynomialRing_general

        for I in ideals:
            I = Ideal(map(numerator, I))
            if isinstance(I.ring(), PolynomialRing_general):
                # univariate case - ideals over univariate rings don't
                # implement dimension() or variety()
                algebraics.extend(I.gen().roots(multiplicities=False))
            elif isinstance(I.ring(), MPolynomialRing_base):
                if I.dimension() == 0:
                    for point in I.variety():
                        algebraics.extend(point.values())

        (constant_field, new_algebraics, nftoQQbar) = number_field_elements_from_algebraics(algebraics)

        # If we decided to expand our number field, then duplicate constant_field,
        # but with an embedding attached, and set it as our new working constant field.

        if constant_field.degree() != self.working_constant_field.degree():
            import sage.rings.number_field.number_field as number_field
            constant_field = number_field.NumberField(constant_field.polynomial(), constant_field.gen(),
                                         embedding = nftoQQbar(constant_field.gen()))
            self._set_working_constant_field(constant_field)

    def kash(self):
        return self._kash_

    def to_kash(self, c):
        x = self.kash().gen(1)
        # c will be an element of a fraction field with polynomial coefficients
        c = c.element()
        # This isn't just an optimization.  Kash won't build multivariate polynomials
        # if a constant coefficient is in a polynomial ring / fraction field.  The
        # main line code below does that; this code makes sure that constants
        # are presented as constants, and not as zero-degree polynomials.
        # XXX probably still has problems with algebraic constants like sqrt(2),
        # which fall through into the code below
        try:
            if self.constant_field() is QQbar:
                return QQ(c)
            else:
                return self.constant_field()(c)
        except:
            pass
        # Sage considers QQ to be a number field, which is why I check
        # to see if there's a reverse_map instead of just checking
        # to see if self.constant_field() is a NumberField.
        if len(self.reverse_map) > 0:
            # this is the number field case
            if c.numerator().parent().base_ring() is QQbar:
                # XXX assumes that the polynomial coefficients are in working_constant_field
                # If they're not, we should extend working_constant_field, but instead
                # we just throw an exception.  This isn't a problem if we work with
                # divisors and ideals generated by the code itself, but if the user tries to
                # create an ideal himself with a previously unseen algebraic number,
                # we'll get an exception here.
                n = c.numerator().change_ring(self.working_constant_field)
                d = c.denominator().change_ring(self.working_constant_field)
            else:
                n = c.numerator()
                d = c.denominator()
            # P will be a polynomial ring with coeffs in a number field
            P = c.numerator().parent()
            # create a bivariate ring, one variable for the polynomial var, and one for the number field generator
            R2 = QQ[P.gens() + ('a',)]
            # the number field generator will map to 'a'
            a = R2.gen(1)
            # for each coefficient, convert it to a polynomial in the number field generator, then map it into R2
            n = n.map_coefficients(lambda v: v.polynomial().subs({v.polynomial().parent().gen(0) : a}), new_base_ring = R2)
            # now map n's polynomial variable to R2's polynomial variable (n is now in R2)
            n = n.subs({n.parent().gen(0): R2.gen(0)})
            # now map n into kash
            n = n.subs({R2.gen(0) : x, R2.gen(1) : self.kash_constant_field.gen(1)})
            # ditto for the denominator
            d = d.map_coefficients(lambda v: v.polynomial().subs({v.polynomial().parent().gen(0) : a}), new_base_ring = R2)
            d = d.subs({d.parent().gen(0): R2.gen(0)})
            d = d.subs({R2.gen(0) : x, R2.gen(1) : self.kash_constant_field.gen(1)})
            return n/d
        else:
            # P will be a polynomial ring with coeffs in either QQ or QQbar
            # XXX this code assumes that even if we're in QQbar, actual coeffs are only in QQ
            # (same issue as the comment above)
            P = c.numerator().parent()
            R2 = QQ[P.gens()]
            n = c.numerator().change_ring(QQ)(x)
            d = c.denominator().change_ring(QQ)(x)
            return n/d

    def extension(self, f, names=None):
        """
        Create an extension `L = K[y]/(f(y))` of the rational function field.

        INPUT:

        - ``f`` -- univariate polynomial over self

        - ``names`` -- string or length-1 tuple

        OUTPUT:

        - a function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: K.extension(y^5 - x^3 - 3*x + x*y)          # optional - kash
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        A nonintegral defining polynomial::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: K.extension(y^3 + (1/t)*y + t^3/(t+1))      # optional - kash
            Function field in y defined by y^3 + 1/t*y + t^3/(t + 1)

        The defining polynomial need not be monic or integral::

            sage: K.extension(t*y^3 + (1/t)*y + t^3/(t+1))    # optional - kash
            Function field in y defined by t*y^3 + 1/t*y + t^3/(t + 1)
        """
        from . import constructor
        return constructor.FunctionFieldExtension(f, names, implementation='kash')

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: K.place_set()                               # optional - kash
            Set of places of Rational function field in t over Rational Field
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5), implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2)) # optional - kash
            sage: L.space_of_differentials()                  # optional - kash
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

            sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: K.space_of_differentials()                  # optional - kash
            Space of differentials of Rational function field in t over Rational Field
        """
        return DifferentialsSpace_kash(self)

def divisor_decorator(f):

    def wrapper(self):

        if self.is_zero():
            raise ValueError("divisor not defined for zero")

        # if we're working over QQbar, this divisor is over the current constant field
        D = getattr(super(FunctionFieldElement_polymod_kash, self), f.__name__)()

        if self.parent().constant_base_field() is not QQbar:
            return D

        # QQbar case - compute a number field that can factor all of our support
        # polynomials in the base function field C(x).

        # If we decided to expand our number field, then redo the divisor
        # computation in a new function field

        self.parent().base_field()._extend_constant_field(D)

        D = getattr(super(FunctionFieldElement_polymod_kash, self), f.__name__)()

        return D

    return wrapper

class FunctionFieldDifferential_kash(FunctionFieldDifferential_global):

    def divisor_of_poles(self):
        """
        Return the divisor of poles of the differential.

        NOTE::
            To avoid extending the constant field unnecessarily, we first
            compute the full divisor in the current constant field, then
            discard any positive places to get the divisor of poles,
            then extend the constant field so that all of of those places
            factor, then repeat the calculation for the final result.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y>=K[] # optional - kash
            sage: L.<y> = K.extension(Y^3+x+x^3*Y)            # optional - kash
            sage: w = (1/y) * y.differential()                # optional - kash
            sage: w.divisor_of_poles()                        # optional - kash
            Place (1/x, -1/x^3*y^2 + 1/x)
             + Place (1/x, 1/x^3*y^2 + 1/x^2*y + 1)
             + Place (x, y)

            sage: R.<x> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: root = x^4+4*x^3+2*x^2+1                    # optional - kash
            sage: F.<y> = R.extension(y^2 - root)             # optional - kash
            sage: num = 6*x^2 + 5*x +7                        # optional - kash
            sage: den = 2*x^6 + 8*x^5 + 3*x^4 + - 4*x^3 - 1   # optional - kash
            sage: integrand = y*num/den * x.differential();   # optional - kash
            sage: integrand.divisor_of_poles()                # optional - kash
            Place (x^2 - 1/2, y - 2*x - 1/2)
             + Place (x^2 - 1/2, y + 2*x + 1/2)

            sage: QQbar.options.display_format = 'radical'    # optional - kash
            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: root = x^4+4*x^3+2*x^2+1                    # optional - kash
            sage: F.<y> = R.extension(y^2 - root)             # optional - kash
            sage: num = 6*x^2 + 5*x +7                        # optional - kash
            sage: den = 2*x^6 + 8*x^5 + 3*x^4 + - 4*x^3 - 1   # optional - kash
            sage: integrand = y*num/den * x.differential();   # optional - kash
            sage: integrand.divisor_of_poles()                # optional - kash
            Place (x - 1/2*sqrt(2), y - sqrt(2) - 1/2)
             + Place (x - 1/2*sqrt(2), y + sqrt(2) + 1/2)
             + Place (x + 1/2*sqrt(2), y - sqrt(2) + 1/2)
             + Place (x + 1/2*sqrt(2), y + sqrt(2) - 1/2)

        """
        from .divisor import FunctionFieldDivisor
        F = self.base_ring()
        x = F.base_field().gen()
        D = self._f._divisor() + (-2) * F(x)._divisor_of_poles() + F.different()
        D = FunctionFieldDivisor(D.parent().function_field(), {p:-m for p,m in D.list() if m < 0})
        if F.constant_base_field() is QQbar:
            F.base_field()._extend_constant_field(D)
            D = self._f._divisor() + (-2) * F(x)._divisor_of_poles() + F.different()
            D = FunctionFieldDivisor(D.parent().function_field(), {p:-m for p,m in D.list() if m < 0})
        return D

    def residue(self, place):
        return self._field.completion(place, prec=0)(self)[-1]

class DifferentialsSpace_kash(DifferentialsSpace):
    """
    Space of differentials of a function field.
    """
    Element = FunctionFieldDifferential_kash

class FunctionFieldElement_polymod_kash(FunctionFieldElement_polymod):

    def valuation(self, place):
        """
        Return the valuation of the element at the place.

        INPUT:

        - ``place`` -- a place of the function field

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 - (x^2+1))          # optional - kash
            sage: pl = F.maximal_order().ideal(x-QQbar(sqrt(-1)),y).place() # optional - kash
            sage: pl                                          # optional - kash
            Place (x - I, y)
            sage: x.valuation(pl)                             # optional - kash
            0
            sage: y.valuation(pl)                             # optional - kash
            1
            sage: (x-QQbar(sqrt(-1))).valuation(pl)           # optional - kash
            2
        """
        prime = place.prime_ideal()
        ideal = prime.ring().ideal(self)
        return prime.valuation(ideal)

    @divisor_decorator
    def divisor(self):
        pass

    @divisor_decorator
    def divisor_of_zeros(self):
        pass

    @divisor_decorator
    def divisor_of_poles(self):
        pass

    def _divisor(self):
        """
        Return the divisor of the element, without extending the constant subfield.
        """
        return super(FunctionFieldElement_polymod_kash, self).divisor()

    def _divisor_of_poles(self):
        """
        Return the divisor of poles of the element, without extending the constant
        subfield.
        """
        return super(FunctionFieldElement_polymod_kash, self).divisor_of_poles()


class FunctionField_polymod_kash(FunctionField_polymod):
    """
    Function fields defined by a univariate polynomial, as an extension of the
    base field, implemented using kash.

    EXAMPLES:

    We make a function field defined by a degree 5 polynomial over the
    rational function field over the rational numbers::

        sage: K.<x> = FunctionField(QQ, implementation='kash') # optional - kash
        sage: R.<y> = K[]                                     # optional - kash
        sage: L.<y> = K.extension(y^5 - (x^3 + 2*x*y + 1/x)); L # optional - kash
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x

    We next make a function field over the above nontrivial function
    field L::

        sage: S.<z> = L[]                                     # optional - kash
        sage: M.<z> = L.extension(z^2 + y*z + y); M           # optional - kash
        Function field in z defined by z^2 + y*z + y
        sage: 1/z                                             # optional - kash
        ((-x/(x^4 + 1))*y^4 + 2*x^2/(x^4 + 1))*z - 1
        sage: z * (1/z)                                       # optional - kash
        1

    We drill down the tower of function fields::

        sage: M.base_field()                                  # optional - kash
        Function field in y defined by y^5 - 2*x*y + (-x^4 - 1)/x
        sage: M.base_field().base_field()                     # optional - kash
        Rational function field in x over Rational Field
        sage: M.base_field().base_field().constant_field()    # optional - kash
        Rational Field
        sage: M.constant_base_field()                         # optional - kash
        Rational Field

    The polynomial must be irreducible::

        sage: K.<x>=FunctionField(QQ, implementation='kash')  # optional - kash
        sage: R.<y> = K[]                                     # optional - kash
        sage: L.<y>=K.extension(x^2-y^2)                      # optional - kash
        Traceback (most recent call last):
        ...
        ValueError: The polynomial must be irreducible
    """
    Element = FunctionFieldElement_polymod_kash

    def __init__(self, polynomial, names, category=None):
        """
        Create a function field defined as an extension of another function
        field by adjoining a root of a univariate polynomial.

        INPUT:

        - ``polynomial`` -- univariate polynomial over a function field

        - ``names`` -- tuple of length 1 or string; variable names

        - ``category`` -- category (default: category of function fields)

        EXAMPLES:

        We create an extension of a function field::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: L = K.extension(y^5 - x^3 - 3*x + x*y); L   # optional - kash
            Function field in y defined by y^5 + x*y - x^3 - 3*x

        Note the type::

            sage: type(L)                                     # optional - kash
            <class 'sage.rings.function_field.function_field_kash.FunctionField_polymod_kash_with_category'>

        We can set the variable name, which doesn't have to be y::

            sage: L.<w> = K.extension(y^5 - x^3 - 3*x + x*y); L # optional - kash
            Function field in w defined by w^5 + x*w - x^3 - 3*x

        TESTS:

        Test that :trac:`17033` is fixed::

            sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: R.<x> = QQ[]                                # optional - kash
            sage: M.<z> = K.extension(x^7-x-t)                # optional - kash
            sage: M(x)                                        # optional - kash
            z
            sage: M('z')                                      # optional - kash
            z
            sage: M('x')                                      # optional - kash
            Traceback (most recent call last):
            ...
            TypeError: unable to evaluate 'x' in Fraction Field of Univariate
            Polynomial Ring in t over Rational Field
        """

        FunctionField_polymod.__init__(self, polynomial, names, category)

        assert isinstance(polynomial.base_ring(), RationalFunctionField_kash)

        # the field, without the function field structure (referenced by divisor.py)

        # don't set these here; they were set in the superclass constructor
        # self._base_field = polynomial.base_ring()
        # self._polynomial = polynomial

        # self._field = self._base_field[self._gen].extension(polynomial)
        # self._field = self._base_field[polynomial.parent().gen(0)].extension(polynomial)
        self._field = self

        self._place_class = FunctionFieldPlace_kash

        self.kash_constant_field = None

        # compute the kash field, so we check for irreducibility
        # XXX for QQbar, need to check for irreducibility in Sage, over QQbar
        self.kash()

    def kash(self):

        if self.kash_constant_field != self.base_field().kash_constant_field:

            self.kash_constant_field = self.base_field().kash_constant_field

            polynomial = self.polynomial()

            kash_base_field = self.base_field().kash()
            kTy = kash_base_field.PolynomialAlgebra()

            x = kash_base_field.gen(1)
            y = kTy.gen(1)

            # Goofiness because 'subs' doesn't recurse into coefficient
            # rings, so we convert our polynomial to a list of
            # coefficients, map them to Kash, and pass the list to Kash's
            # Element constructor, which builds a polynomial from a list.

            # Also, FunctionFieldElement's aren't callable, and thus the
            # extra call to element() to convert the FunctionFieldElement
            # to a FractionFieldElement

            #kash_polynomial = kTy.Element(map(lambda c: c.element()(x), list(polynomial)))
            kash_polynomial = kTy.Element([self.base_field().to_kash(c) for c in polynomial])

            try:
                self._kash_ = kash_polynomial.FunctionField()
            except TypeError as ex:
                if 'Polynomial must be irreducible' in ex.args[0]:
                    raise ValueError("The polynomial must be irreducible")
                raise ex

            ybar = self._kash_.gen(1)

            # Set up reverse map to convert elements back from kash

            self.reverse_map = {x : polynomial.base_ring().gen(0), ybar : self.gen(0)}
            self.reverse_map.update(polynomial.base_ring().reverse_map)

        return self._kash_

    def to_kash(self, polynomial):

        # Same kind of goofiness as above.  When building an element
        # of an algebraic extension, Kash requires the list to be the
        # same length as the degree of the extension.  Also, Kash's
        # Element function doesn't accept constant polynomials - they
        # have to actually be constants.

        x = self.base_ring().kash().gen(1)
        coeffs = list(self(polynomial).element())
        coeffs += [self.base_ring().zero()] * (self.degree() - len(coeffs))
        coeffs = [self.base_ring().to_kash(c) for c in coeffs]
        return self.kash().Element(coeffs)

    @cached_method
    def place_set(self):
        """
        Return the set of all places of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: L.place_set()                               # optional - kash
            Set of places of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        from .place import PlaceSet
        return PlaceSet(self)

    @cached_method
    def space_of_differentials(self):
        """
        Return the space of differentials attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(5), implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2)) # optional - kash
            sage: L.space_of_differentials()                  # optional - kash
            Space of differentials of Function field in y defined by y^3 + (4*x^3 + 1)/(x^3 + 3)

            sage: K.<t> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: K.space_of_differentials()                  # optional - kash
            Space of differentials of Rational function field in t over Rational Field
        """
        return DifferentialsSpace_kash(self)

    @cached_method
    def divisor_group(self):
        """
        Return the group of divisors attached to the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^3 - (x^3 - 1)/(x^3 - 2)) # optional - kash
            sage: L.divisor_group()                           # optional - kash
            Divisor group of Function field in y defined by y^3 + (-x^3 + 1)/(x^3 - 2)
        """
        from .divisor import DivisorGroup
        return DivisorGroup(self)

    @cached_method
    def maximal_order(self):
        """
        Return the maximal order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); # optional - kash
            sage: R.<t> = PolynomialRing(K);                  # optional - kash
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18); # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: O.basis()                                   # optional - kash
            (1, 1/x^4*y, 1/x^9*y^2, 1/x^13*y^3)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); # optional - kash
            sage: R.<t> = PolynomialRing(K);                  # optional - kash
            sage: F.<y> = K.extension(t^4 + x^12*t^2 + x^18*t + x^21 + x^18); # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: O.basis()                                   # optional - kash
            (1, 1/x^4*y, 1/x^11*y^2 + 1/x^2, 1/x^15*y^3 + 1/x^6*y)

        The basis of the maximal order *always* starts with 1. This is assumed
        in some algorithms.
        """
        return FunctionFieldMaximalOrder_kash(self)

    @cached_method
    def maximal_order_infinite(self):
        """
        Return the maximal infinite order of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[] # optional - kash
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2) # optional - kash
            sage: F.maximal_order_infinite()                  # optional - kash
            Maximal infinite order of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: L.maximal_order_infinite()                  # optional - kash
            Maximal infinite order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """
        return FunctionFieldMaximalOrderInfinite_kash(self)

    def different(self):
        """
        Return the different divisor of the function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); R.<t> = PolynomialRing(K) # optional - kash
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)    # optional - kash
            sage: F.different()                               # optional - kash
            2*Place (x, (1/(x^3 + x^2 + x))*y^2)
             + 2*Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)
        """
        D = self.kash().DifferentDivisor()
        support = D.Support()
        data = {FunctionFieldPlace_kash(self, place) : D.Valuation(place).sage() for place in support}
        return FunctionFieldDivisor(self.divisor_group(), data)

    def completion(self, place, name=None, prec=None, uvar=None):
        """
        Return the completion of the function field at the place.

        INPUT:

        - ``place`` -- place

        - ``name`` -- string; name of the series variable

        - ``prec`` -- integer; default precision

        - ``uvar`` -- uniformizing variable.  If ``None``, then an
          arbitrary uniformizing variable is selected.  Must be either
          an element of the function field, or an element of the
          Symbolic Ring constructed from a root of such an element.
          In either case, the specified variable is checked to see
          if it has valuation one at the specified place, but in
          the later case (element of the Symbolic Ring), no check
          is made to see if such a variable actually exists in the
          function field.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 - (x^2+1))          # optional - kash
            sage: D = (y/x).divisor()                         # optional - kash
            sage: p = D.support()[0]                          # optional - kash
            sage: m = F.completion(p)                         # optional - kash
            sage: m                                           # optional - kash
            Completion map:
              From: Function field in y defined by y^2 - x^2 - 1
              To:   Laurent Series Ring in s over Algebraic Field
            sage: m(x, 10)                                    # optional - kash
            I + 2*I*s^2 + 2*I*s^4 + 2*I*s^6 + 2*I*s^8 + O(s^10)
            sage: m(y, 10)                                    # optional - kash
            2*I*s + 2*I*s^3 + 2*I*s^5 + 2*I*s^7 + 2*I*s^9 + O(s^10)

            sage: I = sqrt(QQbar(-1))                         # optional - kash
            sage: QQbar.options.display_format = 'radical'    # optional - kash
            sage: m2 = F.completion(p, uvar=(x-I))            # optional - kash
            Traceback (most recent call last):
            ...
            ValueError: x - I is not a uniformizing variable at Place (x - I, y)
            sage: m2 = F.completion(p, uvar=SR(x-I)^(1/3))    # optional - kash
            Traceback (most recent call last):
            ...
            ValueError: (x - I)^(1/3) is not a uniformizing variable at Place (x - I, y)
            sage: m2 = F.completion(p, uvar=sqrt(x-I))        # optional - kash
            sage: m2(x, 10)                                   # optional - kash
            I + s^2
            sage: m2(y, 10)                                   # optional - kash
            (I + 1)*s + (-1/4*I + 1/4)*s^3 + (1/32*I + 1/32)*s^5 + (1/128*I - 1/128)*s^7 + (-5/2048*I - 5/2048)*s^9 + O(s^10)

        .. TODO:
            Remove the need to explicitly cast into SR when constructing
            roots of uniformizing variables.
        """
        return FunctionFieldCompletion_kash(self, place, name=name, prec=prec, uvar=uvar)

    def places_infinite(self, degree=1):
        """
        Return a list of the infinite places of degree ``degree``.

        INPUT:

        - ``degree`` -- positive integer (default: `1`)

        EXAMPLES::

            sage: F.<a>=GF(2)                                 # optional - kash
            sage: K.<x>=FunctionField(F, implementation='kash') # optional - kash
            sage: R.<t>=PolynomialRing(K)                     # optional - kash
            sage: L.<y>=K.extension(t^4+t-x^5)                # optional - kash
            sage: L.places_infinite(1)                        # optional - kash
            [Place (1/x, 1/x^4*y^3)]
        """
        return [place for place in self._places_infinite(degree)]

    def _places_infinite(self, degree):
        """
        Return a generator of *infinite* places of the function field of the degree.

        INPUT:

        - ``degree`` -- positive integer

        EXAMPLES::

            sage: F.<a>=GF(2)                                 # optional - kash
            sage: K.<x>=FunctionField(F, implementation='kash') # optional - kash
            sage: R.<t>=PolynomialRing(K)                     # optional - kash
            sage: L.<y>=K.extension(t^4+t-x^5)                # optional - kash
            sage: L._places_infinite(1)                       # optional - kash
            <generator object ...>
        """
        Oinf = self.maximal_order_infinite()
        for prime,_,_ in Oinf.decomposition():
            place = prime.place()
            if place.degree() == degree:
                yield place

from sage.rings.laurent_series_ring import LaurentSeriesRing
from sage.rings.laurent_series_ring_element import LaurentSeries
from sage.rings.infinity import Infinity

class LaurentSeriesDifferential(LaurentSeries):
    """
    A slight variant on LaurentSeries that prints the series in brackets
    followed by a differential.

    No attempt is made to enforce consistency of operations, i.e, you
    can just add this to a normal LaurentSeries and the differential will
    disappear.
    """
    def _repr_(self):
        differential = "d" + self.variable()
        if self == 1:
            return differential
        if (self.prec() == Infinity) and (len(self.coefficients()) == 1):
            return super(LaurentSeriesDifferential, self)._repr_() + " " + differential
        return "[" + super(LaurentSeriesDifferential, self)._repr_() + "] " + differential

    def _latex_(self, **kwds):
        differential = "d" + latex(self.parent().gen(0))
        if self == 1:
            return differential
        if (self.prec() == Infinity) and (len(self.coefficients()) == 1):
            return super(LaurentSeriesDifferential, self)._latex_(**kwds) + "\\," + differential

        # Sizing the brackets is problematic if there are alignment
        # characters in the LaTeX, because we can't match the brackets
        # across multiple cells in a table.  In this case, we put a
        # non-aligned version of the LaTeX in a "vphantom" to size
        # the brackets correctly.
        #
        # https://www.giss.nasa.gov/tools/latex/ltx-403.html

        latex_series = super(LaurentSeriesDifferential, self)._latex_(**kwds)
        if kwds.pop('table', None):
            latex_series_no_alignment = super(LaurentSeriesDifferential, self)._latex_(**kwds)
        else:
            latex_series_no_alignment = latex_series_no_alignment
        if latex_series == latex_series_no_alignment:
            return "\\left[" + latex_series + "\\right]\\," + differential
        else:
            return "\\left[\\vphantom{{{0}}}\\right.{1}\\left.\\vphantom{{{0}}}\\right]\\,".format(latex_series_no_alignment, latex_series) + differential

class FunctionFieldCompletion_kash(FunctionFieldCompletion):
    """
    Completions on kash function fields.  Currently only supports
    QQbar as the field of constants.

    EXAMPLES::

        sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
        sage: L.<y> = R[]                                     # optional - kash
        sage: F.<y> = R.extension(y^2 - (x^2+1))              # optional - kash
        sage: D = (y/x).divisor()                             # optional - kash
        sage: p = D.support()[0]                              # optional - kash
        sage: m = F.completion(p)                             # optional - kash
        sage: m                                               # optional - kash
        Completion map:
          From: Function field in y defined by y^2 - x^2 - 1
          To:   Laurent Series Ring in s over Algebraic Field
        sage: m(x)                                            # optional - kash
        I + 2*I*s^2 + 2*I*s^4 + 2*I*s^6 + 2*I*s^8 + 2*I*s^10 + 2*I*s^12 + 2*I*s^14 + 2*I*s^16 + 2*I*s^18 + O(s^20)
        sage: m(y)                                            # optional - kash
        2*I*s + 2*I*s^3 + 2*I*s^5 + 2*I*s^7 + 2*I*s^9 + 2*I*s^11 + 2*I*s^13 + 2*I*s^15 + 2*I*s^17 + 2*I*s^19 + O(s^20)
        sage: m(x*y) == m(x) * m(y)                           # optional - kash
        True
        sage: m(x+y) == m(x) + m(y)                           # optional - kash
        True
        sage: m(y)^2 == m(x)^2 + 1                            # optional - kash
        True

    The variable name of the series can be supplied, as can the default precision.

        sage: p2 = D.support()[1]                             # optional - kash
        sage: p2                                              # optional - kash
        Place (x, y + x - 1)
        sage: m2 = F.completion(p2, 't', prec=10)             # optional - kash
        sage: m2(x)                                           # optional - kash
        t + 1/4*t^3 + 1/16*t^5 + 1/64*t^7 + 1/256*t^9 + O(t^10)
        sage: m2(y)                                           # optional - kash
        1 + 1/2*t^2 + 1/8*t^4 + 1/32*t^6 + 1/128*t^8 + O(t^10)
    """
    def __init__(self, field, place, name=None, prec=None, uvar=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        - ``place`` -- place of the function field

        - ``name`` -- string for the name of the series variable

        - ``prec`` -- positive integer; default precision

        - ``uvar`` -- uniformizing variable.  If ``None``, then an
          arbitrary uniformizing variable is selected.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 - (x^2+1))          # optional - kash
            sage: D = (y/x).divisor()                         # optional - kash
            sage: p = D.support()[0]                          # optional - kash
            sage: m = F.completion(p)                         # optional - kash
            sage: m                                           # optional - kash
            Completion map:
              From: Function field in y defined by y^2 - x^2 - 1
              To:   Laurent Series Ring in s over Algebraic Field

        """

        if name is None:
            name = 's' # default

        # Currently we only work on places of degree one, where our residue field
        # is simply the constant base field.  In particular, this always works over QQbar.

        if place.degree() != 1:
            raise NotImplementedError("series expansions not implemented at places of degree > 1")

        # if prec is None, the Laurent series ring provides default
        # precision
        codomain = LaurentSeriesRing(field.constant_base_field(), name=name, default_prec=prec)

        FunctionFieldCompletion.__init__(self, field, codomain)

        if uvar is not None:
            if uvar.parent() is SR:
                if uvar.operator() is operator.pow:
                    u,p = uvar.operands()
                    u = field(u)
                else:
                    u = field(uvar)
                    p = 1
            else:
                u = uvar
                p = 1
            if u.valuation(place) != 1/p:
                raise ValueError("{} is not a uniformizing variable at {}".format(uvar, place))

        self._place = place
        self._precision = codomain.default_prec()
        self._uvar = uvar

    def _call_(self, f):
        """
        Call the completion for f

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)      # optional - kash
            sage: D = (x*y).divisor()                         # optional - kash
            sage: p = D.support()[2]                          # optional - kash
            sage: m = F.completion(p, prec=10)                # optional - kash
            sage: m                                           # optional - kash
            Completion map:
              From: Function field in y defined by y^2 + y + (x^2 + 1)/x
              To:   Laurent Series Ring in s over Algebraic Field
            sage: m(x)                                        # optional - kash
            -s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        return self._expand(f, prec=None)

    def _call_with_args(self, f, args=(), kwds={}):
        """
        Call the completion with ``args`` and ``kwds``.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)      # optional - kash
            sage: D = (x*y).divisor()                         # optional - kash
            sage: p = D.support()[2]                          # optional - kash
            sage: m = F.completion(p)                         # optional - kash
            sage: m(x+y, 10)                                  # indirect doctest, optional - kash
            -s^-1 - s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        return self._expand(f, *args, **kwds)

    def _expand(self, f, prec=None, uvar=None):
        """Return the power series representation of f of precision prec.

        INPUT:

        - ``f`` -- element of the function field

        - ``prec`` -- integer; absolute precision of the series

        - ``uvar`` -- uniformizing variable.  If ``None``, then the
          uniformizing variable specified with the completion is
          selected, or if no uniformizing variable was specified with
          the completion, then an arbitrary one is selected.
          ``False`` is used internally to force selection of an
          arbitrary uniformizing variable, even if one was specified
          when the completion was created.

        OUTPUT:

        - a series of precision prec

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)      # optional - kash
            sage: D = (x*y).divisor()                         # optional - kash
            sage: p = D.support()[2]                          # optional - kash
            sage: m = F.completion(p)                         # optional - kash
            sage: m(x, 10)                                    # indirect doctest, optional - kash
            -s^2 - s^3 - s^4 - s^5 - 2*s^6 - 4*s^7 - 7*s^8 - 11*s^9 + O(s^10)
        """
        if prec is None:
            prec = self._precision
        if uvar is None:
            uvar = self._uvar

        place = self._place
        F = place.function_field()

        kash_series = F.to_kash(f).Expand(place.kash(), AbsPrec=prec)

        val = kash_series.Valuation()
        coeffs = [c.sage(F.reverse_map) for c in list(kash_series.Coefficients())]

        s_series = self.codomain()(coeffs, val).add_bigoh(prec)

        if uvar is not None and uvar is not False:
            if uvar.parent() is SR:
                if uvar.operator() is operator.pow:
                    u,p = uvar.operands()
                    u = F(u)
                else:
                    u = F(uvar)
                    p = 1
            else:
                u = uvar
                p = 1
            if u.valuation(place) != 1/p:
                raise ValueError("{} is not a uniformizing variable at {}".format(uvar, place))
            v = s_series.valuation()
            t_series = s_series((self(u, uvar=False, prec=int((prec-v+1)/p))**QQ(p)).reverse())
            # Check to see if expansion is exact.
            # If so, return an exact result.
            try:
                p = t_series.laurent_polynomial()
                if F(p(uvar)) == f:
                    return self.codomain()(p)
            except:
                pass
            return t_series

        return s_series

    def pushforward(self, f, prec=None):
        """
        Allows the map to work on differentials in addition to function field elements.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 - (x^2+1))          # optional - kash
            sage: D = (x/y).divisor()                         # optional - kash
            sage: D                                           # optional - kash
            - Place (x - I, y)
             + Place (x, y + x - 1)
             + Place (x, y + x + 1)
             - Place (x + I, y)
            sage: pl = D.support()[0]                         # optional - kash
            sage: m = F.completion(pl)                        # optional - kash
            sage: m(x.differential())                         # optional - kash
            [4*I*s + 8*I*s^3 + 12*I*s^5 + 16*I*s^7 + 20*I*s^9 + 24*I*s^11 + 28*I*s^13 + 32*I*s^15 + 36*I*s^17 + 40*I*s^19 + O(s^20)] ds
            sage: m(x/y)                                      # optional - kash
            1/2*s^-1 + 1/2*s + O(s^20)
            sage: m(x/y*x.differential())                     # optional - kash
            [2*I + 6*I*s^2 + 10*I*s^4 + 14*I*s^6 + 18*I*s^8 + 22*I*s^10 + 26*I*s^12 + 30*I*s^14 + 34*I*s^16 + 38*I*s^18 + O(s^20)] ds

            sage: K.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = K[]                                 # optional - kash
            sage: F.<y> = K.extension(y^2-x+1)                # optional - kash
            sage: pl = y.divisor().support()[1]               # optional - kash
            sage: m = F.completion(pl, prec=1)                # optional - kash
            sage: m(y*x.differential())                       # optional - kash
            [O(s^1)] ds

        """
        if prec is None:
            prec = self._precision

        t = self.domain().base_field().gen()  # all differentials construction w.r.t this differential
        vt = t.valuation(self._place)         # t's valuation
        vf = f._f.valuation(self._place)      # f's valuation
        st = self._expand(t, max(prec-vf,vt)+1) # t's series
        sdt = st.derivative()                 # dt's series
        sf = self._expand(f._f, prec-(vt-1))  # f's series
        return LaurentSeriesDifferential(sf.parent(), sf * sdt)

    def default_precision(self):
        """
        Return the default precision.

        EXAMPLES::

            sage: R.<x> = FunctionField(QQbar, implementation='kash') # optional - kash
            sage: L.<y> = R[]                                 # optional - kash
            sage: F.<y> = R.extension(y^2 + y + x + 1/x)      # optional - kash
            sage: D = (x*y).divisor()                         # optional - kash
            sage: p = D.support()[2]                          # optional - kash
            sage: m = F.completion(p)                         # optional - kash
            sage: m.default_precision()                       # optional - kash
            20
        """
        return self._precision


class FunctionFieldMaximalOrder_kash(FunctionFieldMaximalOrder):
    """
    Base class of kash-implemented maximal orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K # optional - kash
            Rational function field in t over Rational Field
            sage: R = K.maximal_order(); R                    # optional - kash
            Maximal order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrder.__init__(self, field, category)
        self.kash_constant_field = None

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.function_field().kash().MaximalOrderFinite()
        return self._kash_

    def _element_constructor_(self, f, check=True):
        """
        Construct an element of this order from ``f``.

        INPUT:

        - ``f`` -- element convertible to the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2-x*Y+x^2+1)          # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: y in O                                      # optional - kash
            True
            sage: 1/y in O                                    # optional - kash
            False
            sage: x in O                                      # optional - kash
            True
            sage: 1/x in O                                    # optional - kash
            False
            sage: L.<y>=K.extension(Y^2+Y+x+1/x)              # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: 1 in O                                      # optional - kash
            True
            sage: y in O                                      # optional - kash
            False
            sage: x*y in O                                    # optional - kash
            True
            sage: x^2*y in O                                  # optional - kash
            True
        """
        field = self.function_field()

        #if f.parent() is field:
        #    f = f.element()
        #f = self._field._ring(f)
        if check:
            if not kash._contains(self._field.to_kash(f).name(), self.kash().name()):
                raise TypeError("%r is not an element of %r"%(f,self))
        # return f and not field._element_class(self, f) because
        # 1. that's what FunctionFieldMaximalOrder_global does
        # 2. that's what makes "d in O" work (it tests if O(d) == d)
        # 3. but it's probably not right (O(d).parent() should report the order, not the field)
        # return field._element_class(self, f)
        return f

    def polynomial(self):
        """
        Return the defining polynomial of the function field of which this is an order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)    # optional - kash
            sage: O = L.equation_order()                      # optional - kash
            sage: O.polynomial()                              # optional - kash
            y^4 + x*y + 4*x + 1
        """
        return self._field.polynomial()

    def basis(self):
        """
        Return the basis of the order as a module over the polynomial ring.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: L.<y> = K.extension(y^4 + x*y + 4*x + 1)    # optional - kash
            sage: O = L.equation_order()                      # optional - kash
            sage: O.basis()                                   # optional - kash
            (1, y, y^2, y^3)
        """

        return tuple(self.kash().Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- list of generators or an ideal in a ring which
          coerces to this order

        EXAMPLES::

            sage: K.<y> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: O = K.maximal_order()                       # optional - kash
            sage: O.ideal(y)                                  # optional - kash
            Ideal (y) of Maximal order of Rational function field in y over Rational Field
            sage: O.ideal([y,1/y]) == O.ideal(y,1/y)          # multiple generators may be given as a list, optional - kash
            True

        A fractional ideal of a nontrivial extension::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: O = K.maximal_order()                       # optional - kash
            sage: I = O.ideal(x^2-4)                          # optional - kash
            sage: L.<y> = K.extension(y^2 - x^3 - 1)          # optional - kash
            sage: S = L.equation_order()                      # optional - kash
            sage: S.ideal(1/y)                                # optional - kash
            Ideal (1, (-1/(x^3 + 1))*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 = S.ideal(x^2-4); I2                     # optional - kash
            Ideal (x^2 - 4, (x^2 - 4)*y) of Order in Function field in y defined by y^2 - x^3 - 1
            sage: I2 == S.ideal(I)                            # optional - kash
            True

            sage: K.<x> = FunctionField(GF(7), implementation='kash'); R.<y> = K[] # optional - kash
            sage: O = K.maximal_order()                       # optional - kash
            sage: I = O.ideal(x^2-4)                          # optional - kash
            sage: L.<y> = K.extension(y^2 - x^3 - 1)          # optional - kash
            sage: S = L.equation_order()                      # optional - kash
            sage: S.ideal(1/y)                                # optional - kash
            Ideal (1, (6/(x^3 + 1))*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 = S.ideal(x^2-4); I2                     # optional - kash
            Ideal (x^2 + 3, (x^2 + 3)*y) of Order in Function field in y defined by y^2 + 6*x^3 + 6
            sage: I2 == S.ideal(I)                            # optional - kash
            True
        """

        return FunctionFieldIdeal_kash(self, gens)

class FunctionFieldMaximalOrderInfinite_kash(FunctionFieldMaximalOrderInfinite):
    """
    Base class of kash-implemented maximal infinite orders of function fields.
    """

    def __init__(self, field, category=None):
        """
        Initialize.

        INPUT:

        - ``field`` -- function field

        EXAMPLES::

            sage: K.<t> = FunctionField(QQ, implementation='kash'); K # optional - kash
            Rational function field in t over Rational Field
            sage: R = K.maximal_order_infinite(); R           # optional - kash
            Maximal infinite order of Rational function field in t over Rational Field
        """

        FunctionFieldMaximalOrderInfinite.__init__(self, field, category)
        self.kash_constant_field = None

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.function_field().kash().MaximalOrderInfinite()
        return self._kash_

    def _element_constructor_(self, f):
        """
        Make ``f`` an element of this order.

        EXAMPLES::

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.basis()                                # optional - kash
            (1, 1/x*y)
            sage: 1 in Oinf                                   # optional - kash
            True
            sage: 1/x*y in Oinf                               # optional - kash
            True
            sage: x*y in Oinf                                 # optional - kash
            False
            sage: 1/x in Oinf                                 # optional - kash
            True
        """
        if not f.parent() is self.function_field():
            f = self.function_field()(f)

        if not kash._contains(self._field.to_kash(f).name(), self.kash().name()):
            raise TypeError("%r is not an element of %r"%(f,self))

        return f

    def basis(self):
        """
        Return a basis of this order as a module over the maximal order
        of the base function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[] # optional - kash
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2) # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.basis()                                # optional - kash
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.basis()                                # optional - kash
            (1, 1/x*y)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<t> = K[] # optional - kash
            sage: L.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2) # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.basis()                                # not tested - kash returns a different (but equivalent) basis
            (1, 1/x^2*y, 1/x^4*y^2)

            sage: K.<x> = FunctionField(GF(2), implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.basis()                                # optional - kash
            (1, 1/x*y)
        """

        return tuple(self.kash().Basis().sage(self._field.reverse_map))

    def ideal(self, *gens):
        """
        Return the fractional ideal generated by ``gens``.

        INPUT:

        - ``gens`` -- tuple of elements of the function field

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[] # optional - kash
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2) # optional - kash
            sage: Oinf = F.maximal_order_infinite()           # optional - kash
            sage: I = Oinf.ideal(x,y); I                      # optional - kash
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: I = Oinf.ideal(x,y); I                      # optional - kash
            Ideal (x, y) of Maximal infinite order of Function field
            in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self, gens)

    def decomposition(self):
        """
        Return prime ideal decomposition of `pO_\infty` where `p` is the unique
        prime ideal of the maximal infinite order.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[] # optional - kash
            sage: F.<y> = K.extension(t^3 - x^2*(x^2 + x + 1)^2) # optional - kash
            sage: Oinf = F.maximal_order_infinite()           # optional - kash
            sage: Oinf.decomposition()                        # optional - kash
            [(Ideal (1/x, 1/x^2*y - 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 1, 1),
             (Ideal (1/x, 1/x^4*y^2 + 1/x^2*y + 1) of Maximal infinite order
             of Function field in y defined by y^3 - x^6 - 2*x^5 - 3*x^4 - 2*x^3 - x^2, 2, 1)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: Oinf = L.maximal_order_infinite()           # optional - kash
            sage: Oinf.decomposition()                        # optional - kash
            [(Ideal (1/x, 1/x*y) of Maximal infinite order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x, 1, 2)]
        """

        return [ ( FunctionFieldIdeal_kash(self, i), i.Degree().sage(), i.RamificationIndex().sage() )
                 for i in self.kash().Decomposition()]

class FunctionFieldIdeal_kash(FunctionFieldIdeal):
    """
    Fractional ideal of the maximal order of a kash-implemented function field.
    """

    def __init__(self, ring, gens):
        """
        Initialize.

        INPUT:

        - ``ring`` -- maximal order

        - ``gens``-- a list of generators, or a kash ideal
        """

        FunctionFieldIdeal.__init__(self, ring)

        if isinstance(gens, KashElement):
            self._kash_ = gens
            self._gens = tuple(gens.Generators().sage(ring._field.reverse_map))
        else:
            if ring.function_field().constant_base_field() is QQbar:
                ring.function_field().base_field()._extend_constant_field(gens)
            self._gens = tuple(flatten(gens))
            self._kash_ = ring.kash().Ideal(map(ring._field.to_kash, self._gens))

        self.kash_constant_field = ring.function_field().kash_constant_field

    def kash(self):
        if self.kash_constant_field != self.ring().function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.ring().function_field().base_field().kash_constant_field
            self._kash_ = self.ring().kash().Ideal(map(self.ring().function_field().to_kash, self._gens))
        return self._kash_

    def __hash__(self):
        """
        Return hash computed from the data.
        """

        return hash( (self._ring, self._gens) )

    def _richcmp_(self, other, op):
        """
        Compare the ideal with the other ideal with respect to ``op``.

        INPUT:

        - ``other`` -- ideal

        - ``op`` -- comparison operator

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(1/y)                            # optional - kash
            sage: I == I + I                                  # optional - kash
            True
            sage: I == I * I                                  # optional - kash
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(1/y)                            # optional - kash
            sage: I == I + I                                  # optional - kash
            True
            sage: I == I * I                                  # optional - kash
            False
            sage: I < I * I                                   # optional - kash
            True
            sage: I > I * I                                   # optional - kash
            False
        """

        return richcmp((self.denominator(), self.gens_over_base()), (other.denominator(), other.gens_over_base()), op)

    def __repr__(self):
        """
        Return a string representation of the ideal.
        """
        return "Ideal (%s) of %s"%(', '.join([repr(g) for g in self.gens()]), self.ring())

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)          # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal([y]); I                         # optional - kash
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 - x^3 - 1
            sage: x * y in I                                  # optional - kash
            True
            sage: y / x in I                                  # optional - kash
            False
            sage: y^2 - 2 in I                                # optional - kash
            False

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal([y]); I                         # optional - kash
            Ideal (y) of Maximal order of Function field in y
            defined by y^2 + y + (x^2 + 1)/x
            sage: x * y in I                                  # optional - kash
            True
            sage: y / x in I                                  # optional - kash
            False
            sage: y^2 - 2 in I                                # optional - kash
            False
        """
        kashx = self._ring._field.to_kash(x)
        return kash.eval("%s in %s" % (kashx.name(), self.kash().name())) == 'TRUE'

    def __invert__(self):
        """
        Return the inverse fractional ideal of the ideal.

        EXAMPLES::
            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3 - 1)          # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: ~I                                          # optional - kash
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: I^(-1)                                      # optional - kash
            Ideal (1, (1/(x^3 + 1))*y) of Maximal order of Function field in y defined by y^2 - x^3 - 1
            sage: ~I * I                                      # optional - kash
            Ideal (1) of Maximal order of Function field in y defined by y^2 - x^3 - 1

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y>=K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: ~I                                          # optional - kash
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: I^(-1)                                      # optional - kash
            Ideal (x, (x/(x^2 + 1))*y + x/(x^2 + 1)) of Maximal order
            of Function field in y defined by y^2 + y + (x^2 + 1)/x
            sage: ~I * I                                      # optional - kash
            Ideal (1) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, 1 / self.kash())

    def _add_(self, other):
        """
        Add with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: I + J                                       # optional - kash
            Ideal (x, y) of Maximal order of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: I + J                                       # optional - kash
            Ideal (1, y) of Maximal order of Function field in y defined by y^2 + y + (x^2 + 1)/x
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash() + other.kash())

    def _mul_(self, other):
        """
        Multiply with other ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: I * J                                       # optional - kash
            Ideal (x^4 + x^2 - x, x*y + x^2) of Maximal order
            of Function field in y defined by y^2 - x^3*y - x

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: I * J                                       # optional - kash
            Ideal ((x^5 + x^3 + x^2 + 1)/x,
                   y + (1/2*x^4 + 1/2*x^3 + x^2 + 1/2*x + 1/2)/x)
                of Maximal order of Function field in y
                defined by y^2 + y + (x^2 + 1)/x

        TESTS:

        Verify the examples using Sage's standard ideal operations::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: IJ = I * J                                  # optional - kash

            sage: R = PolynomialRing(QQ, K.gens() + L.gens()) # optional - kash
            sage: Q = R.quo(L.polynomial())                   # optional - kash
            sage: I1 = ideal(Q(g*b) for g in I.gens() for b in O.basis()) # optional - kash
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis()) # optional - kash
            sage: IJ = ideal(Q(g*b) for g in (IJ).gens() for b in O.basis()) # optional - kash
            sage: I1 * J1 == IJ                               # optional - kash
            True

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: J = O.ideal(x+y)                            # optional - kash
            sage: IJ = I * J                                  # optional - kash

            sage: Q = R.quo(L.polynomial().numerator())       # optional - kash
            sage: I1 = ideal(Q(g*b*x) for g in I.gens() for b in O.basis()) # optional - kash
            sage: J1 = ideal(Q(g*b) for g in J.gens() for b in O.basis()) # optional - kash
            sage: IJ = ideal(Q(g*b*x) for g in (IJ).gens() for b in O.basis()) # optional - kash
            sage: I1 * J1 == IJ                               # optional - kash
            True
        """

        return FunctionFieldIdeal_kash(self._ring, self.kash() * other.kash())

    def denominator(self):
        """
        Return the denominator of the fractional ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); R.<y> = K[] # optional - kash
            sage: L.<y> = K.extension(y^2 - x^3 - 1)          # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y/(y+1))                        # optional - kash
            sage: d = I.denominator(); d                      # optional - kash
            x^3
            sage: d in O                                      # optional - kash
            True
        """
        return self.kash().Denominator().sage(self._ring._field.reverse_map)

    def factor(self):
        """
        Return the factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K) # optional - kash
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)    # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: I == I.factor().prod()                      # optional - kash
            True

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: I == I.factor().prod()                      # optional - kash
            True
        """
        factors = self._factor()
        return Factorization(factors, cr=True)

    def _factor(self):
        """
        Return the list of prime and multiplicity pairs of the
        factorization of the ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = K[] # optional - kash
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)    # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: I == I.factor().prod()                      # indirect doctest, optional - kash
            True
        """
        factors = []
        for f,m in self.kash().Factorization():
            factors.append( (FunctionFieldIdeal_kash(self._ring, f), m) )
        return factors

    def gens(self):
        """
        Return the generators of the ideal.

        This provides whatever set of generators as quickly
        as possible.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(x+y)                            # optional - kash
            sage: I.gens()                                    # optional - kash
            (y + x,)

            sage: L.<y> = K.extension(Y^2 +Y + x + 1/x)       # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(x+y)                            # optional - kash
            sage: I.gens()                                    # optional - kash
            (y + x,)
        """
        return self._gens

    def gens_over_base(self):
        """
        Return the generators of the ideal as a module over the
        maximal order of the base rational function field.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 - x^3*Y - x)        # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(x+y)                            # optional - kash
            sage: I.gens_over_base()                          # optional - kash
            (x^4 + x^2 - x, y + x)

            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(x+y)                            # optional - kash
            sage: I.gens_over_base()                          # optional - kash
            (x^3 + 1, y + x)
        """
        return tuple(self.kash().Basis().sage(self._ring._field.reverse_map))

    def is_prime(self):
        """
        Return ``True`` if the ideal is a prime ideal.

        If checked to be a prime ideal, then the ideal can be used
        as a prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K) # optional - kash
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)    # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: I.is_prime()                                # optional - kash
            False
            sage: [f.is_prime() for f,_ in I.factor()]        # optional - kash
            [True, True]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: I.is_prime()                                # optional - kash
            False
            sage: [f.is_prime() for f,_ in I.factor()]        # optional - kash
            [True, True]
        """

        return bool(self.kash().IsPrime())

    def place(self):
        """
        Return the place corresponding to the prime ideal.

        EXAMPLES::

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<t> = PolynomialRing(K) # optional - kash
            sage: F.<y> = K.extension(t^3-x^2*(x^2+x+1)^2)    # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: [f.place() for f,_ in I.factor()]           # optional - kash
            [Place (x, (1/(x^3 + x^2 + x))*y^2),
             Place (x^2 + x + 1, (1/(x^3 + x^2 + x))*y^2)]

            sage: K.<x> = FunctionField(QQ, implementation='kash'); _.<Y> = K[] # optional - kash
            sage: L.<y> = K.extension(Y^2 + Y + x + 1/x)      # optional - kash
            sage: O = L.maximal_order()                       # optional - kash
            sage: I = O.ideal(y)                              # optional - kash
            sage: [f.place() for f,_ in I.factor()]           # optional - kash
            [Place (x, x*y), Place (x^2 + 1, x*y)]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")
        return FunctionFieldPlace_kash(self._ring._field, self)

    def valuation(self, ideal):
        """
        Return the valuation of the ideal at this prime ideal.

        INPUT:

        - ``ideal`` -- fractional ideal

        EXAMPLES::

            sage: F.<x> = FunctionField(QQ, implementation='kash') # optional - kash
            sage: O = F.maximal_order()                       # optional - kash
            sage: I = O.ideal(x^2*(x^2+x+1)^3)                # optional - kash
            sage: [f.valuation(I) for f,_ in I.factor()]      # optional - kash
            [2, 3]
        """
        if not self.is_prime():
            raise TypeError("not a prime ideal")

        return ideal.kash().Valuation(self.kash()).sage()


class FunctionFieldPlace_kash(FunctionFieldPlace):
    """
    Places of kash-implemented function field.
    """

    def __init__(self, field, arg):
        """
        Initialize the place.

        INPUT:

        - ``field`` -- function field

        - ``arg`` -- prime ideal associated with the place, or a kash place

        .. TESTS::

            Check that when creating Sage places from kash places, we set
            ``order`` correctly.  Otherwise, the divisor below isn't
            calculated correctly.

            sage: R.<x> = FunctionField(QQbar, implementation='kash'); # optional - kash
            sage: L.<Y> = R[];                                # optional - kash
            sage: F.<y> = R.extension(Y^4 - (x^2+1)^3);       # optional - kash
            sage: (1/y * x.differential()).divisor()          # indirect doctest, optional - kash
            0

        """

        if isinstance(arg, KashElement):
            ideal = arg.Ideal()
            if arg.IsFinite():
                order = field.maximal_order()
            else:
                order = field.maximal_order_infinite()
            prime = FunctionFieldIdeal_kash(order, order.kash().CoerceIdeal(ideal))
            self._kash_ = arg
        else:
            prime = arg
            self._kash_ = arg.kash().Place()

        FunctionFieldPlace.__init__(self, field, prime)
        self.kash_constant_field = field.base_field().kash_constant_field

    def kash(self):
        if self.kash_constant_field != self.function_field().base_field().kash_constant_field:
            self.kash_constant_field = self.function_field().base_field().kash_constant_field
            self._kash_ = self.prime_ideal().kash().Place()
        return self._kash_

    def is_infinite_place(self):
        """
        Return ``True`` if the place is at infinity.
        """
        F = self.function_field()
        return self.prime_ideal().ring() == F.maximal_order_infinite()

    def degree(self):
        return self.kash().Degree()
