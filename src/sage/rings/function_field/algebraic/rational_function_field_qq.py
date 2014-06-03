"""
Algebraic Function Fields

EXAMPLES:

We create an extension of a rational function fields, and do some
simple arithmetic in it::

    sage: K.<x> = RationalFunctionField(GF(5^2,'a')); K
    Rational function field in x over Finite Field in a of size 5^2
    sage: R.<y> = K[]
    sage: L.<y> = K.extension(y^3 - (x^3 + 2*x*y + 1/x)); L
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: y^2
    y^2
    sage: y^3
    2*x*y + (x^4 + 1)/x
    sage: a = 1/y; a
    (4*x/(4*x^4 + 4))*y^2 + 2*x^2/(4*x^4 + 4)
    sage: a * y
    1

We next make an extension of the above function field, illustrating
that arithmetic with a tower of 3 fields is fully supported::

    sage: S.<t> = L[]
    sage: M.<t> = L.extension(t^2 - x*y)
    sage: M
    Function field in t defined by t^2 + 4*x*y
    sage: t^2
    x*y
    sage: 1/t
    ((1/(x^4 + 1))*y^2 + 2*x/(4*x^4 + 4))*t
    sage: M.base_field()
    Function field in y defined by y^3 + 3*x*y + (4*x^4 + 4)/x
    sage: M.base_field().base_field()
    Rational function field in x over Finite Field in a of size 5^2

TESTS::

    sage: TestSuite(K).run()
    sage: TestSuite(L).run()  # long time (8s on sage.math, 2012)
    sage: TestSuite(M).run()  # long time (52s on sage.math, 2012)

The following two test suites do not pass ``_test_elements`` yet since
``R.an_element()`` has a ``_test_category`` method wich it should not have.
It is not the fault of the function field code so this will
be fixed in another ticket::

    sage: TestSuite(R).run(skip = '_test_elements')
    sage: TestSuite(S).run(skip = '_test_elements')
"""

import collections

from sage.categories.function_fields import FunctionFields

from sage.structure.parent import Parent

from sage.rings.rational_field import RationalField
from sage.rings.function_field.algebraic.rational_function_qq import RationalFunctionQQ

from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class RationalFunctionFieldQQ(Parent):
    Element = RationalFunctionQQ
    def __init__(self, names, category = None):
        if (isinstance(names, collections.Iterable)
                and not isinstance(names, basestring)):
            names = tuple(*names)
        else:
            names = (names, )

        if len(names) > 1:
            raise ValueError("only supports univariant rational functions")

        constant_field = RationalField()
        self._ring = PolynomialRing(constant_field, names)

        if category is None:
            #category = FunctionFields(constant_field)
            category = FunctionFields()

        Parent.__init__(self, base=constant_field, names=names, category=category)

    def __hash__(self):
        return hash((self.base_ring(), self.variable_names()))

    def _repr_(self):
        return "Rational function field in %s over %s"%(
            ', '.join(self.variable_names()), self.base_ring())

    def _element_constructor_(self, *args):
        return self.element_class(self, *args)

    def gens(self):
        return (self.element_class(self, [0, 1]),)

    def genus(self):
        return 0

    def _coerce_map_from_(self, other):
        if self.base_ring().has_coerce_map_from(other):
            return True
        if isinstance(other, PolynomialRing_general):
            return self.base_ring().has_coerce_map_from(other.base_ring())
