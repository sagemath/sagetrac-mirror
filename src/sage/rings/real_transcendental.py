# -*- coding: utf-8 -*-
r"""
Real transcendental field extensions

A real transcendental extension of a base field is an extension of a field by
adjoining one or more reals. The reals must be collectively transcendental over
the base field, i.e., no polynomial equations are satisfied. Then the field
extension is isomorphic to the field of rational functions over the generators.
We can thus do exact arithmetic.

This file also contains classes for working with random lazy real numbers. Note
that if r is taken at random from an interval with respect to Lebesgue measure,
then Q(r) is almost surely a transcendental extension.

AUTHORS:

- Pat Hooper (2018): initial commit (:trac:`26042`)
"""

#*****************************************************************************
#       Copyright (C) 2018 Pat Hooper <wphooper@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import absolute_import

from sage.categories.fields import Fields
from sage.categories.sets_cat import Sets
from sage.misc.cachefunc import cached_method
from sage.structure.element import Element, FieldElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE
from sage.structure.unique_representation import UniqueRepresentation
from .function_field.constructor import FunctionField
from .function_field.function_field_element import FunctionFieldElement
from .integer_ring import ZZ
from .real_double import RDF
from .real_lazy import RLF
from .real_mpfi import RealIntervalField
from .real_mpfr import RealField, RR
from .ring import Field

class LazyRandomNumber(Element):
    r"""
    Represents a random number in [0,1] chosen with respect to Lebesgue measure.

    The random number is computed lazily. You can access a representation in a
    real interval field. The number is improved randomly when more precision
    is requested.
    """
    def __init__(self, parent):
        assert isinstance(parent, LazyRandomNumberSet_class)
        Element.__init__(self, parent)
        self._approx=[]

    def __hash__(self):
        return 43*hash(self.__float__())

    def _real_mpfi_(self, real_interval_field):
        i = (real_interval_field.prec()-1)//53
        if i < len(self._approx):
            return real_interval_field(self._approx[i])
        else:
            from sage.rings.integer_ring import ZZ
            from sage.rings.real_mpfi import RealIntervalField
            if len(self._approx)==0:
                k = ZZ.random_element(0, 2**53)
                self._approx.append(RealIntervalField(53)(k/2**53, (k+1)/2**53))
            for j in xrange(len(self._approx), i+1):
                k = ZZ.random_element(0, 2**53)
                a = RealField( 53*(j+1) )(self._approx[j-1].lower())
                self._approx.append(RealIntervalField(53*(j+1))( a + k / ZZ(2)**(53*(j+1)),
                                                                 a + (k+1) / ZZ(2)**(53*(j+1)) ) )
        return real_interval_field(self._approx[i])

    def eval(self, R):
        r"""
        Convert ``self`` into an element of the field ``R``.
        """
        if R is float:
            return self.__float__()
        from sage.rings.real_mpfi import RealIntervalField
        return R(self._real_mpfi_(RealIntervalField(R.prec())))

    def _mpfr_(self, R):
        return self.eval(R)

    def _real_double_(self, R):
        return self.eval(R)

    def __float__(self):
        return float(self.eval(RR))

    def __repr__(self):
        return repr(self._real_mpfi_(RIF))

class LazyRandomNumberSet_class(Parent, UniqueRepresentation):
    r"""
    The set of random real numbers of type LazyRandomNumber.
    """
    Element = LazyRandomNumber

    def __init__(self):
        Parent.__init__(self, category=Sets())

    def __repr__(self):
        return "Set of lazy random numbers"

    def _an_element_(self):
        return self.element_class(self)

    def random_real(self):
        r"""
        Return a random real number in the interval [0,1] as an element of RealLazyField.

        The number is randomly extended as needed.
        """
        return RLF(self.element_class(self))

LazyRandomNumberSet = LazyRandomNumberSet_class()

def random_real(a=0, b=1):
    r"""
    Return a real number in the interval [a, b] as an element of RealLazyField.
    The real number is determined at random with respect to Lebesgue measure on
    the interval. By default the interval is taken to be [0, 1].

    The number is randomly extended as needed. This has consequences for
    pickling. A restored pickled random number will agree with the original
    random number up to the known accuracy at the time of pickling, but if
    further digits are needed later, then the numbers will diverge in value.
    Consequently a random number is unequal to its recovered pickle.

    EXAMPLES::

    By default we return a number in [0,1]::

        sage: from sage.rings.real_transcendental import random_real
        sage: 0 < random_real() < 1
        True

    We can get any interval::

        sage: from sage.symbolic.constants import e, pi
        sage: RLF(e) < random_real(e, pi) < RLF(pi)
        True

    Two random reals are never the same (unless they are the same object)::

        sage: random_real() == random_real()
        False

    Pickling::

        sage: r = random_real()
        sage: 0 < float(r) < 1
        True
        sage: rr = loads(dumps(r))
        sage: float(r) == float(rr)
        True
        sage: F = RealField(100)
        sage: F(r) == F(rr)
        False
        sage: r == rr
        False
        sage: TestSuite(r).run(skip="_test_pickling")
    """
    if a==0 and b==1:
        return LazyRandomNumberSet.random_real()
    else:
        assert a != b, "Can not construct random number in degenerate interval."
        return RLF(a)+RLF(b-a)*LazyRandomNumberSet.random_real()

class RealTranscendentalExtensionFieldElement(FieldElement):
    r"""
    An element of a RealTranscendentalExtensionField.
    """
    def __init__(self, parent, rf):
        r"""
        Construct an element in the field from a rational function.
        """
        # Presumably rf has already been coerced into the correct function field.
        assert rf.parent() == parent.function_field()
        self._rf = rf

        #try:
        #    # This .factor().expand() stuff is to deal with the fact that
        #    # FunctionField does not automatically cancel common factors:
        #    # sage: K.<x>=FunctionField(QQ)
        #    # sage: (2*x)/(2*(x+1))
        #    # 2*x/(2*x + 2)
        #    self._rf = parent.function_field()(rf.factor().expand())
        #except ArithmeticError:
        #    self._rf = rf
        FieldElement.__init__(self, parent)

    def __hash__(self):
        return hash(self._rf)

    def rational_function(self):
        r"""
        Return the underlying rational function (defined in terms of the generators
        of the parent).
        """
        return self._rf

    @cached_method
    def real_lazy(self):
        r"""
        Return this number as an element of RealLazyField by viewing it as a
        rational function in the generators of the transcendental extension.
        """
        return self.parent()._eval(self._rf)

    def _real_mpfi_(self, real_field):
        return real_field(self.real_lazy())

    def __repr__(self):
        return repr(self._rf)

    def _add_(self, other):
        return self.parent().element_class(self.parent(), self._rf+other._rf)

    def _sub_(self, other):
        return self.parent().element_class(self.parent(), self._rf-other._rf)

    def _neg_(self):
        return self.parent().element_class(self.parent(), -self._rf)

    def _mul_(self, other):
        other =  self.parent()(other)
        return self.parent().element_class(self.parent(), self._rf*other._rf)

    def _invert_(self):
        return self.parent().element_class(self.parent(), ~self._rf)

    def _div_(self, other):
        return self.parent().element_class(self.parent(), self._rf/other._rf)

    def _pow_int(self, k):
        return self.parent().element_class(self.parent(), self._rf**k)

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a numerical approximation of "self" with "prec" bits (or
        decimal "digits") of precision.
        """
        if self == self.parent().zero():
            return self.parent().constant_field().zero().numerical_approx(prec=prec, digits=digits, algorithm=algorithm)
        if digits is None:
            if prec is None:
                return self._mpfr_(RDF)
            return self._mpfr_(RealField(prec))
        else:
            from sage.functions.other import ceil
            return self._mpfr_(RealField( (digits * RDF(10).log2()).ceil() ))

    def _mpfr_(self, R):
        """
        EXAMPLES::

            sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
            sage: RTE.<ee> = RealTranscendentalExtensionField(QQ, e)
            sage: RF = RealField(200)
            sage: RF(e) == RF(ee)
            True

            sage: from sage.rings.real_transcendental import random_real, RealTranscendentalExtensionField
            sage: RTE.<r> = RealTranscendentalExtensionField(QQ, random_real())
            sage: z = 2^100*r - floor(2^100*r)
            sage: # Note the below will fail if r is in RealLazyField
            sage: 0 < z < 1
            True
            sage: 0 < float(z) < 1
            True
        """
        prec = R.precision()
        rl = self.real_lazy()
        guess = RealIntervalField(prec)(rl)
        while guess.relative_diameter() > 2**(-R.precision()):
            prec += 53
            guess = RealIntervalField(prec)(rl)
        return R(guess)

    def _real_double_(self, R):
        """
        EXAMPLES::

            sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
            sage: RTE.<ee> = RealTranscendentalExtensionField(QQ, e)
            sage: RF = RealField(200)
            sage: RDF(e) == RDF(ee)
            True
        """
        return R(self._mpfr_(RDF))

    def __float__(self):
        """
        EXAMPLES::

            sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
            sage: RTE.<ee> = RealTranscendentalExtensionField(QQ, e)
            sage: RF = RealField(200)
            sage: float(e) == float(ee)
            True
        """
        return float(self._mpfr_(RDF))

    def _richcmp_(self, other, op):
        if self._rf == other._rf:
            # Equal!
            return op == op_LE or op == op_EQ or op == op_GE
        if op == op_EQ:
            return False
        if op == op_NE:
            return True
        diff = (self-other).real_lazy()
        i=1
        F = RealIntervalField(53)
        value = F(diff)
        while 0 in value:
            i += 1
            F = RealIntervalField(53*i)
            value = F(diff)
        return richcmp(value,F.zero(), op)

    def euclidean_degree(self):
        r"""
        Return the degree of this element as an element of an Euclidean domain.

        In a field, this returns 0 for all but the zero element (for which
        it is undefined).
        """
        if self == self.parent().zero():
            raise ValueError("euclidean degree not defined for the zero element")
        return 0

class RealTranscendentalExtensionField(Field):
    r"""
    Represents a real transcendental extension of a field in sage.

    The field we are extending is called the "constant field".

    EXAMPLES::

    The field Q(e)::

        sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
        sage: RTE.<ee> = RealTranscendentalExtensionField(QQ, e)
        sage: RF = RealField(200)
        sage: RF(e) == RF(ee)
        True
        sage: TestSuite(RTE).run()

    The field Q(pi)::

        sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
        sage: RTE.<pipi> = RealTranscendentalExtensionField(QQ, pi)
        sage: RF = RealField(200)
        sage: RF(pipi) == RF(pi)
        True
        sage: TestSuite(RTE).run()

    The field Q(r) where r in [0,1] is taken at random:

        sage: from sage.rings.real_transcendental import random_real
        sage: K.<r> = RealTranscendentalExtensionField(QQ, random_real())
        sage: 0 < r < 1
        True
        sage: TestSuite(r).run(skip="_test_pickling")
        sage: TestSuite(K).run(skip=["_test_elements","_test_pickling"])

    Exact arithmetic with a random unit vector::

        sage: # Construct a random angle
        sage: theta = random_real(0, 2*pi)
        sage: # Choose y so that the point (1,y) lies on the line
        sage: # through (-1,0) and (cos(theta), sin(theta))
        sage: y_RLF = 2*sin(theta)/(1+cos(theta))
        sage: from sage.rings.real_transcendental import RealTranscendentalExtensionField
        sage: # Construct Q(y)
        sage: F.<y> = RealTranscendentalExtensionField(QQ, y_RLF)
        sage: # The random unit vector (cos(theta), sin(theta))
        sage: unit_vector = vector(F, [(4-y^2)/(4+y^2), 4*y/(4+y^2)])
        sage: a,b = unit_vector
        sage: a^2 + b^2
        1
        sage: -1 < a < 1 and -1 < b < 1
        True
        sage: # Pickling is okay here because of the way RLF tests for equality:
        sage: TestSuite(F).run()
        sage: TestSuite(a).run()
        sage: TestSuite(b).run()

    Extending by more than one element:

        sage: FF.<r1,r2> = RealTranscendentalExtensionField(QQ, [random_real(), random_real()])
        sage: RDF(r1-r2) != 0
        True
        sage: 0 < r1 < 1 and 0 < r2 < 1
        True
        sage: TestSuite(FF).run(skip=["_test_elements","_test_pickling"])

        sage: F.<sqrt2> = NumberField(x^2-2, embedding=AA(sqrt(2)))
        sage: L.<ee,r> = RealTranscendentalExtensionField(F, [e, random_real()])
        sage: r < ee-sqrt2
        True
        sage: TestSuite(FF).run(skip=["_test_elements","_test_pickling"])
        sage: TestSuite( r*(ee-sqrt2) ).run(skip="_test_pickling")
    """

    Element = RealTranscendentalExtensionFieldElement

    def __init__(self, constant_field, transcendentals, names=None):
        r"""
        Construct the real transcendental extension of the provided
        constant_field by adjoining a transcendental or list of transcendentals.

        INPUT:

        - ``constant_field`` -- a field to extend

        - ``transcendentals`` -- a NestedIntervalReal defined over constant_field which
          represents a transcentental number over the constant_field, or a list of such
          transcentals. If a list is provided, then each number in the list be independently
          transcendental (i.e., no polynomial equation in the transcendentals with coefficients
          in the constant_field should be satisfied).

        - ``names`` -- a string giving a name to the transcendental generator or, in the case
          a list of transcendentals was provided, a list of names for the transcentals. If names
          are not provided, then the name `t` is given to a single transcendental and the names
          `t_0`, `t_1`, ... are given to a list of transcendentals.
        """
        self._constant_field = constant_field
        Field.__init__(self, self._constant_field, category=Fields())
        try:
            # This line is just to test if transcendentals is a single number or a list:
            real_value = RDF(transcendentals)
            # Found just one transcendental
            transcendentals = [transcendentals]
        except TypeError:
            pass
        n = len(transcendentals)
        if n == 0:
            raise ValueError("Must have at least one transcendental to extend.")
        self._values = tuple([RLF(t) for t in transcendentals])
        if not names is None:
            if n == 1:
                if isinstance(names,str):
                    self._names = (names,)
                else:
                    assert len(names)==n, "Wrong number of names provided."
                    assert isinstance(names[0],str)
                    self._names = (names[0],)
            else:
                self._names = tuple([str(names[i]) for i in xrange(n)])
        else:
            if n == 1:
                self._names = ("t",)
            else:
                self._names = tuple(["t_"+str(i) for i in xrange(n)])

        # Construct the underlying function field.
        ff = self._constant_field
        gens = []
        for name in self._names:
            ff = FunctionField(ff,name)
            # Convert the generators that already exist into the current Field.
            gens = [ff(gen) for gen in gens]
            gens.append(ff.gen())

        self._ff = ff

        self._gens = tuple([self.element_class(self,gens[i]) for i in xrange(n)])

    def constant_field(self):
        r"""
        Return the field which is being extended.
        """
        return self._constant_field

    def function_field(self):
        r"""
        Return the underlying field of rational functions over the generators.
        """
        return self._ff

    def gen(self, n=None):
        r"""
        Return the n-th generator. If there is only one generator, then n can be omitted.
        """
        if n is None:
            assert len(self._gens)==1, "Which generator do you want? Use gen(n)."
            return self._gens[0]
        return self._gens[n]

    def gens(self):
        r"""
        Return the list of generators for the extension.
        """
        return self._gens

    def some_elements(self):
        r"""
        Return an iterator over some of the elements of this field.
        """
        for f in self._ff.some_elements():
            yield self.element_class(self, f)

    def transcendental(self, n=None):
        r"""
        Return the n-th transcendental used to construct the extension.

        If there is only one generator then n is not needed.
        """
        if n is None:
            assert len(self._values)==1, "Which transcendental do you want? Use transcendental(n)."
            return self._values[0]
        return self._values[n]

    @cached_method
    def _transcendental_power(self, n, power):
        r"""
        Return an integer power of the $n$-th generator.
        """
        if power==0:
            return RLF.one()
        return self.transcendental(n)*self._transcendental_power(n,power-1)

    def _eval(self, rf, num_variables=None):
        r"""
        Return a NestedIntervalReal determined by evaluating the provided rational function
        (element of the function_field) over the generators considered to be
        NestedIntervalReals.
        """
        if num_variables is None:
            num_variables = len(self._gens)
        if num_variables == 0:
            return rf
        n = num_variables - 1 # The variable we are working with.
        polys = [rf.numerator(), rf.denominator()]
        polys_out = []
        for poly in polys:
            polys_out.append(sum([
                self._transcendental_power(n,power)*self._eval(expression,n)
                for power,expression in poly.dict().items()
            ]))
        return polys_out[0]/polys_out[1]

    def transcendentals(self):
        r"""
        Return the list of transcendental generators.
        """
        return self._values

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, RealTranscendentalExtensionField):
            return False
        st = self.transcendentals()
        ot = other.transcendentals()
        if len(st) != len(ot):
            return False
        for i in xrange(len(st)):
            if st[i] != ot[i]:
                return False
        return True

    def __repr__(self):
        if len(self._names)==1:
            return "Real transcendental extension of %s by %s" % \
                (self._constant_field, self.gen())
        else:
            return "Real transcendental extension of %s by %s" % \
                (self._constant_field, self.gens())

    def _element_constructor_(self, x):
        if isinstance(x,RealTranscendentalExtensionFieldElement):
            if x.parent() is self:
                return x
            assert x.parent() == self
            return self.element_class(self, self._ff(x._rf))
        if isinstance(x, FunctionFieldElement):
            return self.element_class(self, self._ff(x))
        xx = self._constant_field(x) # The quantity x should be in the constant field.
        return self.element_class(self, self._ff(xx))

    def _coerce_map_from_(self, S):
        r"""
        The only thing that we accept coersions from is the constant field.
        """
        if S == self:
            return True
        if S is self._constant_field:
            return True
        if self._constant_field.has_coerce_map_from(S):
            return True
        return False

    def is_exact(self):
        r"""
        Return True; We do exact arithmetic.
        """
        return True

    def is_perfect(self):
        r"""
        Return whether the field is perfect, i.e., its characteristic `p` is zero
        or every element has a `p`-th root.

        Returns True since a real subfield has characteristic zero.
        """
        return True

    def characteristic(self):
        """
        Return the characteristic of the this real subfield which is zero.
        """
        return ZZ.zero()

