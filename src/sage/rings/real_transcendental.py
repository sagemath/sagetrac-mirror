# -*- coding: utf-8 -*-
r"""
Real transcendental field extensions

A real transcendental extension of a base field is an extension of a field by
adjoining one or more reals. The reals must be collectively transcendental over
the base field, i.e., no polynomial equations are satisfied. Then the field
extension is algebraically conjugate to the field of rational functions over
the generators. We can thus do exact arithmetic.

We also provide classes to work with transcentals as real numbers. We consider
real numbers that can be described as the unique point in an intersection of
closed intervals with endpoints in the base field. The parent of such numbers
is the NestedIntervalRealSet over the base field. Our transcendentals must be
given as elements of this set.

EXAMPLES::

Construct $\mathbb Q(e)$:

    sage: from sage.rings.real_transcendental import *
    sage: NIR = NestedIntervalRealSet(QQ)
    sage: K = RealTranscendentalExtensionField(QQ, NIR.e(), "ee")
    sage: ee = K.gen()
    sage: from sage.functions.other import factorial
    sage: ee > sum([1/factorial(i) for i in xrange(25)])
    True
    sage: ee < sum([1/factorial(i) for i in xrange(25)]) + 1/factorial(24)
    True

Construct $\mathbb Q(sqrt2)(e)$:

    sage: F.<sqrt2> = NumberField(x^2-2, embedding=AA(sqrt(2)))
    sage: NIR_F = NestedIntervalRealSet(F)
    sage: L = RealTranscendentalExtensionField(F, NIR_F.e(), "eee")
    sage: eee = L.gen()
    sage: (eee/sqrt2 + 7).parent()
    Real transcendental extension of Number Field in sqrt2 with defining polynomial x^2 - 2 by eee
    sage: (eee-sqrt2).n(digits=5)
    1.3041
    sage: eee < 2*sqrt2
    True

Construct $\mathbb Q(r)$ where $r$ is a random number between zero and one.

    sage: M = RealTranscendentalExtensionField(QQ, NIR.random(0,1), "r")
    sage: r = M.gen()
    sage: 0 < r < 1
    True

Construct $\mathbb Q(\sqrt{2})(e, r_1, r_2)$ where $r_1$ is a random number 
between zero and one and $r_2$ is a random number in $[\sqrt{2},2]$.

    sage: N = RealTranscendentalExtensionField(F, \
        (NIR_F.e(), NIR_F.random(0,1), NIR_F.random(sqrt2,2)), \
        ("eeee", "r_1", "r_2"))
    sage: eeee, r1,r2 = N.gens()
    sage: 0 < r1 < 1
    True
    sage: sqrt2 < r2 < 2
    True
    sage: r2 < eeee
    True

AUTHORS:

- Pat Hooper (2018): initial commit (:trac:`26042`)
"""
# python3
from __future__ import absolute_import

from exceptions import ArithmeticError
from sage.categories.fields import Fields
from sage.categories.sets_cat import Sets
from sage.functions.log import log
from sage.functions.other import factorial, floor
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.misc.prandom import randint
from sage.structure.element import Element, FieldElement
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp, op_LT, op_LE, op_EQ, op_NE, op_GT, op_GE
from sage.structure.unique_representation import UniqueRepresentation
from .function_field.constructor import FunctionField
from .integer_ring import ZZ
from .ring import Field


class ClosedInterval(Element):
    r"""
    Represents a bounded closed interval in the real numbers with endpoints
    in a field (obtainable via .parent().base_ring() )

    The class supports interval arithmetic.

    EXAMPLES::

    sage: from sage.rings.real_transcendental import ClosedIntervalSet
    sage: CIS = ClosedIntervalSet(QQ)
    sage: I = CIS(1,2); I
    [1, 2]
    sage: 1 in I
    True
    sage: 2 in I
    True
    sage: 3/2 in I
    True
    sage: 4 in I
    False
    """

    def __init__(self, parent, a, b):
        r"""
        Initialize the closed and bounded interval [a,b].
        """
        self._a = a
        self._b = b
        Element.__init__(self, parent)

    def lower(self):
        r"""
        Return the lower endpoint of the interval.
        """
        return self._a

    def upper(self):
        r"""
        Return the upper endpoint of the interval.
        """
        return self._b

    def center(self):
        r"""
        Return the midpoint of the interval.
        """
        return (self._a+self._b)/2

    def diameter(self):
        r"""
        Return the distance between the endpoints of the interval.
        """
        return self._b-self._a

    def intersects(self, other):
        r"""
        Return true if this interval intersects other, and False otherwise.

        EXAMPLES::

        sage: from sage.rings.real_transcendental import ClosedIntervalSet
        sage: CIS = ClosedIntervalSet(QQ)
        sage: CIS(1,2).intersects(CIS(2,3))
        True
        sage: CIS(1,2).intersects(CIS(5/2,3))
        False
        """
        return self.parent().base_ring().zero() in (self-other)

    def is_degenerate(self):
        r"""
        Return true if this interval has the form [a,a] for some a.

        EXAMPLES::

        sage: from sage.rings.real_transcendental import ClosedIntervalSet
        sage: CIS = ClosedIntervalSet(QQ)
        sage: CIS(0,1).is_degenerate()
        False
        sage: CIS(0,0).is_degenerate()
        True
        """
        return self._a == self._b

    def _repr_(self):
        return "["+str(self._a)+", "+str(self._b)+"]"

    def __add__(self, other):
        assert self.parent()==other.parent(), "Must have same parents"
        return ClosedInterval(self.parent(), self._a+other._a, self._b+other._b)

    def __sub__(self, other):
        assert self.parent()==other.parent(), "Must have same parents"
        return ClosedInterval(self.parent(), self._a-other._b, self._b-other._a)

    def __neg__(self):
        return ClosedInterval(self.parent(), -self._b, -self._a)

    def __mul__(self, other):
        assert self.parent()==other.parent(), "Must have same parents"
        l = [self._a*other._a, self._a*other._b, self._b*other._a, self._b*other._b]
        return ClosedInterval(self.parent(), min(l), max(l))

    def _pow_int(self, n):
        if n==0:
            return self.parent()(self.parent().base_ring().one())
        if n==1:
            return self
        if n>=2:
            if n%2 == 0:
                if self._a >= self.parent().base_ring().zero():
                    return ClosedInterval(self.parent(), self._a**n, self._b**n)
                if self._b <= self.parent().base_ring().zero():
                    return ClosedInterval(self.parent(), self._b**n, self._a**n)
                else:
                    # interval contains zero.
                    return ClosedInterval(self.parent(), self.parent().base_ring().zero(), max(self._a**n, self._b**n))
            else:
                # Fortunately raising to an odd power is monotonic.
                return ClosedInterval(self.parent(), self._a**n, self._b**n)
        else:
            # Negative power.
            return (~self)**(-n)

    def __contains__(self, value):
        return self._a <= value <= self._b

    def __invert__(self):
        if self.parent().base_ring().zero() in self:
            raise ZeroDivisionError("Inverting interval containing zero.")
        return ClosedInterval(self.parent(), ~self._b, ~self._a)

    def _div_(self, other):
        return self.__mul__(~other)

    def __eq__(self,other):
        if other is None:
            return False
        other = self.parent()(other)
        assert self.parent()==other.parent(), "Must have same parents"
        return self._a == other._a and self._b == other._b

    def __hash__(self):
        return 73*hash(self._a)-19*hash(self._b)

    def n(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a pair consisting of numerical approximations of the endpoints.

        EXAMPLES::

        sage: from sage.rings.real_transcendental import ClosedIntervalSet
        sage: CIS = ClosedIntervalSet(SR)
        sage: CIS(e,pi).n(digits=3)
        (2.72, 3.14)
        """
        return (self._a.n(prec, digits, algorithm), self._b.n(prec, digits, algorithm))

class ClosedIntervalSet(UniqueRepresentation, Parent):
    r"""
    The collection of closed bounded intervals with endpoints in a provided field.

    Interval arithmetic is supported.

    EXAMPLES::

    sage: from sage.rings.real_transcendental import ClosedIntervalSet
    sage: CIS = ClosedIntervalSet(QQ)
    sage: CIS(0)
    [0, 0]
    sage: CIS(-1/2,5)
    [-1/2, 5]

    sage: CIS(1,2)*CIS(2,3)
    [2, 6]
    sage: CIS(-3,2)^2
    [0, 9]
    sage: CIS(1,2)/CIS(-2,-1)
    [-2, -1/2]
    """
    Element = ClosedInterval

    def __init__(self, field):
        self._f = field
        Parent.__init__(self, category=Sets())

    def base_ring(self):
        return self._f

    def _element_constructor_(self, *args, **kwds):
        if len(args)==2:
            a = self._f(args[0])
            b = self._f(args[1])
            assert a <= b, "the first endpoint should be less than the second."
            return self.element_class(self, a, b)
        elif len(args)==1:
            if isinstance(args[0],ClosedInterval):
                a = self._f(args[0].lower())
                b = self._f(args[0].upper())
                return self.element_class(self, a, b)
            else:
                # A single point passed representing a degenerate interval.
                a = self._f(args[0])
                return self.element_class(self,a,a)

    def _an_element_(self):
        """
        Return [0,1].
        """
        return self.element_class(self, self._f.zero(), self._f.one())

    def some_elements(self):
        for a in self._f.some_elements():
            for b in self._f.some_elements():
                if a <= b:
                    yield self.element_class(self, a, b)

    def _repr_(self):
        return "The set of closed intervals with endpoints in "+str(self._f)

    def __contains__(self, x):
        if not isinstance(x, ClosedInterval):
            return False
        return x.parent()==self

class NestedIntervalReal(Element):
    r"""
    An instance of this class represents a real number described by a sequence of nested intervals,
    whose infinite intersection yields the represented real number.
    """
    def __init__(self, parent):
        Element.__init__(self, parent)

    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        pass

    def __repr__(self):
        return "NestedIntervalReal in "+repr(self.interval(0))

    def __add__(self, other):
        return NIR_sum(self.parent(), self, other)

    def __sub__(self, other):
        return NIR_difference(self.parent(), self, other)

    def __neg__(self):
        return NIR_negate(self.parent(), self)

    def __mul__(self, other):
        return NIR_product(self.parent(), self, other)

    def __invert__(self):
        return NIR_invert(self.parent(), self)

    def __div__(self, other):
        return self.__mul__(~(self.parent()(other)))

    def _pow_int(self, k):
        return NIR_power(self.parent(), self, k)

    def n(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a numerical approximation of "self" with "prec" bits (or
        decimal "digits") of precision.
        """
        if prec is None and digits is None:
            prec = 53
        if not prec is None:
            assert digits is None, "Only prec or digits can be provided, not both."
            # Figure out the power to use.
            k = 0
            while 0 in self.interval(k):
                k += 1
            I = self.interval(k)
            ratio = I.upper()/I.lower()
            while ratio > 2:
                k += 1
                I = self.interval(k)
                ratio = I.upper()/I.lower()
            if I.lower() > 0:
                # positive case.
                exp = floor(log(I.lower().n(),2))
            else:
                exp = floor(log((-I.upper()).n(),2))
            diameter_bound = self.parent().base_ring()(2)**(exp-prec)
        elif not digits is None:
            # Figure out the power to use.
            k = 0
            while 0 in self.interval(k):
                k += 1
            I = self.interval(k)
            ratio = I.upper()/I.lower()
            while ratio > 10:
                k += 1
                I = self.interval(k)
                ratio = I.upper()/I.lower()
            if I.lower() > 0:
                # positive case.
                exp = floor(log(I.lower().n(),10))
            else:
                exp = floor(log((-I.upper()).n(),10))
            diameter_bound = self.parent().base_ring()(10)**(exp-digits)
        I = self.interval(k)
        while I.diameter() > diameter_bound:
            k += 1
            I = self.interval(k)
        return I.center().n(prec=prec,digits=digits,algorithm=algorithm)

    def _richcmp_(self, other, op):
        #print "Executing comparison of "+str((self,other))
        if self is other:
            # Equal!
            return op == op_LE or op == op_EQ or op == op_GE
        else:
            # This will go into an infinite loop if other equals this.
            #print "checking with n=0"
            n = 0
            i = self.interval(n)
            j = other.interval(n)
            while i.intersects(j):
                #print "Found that these intersect: "+str((i,j))
                n += 1
                i = self.interval(n)
                j = other.interval(n)
            return richcmp(i.lower(), j.lower(), op)

class NestedIntervalRealSet(UniqueRepresentation, Parent):
    r"""
    The parent of NestedIntervalReal.
    """

    Element = NestedIntervalReal

    def __init__(self, field):
        self._f = field
        Parent.__init__(self, category=Sets())


    @cached_method
    def base_ring(self):
        r"""
        Return the field containing the endpoints of the intervals.
        """
        return self._f

    @cached_method
    def closed_interval_set(self):
        return ClosedIntervalSet(self._f)

    @cached_method
    def zero(self):
        return self(self._f.zero())

    @cached_method
    def one(self):
        return self(self._f.one())

    @cached_method
    def e(self):
        r"""
        Return the constant e as a NestedIntervalReal.

        EXAMPLES::

        sage: from sage.rings.real_transcendental import NestedIntervalRealSet
        sage: NIRS = NestedIntervalRealSet(QQ)
        sage: e = NIRS.e()
        sage: from sage.functions.other import factorial
        sage: n = 10
        sage: lower_bound = sum([1/factorial(i) for i in xrange(n+1)])
        sage: lower_bound < e
        True
        sage: e < lower_bound + 1/factorial(n)
        True
        """

        return E_as_NIR(self)

    def random(self, a, b):
        r"""
        Construct a random number between a and b.
        """
        return Random_NIR(self, self.closed_interval_set()(a,b))

    def _element_constructor_(self, *args, **kwds):
        assert len(args)==1, "Only takes one argument"
        if isinstance(args[0],NestedIntervalReal):
            if self == args[0].parent():
                return args[0]
            else:
                raise ValueError("Currently can not convert NestedIntervalReal to different field.")
        else:
            # Should be a constant
            c = self._f(args[0])
            return NIR_constant(self, c)

    def _coerce_map_from_(self, S):
        """
        The only thing that we except coersions from is the constant field.
        """
        if S is self._f:
            return True
        if self._f.has_coerce_map_from(S):
            return True
        return False


    def __repr__(self):
        return "NestedIntervalRealSet over "+repr(self._f)

class NIR_constant(NestedIntervalReal):
    r"""
    Represents a real number in the base_ring() of the parent NestedIntervalRealSet.
    """
    def __init__(self, parent, constant):
        self._c = parent.base_ring()(constant)
        NestedIntervalReal.__init__(self, parent)

    def interval(self, n):
        return self.parent().closed_interval_set()(self._c)

    def __repr__(self):
        return repr(self._c)

class NIR_sum(NestedIntervalReal):
    r"""
    Represents the sum of two NestedIntervalReals.
    """
    def __init__(self, *args):
        r"""
        Construct the sum of a list of NestedIntervalReals provided as
        arguments.
        """
        assert len(args)>2, "Require at least three arguments."
        parent = args[0]
        assert isinstance(parent,NestedIntervalRealSet)
        self._x = [parent(a) for a in args[1:]]
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return sum([val.interval(n) for val in self._x])

class NIR_difference(NestedIntervalReal):
    r"""
    Represents the difference of two NestedIntervalReals.
    """
    def __init__(self, parent, x, y):
        r"""
        Construct the difference `x-y` of two NestedIntervalReals.
        """
        assert isinstance(parent,NestedIntervalRealSet)
        self._x = parent(x)
        self._y = parent(y)
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return self._x.interval(n) - self._y.interval(n)


class NIR_negate(NestedIntervalReal):
    r"""
    Represents the negation of a NestedIntervalReal.
    """
    def __init__(self, parent, x):
        r"""
        Construct the negation `-x` of a NestedIntervalReal `x`.
        """
        assert isinstance(parent,NestedIntervalRealSet)
        self._x = parent(x)
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return -self._x.interval(n)

class NIR_product(NestedIntervalReal):
    r"""
    Represents the product of a list of NestedIntervalReals.
    """
    def __init__(self, *args):
        r"""
        Construct the products of a list of NestedIntervalReals provided
        as arguments.
        """
        assert len(args)>2, "Require at least three arguments."
        parent = args[0]
        assert isinstance(parent,NestedIntervalRealSet)
        self._x = [parent(a) for a in args[1:]]
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return prod([val.interval(n) for val in self._x])

class NIR_invert(NestedIntervalReal):
    r"""
    Construct the multiplicative inverse of a NestedIntervalReal.
    """
    def __init__(self, parent, x):
        r"""
        Construct `1/x`. Here x must be non-zero or this will result in an infinite loop.
        """
        assert isinstance(parent,NestedIntervalRealSet)
        self._x = parent(x)
        zero = parent.base_ring().zero()
        self._shift = 0
        while zero in self._x.interval(self._shift):
            self._shift += 1
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return ~self._x.interval(n+self._shift)

class NIR_power(NestedIntervalReal):
    r"""
    Construct an integer power of a NestedIntervalReal.
    """
    def __init__(self, parent, x, k):
        r"""
        Construct x^k.
        """
        self._x = x
        self._k = k
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        r"""
        For integers n>=0 return the n-th containing closed interval.
        """
        return (self._x.interval(n))**self._k

class E_as_NIR(NestedIntervalReal):
    r"""
    Represents the real number $e$ as a NestedIntervalReal.
    """
    def __init__(self, parent):
        self._shift = 1 # Some integer at least one.
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        f = self.parent().base_ring()
        total = f.zero()
        for i in xrange(n+1+self._shift):
            total += ~f(factorial(i))
        return self.parent().closed_interval_set()(total, total + ~f(factorial(n+self._shift)))

class Random_NIR(NestedIntervalReal):
    r"""
    Represents a random number in some closed interval as a NestedIntervalReal.

    This real number is taken with respect to Lebesgue measure on the interval.
    """
    def __init__(self, parent, interval, step=4096):
        r"""
        Represents a random number within the provided interval as a NestedIntervalReal.
        """
        self._step = step # Some integer at least one.
        self._n = 0
        interval = parent.closed_interval_set()(interval)
        assert not interval.is_degenerate(), "Degenerate interval provided."
        self._interval = interval
        NestedIntervalReal.__init__(self, parent)

    @cached_method
    def interval(self, n):
        if n==0:
            return self._interval
        assert n in ZZ and n >= 0, "n must be a positive integer"
        #print "Generating interval at level "+str(n)
        i = self.interval(n-1)
        a = i.lower()
        b = i.upper()
        r = randint(0,self._step - 1)
        self._n = n
        return self.parent().closed_interval_set()(a+r*(b-a)/self._step, a+(r+1)*(b-a)/self._step)

    def best_interval(self):
        return self.interval(self._n)

    def __repr__(self):
        i = self.best_interval()
        return "Random number within "+repr(i)

class RealTranscendentalExtensionFieldElement(FieldElement):
    r"""
    An element of a RealTranscendentalExtensionField.
    """
    def __init__(self, parent, rf, nir=None):
        r"""
        Construct an element in the field from a rational function.
        """
        try:
            # This .factor().expand() stuff is to deal with the fact that
            # FunctionField does not automatically cancel common factors:
            # sage: K.<x>=FunctionField(QQ)
            # sage: (2*x)/(2*(x+1))
            # 2*x/(2*x + 2)
            self._rf = parent.function_field()(rf.factor().expand())
        except ArithmeticError:
            self._rf = rf
        self._nir = nir
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
    def nested_interval_real(self):
        r"""
        Return the NestedIntervalReal given by evaluating the underlying rational function.
        """
        if self._nir is None:
            self._nir = self.parent()._eval(self._rf)
        return self._nir

    def __repr__(self):
        return repr(self._rf)

    def _add_(self, other):
        return RealTranscendentalExtensionFieldElement(self.parent(), self._rf+other._rf)

    def _sub_(self, other):
        return RealTranscendentalExtensionFieldElement(self.parent(), self._rf-other._rf)

    def _neg_(self):
        return RealTranscendentalExtensionFieldElement(self.parent(), -self._rf)

    def _mul_(self, other):
        other =  self.parent()(other)
        return RealTranscendentalExtensionFieldElement(self.parent(), self._rf*other._rf)

    def _invert_(self):
        return RealTranscendentalExtensionFieldElement(self.parent(), ~self._rf)

    def _div_(self, other):
        return RealTranscendentalExtensionFieldElement(self.parent(), self._rf/other._rf)

    def _pow_int(self, k):
        return RealTranscendentalExtensionFieldElement(self.parent(), self._rf**k)

    def numerical_approx(self, prec=None, digits=None, algorithm=None):
        r"""
        Return a numerical approximation of "self" with "prec" bits (or
        decimal "digits") of precision.
        """
        if self == self.parent().zero():
            return self.parent().constant_field().zero().numerical_approx(prec=prec, digits=digits, algorithm=algorithm)
        return self.nested_interval_real().n(prec=prec, digits=digits, algorithm=algorithm)

    def _mpfr_(self, R):
        """
        EXAMPLES::

        sage: from sage.rings.real_transcendental import NestedIntervalRealSet, \
            RealTranscendentalExtensionField
        sage: NIRS = NestedIntervalRealSet(QQ)
        sage: RTE = RealTranscendentalExtensionField(QQ, NIRS.e(), "ee")
        sage: ee = RTE.gen()
        sage: RF = RealField(200)
        sage: RF(e) == RF(ee)
        True
        """
        return R(self.numerical_approx(prec=R.prec()))

    def _real_double_(self, R):
        """
        EXAMPLES::

        sage: from sage.rings.real_transcendental import NestedIntervalRealSet, \
            RealTranscendentalExtensionField
        sage: NIRS = NestedIntervalRealSet(QQ)
        sage: RTE = RealTranscendentalExtensionField(QQ, NIRS.e(), "ee")
        sage: ee = RTE.gen()
        sage: RDF(e) == RDF(ee)
        True
        """
        return R(self.numerical_approx(prec=R.prec()))

    def __float__(self):
        """
        EXAMPLES::

        sage: from sage.rings.real_transcendental import NestedIntervalRealSet, \
            RealTranscendentalExtensionField
        sage: NIRS = NestedIntervalRealSet(QQ)
        sage: RTE = RealTranscendentalExtensionField(QQ, NIRS.e(), "ee")
        sage: ee = RTE.gen()
        sage: float(e) == float(ee)
        True
        """
        return float(self.numerical_approx(prec=53))


    def _richcmp_(self, other, op):
        if self._rf == other._rf:
            # Equal!
            return op == op_LE or op == op_EQ or op == op_GE
        else:
            return richcmp(self.nested_interval_real(),other.nested_interval_real(), op)

class RealTranscendentalExtensionField(Field):
    r"""
    Represents a real transcendental extension of a field in sage.

    The field we are extending is called the "constant field".

    EXAMPLES::

    sage: from sage.rings.real_transcendental import NestedIntervalRealSet, \
        RealTranscendentalExtensionField
    sage: NIR = NestedIntervalRealSet(QQ)
    sage: K = RealTranscendentalExtensionField(QQ, NIR.e(), "ee")
    sage: ee = K.gen()
    sage: ee.n(digits=5)
    2.7183

    sage: K = RealTranscendentalExtensionField(QQ,[NIR.e(),NIR.random(0,1)],["ee","r"])
    sage: ee,r = K.gens()
    sage: ee.n(digits=5)
    2.7183
    sage: 0<r<1
    True
    sage: r*(ee+1)/(r-2)^2
    ((ee + 1)*r)/(r^2 - 4*r + 4)
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
        self._NIRS = NestedIntervalRealSet(constant_field)
        try:
            n = len(transcendentals)
        except TypeError:
            # Just one transcendental
            n = -1 # Flag that we have already done the processing.

            self._values = tuple([self._NIRS(transcendentals)])
            if names is None:
                self._names = tuple(["t"])
            else:
                self._names = tuple([str(names)])
            #if latex_names is None:
            #    self._latex_names = self._names
            #else:
            #    self._latex_names = tuple([str(latex_names)])
        if n >= 0:
            if n == 0:
                raise ValueError("Must have at least one transcendental to extend.")
            self._values = tuple([self._NIRS(t) for t in transcendentals])
            if not names is None:
                self._names = tuple([str(names[i]) for i in xrange(n)])
                # We will not deal with latex_names yet.
                #if latex_names is None:
                #    self._latex_names = self._names
                #else:
                #    self._latex_names = tuple([str(latex_names[i]) for i in xrange(n)])
            else:
                self._names = tuple(["t_"+str(i) for i in xrange(n)])
                #if not latex_names is None:
                #    self._latex_names = tuple([str(latex_names[i]) for i in xrange(n)])
                #else:
                #    self._latex_names = tuple(["t_{"+str(i)+"}" for i in xrange(n)])
        else:
            # Fix the "flag"
            n=1

        # Construct the underlying function field.
        ff = self._constant_field
        gens = []
        for name in self._names:
            ff = FunctionField(ff,name)
            # Convert the generators that already exist into the current Field.
            gens = [ff(gen) for gen in gens]
            gens.append(ff.gen())

        self._ff = ff

        self._gens = tuple([self.element_class(self,gens[i], self._values[i]) for i in xrange(n)])

        Field.__init__(self, self._constant_field, category=Fields())

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

    def nested_interval_real_set(self):
        r"""
        Return the NestedIntervalSet over the constant_field.

        The NestedIntervalSet is used to approximate elements in the field and to make
        comparisons.
        """
        return self._NIRS

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
        return self.transcendental(n)**power

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

    def __repr__(self):
        if len(self._names)==1:
            return "Real transcendental extension of "+repr(self._constant_field)+" by "+str(self.gen())
        else:
            return "Real transcendental extension of "+repr(self._constant_field)+" by "+str(self.gens())

    def _element_constructor_(self, x):
        if isinstance(x,RealTranscendentalExtensionFieldElement):
            assert x.parent() == self
            return x
        else:
            xx = self._constant_field(x) # The quantity x should be in the constant field.
            return self.element_class(self, self._ff(xx), self._NIRS(xx))

    def _coerce_map_from_(self, S):
        r"""
        The only thing that we accept coersions from is the constant field.
        """
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

