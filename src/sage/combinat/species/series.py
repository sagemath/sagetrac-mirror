"""
Lazy Power Series

This file provides an implementation of lazy univariate power
series, which uses the stream class for its internal data
structure. The lazy power series keep track of their approximate
order as much as possible without forcing the computation of any
additional coefficients. This is required for recursively defined
power series.

This code is based on the work of Ralf Hemmecke and Martin Rubey's
Aldor-Combinat, which can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/aldor/combinat/index.html.
In particular, the relevant section for this file can be found at
http://www.risc.uni-linz.ac.at/people/hemmecke/AldorCombinat/combinatse9.html.
"""
#*****************************************************************************
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
#from stream import Stream, Stream_class
from new_stream import OldStreamBehavior as Stream, Stream as Stream_class, ListCachedStream, StreamFromList, StreamFromIterator
from series_order import  bounded_decrement, increment, inf, unk
from sage.rings.all import Integer, prod
from functools import partial
from sage.misc.misc import repr_lincomb, is_iterator
from sage.algebras.algebra import Algebra
from sage.algebras.algebra_element import AlgebraElement
import sage.structure.parent_base
from sage.categories.all import Rings

class SeriesStream(ListCachedStream):
    def __init__(self, order=unk, aorder=unk, base_ring=None,
                 aorder_changed=True, **kwds):
        assert base_ring is not None
        self._base_ring = base_ring
        self.aorder = aorder
        self.order = order
        if aorder == inf:
            self.order = inf
        self.aorder_changed = aorder_changed
        self._zero = base_ring(0)
        super(SeriesStream, self).__init__(**kwds)

    def __getitem__(self, n):
        # The following line must not be written n < self.get_aorder()
        # because comparison of Integer and OnfinityOrder is not implemented.
        if self.get_aorder() > n:
            return self._base_ring.zero_element()
        return super(SeriesStream, self).__getitem__(n)

    def children(self):
        return []

    def order_operation(self, *series):
        return 0

    def __mul__(self, other):
        return ProductStream(self, other, base_ring=self._base_ring)

    def stretch(self, k):
        """
        EXAMPLES::
        """
        return StretchedStream(k, self, base_ring=self._base_ring)

    def get_aorder(self):
        """
        Returns the approximate order of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: a._stream.get_aorder()
            1
        """
        if self.order is unk:
            self.refine_aorder()
        return self.aorder

    def get_order(self):
        """
        Returns the order of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: a._stream.get_order()
            1
        """
        if self.order is unk:
            self.refine_aorder()
        return self.order

    def refine_aorder(self):
        """
        Refines the approximate order of self as much as possible without
        computing any coefficients.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L(iter([0,0,0,0,1]))
            sage: a = a._stream
            sage: a.aorder
            Unknown series order
            sage: a[2]
            0
            sage: a.aorder
            0
            sage: a.refine_aorder()
            sage: a.aorder
            3

        ::

            sage: a = L(iter([0,0]))
            sage: a = a._stream
            sage: a.aorder
            Unknown series order
            sage: a[5]
            0
            sage: a.refine_aorder()
            sage: a.aorder
            Infinite series order

        ::

            sage: a = L(iter([0,0,1,0,0,0]))
            sage: a = a._stream
            sage: a[4]
            0
            sage: a.refine_aorder()
            sage: a.aorder
            2
        """
        #If we already know the order, then we don't have
        #to worry about the approximate order
        if self.order != unk:
            return

        #aorder can never be infinity since order would have to
        #be infinity as well
        assert self.aorder != inf

        if self.aorder == unk:
            self.compute_aorder()
        else:
            #Try to improve the approximate order
            ao = self.aorder
            n = self.number_computed()

            if self.aorder < n:
                while self.aorder < n:
                    if self._cache[self.aorder] == 0:
                        self.aorder += 1
                    else:
                        break

            #Try to recognize the zero series
            if self.aorder == n:
                #For non-constant series, we cannot do anything
                if not self.is_constant():
                    return
                elif self.get_constant() == 0:
                    if n >= self.get_constant_position():
                        self.set_approximate_order(inf)
                        return

            if self.aorder < n:
                self.order = self.aorder

    def compute_aorder(self):
        changed = self.aorder_changed
        
        ao = self.order_operation(*[c.aorder for c in self.children()])
        if ao == unk:
            ao = inf
        changed = self.set_approximate_order(ao) or changed
        
        if changed:
            for s in self.children():
                s.compute_aorder()
            ao = self.order_operation(*[c.aorder for c in self.children()])
            self.set_approximate_order(ao)

        if self.aorder == inf:
            self.order = inf

        if hasattr(self, '_reference'):
            self._reference._copy()

    def set_approximate_order(self, new_order):
        """
        Sets the approximate order of self and returns True if the
        approximate order has changed otherwise it will return False.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([0,0,0,3,2,1,0])
            sage: f = f._stream
            sage: f.get_aorder()
            3
            sage: f.set_approximate_order(3)
            False
        """
        self.aorder_changed = ( self.aorder != new_order )
        self.aorder = new_order
        if self.aorder == inf:
            self.order = inf
        return self.aorder_changed

class SeriesStreamFromList(SeriesStream, StreamFromList):
    def __init__(self, list=None, base_ring=None, **kwds):
        assert list is not None
        assert base_ring is not None
        list = map(base_ring, list)
        super(SeriesStreamFromList, self).__init__(list=list, base_ring=base_ring, **kwds)
        self.get_aorder()
        
    def compute_aorder(self):
        for i, value in enumerate(self._cache):
            if value != 0:
                ao = i
                break
        else:
            ao = inf
        self.set_approximate_order(ao)

class SeriesStreamFromIterator(SeriesStream, StreamFromIterator):
    def __init__(self, **kwds):
        super(SeriesStreamFromIterator, self).__init__(**kwds)

    def compute_aorder(self):
        if self._cache:
            for i, value in enumerate(self._cache):
                if value != self._zero:
                    break
            ao = i
        else:
            ao = 0
        self.set_approximate_order(ao)

class TermStream(SeriesStream):
    def __init__(self, n=None, value=None, **kwds):
        kwds['order'] = kwds['aorder'] = n if value != 0 else inf
        kwds['aorder_changed'] = False
        super(TermStream, self).__init__(**kwds)
        self._n = n
        self._value = self._base_ring(value)
        if value == 0:
            self.set_constant(0, self._zero)
        else:
            self.set_constant(n + 1, self._zero)

    def order_operation(self):
        return self.order
    
    def compute(self, n):
        if n == self._n:
            return self._value
        else:
            return self._zero

class ChangeRingStream(SeriesStream):
    def __init__(self, stream=None, new_ring=None, **kwds):
        self._stream = stream
        self._new_ring = new_ring
        super(ChangeRingStream, self).__init__(**kwds)

    def compute(self, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L2 = LazyPowerSeriesRing(RR)
            sage: a = L([1])
            sage: b = L2(a)
            sage: b.parent()
            Lazy Power Series Ring over Real Field with 53 bits of precision
            sage: b.coefficients(3)
            [1.00000000000000, 1.00000000000000, 1.00000000000000]
        """
        return self._new_ring(self._stream[n])

    def children(self):
        return [self._stream]

    def order_operation(self, ao):
        return ao

    def is_constant(self):
        return self._stream.is_constant()

    def get_constant(self):
        return self._stream.get_constant()

    def get_constant_position(self):
        return self._stream.get_constant_position()

class SumStream(SeriesStream):
    def __init__(self, left_summand=None, right_summand=None, **kwds):
        self._left = left_summand
        self._right = right_summand
        super(SumStream, self).__init__(**kwds)

    def children(self):
        return [self._left, self._right]

    order_operation = staticmethod(min)

    def compute(self, n):
        return self._left[n] + self._right[n]

    def is_constant(self):
        return self._left.is_constant() and self._right.is_constant()

    def get_constant(self):
        return self._left.get_constant() + self._right.get_constant()

    def get_constant_position(self):
        return max(self._left.get_constant_position(),
                   self._right.get_constant_position())

class ProductStream(SeriesStream):
    def __init__(self, left_factor=None, right_factor=None, **kwds):
        self._left = left_factor
        self._right = right_factor
        super(ProductStream, self).__init__(**kwds)

    def children(self):
        return [self._left, self._right]

    def order_operation(self, a, b):
        return a + b

    def compute(self, n):
        """
        Returns an iterator for the coefficients of self \* y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])
            sage: g = f * f
            sage: [g._stream.compute(i) for i in range(5)]
            [1, 2, 1, 0, 0]
        """
        low = self._left.aorder
        high = n - self._right.aorder
        res = self._base_ring(0)
        if low == inf or high == inf:
            return res

        for k in range(low, high + 1):
            cx = self._left[k]
            if cx == 0:
                continue
            res += cx * self._right[n - k]

        return res

    def is_constant(self):
        return ((self._left.is_constant() and self._left.is_constant() == 0) and
                (self._right.is_constant() and self._right.is_constant() == 0))

    def get_constant(self):
        assert self.is_constant()
        return self._base_ring(0)

    def get_constant_position(self):
        return (self._left.get_constant_position() *
                self._right.get_constant_position()) - 1
                
    
class TailStream(SeriesStream):
    def __init__(self, stream, **kwds):
        self._stream = stream
        super(TailStream, self).__init__(**kwds)
        
    def compute(self, n):
        return self._stream[n + 1]

    def children(self):
        return [self._stream]

    order_operation = staticmethod(bounded_decrement)

    def is_constant(self):
        return self._stream.is_constant()

    def get_constant(self):
        return self._stream.get_constant()

    def get_constant_position(self):
        return self._stream.get_constant_position() - 1

class DerivativeStream(SeriesStream):
    def __init__(self, stream, **kwds):
        self._stream = stream
        super(DerivativeStream, self).__init__(**kwds)

    def children(self):
        return [self._stream]

    order_operation = staticmethod(bounded_decrement)

    def compute(self, n):
        """
        Returns an iterator for the coefficients of the derivative of
        self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1])
            sage: g = f.derivative()
            sage: [g._stream.compute(i) for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        np1 = self._base_ring(n + 1)
        return np1 * self._stream[n + 1]

    def is_constant(self):
        return self._stream.is_constant() and self._stream.get_constant() == 0

    def get_constant(self):
        if self.is_constant():
            return self._base_ring(0)
        else:
            raise ValueError

    def get_constant_position(self):
        return self._stream.get_constant_position() - 1

class IntegralStream(SeriesStream):
    def __init__(self, stream, integration_constant=0, **kwds):
        if integration_constant != 0:
            kwds['order'] = kwds['aorder'] = 0
            kwds['aorder_changed'] = False
        super(IntegralStream, self).__init__(**kwds)
        self._stream = stream
        self._ic = self._base_ring(integration_constant)

    def children(self):
        return [self._stream]

    def order_operation(self, a):
        if self._ic == 0:
            return a + 1
        else:
            return 0
        
    def compute(self, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L.gen()
            sage: g = s.integral(0)
            sage: [g._stream.compute(i) for i in range(5)]
            [0, 0, 1/2, 0, 0]

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([0,0,4,5,6,0]).derivative()
            sage: g = f.integral(1)
            sage: [g._stream.compute(i) for i in range(5)]
            [1, 0, 4, 5, 6]
        """
        if n == 0:
            return self._ic
        else:
            return (Integer(1) / Integer(n)) * self._stream[n - 1]

    def is_constant(self):
        return self._stream.is_constant() and self._stream.get_constant() == 0

    def get_constant(self):
        if self.is_constant():
            return self._base_ring(0)
        else:
            raise ValueError

    def get_constant_position(self):
        return self._stream.get_constant_position() + 1

class CompositionStream(SeriesStreamFromIterator):
    def __init__(self, outer_stream, inner_stream, **kwds):
        self._outer = outer_stream
        self._inner = inner_stream
        super(CompositionStream, self).__init__(iterator=self.compute_iterator(), **kwds)

    def children(self):
        return [self._outer, self._inner]
    
    def order_operation(self, a, b):
        return a*b

    def compute_iterator(self):
        """
        Returns a iterator for the coefficients of the composition of this
        power series with the power series y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: t = L([0,1])
            sage: g = s(t)
            sage: [g[i] for i in range(10)]
            [1, 1, 2, 4, 8, 16, 32, 64, 128, 256]
        """
        assert self._inner[0] == 0
        yield self._outer[0]

        res = TailStream(self._outer, base_ring=self._base_ring)
        res = CompositionStream(res, self._inner, base_ring=self._base_ring)
        res = ProductStream(res, self._inner, base_ring=self._base_ring)
        self._stream = res

        self._stream[0]
        n = 1
        while True:
            yield self._stream[n]
            n += 1

class RecursiveStream(SeriesStream):
    def define(self, stream):
        self._stream = stream
        self._stream._reference = self
        self._copy()

    def order_operation(self):
        return self._stream.aorder

    def _copy(self):
        self.aorder = self._stream.aorder
        self.order = self._stream.order
        self.aorder_changed = self._stream.aorder_changed

    def refine_aorder(self):
        try:
            return self._stream.refine_aorder()
        except AttributeError:
            raise NotImplementedError
        self._copy()

    def compute_aorder(self):
        try:
            return self._stream.compute_aorder()
        except AttributeError:
            raise NotImplementedError
        self._copy()
    
    def compute(self, n):
        try:
            return self._stream[n]
        except AttributeError:
            raise NotImplementedError

class RestrictedStream(SeriesStream):
    def __init__(self, stream, min=None, max=None, **kwds):
        self._min = min
        self._max = max
        self._stream = stream
        super(RestrictedStream, self).__init__(**kwds)

    def children(self):
        return [self._stream]

    def order_operation(self, a):
        return max(a, self._min)

    def compute(self, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1])
            sage: g = a.restricted(min=2)
            sage: [g._stream.compute(i) for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: g = a.restricted(min=2, max=4)
            sage: [g._stream.compute(i) for i in range(10)]
            [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]

        ::

            sage: g = a.restricted(min=2, max=5)
            sage: [g._stream.compute(i) for i in range(6)]
            [0, 0, 1, 1, 1, 0]
        """
        if self._min is not None and n < self._min:
            return self._zero
        if self._max is not None and n >= self._max:
            return self._zero
        return self._stream[n]

class ListSumStream(SeriesStream):
    def __init__(self, stream_list, **kwds):
        self._stream_list = stream_list
        super(ListSumStream, self).__init__(**kwds)

    def children(self):
        return self._stream_list

    def order_operation(self, *orders):
        return min(orders)
    
    def compute(self, n):
        """
        Returns a generator for the coefficients of the sum the the lazy
        power series in series_list.

        INPUT:


        -  ``series_list`` - a list of lazy power series


        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: series_list = [ L([1]), L([0,1]), L([0,0,1]) ]
            sage: g = L.sum(series_list)
            sage: [g._stream.compute(i) for i in range(5)]
            [1, 2, 3, 3, 3]
        """
        return sum([c[n] for c in self.children()], self._zero)

    def is_constant(self):
        return all(c.is_constant() for c in self.children())

    def get_constant(self):
        return sum(c.get_constant() for c in self.children())

    def get_constant_position(self):
        return max(c.get_constant_position() for c in self.children())

class SumGeneratorStream(SeriesStream):
    def __init__(self, series_stream, **kwds):
        self._series_stream = SeriesStreamFromIterator(iterator=iter(series_stream), **kwds)
        super(SumGeneratorStream, self).__init__(**kwds)

    def compute(self, n):
        s = self._series_stream
        r = s[n][n]
        for i in range(min(n, s.number_computed() - 1)):
            r += s[i][n]
        return r

class ProductGeneratorStream(SeriesStreamFromIterator):
    def __init__(self, series_iter, **kwds):
        self._series_it = iter(series_iter)
        super(ProductGeneratorStream, self).__init__(iterator=self.compute_iterator(), **kwds)

    def compute_iterator(self):
        z = self._series_it.next()
        yield z[0]
        yield z[1]

        n = 2
        for x in self._series_it:
            z = z * x
            yield z[n]
            n += 1

        while True:
            yield z[n]
            n += 1

class PowerStream(StreamFromIterator):
    def __init__(self, stream):
        self._stream = stream
        super(PowerStream, self).__init__(iterator=self.compute_iterator())

    def compute_iterator(self):
        k = 1
        z = self._stream
        while True:
            yield z
            z = ProductStream(z, self._stream, base_ring=z._base_ring)

class StretchedStream(SeriesStream):
    def __init__(self, k, stream, **kwds):
        self._k = k
        self._stream = stream
        super(StretchedStream, self).__init__(**kwds)

    def compute(self, n):
        n = ZZ(n)
        quo, rem = n.quo_rem(self._k)
        if rem == 0:
            return self._stream[quo]
        else:
            return self._zero

class LazyPowerSeriesRing(Algebra):
    def __init__(self, R, element_class = None, names=None):
        """
        TESTS::

            sage: from sage.combinat.species.series import LazyPowerSeriesRing
            sage: L = LazyPowerSeriesRing(QQ)
            sage: loads(dumps(L))
            Lazy Power Series Ring over Rational Field
        """
        #Make sure R is a ring with unit element
        if not R in Rings():
            raise TypeError, "Argument R must be a ring."
        try:
            z = R(Integer(1))
        except StandardError:
            raise ValueError, "R must have a unit element"

        #Take care of the names
        if names is None:
            names = 'x'
        else:
            names = names[0]

        self._element_class = element_class if element_class is not None else LazyPowerSeries
        self._order = None
        self._name = names
        sage.structure.parent_base.ParentWithBase.__init__(self, R)

    def ngens(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ).ngens()
            1
        """
        return 1

    def __repr__(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ)
            Lazy Power Series Ring over Rational Field
        """
        return "Lazy Power Series Ring over %s"%self.base_ring()

    def __cmp__(self, x):
        """
        EXAMPLES::

            sage: LQ = LazyPowerSeriesRing(QQ)
            sage: LZ = LazyPowerSeriesRing(ZZ)
            sage: LQ == LQ
            True
            sage: LZ == LQ
            False
        """
        if self.__class__ is not x.__class__:
            return cmp(self.__class__, x.__class__)
        return cmp(self.base_ring(), x.base_ring())

    def _coerce_impl(self, x):
        """
        EXAMPLES::

            sage: L1 = LazyPowerSeriesRing(QQ)
            sage: L2 = LazyPowerSeriesRing(RR)
            sage: L2.has_coerce_map_from(L1)
            True
            sage: L1.has_coerce_map_from(L2)
            False

        ::

            sage: a = L1([1]) + L2([1])
            sage: a.coefficients(3)
            [2.00000000000000, 2.00000000000000, 2.00000000000000]
        """
        return self(x)

    def _new(self, stream_cls, *args, **kwds):
        kwds['base_ring'] = self.base_ring()
        return self._element_class(self, stream_cls(*args, **kwds))
        
    def __call__(self, x=None, order=unk):
        """
        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: L()
            O(1)
            sage: L(1)
            1
            sage: L(ZZ).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: L(iter(ZZ)).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: L(Stream(ZZ)).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]

        ::

            sage: a = L([1,2,3])
            sage: a.coefficients(3)
            [1, 2, 3]
            sage: L(a) is a
            True
            sage: L_RR = LazyPowerSeriesRing(RR)
            sage: b = L_RR(a)
            sage: b.coefficients(3)
            [1.00000000000000, 2.00000000000000, 3.00000000000000]
            sage: L(b)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to coerce ... into self

        TESTS::

            sage: L(pi)
            Traceback (most recent call last):
            ...
            TypeError: do not know how to coerce pi into self
        """
        cls = self._element_class
        base_ring = self.base_ring()

        if x is None:
            return self._new(RecursiveStream)

        if isinstance(x, LazyPowerSeries):
            x_parent = x.parent()
            if x_parent.__class__ != self.__class__:
                raise ValueError

            if x_parent.base_ring() == base_ring:
                return x
            else:
                if base_ring.has_coerce_map_from(x_parent.base_ring()):
                    return self._new(ChangeRingStream, stream=x._stream,
                                     new_ring=base_ring)


        if hasattr(x, "parent") and base_ring.has_coerce_map_from(x.parent()):
            x = base_ring(x)
            return self.term(x, 0)

        if isinstance(x, (list, tuple)):
            x = SeriesStreamFromList(x, base_ring=base_ring)
        elif hasattr(x, "__iter__") and not isinstance(x, Stream_class):
            x = iter(x)

        if is_iterator(x):
            x = SeriesStreamFromIterator(iterator=x, base_ring=base_ring)

        if isinstance(x, SeriesStream):
            return cls(self, x)
        elif not hasattr(x, "parent"):
            return self.term(x, 0)

        raise TypeError, "do not know how to coerce %s into self"%x

    def zero_element(self):
        """
        Returns the zero power series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.zero_element()
            0
        """
        return self(self.base_ring()(0))

    def identity_element(self):
        """
        Returns the one power series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.identity_element()
            1
        """
        return self(self.base_ring()(1))

    def gen(self, i=0):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.gen().coefficients(5)
            [0, 1, 0, 0, 0]
        """
        res = self._new(TermStream, 1, 1)
        res.set_name(self._name)
        return res

    def term(self, r, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.term(0,0)
            0
            sage: L.term(3,2).coefficients(5)
            [0, 0, 3, 0, 0]
        """
        if n < 0:
            raise ValueError, "n must be non-negative"
        BR = self.base_ring()
        if r == 0:
            res = self._new(TermStream, 0, 0)
            res.set_name("0")
        else:
            res = self._new(TermStream, n, r)

            if n == 0:
                res.set_name(repr(r))
            elif n == 1:
                res.set_name(repr(r) + "*" + self._name)
            else:
                res.set_name("%s*%s^%s"%(repr(r), self._name, n))

        return res

    def sum(self, a):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: l = [L(ZZ)]*3
            sage: L.sum(l).coefficients(10)
            [0, 3, -3, 6, -6, 9, -9, 12, -12, 15]
        """
        return self._new(ListSumStream, [x._stream for x in a])

    #Potentially infinite sum
    def _sum_generator_gen(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: def f():
            ...       while True:
            ...           yield s
            sage: g = L._sum_generator_gen(f())
            sage: [g.next() for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        s = Stream(g)
        n = 0
        while True:
            r = s[n].coefficient(n)
            for i in range(min(n, len(s)-1)):
                r += s[i].coefficient(n)
            yield r
            n += 1

    def sum_generator(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = [L([1])]*6 + [L(0)]
            sage: t = L.sum_generator(g)
            sage: t.coefficients(10)
            [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]

        ::

            sage: s = L([1])
            sage: def g():
            ...       while True:
            ...           yield s
            sage: t = L.sum_generator(g())
            sage: t.coefficients(9)
            [1, 2, 3, 4, 5, 6, 7, 8, 9]
        """
        return self._new(SumGeneratorStream, g)

    #Potentially infinite product
    def product_generator(self, g):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s1 = L([1,1,0])
            sage: s2 = L([1,0,1,0])
            sage: s3 = L([1,0,0,1,0])
            sage: s4 = L([1,0,0,0,1,0])
            sage: s5 = L([1,0,0,0,0,1,0])
            sage: s6 = L([1,0,0,0,0,0,1,0])
            sage: s = [s1, s2, s3, s4, s5, s6]
            sage: def g():
            ...       for a in s:
            ...           yield a
            sage: p = L.product_generator(g())
            sage: p.coefficients(26)
            [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 4, 4, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0]

        ::

            sage: def m(n):
            ...       yield 1
            ...       while True:
            ...           for i in range(n-1):
            ...               yield 0
            ...           yield 1
            ...
            sage: def s(n):
            ...       q = 1/n
            ...       yield 0
            ...       while True:
            ...           for i in range(n-1):
            ...               yield 0
            ...           yield q
            ...

        ::

            sage: def lhs_gen():
            ...       n = 1
            ...       while True:
            ...           yield L(m(n))
            ...           n += 1
            ...

        ::

            sage: def rhs_gen():
            ...       n = 1
            ...       while True:
            ...           yield L(s(n))
            ...           n += 1
            ...
            sage: lhs = L.product_generator(lhs_gen())
            sage: rhs = L.sum_generator(rhs_gen()).exponential()
            sage: lhs.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
            sage: rhs.coefficients(10)
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30]
        """
        return self._new(ProductGeneratorStream, g)

class LazyPowerSeries(AlgebraElement):
    def __init__(self, A, stream=None):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L()
            sage: loads(dumps(f))
            O(1)
        """
        AlgebraElement.__init__(self, A)
        self._stream = stream
        self._zero = A.base_ring().zero_element()
        self._name = None

    def set_name(self, name):
        self._name = name

    @property
    def aorder(self):
        return self._stream.aorder

    def get_aorder(self):
        return self._stream.get_aorder()

    @property
    def order(self):
        return self._stream.order

    def get_order(self):
        return self._stream.get_order()

    def _get_repr_info(self, x):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1,2,3])
            sage: a.compute_coefficients(5)
            sage: a._get_repr_info('x')
            [('1', 1), ('x', 2), ('x^2', 3)]
        """
        n = len(self._stream)
        m = ['1', x]
        m += [x+"^"+str(i) for i in range(2, n)]
        c = [ self._stream[i] for i in range(n) ]
        return [ (m,c) for m,c in zip(m,c) if c != 0]

    def __repr__(self):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L(); s.set_name('s'); s
            s

        ::

            sage: L()
            O(1)

        ::

            sage: a = L(iter([1,2,3]))
            sage: a
            O(1)
            sage: a.compute_coefficients(2)
            sage: a
            1 + 2*x + 3*x^2 + O(x^3)
            sage: a.compute_coefficients(4)
            sage: a
            1 + 2*x + 3*x^2 + 3*x^3 + 3*x^4 + 3*x^5 + ...

        ::

            sage: a = L([1,2,3,0])
            sage: a.compute_coefficients(5)
            sage: a
            1 + 2*x + 3*x^2
        """
        if self._name is not None:
            return self._name

        n = len(self._stream)
        x = self.parent()._name
        baserepr = repr_lincomb(self._get_repr_info(x))
        if (self._stream.is_constant() and
            (n - 1) >= self._stream.get_constant_position()):
            if self._stream[n-1] == 0:
                l = baserepr
            else:
                l = baserepr + " + " + repr_lincomb([(x+"^"+str(i), self._stream[n-1]) for i in range(n, n+3)]) + " + ..."
        else:
            l = baserepr + " + O(x^%s)"%n if n > 0 else "O(1)"
        return l


    def compute_coefficients(self, i):
        """
        Computes all the coefficients of self up to i.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1,2,3])
            sage: a.compute_coefficients(5)
            sage: a
            1 + 2*x + 3*x^2 + 3*x^3 + 3*x^4 + 3*x^5 + ...
        """
        self.coefficient(i)

    def coefficients(self, n):
        """
        Returns the first n coefficients of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,2,3,0])
            sage: f.coefficients(5)
            [1, 2, 3, 0, 0]
        """
        return [self.coefficient(i) for i in range(n)]

    def is_zero(self):
        """
        Returns True if and only if self is zero.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([0,2,3,0])
            sage: s.is_zero()
            False

        ::

            sage: s = L(0)
            sage: s.is_zero()
            True

        ::

            sage: s = L(iter([0]))
            sage: s.is_zero()
            False
            sage: s.coefficient(0)
            0
            sage: s.coefficient(1)
            0
            sage: s.is_zero()
            True
        """
        self._stream.refine_aorder()
        return self._stream.order == inf

    def define(self, x):
        """
        EXAMPLES: Test Recursive 0

        ::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: one = L(1)
            sage: monom = L.gen()
            sage: s = L()
            sage: s.define(one+monom*s)
            sage: s.aorder
            Unknown series order
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(6)]
            [1, 1, 1, 1, 1, 1]

        Test Recursive 1

        ::

            sage: s = L()
            sage: s.define(one+monom*s*s)
            sage: s.aorder
            Unknown series order
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(6)]
            [1, 1, 2, 5, 14, 42]

        Test Recursive 1b

        ::

            sage: s = L()
            sage: s.define(monom + s*s)
            sage: s.aorder
            Unknown series order
            sage: s.order
            Unknown series order
            sage: [s.coefficient(i) for i in range(7)]
            [0, 1, 1, 2, 5, 14, 42]

        Test Recursive 2

        ::

            sage: s = L()
            sage: t = L()
            sage: s.define(one+monom*t*t*t)
            sage: t.define(one+monom*s*s)
            sage: [s.coefficient(i) for i in range(9)]
            [1, 1, 3, 9, 34, 132, 546, 2327, 10191]
            sage: [t.coefficient(i) for i in range(9)]
            [1, 1, 2, 7, 24, 95, 386, 1641, 7150]

        Test Recursive 2b

        ::

            sage: s = L()
            sage: t = L()
            sage: s.define(monom + t*t*t)
            sage: t.define(monom + s*s)
            sage: [s.coefficient(i) for i in range(9)]
            [0, 1, 0, 1, 3, 3, 7, 30, 63]
            sage: [t.coefficient(i) for i in range(9)]
            [0, 1, 1, 0, 2, 6, 7, 20, 75]

        Test Recursive 3

        ::

            sage: s = L()
            sage: s.define(one+monom*s*s*s)
            sage: [s.coefficient(i) for i in range(10)]
            [1, 1, 3, 12, 55, 273, 1428, 7752, 43263, 246675]
        """
        self._stream.define(x._stream)

    def _new(self, cls, *args, **kwds):
        parent = self.parent()
        kwds['base_ring'] = parent.base_ring()
        return self.__class__(parent, cls(*args, **kwds))

    def coefficient(self, n):
        """
        Returns the coefficient of xn in self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(ZZ)
            sage: [f.coefficient(i) for i in range(5)]
            [0, 1, -1, 2, -2]
        """
        return self._stream[n]

    def get_stream(self):
        """
        Returns self's underlying Stream object.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L.gen()
            sage: s = a.get_stream()
            sage: [s[i] for i in range(5)]
            [0, 1, 0, 0, 0]
        """
        self._stream.refine_aorder()
        return self._stream

    def _add_(self, y):
        """
        EXAMPLES: Test Plus 1

        ::

            sage: from sage.combinat.species.series import *
            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: gs0 = L([0])
            sage: gs1 = L([1])
            sage: sum1 = gs0 + gs1
            sage: sum2 = gs1 + gs1
            sage: sum3 = gs1 + gs0
            sage: [gs0.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            sage: [gs1.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [sum1.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: [sum2.coefficient(i) for i in range(11)]
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
            sage: [sum3.coefficient(i) for i in range(11)]
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

        Test Plus 2

        ::

            sage: gs1 = L([1,2,4,8,0])
            sage: gs2 = L([-1, 0,-1,-9,22,0])
            sage: sum = gs1 + gs2
            sage: sum2 = gs2 + gs1
            sage: [ sum.coefficient(i) for i in range(5) ]
            [0,  2, 3, -1, 22]
            sage: [ sum.coefficient(i) for i in range(5, 11) ]
            [0, 0, 0, 0, 0, 0]
            sage: [ sum2.coefficient(i) for i in range(5) ]
            [0,  2, 3, -1, 22]
            sage: [ sum2.coefficient(i) for i in range(5, 11) ]
            [0, 0, 0, 0, 0, 0]
        """
        return self._new(SumStream, self._stream, y._stream)

    add = _add_

    def _mul_(self, y):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: gs0 = L(0)
            sage: gs1 = L([1])

        ::

            sage: prod0 = gs0 * gs1
            sage: [prod0.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        ::

            sage: prod1 = gs1 * gs0
            sage: [prod1.coefficient(i) for i in range(11)]
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        ::

            sage: prod2 = gs1 * gs1
            sage: [prod2.coefficient(i) for i in range(11)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

        ::

            sage: gs1 = L([1,2,4,8,0])
            sage: gs2 = L([-1, 0,-1,-9,22,0])

        ::

            sage: prod1 = gs1 * gs2
            sage: [prod1.coefficient(i) for i in range(11)]
            [-1, -2, -5, -19, 0, 0, 16, 176, 0, 0, 0]

        ::

            sage: prod2 = gs2 * gs1
            sage: [prod2.coefficient(i) for i in range(11)]
            [-1, -2, -5, -19, 0, 0, 16, 176, 0, 0, 0]
        """
        return self._new(ProductStream, self._stream, y._stream)

    times = _mul_

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])  # 1+x
            sage: g = f^3
            sage: g.coefficients(4)
            [1, 3, 3, 1]

        ::

            sage: f^0
            1
        """
        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError, "n must be a nonnegative integer"
        if n == 0:
            return self.parent().term(1, 0)
        
        #Compute power using binary powering
        res = self
        for i in reversed(Integer(n).bits()[:-1]):
            res = res * res
            if i == 1:
                res = res * self
        return res

    CompositionStream = CompositionStream
            
    def __call__(self, y):
        """
        Returns the composition of this power series and the power series
        y.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: t = L([0,0,1])
            sage: u = s(t)
            sage: u.coefficients(11)
            [1, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34]

        Test Compose 2

        ::

            sage: s = L([1])
            sage: t = L([0,0,1,0])
            sage: u = s(t)
            sage: u.aorder
            Unknown series order
            sage: u.order
            Unknown series order
            sage: u.coefficients(10)
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
            sage: u.aorder
            0
            sage: u.order
            0

        Test Compose 3 s = 1/(1-x), t = x/(1-x) s(t) = (1-x)/(1-2x)

        ::

            sage: s = L([1])
            sage: t = L([0,1])
            sage: u = s(t)
            sage: u.coefficients(14)
            [1, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
        """
        return self._new(self.CompositionStream, self._stream, y._stream)

    composition = __call__

    def tail(self):
        """
        Returns the power series whose coefficients obtained by subtracting
        the constant term from this series and then dividing by x.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(range(20))
            sage: g = f.tail()
            sage: g.coefficients(10)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        return self._new(TailStream, self._stream)

    def iterator(self, n=0, initial=None):
        """
        Returns an iterator for the coefficients of self starting at n.

        EXAMPLES::

            sage: from sage.combinat.species.stream import Stream
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(range(10))
            sage: g = f.iterator(2)
            sage: [g.next() for i in range(5)]
            [2, 3, 4, 5, 6]
            sage: g = f.iterator(2, initial=[0,0])
            sage: [g.next() for i in range(5)]
            [0, 0, 2, 3, 4]
        """
        if initial is not None:
            for x in initial:
                yield x
        while True:
            yield self._stream[n]
            n += 1

    compose = __call__

    def _power_gen(self):
        """
        Returns a generator for all the powers self^k starting with k = 1.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,1,0])
            sage: g = f._power_gen()
            sage: g.next().coefficients(5)
            [1, 1, 0, 0, 0]
            sage: g.next().coefficients(5)
            [1, 2, 1, 0, 0]
            sage: g.next().coefficients(5)
            [1, 3, 3, 1, 0]
        """
        k = 1
        z = self
        while True:
            yield z
            z = z*self

    def derivative(self):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: one = L(1)
            sage: monom = L.gen()
            sage: s = L([1])
            sage: u = s.derivative()
            sage: u.coefficients(10)
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

        ::

            sage: s = L()
            sage: s._name = 's'
            sage: s.define(one+monom*s*s)
            sage: u = s.derivative()
            sage: u.coefficients(5) #[1*1, 2*2, 3*5, 4*14, 5*42]
            [1, 4, 15, 56, 210]

        ::

            sage: s = L([1])
            sage: t = L([0,1])
            sage: u = s(t).derivative()
            sage: v = (s.derivative().compose(t))*t.derivative()
            sage: u.coefficients(11)
            [1, 4, 12, 32, 80, 192, 448, 1024, 2304, 5120, 11264]
            sage: v.coefficients(11)
            [1, 4, 12, 32, 80, 192, 448, 1024, 2304, 5120, 11264]

        ::

            sage: s = L();
            sage: t = L();
            sage: s.define(monom+t*t*t)
            sage: t.define(monom+s*s)
            sage: u = (s*t).derivative()
            sage: v = s.derivative()*t + s*t.derivative()
            sage: u.coefficients(10)
            [0, 2, 3, 4, 30, 72, 133, 552, 1791, 4260]
            sage: v.coefficients(10)
            [0, 2, 3, 4, 30, 72, 133, 552, 1791, 4260]
            sage: u.coefficients(10) == v.coefficients(10)
            True

        ::

            sage: f = L([0,0,4,5,6,0])
            sage: d = f.derivative()
            sage: d.get_aorder()
            1
            sage: d.coefficients(5)
            [0, 8, 15, 24, 0]
        """
        return self._new(DerivativeStream, self._stream)

    ###########
    #Integrals#
    ###########
    def integral(self, integration_constant = 0):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: zero = L(0)
            sage: s = zero
            sage: t = s.integral()
            sage: t.is_zero()
            True

        ::

            sage: s = zero
            sage: t = s.integral(1)
            sage: t.coefficients(6)
            [1, 0, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1, 0)
            sage: t = s.integral()
            sage: t.coefficients(6)
            [0, 1, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1,0)
            sage: t = s.integral(1)
            sage: t.coefficients(6)
            [1, 1, 0, 0, 0, 0]
            sage: t._stream.is_constant()
            True

        ::

            sage: s = L.term(1, 4)
            sage: t = s.integral()
            sage: t.coefficients(10)
            [0, 0, 0, 0, 0, 1/5, 0, 0, 0, 0]

        ::

            sage: s = L.term(1,4)
            sage: t = s.integral(1)
            sage: t.coefficients(10)
            [1, 0, 0, 0, 0, 1/5, 0, 0, 0, 0]

        TESTS::

            sage: f = L([0,0,4,5,6,0])
            sage: i = f.derivative().integral()
            sage: i.get_aorder()
            2
            sage: i.coefficients(5)
            [0, 0, 4, 5, 6]
            sage: i = f.derivative().integral(1)
            sage: i.get_aorder()
            0
            sage: i.coefficients(5)
            [1, 0, 4, 5, 6]
        """
        return self._new(IntegralStream, stream=self._stream,
                         integration_constant=integration_constant)

    def is_finite(self, n=None):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L(iter([0,0,1,0,0])); a
            O(1)
            sage: a.is_finite()
            False
            sage: c = a[4]
            sage: a.is_finite()
            False
            sage: a.is_finite(4)
            False
            sage: c = a[5]
            sage: a.is_finite()
            True
            sage: a.is_finite(4)
            True
        """
        if self.order is inf:
            return True

        s = self._stream

        if n is None:
            n = len(s)

        if s.is_constant() and all(s[i] == 0 for i in range(n-1, max(n,len(s)))):
            return True

        return False

    def exponential(self):
        """
        TESTS::

            sage: def inv_factorial():
            ...       q = 1
            ...       yield 0
            ...       yield q
            ...       n = 2
            ...       while True:
            ...           q = q / n
            ...           yield q
            ...           n += 1
            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L(inv_factorial()) #e^(x)-1
            sage: u = f.exponential()
            sage: g = inv_factorial()
            sage: z1 = [1,1,2,5,15,52,203,877,4140,21147,115975]
            sage: l1 = [z*g.next() for z in z1]
            sage: l1 = [1] + l1[1:]
            sage: u.coefficients(11)
            [1, 1, 1, 5/6, 5/8, 13/30, 203/720, 877/5040, 23/224, 1007/17280, 4639/145152]
            sage: l1 == u.coefficients(11)
            True
        """
        base_ring = self.parent().base_ring()
        s = self.parent()()
        s.define( (self.derivative()*s).integral(base_ring(1)) )
        return s

    def __getitem__(self, i):
        """
        Returns the ith coefficient of self.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: f = L([1,2,3,0])
            sage: [f[i] for i in range(5)]
            [1, 2, 3, 0, 0]
        """
        return self.coefficient(i)


    #########################
    #Min and max restriction#
    #########################
    def restricted(self, min=None, max=None):
        """
        Returns the power series restricted to the coefficients starting at
        min and going up to, but not including max. If min is not
        specified, then it is assumed to be zero. If max is not specified,
        then it is assumed to be infinity.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: a = L([1])
            sage: a.restricted().coefficients(10)
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: a.restricted(min=2).coefficients(10)
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: a.restricted(max=5).coefficients(10)
            [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
            sage: a.restricted(min=2, max=6).coefficients(10)
            [0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
        """
        if ((min is None and max is None) or
            (max is None and self._stream.get_aorder() >= min)):
            return self
        return self._new(RestrictedStream, self._stream, min=min, max=max)
