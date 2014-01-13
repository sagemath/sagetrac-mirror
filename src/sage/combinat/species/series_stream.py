"""
Series Streams

"""
#*****************************************************************************
#       Copyright (C) 2013 Mike Hansen <mhansen@gmail.com>,
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
from new_stream import ListCachedStream, StreamFromList, StreamFromIterator
from series_order import bounded_decrement, inf, unk
from sage.rings.all import Integer

class SeriesStream(ListCachedStream):
    def __init__(self, order=unk, aorder=unk, base_ring=None,
                 aorder_changed=True, convert=False, **kwds):
        assert base_ring is not None
        self._base_ring = base_ring
        self.aorder = aorder
        self.order = order
        if aorder == inf:
            self.order = inf
        self.aorder_changed = aorder_changed
        self._zero = base_ring(0)
        self._convert = convert
        self._children = kwds.pop('children', [])
        super(SeriesStream, self).__init__(**kwds)

    def __getitem__(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.  This method
        returns zero without doing further computations if ``n`` is
        less than the approximate order.
        """
        # The following line must not be written n < self.get_aorder()
        # because comparison of Integer and InfinityOrder is not implemented.
        if self.get_aorder() > n:
            return self._zero
        result = super(SeriesStream, self).__getitem__(n)
        if self._convert:
            return self._base_ring(result)
        else:
            return result

    def base_ring(self):
        return self._base_ring

    def children(self):
        return self._children

    def order_operation(self, *series):
        if self.aorder != unk:
            return self.aorder
        else:
            return 0

    def __mul__(self, other):
        """
        Returns the ProductStream of ``self`` and ``other``.

        """
        return ProductStream(self, other, base_ring=self._base_ring)

    def stretch(self, k):
        """
        EXAMPLES::
        """
        return StretchedStream(k, self, base_ring=self._base_ring)

    def get_aorder(self):
        """
        Returns the approximate order for this stream after refining
        it (without computing additional coefficients).

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[0,1,0], base_ring=QQ)
            sage: s.get_aorder()
            1
        """
        if self.order is unk:
            self.refine_aorder()
        return self.aorder

    def get_order(self):
        """
        Returns the order of self.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[0,1,0], base_ring=QQ)
            sage: s.get_order()
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

            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator
            sage: a = SeriesStreamFromIterator(iterator=iter([0,0,0,0,1]), base_ring=QQ)
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

            sage: a = SeriesStreamFromIterator(iterator=iter([0,0]), base_ring=QQ)
            sage: a.aorder
            Unknown series order
            sage: a[5]
            0
            sage: a.refine_aorder()
            sage: a.aorder
            Infinite series order

        ::

            sage: a = SeriesStreamFromIterator(iterator=iter([0,0,1,0,0,0]), base_ring=QQ)
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
        """
        Computes the approximate order of this stream based on its
        children.

        EXAMPLES:

        """
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

    def set_approximate_order(self, new_order):
        """
        Sets the approximate order of self and returns True if the
        approximate order has changed otherwise it will return False.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: f = SeriesStreamFromList(list=[0,0,0,3,2,1,0], base_ring=QQ)
            sage: f.get_aorder()
            3
            sage: f.set_approximate_order(3)
            False
            sage: f.set_approximate_order(2)
            True
            sage: f.set_approximate_order(3)
            True
            
        """
        self.aorder_changed = ( self.aorder != new_order )
        self.aorder = new_order
        if self.aorder == inf:
            self.order = inf
        return self.aorder_changed

class SeriesStreamFromList(SeriesStream, StreamFromList):
    def __init__(self, **kwds):
        if 'list' in kwds:
            kwds['list'] = map(kwds['base_ring'], kwds['list'])
        super(SeriesStreamFromList, self).__init__(**kwds)
        self.get_aorder() # compute the approximate order right away

    def order_operation(self):
        """
        Returns the order of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
        """
        for i, value in enumerate(self._cache):
            if value != 0:
                return i
        else:
            return inf

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
        super(ChangeRingStream, self).__init__(children=[stream], **kwds)

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

    def order_operation(self, ao):
        return ao

    def is_constant(self):
        return self._stream.is_constant()

    def get_constant(self):
        return self._new_ring(self._stream.get_constant())

    def get_constant_position(self):
        return self._stream.get_constant_position()

class SumStream(SeriesStream):
    def __init__(self, left_summand=None, right_summand=None, **kwds):
        self._left = left_summand
        self._right = right_summand
        super(SumStream, self).__init__(children=[self._left, self._right], **kwds)

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
        super(ProductStream, self).__init__(children=[self._left, self._right],
                                            **kwds)

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
        return ((self._left.is_constant() and self._left.get_constant() == 0) and
                (self._right.is_constant() and self._right.get_constant() == 0))

    def get_constant(self):
        assert self.is_constant()
        return self._base_ring(0)

    def get_constant_position(self):
        return (self._left.get_constant_position() *
                self._right.get_constant_position())
                
    
class TailStream(SeriesStream):
    def __init__(self, stream, **kwds):
        self._stream = stream
        super(TailStream, self).__init__(children=[stream], **kwds)
        
    def compute(self, n):
        return self._stream[n + 1]

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
        super(DerivativeStream, self).__init__(children=[stream], **kwds)

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
        super(IntegralStream, self).__init__(children=[stream], **kwds)
        self._stream = stream
        self._ic = self._base_ring(integration_constant)

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

class CompositionStream(SeriesStream):
    def __init__(self, outer_stream, inner_stream, **kwds):
        """
        A stream for thethe coefficients of the composition of
        ``outer_stream`` with ``inner_stream``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([1])
            sage: t = L([0,1])
            sage: g = s(t)
            sage: stream = g._stream; stream
            <class 'sage.combinat.species.series_stream.CompositionStream'>
            sage: [stream[i] for i in range(10)]
            [1, 1, 2, 4, 8, 16, 32, 64, 128, 256]
        """

        self._outer = outer_stream
        self._inner = inner_stream
        super(CompositionStream, self).__init__(children=[self._outer, self._inner],
                                                **kwds)
        
    def order_operation(self, a, b):
        """
        The order of a composition of two streams is the product of
        their orders.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([0,0,0,1])
            sage: t = L([0,0,1])
            sage: g = s(t)
            sage: stream = g._stream; stream
            <class 'sage.combinat.species.series_stream.CompositionStream'>
            sage: stream.order_operation(s.get_order(), t.get_order())
            6
        """
        return a*b

    def recursive_stream(self):
        """
        Returns a stream whose tail is equal to the tail of this
        composition stream.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([2])
            sage: t = L([0,3])
            sage: g = s(t)
            sage: stream = g._stream; stream
            <class 'sage.combinat.species.series_stream.CompositionStream'>
            sage: [stream[i] for i in range(1, 5)]
            [6, 24, 96, 384]
            sage: f = stream.recursive_stream()
            sage: [f[i] for i in range(1, 5)]
            [6, 24, 96, 384]

        """
        res = TailStream(self._outer, base_ring=self._base_ring)
        res = CompositionStream(res, self._inner, base_ring=self._base_ring)
        res = ProductStream(res, self._inner, base_ring=self._base_ring)
        return res

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this composition stream.
        If $n = 0$, then it will return the $0^{th}$ coefficient of
        ``self._outer``; otherwise, it will return the $n^th$ coefficient
        of the stream defined by :meth:`recursive_stream`.

        EXAMPLES::
        
            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L([2])
            sage: t = L([0,3])
            sage: g = s(t)
            sage: stream = g._stream; stream
            <class 'sage.combinat.species.series_stream.CompositionStream'>
            sage: [stream.compute(i) for i in range(1, 5)]
            [6, 24, 96, 384]
        """
        if n == 0:
            assert self._inner[0] == 0
            return self._outer[0]
        
        try:
            res = self._res
        except AttributeError:
            res = self._res = self.recursive_stream()
            res[0]
        return res[n]

class RecursiveStream(SeriesStream):
    """
    A stream used for recursively defined series.
    """
    def define(self, stream):
        """
        Defines this stream to be equal to ``stream``.  The stream
        ``stream`` will be referred to as the definition stream.
        """
        self._stream = stream

    def order_operation(self):
        """
        The order of this recursive stream is the order of the
        definition stream.
        """
        return self.aorder

    def not_defined_check(func):
        """
        A decorator which checks to see if the definition stream as
        been set.  If so, the decorated function executes as normal;
        otherwise, a ``NotImplementedError`` is raised.
        """
        from functools import wraps
        @wraps(func)
        def wrapper(self, *args, **kwds):
            try:
                self._stream
            except AttributeError:
                raise NotImplementedError, "must call define() first'"
            return func(self, *args, **kwds)
        return wrapper

    @property
    @not_defined_check
    def aorder(self):
        """
        Returns the approximate order of the definition stream.
        """
        return self._stream.aorder

    @aorder.setter
    def aorder(self, value):
        """
        Setting the approximate order on this stream does nothing.
        """
        pass

    @property
    @not_defined_check
    def order(self):
        """
        Returns the order of the definition stream.
        """
        return self._stream.order

    @order.setter
    def order(self, value):
        """
        Setting the order on this stream does nothing.
        """
        pass

    @property
    @not_defined_check
    def aorder_changed(self):
        """
        Returns whether or not the approximate order on the definition
        stream has changed.
        """
        return self._stream.aorder_changed

    @aorder_changed.setter
    def aorder_changed(self, value):
        """
        Setting the ``aorder_changed`` flag on this stream does nothing.
        """
        pass

    @not_defined_check
    def refine_aorder(self):
        """
        Refines the approximate order of the definition stream.
        """
        return self._stream.refine_aorder()

    @not_defined_check
    def compute_aorder(self):
        """
        Computes the approximate order of the definition stream.
        """
        return self._stream.compute_aorder()

    @not_defined_check
    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream, which is the
        same as the $n^{th}$ coefficient of the definition stream.
        """
        return self._stream[n]

class RestrictedStream(SeriesStream):
    def __init__(self, stream, min=None, max=None, **kwds):
        """
        Returns a stream whose $n^{th}$ coefficient is the $n^{th}$
        coefficient of ``stream`` unless ``n < min`` or ``n >= max``
        in which case it is $0$.  If either ``min`` or ``max`` are
        ``None`` then the previous corresponding conditions are
        considered ``False``.
        """
        self._min = min
        self._max = max
        self._stream = stream
        super(RestrictedStream, self).__init__(children=[stream], **kwds)

    def order_operation(self, a):
        """
        The order of the a restricted stream is the maximum of
        ``self._min`` and the order of the underlying stream.
        """
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
        """
        Returns a stream whose $n^{th}$ coefficient is the sum of the
        $n^{th}$ coefficients of the streams in ``stream_list``.
        These streams will be referred to as "child streams".
        """
        self._stream_list = stream_list
        super(ListSumStream, self).__init__(children=stream_list, **kwds)

    def order_operation(self, *orders):
        """
        The order of a :class:`ListSumStream` is the minimum of all
        the child streams.
        """
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
        """
        Returns whether or not this stream is constant.
        """
        return all(c.is_constant() for c in self.children())

    def get_constant(self):
        """
        """
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
        n = Integer(n)
        quo, rem = n.quo_rem(self._k)
        if rem == 0:
            return self._stream[quo]
        else:
            return self._zero
