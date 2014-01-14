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
        """
        A class for streams that represent power series.  The class
        keeps track of the (approximate) order of the power series
        which allows for recursively defined streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStream
            sage: class MyStream(SeriesStream):
            ....:     def compute(self, n):
            ....:         return 4
            sage: s = MyStream(order=3, aorder=3, base_ring=QQ, aorder_changed=False, convert=True)
            sage: [s[i] for i in range(6)]
            [0, 0, 0, 4, 4, 4]
            sage: parent(_[0])
            Rational Field
        """
        assert base_ring is not None
        self._base_ring = base_ring
        self.aorder = aorder
        self.order = order
        if aorder == inf:
            self.order = inf
        self.aorder_changed = aorder_changed
        self._zero = base_ring(0)
        self._convert = convert
        self._children = tuple(kwds.pop('children', ()))
        super(SeriesStream, self).__init__(**kwds)

    def __getitem__(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.  This method
        returns zero without doing further computations if ``n`` is
        less than the approximate order.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: [s[i] for i in range(10)]
            [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
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
        """
        Returns the base ring of this series stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[1,2,3], base_ring=ZZ)
            sage: s.base_ring()
            Integer Ring
        """
        return self._base_ring

    def children(self):
        """
        Returns the children of this series stream.  If this series
        stream depends on other streams, then those streams are
        considered its "children".

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[1,2,3], base_ring=ZZ)
            sage: s.children()
            ()
            sage: s2 = s * s
            sage: s2.children() == (s, s)
            True
        """
        return self._children

    def order_operation(self, *orders):
        """
        Returns the approximate order of this stream give the
        approximate orders of its children.  For streams with no
        children and unknown approximate order, it returns ``0``.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[1,2,3], base_ring=ZZ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            0
        """
        if self.aorder != unk:
            return self.aorder
        else:
            return 0

    def __mul__(self, other):
        """
        Returns the product of two streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[1,2,3,0], base_ring=ZZ)
            sage: s2 = s * s
            sage: [s2[i] for i in range(6)]
            [1, 4, 10, 12, 9, 0]
        """
        return ProductStream(self, other, base_ring=self._base_ring)

    def __add__(self, other):
        """
        Returns the sum of two streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[1,2,3,0], base_ring=ZZ)
            sage: s2 = s + s
            sage: [s2[i] for i in range(6)]
            [2, 4, 6, 0, 0, 0]
        """
        return SumStream(self, other, base_ring=self._base_ring)

    def stretch(self, k):
        """
        Returns a this stream stretched by $k$.  See
        :class:`StretchedStream` for a defintion.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1,2,3,], base_ring=ZZ)
            sage: s = t.stretch(3)
            sage: [s[i] for i in range(12)]
            [1, 0, 0, 2, 0, 0, 3, 0, 0, 3, 0, 0]
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
        Returns the order of this stream after refining it (without
        computing additional coefficients).

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
                        self.set_order(inf)
                        return

            if self.aorder < n:
                self.order = self.aorder

    def compute_aorder(self):
        """
        Computes the approximate order of this stream based on its
        children.

        EXAMPLES:

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.compute_aorder()
            sage: s.aorder
            1
        """
        changed = self.aorder_changed
        
        # Compute the approximate order from the children
        ao = self.order_operation(*[c.aorder for c in self.children()])

        # If we still don't know the approximate order, set the
        # approximate order to infinity and recompute from there.
        if ao == unk:
            ao = inf

        # If the approximate order of this series has changed (which
        # includes when we set it to infinity), then recompute the
        # approximate orders of the children and update our
        # approximate order
        changed = self.set_approximate_order(ao) or changed
        if changed:
            for s in self.children():
                s.compute_aorder()
            ao = self.order_operation(*[c.aorder for c in self.children()])
            self.set_approximate_order(ao)

        if self.aorder == inf:
            self.order = inf

    def set_order(self, order):
        """
        Sets the order of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: f = SeriesStreamFromList(list=[3,2,1,0], base_ring=QQ)
            sage: f.set_order(1)
            sage: [f[i] for i in range(5)]
            [0, 2, 1, 0, 0]
        """
        self.aorder = order
        self.aorder_changed = False
        self.order = order

    def set_approximate_order(self, new_order):
        """
        Sets the approximate order of this stream and returns ``True``
        if the approximate order has changed otherwise it will return
        ``False``.

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
        return self.aorder_changed

class SeriesStreamFromList(SeriesStream, StreamFromList):
    def __init__(self, **kwds):
        """
        A class for streams where the coefficients of the stream are
        given by a list.  Either the ``list`` keyword can be specified
        in ``__init__``, or a subclass can define a ``list`` method to
        determine the list of coefficients.   See
        :class:`sage.combinat.species.new_stream.StreamFromList`.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: f = SeriesStreamFromList(list=[0,0,0,3,2,1,0], base_ring=QQ)
            sage: [f[i] for i in range(10)]
            [0, 0, 0, 3, 2, 1, 0, 0, 0, 0]

            sage: class MyStream(SeriesStreamFromList):
            ....:     def list(self):
            ....:         return [1, 2, 3]
            sage: f = MyStream(base_ring=QQ)
            sage: [f[i] for i in range(5)]
            [1, 2, 3, 3, 3]
            sage: parent(_[0])
            Rational Field
        """
        super(SeriesStreamFromList, self).__init__(**kwds)
        self._cache = map(self._base_ring, self._cache)
        self.get_aorder() # compute the approximate order right away

    def set_constant(self, n, value):
        """
        Sets this stream to be constant at position $n$. We override
        :func:`StreamFromList.set_constant` in order to make sure that
        the constant is in the proper ring.


        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: f = SeriesStreamFromList(list=[0,0,0,3,2,1,0], base_ring=QQ)
            sage: f.set_constant(f.get_constant(), 3.0)
            sage: a = f[10]; a
            3
            sage: parent(a)
            Rational Field
        """
        super(SeriesStreamFromList, self).set_constant(n, self._base_ring(value))

    def order_operation(self):
        """
        Returns the order of this stream.  Since all of the
        coefficients are known at creation time, we can go through
        them to determine the order.
        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: s = SeriesStreamFromList(list=[0,1,2,3], base_ring=QQ)
            sage: s.order_operation()
            1
        """
        for i, value in enumerate(self._cache):
            if value != self._zero:
                return i
        else:
            return inf

class SeriesStreamFromIterator(SeriesStream, StreamFromIterator):
    """
    A class for streams where the coefficients of the stream are
    given by an iterator.  See
    :class:`sage.combinat.species.new_stream.StreamFromIterator`.

    EXAMPLES::

        sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator
        sage: s = SeriesStreamFromIterator(iterator=iter([1,2,3]), base_ring=QQ, convert=True)
        sage: [s[i] for i in range(5)]
        [1, 2, 3, 3, 3]
        sage: parent(_[0])
        Rational Field

    """
    def compute_aorder(self):
        """
        Computes the approximate order of this stream by examining the
        coefficients already produced by the iterator.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator
            sage: s = SeriesStreamFromIterator(iterator=iter([0,1,2,3]), base_ring=QQ, convert=True)
            sage: s[1]
            1
            sage: s.aorder
            0
            sage: s.compute_aorder()
            sage: s.aorder
            1
        """
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
        """
        A class for streams with at most one non-zero coefficient.

        INPUT::

        - ``n`` (int) the location of the nonzero coefficient.

        - ``value`` - the value of the nonzero coefficient.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: t.order
            1
            sage: [t[i] for i in range(4)]
            [0, 10, 0, 0]
            sage: t = TermStream(n=1, value=0, base_ring=QQ)
            sage: t.order
            Infinite series order
            sage: [t[i] for i in range(4)]
            [0, 0, 0, 0]
        """
        kwds['order'] = kwds['aorder'] = n if value != 0 else inf
        kwds['aorder_changed'] = False
        super(TermStream, self).__init__(**kwds)
        self._n = n
        self._value = self._base_ring(value)
        if value == self._zero:
            self.set_constant(0, self._zero)
        else:
            self.set_constant(n + 1, self._zero)

    def order_operation(self):
        """
        Returns the order of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: t.order_operation()
            1
        """
        return self.order
    
    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: t.compute(0)
            0
            sage: t.compute(1)
            10
        """
        if n == self._n:
            return self._value
        else:
            return self._zero

class ChangeRingStream(SeriesStream):
    def __init__(self, stream=None, new_ring=None, **kwds):
        """
        A class for streams which changes the base ring of another
        stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: [s[i] for i in range(4)]
            [0.000000000000000, 10.0000000000000, 0.000000000000000, 0.000000000000000]
        """
        self._stream = stream
        self._new_ring = new_ring
        super(ChangeRingStream, self).__init__(children=[stream], **kwds)

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.
        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: [s.compute(i) for i in range(4)]
            [0.000000000000000, 10.0000000000000, 0.000000000000000, 0.000000000000000]
        """
        return self._new_ring(self._stream[n])

    def order_operation(self, ao):
        """
        Returns the approximate order of this stream, which is the
        approximate order of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            1
        """
        return ao

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.  It
        is only eventually constant if the underlying stream is.

        EXAMPLES::
        
            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream, SeriesStreamFromIterator
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: s.is_constant() and t.is_constant()
            True
            sage: t = SeriesStreamFromIterator(iterator=ZZ, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: s.is_constant()
            False
        """
        return self._stream.is_constant()

    def get_constant(self):
        """
        If the stream is eventually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: s.get_constant()
            0.000000000000000
        """
        return self._new_ring(self._stream.get_constant())

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ChangeRingStream
            sage: t = TermStream(n=1, value=10, base_ring=QQ)
            sage: s = ChangeRingStream(stream=t, new_ring=RR, base_ring=RR)
            sage: s.get_constant_position()
            2
        """
        return self._stream.get_constant_position()

class SumStream(SeriesStream):
    def __init__(self, left_summand=None, right_summand=None, **kwds):
        """
        A class for a stream whose $n^{th}$ coefficient is the sum of
        the $n^{th}$ coefficient of two other streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, SumStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t2, base_ring=QQ)
            sage: [s[i] for i in range(6)]
            [0, 10, 5, 0, 0, 0]
        """
        self._left = left_summand
        self._right = right_summand
        super(SumStream, self).__init__(children=[self._left, self._right], **kwds)

    order_operation = staticmethod(min)

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.

        EXAMPLES::
        
            sage: from sage.combinat.species.series_stream import TermStream, SumStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t2, base_ring=QQ)
            sage: [s.compute(i) for i in range(6)]
            [0, 10, 5, 0, 0, 0]
        """
        return self._left[n] + self._right[n]

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::
        
            sage: from sage.combinat.species.series_stream import TermStream, SumStream, SeriesStreamFromIterator
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t2, base_ring=QQ)
            sage: s.is_constant()
            True
            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator
            sage: t3 = SeriesStreamFromIterator(iterator=ZZ, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t3, base_ring=QQ)
            sage: s.is_constant()
            False
        """
        return all(c.is_constant() for c in self.children())

    def get_constant(self):
        """
        If the stream is eventually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, SumStream, SeriesStreamFromIterator
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t2, base_ring=QQ)
            sage: s.get_constant()
            0
        """
        return self._left.get_constant() + self._right.get_constant()

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, SumStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = SumStream(left_summand=t1, right_summand=t2, base_ring=QQ)
            sage: s.get_constant_position()
            3
        """
        return max(self._left.get_constant_position(),
                   self._right.get_constant_position())

class ProductStream(SeriesStream):
    def __init__(self, left_factor=None, right_factor=None, **kwds):
        """
        A class for a stream for the which represents the product of
        two power series.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import ProductStream, SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1], base_ring=QQ)
            sage: s = ProductStream(left_factor=t, right_factor=t, base_ring=QQ)
            sage: [s[i] for i in range(6)]
            [1, 2, 3, 4, 5, 6]
        """
        self._left = left_factor
        self._right = right_factor
        super(ProductStream, self).__init__(children=[self._left, self._right],
                                            **kwds)

    def order_operation(self, a, b):
        """
        Returns the approximate order of this stream, which is the sum
        of the approximate orders of the child streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ProductStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = ProductStream(left_factor=t1, right_factor=t2, base_ring=QQ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            3
        """
        return a + b

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this product stream.
        Currently, this just performs the naive O(n^2) multiplication.

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
        """
        Returns whether or not this stream is eventually constant.

        ..note ::

          :class:`ProductStream` only recognizes if the stream is
          eventually constant when the constant is zero.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ProductStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = ProductStream(left_factor=t1, right_factor=t2, base_ring=QQ)
            sage: s.is_constant()
            True
        """
        return ((self._left.is_constant() and self._left.get_constant() == 0) and
                (self._right.is_constant() and self._right.get_constant() == 0))

    def get_constant(self):
        """
        If the stream is eventually constant, returns the constant
        value.

        ..note ::

          :class:`ProductStream` only recognizes if the stream is
          eventually constant when the constant is zero.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ProductStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = ProductStream(left_factor=t1, right_factor=t2, base_ring=QQ)
            sage: s.get_constant()
            0        
        """
        assert self.is_constant()
        return self._base_ring(0)

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, ProductStream
            sage: t1 = TermStream(n=1, value=10, base_ring=QQ)
            sage: t2 = TermStream(n=2, value=5, base_ring=QQ)
            sage: s = ProductStream(left_factor=t1, right_factor=t2, base_ring=QQ)
            sage: s.get_constant_position()
            6
        """
        return (self._left.get_constant_position() *
                self._right.get_constant_position())
                
    
class TailStream(SeriesStream):
    def __init__(self, stream, **kwds):
        """
        A class for a stream which is the tail of another stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator, TailStream
            sage: t = SeriesStreamFromIterator(iterator=NN, base_ring=QQ)
            sage: s = TailStream(t, base_ring=QQ)
            sage: [t[i] for i in range(7)]
            [0, 1, 2, 3, 4, 5, 6]
            sage: [s[i] for i in range(6)]
            [1, 2, 3, 4, 5, 6]
        """
        self._stream = stream
        super(TailStream, self).__init__(children=[stream], **kwds)
        
    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream which is the
        $(n+1)^{th}$ coefficient of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator, TailStream
            sage: t = SeriesStreamFromIterator(iterator=NN, base_ring=QQ)
            sage: s = TailStream(t, base_ring=QQ)
            sage: [s.compute(i) for i in range(5)]
            [1, 2, 3, 4, 5]
        """
        return self._stream[n + 1]

    order_operation = staticmethod(bounded_decrement)

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, TailStream
            sage: t = TermStream(n=2, value=3, base_ring=QQ)
            sage: s = TailStream(t, base_ring=QQ)
            sage: s.is_constant()
            True
        """
        return self._stream.is_constant()

    def get_constant(self):
        """
        If the stream is eventually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, TailStream
            sage: t = TermStream(n=2, value=3, base_ring=QQ)
            sage: s = TailStream(t, base_ring=QQ)
            sage: s.get_constant()
            0
        """
        return self._stream.get_constant()

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, TailStream
            sage: t = TermStream(n=2, value=3, base_ring=QQ)
            sage: s = TailStream(t, base_ring=QQ)
            sage: s.get_constant_position()
            2
        """
        return self._stream.get_constant_position() - 1

class DerivativeStream(SeriesStream):
    def __init__(self, stream, **kwds):
        """
        A class for a stream whose coefficients represent the
        derivative of a power series.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, DerivativeStream
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        self._stream = stream
        super(DerivativeStream, self).__init__(children=[stream], **kwds)

    order_operation = staticmethod(bounded_decrement)

    def compute(self, n):
        """
        Returns an iterator for the coefficients of the derivative of
        self.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, DerivativeStream
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        """
        np1 = self._base_ring(n + 1)
        return np1 * self._stream[n + 1]

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, DerivativeStream
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: s.is_constant()
            False
            sage: t = SeriesStreamFromList(list=[1,2,3,0], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: s.is_constant()
            True
        """
        return self._stream.is_constant() and self._stream.get_constant() == 0

    def get_constant(self):
        """
        If this stream is eventually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, DerivativeStream
            sage: t = SeriesStreamFromList(list=[1, 2, 3, 0], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: s.get_constant()
            0        
        """
        assert self.is_constant()
        return self._zero

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, DerivativeStream
            sage: t = SeriesStreamFromList(list=[1, 2, 3, 0], base_ring=ZZ)
            sage: s = DerivativeStream(t, base_ring=ZZ)
            sage: s.get_constant_position()
            2
        """
        return self._stream.get_constant_position() - 1

class IntegralStream(SeriesStream):
    def __init__(self, stream, integration_constant=0, **kwds):
        """
        A class for a stream whose coefficients represent the
        derivative of a power series.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [0, 1, 1, 1, 0, 0, 0, 0, 0, 0]
            sage: s = IntegralStream(t, integration_constant=3, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [3, 1, 1, 1, 0, 0, 0, 0, 0, 0]
        """
        if integration_constant != 0:
            kwds['order'] = kwds['aorder'] = 0
            kwds['aorder_changed'] = False
        super(IntegralStream, self).__init__(children=[stream], **kwds)
        self._stream = stream
        self._ic = self._base_ring(integration_constant)

    def order_operation(self, a):
        """
        Returns the approximate order of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[0, 1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            2
       """
        if self._ic == 0:
            return a + 1
        else:
            return 0
        
    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.
        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[0, 1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [0, 0, 1/2, 2/3, 3/4, 0, 0, 0, 0, 0]
        """
        if n == 0:
            return self._ic
        else:
            return (Integer(1) / Integer(n)) * self._stream[n - 1]

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[0, 1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: s.is_constant()
            True
        """
        return (self._stream.is_constant() and
                self._stream.get_constant() == 0)

    def get_constant(self):
        """
        If this stream is evetnually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[0, 1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: s.get_constant()
            0
        """
        assert self.is_constant()
        return self._zero

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, IntegralStream
            sage: t = SeriesStreamFromList(list=[0, 1, 2, 3, 0], base_ring=ZZ)
            sage: s = IntegralStream(t, base_ring=ZZ)
            sage: s.get_constant_position()
            5
        """
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

    def tail_stream(self):
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
            sage: f = stream.tail_stream()
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
        of the stream defined by :meth:`tail_stream`.

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
            res = self._res = self.tail_stream()
            res[0]
        return res[n]

class RecursiveStream(SeriesStream):
    """
    A stream used for recursively defined series.

    EXAMPLES::

        sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
        sage: t = TermStream(n=1, value=1, base_ring=ZZ)
        sage: s = RecursiveStream(base_ring=ZZ)
        sage: s.define(t + s*s)
        sage: [s[i] for i in range(10)]
        [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
    """
    def define(self, stream):
        """
        Defines this stream to be equal to ``stream``.  The stream
        ``stream`` will be referred to as the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: [s[i] for i in range(10)]
            [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
        """
        self._stream = stream

    def order_operation(self):
        """
        The order of this recursive stream is the order of the
        definition stream.
        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s[2]
            1
            sage: s.order_operation(*[c for c in s.children()])
            1
        """
        return self.aorder

    def not_defined_check(func):
        """
        A decorator which checks to see if the definition stream as
        been set.  If so, the decorated function executes as normal;
        otherwise, a ``NotImplementedError`` is raised.

        EXAMPLES::
        
            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.aorder # indirect doctest
            Traceback (most recent call last):
            ...
            ValueError: must call define() first
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.aorder
            Unknown series order
            sage: s[5]
            14
            sage: s.aorder
            1            
        """
        from functools import wraps
        @wraps(func)
        def wrapper(self, *args, **kwds):
            try:
                self._stream
            except AttributeError:
                raise ValueError, "must call define() first"
            return func(self, *args, **kwds)
        return wrapper

    @property
    @not_defined_check
    def aorder(self):
        """
        Returns the approximate order of the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s[5]
            14
            sage: s.aorder
            1
        """
        return self._stream.aorder

    @aorder.setter
    def aorder(self, value):
        """
        Setting the approximate order on this stream does nothing.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s[5]
            14
            sage: s.aorder = 10
            sage: s.aorder
            1
        """
        pass

    @property
    @not_defined_check
    def order(self):
        """
        Returns the order of the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s[5]
            14
            sage: s.order
            1
        """
        return self._stream.order

    @order.setter
    def order(self, value):
        """
        Setting the order on this stream does nothing.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s[5]
            14
            sage: s.order = 10
            sage: s.order
            1        
        """
        pass

    @property
    @not_defined_check
    def aorder_changed(self):
        """
        Returns whether or not the approximate order on the definition
        stream has changed.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.aorder_changed
            True
            sage: s[5]
            14
            sage: s.aorder_changed
            False        
        """
        return self._stream.aorder_changed

    @aorder_changed.setter
    def aorder_changed(self, value):
        """
        Setting the ``aorder_changed`` flag on this stream does nothing.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.aorder_changed = False
            sage: s.aorder_changed
            True
        """
        pass

    @not_defined_check
    def refine_aorder(self):
        """
        Refines the approximate order of the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.aorder
            Unknown series order
            sage: s.refine_aorder()
            sage: s.aorder
            1
        """
        return self._stream.refine_aorder()

    @not_defined_check
    def compute_aorder(self):
        """
        Computes the approximate order of the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: s.compute_aorder()
            sage: s.aorder
            1
        """
        return self._stream.compute_aorder()

    @not_defined_check
    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream, which is the
        same as the $n^{th}$ coefficient of the definition stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream, RecursiveStream
            sage: t = TermStream(n=1, value=1, base_ring=ZZ)
            sage: s = RecursiveStream(base_ring=ZZ)
            sage: s.define(t + s*s)
            sage: [s.compute(i) for i in range(10)]
            [0, 1, 1, 2, 5, 14, 42, 132, 429, 1430]
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

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, RestrictedStream
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = RestrictedStream(t, min=2, max=6, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [0, 0, 1, 1, 1, 1, 0, 0, 0, 0]
        """
        self._min = min
        self._max = max
        self._stream = stream
        super(RestrictedStream, self).__init__(children=[stream], **kwds)

    def order_operation(self, a):
        """
        The order of the a restricted stream is the maximum of
        ``self._min`` and the order of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, RestrictedStream
            sage: t = SeriesStreamFromList(list=[0,0,0,1], base_ring=ZZ)
            sage: s = RestrictedStream(t, min=2, max=6, base_ring=ZZ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            3
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = RestrictedStream(t, min=2, max=6, base_ring=ZZ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            2
        """
        return max(a, self._min)

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.
        
        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, RestrictedStream
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = RestrictedStream(t, min=2, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [0, 0, 1, 1, 1, 1, 1, 1, 1, 1]
            sage: s = RestrictedStream(t, min=2, max=4, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [0, 0, 1, 1, 0, 0, 0, 0, 0, 0]
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

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [0, 1, 3, 6, 6, 6, 6, 6, 6, 6]
        """
        assert stream_list != []
        self._stream_list = stream_list
        super(ListSumStream, self).__init__(children=stream_list, **kwds)

    def order_operation(self, *orders):
        """
        The order of a :class:`ListSumStream` is the minimum of all
        the child streams.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            1
        """
        return min(orders)
    
    def compute(self, n):
        """
        Computes the $n^{th}$ coefficient of this sum.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [0, 1, 3, 6, 6, 6, 6, 6, 6, 6]
        """
        return sum([c[n] for c in self.children()], self._zero)

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream, SeriesStreamFromIterator
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: s.is_constant()
            True
            sage: z = SeriesStreamFromIterator(iterator=ZZ, base_ring=ZZ)
            sage: s = ListSumStream(t + [z], base_ring=ZZ)
            sage: s.is_constant()
            False
        """
        return all(c.is_constant() for c in self.children())

    def get_constant(self):
        """
        If the stream is eventually constant, returns the constant
        value.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: s.get_constant()
            6
        
        """
        return sum(c.get_constant() for c in self.children())

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ListSumStream
            sage: t = [SeriesStreamFromList(list=[0]*i + [i], base_ring=ZZ) for i in range(4)]
            sage: s = ListSumStream(t, base_ring=ZZ)
            sage: s.get_constant_position()
            3    
        """
        return max(c.get_constant_position() for c in self.children())

class SumGeneratorStream(SeriesStream):
    def __init__(self, series_stream, **kwds):
        """
        A class for a stream whose $n^{th}$ coefficient is given by
        the sum of $n^{th}$ coefficients for the first $n$ series in
        ``series_stream``.

        INPUT:

        - ``series_stream`` - a :class:`Stream` whose values are
          themselves series or instances of :class:`SeriesStream`.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, StreamFromList, SumGeneratorStream
            sage: one = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: zero = SeriesStreamFromList(list=[0], base_ring=ZZ)
            sage: t = StreamFromList(list=[one]*6 + [zero])
            sage: s = SumGeneratorStream(t, base_ring=ZZ)
            sage: [s[i] for i in range(10)]
            [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]
        """
        self._series_stream = SeriesStreamFromIterator(iterator=iter(series_stream), **kwds)
        super(SumGeneratorStream, self).__init__(**kwds)

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, StreamFromList, SumGeneratorStream
            sage: one = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: zero = SeriesStreamFromList(list=[0], base_ring=ZZ)
            sage: t = StreamFromList(list=[one]*6 + [zero])
            sage: s = SumGeneratorStream(t, base_ring=ZZ)
            sage: [s.compute(i) for i in range(10)]
            [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]
        """
        s = self._series_stream
        r = s[n][n]
        for i in range(min(n, s.number_computed() - 1)):
            r += s[i][n]
        return r

class ProductGeneratorStream(SeriesStreamFromIterator):
    def __init__(self, series_iter, **kwds):
        r"""
        A class for streams representing the product of a potentially
        infinite number of power series.  Let
        $g_0, g_1, \ldots, g_n, \ldots$ represent the power series given by
        ``series_iter``.  In order to do so, we place restrictions
        on the input series.  We require that

        .. math::

           g_n = 1 + \sum_{i=n}^{\infty} a_{n,i} x^i.

        With that restriction, we can compute the $n^{th}$ coefficient
        by multiplying the first $n$ series.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ProductGeneratorStream
            sage: g = [SeriesStreamFromList(list=[1] + [0]*i + [1, 0], base_ring=ZZ) for i in range(6)]
            sage: s = ProductGeneratorStream(g, base_ring=ZZ)
            sage: [s[i] for i in range(26)]
            [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 4, 4, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0]
        """
        self._series_it = iter(series_iter)
        super(ProductGeneratorStream, self).__init__(iterator=self.compute_iterator(), **kwds)

    def compute_iterator(self):
        """
        A generator for the coefficients of this potentially infinite
        product.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList, ProductGeneratorStream
            sage: g = [SeriesStreamFromList(list=[1] + [0]*i + [1, 0], base_ring=ZZ) for i in range(6)]
            sage: s = ProductGeneratorStream(g, base_ring=ZZ)
            sage: it = s.compute_iterator()
            sage: [it.next() for i in range(26)]
            [1, 1, 1, 2, 2, 3, 4, 4, 4, 5, 5, 5, 5, 4, 4, 4, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0]
        """
        z = self._series_it.next()
        yield z[0]
        yield z[1]

        n = 2
        for x in self._series_it:
            z = z * x
            yield z[n]
            n += 1

        # In case there are only finitely many series in the iterator
        while True:
            yield z[n]
            n += 1

class PowerStream(StreamFromIterator):
    def __init__(self, stream, **kwds):
        """
        A stream for the powers of a series stream ``stream`` starting
        with ``stream^1``.

        .. note::

           This is not a :class:`SeriesStream`.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import PowerStream, SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = PowerStream(t)
            sage: [s[0][i] for i in range(6)] # t^1
            [1, 1, 1, 1, 1, 1]
            sage: [s[1][i] for i in range(6)] # t^2
            [1, 2, 3, 4, 5, 6]
            sage: [s[2][i] for i in range(6)] # t^3
            [1, 3, 6, 10, 15, 21]
        """
        self._stream = stream
        super(PowerStream, self).__init__(iterator=self.compute_iterator(), **kwds)

    def compute_iterator(self):
        """
        A generator for the powers of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import PowerStream, SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1], base_ring=ZZ)
            sage: s = PowerStream(t)
            sage: g = s.compute_iterator()
            sage: t1 = g.next()
            sage: [t1[i] for i in range(6)] # t^1
            [1, 1, 1, 1, 1, 1]
            sage: t2 = g.next()
            sage: [t2[i] for i in range(6)] # t^2
            [1, 2, 3, 4, 5, 6]
            sage: t3 = g.next()
            sage: [t3[i] for i in range(6)] # t^3
            [1, 3, 6, 10, 15, 21]
        """
        z = self._stream
        while True:
            yield z
            z = ProductStream(z, self._stream, base_ring=z._base_ring)

class StretchedStream(SeriesStream):
    def __init__(self, k, stream, **kwds):
        """
        A class for a stretched series stream.  The $n^{th}$
        coefficient of a stretched stream is $0$ if $n$ is not
        divisible by $k$; otherwise, it is the $(n/k)^{th}$
        coefficient of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1, 2, 3], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: [s[i] for i in range(10)]
            [1, 0, 2, 0, 3, 0, 3, 0, 3, 0]
        """
        self._k = k
        self._stream = stream
        super(StretchedStream, self).__init__(children=[stream], **kwds)

    def order_operation(self, a):
        """
        Returns the approximate order of a stretched stream, which is
        $k$ times the approximate order of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[0,1,1], base_ring=ZZ)
            sage: s = t.stretch(3)
            sage: s.order_operation(*[c.aorder for c in s.children()])
            3
        """
        return self._k * a

    def compute(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream, which is $0$
        if $n$ is not divisible by $k$; otherwise, it returns the
        $(n/k)^{th}$ coefficient of the underlying stream.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1, 2, 3], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: [s.compute(i) for i in range(10)]
            [1, 0, 2, 0, 3, 0, 3, 0, 3, 0]
        """
        n = Integer(n)
        quo, rem = n.quo_rem(self._k)
        if rem == 0:
            return self._stream[quo]
        else:
            return self._zero

    def is_constant(self):
        """
        Returns whether or not this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1, 2, 0], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: s.is_constant()
            True
            sage: t = SeriesStreamFromList(list=[1, 2, 3], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: s.is_constant()
            False        
        """
        return (self._stream.is_constant() and
                self._stream.get_constant() == 0)

    def get_constant(self):
        """
        If this stream is eventually constant, returns the constant
        value.  For stretched streams, this always has to be zero.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1, 2, 0], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: s.get_constant()
            0
        """
        assert self.is_constant()
        return self._zero

    def get_constant_position(self):
        """
        If the stream is eventually constant, returns the position
        where the stream becomes constant.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromList
            sage: t = SeriesStreamFromList(list=[1, 2, 0], base_ring=ZZ)
            sage: s = t.stretch(2)
            sage: s.get_constant_position()
            3
            sage: [s[i] for i in range(10)]
            [1, 0, 2, 0, 0, 0, 0, 0, 0, 0]
        """
        p = self._stream.get_constant_position()
        return (p - 1) * self._k + 1
