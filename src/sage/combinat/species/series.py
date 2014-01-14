"""
Lazy Power Series

This file provides an implementation of lazy univariate power
series, which uses the stream class for its internal data
structure.

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
from series_order import  inf, unk
from sage.rings.all import Integer
from sage.misc.misc import repr_lincomb, is_iterator
from sage.misc.cachefunc import cached_method
import sage.structure.parent_base
from sage.categories.all import Rings
from series_stream import (SeriesStream, SeriesStreamFromList, SeriesStreamFromIterator, TermStream,
                           ChangeRingStream, SumStream, ProductStream, TailStream, DerivativeStream,
                           IntegralStream, CompositionStream, RecursiveStream, RestrictedStream,
                           ListSumStream, SumGeneratorStream, ProductGeneratorStream)
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element, AlgebraElement
from sage.algebras.algebra import Algebra
from sage.categories.algebras import Algebras

class LazyPowerSeriesRing(UniqueRepresentation, Algebra):
    def __init__(self, R, element_class=None, category=None, names=None):
        """
        TESTS::

            sage: from sage.combinat.species.series import LazyPowerSeriesRing
            sage: L = LazyPowerSeriesRing(QQ)
            sage: TestSuite(L).run()
        """
        #Make sure R is a ring with unit element
        if not R in Rings():
            raise TypeError, "Argument R must be a ring."
        try:
            R(Integer(1))
        except StandardError:
            raise ValueError, "R must have a unit element"

        if category is None:
            category = Algebras(R)

        #Take care of the names
        if names is None:
            names = 'x'
        else:
            names = names[0]

        self.Element = element_class if element_class is not None else LazyPowerSeries
        self._name = names

        Parent.__init__(self, base=R, category=category)

    def ngens(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ).ngens()
            1
        """
        return 1

    def _repr_(self):
        """
        EXAMPLES::

            sage: LazyPowerSeriesRing(QQ)
            Lazy Power Series Ring over Rational Field
        """
        return "Lazy Power Series Ring over %s"%self.base_ring()

    def _new(self, stream_cls, *args, **kwds):
        """
        Returns a new lazy power series with class ``stream_cls``.
        This function automatically sets the ``base_ring`` parameter.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream
            sage: L1 = LazyPowerSeriesRing(QQ)
            sage: g = L1._new(TermStream, n=2, value=3)
            sage: g.coefficients(5)
            [0, 0, 3, 0, 0]

        """
        kwds['base_ring'] = self.base_ring()
        return self.element_class(self, stream_cls(*args, **kwds))

    def _coerce_map_from_(self, ring):
        """
        EXAMPLES::

            sage: L1 = LazyPowerSeriesRing(QQ)
            sage: L2 = LazyPowerSeriesRing(RR)
            sage: L1.has_coerce_map_from(ZZ)
            True
            sage: L2.has_coerce_map_from(L1)
            True
            sage: L1.has_coerce_map_from(L2)
            False

        ::

            sage: a = L1([1]) + L2([1])
            sage: a.coefficients(3)
            [2.00000000000000, 2.00000000000000, 2.00000000000000]
        """
        if self.base_ring().has_coerce_map_from(ring):
            return True
        if (self.__class__.__base__ is ring.__class__.__base__ and
            self.base_ring().has_coerce_map_from(ring.base_ring())):
            return True

    def __call__(self, x=None, *args, **kwds):
        """
        We override the default call method of ``Parent`` to make the
        default value ``None`` instead of ``x``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L()
            O(1)
            sage: L(0)
            0        
        """
        return super(LazyPowerSeriesRing, self).__call__(x=x, *args, **kwds)
        
    def _element_constructor_(self, x=None, order=unk):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L()
            O(1)
            sage: L(1)
            1
            sage: L(ZZ).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: L(iter(ZZ)).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]

        ::

            sage: from sage.combinat.species.series_stream import SeriesStreamFromIterator
            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = SeriesStreamFromIterator(iterator=ZZ, base_ring=QQ, convert=True)
            sage: L(s).coefficients(10)
            [0, 1, -1, 2, -2, 3, -3, 4, -4, 5]
            sage: _[0].parent()
            Rational Field

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
        cls = self.element_class
        base_ring = self.base_ring()

        if x is None:
            return self._new(RecursiveStream)

        if hasattr(x, "parent") and base_ring.has_coerce_map_from(x.parent()):
            x = base_ring(x)
            return self.term(x, 0)

        if isinstance(x, LazyPowerSeries):
            x_parent = x.parent()
            if x_parent.__class__.__base__ != self.__class__.__base__:
                raise ValueError

            if x_parent.base_ring() == base_ring:
                return x
            else:
                if base_ring.has_coerce_map_from(x_parent.base_ring()):
                    return self._new(ChangeRingStream, stream=x._stream,
                                     new_ring=base_ring)


        if isinstance(x, SeriesStream):
            assert x.base_ring() == base_ring
            return cls(self, x)
        
        if isinstance(x, (list, tuple)):
            x = SeriesStreamFromList(list=x, base_ring=base_ring)
        elif hasattr(x, "__iter__"):
            x = iter(x)

        if is_iterator(x):
            x = SeriesStreamFromIterator(iterator=x, base_ring=base_ring, convert=True)

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

    @cached_method
    def gen(self, i=0):
        """
        Returns the generator for this algebra.
        
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L.gen().coefficients(5)
            [0, 1, 0, 0, 0]
        """
        assert i == 0
        res = self._new(TermStream, 1, 1)
        res.rename(self._name)
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
        if r == 0:
            res = self._new(TermStream, 0, 0)
            res.rename("0")
        else:
            res = self._new(TermStream, n, r)

            if n == 0:
                res.rename(repr(r))
            elif n == 1:
                res.rename(repr(r) + "*" + self._name)
            else:
                res.rename("%s*%s^%s"%(repr(r), self._name, n))

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

    def sum_generator(self, g):
        """
        Returns a lazy power series whose $n^{th}$ coefficients is
        the sum of $n^{th}$ coefficients for the first $n$ series in
        ``g``.

        INPUT:

        - ``g`` - an iterable of lazy power series
        
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
        r"""
        Returns a lazy power series which is the product of the
        (potentially infinite) series in the iterable ``g``.  In order
        to do this, we place restrictions on the form of the series
        $g_n$.  In particular we require that

        .. math::

           g_n = 1 + \sum_{i=n}^{\infty} a_{n,i} x^i.

        Then, the $n^{th}$ coefficient can be obtained from the first
        $n$ series.

        INPUT:

        - ``g`` - an iterable of lazy power series
        
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
            sage: f = L(2)
            sage: TestSuite(f).run()
            sage: f = L(); f
            O(1)
        """
        super(LazyPowerSeries, self).__init__(A)
        self._stream = stream
        self._zero = A.base_ring().zero_element()

    def __nonzero__(self):
        """
        Returns ``True`` if this series is nonzero; otherwise, it
        returns ``False``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: bool(L(2))
            True
            sage: bool(L(0))
            False
        """
        return self != self.parent().zero_element()

    def __ne__(self, other):
        """
        Returns ``True`` if ``self != other``.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: L(2) != 2.3
            True
            sage: L(3) != L(3)
            False
            sage: L([1,2,3]) != L([2,3,4])
            True
            sage: L([1,2,3]) != L([1,2,3])
            False
        """
        return not (self == other)

    def __eq__(self, other):
        """
        Returns ``True`` if ``self == other``.

        EXAMPLES::

            sage: L.<x> = LazyPowerSeriesRing(QQ)
            sage: x^2 == L([0,0,1,0])
            True
            sage: L(3) == 2.3
            False
        """
        if type(self) != type(other):
            return False
        if not (self._stream.is_constant() and other._stream.is_constant()):
            return False
        for i in range(max(self._stream.get_constant_position(),
                           other._stream.get_constant_position())):
            if self[i] != other[i]:
                return False
        else:
            return True
        
    @property
    def aorder(self):
        """
        Returns the (currently computed) approximate order for this
        series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = L(iter([0,3,0]))
            sage: g.aorder
            Unknown series order
            sage: g[1]
            3
            sage: g.aorder
            0
            sage: g.get_aorder()
            1
            sage: g.aorder
            1
        """
        return self._stream.aorder

    def get_aorder(self):
        """
        Returns the approximate order for this series after refining
        it (without computing additional coefficients).

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = L(iter([0,3,0]))
            sage: g.aorder
            Unknown series order
            sage: g[1]
            3
            sage: g.get_aorder()
            1
        """
        return self._stream.get_aorder()

    @property
    def order(self):
        """
        Returns the (currently computed) approximate order for this
        series.

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = L(iter([0,0,3,0]))
            sage: g.order
            Unknown series order
            sage: g[1]
            0
            sage: g.order
            Unknown series order
            sage: g.get_order()
            Unknown series order
            sage: g[2]
            3
            sage: g.order
            Unknown series order
            sage: g.get_order()
            2
            sage: g.order
            2
        """
        return self._stream.order

    def get_order(self):
        """
        Returns the order for this series after refining
        it (without computing additional coefficients).

        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: g = L(iter([0,3,0]))
            sage: g.order
            Unknown series order
            sage: g[1]
            3
            sage: g.get_order()
            1
        """
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
        return [(_m, _c) for _m, _c in zip(m, c) if c != 0]

    def _repr_(self):
        """
        EXAMPLES::

            sage: L = LazyPowerSeriesRing(QQ)
            sage: s = L(); s.rename('s'); s
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
        """
        Returns a new lazy power series with class ``stream_cls``.
        This function automatically sets the ``base_ring`` parameter.

        EXAMPLES::

            sage: from sage.combinat.species.series_stream import TermStream
            sage: L1 = LazyPowerSeriesRing(QQ)
            sage: x = L1.gen()
            sage: g = x._new(TermStream, n=2, value=3)
            sage: g.coefficients(5)
            [0, 0, 3, 0, 0]

        """
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
        
        return (self._stream.is_constant() and
                self._stream.get_constant() == 0)

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
