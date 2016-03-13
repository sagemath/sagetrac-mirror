r"""
Lazy lists

A lazy list is an iterator that behaves like a list and possesses a cache
mechanism. A lazy list is potentially infinite and speed performances of the
cache is comparable with Python lists. One major difference with original
Python list is that lazy list are immutable. The advantage is that slices
share memory.

EXAMPLES::

    sage: from sage.misc.lazy_list import lazy_list
    sage: P = lazy_list(Primes())
    sage: P[100]
    547
    sage: P[10:34]
    lazy list [31, 37, 41, ...]
    sage: P[12:23].list()
    [41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83]

    sage: f = lazy_list((i**2-3*i for i in xrange(10)))
    sage: for i in f: print i,
    0 -2 -2 0 4 10 18 28 40 54
    sage: i1 = iter(f)
    sage: i2 = iter(f)
    sage: print next(i1), next(i1)
    0 -2
    sage: print next(i2), next(i2)
    0 -2
    sage: print next(i1), next(i2)
    -2 -2

It is possible to prepend a list to a lazy list::

    sage: from itertools import count
    sage: l = [3,7] + lazy_list(i**2 for i in count())
    sage: l
    lazy list [3, 7, 0, ...]

But, naturally, not the other way around::

    sage: lazy_list(i-1 for i in count()) + [3,2,5]
    Traceback (most recent call last):
    ...
    TypeError: can only add list to lazy_list

You can easily create your own class inheriting from :class:`lazy_list_generic`. You
should call the :class:`lazy_list_generic` constructor (optionally with some
precomputed values for the cache) and implement the method ``_new_slice`` that
returns a new chunk of data at each call. Here is an example of implementation
of the Thue--Morse word that is obtained as the fixed point of the substitution
`0 \to 01` and `1 \to 10`::

    sage: from sage.misc.lazy_list import lazy_list_generic
    sage: class MyThueMorseWord(lazy_list_generic):
    ....:     def __init__(self):
    ....:         self.i = 1
    ....:         lazy_list_generic.__init__(self, cache=[0,1])
    ....:     def _new_slice(self):
    ....:         letter = self.get(self.i)
    ....:         self.i += 1
    ....:         return [0,1] if letter == 0 else [1,0]
    sage: w = MyThueMorseWord()
    sage: w
    lazy list [0, 1, 1, ...]
    sage: all(w[i] == ZZ(i).popcount()%2 for i in range(100))
    True
    sage: w[:500].list() == w[:1000:2].list()
    True

Alternatively, you can create the lazy list from an update function::

    sage: def thue_morse_update(values):
    ....:     n = len(values)
    ....:     if n == 0:
    ....:         letter = 0
    ....:     else:
    ....:         assert n%2 == 0
    ....:         letter = values[n//2]
    ....:     values.append(letter)
    ....:     values.append(1-letter)
    sage: w2 = lazy_list(update_function=thue_morse_update)
    sage: w2
    lazy list [0, 1, 1, ...]
    sage: w2[:500].list() == w[:500].list()
    True

You can also create extension type inheriting from :class:`lazy_list_generic`
(with Cython). In that case you would better implement directly the method
:meth:`~lazy_list_generic.update_cache_up_to`. See the examples in this file with the classes
:class:`lazy_list_from_iterator` and :class:`lazy_list_from_function`.

Classes and Methods
===================
"""
#*****************************************************************************
#       Copyright (C) 2015 Vincent Delecroix <20100.delecroix@gmail.com>
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from "Python.h":
    Py_ssize_t PY_SSIZE_T_MAX

# make a unique instance of empty lazy lists
cdef lazy_list_generic empty_lazy_list
empty_lazy_list = lazy_list_generic.__new__(lazy_list_generic)
empty_lazy_list.start = 0
empty_lazy_list.stop = 0
empty_lazy_list.step = 1
empty_lazy_list.cache = []
empty_lazy_list.name = 'lazy list'
empty_lazy_list.opening_delimiter = '['
empty_lazy_list.separator = ', '
empty_lazy_list.more = '...'
empty_lazy_list.closing_delimiter = ']'
empty_lazy_list.preview = 3
empty_lazy_list.cls = lazy_list_generic
empty_lazy_list.cls_kwds = {}
#empty_lazy_list = lazy_list_generic(initial_values=[],
#                                    start=0, stop=0, step=1)  # ... does not work


def lazy_list(data=None, initial_values=None, start=None, stop=None, step=None,
              update_function=None, **kwds):
    r"""
    Return a lazy list.

    INPUT:

    - ``data`` -- data to create a lazy list from. This can be

      #. a (possibly infinite) iterable,
      #. a function (that takes as input an integer ``n`` and return
         the ``n``-th term of the list),
      #. or a standard Python container ``list`` or ``tuple``.

    - ``initial_values`` -- the beginning of the sequence that will not be computed from
      the ``data`` provided.

    - ``update_function`` -- you can also construct a lazy list from a function
      that takes as input a list of precomputed values and updates it with some
      more values.

    Representation string options:

    - ``name`` -- (default: ``'lazy list'``) a string.

    - ``opening_delimiter`` -- (default: ``'['``) a string.

    - ``closing_delimiter`` -- (default: ``']'``) a string.

    - ``separator`` -- (default: ``', '``) a string.

    - ``more`` -- (default: ``'...'``) a string.

    - ``preview`` -- (default: ``3``) an integer specifying the number of
      elements shown in the representation string.

    Other options:

    - ``cls`` -- (default: ``None``) a class. If specified, then when
      creating sublists (like slices, or
      :meth:`~lazy_list_generic.dropwhile` or
      :meth:`~lazy_list_generic.takewhile`) the result is assured to be an
      instance of ``cls``.

    - ``cls_kwds`` -- (default: ``None``) a dictionary with keyword
      arguments which are passed to an instantiating of ``cls``.

    - ``start``, ``stop``, ``step`` -- deprecated arguments

    .. NOTE::

        If you want finer tuning of the constructor you can directly instantiate
        the classes associated to lazy lists that are
        :class:`lazy_list_generic`, :class:`lazy_list_from_iterator`,
        :class:`lazy_list_from_function`.

    EXAMPLES:

    The basic construction of lazy lists.
    ::

        sage: from sage.misc.lazy_list import lazy_list

    #. Iterators::

        sage: from itertools import count
        sage: lazy_list(count())
        lazy list [0, 1, 2, ...]

    #. Functions::

        sage: lazy_list(lambda n: (n**2)%17)
        lazy list [0, 1, 4, ...]

    #. Plain lists::

        sage: lazy_list([1,5,7,2])
        lazy list [1, 5, 7, ...]

    If a function is only defined for large values, you can provide the beginning
    of the sequence manually::

        sage: l = lazy_list(divisors, [None])
        sage: l
        lazy list [None, [1], [1, 2], ...]

    Lazy lists behave like lists except that they are immutable::

        sage: l[3::5]
        lazy list [[1, 3], [1, 2, 4, 8], [1, 13], ...]

    If your lazy list is finite, you can obtain the underlying list with the
    method `.list()`::

        sage: l[30:50:5].list()
        [[1, 2, 3, 5, 6, 10, 15, 30],
         [1, 5, 7, 35],
         [1, 2, 4, 5, 8, 10, 20, 40],
         [1, 3, 5, 9, 15, 45]]

    TESTS::

        sage: lazy_list(count(), start=5)
        doctest:...: DeprecationWarning: The arguments start, stop, step are deprecated. Use
        direct slicing as in my_data[start:stop:step]
        See http://trac.sagemath.org/16137 for details.
        lazy list [5, 6, 7, ...]

        sage: lazy_list()
        lazy list []
        sage: lazy_list(data='hey', update_function='hello')
        Traceback (most recent call last):
        ...
        ValueError: only one of the arguments 'data' or 'update_function'
        can be used

    Lazy lists created from each other share their cache::

        sage: C = lazy_list(count())
        sage: C[4]
        4
        sage: D = lazy_list(C)
        sage: D._info()
        cache length 5
        start        0
        stop         9223372036854775807
        step         1
        sage: D[6]
        6
        sage: C._info()
        cache length 7
        start        0
        stop         9223372036854775807
        step         1

    ::

        sage: from sage.misc.lazy_list import lazy_list_generic
        sage: class Z(lazy_list_generic):
        ....:     pass
        sage: P = lazy_list(Primes(), cls=Z, name='primes'); P
        primes [2, 3, 5, ...]
        sage: type(P)
        <class '__main__.Z'>
    """
    cdef lazy_list_generic l

    if data is None and update_function is None:
        return empty_lazy_list
    elif data is not None and update_function is not None:
        raise ValueError("only one of the arguments 'data' or 'update_function' can be used")

    if initial_values is None:
        cache = []
    else:
        cache = list(initial_values)

    if update_function is not None:
        assert callable(update_function)
        return lazy_list_from_update_function(update_function, cache)

    if isinstance(data, (tuple,list)):
        data = cache + list(data)
        l = lazy_list_generic(data, start=0, stop=len(data), step=1, **kwds)
    elif isinstance(data, lazy_list_generic):
        l = data.make_linked_copy(**kwds)  # to share cache
    else:
        # the code below is not very clean
        # we just want to differentiate on the one hand iterable (= object with a
        # proper __iter__ method)/iterator (= object with a next method) and on the
        # other hand callable (= object with __call__)
        try:
            data = iter(data)
        except TypeError:
            pass

        from sage.misc.misc import is_iterator
        if is_iterator(data):
            l = lazy_list_from_iterator(iter(data), cache, **kwds)
        elif callable(data):
            l = lazy_list_from_function(data, cache, **kwds)
        else:
            raise ValueError("not able to build a lazy list from {}".format(type(data)))

    l = l.change_class()

    if start is not None or stop is not None or step is not None:
        from sage.misc.superseded import deprecation
        deprecation(16137, "The arguments start, stop, step are deprecated. "
                           "Use direct slicing as in my_data[start:stop:step]")
        return l[start:stop:step]
    else:
        return l


def slice_unpickle(master, start, stop, step):
    r"""
    Unpickle helper

    TESTS::

        sage: from sage.misc.lazy_list import slice_unpickle
        sage: slice_unpickle(range(35), 1, 3, 7) == range(35)[1:3:7]
        True
    """
    return master[start:stop:step]


cdef class lazy_list_generic(object):
    def __init__(self, cache=None, master=None,
                 start=None, stop=None, step=None,
                 name=None, separator=None, more=None,
                 opening_delimiter=None, closing_delimiter=None,
                 preview=None,
                 cls=None, cls_kwds=None):
        r"""
        A lazy list (generic base class).

        INPUT:

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy.

        Representation string options:

        - ``name`` -- (default: ``'lazy list'``) a string.

        - ``opening_delimiter`` -- (default: ``'['``) a string.

        - ``closing_delimiter`` -- (default: ``']'``) a string.

        - ``separator`` -- (default: ``', '``) a string.

        - ``more`` -- (default: ``'...'``) a string.

        - ``preview`` -- (default: ``3``) an integer specifying the number of
          elements shown in the representation string.

        Other options:

        - ``cls`` -- (default: ``None``) a class. If specified, then when
          creating sublists (like slices, or
          :meth:`~lazy_list_generic.dropwhile` or
          :meth:`~lazy_list_generic.takewhile`) the result is assured to be an
          instance of ``cls``.

        - ``cls_kwds`` -- (default: ``None``) a dictionary with keyword
          arguments which are passed to an instantiating of ``cls``.

        - ``start``, ``stop``, ``step`` -- for slices

        .. NOTE::

            Everywhere the constant ``PY_SSIZE_T_MAX`` plays the role of infinity

        .. WARNING::

            No check is performed on input and bad input can result in a Sage crash.
            You are advised to use the function :func:`lazy_list` instead. The only
            case where you might want to use directly this constructor is if you
            have a list that you want to wrap (without copy) into a lazy list.
            See also the example below.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: l = lazy_list(Primes())
            sage: l
            lazy list [2, 3, 5, ...]
            sage: l[200]
            1229

        ::

            sage: from sage.misc.lazy_list import lazy_list_generic
            sage: l = [0, 1, 2]
            sage: ll = lazy_list_generic(l, start=0, stop=2, step=None)
            sage: ll
            lazy list [0, 1]

        The above code may be dangerous since the lazy list holds a reference
        to the initial list::

            sage: l[0] = 'haha'
            sage: ll
            lazy list ['haha', 1]
        """
        self.cache = [] if cache is None else cache
        self.master = master
        self._start_master = 0 if self.master is None else self.master.start

        self.start = 0 if start is None else start
        self.stop = PY_SSIZE_T_MAX if stop is None else stop
        self.step = 1 if step is None else step

        self.name = 'lazy list' if name is None else name
        self.opening_delimiter = '[' if opening_delimiter is None else opening_delimiter
        self.separator = ', ' if separator is None else separator
        self.more = '...' if more is None else more
        self.closing_delimiter = ']' if closing_delimiter is None else closing_delimiter
        self.preview = 3 if preview is None else preview

        self.cls = lazy_list_generic if cls is None else cls
        self.cls_kwds = {} if cls_kwds is None else cls_kwds

    def start_stop_step(self):
        r"""
        Return the triple ``(start, stop, step)`` of reference points of the
        original lazy list.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: p = lazy_list(Primes())[100:1042240:12]
            sage: p.start_stop_step()
            doctest:...: DeprecationWarning: The method start_stop_step is deprecated. Consider using _info() instead.
            See http://trac.sagemath.org/16137 for details.
            (100, 1042240, 12)
        """
        from sage.misc.superseded import deprecation
        deprecation(16137, "The method start_stop_step is deprecated. Consider using _info() instead.")
        return (self.start, self.stop, self.step)


    def list(self):
        r"""
        Return the list made of the elements of ``self``.

        .. NOTE::

            If the iterator is sufficiently large, this will build a list
            of length ``PY_SSIZE_T_MAX`` which should be beyond the capacity of
            your RAM!

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: P = lazy_list(Primes())
            sage: P[2:143:5].list()
            [5, 19, 41, 61, 83, 107, 137, 163, 191, 223, 241, 271, 307, 337, 367, 397, 431, 457, 487, 521, 563, 593, 617, 647, 677, 719, 751, 787, 823]
            sage: P = lazy_list(iter([1,2,3]))
            sage: P.list()
            [1, 2, 3]
            sage: P[:100000].list()
            [1, 2, 3]
            sage: P[1:7:2].list()
            [2]

        TESTS:

        Check that the cache is immutable::

            sage: lazy = lazy_list(iter(Primes()))[:5]
            sage: l = lazy.list(); l
            [2, 3, 5, 7, 11]
            sage: l[0] = -1; l
            [-1, 3, 5, 7, 11]
            sage: lazy.list()
            [2, 3, 5, 7, 11]
        """
        self._fit(self.stop - self.step)
        return self.cache[self.start:self.stop:self.step]


    def info(self):
        r"""
        Deprecated method

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: lazy_list([0]).info()
            doctest:...: DeprecationWarning: info is deprecated in favor of a private method.
            Use _info() instead
            See http://trac.sagemath.org/19428 for details.
            cache length 1
            start        0
            stop         1
            step         1
        """
        from sage.misc.superseded import deprecation
        deprecation(19428, "info is deprecated in favor of a private method. Use _info() instead")
        return self._info()


    def _info(self):
        r"""
        Print information about ``self`` on standard output.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: P = lazy_list(iter(Primes()))[10:21474838:4]
            sage: P._info()
            cache length 0
            start        10
            stop         21474838
            step         4
            sage: P[0]
            31
            sage: P._info()
            cache length 11
            start        10
            stop         21474838
            step         4
        """
        print "cache length", len(self.cache)
        print "start       ", self.start
        print "stop        ", self.stop
        print "step        ", self.step


    def __add__(self, other):
        r"""
        If ``self`` is a list then return the lazy_list that consists of the
        concatenation of ``self`` and ``other``.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import count
            sage: l = lazy_list(i**3 - i + 1 for i in count()); l
            lazy list [1, 1, 7, ...]
            sage: p = ['huit', 'douze']
            sage: ll = p + l; ll
            lazy list ['huit', 'douze', 1, ...]
            sage: l[:10].list() == ll[2:12].list()
            True
            sage: p
            ['huit', 'douze']
            sage: ([0,2] + lazy_list([0,1])).list()
            [0, 2, 0, 1]
        """
        if not isinstance(self, list):
            raise TypeError("can only add list to lazy_list")

        return lazy_list_from_iterator(iter(other), cache=self[:])


    def __repr__(self):
        r"""
        Return a string representation.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import count
            sage: r = lazy_list(count()); r  # indirect doctest
            lazy list [0, 1, 2, ...]
            sage: r[:0]
            lazy list []
            sage: r[:1]
            lazy list [0]
            sage: r[:2]
            lazy list [0, 1]
            sage: r[:3]
            lazy list [0, 1, 2]
            sage: r[:4]
            lazy list [0, 1, 2, ...]
            sage: lazy_list([0,1])
            lazy list [0, 1]
            sage: lazy_list([0,1,2,3])
            lazy list [0, 1, 2, ...]
        """
        cdef str s = self.name
        if s:
            s += ' '
        s += self.opening_delimiter
        cdef list P = list(self[:self.preview+1])
        cdef list E = list('{!r}'.format(e)
                           for e in P[:self.preview])
        cdef Py_ssize_t num_elts = 1 + (self.stop-self.start-1) / self.step

        if num_elts > self.preview:
            E.append(self.more)
        s += self.separator.join(E)
        s += self.closing_delimiter
        return s


    def __reduce__(self):
        r"""
        Pickling support

        EXAMPLES::

            sage: from itertools import count
            sage: from sage.misc.lazy_list import lazy_list
            sage: m = lazy_list(count())
            sage: x = loads(dumps(m))
            sage: y = iter(x)
            sage: print next(y), next(y), next(y)
            0 1 2
            sage: m2 = m[3::2]
            sage: loads(dumps(m2))
            lazy list [3, 5, 7, ...]
        """
        if self.master is None:
            raise NotImplementedError
        return slice_unpickle, (self.master, self.start, self.stop, self.step)


    cpdef int _fit(self, Py_ssize_t n) except -1:
        r"""
        Fill the cache making the term at index ``n`` available.

        You can access the term at position ``n`` from the cache when it returns
        ``0``.

        OUTPUT:

        - ``1`` -- the lazy list is actually finite and shorter than ``n``

        - ``0`` -- you can safely access the term at position ``n`` after this call

        - ``-1`` -- to handle Python errors (you can ignore it in Python code)

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: l = lazy_list(iter([0,1,2,-34,3,2,-5,12,1,4,-18,5,-12]))[2::3]
            sage: l._info()
            cache length 0
            start        2
            stop         9223372036854775807    # 64-bit
            stop         2147483647             # 32-bit
            step         3
            sage: l._fit(13)
            1
            sage: l._info()
            cache length 13
            start        2
            stop         14
            step         3

            sage: l = lazy_list([0]*12)[1::2]
            sage: l._fit(100)
            1
            sage: l._info()
            cache length 12
            start        1
            stop         13
            step         2
            sage: l._fit(100)
            1
        """
        if n > self.stop - self.step:
            return 1

        cdef int result = self.update_cache_up_to(n)

        if self.master is not None:
            # We have to update start in case the tracked master list
            # changed its start (e.g. dropwhile).
            if self.master.start != self._start_master:
                self.start += self.master.start - self._start_master
                s = PY_SSIZE_T_MAX - self.stop
                s = min(s, self.master.start - self._start_master)
                self.stop += s
                self._start_master = self.master.start

        if result:
            self.stop = min(self.stop, len(self.cache))
            if self.master is not None:
                self.stop = min(self.stop, self.master.stop)
            if self.stop <= self.start:
                self.start = self.stop = 0
                self.step = 1
            if (self.start - self.stop) % self.step:
                self.stop += self.step + (self.start - self.stop) % self.step
            return 1

        return 0


    cpdef get(self, Py_ssize_t i):
        r"""
        Return the element at position ``i``.

        If the index is not an integer, then raise a ``TypeError``.  If the
        argument is negative then raise a ``ValueError``.  Finally, if the
        argument is beyond the size of that lazy list it raises a
        ``IndexError``.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import chain, repeat
            sage: f = lazy_list(chain(iter([1,2,3]), repeat('a')))
            sage: f.get(0)
            1
            sage: f.get(3)
            'a'
            sage: f.get(0)
            1
            sage: f.get(4)
            'a'

            sage: g = f[:10]
            sage: g.get(5)
            'a'
            sage: g.get(10)
            Traceback (most recent call last):
            ...
            IndexError: lazy list index out of range
            sage: g.get(1/2)
            Traceback (most recent call last):
            ...
            TypeError: rational is not an integer

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_generic, lazy_list_from_iterator
            sage: class Z(lazy_list_generic):
            ....:     pass
            sage: lazy_list(Primes(), cls=Z).dropwhile(lambda x: x <= 5)[1]
            11
        """
        if i < 0:
            raise ValueError("indices must be non-negative")

        n = self.start + i*self.step
        if self._fit(n):
            raise IndexError("lazy list index out of range")
        n = self.start + i*self.step  # update since _fit may have
                                      # changed self.start
        return self.cache[n]


    def __call__(self, i):
        r"""
        An alias for :meth:`get`

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import chain, repeat
            sage: f = lazy_list(chain(iter([1,2,3]), repeat('a')))
            sage: f(2)
            3
            sage: f(3)
            'a'
        """
        return self.get(i)


    def __iter__(self):
        r"""
        Return an iterator.

        TESTS::

            sage: from itertools import count
            sage: from sage.misc.lazy_list import lazy_list
            sage: iter(lazy_list(count()))
            <generator object at 0x...>

        ::

            sage: l = lazy_list(i^2 for i in xrange(5))
            sage: list(l)
            [0, 1, 4, 9, 16]
            sage: l._info()
            cache length 5
            start        0
            stop         5
            step         1

        ::

            sage: tuple(lazy_list(Primes()).dropwhile(lambda x: x < 10)[:5])
            (11, 13, 17, 19, 23)
        """
        cdef Py_ssize_t i

        if self._fit(self.start):
            return
        i = self.start
        while i < self.stop:
            if self._fit(i):
                return
            yield self.cache[i]
            i += self.step


    def __getitem__(self, key):
        r"""
        Return a lazy list which shares the same cache.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: f = lazy_list(iter([1,2,3]))
            sage: f0 = f[0:]
            sage: print f.get(0), f.get(1), f.get(2)
            1 2 3
            sage: f1 = f[1:]
            sage: print f1.get(0), f1.get(1)
            2 3
            sage: f2 = f[2:]
            sage: print f2.get(0)
            3
            sage: f3 = f[3:]
            sage: print f3.get(0)
            Traceback (most recent call last):
            ...
            IndexError: lazy list index out of range

            sage: l = lazy_list([0]*12)[1::2]
            sage: l[2::3]
            lazy list [0, 0]
            sage: l[3::2]
            lazy list [0, 0]

        A lazy list automatically adjusts the indices in order that start and
        stop are congruent modulo step::

            sage: P = lazy_list(iter(Primes()))
            sage: P[1:12:4]._info()
            cache length 0
            start        1
            stop         13
            step         4
            sage: P[1:13:4]._info()
            cache length 0
            start        1
            stop         13
            step         4
            sage: P[1:14:4]._info()
            cache length 0
            start        1
            stop         17
            step         4
            sage: Q = P[100:1042233:12]
            sage: Q._info()
            cache length 0
            start        100
            stop         1042240
            step         12
            sage: R = Q[233::3]
            sage: R._info()
            cache length 0
            start        2896
            stop         1042252
            step         36
            sage: 1042252%36 == 2896%36
            True

        We check commutation::

            sage: l = lazy_list(iter(xrange(10000)))
            sage: l1 = l[::2][:3001]
            sage: l2 = l[:6002][::2]
            sage: l1._info()
            cache length 0
            start        0
            stop         6002
            step         2
            sage: l2._info()
            cache length 0
            start        0
            stop         6002
            step         2
            sage: l3 = l1[13::2][:50:2]
            sage: l4 = l1[:200][13:113:4]
            sage: l3._info()
            cache length 0
            start        26
            stop         226
            step         8
            sage: l4._info()
            cache length 0
            start        26
            stop         226
            step         8

        Further tests::

            sage: l = lazy_list(iter([0]*25))
            sage: l[2::3][2::3][4::5]
            lazy list []
            sage: l[2::5][3::][1::]
            lazy list [0]
            sage: l[3:24:2][1::][1:7:3]
            lazy list [0, 0]
            sage: l[::2][2::][2::3]
            lazy list [0, 0, 0]
            sage: l[4:3][:] is l[18:2]   # *the* empty_lazy_list
            True
        """
        if not isinstance(key, slice):
            return self.get(key)

        # the following make all terms > 0
        cdef Py_ssize_t start, stop, step
        start = 0 if key.start is None else key.start
        stop = PY_SSIZE_T_MAX if key.stop is None else key.stop
        step = 1 if key.step is None else key.step

        if step == 0:
            raise TypeError("step may not be 0")
        if step < 0 or start < 0 or stop < 0:
            raise ValueError("slice indices must be non negative")

        step = step * self.step
        start = self.start + start * self.step
        if stop != PY_SSIZE_T_MAX:
            stop = self.start + stop * self.step
        if stop > self.stop:
            stop = self.stop
        if stop != PY_SSIZE_T_MAX and stop%step != start%step:
            stop = stop - (stop-start)%step + step

        if stop <= start:
            return empty_lazy_list

        # here we return a slice of self. That is to say, a lazy list which
        # shares the same cache of values
        return self.make_linked_copy(
            start=start, stop=stop, step=step)


    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        This is the default implementation that calls ``_new_slice``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- the cache has now size larger than ``i``

        - ``1`` -- the lazy list is actually finite and shorter than ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_generic
            sage: L = lazy_list(Primes())[2:]
            sage: L.update_cache_up_to(4)
            0
            sage: L._info()
            cache length 5
            start        2
            stop         9223372036854775807
            step         1
        """
        if self.master is not None:    # this is a slice or other sublist
            return self.master._fit(i)

        cdef list l
        while len(self.cache) <= i:
            l = self._new_slice()
            if not l:
                return 1
            self.cache.extend(l)
        return 0


    cpdef list _get_cache_(self):
        r"""
        Return the internal cache.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: L = lazy_list(Primes()); L
            lazy list [2, 3, 5, ...]
            sage: L._get_cache_()
            [2, 3, 5, 7]
        """
        # is there a more direct way to get the cache in a Python
        # function (needed in lazy_list_takewhile.__init__) ?
        return self.cache


    def make_linked_copy(self, **kwds):
        r"""
        Return a lazy list which tracks this (``self``) lazy list, but
        shares the same cache.

        INPUT:

        - ``cls`` -- (default: ``None``) a class specifying the output.
          If ``None``, then :class:`lazy_list_generic` is used.

        Keyword arguments are passed on to the constructor of the
        new lazy list.

        OUTPUT:

        An instance of ``cls``.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: P = lazy_list(Primes())
            sage: Q = P.make_linked_copy()  # indirect doctest
            sage: P[5]
            13
            sage: Q._info()
            cache length 6
            start        0
            stop         9223372036854775807
            step         1
            sage: Q[10]
            31
            sage: P._info()
            cache length 11
            start        0
            stop         9223372036854775807
            step         1
            sage: Q._info()
            cache length 11
            start        0
            stop         9223372036854775807
            step         1

        ::

            sage: S = P[::2]; S
            lazy list [2, 5, 11, ...]
            sage: T = S.make_linked_copy(); T  # indirect doctest
            lazy list [2, 5, 11, ...]
            sage: T[::2]
            lazy list [2, 11, 23, ...]

        ::

            sage: from sage.misc.lazy_list import lazy_list_generic, lazy_list_from_iterator
            sage: class Z(lazy_list_generic):
            ....:     pass
            sage: i = lazy_list_from_iterator(iter(Primes()), cls=Z)
            sage: type(i)
            <type 'sage.misc.lazy_list.lazy_list_from_iterator'>
            sage: j = i.make_linked_copy()
            sage: type(j)
            <class '__main__.Z'>
            sage: j
            lazy list [2, 3, 5, ...]

        ::

            sage: class Y(lazy_list_generic):
            ....:     def __init__(self, a=0, **kwds):
            ....:         self.a = a
            ....:         lazy_list_generic.__init__(self, **kwds)
            sage: y = lazy_list(Primes(), cls=Y, cls_kwds={'a': 42})
            sage: y.a
            42
            sage: y._properties_()
            {...
             'cls': <class '__main__.Y'>,
             'cls_kwds': {'a': 42},
             ...}
        """
        properties = self._properties_()

        # We track self by using it as master in our new lazy_list and 
        # both share the same cache.
        properties['cache'] = self.cache
        properties['master'] = self

        cls_kwds = kwds['cls_kwds'] if kwds.has_key('cls_kwds') else properties['cls_kwds']
        properties.update(cls_kwds)
        properties.update(kwds)

        cls = properties['cls']
        return cls(**properties)


    def _properties_(self):
        r"""
        Return the values of all the keyword arguments used in the
        constructor.

        OUTPUT:

        A dictionary.

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list
            sage: lazy_list(Primes())._properties_()
            {'closing_delimiter': ']',
             'cls': <type 'sage.misc.lazy_list.lazy_list_generic'>,
             'cls_kwds': {},
             'more': '...',
             'name': 'lazy list',
             'opening_delimiter': '[',
             'preview': 3,
             'separator': ', ',
             'start': 0,
             'step': 1,
             'stop': 9223372036854775807}
        """
        return {
            'start': self.start,
            'stop': self.stop,
            'step': self.step,
            'name': self.name,
            'separator': self.separator,
            'more': self.more,
            'opening_delimiter': self.opening_delimiter,
            'closing_delimiter': self.closing_delimiter,
            'preview': self.preview,
            'cls': self.cls,
            'cls_kwds': self.cls_kwds}


    def change_class(self, cls=None):
        r"""
        Return a copy of this lazy list (sharing caches), which is an instance
        if the specified class ``cls``.

        INPUT:

        - ``cls`` -- (default: ``None``) a class specifying the
          output. If ``None``, then the default value of this lazy
          list (specified at construction and also called ``cls``) is
          used.

        OUTPUT:

        An instance of ``cls``.

        Note that if this lazy list is already an instance of ``cls``, then
        no copy is made.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_generic
            sage: P = lazy_list(Primes())
            sage: P.change_class() is P
            True
            sage: class Z(lazy_list_generic):
            ....:     pass
            sage: type(P.change_class(cls=Z))
            <class '__main__.Z'>
        """
        if cls is None:
            cls = self.cls
        if isinstance(self, cls):
            return self
        return self.make_linked_copy(cls=cls)


    def dropwhile(self, predicate):
        r"""
        Return a sublist (slice) but drop elements as long as the
        ``predicate`` is ``True``; afterwards take every element.

        INPUT:

        - ``predicate`` -- a callable returning ``True`` or ``False``.

        OUTPUT:

        A lazy list.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: lazy_list(Primes()).dropwhile(lambda x: x <= 5)
            lazy list [7, 11, 13, ...]

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_generic, lazy_list_from_iterator
            sage: class Z(lazy_list_generic):
            ....:     pass
            sage: lazy_list(Primes(), cls=Z).dropwhile(lambda x: x <= 5)
            lazy list [7, 11, 13, ...]
            sage: type(_)
            <class '__main__.Z'>
        """
        return lazy_list_dropwhile(
            self, predicate, **self._properties_()).change_class()


    def takewhile(self, predicate):
        r"""
        Return a sublist (slice) and take elements as long as
        the ``predicate`` is ``True``.

        INPUT:

        - ``predicate`` -- a callable returning ``True`` or ``False``.

        OUTPUT:

        A lazy list.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: lazy_list(Primes()).takewhile(lambda x: x <= 3)
            lazy list [2, 3]

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_generic, lazy_list_from_iterator
            sage: class Z(lazy_list_generic):
            ....:     pass
            sage: lazy_list(Primes(), cls=Z).takewhile(lambda x: x <= 3)
            lazy list [2, 3]
            sage: type(_)
            <class '__main__.Z'>
            sage: lazy_list(iter([4,3])).make_linked_copy(cls=Z)
            lazy list [4, 3]
            sage: type(_)
            <class '__main__.Z'>
            sage: from sage.misc.lazy_list import lazy_list_takewhile
        """
        return lazy_list_takewhile(
            self, predicate, **self._properties_()).change_class()


cdef class lazy_list_from_iterator(lazy_list_generic):

    def __init__(self, iterator, cache=None, stop=None, **kwds):
        r"""
        Lazy list built from an iterator.

        INPUT:

        - ``iterator`` -- an iterator

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy.

        - ``stop`` -- an optional stop point

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list
            sage: from itertools import count
            sage: m = lazy_list(count()); m
            lazy list [0, 1, 2, ...]

            sage: m2 = lazy_list(count(), start=8, stop=20551, step=2)
            sage: m2
            lazy list [8, 10, 12, ...]

            sage: x = iter(m)
            sage: print next(x), next(x), next(x)
            0 1 2
            sage: y = iter(m)
            sage: print next(y), next(y), next(y)
            0 1 2
            sage: print next(x), next(y)
            3 3
            sage: loads(dumps(m))
            lazy list [0, 1, 2, ...]

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_iterator
            sage: from itertools import count
            sage: lazy_list_from_iterator(count())
            lazy list [0, 1, 2, ...]
            sage: lazy_list_from_iterator(count(), ['a'], 10)
            lazy list ['a', 0, 1, ...]
            sage: _.info()
            cache length 4
            start        0
            stop         10
            step         1
        """
        self.iterator = iterator
        lazy_list_generic.__init__(self, cache, start=None, stop=stop, step=None, **kwds)


    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_iterator
            sage: L = lazy_list_from_iterator(iter(Primes()))
            sage: L.update_cache_up_to(10)
            0
            sage: L._info()
            cache length 11
            start        0
            stop         9223372036854775807
            step         1
        """
        while len(self.cache) <= i:
            try:
                o = next(self.iterator)
            except StopIteration:
                return 1
            self.cache.append(o)
        return 0


    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_iterator
            sage: from itertools import count
            sage: loads(dumps(lazy_list_from_iterator(count())))
            lazy list [0, 1, 2, ...]
            sage: loads(dumps(lazy_list_from_iterator(count(), ['a'])))
            lazy list ['a', 0, 1, ...]
        """
        return lazy_list_from_iterator, (self.iterator, self.cache, self.stop)


cdef class lazy_list_from_function(lazy_list_generic):

    def __init__(self, function, cache=None, stop=None, **kwds):
        r"""
        Lazy list built from a function.

        INPUT:

        - ``function`` -- a function that maps ``n`` to the element
          at position ``n``. (This
          function only needs to be defined for length larger than the length of
          the cache.)

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy.

        - ``stop`` -- an optional integer to specify the length of this lazy list.
          (Otherwise it is considered infinite).

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list_from_function
            sage: lazy_list_from_function(euler_phi)
            lazy list [0, 1, 1, ...]
            sage: lazy_list_from_function(divisors, [None])
            lazy list [None, [1], [1, 2], ...]

        TESTS::

            sage: def f(n):
            ....:     if n >= 5: raise StopIteration
            ....:     return 5 - n
            sage: list(lazy_list_from_function(f))
            [5, 4, 3, 2, 1]
        """
        self.callable = function
        lazy_list_generic.__init__(self, cache, stop=stop, **kwds)


    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_function
            sage: L = lazy_list_from_function(lambda x: 2*x); L
            lazy list [0, 2, 4, ...]
            sage: L.update_cache_up_to(10)
            0
            sage: L._info()
            cache length 11
            start        0
            stop         9223372036854775807
            step         1
        """
        while len(self.cache) <= i:
            self.cache.append(self.callable(len(self.cache)))


    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_function
            sage: loads(dumps(lazy_list_from_function(euler_phi)))
            lazy list [0, 1, 1, ...]
            sage: loads(dumps(lazy_list_from_function(divisors, [None])))
            lazy list [None, [1], [1, 2], ...]
        """
        if self.start != 0 or self.step != 1:
            raise RuntimeError
        return lazy_list_from_function, (self.callable, self.cache, self.stop)


cdef class lazy_list_from_update_function(lazy_list_generic):

    def __init__(self, function, cache=None, stop=None, **kwds):
        r"""
        Lazy list built from an update function.

        INPUT:

        - ``function`` -- a function that updates a list of precomputed values.
          The update function should take as input a list and make it longer
          (using either the methods ``append`` or ``extend``). If after a call
          to the update function the list of values is shorter a
          ``RuntimeError`` will occurr. If no value is added then the lazy list
          is considered finite.

        - ``cache`` -- an optional list to be used as the cache. Be careful that
          there is no copy.

        - ``stop`` -- an optional integer to specify the length of this lazy list
          (otherwise it is considered infinite)

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_update_function
            sage: def update_function(values):
            ....:     n = len(values)+1
            ....:     values.extend([n]*n)
            sage: l = lazy_list_from_update_function(update_function)
            sage: l[:20].list()
            [1, 2, 2, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 16, 16, 16, 16, 16]
        """
        self.update_function = function
        lazy_list_generic.__init__(self, cache,
                                   start=None, stop=stop, step=None,
                                   **kwds)


    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_from_update_function
            sage: L = lazy_list_from_update_function(lambda V: V.extend([1, 2])); L
            lazy list [1, 2, 1, ...]
            sage: L.update_cache_up_to(10)
            0
            sage: L._info()
            cache length 12
            start        0
            stop         9223372036854775807
            step         1
        """
        cdef Py_ssize_t l,ll
        l = len(self.cache)
        while l <= i:
            self.update_function(self.cache)
            ll = len(self.cache)
            if ll < l:
                raise RuntimeError("the update function made the cache shorter")
            elif l == ll:
                return 1
            l = ll
        return 0


    def __reduce__(self):
        r"""
        TESTS::

            sage: from sage.misc.lazy_list import lazy_list

            sage: def my_update_function(values): values.append(ZZ(len(values)).is_prime())
            sage: l = lazy_list(update_function=my_update_function)
            sage: l[4]
            False
            sage: loads(dumps(l))   # not tested (works in console though)
            lazy list [False, False, True, ...]

            sage: def say_hey(cache): print "hey"
            sage: l = lazy_list(update_function=say_hey, initial_values=range(10))
            sage: l._fit(10)
            hey
            1
            sage: l._info()
            cache length 10
            start        0
            stop         10
            step         1
            sage: l2 = loads(dumps(l))   # not tested
            sage: l2._info()             # not tested
            sage: l2._info()             # not tested
            cache length 10
            start        0
            stop         10
            step         1
            sage: l.list() == l2.list()  # not tested
            True
        """
        if self.start != 0 or self.step != 1:
            raise RuntimeError
        return lazy_list_from_update_function, (self.update_function, self.cache, self.stop)


cdef class lazy_list_takewhile(lazy_list_generic):

    def __init__(self, master, predicate, **kwds):
        r"""
        A lazy sublist (slice), which takes elements as long as
        the given ``predicate`` is ``True``.

        INPUT:

        - ``master`` -- a lazy list.

        - ``predicate`` -- a callable returning ``True`` or ``False``.

        Keyword arguments are passed on to the constructor of
        :class:`lazy_list_generic`.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_takewhile
            sage: L = lazy_list(Primes())
            sage: lazy_list_takewhile(L, lambda x: x <= 3)
            lazy list [2, 3]

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list_generic
            sage: class Z(lazy_list_generic):
            ....:     def __init__(self):
            ....:         self.i = -2
            ....:         lazy_list_generic.__init__(self)
            ....:     def _new_slice(self):
            ....:         self.i += 2
            ....:         return range(self.i, -1, -1)
            sage: list(Z()[:16])
            [0, 2, 1, 0, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0]
            sage: list(lazy_list_takewhile(Z(), lambda x: x <= 3))
            [0, 2, 1, 0]
            sage: list(lazy_list_takewhile(Z(), lambda x: x <= 2))
            [0, 2, 1, 0]
            sage: lazy_list_takewhile(Z(), lambda x: x != 3)[4]
            4
        """
        lazy_list_generic.__init__(
            self, master=master, cache=master._get_cache_(), **kwds)
        self.predicate = predicate
        self.taking = True
        self.to_test = self.start


    cpdef int update_cache_up_to(self, Py_ssize_t i) except -1:
        r"""
        Update the cache up to ``i``.

        OUTPUT:

        - ``-1`` -- a Python error occurred

        - ``0`` -- everything went fine

        - ``1`` -- the iterator stopped before ``i``

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_takewhile
            sage: L = lazy_list_takewhile(
            ....:    lazy_list(Primes()), lambda x: x <= 10); L
            lazy list [2, 3, 5, ...]
            sage: L._info()
            cache length 4
            start        0
            stop         9223372036854775807
            step         1
            sage: L.update_cache_up_to(10)
            1
            sage: tuple(L)
            (2, 3, 5, 7)
            sage: L._info()
            cache length 5
            start        0
            stop         4
            step         1

        ::

            sage: P = lazy_list(Primes())
            sage: a = tuple(P.dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20)); a
            (11, 13, 17, 19)
            sage: b = tuple(P.dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20))
            sage: a == b
            True

        ::

            sage: F = lazy_list(lambda x: fibonacci(x))
            sage: tuple(F.dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20))
            (13,)

        ::

            sage: P = lazy_list(srange(1,30,2))
            sage: a = tuple(P.dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20)); a
            (11, 13, 15, 17, 19)
            sage: b = tuple(P.dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20))
            sage: a == b
            True

        ::

            sage: P = lazy_list(Primes())
            sage: Q = P.takewhile(lambda x: x < 20); c = tuple(Q); c
            (2, 3, 5, 7, 11, 13, 17, 19)
            sage: Q._info()
            cache length 9
            start        0
            stop         8
            step         1
            sage: Q = P.takewhile(lambda x: x < 20); d = tuple(Q)
            sage: Q._info()
            cache length 9
            start        0
            stop         8
            step         1
            sage: c == d
            True
        """
        if not self.taking:
            return 0

        if self.master._fit(self.to_test):
            self.stop = self.to_test
            self.taking = False
            return 1
        self.to_test = max(self.to_test, self.start)

        while self.to_test <= i:
            if self.master._fit(self.to_test):
                self.stop = self.to_test
                self.taking = False
                return 1
            if not (self.to_test < len(self.cache) and
                    self.predicate(self.cache[self.to_test])):
                break
            self.to_test += 1

        if self.to_test <= i:
            self.stop = self.to_test
            self.taking = False
            return 1

        return 0


cdef class lazy_list_dropwhile(lazy_list_generic):

    def __init__(self, master, predicate, **kwds):
        r"""
        A lazy sublist (slice) which drops elements as long as the
        ``predicate`` is ``True``; afterwards takes every element.

        INPUT:

        - ``master`` -- a lazy list.

        - ``predicate`` -- a callable returning ``True`` or ``False``.

        Keyword arguments are passed on to the constructor of
        :class:`lazy_list_generic`.

        EXAMPLES::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_dropwhile
            sage: L = lazy_list(Primes())
            sage: lazy_list_dropwhile(L, lambda x: x <= 3)
            lazy list [5, 7, 11, ...]

        TESTS::

            sage: lazy_list_dropwhile(lazy_list(Primes())[3:], lambda x: x <= 13)
            lazy list [17, 19, 23, ...]
            sage: A = lazy_list_dropwhile(lazy_list(Primes())[3:], lambda x: x <= 13)
            sage: A[:10].list()
            [17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
            sage: A
            lazy list [17, 19, 23, ...]
            sage: A[:10].list()
            [17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
            sage: A._info()
            cache length 16
            start        6
            stop         9223372036854775807
            step         1

        ::

            sage: lazy_list(Primes()).dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20)
            lazy list [11, 13, 17, ...]
            sage: tuple(lazy_list(Primes()).dropwhile(lambda x: x < 10).takewhile(lambda x: x < 20))
            (11, 13, 17, 19)
        """
        lazy_list_generic.__init__(
            self, master=master, cache=master._get_cache_(), **kwds)
        self.predicate = predicate
        self.dropping = True


    cpdef int _fit(self, Py_ssize_t n) except -1:
        r"""
        Fill the cache making the term at index ``n`` available.

        You can access the term at position ``n`` from the cache when it returns
        ``0``.

        OUTPUT:

        - ``1`` -- the lazy list is actually finite and shorter than ``n``

        - ``0`` -- you can safely access the term at position ``n`` after this call

        - ``-1`` -- to handle Python errors (you can ignore it in Python code)

        TESTS::

            sage: from sage.misc.lazy_list import lazy_list, lazy_list_dropwhile
            sage: L = lazy_list_dropwhile(
            ....:    lazy_list(Primes()), lambda x: x <= 10); L
            lazy list [11, 13, 17, ...]
            sage: L._fit(10)
            0
            sage: L._info()
            cache length 11
            start        4
            stop         9223372036854775807
            step         1
        """
        cdef Py_ssize_t m
        if self.dropping:
            m = self.start
            lazy_list_generic._fit(self, m)
            while m < len(self.cache) and self.predicate(self.cache[m]):
                m += 1
                lazy_list_generic._fit(self, m)
            self.start = m
            self.dropping = False
            n += m

        return lazy_list_generic._fit(self, n)
