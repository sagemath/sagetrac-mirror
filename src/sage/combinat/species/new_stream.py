from sage.structure.sage_object import SageObject

# TODO:
# 1. __len__ / number_computed / max_computed??


def check_constant_decorator(func):
    """
    A method decorator for ``__getitem__`` which checks is the stream
    is (eventually) constant before computing the $n^{th}$
    coefficient.

    EXAMPLES::

        sage: from sage.combinat.species.new_stream import Stream, check_constant_decorator
        sage: s = Stream()
        sage: def foo(self, n):
        ....:     return self.compute(n)
        sage: import types
        sage: s.foo = types.MethodType(check_constant_decorator(foo), s)
        sage: s.set_constant(1, 3)
        sage: s.compute(5)
        Traceback (most recent call last):
        ...
        NotImplementedError
        sage: s.foo(5)
        3

    """
    from sage.misc.all import sage_wraps
    @sage_wraps(func)
    def wrapper(self, n):
        if self._constant is not False or self.is_constant():
            if self._constant is False:
                self.set_constant(self.get_constant_position(),
                                  self.get_constant())
            pos, value = self._constant
            if n >= pos:
                return value
        return func(self, n)
    return wrapper


class Stream(SageObject):
    def __init__(self):
        """
        A base class for streams.  This class is typically subclassed.
        """
        self._constant = False

    def compute(self, n):
        """
        Computes the $n^{th}$ coefficient of this stream.  This should
        be overridden by subclasses.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.compute(2)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __setitem__(self, n, value):
        """
        Sets the $n^{th}$ coefficient of this stream to ``value``.
        This should be overridden by subclasses.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s[0] = 2
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    @check_constant_decorator
    def __getitem__(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import Stream
            sage: class MyStream(Stream):
            ....:     def compute(self, n):
            ....:         return n
            sage: s = MyStream()
            sage: s[10]
            10
        """
        return self.compute(n)

    def is_constant(self):
        """
        Returns True if this stream is eventually constant.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.is_constant()
            False
            sage: s.set_constant(2, 4)
            sage: s.is_constant()
            True
        """
        return self._constant is not False

    def set_constant(self, n, value):
        """
        Sets this stream to be eventually constant at coefficient
        ``n`` with value ``value``.

        EXAMPLES::
        
            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.set_constant(0, 2)
            sage: s.get_constant()
            2
            sage: s.get_constant_position()
            0
            sage: s[3]
            2
        """
        self._constant = (n, value)

    def get_constant(self):
        """
        Returns the constant value if this stream is eventually
        constant.

        EXAMPLES::
        
            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.set_constant(0, 2)
            sage: s.get_constant()
            2
        """
        assert self._constant is not False
        return self._constant[1]

    def get_constant_position(self):
        """
        Returns the position where this stream is eventually constant.

        EXAMPLES::
        
            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.set_constant(1, 2)
            sage: s.get_constant_position()
            1
        """
        assert self._constant is not False
        return self._constant[0]

    def __iter__(self):
        """
        Returns an iterator for this stream.

        EXAMPLES::
        
            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream()
            sage: s.set_constant(0, 2)
            sage: it = iter(s)
            sage: [it.next() for i in range(5)]
            [2, 2, 2, 2, 2]

        """
        i = 0
        while True:
            try:
                yield self[i]
            except IndexError:
                break
            i += 1
        raise StopIteration

class ListCachedStream(Stream):
    def __init__(self, **kwds):
        """
        Returns a stream whose computed values are cached in a list.
        Additionally, when the $n^{th}$ coefficient is requested, it
        guarantees that all of the coefficients up to $n$ have been
        computed.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc, ListCachedStream
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: isinstance(s, ListCachedStream)
            True
            sage: s[5]
            8
            sage: s._cache
            [1, 1, 2, 3, 5, 8]

        """
        self._cache = []
        super(ListCachedStream, self).__init__(**kwds)

    def __setitem__(self, n, value):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc, ListCachedStream
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: isinstance(s, ListCachedStream)
            True
            sage: s[1]
            1
            sage: s[0] = 2
            sage: s[0]
            2
            sage: s[10]
            123
            sage: s[15] = 100
            sage: s[15]
            100
        """
        pos = len(self._cache)
        while n >= pos:
            self[pos]
            pos += 1
        self._cache[n] = value

    def number_computed(self):
        """
        Returns the number of coefficients that have been computed so
        far.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc, ListCachedStream
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: isinstance(s, ListCachedStream)
            True
            sage: s[5]
            8
            sage: s.number_computed()
            6
        """
        return len(self._cache)

    __len__ = number_computed

    @check_constant_decorator
    def __getitem__(self, n):
        """
        Returns the $n^{th}$ coefficient of this stream, checking the
        cache before trying to compute the value.  This method
        guarantees that all of the coefficients up to $n$ have been
        computed first.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc, ListCachedStream
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: isinstance(s, ListCachedStream)
            True
            sage: s[5]
            8
            sage: s._cache
            [1, 1, 2, 3, 5, 8]

        We check to see that values are indeed returned from the cache
        if already computed::

            sage: def foo(self, n):
            ....:    raise NotImplementedError
            sage: s.compute = foo
            sage: s[2]
            2
        """
        pos = len(self._cache)
        while pos <= n:
            value = self.compute(pos)
            self._cache.append(value)
            pos += 1
        return self._cache[n]

class StreamFromIterator(ListCachedStream):
    def __init__(self, iterator=None, **kwds):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromIterator
            sage: from sage.combinat.species.stream import _integers_from
            sage: s = StreamFromIterator(_integers_from(1))
            sage: s[0]
            1
            sage: s[10]
            11
        """
        self._it = iterator
        super(StreamFromIterator, self).__init__(**kwds)            
    
    def compute(self, n):
        """
        EXAMPLES:

        We test to make sure that iterator which finish iterating are
        constnat for the rest of the values::

            sage: from sage.combinat.species.new_stream import StreamFromIterator
            sage: s = StreamFromIterator(iter([1,2,3]))
            sage: s[0]
            1
            sage: s.is_constant()
            False
            sage: s[2], s[10]
            (3, 3)
            sage: s.is_constant()
            True
        """
        # ListCachedStream verifies that compute will be called with n in order
        assert n == len(self._cache), "compute called out of order"

        try:
            return self._it.next()
        except StopIteration:
            value = self._cache[-1]
            self.set_constant(len(self._cache) - 1, value)
            return value
            

class StreamFromFunc(ListCachedStream):
    def __init__(self, func=None, **kwds):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: s[0]
            1
            sage: s[1]
            1
            sage: s[2]
            2
            sage: s[10]
            89
        """
        self._func = func
        super(StreamFromFunc, self).__init__(**kwds)

    def compute(self, n):
        """
        .. note::
       
           This should not be called directly.  Instead, you should
           use :meth:`__getitem__`.

        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromFunc
            sage: h = lambda l: 1 if len(l) < 2 else l[-1] + l[-2]
            sage: s = StreamFromFunc(h)
            sage: s.compute(0)
            1
            sage: s.compute(2)
            Traceback (most recent call last):
            ...
            AssertionError: compute called out of order
        """
        # ListCachedStream verifies that compute will be called with n in order
        assert n == len(self._cache), "compute called out of order"
        return self._func(self._cache)
        

class StreamFromList(ListCachedStream):
    def __init__(self, list=None, **kwds):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_stream import StreamFromList
            sage: s = StreamFromList([1,2,3])
            sage: s[0]
            1
            sage: s[5]
            3

        """
        super(StreamFromList, self).__init__(**kwds)
        if list is None:
            list = self.list()
        if len(list) < 0:
            raise ValueError, "list cannot be empty"
        self._cache = list[:]
        self.set_constant(len(list) - 1, list[-1])


def OldStreamBehavior(x=None, const=None):
    import types
    if const is not None:
        return ConstantStream(const)
    elif isinstance(x, list):
        return StreamFromList(x)
    elif hasattr(x, '__iter__'):
        return StreamFromIterator(iter(x))
    elif isinstance(x, (types.FunctionType, types.LambdaType)):
        return StreamFromFunc(x)
    else:
        return StreamFromIterator(iter([x,0]))
