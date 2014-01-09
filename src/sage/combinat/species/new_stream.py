from sage.rings.all import ZZ
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

# TODO:
# 1. __len__ / number_computed / max_computed??


class Stream(SageObject):
    def __init__(self, compute=None):
        """
        """
        if compute is not None:
            self.compute = compute

        self._constant = False

    def compute(self, n):
        raise NotImplementedError

    def __setitem__(self, n, value):
        raise NotImplementedError

    # Will be changed to staticmethod later
    def getitem_decorator(func):
        from sage.misc.all import sage_wraps
        @sage_wraps(func)
        def wrapper(self, n):
            # Handle constant
            if self._constant is not False:
                pos, value = self._constant
                if n >= pos:
                    return value
            return func(self, n)
        return wrapper

    @getitem_decorator
    def __getitem__(self, n):
        return self.compute(n)

    getitem_decorator = staticmethod(getitem_decorator)

    def is_constant(self):
        """
        Returns True if this stream is eventually constant.
        """
        return self._constant is not False

    def set_constant(self, n, value):
        self._constant = (n, value)

    def get_constant(self):
        assert self._constant is not False
        return self._constant[1]

    def get_constant_position(self):
        assert self._constant is not False
        return self._constant[0]

    def __iter__(self):
        i = 0
        while True:
            try:
                yield self[i]
            except IndexError:
                break
            i += 1
        raise StopIteration

    def map(self, func):
        """
        EXAMPLES::

            sage: from sage.combinat.species.new_stream import Stream
            sage: s = Stream(compute=lambda x: x)
            sage: ss = s.map(lambda x: x*x)
            sage: [ss[i] for i in range(10)]
            [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]
        """
        return MappedStream(func, self)

class ConstantStream(Stream):
    def __init__(self, constant):
        super(ConstantStream, self).__init__()
        self.set_constant(0, constant)

class DictCachedStream(Stream):
    def __init__(self, **kwds):
        self._cache = {}
        super(DictCachedStream, self).__init__(**kwds)

    def __setitem__(self, n, value):
        self._cache[n] = value

    @Stream.getitem_decorator
    def __getitem__(self, n):
        value = self._cache.get(n, None)
        if value is None:
            value = self.compute(n)
            self._cache[n] = value
            return value
        else:
            return value
            
class MappedStream(DictCachedStream):
    def __init__(self, func, stream):
        self._func = func
        self._stream = stream
        super(MappedStream, self).__init__()

    def compute(self, n):
        value = self._func(self._stream[n])
        if self._constant is False and self._stream.is_constant():
            self._constant = (n, value)
        return value


class ListCachedStream(Stream):
    def __init__(self, **kwds):
        self._cache = []
        self._cached = kwds.pop('cached', True)
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
            _ = self[pos]
            pos += 1
        self._cache[n] = value

    def __len__(self):
        return len(self._cache)

    number_computed = __len__

    @Stream.getitem_decorator
    def __getitem__(self, n):
        if not self._cached:
            return self.compute(n)
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
