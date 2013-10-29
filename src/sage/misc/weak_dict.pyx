"""
Fast and safe weak value dictionary

AUTHORS:

- Simon King (2013-10)

Python's :mod:`weakref` module provides
:class:`~weakref.WeakValueDictionary`. This behaves similar to a dictionary,
but it does not prevent its values from garbage collection. Hence, it stores
the values by weak references with callback functions: The callback function
deletes a key-value pair from the dictionary, as soon as the value becomes
subject to garbage collection.

However, a problem arises if hash and comparison of the key depend on the
value that is being garbage collected::

    sage: import weakref
    sage: class Vals(object): pass
    sage: class Keys:
    ....:     def __init__(self, val):
    ....:         self.val = weakref.ref(val)
    ....:     def __hash__(self):
    ....:         return hash(self.val())
    ....:     def __eq__(self, other):
    ....:         return self.val() == other.val()
    sage: ValList = [Vals() for _ in range(10)]
    sage: D = weakref.WeakValueDictionary()
    sage: for v in ValList:
    ....:     D[Keys(v)] = v
    sage: len(D)
    10
    sage: del ValList, v
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    Exception KeyError: (<__main__.Keys instance at ...>,) in <function remove at ...> ignored
    ...
    sage: len(D) > 1
    True

Hence, there are scary error messages, and moreover the defunct items have not
been removed from the dictionary.

Therefore, Sage provides an alternative implementation
:class:`sage.misc.weak_dict.WeakValueDictionary`, using a callback that
removes the defunct item not based on hash and equality check of the key (this
is what fails in the example above), but based on comparison by identity. This
is possible, since references with callback function are distinct even if they
point to the same object. Hence, even if the same object ``O`` occurs as value
for several keys, each reference to ``O`` corresponds to a unique key. We see
no error messages, and the items get correctly removed::

    sage: ValList = [Vals() for _ in range(10)]
    sage: import sage.misc.weak_dict
    sage: D = sage.misc.weak_dict.WeakValueDictionary()
    sage: for v in ValList:
    ....:     D[Keys(v)] = v
    sage: len(D)
    10
    sage: del ValList
    sage: len(D)
    1
    sage: del v
    sage: len(D)
    0

Another proplem arises when iterating over the items of a dictionary: If
garbage collection occurs during iteration, then the content of the dictionary
changes, and the iteration breaks for :class:`weakref.WeakValueDictionary`::

    sage: class Cycle:
    ....:     def __init__(self):
    ....:         self.selfref = self
    sage: C = [Cycle() for n in range(10)]
    sage: D = weakref.WeakValueDictionary(enumerate(C))
    sage: import gc
    sage: gc.disable()
    sage: del C[:5]
    sage: len(D)
    10
    sage: for k in D.iterkeys():
    ....:     gc.enable()
    ....:     _ = gc.collect()
    Traceback (most recent call last):
    ...
    RuntimeError: dictionary changed size during iteration

With :class:`~sage.misc.weak_dict.WeakValueDictionary`, the behaviour is
safer. During iteration, we can delete one item from the list ``C`` and make
sure that the cyclic garbage collector removes it, we can directly remove a
not yet visited item from the dictionary, we can pop a prescribed not yet
visited item from the dictionary, and we can let the dictionary pop an
arbitrary item. All the time, the length of the dictionary is correctly
reported, and it does not happen that we visit the key of an item that is
subject to deletion::

    sage: C = [Cycle() for n in range(5)]
    sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(C))
    sage: remaining_keys = set(range(5))
    sage: for k in D:
    ....:     print "len",len(D)
    ....:     remaining_keys.remove(k)
    ....:     if len(D)==5:
    ....:         C[k] = None
    ....:         C[remaining_keys.pop()] = None
    ....:         _ = gc.collect()
    ....:     elif len(D)==3:
    ....:         del D[remaining_keys.pop()]
    ....:     elif len(D)==2:
    ....:         l,v = D.popitem()
    len 5
    len 3
    len 2

Note that Sage's weak value dictionary is actually an instance of
:class:`dict`, in contrast to :mod:`weakref`'s weak value dictionary::

    sage: issubclass(weakref.WeakValueDictionary, dict)
    False
    sage: issubclass(sage.misc.weak_dict.WeakValueDictionary, dict)
    True

"""
########################################################################
#       Copyright (C) 2013 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################

import weakref
from weakref import KeyedRef

from cpython.dict cimport *
from cpython.weakref cimport *
from cpython.list cimport *

cdef PyObject* Py_None = <PyObject*>None

cdef class WeakValueDictionary(dict):
    """
    IMPLEMENTATION:

    The :class:`WeakValueDictionary` inherits from :class:`dict`. As a
    :class:`dict`, it stores lists that are indexed by the hash of the given
    key. These lists (also known as "hash buckets") are organised as ``[key_1,
    ref_value_1, key_2, ref_value_2, ... ]``.

    ``ref_value_n`` is a keyed weak reference to the `n`-th value. The key of
    this weak reference is the hash of ``key_n``. The callback of the weak
    reference thus has enough information to find the correct hash bucket. It
    will then search the hash bucket for the weak reference that is subject to
    the callback, and will delete the corresponding key-reference pair from
    the hash bucket.

    NOTE:

    In contrast to :class:`weakref.WeakValueDictionary` in Python's
    :mod:`weakref` module, the callback does not need to assume that the
    dictionary key is a valid Python object when it is called. There is no
    need to compute the hash or compare the dictionary keys. This is why
    the example below would not work with
    :class:`weakref.WeakValueDictionary`, but does work with
    :class:`sage.misc.weak_dict.WeakValueDictionary`.

    EXAMPLES::

        sage: import weakref
        sage: class Vals(object): pass
        sage: class Keys:
        ....:     def __init__(self, val):
        ....:         self.val = weakref.ref(val)
        ....:     def __hash__(self):
        ....:         return hash(self.val())
        ....:     def __eq__(self, other):
        ....:         return self.val() == other.val()
        sage: ValList = [Vals() for _ in range(10)]
        sage: import sage.misc.weak_dict
        sage: D = sage.misc.weak_dict.WeakValueDictionary()
        sage: for v in ValList:
        ....:     D[Keys(v)] = v
        sage: len(D)
        10
        sage: del ValList
        sage: len(D)
        1
        sage: del v
        sage: len(D)
        0

    """
    cdef callback
    cdef int _guard_level
    cdef dict _pending_removals
    cdef _iteration_context

    def __init__(self, data=()):
        """
        Create a :class:`WeakValueDictionary` with given initial data.

        INPUT:

        Optional iterable of key-value pairs

        EXAMPLES::

            sage: L = [(p,GF(p)) for p in prime_range(10)]
            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: len(D)
            0
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)
            sage: len(D) == len(L)
            True

        """
        dict.__init__(self)
        # Define a callback function. In contrast to what is done in Python's
        # weakref.WeakValueDictionary, we use a closure and not a weak
        # reference to refer to self. Reason: We trust that we will not create
        # sub-classes of keyed references or weak value dictionaries providing
        # a __del__ method.
        def callback(r):
            cdef WeakValueDictionary cself = self
            pending = False
            if cself._guard_level:
                cself._pending_removals[<Py_ssize_t><void*>r] = r
                pending = True
                #return
            cdef int hashk = r.key
            cdef Py_ssize_t idr = <Py_ssize_t><void *>r
            cdef void * buckref = PyDict_GetItem(cself, hashk)
            if buckref==NULL:
                return
            cdef list bucket = <list>buckref
            cdef Py_ssize_t i,l
            l = len(bucket)
            for i from 0 <= i < l by 2:
                if <Py_ssize_t>PyList_GetItem(bucket,i+1) == idr:
                    if pending:
                        return
                    del bucket[i:i+2]
                    if not bucket:
                        PyDict_DelItem(cself,hashk)
                    return
        self.callback = callback
        self._guard_level = 0
        self._pending_removals = {}
        self._iteration_context = _IterationContext(self)
        for k,v in data:
            self[k] = v

    def __repr__(self):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: repr(sage.misc.weak_dict.WeakValueDictionary([(1,ZZ),(2,QQ)]))  # indirect doctest
            '<WeakValueDictionary at 0x...>'

        """
        return "<WeakValueDictionary at 0x%x>" % id(self)

    def __str__(self):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: str(sage.misc.weak_dict.WeakValueDictionary([(1,ZZ),(2,QQ)]))  # indirect doctest
            '<WeakValueDictionary at 0x...>'

        """
        return "<WeakValueDictionary at 0x%x>" % id(self)

    def setdefault(self, k, default=None):
        """
        Return the stored value for a given key; return and store a default
        value if no previous value is stored.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [(p,GF(p)) for p in prime_range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)
            sage: len(D)
            4

        The value for an existing key is returned and not overridden::

            sage: D.setdefault(5, ZZ)
            Finite Field of size 5
            sage: D[5]
            Finite Field of size 5

        For a non-existing key, the default value is stored and returned::

            sage: D.has_key(4)
            False
            sage: D.setdefault(4, ZZ)
            Integer Ring
            sage: D.has_key(4)
            True
            sage: D[4]
            Integer Ring
            sage: len(D)
            5

        """
        cdef int hashk = hash(k)
        cdef void * buckref = PyDict_GetItem(self, hashk)
        cdef list bucket
        cdef Py_ssize_t i,l
        cdef PyObject* out
        if buckref==NULL:
            bucket = [k, KeyedRef(default, self.callback, hashk)]
            PyDict_SetItem(self,hashk,bucket)
            return default
        bucket = <list>buckref
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket, i) == k:
                out = PyWeakref_GetObject(<object>PyList_GetItem(bucket, i+1))
                if out!=Py_None:
                    return <object>out
                bucket[i+1] = KeyedRef(default, self.callback, hashk)
                return default
        if self._guard_level:
            raise RuntimeError("setdefaul(): Can not add items while iterating over the dictionary")
        bucket.extend([k, KeyedRef(default, self.callback, hashk)])
        return default

    def __setitem__(self, k, v):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: ZZ in D
            False

        One can set new items::

            sage: D[ZZ] = QQ   # indirect doctest
            sage: D[ZZ]
            Rational Field
            sage: len(D)
            1
            sage: ZZ in D
            True

       One can also override existing items::

           sage: D[ZZ] = RLF
           sage: ZZ in D
           True
           sage: D[ZZ]
           Real Lazy Field
           sage: len(D)
           1

        TESTS:

        One may wonder whether it causes problems when garbage collection for
        a previously exististing item happens *after* overriding the
        item. The example shows that it is not a problem::

            sage: class Cycle:
            ....:     def __init__(self):
            ....:         self.selfref = self
            sage: L = [Cycle() for _ in range(5)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: len(D)
            5
            sage: import gc
            sage: gc.disable()
            sage: del L
            sage: len(D)
            5
            sage: D[2] = ZZ
            sage: len(D)
            5
            sage: gc.enable()
            sage: _ = gc.collect()
            sage: len(D)
            1
            sage: D.items()
            [(2, Integer Ring)]

        """
        cdef int hashk = hash(k)
        cdef void * buckref = PyDict_GetItem(self, hashk)
        cdef list bucket
        if buckref==NULL:
            bucket = []
            PyDict_SetItem(self,hashk,bucket)
        else:
            bucket = <list>buckref
        cdef object k0
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket, i) == k:
                bucket[i+1] = KeyedRef(v, self.callback, hashk)
                return
        if self._guard_level:
            raise RuntimeError("Can not add items while iterating over the dictionary")
        bucket.extend((k, KeyedRef(v, self.callback, hashk)))

    def __delitem__(self, k):
        """
        TESTS::

            sage: import sage.misc.weak_dict
            sage: L = [GF(p) for p in prime_range(10^3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 4 in D
            True
            sage: D[4]
            Finite Field of size 11
            sage: del D[4]
            sage: 4 in D
            False
            sage: D[4]
            Traceback (most recent call last):
            ...
            KeyError: 4

        """
        cdef int hashk = hash(k)
        cdef void * buckref = PyDict_GetItem(self, hashk)
        if buckref==NULL:
            raise KeyError(k)
        cdef list bucket = <list>buckref
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i) == k:
                if self._guard_level:
                    v = <object>PyList_GetItem(bucket,i+1)
                    if <Py_ssize_t><void*>v in self._pending_removals:
                        raise KeyError(k)
                    self._pending_removals[<Py_ssize_t><void*>v] = v
                    return
                del bucket[i:i+2]
                if not bucket:
                    PyDict_DelItem(self,hashk)
                return
        raise KeyError(k)

    def pop(self, k):
        """
        Return the value for a given key, and delete it from the dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [GF(p) for p in prime_range(10^3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 20 in D
            True
            sage: D.pop(20)
            Finite Field of size 73
            sage: 20 in D
            False
            sage: D.pop(20)
            Traceback (most recent call last):
            ...
            KeyError: 20

        """
        cdef int hashk = hash(k)
        cdef void * buckref = PyDict_GetItem(self, hashk)
        if buckref==NULL:
            raise KeyError(k)
        cdef list bucket = <list>buckref
        cdef Py_ssize_t i,l
        cdef PyObject * out
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket, i) == k:
                v = <object>PyList_GetItem(bucket, i+1)
                if not self._guard_level:
                    del bucket[i:i+2]
                elif <Py_ssize_t><void*>v in self._pending_removals:
                    break
                else:
                    self._pending_removals[<Py_ssize_t><void*>v] = v
                out = PyWeakref_GetObject(v)
                if not bucket:
                    PyDict_DelItem(self,hashk)
                if out != Py_None:
                    return <object>out
                break
        raise KeyError(k)

    def popitem(self):
        """
        Return and delete some item from the dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[1] = ZZ

        The dictionary only contains a single item, hence, it is clear which
        one will be returned::

            sage: D.popitem()
            (1, Integer Ring)

        Now, the dictionary is empty, and hence the next attempt to pop an
        item will fail with a ``KeyError``::

            sage: D.popitem()
            Traceback (most recent call last):
            ...
            KeyError: 'popitem(): weak value dictionary is empty'

        """
        for k,v in self.iteritems():
            del self[k]
            return k, v
        raise KeyError('popitem(): weak value dictionary is empty')
#        if self._guard_level:
#            raise RuntimeError("popitem(): not allowed while iterating over the dictionary")
#        cdef Py_ssize_t hashk
#        cdef list bucket
#        cdef PyObject *basekey, *bucketref
#        cdef Py_ssize_t pos = 0
#        if not PyDict_Next(self, &pos, &basekey, &bucketref):
#            raise KeyError('popitem(): weak value dictionary is empty')
#        bucket = <list>bucketref
#        if not bucket:
#            raise KeyError('popitem(): weak value dictionary is empty')
#        k = bucket.pop(0)
#        v = bucket.pop(0)
#        if not bucket:
#            PyDict_DelItem(self,<object>basekey)
#        cdef PyObject* out = PyWeakref_GetObject(v)
#        if out!=Py_None:
#            return k, <object>out
#        return self.popitem()

    def get(self, k, d=None):
        """
        Return the stored value for a key, or a default value for unkown keys.

        The default value defaults to None.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: L = [GF(p) for p in prime_range(10^3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 100 in D
            True
            sage: 200 in D
            False
            sage: D.get(100, "not found")
            Finite Field of size 547
            sage: D.get(200, "not found")
            'not found'
            sage: D.get(200) is None
            True

        """
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            return d
        cdef list bucket = <list>buckref
        cdef PyObject* out
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i)==k:
                v = <object>PyList_GetItem(bucket,i+1)
                if <Py_ssize_t><void*>v in self._pending_removals:
                    return d
                out = PyWeakref_GetObject(v)
                if out!=Py_None:
                    return <object>out
                break
        return d

    def __getitem__(self, k):
        """
        TESTS::

            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: D[ZZ] = QQ
            sage: D[QQ]
            Traceback (most recent call last):
            ...
            KeyError: Rational Field
            sage: D[ZZ]     # indirect doctest
            Rational Field

        As usual, the dictionary keys are compared by `==` and not by
        identity::

            sage: D[10] = ZZ
            sage: D[int(10)]
            Integer Ring

        """
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            raise KeyError(k)
        cdef list bucket = <list>buckref
        cdef PyObject* out
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i)==k:
                v = <object>PyList_GetItem(bucket,i+1)
                if <Py_ssize_t><void*>v in self._pending_removals:
                    raise KeyError(k)
                out = PyWeakref_GetObject(v)
                if out!=Py_None:
                    return <object>out
                break
        raise KeyError(k)

    def has_key(self, k):
        """
        Returns True, if the key is known to the dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: D.has_key(3)
            True

        As usual, keys are compared by equality and not by identity::

            sage: D.has_key(int(3))
            True

        This is a weak value dictionary. Hence, the existence of the
        dictionary does not prevent the values from garbage collection,
        thereby removing the corresponding key-value pairs::

            sage: del L[3]
            sage: D.has_key(3)
            False

        """
        return k in self

    def __contains__(self, k):
        """
        Containment in the set of keys.

        TESTS::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: 3 in D     # indirect doctest
            True

        As usual, keys are compared by equality and not by identity::

            sage: int(3) in D
            True

        This is a weak value dictionary. Hence, the existence of the
        dictionary does not prevent the values from garbage collection,
        thereby removing the corresponding key-value pairs::

            sage: del L[3]
            sage: 3 in D
            False

        """
        cdef list bucket
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            return False
        bucket = <list>buckref
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i) == k:
                v = <object>PyList_GetItem(bucket,i+1)
                if <Py_ssize_t><void*>v in self._pending_removals:
                    return False
                if PyWeakref_GetObject(v)==Py_None:
                    return False
                else:
                    return True
        return False

    def __len__(self):
        """
        TESTS::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: len(D)
            10
            sage: del D[2]
            sage: len(D)
            9
            sage: del L[4]
            sage: len(D)
            8
            sage: D[1] = ZZ
            sage: len(D)
            8
            sage: D[4] = ZZ
            sage: len(D)
            9
            sage: len(D) == len(D.items()) == len(D.keys()) == len(D.values())
            True

        """
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t length = 0
        cdef list bucket
        with self._iteration_context:
            while PyDict_Next(self, &pos, &basekey, &bucketref):
                bucket = <list>bucketref
                length += PyList_Size(bucket)
        return length//2-len(self._pending_removals)

    def iterkeys(self):
        """
        Iterate over the keys of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: list(sorted(D.iterkeys()))
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        with self._iteration_context:
            while PyDict_Next(self, &pos, &basekey, &bucketref):
                bucket = <list>bucketref
                l = len(bucket)
                for i from 0 <= i < l by 2:
                    v = <object>PyList_GetItem(bucket,i+1)
                    if <Py_ssize_t><void*>v in self._pending_removals:
                        continue
                    if PyWeakref_GetObject(v)!=Py_None:
                        yield <object>PyList_GetItem(bucket,i)

    def __iter__(self):
        """
        Iterate over the keys of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: sorted(list(D))    # indirect doctest
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        return self.iterkeys()

    def keys(self):
        """
        The list of keys.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object): pass
            sage: L = [Vals() for _ in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))
            sage: del L[4]

        One item got deleted from the list ``L`` and hence the corresponding
        item in the dictionary got deleted as well. Therefore, the
        corresponding key 4 is missing in the list of keys::

            sage: sorted(D.keys())
            [0, 1, 2, 3, 5, 6, 7, 8, 9]

        """
        return list(self.iterkeys())

    def itervalues(self):
        """
        Iterate over the values of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [Vals(n) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))

        We delete one item from ``D`` and we delete one item from the list
        ``L``. The latter implies that the corresponding item from ``D`` gets
        deleted as well. Hence, there remain eight values::

            sage: del D[2]
            sage: del L[5]
            sage: for v in sorted(D.itervalues()):
            ....:     print v
            <0>
            <1>
            <3>
            <4>
            <6>
            <7>
            <8>
            <9>

        """
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        with self._iteration_context:
            while PyDict_Next(self, &pos, &basekey, &bucketref):
                bucket = <list>bucketref
                l = len(bucket)
                for i from 0 < i <= l by 2:
                    v = <object>PyList_GetItem(bucket,i)
                    if <Py_ssize_t><void*>v in self._pending_removals:
                        continue
                    obj = PyWeakref_GetObject(v)
                    if obj != Py_None:
                        yield <object>obj

    def values(self):
        """
        Return the list of values.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [Vals(n) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(L))

        We delete one item from ``D`` and we delete one item from the list
        ``L``. The latter implies that the corresponding item from ``D`` gets
        deleted as well. Hence, there remain eight values::

            sage: del D[2]
            sage: del L[5]
            sage: sorted(D.values())
            [<0>, <1>, <3>, <4>, <6>, <7>, <8>, <9>]

        """
        return list(self.itervalues())

    def iteritems(self):
        """
        Iterate over the items of this dictionary.

        .. WARNING::

            Iteration is unsafe, if the length of the dictionary changes
            during the iteration! This can also happen by garbage collection.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: class Keys(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __hash__(self):
            ....:         if self.n%2:
            ....:             return 5
            ....:         return 3
            ....:     def __repr__(self):
            ....:         return "[%s]"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [(Keys(n), Vals(n)) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)

        We remove one dictionary item directly. Another item is removed by
        means of garbage collection. By consequence, there remain eight
        items in the dictionary::

            sage: del D[Keys(2)]
            sage: del L[5]
            sage: for k,v in sorted(D.iteritems()):
            ....:     print k, v
            [0] <0>
            [1] <1>
            [3] <3>
            [4] <4>
            [6] <6>
            [7] <7>
            [8] <8>
            [9] <9>

        """
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        with self._iteration_context:
            while PyDict_Next(self, &pos, &basekey, &bucketref):
                bucket = <list>bucketref
                l = len(bucket)
                for i from 0 <= i < l by 2:
                    v = <object>PyList_GetItem(bucket,i+1)
                    if <Py_ssize_t><void*>v in self._pending_removals:
                        continue
                    obj = PyWeakref_GetObject(v)
                    if obj != Py_None:
                        yield <object>PyList_GetItem(bucket,i), <object>obj

    def items(self):
        """
        The key-value pairs of this dictionary.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals:
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: class Keys(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __hash__(self):
            ....:         if self.n%2:
            ....:             return 5
            ....:         return 3
            ....:     def __repr__(self):
            ....:         return "[%s]"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [(Keys(n), Vals(n)) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)

        We remove one dictionary item directly. Another item is removed by
        means of garbage collection. By consequence, there remain eight
        items in the dictionary::

            sage: del D[Keys(2)]
            sage: del L[5]
            sage: sorted(D.items())
            [([0], <0>),
             ([1], <1>),
             ([3], <3>),
             ([4], <4>),
             ([6], <6>),
             ([7], <7>),
             ([8], <8>),
             ([9], <9>)]

        """
        return list(self.iteritems())

    def _buckets_(self):
        """
        Returns the internal data structure (DEBUGGING ONLY).

        NOTE:

        We have a list for each occuring hash value of the dictionary
        keys. Each list alternatingly comprises the key and a weak reference
        to the corresponding value.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Vals(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __repr__(self):
            ....:         return "<%s>"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: class Keys(object):
            ....:     def __init__(self, n):
            ....:         self.n = n
            ....:     def __hash__(self):
            ....:         if self.n%2:
            ....:             return 5
            ....:         return 3
            ....:     def __repr__(self):
            ....:         return "[%s]"%self.n
            ....:     def __cmp__(self, other):
            ....:         c = cmp(type(self),type(other))
            ....:         if c:
            ....:             return c
            ....:         return cmp(self.n,other.n)
            sage: L = [(Keys(n), Vals(n)) for n in range(10)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(L)
            sage: sorted(D._buckets_())
            [(3,
              [[0],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [2],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [4],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [6],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [8],
               <weakref at 0x...; to 'Vals' at 0x...>]),
             (5,
              [[1],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [3],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [5],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [7],
               <weakref at 0x...; to 'Vals' at 0x...>,
               [9],
               <weakref at 0x...; to 'Vals' at 0x...>])]
            sage: sorted(D.items())
            [([0], <0>),
             ([1], <1>),
             ([2], <2>),
             ([3], <3>),
             ([4], <4>),
             ([5], <5>),
             ([6], <6>),
             ([7], <7>),
             ([8], <8>),
             ([9], <9>)]

        """
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        cdef list buckets = []
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            buckets.append((<object>basekey, bucket))
        return buckets

cdef class _IterationContext:
    """
    An iterator that protects :class:`WeakValueDictionary` from negative
    effects of deletions that occur during iteration.

    EXAMPLES::

        sage: import sage.misc.weak_dict, gc
        sage: class Cycle:
        ....:     def __init__(self):
        ....:         self.selfref = self
        sage: C = [Cycle() for n in range(5)]
        sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(C))
        sage: remaining_keys = set(range(5))
        sage: for k in D:
        ....:     print "len",len(D)
        ....:     remaining_keys.remove(k)
        ....:     if len(D)==5:
        ....:         C[k] = None
        ....:         C[remaining_keys.pop()] = None
        ....:         _ = gc.collect()
        ....:     elif len(D)==3:
        ....:         del D[remaining_keys.pop()]
        ....:     elif len(D)==2:
        ....:         l,v = D.popitem()      # indirect doctest
        len 5
        len 3
        len 2

    This is implemented by putting all items that are deleted during iteration
    are not actually deleted, but only *marked* for deletion. Note, however,
    that they behave like items that actually *are* deleted. E.g., it is not
    possible to delete the same item twice::

        sage: class Val: pass
        sage: V = [Val() for n in range(3)]
        sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(V))
        sage: for k,v in D.iteritems():
        ....:     print k in D.keys()
        ....:     del D[k]
        ....:     print k in D.keys()
        ....:     try:
        ....:         del D[k]
        ....:     except KeyError:
        ....:         print "due error"
        True
        False
        due error
        True
        False
        due error
        True
        False
        due error

    """
    cdef WeakValueDictionary Dict
    def __init__(self, Dict):
        """
        INPUT:

        A :class:`WeakValueDictionary`

        This context manager is used during iteration over the given
        dictionary, and prevents negative side-effects of item deletions
        during iteration.

        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: class Val: pass
            sage: V = [Val() for n in range(3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(V))
            sage: for k,v in D.iteritems():  # indirect doctest
            ....:     k in D.keys()
            ....:     del D[k]
            ....:     k in D.keys()
            True
            False
            True
            False
            True
            False

        """
        self.Dict = Dict

    def __enter__(self):
        """
        Make sure that items of a weak value dictionary are not actually
        deleted, but only *marked* for deletion.

        TESTS::

            sage: import sage.misc.weak_dict
            sage: class Val: pass
            sage: V = [Val() for n in range(3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(V))
            sage: for k,v in D.iteritems():        # indirect doctest
            ....:     k in D.keys()
            ....:     del D[k]
            ....:     k in D.keys()
            True
            False
            True
            False
            True
            False

        """
        self.Dict._guard_level += 1
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Make sure that all items of a weak value dictionary that are marked
        for deletion are actually deleted, as soon as there is no iteration
        over the dictionary.

        TESTS::

            sage: import sage.misc.weak_dict
            sage: class Val: pass
            sage: V = [Val() for n in range(3)]
            sage: D = sage.misc.weak_dict.WeakValueDictionary(enumerate(V))
            sage: for k,v in D.iteritems():      # indirect doctest
            ....:     k in D.keys()
            ....:     del D[k]
            ....:     k in D.keys()
            True
            False
            True
            False
            True
            False

        """
        # Propagate errors by returning "False"
        self.Dict._guard_level -= 1
        if not self.Dict._guard_level:
            while self.Dict._pending_removals:
                self.Dict.callback(self.Dict._pending_removals.popitem()[1])
        return False
