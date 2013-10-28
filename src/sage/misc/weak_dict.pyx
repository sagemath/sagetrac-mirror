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
    ....: 
    sage: ValList = [Vals() for _ in range(10)]
    sage: D = weakref.WeakValueDictionary()
    sage: for v in ValList:
    ....:     D[Keys(v)] = v
    ....:     
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
    ....:     
    sage: len(D)
    10
    sage: del ValList
    sage: len(D)
    1
    sage: del v
    sage: len(D)
    0

Note that Sage's weak value dictionary is actually an instance of
:class:`dict`, in contrast to :mod:`weakref`'s weak value dictionary::

    sage: issubclass(weakref.WeakValueDictionary, dict)
    False
    sage: issubclass(sage.misc.weak_dict.WeakValueDictionary, dict)
    True

In addition, Sage's implementation has a better performance.

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
import sage
from sage.all import add

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

    """
    def __init__(self, data=()):
        """
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
        for k,v in data:
            self[k] = v

    def __repr__(self):
        """
        EXAMPLES::

            sage: import sage.misc.weak_dict
            sage: sage.misc.weak_dict.WeakValueDictionary()
            <WeakValueDictionary at 0x...>

        """
        return "<WeakValueDictionary at 0x%x>" % id(self)

    def callback(self, r):
        """
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
            ....: 
            sage: ValList = [Vals() for _ in range(10)]
            sage: import sage.misc.weak_dict
            sage: D = sage.misc.weak_dict.WeakValueDictionary()
            sage: for v in ValList:
            ....:     D[Keys(v)] = v
            ....:     
            sage: len(D)
            10
            sage: del ValList
            sage: len(D)
            1
            sage: del v
            sage: len(D)
            0

        """
        cdef int hashk = r.key
        cdef Py_ssize_t idr = <Py_ssize_t><void *>r
        cdef void * buckref = PyDict_GetItem(self, hashk)
        if buckref==NULL:
            return
        cdef list bucket = <list>buckref
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <Py_ssize_t>PyList_GetItem(bucket,i+1) == idr:
                del bucket[i:i+2]
                if not bucket:
                    PyDict_DelItem(self,hashk)
                return

    def setdefault(self, k, default=None):
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
        bucket.extend([k, KeyedRef(default, self.callback, hashk)])
        return default

    def __setitem__(self, k, v):
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
        bucket.extend((k, KeyedRef(v, self.callback, hashk)))

    def __delitem__(self, k):
        cdef int hashk = hash(k)
        cdef void * buckref = PyDict_GetItem(self, hashk)
        if buckref==NULL:
            raise KeyError(k)
        cdef list bucket = <list>buckref
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i) == k:
                del bucket[i:i+2]
                if not bucket:
                    PyDict_DelItem(self,hashk)
                return

    def pop(self, k):
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
                bucket.pop(i)
                out = PyWeakref_GetObject(bucket.pop(i))
                if not bucket:
                    PyDict_DelItem(self,hashk)
                if out != Py_None:
                    return <object>out
                break
        raise KeyError(k)
        
    def popitem(self):
        cdef Py_ssize_t hashk
        cdef list bucket = []
        try:
            while not bucket:
                hashk, bucket = dict.popitem(self)
        except KeyError:
            raise KeyError('popitem(): weak value dictionary is empty')
        PyDict_SetItem(self, hashk, bucket)
        k = bucket.pop(0)
        v = bucket.pop(0)
        if not bucket:
            PyDict_DelItem(self,hashk)
        cdef PyObject* out = PyWeakref_GetObject(v)
        if out!=Py_None:
            return k, <object>out
        return self.popitem()

    def get(self, k, default=None):
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            return default
        cdef list bucket = <list>buckref
        cdef PyObject* out
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i)==k:
                out = PyWeakref_GetObject(<object>PyList_GetItem(bucket,i+1))
                if out!=Py_None:
                    return <object>out
                break
        return default

    def __getitem__(self, k):
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            raise KeyError(k)
        cdef list bucket = <list>buckref
        cdef PyObject* out
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i)==k:
                out = PyWeakref_GetObject(<object>PyList_GetItem(bucket,i+1))
                if out!=Py_None:
                    return <object>out
                break
        raise KeyError(k)

    def has_key(self, k):
        return k in self

    def __contains__(self, k):
        cdef list bucket
        cdef void * buckref = PyDict_GetItem(self, hash(k))
        if buckref==NULL:
            return False
        bucket = <list>buckref
        cdef Py_ssize_t i,l
        l = len(bucket)
        for i from 0 <= i < l by 2:
            if <object>PyList_GetItem(bucket,i) == k:
                if PyWeakref_GetObject(<object>PyList_GetItem(bucket,i+1))==Py_None:
                    return False
                else:
                    return True
        return False

    def __len__(self):
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t length = 0
        cdef list bucket
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            length += PyList_Size(bucket)
        return length//2

    def iterkeys(self):
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            l = len(bucket)
            for i from 0 <= i < l by 2:
                if PyWeakref_GetObject(<object>PyList_GetItem(bucket,i+1))!=Py_None:
                    yield <object>PyList_GetItem(bucket,i)

    def __iter__(self):
        return self.iterkeys()

    def keys(self):
        return list(self.iterkeys())

    def itervalues(self):
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            l = len(bucket)
            for i from 0 < i <= l by 2:
                obj = PyWeakref_GetObject(<object>PyList_GetItem(bucket,i))
                if obj != Py_None:
                    yield <object>obj

    def values(self):
        return list(self.itervalues())

    def iteritems(self):
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            l = len(bucket)
            for i from 0 <= i < l by 2:
                obj = PyWeakref_GetObject(<object>PyList_GetItem(bucket,i+1))
                if obj != Py_None:
                    yield <object>PyList_GetItem(bucket,i), <object>obj

    def buckets(self):
        cdef PyObject * obj
        cdef list bucket
        cdef PyObject *basekey, *bucketref
        cdef Py_ssize_t pos = 0
        cdef Py_ssize_t i,l
        cdef list buckets = []
        while PyDict_Next(self, &pos, &basekey, &bucketref):
            bucket = <list>bucketref
            buckets.append(bucket)
        return buckets

    def items(self):
        return list(self.iteritems())
