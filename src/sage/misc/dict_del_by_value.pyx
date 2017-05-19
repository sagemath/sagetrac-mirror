"""
Delete item from PyDict by exact value and hash

Beware that the implementation of the routine here relies on implementation
details of CPython's dict that go beyond the published API. This file depends
on python version when cythonized. It expects PY_VERSION_HEX to be available
in the cythonization and the result depends on it (and needs to match the
python version the C-compiler compiles it against). Usage should do something
along the lines of

    cythonize("dict_del_by_value.pyx",
        compile_time_env({"PY_VERSION_HEX": sys.hexversion}))

AUTHORS:

- Nils Bruin (2017-05)
"""
########################################################################
#       Copyright (C) 2017 Nils Bruin <nbruin@sfu.ca>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
########################################################################
from __future__ import print_function

import weakref
from weakref import KeyedRef
from copy import deepcopy

from cpython.list cimport PyList_New
from cpython cimport Py_XINCREF, Py_XDECREF

IF PY_VERSION_HEX<=0x02ffffff:
    cdef extern from "Python.h":
        ctypedef struct PyDictEntry:
            Py_ssize_t me_hash
            PyObject* me_key
            PyObject* me_value
        ctypedef struct PyDictObject:
            Py_ssize_t ma_fill
            Py_ssize_t ma_used
            Py_ssize_t ma_mask
            PyDictEntry* ma_table
            PyDictEntry* (*ma_lookup)(PyDictObject *mp, PyObject *key, Py_hash_t hash) except NULL

        PyObject* Py_None
        #we need this redefinition because we want to be able to call
        #PyWeakref_GetObject with borrowed references. This is the recommended
        #strategy according to Cython/Includes/cpython/__init__.pxd
        PyObject* PyWeakref_GetObject(PyObject * wr)
        int PyList_SetItem(object list, Py_ssize_t index, PyObject * item) except -1

    from cpython.object cimport PyObject_Hash
    cdef PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL:
        cdef PyDictEntry* ep
        cdef PyDictObject* mp = <PyDictObject*><void *>op
        ep = mp.ma_lookup(mp, <PyObject*><void*>key, PyObject_Hash(key))
        if ep:
            return ep.me_value
        else:
            return NULL

    #this routine extracts the "dummy" sentinel value that is used in dicts to mark
    #"freed" slots. We need that to delete things ourselves.

    cdef PyObject* init_dummy() except NULL:
        cdef dict D = dict()
        cdef PyDictObject* mp = <PyDictObject *><void *>D
        cdef size_t mask
        cdef PyDictEntry* ep0 = mp.ma_table
        cdef PyDictEntry* ep
        cdef size_t i
        global dummy

        D[0]=0; del D[0] #ensure that there is a "deleted" entry in the dict
        mask = mp.ma_mask
        #since our entry had hash 0, we should really succeed on the first iteration
        for i in range(mask+1):
            ep = &(ep0[i])
            if ep.me_key != NULL:
                return ep.me_key
        raise RuntimeError("Problem initializing dummy")

    #note that dummy here is a borrowed reference. That's not a problem because
    #we're never giving it up and dictobject.c is also holding a permanent reference
    #to this object
    cdef PyObject* dummy = init_dummy()

    #this routine looks for the first entry in dict D with given hash of the
    #key and given (identical!) value and deletes that entry.
    cdef del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash):
        """
        This is used in callbacks for the weak values of :class:`WeakValueDictionary`.

        INPUT:

        - ``PyDictObject *mp`` -- pointer to a dict
        - ``PyObject *value``  -- pointer to a value of the dictionary
        - ``Py_hash_t hash``        -- hash of the key by which the value is stored in the dict

        The hash bucket determined by the given hash is searched for the item
        containing the given value. If this item cannot be found, the function is
        silently returning. Otherwise, the item is removed from the dict.

        TESTS:

        The following is an indirect doctest, as discussed on :trac:`13394`.
        ::

            sage: from sage.misc.weak_dict import WeakValueDictionary
            sage: V = [set(range(n)) for n in range(5)]
            sage: D = WeakValueDictionary(enumerate(V))

        The line ``V[k] = None`` triggers execution of the callback functions of
        the dict values. However, the actual deletion is postponed till after the
        iteration over the dictionary has finished. Hence, when the callbacks are
        executed, the values which the callback belongs to has already been
        overridded by a new value. Therefore, the callback does not delete the
        item::

            sage: for k in D:    # indirect doctest
            ....:     V[k] = None
            ....:     D[k] = ZZ
            sage: len(D)
            5
            sage: D[1]
            Integer Ring

        TESTS:

        The following shows that the deletion of deeply nested structures does not
        result in an error, by :trac:`15506`::

            sage: class A: pass
            sage: a = A(); prev = a
            sage: M = WeakValueDictionary()
            sage: for i in range(10^3+10): newA = A(); M[newA] = prev; prev = newA
            sage: del a

        """
        cdef size_t i
        cdef size_t perturb
        cdef size_t mask = <size_t> mp.ma_mask
        cdef PyDictEntry* ep0 = mp.ma_table
        cdef PyDictEntry* ep
        i = hash & mask
        ep = &(ep0[i])

        if ep.me_key == NULL:
            # key not found
            return

        perturb = hash
        while (<PyObject *>(ep.me_value) != value or ep.me_hash != hash):
            i = (i << 2) + i + perturb +1
            ep = &ep0[i & mask]
            if ep.me_key == NULL:
                # key not found
                return
            perturb = perturb >> 5 #this is the value of PERTURB_SHIFT

        T=PyList_New(2)
        PyList_SetItem(T, 0, ep.me_key)
        if dummy == NULL:
            raise RuntimeError("dummy needs to be initialized")
        Py_XINCREF(dummy)
        ep.me_key = dummy
        PyList_SetItem(T, 1, ep.me_value)
        ep.me_value = NULL
        mp.ma_used -= 1
        #We have transferred the to-be-deleted references to the list T
        #we now delete the list so that the actual decref happens through a
        #deallocation routine that uses the Python Trashcan macros to
        #avoid stack overflow in deleting deep structures.
        del T

ELIF PY_VERSION_HEX>=0x03060000:

    from libc.stdint cimport int8_t,int16_t,int32_t,int64_t
    cdef extern from "Python.h":
        ctypedef struct PyDictObject:
            Py_ssize_t ma_used
            PyDictKeysObject * ma_keys
            PyObject ** ma_values

        PyObject* Py_None
        #we need this redefinition because we want to be able to call
        #PyWeakref_GetObject with borrowed references. This is the recommended
        #strategy according to Cython/Includes/cpython/__init__.pxd
        PyObject* PyWeakref_GetObject(PyObject * wr)
        int PyList_SetItem(object list, Py_ssize_t index, PyObject * item) except -1
        PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
        int PyWeakref_Check(PyObject * ob)

        ctypedef struct PyDictKeysObject
    ####
    #definitions replicated from CPython's Objects/dict-common.h
    #(this file is not exported fron CPython, so we need to be
    #careful the definitions are in step with what happens there.
    ctypedef union IndexBlock:
        int8_t as_1[8]
        int16_t as_2[4]
        int32_t as_4[2]
        int64_t as_8[1]

    ctypedef struct MyPyDictKeysObject:
        Py_ssize_t dk_refcnt
        Py_ssize_t dk_size
        Py_ssize_t (*dk_lookup)(PyDictObject *mp, PyObject *key, Py_hash_t hash, PyObject ** value_addr, Py_ssize_t *hashpos)
        Py_ssize_t dk_usable
        Py_ssize_t dk_nentries
        IndexBlock dk_indices

    ctypedef struct PyDictKeyEntry:
        Py_hash_t me_hash
        PyObject * me_key
        PyObject * me_value

    cdef Py_ssize_t DKIX_EMPTY = -1
    cdef Py_ssize_t DKIX_DUMMY = -2
    cdef Py_ssize_t DKIX_ERROR = -3

    #####
    #These routines are copied from CPython's Object/dictobject.c
    #in order to access PyDictKeysObject fields

    cdef inline int DK_IXSIZE(MyPyDictKeysObject *keys):
        cdef Py_ssize_t s = keys.dk_size
        cdef int i
        if s <= 0xff:
            i = 1
        elif s <= 0xffff:
            i = 2
        elif s <= 0xffffffff:
            i = 4
        else:
            i = 8
        return i
        
    cdef inline PyDictKeyEntry * DK_ENTRIES(MyPyDictKeysObject *keys):
        return <PyDictKeyEntry*> &(keys.dk_indices.as_1[keys.dk_size * DK_IXSIZE(keys)])

    cdef inline Py_ssize_t dk_get_index(MyPyDictKeysObject *keys, Py_ssize_t i):
        cdef Py_ssize_t s = keys.dk_size
        cdef Py_ssize_t ix
        if s <= 0xff:
            ix = keys.dk_indices.as_1[i]
        elif s <= 0xffff:
            ix = keys.dk_indices.as_2[i]
        elif s <= 0xffffffff:
            ix = keys.dk_indices.as_4[i]
        else:
            ix = keys.dk_indices.as_8[i]
        return ix

    cdef inline dk_set_index(MyPyDictKeysObject *keys, Py_ssize_t i, Py_ssize_t ix):
        cdef Py_ssize_t s = keys.dk_size
        if s <= 0xff:
            keys.dk_indices.as_1[i] = ix
        elif s <= 0xffff:
            keys.dk_indices.as_2[i] = ix
        elif s <= 0xffffffff:
            keys.dk_indices.as_4[i] = ix
        else:
            keys.dk_indices.as_8[i] = ix
        return ix
    #End of replication of Object/dictobject.c
    ######

    cdef void * lookdict

    cdef void * DK_LOOKUP(PyDictObject *mp):
        return <void *>((<MyPyDictKeysObject *>(mp.ma_keys)).dk_lookup)

    def init_lookdict():
        global lookdict
        cdef dict D={}
        #this should trigger the initialization of the general
        #lookup function on the dict.
        PyDict_GetItemWithError(D,None)
        #which we store in a global variable
        lookdict=DK_LOOKUP(<PyDictObject *>D)
        
    init_lookdict()

    cdef ensure_allows_deletions(PyDictObject *mp):
        if DK_LOOKUP(mp) != lookdict:
            #on normal dictionaries (non-split table), looking up a key
            #that is not a unicode object triggers installation of the
            #general lookup function (which can deal with DKIX_DUMMY)
            PyDict_GetItemWithError(<dict>mp,None)
            #this can actually fail if mp is a dictionary with split table
            assert DK_LOOKUP(mp) == lookdict
                
    cdef del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash):
        """
        This is used in callbacks for the weak values of :class:`WeakValueDictionary`.

        INPUT:

        - ``PyDictObject *mp`` -- pointer to a dict
        - ``PyObject *value``  -- pointer to a value of the dictionary
        - ``Py_hash_t hash``        -- hash of the key by which the value is stored in the dict

        The hash bucket determined by the given hash is searched for the item
        containing the given value. If this item cannot be found, the function is
        silently returning. Otherwise, the item is removed from the dict.

        TESTS:

        The following is an indirect doctest, as discussed on :trac:`13394`.
        ::

            sage: from sage.misc.weak_dict import WeakValueDictionary
            sage: V = [set(range(n)) for n in range(5)]
            sage: D = WeakValueDictionary(enumerate(V))

        The line ``V[k] = None`` triggers execution of the callback functions of
        the dict values. However, the actual deletion is postponed till after the
        iteration over the dictionary has finished. Hence, when the callbacks are
        executed, the values which the callback belongs to has already been
        overridded by a new value. Therefore, the callback does not delete the
        item::

            sage: for k in D:    # indirect doctest
            ....:     V[k] = None
            ....:     D[k] = ZZ
            sage: len(D)
            5
            sage: D[1]
            Integer Ring

        TESTS:

        The following shows that the deletion of deeply nested structures does not
        result in an error, by :trac:`15506`::

            sage: class A: pass
            sage: a = A(); prev = a
            sage: M = WeakValueDictionary()
            sage: for i in range(10^3+10): newA = A(); M[newA] = prev; prev = newA
            sage: del a

        """
        cdef size_t i
        cdef MyPyDictKeysObject * keys = <MyPyDictKeysObject *>(mp.ma_keys)
        cdef size_t perturb
        cdef size_t mask = <size_t> keys.dk_size-1
        cdef PyDictKeyEntry *entries, *ep
        entries = DK_ENTRIES(keys)
        
        if mp.ma_values != NULL:
            print ("del_dictitem_by_exact_value cannot be applied to a shared key dict")
            return
        
        i = <size_t>hash & mask
        ix = dk_get_index(keys, i)

        if ix == DKIX_EMPTY:
            # key not found
            return

        ep = &(entries[ix])
        perturb = hash
        while (ep.me_value != value or ep.me_hash != hash):
            perturb = perturb >> 5 #this is the value of PERTURB_SHIFT
            i = ((i << 2) + i + perturb +1) & mask
            ix = dk_get_index(keys, i)
            if ix == DKIX_EMPTY:
                # key not found
                return
            ep = &(entries[ix])

        ensure_allows_deletions(mp)
        T=PyList_New(2)
        PyList_SetItem(T, 0, ep.me_key)
        PyList_SetItem(T, 1, ep.me_value)
        ep.me_key = NULL
        ep.me_value = NULL
        mp.ma_used -= 1
        dk_set_index(keys, i, DKIX_DUMMY)
        #We have transferred the to-be-deleted references to the list T
        #we now delete the list so that the actual decref happens through a
        #deallocation routine that uses the Python Trashcan macros to
        #avoid stack overflow in deleting deep structures.
        del T

#END IF PY_VERSION_HEX

def test_del_dictitem_by_exact_value(D, value, h):
    """
    This function helps testing some cdef function used to delete dictionary items.

    INPUT:

    - ``D`` -- a Python ``<dict>``.
    - ``value`` -- an object that is value ``D``.
    - ``h`` -- the hash of the key under which to find ``value`` in ``D``.

    The underlying cdef function deletes an item from ``D`` that is in the
    hash bucket determined by ``h`` and whose value is identic with
    ``value``. Of course, this only makes sense if the pairs ``(h, value)``
    corresponding to items in ``D`` are pair-wise distinct.

    If a matching item can not be found, the function does nothing and
    silently returns.

    TESTS:

    See :trac:`13394` for a discussion.
    ::

        sage: from sage.misc.dict_del_by_value import test_del_dictitem_by_exact_value
        sage: B=1000
        sage: L=list(range(B))
        sage: D1=dict()
        sage: D2=dict()
        sage: for i in range(100000):        # long time
        ....:     ki=L[floor(random()*B)]
        ....:     vi=L[floor(random()*B)]
        ....:     D1[ki]=vi
        ....:     D2[ki]=vi
        ....:     ko=L[floor(random()*B)]
        ....:     if ko in D1:
        ....:         vo=D1[ko]
        ....:         del D1[ko]
        ....:         test_del_dictitem_by_exact_value(D2,vo,hash(ko))
        ....:     assert D1 == D2

    No action is taken if the item prescribed by key hash and value does not
    exist in the dictionary::

        sage: D = {1: ZZ}
        sage: test_del_dictitem_by_exact_value(D, ZZ, 2)
        sage: D
        {1: Integer Ring}
        sage: test_del_dictitem_by_exact_value(D, QQ, 1)
        sage: D
        {1: Integer Ring}

    """
    return del_dictitem_by_exact_value(<PyDictObject *>D, <PyObject *>value, h)
