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

TESTS:

See :trac:`13394` for a discussion. ::

    sage: from sage.cpython.dict_del_by_value import del_dictitem_by_exact_value
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
    ....:         del_dictitem_by_exact_value(D2,vo,hash(ko))
    ....:     assert D1 == D2

No action is taken if the item prescribed by key hash and value does not
exist in the dictionary::

    sage: D = {1: ZZ}
    sage: del_dictitem_by_exact_value(D, ZZ, 2)
    sage: D
    {1: Integer Ring}
    sage: del_dictitem_by_exact_value(D, QQ, 1)
    sage: D
    {1: Integer Ring}

The following is an indirect doctest, as discussed on :trac:`13394`. ::

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

The following shows that the deletion of deeply nested structures does not
result in an error, by :trac:`15506`::

    sage: class A: pass
    sage: a = A(); prev = a
    sage: M = WeakValueDictionary()
    sage: for i in range(10^3+10): newA = A(); M[newA] = prev; prev = newA
    sage: del a

AUTHORS:

- Nils Bruin (2017-05)
"""

#*****************************************************************************
#       Copyright (C) 2017 Nils Bruin <nbruin@sfu.ca>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function, absolute_import

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t
from cpython.list cimport PyList_New
from cpython.ref cimport PyObject, Py_INCREF
from cpython.dict cimport PyDictEntry
from cpython.tuple cimport PyTuple_New

cdef extern from "Python.h":
    #we need this redefinition because we want to be able to call
    #PyWeakref_GetObject with borrowed references. This is the recommended
    #strategy according to Cython/Includes/cpython/__init__.pxd
    PyObject* PyWeakref_GetObject(PyObject * wr)
    int PyList_SetItem(object list, Py_ssize_t index, PyObject * item) except -1
    void PyTuple_SET_ITEM(object tuple, Py_ssize_t index, PyObject* item)


IF PY_VERSION_HEX < 0x03060000:
    # This routine extracts the "dummy" sentinel value that is used in
    # dicts to mark "freed" slots. We need that to delete things
    # ourselves.
    cdef get_dummy():

        cdef dict D = {0: 0}
        del D[0] #ensure that there is a "deleted" entry in the dict
        cdef size_t mask = D.ma_mask
        cdef size_t i
        # Since our entry had hash 0, we should really succeed on the
        # first iteration
        for i in range(mask+1):
            ep = &(D.ma_table[i])
            if ep.me_key is not NULL:
                return <object>ep.me_key
        raise RuntimeError("Problem initializing dummy")

    cdef object dummy = get_dummy()

    cpdef del_dictitem_by_exact_value(dict mp, value, Py_hash_t hash):
        """
        Delete an item in a dict with a given value. The value is looked
        up by identity (using the ``is`` operator).

        INPUT:

        - ``mp`` -- dict
        - ``value`` -- a value of the dictionary
        - ``hash`` -- hash of the key by which the value is stored in the dict

        OUTPUT:

        - ``None`` if the value was not found
        - The deleted ``(key, value)`` item otherwise

        The hash bucket determined by the given hash is searched for the
        item containing the given value (by identity). Of course, this
        only makes sense if the pairs ``(h, value)`` corresponding to
        items in ``D`` are pair-wise distinct.

        If a matching item can not be found, the function does nothing
        and silently returns. Otherwise, the item is removed from the
        dict.

        This is used in callbacks for the weak values of
        :class:`WeakValueDictionary`.
        """
        cdef size_t i
        cdef size_t perturb
        cdef size_t mask = <size_t> mp.ma_mask
        cdef PyDictEntry* ep0 = mp.ma_table
        cdef PyDictEntry* ep
        i = hash & mask
        ep = &(ep0[i])

        if ep.me_key is NULL:
            # key not found
            return None

        perturb = hash
        while (ep.me_value is not <PyObject*>value or ep.me_hash != hash):
            i = (i << 2) + i + perturb +1
            ep = &ep0[i & mask]
            if ep.me_key is NULL:
                # key not found
                return None
            perturb = perturb >> 5 #this is the value of PERTURB_SHIFT

        T = PyTuple_New(2)
        PyTuple_SET_ITEM(T, 0, ep.me_key)
        PyTuple_SET_ITEM(T, 1, ep.me_value)
        Py_INCREF(dummy)
        ep.me_key = <PyObject*>dummy
        ep.me_value = NULL
        mp.ma_used -= 1
        # We have transferred the to-be-deleted references to the tuple
        # T. We now return this tuple so that it will be deallocated by
        # the caller. Since this is a tuple, the deallocation routine
        # uses the Python Trashcan macros to avoid stack overflow in
        # deleting deep structures.
        return T

ELSE:  # Python >= 3.6

    ####
    #definitions replicated from CPython's Objects/dict-common.h
    #(this file is not exported from CPython, so we need to be
    #careful the definitions are in step with what happens there.
    ctypedef union IndexBlock:
        int8_t as_1[8]
        int16_t as_2[4]
        int32_t as_4[2]
        int64_t as_8[1]

    ctypedef struct PyDictKeysObject:
        Py_ssize_t dk_refcnt
        Py_ssize_t dk_size
        void* dk_lookup
        Py_ssize_t dk_usable
        Py_ssize_t dk_nentries
        IndexBlock dk_indices

    ctypedef struct PyDictKeyEntry:
        Py_hash_t me_hash
        PyObject* me_key
        PyObject* me_value

    cdef Py_ssize_t DKIX_EMPTY = -1
    cdef Py_ssize_t DKIX_DUMMY = -2
    cdef Py_ssize_t DKIX_ERROR = -3

    #####
    #These routines are copied from CPython's Object/dictobject.c
    #in order to access PyDictKeysObject fields

    cdef inline int DK_IXSIZE(PyDictKeysObject *keys):
        cdef Py_ssize_t s = keys.dk_size
        if s <= 0xff:
            return 1
        elif s <= 0xffff:
            return 2
        elif s <= 0xffffffff:
            return 4
        else:
            return 8

    cdef inline PyDictKeyEntry * DK_ENTRIES(PyDictKeysObject *keys):
        return <PyDictKeyEntry*> &(keys.dk_indices.as_1[keys.dk_size * DK_IXSIZE(keys)])

    cdef inline Py_ssize_t dk_get_index(PyDictKeysObject *keys, Py_ssize_t i):
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

    cdef inline dk_set_index(PyDictKeysObject *keys, Py_ssize_t i, Py_ssize_t ix):
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

    cdef void* DK_LOOKUP(dict mp):
        return <void *>((<PyDictKeysObject *>(mp.ma_keys)).dk_lookup)

    cdef void* get_lookdict():
        cdef dict D = {}
        # This should trigger the initialization of the general
        # lookup function on the dict.
        PyDict_GetItemWithError(D, None)
        return DK_LOOKUP(D)

    cdef void* lookdict = get_lookdict()

    cdef int ensure_allows_deletions(dict mp) except -1:
        if DK_LOOKUP(mp) is not lookdict:
            # On normal dictionaries (non-split table), looking up a key
            # that is not a unicode object triggers installation of the
            # general lookup function (which can deal with DKIX_DUMMY)
            PyDict_GetItemWithError(mp, None)
            # this can actually fail if mp is a dictionary with split table
            if DK_LOOKUP(mp) is not lookdict:
                raise TypeError("cannot deal with split table dictionaries")

    cpdef del_dictitem_by_exact_value(dict mp, value, Py_hash_t hash):
        """
        Delete an item in a dict with a given value. The value is looked
        up by identity (using the ``is`` operator).

        INPUT:

        - ``mp`` -- dict
        - ``value`` -- a value of the dictionary
        - ``hash`` -- hash of the key by which the value is stored in the dict

        OUTPUT:

        - ``None`` if the value was not found
        - The deleted ``(key, value)`` item otherwise

        The hash bucket determined by the given hash is searched for the
        item containing the given value (by identity). Of course, this
        only makes sense if the pairs ``(h, value)`` corresponding to
        items in ``D`` are pair-wise distinct.

        If a matching item can not be found, the function does nothing
        and silently returns. Otherwise, the item is removed from the
        dict.

        This is used in callbacks for the weak values of
        :class:`WeakValueDictionary`.
        """
        if mp.ma_values is not NULL:
            raise NotImplementedError("del_dictitem_by_exact_value cannot be applied to a shared key dict")

        cdef size_t i
        cdef PyDictKeysObject * keys = <PyDictKeysObject *>(mp.ma_keys)
        cdef size_t perturb
        cdef size_t mask = <size_t> keys.dk_size-1
        cdef PyDictKeyEntry* entries = DK_ENTRIES(keys)
        cdef PyDictKeyEntry* ep

        i = <size_t>hash & mask
        ix = dk_get_index(keys, i)
        if ix == DKIX_EMPTY:
            # key not found
            return None

        ep = &(entries[ix])
        perturb = hash
        while (ep.me_value != value or ep.me_hash != hash):
            perturb = perturb >> 5 #this is the value of PERTURB_SHIFT
            i = ((i << 2) + i + perturb +1) & mask
            ix = dk_get_index(keys, i)
            if ix == DKIX_EMPTY:
                # key not found
                return None
            ep = &(entries[ix])

        ensure_allows_deletions(mp)

        T = PyTuple_New(2)
        PyTuple_SET_ITEM(T, 0, ep.me_key)
        PyTuple_SET_ITEM(T, 1, ep.me_value)
        ep.me_key = NULL
        ep.me_value = NULL
        mp.ma_used -= 1
        dk_set_index(keys, i, DKIX_DUMMY)
        # We have transferred the to-be-deleted references to the tuple
        # T. We now return this tuple so that it will be deallocated by
        # the caller. Since this is a tuple, the deallocation routine
        # uses the Python Trashcan macros to avoid stack overflow in
        # deleting deep structures.
        return T
