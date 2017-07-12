from cpython.object cimport PyObject

cdef extern from *:
    ctypedef struct PyDictEntry:     # Py 2
        Py_ssize_t me_hash
        PyObject* me_key
        PyObject* me_value

    ctypedef struct PyDictObject:
        Py_ssize_t ma_used
        Py_ssize_t ma_fill           # Py 2
        Py_ssize_t ma_mask           # Py 2
        PyDictEntry* ma_table        # Py 2
        PyDictEntry* (*ma_lookup)(   # Py 2
                PyDictObject* mp, object key, Py_hash_t hash) except NULL
        # ma_keys is really a PyDictKeysObject* but that struct is
        # not part of the public API
        void* ma_keys                # Py 3
        PyObject ** ma_values        # Py 3


cdef del_dictitem_by_exact_value(PyDictObject *mp, PyObject *value, Py_hash_t hash)

# This function is available in Py3; we implement it ourselves in Py2.
IF PY_VERSION_HEX < 0x03000000:
    cdef PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
ELSE:
    cdef extern from "Python.h":
        PyObject* PyDict_GetItemWithError(dict op, object key) except? NULL
