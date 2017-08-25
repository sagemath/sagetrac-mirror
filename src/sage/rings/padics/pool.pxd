from cpython.object cimport PyObject, PyTypeObject, destructor

cdef class Pool:
    cdef object __weakref__

    cdef PyObject** elements
    cdef bint enabled
    cdef long size
    cdef long allocated

    cdef PyTypeObject* type
    cdef destructor save_tp_dealloc

    cdef Pool _new(self)


cdef PY_NEW_WITH_POOL(type t, Pool pool)

cdef get_pool(t, size)
  
