import _weakref as wr

from cpython.object cimport PyObject, PyTypeObject, destructor
from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from sage.rings.padics.local_generic cimport LocalGeneric
from sage.rings.padics.local_generic_element cimport LocalGenericElement


cdef class Pool:
    def __cinit__(self, t):
        self.type = <PyTypeObject*>t
        self.save_tp_dealloc = self.type.tp_dealloc
        self.type.tp_dealloc = &tp_dealloc
        self.enabled = 0

    def __repr__(self):
        if self.enabled:
            if self.allocated > 1:
                s = "Pool of %s elements of type " % self.allocated
            else:
                s = "Pool of %s element of type " % self.allocated
        else:
            s = "Disabled pool of type "
        s += str(<type>self.type)[7:-2]
        return s

    cdef Pool _new(self):
        cdef Pool pool = Pool(<type>self.type)
        pool.save_tp_dealloc = self.save_tp_dealloc
        pool.enabled = 1
        return pool

    def resize(self, long size):
        if not self.enabled:
            raise ValueError("cannot resize a pool which is not enabled")
        self.clear(size)
        self.size = size
        self.elements = <PyObject**> sig_realloc(self.elements, size*sizeof(PyObject*))    # should we test self.elements == NULL?        

    def clear(self, start=0):
        cdef long s = <long?>start
        cdef long i
        if s < self.allocated:
            for i from s <= i < self.allocated:
                self.save_tp_dealloc(self.elements[i])
            self.allocated = s

    def __dealloc__(self):
        self.clear()
        if self.elements != NULL:
            sig_free(self.elements)


cdef inline PY_NEW_WITH_POOL(type t, Pool pool):
    cdef PyObject* o
    if pool.allocated > 0:
        pool.allocated -= 1
        o = pool.elements[pool.allocated]
        o.ob_refcnt = 0
        return <object>o
    else:
        return (<PyTypeObject*>t).tp_new(t, <object>NULL, <object>NULL)


cdef void tp_dealloc(PyObject* o):
    cdef Pool pool = (<LocalGeneric>(<LocalGenericElement>o)._parent)._pool
    if pool.allocated < pool.size:
        o.ob_refcnt = 1
        pool.elements[pool.allocated] = o
        pool.allocated += 1
    else:
        pool.save_tp_dealloc(o)




cdef long DEFAULT_POOL_SIZE = 100
cdef dict pools = {}

cdef get_pool(t, size):
    cdef Pool pool_empty = None
    cdef Pool pool = None
    if not isinstance(t, type):
        raise TypeError("t must be a type")
    if pools.has_key(t):
        wr_pool_empty, wr_pool = pools[t]
        pool_empty, pool = wr_pool_empty(), wr_pool()
    if pool_empty is None:
        if pool is not None:
            raise RuntimeError
        pool_empty = Pool(t)
    if pool is None:
        pool = pool_empty._new()
        if size is None:
            size = DEFAULT_POOL_SIZE
    if size is not None:
        pool.resize(size)
    pools[t] = wr.ref(pool_empty), wr.ref(pool)
    return pool_empty, pool
