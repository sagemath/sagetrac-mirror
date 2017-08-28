import _weakref

from cpython.object cimport PyObject, PyTypeObject, destructor
from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from sage.rings.padics.local_generic cimport LocalGeneric
from sage.rings.padics.local_generic_element cimport LocalGenericElement

from sage.rings.integer_ring import ZZ


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

    cdef Pool new_enabled(self):
        cdef Pool pool = Pool(<type>self.type)
        pool.save_tp_dealloc = self.save_tp_dealloc
        pool.enabled = 1
        return pool

    def automatic_resize(self, enable=True):
        if enable:
            if self.enabled:
                self.type.tp_dealloc = &tp_dealloc_with_resize
            else:
                raise ValueError("this pool is disabled")
        else:
            self.type.tp_dealloc = &tp_dealloc

    def length(self):
        return ZZ(self.size)

    def usage(self):
        return ZZ(self.allocated)

    def maximal_usage(self):
        return ZZ(self.max_allocated)

    def resize(self, size=None):
        if not self.enabled:
            raise ValueError("cannot resize a pool which is not enabled")
        if size is None:
            s = 2 * self.size 
        else:
            s = <long?>size
        self.clear(s)
        self.size = s
        self.elements = <PyObject**> sig_realloc(self.elements, s*sizeof(PyObject*))    # should we test self.elements == NULL?        

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
        #print("reuse")
        pool.allocated -= 1
        o = pool.elements[pool.allocated]
        o.ob_refcnt = 0
        return <object>o
    else:
        #print("create")
        return (<PyTypeObject*>t).tp_new(t, <object>NULL, <object>NULL)


cdef void tp_dealloc(PyObject* o):
    cdef Pool pool = (<LocalGeneric>(<LocalGenericElement>o)._parent)._pool
    if pool.allocated < pool.size:
        #print("add to pool")
        o.ob_refcnt = 1
        pool.elements[pool.allocated] = o
        pool.allocated += 1
    else:
        #print("dealloc")
        pool.save_tp_dealloc(o)

cdef void tp_dealloc_with_resize(PyObject* o):
    cdef Pool pool = (<LocalGeneric>(<LocalGenericElement>o)._parent)._pool
    if pool.allocated >= pool.size:
        pool.resize()
    #print("add to pool")
    o.ob_refcnt = 1
    pool.elements[pool.allocated] = o
    pool.allocated += 1




cdef long DEFAULT_POOL_SIZE = 100
cdef dict pools_disabled = {}
cdef dict pools_enabled = {}


cdef pool_disabled(type t):
    cdef Pool pool = None
    if pools_disabled.has_key(t):
        wr_pool = pools_disabled[t]
        pool = wr_pool()
    if pool is None:
        pool = Pool(t)
    pools_disabled[t] = _weakref.ref(pool)
    return pool


cdef pool_enabled(Pool pool_dis, long size, bint is_global):
    cdef type t = <type>pool_dis.type
    cdef Pool pool = None
    if is_global:
        if pools_enabled.has_key(t):
            wr_pool = pools_enabled[t]
            pool = wr_pool()
        if pool is None:
            pool = pool_dis.new_enabled()
            if size == 0:
                size = DEFAULT_POOL_SIZE
        pools_enabled[t] = _weakref.ref(pool)
    else:
        pool = pool_dis.new_enabled()
        if size == 0:
            size = DEFAULT_POOL_SIZE
    if size > 0:
        pool.resize(size)
    return pool
