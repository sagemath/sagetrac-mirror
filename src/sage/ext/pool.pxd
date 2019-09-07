cimport cython
from cpython.object cimport (PyObject, PyTypeObject,
        newfunc, destructor, Py_TPFLAGS_HAVE_GC,
        Py_TPFLAGS_HEAPTYPE, Py_TPFLAGS_HAVE_VERSION_TAG)
from cpython.type cimport type


cdef extern from *:
    void PyType_Modified(type)


cdef class ObjectPool:
    cdef readonly type cls
    cdef newfunc real_tp_new
    cdef destructor real_tp_dealloc

    cdef readonly Py_ssize_t size
    cdef readonly Py_ssize_t used
    cdef PyObject** pool

    cdef tp_new(self, cls, args, kwds)
    cdef void tp_dealloc(self, PyObject* obj)

    cdef inline new(self):
        """
        Return a newly object, using the pool if possible.
        This is equivalent to ``(self.cls).__new__(self.cls)``.
        """
        return self.tp_new(self.cls, <tuple>NULL, <dict>NULL)

    cdef inline PyTypeObject* tp(self):
        return <PyTypeObject*>self.cls

    cdef inline bint is_enabled(self):
        return self.tp().tp_new is not self.real_tp_new

    cdef int pool_enable(self) except -1
    cdef void pool_in(self, obj)
    cdef int pool_out(self, obj) except -1
    cdef pool_new(self)
    cdef void pool_dealloc(self, obj)

    cdef inline long enable(self) except -1:
        """
        Changes ``self.cls`` to use this pool.
        """
        # Implementation detail: this must be an inline function because
        # we must call _attach_pool() in the Cython module where the
        # type self.cls is implemented, not in the module where the
        # ObjectPool is implemented.
        #
        # Because of this, enable() cannot be changed in subclasses.
        # To overcome this, we provide a hook pool_enable() which can
        # be customized.
        cdef long n = _attach_pool(self.cls, self)
        self.pool_enable()
        PyType_Modified(self.cls)
        return n


cdef class ObjectPool_NoGC(ObjectPool):
    pass


cdef class ObjectPool_Cloning(ObjectPool_NoGC):
    cdef PyObject* template

    cdef int clone_out(self, obj) except -1
    cdef void clone_in(self, obj)


cdef object (*ObjectPool_new "ObjectPool_new")(ObjectPool, cls, args, kwds)
cdef void (*ObjectPool_dealloc "ObjectPool_dealloc")(ObjectPool, PyObject* self)


cdef extern from "sage/ext/pool_impl.h":
    long _attach_pool(type, ObjectPool)

    # Decrease the refcount of an object without deallocating it,
    # even if the refcount would become zero.
    void Py_DECREF_KEEP(PyObject*)
