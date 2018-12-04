cimport cython
from cpython.object cimport PyObject, PyTypeObject, traverseproc, inquiry, visitproc

cdef extern from *:
    int unlikely(int) nogil  # Defined by Cython


ctypedef unsigned int T_t


@cython.no_gc
cdef class MultiRefManager:
    # Object which is being referenced. This is not a strong reference,
    # it is only valid as long as numstrong > 0. When numstrong becomes
    # zero, this is cleared to NULL. That does not mean that the object
    # is dead, there may still be other unrelated references to it.
    cdef PyObject* object

    # Total number of strong/weak references to object.
    cdef int numstrong, numweak

    # Time when curweak is valid.
    # When T_traverse != T_cur, curweak is assumed to be 0.
    cdef T_t T_cur

    # Current number of known weak references. This is only meaningful
    # during a single traverse loop.
    cdef int curweak

    cdef MultiRef new_ref(self, bint strong)


cdef class MultiRef:
    cdef MultiRefManager manager

    # If T_weak == manager.T_cur == T_traverse, then this reference is
    # known to be weak. In all other cases, its status is undecided.
    cdef T_t T_weak

    cdef inline bint alive(self):
        M = self.manager
        if M is None:
            return False
        if unlikely(M.object is NULL):
            # If the object is dead, kill the manager too
            self.manager = None
            return False
        return True

    cdef inline MultiRefManager get_manager(self):
        if not self.alive():
            raise RuntimeError("dead MultiRef object")
        return self.manager
