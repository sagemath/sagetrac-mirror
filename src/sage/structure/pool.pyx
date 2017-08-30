"""
This module implements pools for Sage elements.

The purpose of a pool is to speed up creation and deletion of 
elements. This works as follows. When an element is collected
by the garbage collected, it is not deleted from memory but
pushed to the pool (except if the pool is full). On the other
hand, when a new element (with the same parent) is requested
for creation, the system checks if there are some available
elements in the pool. If there is one, it is returned and,
of course, removed from the pool. Otherwise, a new element
is created.

When elements are quite simple (e.g. small integers, p-adic 
numbers), the bottleneck is often the creation/deletion of
instances of the corresponding classes. Therefore, having a
pool can improve drastically the performances.

Parents that wants to benefit of pool features should derive
from the class :class:`ParentWithPool`.

EXAMPLES::

    sage: R = Zp(3)
    sage: x = R.random_element()
    sage: y = R.random_element()

    sage: R.pool_disable()
    sage: timeit("z = x + y")   # somewhat random
    625 loops, best of 3: 390 ns per loop

    sage: R.pool_enable()
    sage: timeit("z = x + y")   # somewhat random
    625 loops, best of 3: 193 ns per loop

"""

import _weakref

from cpython.object cimport PyObject, PyTypeObject, destructor
from cysignals.memory cimport sig_malloc, sig_realloc, sig_free

from sage.structure.parent cimport Parent
from sage.structure.element cimport Element


cdef long DEFAULT_POOL_LENGTH = 100
cdef dict pools_disabled = {}
cdef dict pools_enabled = {}


cdef class Pool:
    """
    A class for handling pools in Sage.
    """
    def __cinit__(self, t):
        """
        Initialize this pool.

        By default the pool is empty and disabled.

        This method should not be called directly.
        Instead use the functions :func:`pool_enabled` and
        :func:`pool_disabled` which serve as factory. They
        in particular prevent the creation of many disabled
        pools for the same type.
        """
        self.type = <PyTypeObject*>t
        self.save_tp_dealloc = self.type.tp_dealloc
        self.type.tp_dealloc = &tp_dealloc
        self.enabled = 0

    def __repr__(self):
        """
        Return a string representation of this pool

            sage: R = Zp(3)
            sage: R.pool_enable()
            sage: R.pool()   # indirect doctest
            Pool of ... elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement
        """
        if self.enabled:
            s = "Pool of %s elements of type " % self.allocated
        else:
            s = "Disabled pool of type "
        s += str(<type>self.type)[7:-2]
        return s

    cdef Pool new_enabled(self):
        """
        Return a new enabled pool with the same type

        This method should not be called directly but through
        the function :func:`pool_enabled`.
        """
        cdef Pool pool = Pool(<type>self.type)
        pool.save_tp_dealloc = self.save_tp_dealloc
        pool.enabled = 1
        return pool

    def automatic_resize(self, enable=True):
        """
        Enable or disable automatic resizing for this pool

        When enabled, the length of the pool doubles when the
        pool is full and a new element needs to be added to
        the pool. In this case, deallocations never occur.

        INPUT:

        - `enable` -- boolean (default: True)
          whether we should enable or disable automatic
          resizing for this pool.

        NOTE::

        A disabled pool cannot be resized

        EXAMPLES::

            sage: R = Zp(3)
            sage: R.pool_enable()
            sage: pool = R.pool()
            sage: pool.length()
            100

            sage: pool.automatic_resize()
            sage: M = identity_matrix(R, 50)
            sage: M.determinant()
            1 + O(3^20)
            sage: pool.length()
            3200
            sage: pool.usage()
            2452

            sage: pool.automatic_resize(False)
            sage: pool.clear(); pool.resize(100)
            sage: M = identity_matrix(R, 50)
            sage: M.determinant()
            1 + O(3^20)
            sage: pool.length()
            100
            sage: pool.usage()
            99


        TESTS::

            sage: R = Zp(5)
            sage: R.pool_disable()
            sage: R.pool()
            Disabled pool of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

            sage: R.pool().automatic_resize()
            Traceback (most recent call last):
            ...
            ValueError: this pool is disabled
        """
        if enable:
            if self.enabled:
                self.type.tp_dealloc = &tp_dealloc_with_resize
            else:
                raise ValueError("this pool is disabled")
        else:
            self.type.tp_dealloc = &tp_dealloc

    def length(self):
        """
        Return the length of this pool (that is the total
        number of elements that can be put in this pool)

        EXAMPLES::

            sage: R = Zp(3)
            sage: R.pool_enable()
            sage: R.pool().length()
            100
        """
        from sage.rings.integer_ring import ZZ
        return ZZ(self.size)

    def usage(self):
        """
        Return the number of elements stored in this pool

        EXAMPLES::

            sage: R = Zp(3)
            sage: R.pool_enable()

            sage: R.pool().clear()
            sage: R.pool().usage()
            0

            sage: x = R.random_element()
            sage: del x   # here, x is added to the pool
            sage: R.pool().usage()
            1
        """
        from sage.rings.integer_ring import ZZ
        return ZZ(self.allocated)

    def resize(self, length=None):
        """
        Resize this pool

        INPUT:

        - `length` -- an integer or `None` (default: `None`)
          The new length of the pool. If none, the length is doubled.

        EXAMPLES::

            sage: R = Zp(3)
            sage: R.pool_enable()

            sage: pool = R.pool()
            sage: pool.length()
            100

            sage: pool.resize(500)
            sage: pool.length()
            500

        Without argument, the method :meth:`resize` doubles the length::

            sage: pool.resize()
            sage: pool.length()
            1000

        Of course, it is possible to lower the length as well.
        In this case, if there are elements above the new length, they
        are deallocated::

            sage: pool.resize(100)
            sage: pool.length()
            100

        TESTS::

            sage: R = Zp(5)
            sage: R.pool_disable()
            sage: R.pool()
            Disabled pool of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

            sage: R.pool().resize(1000)
            Traceback (most recent call last):
            ...
            ValueError: this pool is disabled
        """
        cdef long size
        if not self.enabled:
            raise ValueError("this pool is disabled")
        if length is None:
            size = 2 * self.size 
        else:
            size = <long?>length
        self.clear(size)
        self.elements = <PyObject**> sig_realloc(self.elements, size*sizeof(PyObject*))
        # should we test self.elements == NULL?
        self.size = size

    def clear(self, start=0):
        """
        Deallocated elements of the pool

        INPUT:

        - `start` -- an integer (default: 0)
          the position from which elements are deallocated

        EXAMPLES::

            sage: R = Zp(2)
            sage: R.pool_enable()

            sage: pool = R.pool();
            sage: pool.clear()
            sage: pool
            Pool of 0 elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

            sage: L = [ R.random_element() for _ in range(30) ]
            sage: del L
            sage: pool
            Pool of 30 elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

            sage: pool.clear(10)
            sage: pool
            Pool of 10 elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

            sage: pool.clear()
            sage: pool
            Pool of 0 elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement
        """
        cdef long s = <long?>start
        cdef long i
        if s < self.allocated:
            for i from s <= i < self.allocated:
                self.save_tp_dealloc(self.elements[i])
            self.allocated = s

    def __dealloc__(self):
        """
        Deallocate this pool
        """
        self.clear()
        if self.elements != NULL:
            sig_free(self.elements)


cdef inline PY_NEW_FROM_POOL(Pool pool):
    """
    Return a new object using the Pool `pool`

    INPUT:

    - `pool` -- a pool

    NOTE:

    If `pool` is not empty, an element is pulled for it.
    Otherwise, a new element with the type of the pool is
    created and returned.
    """
    cdef PyObject* o
    cdef PyTypeObject* t
    if pool.allocated > 0:
        #print("reuse")
        pool.allocated -= 1
        o = pool.elements[pool.allocated]
        o.ob_refcnt = 0
        return <object>o
    else:
        t = pool.type
        return t.tp_new(<type>t, <object>NULL, <object>NULL)


cdef void tp_dealloc(PyObject* o):
    """
    Add an element to a pool

    INPUT:

    - `o` -- an object

    NOTES:

    The pool is infered from the element using `o.parent().pool()` 
    (or more precisely, a Cython fast equivalent of this).

    If the pool is full, the element is deallocated.

    This function must not be called manually (even in Cython code).
    It is called automatically by Python when the object `o` is collected.
    """
    cdef Pool pool = (<Parent>(<Element>o)._parent)._pool
    if pool.allocated < pool.size:
        #print("add to pool")
        o.ob_refcnt = 1
        pool.elements[pool.allocated] = o
        pool.allocated += 1
    else:
        #print("dealloc")
        pool.save_tp_dealloc(o)


cdef void tp_dealloc_with_resize(PyObject* o):
    """
    Add an element to a pool

    INPUT:

    - `o` -- an object

    NOTES:

    The pool is infered from the element using `o.parent().pool()` 
    (or more precisely, a Cython fast equivalent of this).

    If the pool is full, it is resized.

    This function must not be called manually. 
    It is called automatically by Python when the object `o` is collected.
    """
    cdef Pool pool = (<Parent>(<Element>o)._parent)._pool
    if pool.allocated >= pool.size:
        pool.resize()
    #print("add to pool")
    o.ob_refcnt = 1
    pool.elements[pool.allocated] = o
    pool.allocated += 1


cdef pool_disabled(type t):
    """
    Return a disabled pool of requested type

    The result is cached. Therefore, two different calls
    to this function with the same argument return the same
    pool

    INPUT:

    - `t` -- a Python type

    EXAMPLES::

        sage: R = Zp(3)
        sage: R.pool_disable()
        sage: R.pool()   # indirect doctest
        Disabled pool of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

        sage: S = Zp(5)
        sage: S.pool_disable()
        sage: S.pool()   # indirect doctest
        Disabled pool of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

        sage: R.pool() is S.pool()
        True
    """
    cdef Pool pool = None
    if pools_disabled.has_key(t):
        wr_pool = pools_disabled[t]
        pool = wr_pool()
    if pool is None:
        pool = Pool(t)
    pools_disabled[t] = _weakref.ref(pool)
    return pool


cdef pool_enabled(Pool pool_dis, length, bint is_global):
    """
    Return an enabled pool

    INPUT:

    - `pool_dis` -- a pool
      the disabled pool from which the returned pool is
      constructed (in particular, the associated Python
      type is infered from here)

    - `length` -- an integer or `None`
      the length of the pool

    - `is_global` -- a boolean
      whether this pool should be marked as global (a global
      enabled pool is unique for a given type)

    NOTE:

    If `length` is specified and `is_global` is true, the
    global corresponding pool, if it already exists, is resized.

    If `length` is not specified and a new pool needs to be
    created, it is sized to the default length stored in 
    the global variable `DEFAULT_POOL_LENGTH`.

    EXAMPLES::

        sage: R = Zp(3)
        sage: R.pool_enable()

        sage: R.pool()   # indirect doctest
        Pool of ... elements of type sage.rings.padics.padic_capped_relative_element.pAdicCappedRelativeElement

    In this particular example, the pool is global. Indeed::

        sage: S = Zp(5)
        sage: S.pool_enable()
        sage: R.pool() is S.pool()
        True
    """
    cdef type t = <type>pool_dis.type
    cdef Pool pool = None
    if is_global:
        if pools_enabled.has_key(t):
            wr_pool = pools_enabled[t]
            pool = wr_pool()
        if pool is None:
            pool = pool_dis.new_enabled()
            if length is None:
                length = DEFAULT_POOL_LENGTH
        pools_enabled[t] = _weakref.ref(pool)
    else:
        pool = pool_dis.new_enabled()
        if length is None:
            length = DEFAULT_POOL_LENGTH
    if length is not None:
        pool.resize(length)
    return pool
