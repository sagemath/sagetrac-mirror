"""
Object pools

An object pool is an array which keeps deleted objects of a certain
type instead of deallocating them. When a new object is needed, an
object from the pool is returned instead of allocating a new object.
The goal is faster creation and deletion of objects.

In this implementation, a pool is attached to a specific type, which
must be an extension type (for technical reasons, but also because it
makes little sense to optimize Python classes).

Using a pool in the most naive way would simply optimize the allocation
and deallocation of the Python object itself. However, we can try to do
better and also optimize other initialization (and deinitialization)
steps of the object. For example, an :class:`Integer` is essentially a
wrapper around an ``mpz_t``. When using a pool for integers, we do not
need to deallocate and reallocate this ``mpz_t`` structure.

It should already be clear from this example that it is hard to have one
pool implementation which works for every type. Therefore, this is
implemented as a base class :class:`ObjectPool` with several hooks to
customize the behaviour of the pool. We also provide two derived
classes: :class:`ObjectPool_NoGC` which disables cyclic garbage
collection on the type. Deriving from this, :class:`ObjectPool_Cloning`
creates entirely new objects from a "template" object. This way, we can
also optimize the case of creating new objects when the pool is empty.

One restriction to the use of pools is that any arguments to ``__new__``
are ignored. In other words, the object taken from the pool cannot
depend on the arguments to ``__new__``. This is not a severe restriction
since the arguments can still be handled in ``__init__``.

It is good to realize that there are two ways how an object can be
deleted: the first is that the refcount explicitly goes to zero. At that
point, the object is deleted immediately. A second possibility is that
the cyclic garbage collector collects the object. In the latter case,
parts of the object may already be deallocated. When testing, be sure to
test both of these scenarios. Alternatively, use :class:`ObjectPool_NoGC`
to disable garbage collection for the type.

EXAMPLES::

    sage: cython('''
    ....: from sage.ext.pool cimport ObjectPool
    ....: cdef class X:
    ....:     def __cinit__(self):
    ....:         print("calling __cinit__")
    ....:     def __dealloc__(self):
    ....:         print("calling __dealloc__")
    ....: def start_pool():
    ....:     # Pool of size 1
    ....:     ObjectPool(X, 1).enable()
    ....: ''')
    sage: x0 = X()
    calling __cinit__
    sage: start_pool()
    sage: x1 = X()
    calling __cinit__
    sage: prev_id = id(x1)
    sage: del x1
    sage: x1 = X()
    sage: id(x1) == prev_id
    True
    sage: del x0
    sage: del x1
    calling __dealloc__
"""

#*****************************************************************************
#       Copyright (C) 2017 Jeroen Demeyer <J.Demeyer@UGent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cdef extern from *:
    void Py_INCREF(PyObject*)
    void _Py_NewReference(PyObject*)
    void _Py_ForgetReference(PyObject*)
    PyTypeObject* Py_TYPE(PyObject*)
    bint PyType_HasFeature(PyTypeObject*, long)


from libc.string cimport memcpy
from cysignals.memory cimport check_calloc, check_malloc, sig_free

cdef tuple empty_tuple = ()


cdef class ObjectPool:
    """
    A pool of objects of a certain type for fast
    allocation/deallocation.

    INPUT:

    - ``cls`` -- the type of the objects. This must be an extension type
      (not a Python ``class``).

    - ``size`` -- (integer, default: 128) The maximal number of objects
      to store in the pool.

    Simply creating the ``ObjectPool`` object does not change anything
    yet. The pool is only used after calling the ``enable()`` method.

    The methods ``pool_enable``, ``pool_in``, ``pool_out``, ``pool_new``,
    ``pool_dealloc`` should be seen as hooks which can (and in many
    cases must) be customized for a specific type.

    EXAMPLES:

    We demonstrate all hook functions::

        sage: cython('''
        ....: from sage.ext.pool cimport ObjectPool
        ....: cdef class DebugPool(ObjectPool):
        ....:     cdef int pool_enable(self) except -1:
        ....:         print("pool_enable")
        ....:     cdef void pool_in(self, obj):
        ....:         print("pool_in")
        ....:     cdef int pool_out(self, obj) except -1:
        ....:         print("pool_out")
        ....:     cdef pool_new(self):
        ....:         print("pool_new")
        ....:         return ObjectPool.pool_new(self)
        ....:     cdef void pool_dealloc(self, obj):
        ....:         print("pool_dealloc")
        ....:         ObjectPool.pool_dealloc(self, obj)
        ....: cdef class X:
        ....:     def __cinit__(self):
        ....:         print("__cinit__")
        ....:     def __dealloc__(self):
        ....:         print("__dealloc__")
        ....: DebugPool(X, 1).enable()
        ....: ''')
        pool_enable
        sage: x0 = X()
        pool_new
        __cinit__
        sage: del x0
        pool_in
        sage: x0 = X()
        pool_out
        sage: x1 = X()
        pool_new
        __cinit__
        sage: del x0
        pool_in
        sage: del x1
        pool_dealloc
        __dealloc__

    Subclasses bypass the pool::

        sage: class Y(X): pass
        sage: y = Y()
        __cinit__
        sage: del y
        __dealloc__
    """
    def __cinit__(self, cls, Py_ssize_t size=128):
        """
        EXAMPLES::

            sage: from sage.ext.pool import ObjectPool
            sage: ObjectPool(Rational, size=1234)
            ObjectPool for Rational objects (disabled; 0/1234 filled)

        TESTS::

            sage: ObjectPool(Rational, -1)
            Traceback (most recent call last):
            ...
            ValueError: pool size must be a positive number
            sage: ObjectPool("hello")
            Traceback (most recent call last):
            ...
            TypeError: expected type, got 'hello'
            sage: ObjectPool(Graphics)
            Traceback (most recent call last):
            ...
            TypeError: <class 'sage.plot.graphics.Graphics'> is not an extension type, making it incompatible with ObjectPool
            sage: ObjectPool(tuple)
            Traceback (most recent call last):
            ...
            TypeError: <... 'tuple'> has variable size, making it incompatible with ObjectPool
            sage: ObjectPool(type(None))
            Traceback (most recent call last):
            ...
            TypeError: <... 'NoneType'> is missing tp_version_tag
        """
        if not isinstance(cls, type):
            raise TypeError(f"expected type, got {cls!r}")
        cdef PyTypeObject* t = <PyTypeObject*>cls

        if PyType_HasFeature(t, Py_TPFLAGS_HEAPTYPE):
            raise TypeError(f"{cls!r} is not an extension type, making it incompatible with ObjectPool")
        if t.tp_itemsize:
            raise TypeError(f"{cls!r} has variable size, making it incompatible with ObjectPool")
        if not PyType_HasFeature(t, Py_TPFLAGS_HAVE_VERSION_TAG):
            # Note: it's not clear whether we really require the version
            # tag, but we use this as sanity check that the type is
            # sufficiently implemented. Cython enables this flag on all
            # extension types, so we should be fine.
            raise TypeError(f"{cls!r} is missing tp_version_tag")

        if size <= 0:
            raise ValueError("pool size must be a positive number")

        self.cls = <type>cls
        self.real_tp_new = t.tp_new
        self.real_tp_dealloc = t.tp_dealloc

        self.size = size
        self.used = 0
        self.pool = <PyObject**>check_calloc(self.size, sizeof(PyObject*))

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.ext.pool import ObjectPool
            sage: ObjectPool(Parent)
            ObjectPool for Parent objects (disabled; 0/128 filled)
            sage: Integer.__ObjectPool
            IntegerPool for Integer objects (enabled; ... filled)
        """
        enabled = "enabled" if self.is_enabled() else "disabled"
        return f"{type(self).__name__} for {self.cls.__name__} objects ({enabled}; {self.used}/{self.size} filled)"

    cdef tp_new(self, cls, args, kwds):
        """
        Implementation of ``(self.cls).__new__(cls, *args, **kwds)``.

        This is called by the actual ``tp_new`` function of the type
        ``self.cls``, which is implemented in ``pool_impl.h``

        This may not be overridden in subclasses.
        """
        if cls is not self.cls:
            # We are just a subclass, bypass the pool
            return self.real_tp_new(<type>cls, args, kwds)
        if not self.used:
            return self.pool_new()
        self.used -= 1
        cdef object obj = <object>(self.pool[self.used])
        self.pool_out(obj)
        return obj

    cdef void tp_dealloc(self, PyObject* obj):
        """
        Deallocate ``obj``.

        This is called by the actual ``tp_dealloc`` function of the
        type ``self.cls``, which is implemented in ``pool_impl.h``

        This may not be overridden in subclasses.
        """
        if Py_TYPE(obj) is not self.tp():
            # We are just a subclass, bypass the pool
            self.real_tp_dealloc(obj)
            return

        # The refcount of obj should be zero at this point. We
        # temporarily increase its refcount to prevent an accidental
        # DECREF to zero refcount again.
        Py_INCREF(obj)

        if self.used == self.size:
            self.pool_dealloc(<object>obj)
            return

        self.pool_in(<object>obj)

        self.pool[self.used] = obj
        self.used += 1

        Py_DECREF_KEEP(obj)

    cdef int pool_enable(self) except -1:
        """
        Called during ``enable()``.

        This is a hook for subclasses. Subclasses which override this
        would typically call the ``pool_enable`` method of the base
        class followed by their own customizations.
        """
        pass

    cdef void pool_in(self, obj):
        """
        Take a live object and process it to store it in the pool.

        This is a hook for subclasses. Typically, the implementation
        would be similar to ``__dealloc__`` in Cython: deallocating
        memory which is no longer needed. Also, attributes should be
        set to None such that those referenced objects can be deleted.

        This default implementation clears the object in the same way
        that the garbage collector would do. This makes sense, since the
        object might have been cleared anyway by the garbage collector.
        In Cython, this means that all Python attributes are set to
        ``None``. Other attributes are unchanged.

        EXAMPLES::

            sage: cython('''
            ....: from sage.ext.pool cimport ObjectPool
            ....: cdef class X:
            ....:     cdef public attr
            ....:     cdef public int num
            ....: ObjectPool(X).enable()
            ....: ''')
            sage: x = X()
            sage: x.attr = "hello"
            sage: x.num = 42
            sage: del x
            sage: x = X()
            sage: print(x.attr)
            None
            sage: print(x.num)
            42
        """
        tp_clear = self.tp().tp_clear
        if tp_clear:
            tp_clear(obj)

    cdef int pool_out(self, obj) except -1:
        """
        Take an object from the pool and process it to make it live.

        This is a hook for subclasses. Typically, this is similar to
        ``__cinit__`` in Cython. One difference is that ``__cinit__``
        starts from an "empty" object (all C attributes are 0 or NULL,
        all Python attributes are None). On the other hand, the input to
        ``pool_out`` is the result of an earlier ``pool_in`` or
        ``pool_clone`` call.
        """
        pass

    cdef pool_new(self):
        """
        Create a completely new object, bypassing the pool and return
        that object.

        This is a hook for subclasses.
        """
        return self.real_tp_new(self.cls, empty_tuple, <dict>NULL)

    cdef void pool_dealloc(self, obj):
        """
        Completely deallocate and free an object, bypassing the pool.

        This is a hook for subclasses.
        """
        self.real_tp_dealloc(<PyObject*>obj)


cdef class ObjectPool_NoGC(ObjectPool):
    """
    Variant of :class:`ObjectPool` which disables cyclic garbage
    collection on the type.

    If the type already doesn't enable GC, there is no need for this:
    just use a plain :class:`ObjectPool`.

    .. WARNING::

        This comes with many pitfalls and can easily cause segfaults
        if used badly.

    EXAMPLES:

    Check that instances are indeed not tracked by the garbage
    collector. Subclasses are::

        sage: cython('''
        ....: from sage.ext.pool cimport ObjectPool_NoGC
        ....: cdef class X:
        ....:     cdef attr
        ....: ObjectPool_NoGC(X).enable()
        ....: cdef class Y(X):
        ....:     pass
        ....: ''')
        sage: x = X()
        sage: y = Y()
        sage: from gc import is_tracked
        sage: is_tracked(x)
        False
        sage: is_tracked(y)
        True
    """
    cdef int pool_enable(self) except -1:
        """
        Disable garbage collection.
        """
        self.tp().tp_flags &= ~<unsigned long>Py_TPFLAGS_HAVE_GC

    cdef void pool_in(self, obj):
        # Override ObjectPool.pool_in() which would call tp_clear()
        pass

    cdef void pool_dealloc(self, obj):
        # We cannot call real_tp_dealloc as that will almost certainly
        # do the wrong thing if that was implemented with GC in mind.
        #
        # Note that this exception cannot actually be raised since we
        # are in a dealloc function returning void. Still, we should get
        # a nice traceback.
        raise NotImplementedError("you must implement pool_dealloc")


cdef class ObjectPool_Cloning(ObjectPool_NoGC):
    """
    Variant of :class:`ObjectPool_NoGC` which creates new objects by
    cloning a template object.

    .. WARNING::

        This comes with huge pitfalls and can easily cause segfaults
        if used badly. In particular, cloning an object can cause
        reference counts to be wrong if the template contains Python
        objects which are blindly copied/freed without updating
        refcounts.

    EXAMPLES::

        sage: cython('''
        ....: from sage.ext.pool cimport ObjectPool_Cloning
        ....: cdef class X:
        ....:     def __cinit__(self):
        ....:         print("__cinit__")
        ....:     def __dealloc__(self):
        ....:         print("__dealloc__")
        ....: ObjectPool_Cloning(X).enable()
        ....: ''')
        __cinit__
        sage: x0 = X()
        sage: x1 = X()
        sage: x2 = X()
        sage: del x0
        sage: del x1
        sage: del x2
    """
    cdef int pool_enable(self) except -1:
        """
        Create the template object to be cloned.
        """
        ObjectPool_NoGC.pool_enable(self)
        cdef object obj = ObjectPool_NoGC.pool_new(self)
        self.clone_in(obj)
        self.template = <PyObject*>obj
        # Safely delete obj without deallocating
        Py_INCREF(self.template)
        del obj
        Py_DECREF_KEEP(self.template)

    cdef pool_new(self):
        # Copy template object
        cdef size_t size = self.tp().tp_basicsize
        cdef PyObject* x = <PyObject*>check_malloc(size)
        memcpy(x, self.template, size)
        # Set recount to 1 and do some debugging stuff if Python was
        # compiled with debugging enabled. It is important that we call
        # _Py_NewReference(), otherwise we would get crashes in debug
        # builds.
        _Py_NewReference(x)
        # The <object> cast below already increases the refcount, so we
        # decrease the refcount right away.
        Py_DECREF_KEEP(x)
        cdef object obj = <object>x
        self.clone_out(obj)
        return obj

    cdef void pool_dealloc(self, obj):
        self.clone_in(obj)
        _Py_ForgetReference(<PyObject*>obj)  # Needed for debug builds
        sig_free(<void*>obj)

    cdef int clone_out(self, obj) except -1:
        """
        Take a copy from the template object and prepare it to make
        it live.

        Typically, this means allocating any needed memory.
        """
        pass

    cdef void clone_in(self, obj):
        """
        Called once to create the template. This method receives a new
        object (created with the real ``__new__`` of the type) and
        prepares it for storing as template.

        Typically, this means freeing any memory which was allocated.
        """
        pass


ObjectPool_new = ObjectPool.tp_new
ObjectPool_dealloc = ObjectPool.tp_dealloc
