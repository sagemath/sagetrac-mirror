"""
Manage strong/weak references

TESTS:

While testing this file, we want deterministic garbage collection,
so we disable automatic GC::

    sage: import gc
    sage: gc.disable()
"""

cimport cython
from cpython.object cimport PyObject, PyTypeObject, traverseproc, inquiry, visitproc
from cpython.ref cimport Py_INCREF, Py_XDECREF


# Global "timer" (really a counter) for the number of traverse loops.
# Every time that tp_traverse is called with a *different* visitproc,
# this is increased.
cdef T_t T_traverse
cdef visitproc last_visit


@cython.no_gc
cdef class MultiRefManager:
    """
    Internal class for managing multiple :class:`MultiRef` instances
    to the same object.

    This should never be created directly by the user.
    """

    cdef MultiRef new_ref(self, bint strong):
        """
        Create a new strong/weak reference to the object.

        INPUT:

        - ``strong`` -- True for adding a strong reference, False for
          adding a weak reference.

        This assumes without checking that the object is still alive!
        """
        obj = <object>self.object
        r = <MultiRef>MultiRef.__new__(MultiRef)
        r.manager = self
        if strong:
            Py_INCREF(obj)
            self.numstrong += 1
        else:
            self.numweak += 1
        return r


cdef class MultiRef:
    """
    Manage multiple references to the same object, some of which are
    strong and some of which are weak. However, only the number of
    strong and weak references is specified: any given reference is
    dynamically considered strong or weak.

    INPUT:

    - ``object`` -- object to reference

    When constructing the initial ``MultiRef`` reference to an object
    using ``MultiRef(object)``, the single reference is considered
    strong. Additional strong/weak references to the same object can be
    created using the methods :meth:`strong_ref` and :meth:`weak_ref`
    on the initial ``MultiRef``. However, the newly created references
    are indistinguishable from the initial reference: they all play the
    same role.

    Whenever a ``MultiRef`` instance is deleted, it will delete a
    strong reference. If all strong references are deleted, the
    remaining ``MultiRef`` instances are all considered dead weak
    references. This is different from ``weakref.ref`` where a reference
    only dies when *all* existing references to the object are deleted.
    But ``MultiRef`` only cares about the references that it is managing
    itself (which is much simpler to implement).

    EXAMPLES::

        sage: from sage.cpython.multiref import MultiRef
        sage: A = MultiRef(object()); A
        <MultiRef: 1 strong, 0 weak, 0 external; to ...>
        sage: B = A.weak_ref(); B
        <MultiRef: 1 strong, 1 weak, 0 external; to ...>
        sage: C = A.strong_ref(); C
        <MultiRef: 2 strong, 1 weak, 0 external; to ...>
        sage: del B; A
        <MultiRef: 1 strong, 1 weak, 0 external; to ...>
        sage: del C; A
        <MultiRef: dead>

    Calling it returns the referenced object. An exception is raised
    for a dead reference (the latter is different from ``weakref.ref``
    which returns ``None`` in that case)::

        sage: from sage.cpython.multiref import MultiRef
        sage: A = MultiRef(ZZ)
        sage: B = A.weak_ref(); B
        <MultiRef: 1 strong, 1 weak, ... external; to Integer Ring>
        sage: B()
        Integer Ring
        sage: del A; B()
        Traceback (most recent call last):
        ...
        RuntimeError: dead MultiRef object

    It is possible to initialize multiple ``MultiRef`` instances to the
    same object. However, if these are not created using the
    :meth:`strong_ref` or :meth:`weak_ref` methods, those are
    unrelated::

        sage: from sage.cpython.multiref import MultiRef
        sage: A = MultiRef(object()); A2 = A.strong_ref()
        sage: B = MultiRef(A()); B2 = B.weak_ref()
        sage: A() is B()
        True
        sage: print(A); print(B)
        <MultiRef: 2 strong, 0 weak, 1 external; to ...>
        <MultiRef: 1 strong, 1 weak, 2 external; to ...>
        sage: del B2
        sage: print(A); print(B)
        <MultiRef: 2 strong, 0 weak, 0 external; to ...>
        <MultiRef: dead>

    **GARBAGE COLLECTION**

    As long as we are not collecting garbage, it does not matter which
    of the references is strong and which is weak: only the number of
    strong and weak references matter, since that determines the
    reference count of the object.

    During garbage collection, the ``tp_traverse`` algorithm of
    :class:`MultiRef` decides dynamically which references to
    consider strong and which ones to consider weak. This is done in
    such a way to maximize the amount of garbage that can be collected::

        sage: from sage.cpython.multiref import MultiRef
        sage: from weakref import ref
        sage: class X(object): pass
        sage: C = X()
        sage: A = MultiRef(C)
        sage: B = A.weak_ref()
        sage: C.refs = [A, B]
        sage: W = ref(C); del C

    At this point, we have global strong references to ``A`` and ``B``.
    There are strong references from ``C`` to ``A`` and ``B`` and two
    references (one strong, one weak, managed by ``MultiRef``)
    from ``A`` and ``B`` to ``C``. We also have an ordinary weak
    reference ``W`` to ``C``. When we collect garbage, everything
    remains alive::

        sage: import gc
        sage: _ = gc.collect()
        sage: W
        <weakref ...; to ...>

    Now we delete ``A`` and collect garbage. The reference from ``C``
    to ``A`` will be considered strong and the reference from ``C`` to
    ``B`` weak. As a consequence, the reference cycle between ``A`` and
    ``C`` will be deleted::

        sage: del A
        sage: _ = gc.collect()
        sage: W
        <weakref ...; dead>

    The reference from ``B`` to ``C`` is now a broken weak reference::

        sage: B
        <MultiRef: dead>

    **ALGORTIHM**

    As long as the garbage collector is not running, this class just
    needs to track reference counts. So the only non-trivial part is
    the ``tp_traverse`` implementation.

    We define a *traverse loop* as a loop of the form::

        for obj in set_of_objects:
            type(obj)->tp_traverse(obj, visit, arg)

    where ``visit`` is constant during the loop (``arg`` does not
    matter). More precisely, a traverse loop is a sequence of
    consecutive ``tp_traverse`` calls with the same ``visit``.
    Traverse loops can be detected by comparing the ``visit`` argument
    to the last ``visit``.
    
    During every traverse loop, the first references which are visited
    are considered weak and the last ones are considered strong. This
    targets the ``move_unreachable`` loop from CPython, which starts by
    visiting the known-reachable objects. We want to put weak references
    in the reachable ``MultiRef`` instances. If the number of
    known-reachable ``MultiRef`` instances is at most the number of
    allowed weak references, we will not see any strong reference and
    the object can be garbage collected (assuming that there are no
    other references).

    **ASSUMPTIONS**

    This algorithm makes 3 important assumptions on how the garbage
    collector from CPython works:

    1. Any garbage collection involves at least two traverse loops.
       Thefore, a traverse loop cannot span multiple garbage
       collections.

    2. The precise reference graph is allowed to change between traverse
       loops, as long as refcounts do not change and the reference graph
       is consistent within each individual traverse loop.

    3. Considering the first-visited references in a traverse loop as
       weak references maximizes the amount of garbage that can be
       collected.

    Assumption 1 is true because there are two traverse loops in
    CPython's garbage collector: one in ``subtract_refs`` and one in
    ``move_unreachable``.

    Assumptions 2 and 3 can be justified by noting that CPython's
    garbage collector does not keep track of references A -> B but only
    of the number of times that B is referenced. This can be seen in the
    signature of ``visit``: that takes only the referenced object (B)
    as argument, not the referencing object (A).
    """

    def __cinit__(self):
        # Initially not known to be weak
        self.T_weak = T_traverse - 1

    def __init__(self, obj):
        """
        Create a single strong reference to ``obj``.

        TESTS::

            sage: from sage.cpython.multiref import MultiRef
            sage: A = MultiRef.__new__(MultiRef)
            sage: A.__init__(object())
            sage: A
            <MultiRef: 1 strong, 0 weak, 0 external; to ...>
            sage: A = MultiRef(object())
            sage: A.__init__(object())
            Traceback (most recent call last):
            ...
            RuntimeError: already initialized
        """
        # Only used when creating a new MultiRef instance without a
        # prior manager. A single strong reference to "obj" is added.
        if self.manager is not None:
            raise RuntimeError("already initialized")
        manager = <MultiRefManager>MultiRefManager.__new__(MultiRefManager)
        Py_INCREF(obj)
        manager.object = <PyObject*>obj
        manager.numstrong = 1
        self.manager = manager

    def __dealloc__(self):
        clear(self)

    def strong_ref(self):
        """
        Create a new strong reference to the object and return a new
        :class:`MultiRef` instance pointing to the object.

        This does not mean that this specific reference is strong:
        all ``MultiRef`` references created from a given one are
        considered equivalent.

        EXAMPLES::

            sage: from sage.cpython.multiref import MultiRef
            sage: A = MultiRef(12345)
            sage: B = A.strong_ref(); B
            <MultiRef: 2 strong, 0 weak, 0 external; to 12345>
        """
        return self.get_manager().new_ref(True)

    def weak_ref(self):
        """
        Create a new weak reference to the object and return a new
        :class:`MultiRef` instance pointing to the object.

        This does not mean that this specific reference is weak:
        all ``MultiRef`` references created from a given one are
        considered equivalent.

        EXAMPLES::

            sage: from sage.cpython.multiref import MultiRef
            sage: A = MultiRef(12345)
            sage: B = A.weak_ref(); B
            <MultiRef: 1 strong, 1 weak, 0 external; to 12345>
        """
        return self.get_manager().new_ref(False)

    def __call__(self):
        """
        Return the managed object. Raise an exception if the reference
        is dead.

        EXAMPLES::

            sage: from sage.cpython.multiref import MultiRef
            sage: A = MultiRef(ZZ)
            sage: B = A.weak_ref()
            sage: A()
            Integer Ring
            sage: del B; A()
            Traceback (most recent call last):
            ...
            RuntimeError: dead MultiRef object
        """
        return <object>self.get_manager().object

    def is_alive(self):
        """
        Return whether the reference to the object is still alive.

        EXAMPLES::

            sage: from sage.cpython.multiref import MultiRef
            sage: A = MultiRef(object())
            sage: B = A.weak_ref()
            sage: A.is_alive()
            True
            sage: del B
            sage: A.is_alive()
            False
        """
        return self.alive()

    def __repr__(self):
        """
        TESTS::

            sage: from sage.cpython.multiref import MultiRef
            sage: MultiRef.__new__(MultiRef)
            <MultiRef: dead>
            sage: MultiRef(ZZ)
            <MultiRef: 1 strong, 0 weak, ... external; to Integer Ring>
            sage: _.weak_ref()
            <MultiRef: dead>
        """
        if not self.alive():
            desc = "dead"
        else:
            M = self.manager
            obj = M.object
            numext = obj.ob_refcnt - M.numstrong
            r = repr(<object>obj)
            desc = f"{M.numstrong} strong, {M.numweak} weak, {numext} external; to {r}"
        return f"<{type(self).__name__}: {desc}>"


cdef int traverse(MultiRef self, visitproc visit, void* arg):
    if not self.alive():
        return 0

    # Update T_traverse and self.T_cur
    global T_traverse, last_visit
    if visit is not last_visit:
        # New tp_traverse() loop
        last_visit = visit
        T_traverse += 1

    M = self.manager
    if M.T_cur != T_traverse:
        # M was not yet visited in this tp_traverse() loop
        # Reset all known weak references
        M.T_cur = T_traverse
        M.curweak = 0
    elif self.T_weak == T_traverse:
        # Already a weak reference
        return 0

    # If we can make it a weak reference, do so
    if M.curweak < M.numweak:
        self.T_weak = T_traverse
        M.curweak += 1
        return 0

    # Must be strong reference
    self.T_weak = T_traverse - 1
    return visit(M.object, arg)


cdef int clear(MultiRef self) except -1:
    # Common implementation for tp_clear and tp_dealloc.

    # We check self.manager and set it to None as protection against
    # double frees. We want to run clear() exactly once.
    M = self.manager
    if M is None:
        return 0
    self.manager = None

    # Invalidate last_visit to ensure a consistent state
    global last_visit
    last_visit = NULL

    # Clear strong references first
    obj = M.object
    if obj is not NULL:
        # It is important to ensure that M is in a consistent state
        # before the DECREF, since the latter might cause nasty
        # recursions.
        M.numstrong -= 1
        if M.numstrong == 0:
            M.object = NULL
        Py_XDECREF(obj)
    else:
        M.numweak -= 1


(<PyTypeObject*>MultiRef).tp_traverse = <traverseproc>traverse
(<PyTypeObject*>MultiRef).tp_clear = <inquiry>clear
