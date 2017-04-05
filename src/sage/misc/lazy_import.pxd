cdef int lazy_import_impl(module, names, as_, at_startup, namespace, deprecation) except -1

cdef inline int lazy_import(module, names, as_=None,
    at_startup=None, namespace=None, deprecation=None) except -1:
    """
    Lazy import function for Cython. This must be cimported.

    INPUT: see :func:`lazy_import`

    EXAMPLES::

        sage: cython('''
        ....: from sage.misc.lazy_import cimport lazy_import
        ....: global ZZ  # Avoid "undeclared name not builtin" errors
        ....: lazy_import('sage.rings.all', 'ZZ')
        ....:
        ....: print(type(ZZ))
        ....: x = ZZ(0)
        ....: print(type(ZZ))
        ....: ''')
        <type 'sage.misc.lazy_import.LazyImport'>
        <type 'sage.rings.integer_ring.IntegerRing_class'>
    """
    if namespace is None:
        namespace = globals()
    lazy_import_impl(module, names, as_, at_startup, namespace, deprecation)
