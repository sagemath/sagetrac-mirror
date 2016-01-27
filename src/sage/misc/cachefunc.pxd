from function_mangling cimport ArgumentFixer

cpdef dict_key(o)
cpdef cache_key(o)

cdef class CachedFunction(object):
    cdef public str __name__
    cdef public str __module__
    cdef ArgumentFixer _argument_fixer
    cdef public f
    cdef public cache  # not always of type <dict>
    cdef bint is_classmethod
    cdef int argfix_init(self) except -1
    cdef get_key_args_kwds(self, tuple args, dict kwds)
    cdef fix_args_kwds(self, tuple args, dict kwds)
    cdef empty_key
    cdef key

cdef class CachedMethod(object):
    cdef str _cache_name
    cdef public str __name__
    cdef public str __module__
    cdef CachedFunction _cachedfunc
    cdef Py_ssize_t nargs
    cpdef dict _get_instance_cache(self, inst)

cdef class CachedInParentMethod(CachedMethod):
    pass

cdef class CachedMethodCaller(CachedFunction):
    cdef public _instance
    cdef public CachedMethod _cachedmethod

cdef class CachedMethodCallerNoArgs(CachedFunction):
    cdef public _instance

cdef class GloballyCachedMethodCaller(CachedMethodCaller):
    pass
