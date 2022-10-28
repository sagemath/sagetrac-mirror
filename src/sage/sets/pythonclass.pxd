from sage.structure.parent cimport Set_generic


cdef class Set_PythonType_class(Set_generic):
    cdef type _type


cdef Set_PythonType(typ)
