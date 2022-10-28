from sage.structure.list_clone cimport ClonableArray

cdef class ImmutableListWithParent(ClonableArray):
    cdef _set_index(self, k, value)

cdef class TensorProductOfCrystalsElement(ImmutableListWithParent):
    pass

cdef class TensorProductOfRegularCrystalsElement(TensorProductOfCrystalsElement):
    cdef position_of_last_unmatched_minus(self, i)
    cdef position_of_first_unmatched_plus(self, i)

cdef class CrystalOfTableauxElement(TensorProductOfRegularCrystalsElement):
    pass

cdef class InfinityCrystalOfTableauxElement(CrystalOfTableauxElement):
    pass

cdef class InfinityCrystalOfTableauxElementTypeD(InfinityCrystalOfTableauxElement):
    pass

cdef class TensorProductOfSuperCrystalsElement(TensorProductOfRegularCrystalsElement):
    pass

cdef class CrystalOfBKKTableauxElement(TensorProductOfSuperCrystalsElement):
    pass

cdef class TensorProductOfQueerSuperCrystalsElement(TensorProductOfRegularCrystalsElement):
    pass

cdef class InfinityQueerCrystalOfTableauxElement(TensorProductOfQueerSuperCrystalsElement):
    cdef list _row_lengths

cdef Py_ssize_t count_leading(list row, letter)

