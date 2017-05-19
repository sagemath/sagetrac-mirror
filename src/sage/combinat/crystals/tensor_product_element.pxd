from sage.structure.list_clone cimport ClonableArray

cdef class ImmutableListWithParent(ClonableArray):
    cpdef _set_index(self, k, value)

cdef class TensorProductOfCrystalsElement(ImmutableListWithParent):
    pass

cdef class TensorProductOfRegularCrystalsElement(TensorProductOfCrystalsElement):
    cpdef position_of_last_unmatched_minus(self, i)
    cpdef position_of_first_unmatched_plus(self, i)

cdef class CrystalOfTableauxElement(TensorProductOfRegularCrystalsElement):
    pass

cdef class InfinityCrystalOfTableauxElement(CrystalOfTableauxElement):
    pass

cdef class InfinityCrystalOfTableauxElementTypeD(InfinityCrystalOfTableauxElement):
    pass
