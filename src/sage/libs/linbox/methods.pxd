cdef extern from "linbox/solutions/methods.h":
    cdef struct LinBoxHybridSpecifier "LinBox::HybridSpecifier":
        pass
    cdef struct LinBoxBlackboxSpecifier "LinBox::BlackboxSpecifier":
        pass
    cdef struct LinBoxEliminationSpecifier "LinBox::LinBoxEliminationSpecifier":
        pass
    cdef struct LinBoxWiedemannTraits "LinBox::LinBoxWiedemannTraits":
        pass
    cdef struct LinBoxBlasEliminationTraits "LinBox::LinBoxBlasEliminationTraits":
        pass
    cdef struct LinBoxSparseEliminationTraits "LinBox::SparseEliminationTraits":
        pass

    cdef cppclass LinBoxMethod "LinBox::Method":
        ctypedef LinBoxHybridSpecifier Hybrid
        ctypedef LinBoxBlackboxSpecifier Blackbox
        ctypedef LinBoxEliminationSpecifier Elimination
        ctypedef LinBoxWiedemannTraits Wiedemann
        ctypedef LinBoxBlasEliminationTraits BlasElimination
        ctypedef LinBoxSparseEliminationTraits SparseElimination

