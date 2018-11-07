
cdef extern from "Automaton.h":
    ctypedef char bool
    cdef cppclass State:
        int* f
        bool final

    cdef cppclass Automaton:
        State* e # states
        int n   # number of states
        int na  # number of letters
        int i # initial state

from sage.combinat.words.cautomata cimport DetAutomaton

cdef class BetaAdicSet:
    cdef b
    cdef K
    cdef DetAutomaton a
    
    #cdef draw_zoom(self, n, int sx, int sy, bool ajust, int prec, color, int method, float coeff, bool verb)
    