
cdef extern from "Automaton.h":
    ctypedef char bool
    cdef cppclass Etat:
        int* f
        bool final

    cdef cppclass Automaton:
        Etat* e # états
        int n   # nombre d'états
        int na  # nombre de lettres
        int i # état initial

    cdef cppclass Arete:
        int l # label
        int e # état d'arrivée

    cdef cppclass NEtat:
        Arete* a
        int n
        bool final
        bool initial

    cdef cppclass NAutomaton:
        NEtat* e # états
        int n   # nombre d'états
        int na  # nombre de lettres

cdef extern from "automataC.h":
    cdef cppclass Dict:
        int *e
        int n
    Automaton NewAutomaton(int n, int na)
    void FreeAutomaton(Automaton *a)
    NAutomaton NewNAutomaton(int n, int na)
    void FreeNAutomaton(NAutomaton *a)

cdef class FastAutomaton:
    cdef Automaton* a
    cdef list A	 # alphabet
    cdef dict dA  # dictionnary giving the index in A
    cdef list S	 # states
    cdef dict dS  # dictionnary giving the index in S
    # cdef set_a(self, Automaton a)


cdef class NFastAutomaton:
    cdef NAutomaton* a
    cdef list A
