
#import os
#os.chdir(os.path.dirname(__file__))

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

    cdef cppclass Transition:
        int l # label
        int e # arrival state

    cdef cppclass NState:
        Transition* a
        int n
        bool final
        bool initial

    cdef cppclass NAutomaton:
        NState* e # states
        int n   # number of states
        int na  # number of letters



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
    cdef list S	 # states

#cdef Automaton getAutomaton (a, initial=?, F=?, A=?)
#cdef Dict list_to_Dict(list l)
