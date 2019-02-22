# coding=utf8
"""
Deterministic automaton of Finite state machines using C
DeterministicAutomaton for determinist automata and CAutomaton for non necessarily determinist automaton


AUTHORS:

- Paul Mercat (2013)- I2M AMU Aix-Marseille Universite - initial version
- Dominique Benielli (2018) Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in -SageMath

REFERENCES:

.. [Hopcroft] "Around Hopcroft’s Algorithm"  Manuel of BACLET and
   Claire PAGETTI.

"""

# *****************************************************************************
#       Copyright (C) 2014 Paul Mercat <paul.mercat@univ-amu.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from __future__ import print_function
from libc.stdlib cimport malloc, free
from sage.graphs.digraph import DiGraph
from cysignals.signals cimport sig_on, sig_off, sig_check
from cpython cimport bool as c_bool
from sage.combinat.words.cautomata_generators import DeterministicAutomatonGenerators
cimport sage.combinat.words.cautomata

dag = DeterministicAutomatonGenerators()

cdef extern from "Automaton.h":

    cdef cppclass Transition:
        int l  # label (-1 : epsilon-transition)
        int e  # arrival state

cdef extern from "automataC.h":

    cdef cppclass Dict:
        int* e
        int n
    cdef cppclass InvertDict:
        Dict* d
        int n

    bool DotExists()
    #    Automaton NewAutomaton (int n, int na)
    #    void FreeAutomaton (Automaton *a)
    int hashAutomaton(Automaton a)
    Automaton CopyAutomaton(Automaton a, int nalloc, int naalloc)
    NAutomaton CopyNAutomaton(NAutomaton a, int nalloc, int naalloc)
    Automaton PieceAutomaton(Automaton a, int *w, int n, int e)
    void initAutomaton(Automaton *a)
    void initNAutomaton(NAutomaton *a)
    void printAutomaton(Automaton a)
    void plotDot(const char *file, Automaton a, const char **labels, const char *graph_name, double sx, double sy, const char **vlabels, bool html, bool verb, bool run_dot)
    void NplotDot(const char *file, NAutomaton a, const char **labels, const char *graph_name, double sx, double sy, bool run_dot)
    Automaton Product(Automaton a1, Automaton a2, Dict d, bool verb)
    Automaton Determinize(Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
    Automaton DeterminizeN(NAutomaton a, bool sink, int verb)
    NAutomaton Concat(Automaton a, Automaton b, bool verb)
    NAutomaton CopyN(Automaton a, bool verb)
    void AddTransitionN(NAutomaton *a, int e, int f, int l)
    void AddPathN(NAutomaton *a, int e, int f, int *l, int len, bool verb)
    NAutomaton Proj(Automaton a, Dict d, bool verb)
    void ZeroComplete(Automaton *a, int l0, bool verb)
    Automaton ZeroComplete2(Automaton *a, int l0, bool sink_state, bool verb)
    Automaton ZeroInv(Automaton *a, int l0)
    Automaton prune_inf(Automaton a, bool verb)
    Automaton prune(Automaton a, bool verb)
    Automaton pruneI(Automaton a, bool verb)
    void AccCoAcc(Automaton *a, int *coa)
    void CoAcc(Automaton *a, int *coa)
    bool equalsAutomaton(Automaton a1, Automaton a2, bool verb)
    Dict NewDict(int n)
    void FreeDict(Dict *d)
    void printDict(Dict d)
    InvertDict NewInvertDict(int n)
    void FreeInvertDict(InvertDict id)
    void printInvertDict(InvertDict id)
    Automaton Duplicate(Automaton a, InvertDict id, int na2, bool verb)
    Automaton MirrorDet(Automaton a)
    NAutomaton Mirror(Automaton a)
    int StronglyConnectedComponents(Automaton a, int *res)
    Automaton SubAutomaton(Automaton a, Dict d, bool verb)
    Automaton Permut(Automaton a, int *l, int na, bool verb)
    void PermutOP(Automaton a, int *l, int na, bool verb)
    Automaton Minimise(Automaton a, bool verb)
    void DeleteVertexOP(Automaton* a, int e)
    Automaton DeleteVertex(Automaton a, int e)
    bool equalsLanguages(Automaton *a1, Automaton *a2, Dict a1toa2, bool minimized, bool pruned, bool verb)
    bool Intersect(Automaton a1, Automaton a2, bool verb)
    bool Included(Automaton a1, Automaton a2, bool pruned, bool verb)
    # bool intersectLanguage (Automaton *a1, Automaton *a2, Dict a1toa2, bool pruned, bool verb)
    bool emptyLanguage(Automaton a)
    void AddState(Automaton *a, bool final)
    bool IsCompleteAutomaton(Automaton a)
    bool CompleteAutomaton(Automaton *a)
    Automaton BiggerAlphabet(Automaton a, Dict d, int nna) #copy the automaton with a new bigger alphabet
    bool findWord(Automaton a, Dict *w, bool verb)
    bool shortestWord(Automaton a, Dict *w, int i, int f, bool verb)
    bool shortestWords(Automaton a, Dict *w, int i, bool verb)
    bool recognized_word(Automaton a, Dict d)
    void Test()

cdef imagDict(dict d, list A, list A2=[]):
    """
    Dictionary which is numbering projected alphabet
    """
    d1 = {}
    i = 0
    for a in A:
        if d.has_key(a):
            if not d1.has_key(d[a]):
                d1[d[a]] = i
                A2.append(d[a])
                i += 1
    return d1

cdef imagDict2(dict d, list A, list A2=[]):
    """
    Dictionary which is numbering a new alphabet
    """
    # print("d=%s, A=%s"%(d,A))
    d1 = {}
    i = 0
    for a in A:
        if d.has_key(a):
            for v in d[a]:
                if not d1.has_key(v):
                    d1[v] = i
                    A2.append(v)
                    i += 1
    return d1

cdef Dict getDict(dict d, list A, dict d1=None):
    A = list(A)
    cdef Dict r
    sig_on()
    r = NewDict(len(A))
    sig_off()
    cdef int i
    if d1 is None:
        d1 = imagDict(d, A)
    # print d1
    for i in range(r.n):
        if d.has_key(A[i]):
            r.e[i] = d1[d[A[i]]]
        else:
            r.e[i] = -1
    return r

cdef Dict list_to_Dict(list l):
    cdef Dict d = NewDict(len(l))
    cdef int i
    for i in range(len(l)):
        d.e[i] = l[i]
    return d

cdef InvertDict getDict2(dict d, list A, dict d1=None):
    A = list(A)
    cdef InvertDict r
    sig_on()
    r = NewInvertDict(len(A))
    sig_off()
    cdef int i
    if d1 is None:
        d1 = imagDict2(d, A)
    # print(d1)
    for i in range(r.n):
        if d.has_key(A[i]):
            sig_on()
            r.d[i] = NewDict(len(d[A[i]]))
            sig_off()
            for j in range(r.d[i].n):
                r.d[i].e[j] = d1[d[A[i]][j]]
        else:
            r.d[i].n = 0
    return r

cdef imagProductDict(dict d, list A1, list A2, list Av=[]):
    """
    Dictionary numbering the projected alphabet
    """
    dv = {}
    i = 0
    for a1 in A1:
        for a2 in A2:
            if d.has_key((a1, a2)):
                if not dv.has_key(d[(a1, a2)]):
                    dv[d[(a1, a2)]] = i
                    Av.append(d[(a1, a2)])
                    i += 1
    return dv

cdef Dict getProductDict(dict d, list A1, list A2, dict dv=None, verb=True):
    cdef Dict r
    d1 = {}
    d2 = {}
    cdef int i, n1, n2
    n1 = len(A1)
    for i in range(n1):
        d1[A1[i]] = i
    if verb:
        print(d1)
    n2 = len(A2)
    for i in range(n2):
        d2[A2[i]] = i
    if verb:
        print(d2)
    if dv is None:
        dv = imagProductDict(d, A1, A2)
    sig_on()
    r = NewDict(n1*n2)
    sig_off()
    Keys = d.keys()
    if verb:
        print("Keys=%s" % Keys)
    for (a1, a2) in Keys:
        if d1.has_key(a1) and d2.has_key(a2):
            r.e[d1[a1]+d2[a2]*n1] = dv[d[(a1, a2)]]
    return r

cdef Automaton getAutomaton(a, initial=None, F=None, A=None):
    d = {}
    da = {}
    if F is None:
        if not hasattr(a, 'F'):
            try:
                F = a.vertices()
            except AttributeError:
                    raise AttributeError("No vertices() method")
        else:
            F = a.F
    cdef Automaton r

    if A is None:
        A = list(set(a.edge_labels()))
    V = list(a.vertices())

    cdef int n = len(V)
    cdef int na = len(A)

    sig_on()
    r = NewAutomaton(n, na)
    initAutomaton(&r)
    sig_off()
    for i in range(na):
        da[A[i]] = i
    for i in range(n):
        r.e[i].final = 0
        d[V[i]] = i
    for v in F:
        #if not d.has_key(v):
        #    sig_on()
        #    FreeAutomaton(&r)
        #    r = NewAutomaton(0, 0)
        #    sig_off()
        #    raise ValueError("Incorrect set of final states.")
        if d.has_key(v):
            r.e[d[v]].final = 1

    if initial is None:
        if not hasattr(a, 'I'):
            I = []
        else:
            I = list(a.I)
        if len(I) > 1:
            raise ValueError("The automata must be determistic (I=%s)" % I)
        if len(I) >= 1:
            r.i = d[I[0]]
        else:
            r.i = -1
    else:
        if d.has_key(initial):
            r.i = d[initial]
        else:
            r.i = -1

    w = False
    for e, f, l in a.edges(sort=False):
        if r.e[d[e]].f[da[l]] != -1:
            if not w:
                print("Warning: the automaton was not deterministic! (edge %s -%s-> %s and edge %s -%s-> %s)\nThe result lost some informations."% (d[e], l, r.e[d[e]].f[da[l]], d[e], l, d[f]))
            w = True
        r.e[d[e]].f[da[l]] = d[f]
    #a.dA = da
    a.S = V
    #a.dS = d
    return r

cdef NAutomaton getNAutomaton(a, I=None, F=None, A=None, verb=False):
    d = {}
    da = {}
    dt = {}
    if F is None:
        if not hasattr(a, 'F'):
            F = a.vertices()
        else:
            F = a.F
    if I is None:
        if not hasattr(a, 'I'):
            I = a.vertices()
        else:
            I = a.I
    cdef NAutomaton r

    if A is None:
        A = list(a.Alphabet)

    V = list(a.vertices())
    cdef int n = len(V)
    cdef int na = len(A)

    if verb:
        print("vertices=%s" % V)

    if verb:
        print("Alloc(%s, %s) and init..." % (n, na))
    sig_on()
    r = NewNAutomaton(n, na)
    # initNAutomaton(&r)
    sig_off()
    if verb:
        print("init...")
    for i in range(n):
        d[V[i]] = i
        dt[V[i]] = []
    for i in range(na):
        da[A[i]] = i
    if verb:
        print("final states...")
    for v in F:
        if not d.has_key(v):
            sig_on()
            FreeNAutomaton(&r)
            sig_off()
            raise ValueError("Incorrect set of final states.")
        if verb:
            print("d[%s]=%s" % (v, d[v]))
        r.e[d[v]].final = 1
    if verb:
        print("initial states...")
    for v in I:
        if not d.has_key(v):
            sig_on()
            FreeNAutomaton(&r)
            sig_off()
            raise ValueError("Incorrect set of final states.")
        if verb:
            print("d[%s]=%s" % (v, d[v]))
        r.e[d[v]].initial = 1

    if verb:
        print(dt)
        print("browse transitions...")

    for e, f, l in a.edges():
        if dt.has_key(e):
            dt[e].append((f, l))
        else:
            raise MemoryError("This transition %s -%s-> %s starts from a state not in the set of states."%(e, l, f))
    if verb:
        print(dt)
        print("Alloc...")
    for e in V:
        r.e[d[e]].n = len(dt[e])
        if verb:
            print("alloc r.e[%s].a (%s elements, %so total)" % (d[e], r.e[d[e]].n, sizeof(Transition)*r.e[d[e]].n))
        if r.e[d[e]].n > 0:
            r.e[d[e]].a = <Transition *>malloc(sizeof(Transition)*r.e[d[e]].n)
            if r.e[d[e]].a is NULL:
                raise MemoryError("Could not allocate the transition (size %s)"%(sizeof(Transition)*len(dt[e])))
            for i, (f, l) in enumerate(dt[e]):
                if verb:
                    print("r.e[%s].a[%s] = %s,%s" % (d[e], i, l, f))
                r.e[d[e]].a[i].l = da[l]
                r.e[d[e]].a[i].e = d[f]
        else:
            r.e[d[e]].a = NULL
            r.e[d[e]].n = 0
    a.S = V
    return r

#cdef AutomatonGet(Automaton a, A):
#    """
#    Transform an Automaton a with an alphabet A to a DiGraph
#    """
#    from sage.graphs.digraph import DiGraph
#    r = DiGraph(multiedges=True, loops=True)
#    cdef int i, j
#    r.F = []
#    for i in range(a.n):
#        for j in range(a.na):
#            if a.e[i].f[j] != -1:
#                r.add_transition((i, a.e[i].f[j], A[j]))
#        if a.e[i].final:
#            r.F.append(i)
#    r.I = [a.i]
#    return r

cdef AutomatonToSageAutomaton(Automaton a, A):
    from sage.combinat.finite_state_machine import Automaton as SageAutomaton
    L = []
    if a.i == -1:
        I = []
    else:
        I = [a.i]
    F = []
    cdef int i, j
    for i in range(a.n):
        for j in range(a.na):
            if a.e[i].f[j] != -1:
                L.append((i, a.e[i].f[j], A[j]))
        if a.e[i].final:
            F.append(i)
    return SageAutomaton(L, initial_states=I, final_states=F)

cdef AutomatonToDiGraph(Automaton a, A, keep_edges_labels=True):
    from sage.graphs.digraph import DiGraph
    cdef int i, j
    cdef list L
    L = []
    for i in range(a.n):
        for j in range(a.na):
            if a.e[i].f[j] != -1:
                if keep_edges_labels:
                    L.append((i, a.e[i].f[j], A[j]))
                else:
                    L.append((i, a.e[i].f[j]))
    return DiGraph([range(a.n), L], loops=True, multiedges=True)


cdef class CAutomaton:
    """
    Class :class:`CAutomaton`, this class encapsulates a C structure
    for Automata and implement methods to manipulate them.

    INPUT:

    - ``a`` -- an automaton given as a list, a DiGraph, an Automaton, a DeterministicAutomaton or a CAutomaton

    - ``I`` -- list (default: None)
        The set of initial states.

    - ``F`` -- list (default: None)
        The set of final states.

    - ``A``-- list (default: None)
        The alphabet.

    - ``keep_S`` -- bool (default: True)
        Keep the labels of the states.

    - ``verb`` -- bool (default: False)
        Display informations for debugging.

    OUTPUT:

    Return a instance of :class:`CAutomaton`.

    EXAMPLES::

        sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
        sage: b = CAutomaton(a)
        sage: b
        CAutomaton with 4 states and an alphabet of 2 letters
        sage: c = Automaton({0:{1:'x',2:'z',3:'a'}, 2:{5:'o'}},initial_states=[0])
        sage: b = CAutomaton(c)
        sage: b
        CAutomaton with 5 states and an alphabet of 4 letters
        sage: d = DiGraph({0: [1,2,3], 1: [0,2]})
        sage: b = CAutomaton(d)
        sage: b
        CAutomaton with 4 states and an alphabet of 1 letter
        sage: g =  CAutomaton([(0,1,'a') ,(2,3,'b')])
        sage: g
        CAutomaton with 4 states and an alphabet of 2 letters

    """
    def __cinit__(self):
        sig_on()
        self.a = <NAutomaton *>malloc(sizeof(NAutomaton))
        sig_off()
        if self.a is NULL:
            raise MemoryError("Failed to allocate memory for "
                              "C initialization of CAutomaton.")
        self.a.e = NULL
        self.a.n = 0
        self.a.na = 0
        self.A = []
        self.S = None

    def __init__(self, a, I=None, F=None, A=None, keep_S=True, verb=False):
        """
        TESTS::

            sage: CAutomaton([(0, 1, 'a'), (2, 3, 'b')], I=[1])
            CAutomaton with 4 states and an alphabet of 2 letters
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=1)
            sage: a = CAutomaton(a)
            sage: a
            CAutomaton with 4 states and an alphabet of 2 letters
            sage: CAutomaton(a)
            CAutomaton with 4 states and an alphabet of 2 letters

        """
        cdef DeterministicAutomaton da
        cdef CAutomaton na
        if a is None:
            if verb:
                print("a is None")
            return
        from sage.graphs.digraph import DiGraph
        from sage.combinat.finite_state_machine import Automaton as SageAutomaton
        if isinstance(a, SageAutomaton):
            if verb:
                print("Sage automaton...")
            L = []
            for t in a.transitions():
                if len(t.word_in) != 1 or len(t.word_out) != 0:
                    print("Warning: informations lost during the conversion.")
                L.append((t.from_state, t.to_state, t.word_in[0]))
            if I is None:
                I = a.initial_states()
            if F is None:
                F = a.final_states()
            a = L
        if isinstance(a, list):
            if verb:
                print("List...")
            a = DiGraph(a, multiedges=True, loops=True)
        if isinstance(a, DiGraph):
            if verb:
                print("DiGraph...")
            SL = set(a.edge_labels())
            if A is None:
                A = list(SL)
            else:
                # test that A contains the labels
                if not SL.issubset(A):
                    raise ValueError("A label of a transition is not in the alphabet %s"%A)
            self.A = A
            if verb:
                print("getNAutomaton(%s, I=%s, F=%s, A=%s)" % (a, I, F, self.A))
            self.a[0] = getNAutomaton(a, I=I, F=F, A=self.A, verb=verb)
            if keep_S:
                self.S = a.S
        elif isinstance(a, CAutomaton):
            if verb:
                print("CAutomaton")
            na = a
            sig_on()
            self.a[0] = CopyNAutomaton(na.a[0], na.a.n, na.a.na)
            sig_off()
            self.A = na.A
            self.S = na.S
        elif isinstance(a, DeterministicAutomaton):
            if verb:
                print("DeterministicAutomaton")
            da = a
            self.a[0] = CopyN(da.a[0], verb=verb)
            self.A = da.A
        else:
            raise ValueError("Cannot convert the input to CAutomaton.")

    def __dealloc__(self):
        sig_on()
        if self.a is not NULL:
            FreeNAutomaton(self.a)
            free(self.a)
        sig_off()

    def __repr__(self):
        r"""
        Return a representation of the automaton.

        OUTPUT:

        A string.

        EXAMPLES::

            sage: CAutomaton(dag.AnyWord([0,1]))
            CAutomaton with 1 state and an alphabet of 2 letters

        """
        str = "CAutomaton with %d state"%self.a.n
        if self.a.n > 1:
            str += 's'
        str += " and an alphabet of %d letter"%self.a.na
        if self.a.na > 1:
            str += 's'
        return str

    def _latex_(self):
        r"""
        Return a latex representation of the automaton.

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = CAutomaton(a)
            sage: latex(b)     # optional -  dot2tex
            \documentclass{article}
            \usepackage[x11names, rgb]{xcolor}
            \usepackage[utf8]{inputenc}
            \usepackage{tikz}
            \usetikzlibrary{snakes,arrows,shapes}
            \usepackage{amsmath}
            %
            %
            <BLANKLINE>
            %
            <BLANKLINE>
            %
            <BLANKLINE>
            \begin{document}
            \pagestyle{empty}
            %
            %
            %
            <BLANKLINE>
            \enlargethispage{100cm}
            % Start of code
            % \begin{tikzpicture}[anchor=mid,>=latex',line join=bevel,]
            \begin{tikzpicture}[>=latex',line join=bevel,]
              \pgfsetlinewidth{1bp}
            %%
            \pgfsetcolor{black}
              % Edge: 0 -> 1
              \draw [-stealth'] (44.247bp,22.0bp) .. controls (54.848bp,22.0bp) and (67.736bp,22.0bp)  .. (89.697bp,22.0bp);
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw (67.0bp,33.0bp) node {a};
              % Edge: 2 -> 3
              \draw [-stealth'] (44.247bp,80.0bp) .. controls (54.848bp,80.0bp) and (67.736bp,80.0bp)  .. (89.697bp,80.0bp);
              \draw (67.0bp,91.0bp) node {b};
              % Node: 1
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (112.0bp,22.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (112.0bp,22.0bp) ellipse (22.0bp and 22.0bp);
              \draw (112.0bp,22.0bp) node {1};
            \end{scope}
              % Node: 0
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (22.0bp,22.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (22.0bp,22.0bp) ellipse (22.0bp and 22.0bp);
              \draw (22.0bp,22.0bp) node {0};
            \end{scope}
              % Node: 3
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (112.0bp,80.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (112.0bp,80.0bp) ellipse (22.0bp and 22.0bp);
              \draw (112.0bp,80.0bp) node {3};
            \end{scope}
              % Node: 2
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (22.0bp,80.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (22.0bp,80.0bp) ellipse (22.0bp and 22.0bp);
              \draw (22.0bp,80.0bp) node {2};
            \end{scope}
            %
            \end{tikzpicture}
            % End of code
            <BLANKLINE>
            %
            \end{document}
            %
            <BLANKLINE>
            <BLANKLINE>
            <BLANKLINE>
        """
        sx = 800
        sy = 600
        from sage.misc.latex import LatexExpr
        cdef char *file
        from sage.misc.temporary_file import tmp_filename
        file_name = tmp_filename()+".dot"
        file = file_name
        try:
            from dot2tex import dot2tex
        except ImportError:
            raise ModuleNotFoundError("dot2tex must be installed in order to have the LaTeX representation of the CAutomaton.\n\
                You can install it by doing './sage -i dot2tex' in a shell in the sage directory, or by doing 'install_package(package='dot2tex')' in the notebook.")
        cdef char** ll
        ll = <char **> malloc(sizeof(char*) * self.a.na)
        if ll is NULL:
            raise MemoryError("Failed to allocate memory for ll in "
                              "_latex_")
        cdef int i
        strA = []
        for i in range(self.a.na):
            strA.append(str(self.A[i]))
            ll[i] = strA[i]
        sig_on()
        NplotDot(file, self.a[0], ll, "Automaton", sx, sy, False)
        free(ll)
        sig_off()
        dotfile = open(file_name)
        return LatexExpr(dot2tex(dotfile.read()))

    def copy(self):
        """
        Do a copy of the :class:`CAutomaton`.

        OUTPUT:

        Return a copy of the :class:`CAutomaton`

        EXAMPLES::

            sage: a = CAutomaton([(0, 1, 'a'), (2, 3, 'b')], I=[0])
            sage: a.copy()
            CAutomaton with 4 states and an alphabet of 2 letters

        """
        r = CAutomaton(None)
        sig_on()
        r.a[0] = CopyNAutomaton(self.a[0], self.a.n, self.a.na)
        sig_off()
        r.A = self.A[:]
        return r

    @property
    def n_states(self):
        """
        return the numbers of states

        OUTPUT:

        return the numbers of states

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.n_states
            4
        """
        return self.a.n

    def n_successors(self, int i):
        """
        INPUT:

        - ``i`` -- int successor number

        OUTPUT:

        return the numbers of successor of state ``i``
        (i.e. the number of leaving transitions from state ``i``)

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.n_successors(0)
            1
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        return self.a.e[i].n

    def successor(self, int i, int j):
        """
        Give the state at end of the ``j`` th edge from the state ``i``

        INPUT:

        - ``i`` int state number
        - ``j`` int edge number

        OUTPUT:

        return the state at end of the ``j``th edge of the state ``i``

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.successor(1, 0)
            2
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        if j >= self.a.e[i].n or j < 0:
            raise ValueError("The state %s has no edge number %s !" % (i, j))
        return self.a.e[i].a[j].e

    def label(self, int i, int j):
        """
        Give the label of the ``j``th edge of the state ``i``

        INPUT:

        -``i`` -- int state number
        -``j`` -- int edge number

        OUTPUT:

        return the label index of the ``j``th edge of the state ``i``

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.label(1, 0)
            1
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        if j >= self.a.e[i].n or j < 0:
            raise ValueError("The state %s has no edge number %s !" % (i, j))
        return self.a.e[i].a[j].l

    def is_final(self, int i):
        """
        Return ``True`` if the state``i`` is  final, ``False`` otherwise.

        INPUT:

        -``i`` -- int state number

        OUTPUT:

        Return ``True`` if the state``i`` is  final, ``False`` otherwise.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.is_final(1)
            True
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        sig_on()
        ans = c_bool(self.a.e[i].final)
        sig_off()
        return ans

    def is_initial(self, int i):
        """
        Return `True`` if state ``i`` is initial, ``False`` otherwise.

        INPUT:

        -``i`` -- int state number

        OUTPUT:

        Return `True`` if state ``i`` is initial, ``False`` otherwise.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.is_initial(1)
            False

        TESTS::

            sage: a = CAutomaton(dag.AnyWord([0,1]))
            sage: a.is_initial(2)
            Traceback (most recent call last):
            ...
            ValueError: There is no state 2 !
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        sig_on()
        answ = c_bool(self.a.e[i].initial)
        sig_off()
        return answ

    @property
    def initial_states(self):
        """
        Return the set of initial states of the :class:`CAutomaton`

        OUTPUT:

        list

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = CAutomaton(a)
            sage: b.initial_states
            []
            sage: a = DeterministicAutomaton([(10,11,'a') ,(12,13,'b')], i=12)
            sage: b = CAutomaton(a)
            sage: b.initial_states
            [2]
            sage: a = DeterministicAutomaton([(10,11,'a') ,(12,13,'b')], i=0).mirror()
            sage: a.initial_states
            [0, 1, 2, 3]
        """
        cdef int i
        l = []
        for i in range(self.a.n):
            if self.a.e[i].initial:
                l.append(i)
        return l

    @property
    def final_states(self):
        """
        Return the set of final states

        OUTPUT:

        list

        EXAMPLES::

            sage: a = DeterministicAutomaton([(10,11,'a') ,(12,13,'b')])
            sage: b = CAutomaton(a)
            sage: b.final_states
            [0, 1, 2, 3]
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')], )
            sage: b = CAutomaton(a)
            sage: b.final_states
            [0, 1, 2, 3]
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
            sage: b = CAutomaton(a)
            sage: b.final_states
            [0, 3]
            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                (10,30,'z')], i=10)
            sage: a.final_states
            [0, 1, 2]
        """
        l = []
        for i in range(self.a.n):
            if self.a.e[i].final:
                l.append(i)
        return l

    @property
    def alphabet(self):
        """
        To get the :class:`CAutomaton` attribut alphabet
        Return the alphabet of the :class:`CAutomaton`

        OUTPUT:

        list

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = CAutomaton(a)
            sage: b.alphabet
            ['a', 'b']
            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                (10,30,'z')], i=10)
            sage: a.alphabet
            ['y', 'x', 'z']
        """
        return self.A

    def set_final(self, int i, bool final=True):
        """
        Set the state as a final/non-final state.

        INPUT:

        - ``i`` -- int - the index of the state
        - ``final`` -- (default: ``True``) - if True set this state as final set, otherwise set this state as non-final.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(10,11,'a') ,(12,13,'b')], i=10, final_states=[])
            sage: b = CAutomaton(a)
            sage: b.set_final(2)
            sage: b.final_states
            [2]

        TESTS::

            sage: b.set_final(6)
            Traceback (most recent call last):
            ...
            ValueError: 6 is not an index of a state (must be between 0 and 3)

        """
        if i < 0 or i >= self.a.n:
            raise ValueError("%s is not an index of a state (must be between "%i +
                             "0 and %s)" % (self.a.n - 1))
        self.a.e[i].final = final

    def set_initial_state(self, int i, bool initial=True):
        """
        Set the state as an initial/non initial state.

        INPUT:

        - ``i`` -- int - the index of the state
        - ``initial`` -- (default: ``True``) - if True set this state as initial set,
          otherwise set this state as non-initial.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = CAutomaton(a)
            sage: b.set_initial_state(2)
            sage: b.initial_states
            [2]
            sage: b.set_initial_state(6)
            Traceback (most recent call last):
            ...
            ValueError: i=6 is not the index of a state (i.e. between 0 and 3).

        """
        if i < 0 or i >= self.a.n:
            raise ValueError("i=%s is not the index of a state (i.e. between 0 and %s)."
                             % (i, self.a.n - 1))
        self.a.e[i].initial = initial

    def add_transition(self, int i, l, int f):
        """
        Add a transition from ``i`` to ``f`` labeled by ``l`` in the automaton.
        The states ``i`` and ``f`` must exist, and the label ``l`` must be in the alphabet.
        To add new states, use add_state(). To add letters in the alphabet use bigger_alphabet().

        INPUT:

        - ``i`` -- the first state
        - ``l`` -- the label of the edge
        - ``f`` -- the second state

        EXAMPLES::

            sage: a = DeterministicAutomaton([(10, 11, 'a'), (12, 13, 'b')], i=10)
            sage: b = CAutomaton(a)
            sage: b.add_transition(2,'a',1)
            sage: b.plot()  # random
            
        TESTS::
            
            sage: b.add_transition(2,'v',1)
            Traceback (most recent call last):
            ...
            ValueError: The letter v doesn't exist.
            sage: b.add_transition(2,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state  6 doesn't exist.
            sage: b.add_transition(5,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state  5 doesn't exist.

        """
        if i >= self.a.n or i < 0:
            raise ValueError("The state %s doesn't exist." % i)
        if f >= self.a.n or f < 0:
            raise ValueError("The state  %s doesn't exist." % f)
        try:
            k = self.A.index(l)
        except Exception:
            # La lettre %s n'existe pas.
            raise ValueError("The letter %s doesn't exist." % l)

        sig_on()
        AddTransitionN(self.a, i, f, k)
        sig_off()

    def add_state(self, bool final):
        """
        Add a state in the automaton

        INPUT:

        - ``final`` -- boolean indicate if the added state is final

        OUTPUT:

        return the numbers of states

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.add_state(True)  # not implemented
            TypeError                                 Traceback (most recent call last)
            ...
            TypeError: 'NotImplementedType' object is not callable
        """
        raise NotImplemented()

    def bigger_alphabet(self, list nA):
        """
        Increase the size of the alphabet

        INPUT:

        - ``nA`` -- list - new alphabet or new letters to add

        EXAMPLES::

            sage: a = CAutomaton(dag.AnyLetter([0, 1, 'a']))
            sage: a.alphabet
            ['a', 1, 0]
            sage: a.bigger_alphabet(['b', 'c'])
            sage: a.alphabet
            ['a', 1, 0, 'c', 'b']
        """
        self.A = self.A+list(set(nA).difference(self.A))

    def add_path(self, int e, int f, list li, verb=False):
        """
        Add a path labeled by the list ``li`` from
        state ``e`` to state ``f`` of :class:`CAutomaton`

        INPUT:

        - ``e`` -- int the starting state
        - ``f`` -- int the ending state
        - ``li`` -- list of labels
        - ``verb`` -- boolean (default: ``False``) if True, print debugging informations


        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a'), (2,3,'b')], i=2)
            sage: b = CAutomaton(a)
            sage: b.add_path(1, 2, ['a'])
            sage: b.plot()  # random

        """
        cdef int *l = <int *>malloc(sizeof(int)*len(li))
        if l is NULL:
            raise MemoryError("Failed to allocate memory for l in "
                              "add_path")
        for i in range(len(li)):
            try:
                l[i] = self.A.index(li[i])
            except:
                raise ValueError("The label %s is not in the alphabet! Use bigger_alphabet() if you want to add it"%li[i])
        sig_on()
        AddPathN(self.a, e, f, l, len(li), verb)
        free(l)
        sig_off()

    def determinize(self, sink=False, verb=0):
        """
        Compute a deterministic automaton that recognizes the same language as self.

        INPUT:

        - ``sink``  -- (default: ``False``) - give a complete automaton, with a sink state
        - ``verb`` -- boolean (default: ``False``) if True, print debugging informations

        OUTPUT:

        Return a non determinist automaton  :class:`CAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.determinize()
            DeterministicAutomaton with 2 states and an alphabet of 2 letters

            sage: a = CAutomaton([(0,0,0),(0,0,1),(0,1,1),(1,2,0),(1,2,1),(2,3,0),(2,3,1),(3,4,0),(3,4,1)], I=[0], F=[4])
            sage: a.determinize()
            DeterministicAutomaton with 16 states and an alphabet of 2 letters

        """
        cdef Automaton a
        r = DeterministicAutomaton(None)
        sig_on()
        a = DeterminizeN(self.a[0], sink, verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        return r

    def plot(self, file=None, int sx=10, int sy=8, verb=False):
        """
        plot a representation of the :class:`CAutomaton`.

        INPUT:

        - ``sx`` -- int (default: 10) - width given to dot command
        - ``sy`` -- int (default: 8) - height given to dot command
        - ``verb`` -- boolean (default: ``False``) if True, print debugging informations

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = CAutomaton(a)
            sage: b.plot()   # random

            sage: a = CAutomaton([(0,0,0),(0,0,1),(0,1,1),(1,2,0),(1,2,1),(2,3,0),(2,3,1),(3,4,0),(3,4,1)], I=[0], F=[4])
            sage: a.plot()   # random
            <PIL.PngImagePlugin.PngImageFile image mode=RGBA size=189x147 at 0x7F6711E442D0>
            
        .. PLOT::
           :width: 50%

            a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            b = CAutomaton(a)
            sphinx_plot(b)

        """
        cdef char** ll
        cdef int i
        cdef char *cfile
        cdef bool de
        if verb:
            print("Test if dot exists...")
        sig_on()
        de = DotExists ()
        sig_off()
        if de:
            if verb:
                print(" -> yes !")
                print("alloc ll (size=%s)..." % self.a.na)
            ll = <char **>malloc(sizeof(char*) * self.a.na)
            if ll is NULL:
                raise MemoryError("Failed to allocate memory for ll in "
                                  "plot")
            strA = []
            for i in range(self.a.na):
                strA.append(str(self.A[i]))
                ll[i] = strA[i]
            if file is None:
                from sage.misc.temporary_file import tmp_filename
                file = tmp_filename()
            cfile = file
            if verb:
                print("file=%s" % file)
            sig_on()
            NplotDot(cfile, self.a[0], ll, "Automaton", sx, sy, True)
            free(ll)
            sig_off()
            if verb:
                print("Ouvre l'image produite...")
            from PIL import Image
            return Image.open(file+'.png')
#            import ipywidgets as widgets
#            file = open(file+'.png', "rb")
#            image = file.read()
#            w = widgets.Image(
#                value=image,
#                format='png',
##                width=600,
##                height=400,
#            )
#            return w
        else:
            if verb:
                print(" -> no.")
            raise NotImplementedError("You cannot plot the CAutomaton without dot. Install the dot command of the GraphViz package.")


cdef class DeterministicAutomaton:
    r"""
    Class :class:`DeterministicAutomaton`, this class encapsulates a C structure
    for deterministic automata and implement methods to manipulate them.

    EXAMPLES::

        sage: DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
        DeterministicAutomaton with 4 states and an alphabet of 2 letters
        sage: d = DiGraph({0: [1,2,3], 1: [0,2]})
        sage: DeterministicAutomaton(d)
        Warning: the automaton was not deterministic! (edge 0 -None-> 1 and edge 0 -None-> 2)
        The result lost some informations.
        DeterministicAutomaton with 4 states and an alphabet of 1 letter
        sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'o'}})
        sage: DeterministicAutomaton(g)
        DeterministicAutomaton with 5 states and an alphabet of 4 letters
        sage: a = DeterministicAutomaton([(0, 1,'a') ,(2, 3,'b')], i = 2)
        sage: a
        DeterministicAutomaton with 4 states and an alphabet of 2 letters
        sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
        sage: a
        DeterministicAutomaton with 4 states and an alphabet of 2 letters
        sage: b = DeterministicAutomaton(a)
        sage: b
        DeterministicAutomaton with 4 states and an alphabet of 2 letters
        sage: c = Automaton({0:{1:'x',2:'z',3:'a'}, 2:{5:'o'}},initial_states=[0])
        sage: b = DeterministicAutomaton(c)
        sage: b
        DeterministicAutomaton with 5 states and an alphabet of 4 letters
        sage: dag.AnyWord([0, 1, 'a'])
        DeterministicAutomaton with 1 state and an alphabet of 3 letters
        sage: dag.Random(n=40, A=[None, -1, 1,2,3,'x','y','z'])  # random

    """

#    cdef Automaton* a
#    cdef list A

    def __cinit__(self):
        """

        """
        # print("cinit")
        self.a = <Automaton *>malloc(sizeof(Automaton))
        if self.a is NULL:
            raise MemoryError("Failed to allocate memory for "
                              "C initialization of DeterministicAutomaton.")
        # initialise
        self.a.e = NULL
        self.a.n = 0
        self.a.na = 0
        self.a.i = -1
        self.A = []
        self.S = None
        #self.dA = None
        #self.dS = None

    def __init__(self, a, i=None, final_states=None,
                 A=None, S=None, keep_S=True, verb=False):
        r"""
        INPUT:

        -``a`` - a list, a DiGraph, an Automaton, or a ``DeterministicAutomaton``

        - ``i`` - (default: None) - initial state

        - ``final_states`` - (default: None) - list of final states

        - ``A`` - (default: None) - alphabet

        - ``keep_S``- (default: `True`) - Keep labels of states or not

        OUTPUT:

        Return a instance of :class:`DeterministicAutomaton`.

        TESTS::

            sage: DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            DeterministicAutomaton with 4 states and an alphabet of 2 letters
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: DeterministicAutomaton(a)
            DeterministicAutomaton with 4 states and an alphabet of 2 letters
        """
        cdef DeterministicAutomaton da
        # print("init")
        if a is None:
            if A is not None:
                self.A = A
                self.a.na = len(A)
            if i is not None:
                print("Warning: Ignore i since there is no state.")
            if final_states is not None:
                print("Warning: Ignore final_states since there is no state.")
            if S is not None:
                print("Warning: Ignore S since there is no state.")
            return
        from sage.graphs.digraph import DiGraph
        from sage.combinat.finite_state_machine import Automaton as SageAutomaton
        if isinstance(a, CAutomaton):
            a = a.determinize().minimize()
        if isinstance(a, SageAutomaton):
            L = []
            for t in a.transitions():
                if len(t.word_in) != 1 or len(t.word_out) != 0:
                    print("Warning: informations lost during the conversion.")
                L.append((t.from_state, t.to_state, t.word_in[0]))
            if i is None:
                I = a.initial_states()
                if len(I) == 1:
                    i = I[0]
                else:
                    print("Warning: the automaton has several initial states.")
            if final_states is None:
                final_states = a.final_states()
            a = L
        if isinstance(a, list):
            if verb:
                print("list...")
            if S is None:
                a = DiGraph(a, multiedges=True, loops=True)
            else:
                a = DiGraph([S, a],
                            format='vertices_and_edges',
                            multiedges=True, loops=True)
        if isinstance(a, DiGraph):
            if verb:
                print("DiGraph...")
            SL = set(a.edge_labels())
            if A is None:
                A = list(SL)
            else:
                # test that A contains the labels
                if not SL.issubset(A):
                    raise ValueError("A label of a transition is not in the alphabet %s"%A)
            self.A = A
            self.a[0] = getAutomaton(a, initial=i, F=final_states, A=self.A)
            #self.dA = a.dA
            if keep_S:
                self.S = a.vertices()
                #self.dS = a.dS
        elif isinstance(a, DeterministicAutomaton):
            if verb:
                print("DeterministicAutomaton...")
            da = a
            sig_on()
            self.a[0] = CopyAutomaton(da.a[0], da.a.n, da.a.na)
            sig_off()
            self.A = da.A
            self.S = da.S
        else:
            raise ValueError("Cannot convert the input to DeterministicAutomaton.")

    def __dealloc__(self):
        """
        Desalloc  Automaton  Overwrite built-in function

        TESTS::

        """
        # print("free (%s States) "%self.a.n)
        sig_on()
        FreeAutomaton(self.a)
        # print("free self.a")
        free(self.a)
        sig_off()

    def copy(self):
        """
        Do a copy of the :class:`DeterministicAutomaton`.

        OUTPUT:

        Return a copy of the :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.copy()
            DeterministicAutomaton with 4 states and an alphabet of 2 letters
        
        TESTS::

            sage: a = dag.AnyLetter(['a', 'b'])
            sage: b = a.copy()
            sage: a == b
            True
            sage: a.alphabet.append(0)
            sage: a == b
            False
            sage: a = b.copy()
            sage: a.delete_state_in_place(0)
            sage: a == b
            False

        """
        cdef DeterministicAutomaton r 
        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = CopyAutomaton(self.a[0], self.a.n, self.a.na)
        sig_off()
        if self.S is not None:
            r.S = self.S[:]
        r.A = self.A[:]
        return r

    def __copy__(self):
        """
        Return a copy of the :class:`DeterministicAutomaton`.
        Overwrite builtin function.
        
        See copy() for more details.
        """
        return self.copy()

    def __repr__(self):
        """
        Return a representation of automaton,  Overwrite built-in function

        OUTPUT:

        Return a representation of automaton

        TESTS:

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: repr(a)
            'DeterministicAutomaton with 4 states and an alphabet of 2 letters'

        """
        str = "DeterministicAutomaton with %d state"%self.a.n
        if self.a.n > 1:
            str += 's'
        str += " and an alphabet of %d letter"%self.a.na
        if self.a.na > 1:
            str += 's'
        return str

    def _latex_(self):
        r"""
        Return a latex representation of the automaton.

        OUTPUT:

        string - latex representation of the automaton.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: latex(a)      # optional -  dot2tex
            \documentclass{article}
            \usepackage[x11names, rgb]{xcolor}
            \usepackage[utf8]{inputenc}
            \usepackage{tikz}
            \usetikzlibrary{snakes,arrows,shapes}
            \usepackage{amsmath}
            %
            %
            <BLANKLINE>
            %
            <BLANKLINE>
            %
            <BLANKLINE>
            \begin{document}
            \pagestyle{empty}
            %
            %
            %
            <BLANKLINE>
            \enlargethispage{100cm}
            % Start of code
            % \begin{tikzpicture}[anchor=mid,>=latex',line join=bevel,]
            \begin{tikzpicture}[>=latex',line join=bevel,]
              \pgfsetlinewidth{1bp}
            %%
            \pgfsetcolor{black}
              % Edge: 0 -> 1
              \draw [-stealth'] (44.247bp,22.0bp) .. controls (54.848bp,22.0bp) and (67.736bp,22.0bp)  .. (89.697bp,22.0bp);
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw (67.0bp,33.0bp) node {a};
              % Edge: 2 -> 3
              \draw [-stealth'] (44.247bp,80.0bp) .. controls (54.848bp,80.0bp) and (67.736bp,80.0bp)  .. (89.697bp,80.0bp);
              \draw (67.0bp,91.0bp) node {b};
              % Node: 1
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (112.0bp,22.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (112.0bp,22.0bp) ellipse (22.0bp and 22.0bp);
              \draw (112.0bp,22.0bp) node {1};
             \end{scope}
              % Node: 0
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (22.0bp,22.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (22.0bp,22.0bp) ellipse (22.0bp and 22.0bp);
              \draw (22.0bp,22.0bp) node {0};
            \end{scope}
              % Node: 3
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (112.0bp,80.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (112.0bp,80.0bp) ellipse (22.0bp and 22.0bp);
              \draw (112.0bp,80.0bp) node {3};
            \end{scope}
              % Node: 2
            \begin{scope}
              \definecolor{strokecol}{rgb}{0.0,0.0,0.0};
              \pgfsetstrokecolor{strokecol}
              \draw [solid] (22.0bp,80.0bp) ellipse (18.0bp and 18.0bp);
              \draw [solid] (22.0bp,80.0bp) ellipse (22.0bp and 22.0bp);
              \draw (22.0bp,80.0bp) node {2};
            \end{scope}
            %
            \end{tikzpicture}
            % End of code
            <BLANKLINE>
            %
            \end{document}
            %
            <BLANKLINE>
            <BLANKLINE>
            <BLANKLINE>

        """
        sx = 800
        sy = 600
        vlabels = None
        html = False
        verb = False
        from sage.misc.latex import LatexExpr
        cdef char *file
        from sage.misc.temporary_file import tmp_filename
        file_name = tmp_filename()+".dot"
        file = file_name
        try:
            from dot2tex import dot2tex
        except ImportError:
            print("dot2tex must be installed in order to have the LaTeX representation of the DeterministicAutomaton.")
            print("You can install it by doing './sage -i dot2tex' in a shell in the sage directory, or by doing 'install_package(package='dot2tex')' in the notebook.")
            return None
        cdef char** ll # labels of edges
        cdef char** vl # labels of vertices
        cdef int i
        ll = <char **>malloc(sizeof(char*) * self.a.na)
        if ll is NULL:
            raise MemoryError("Failed to allocate memory for ll in "
                              "_latex_")
        if vlabels is None:
            vl = NULL
        else:
            if verb:
                print("alloc %s..."%self.a.n)
            vl = <char **>malloc(sizeof(char*) * self.a.n)
            if vl is NULL:
                raise MemoryError("Failed to allocate memory for vl in "
                                  "_latex_")
            strV = []
            if verb:
                print("len %s %s" % (self.a.n, len(vlabels)))
            for i in range(self.a.n):
                if html:
                    strV.append("<" + vlabels[i] + ">")
                else:
                    strV.append("\"" + vlabels[i] + "\"")
                if verb:
                    print(strV[i])
                vl[i] = strV[i]
                if verb:
                    print("i=%s : %s" % (i, vl[i]))
        strA = []
        for i in range(self.a.na):
            strA.append(str(self.A[i]))
            ll[i] = strA[i]
        if verb:
            for i in range(self.a.n):
                print("i=%s : %s" % (i, vl[i]))
        if verb:
            print("plot...")
        sig_on()
        plotDot(file, self.a[0], ll, "Automaton", sx, sy, vl, html, verb, False)
        sig_off()
        if verb:
            print("free...plot")
        sig_on()
        free(ll)
        if vlabels is not None:
            free(vl)
        sig_off()
        dotfile = open(file_name)
        return LatexExpr(dot2tex(dotfile.read()))

    def __hash__(self):
        r"""
        Hash the automaton

        OUTPUT:

        Return a ``int``, hash code of the :class:`DeterministicAutomaton`

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: hash(a)
            761033288

            sage: a = dag.AnyLetter([0,1])
            sage: hash(a)
            1519588028
        """
        cdef int h
        h = hash(tuple(self.A))
        sig_on()
        h += hashAutomaton(self.a[0])
        sig_off()
        return h

    def string(self):
        r"""
        Return a ``string`` that can be evaluated to recover the DeterministicAutomaton.
        """
        S = self.states
        r = "DeterministicAutomaton([%s, [" % S
        c = 0
        for i in range(self.a.n):
            for j in range(self.a.na):
                k = self.a.e[i].f[j]
                if k != -1:
                    if c != 0:
                        r += ", "
                    r += "(%r, %r, %r)" % (S[i], S[k], self.A[j])
                    c += 1
        r += "]], A=%s, i=%r, final_states=%s)" % (self.A, S[self.a.i], self.final_states)
        return r

    def is_equal_to(self, DeterministicAutomaton other, verb=False):
        r"""
        Test if the two DetAutomata are equal.
        To be equal, the automata must have the same alphabet
        and exactly the same transitions
        (a permutation of the indices is not allowed).

        INPUT:

        - ``other`` -- other :class:`DeterministicAutomaton` to compare

        OUTPUT:

        Return the result of test (``True`` or ``False``)

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = DeterministicAutomaton([(0, 1, 'a'),(2,3,'b')], i=0)
            sage: a.is_equal_to(b)
            False
            sage: a.set_initial_state(0)
            sage: a.is_equal_to(b)
            True
        """
        if self.A != other.A:
            if verb:
                print("The alphabets are differents.")
            return False
        sig_on()
        r = equalsAutomaton(self.a[0], other.a[0], verb)
        sig_off()
        return c_bool(r)

    def __richcmp__(self, DeterministicAutomaton other, int op):
        r"""
        Compare function, Overwrite built-in function

        INPUT:

        - ``other`` -- other :class:`DeterministicAutomaton` to compare

        OUTPUT:

        Return the result of test (``True`` or ``False``)

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = DeterministicAutomaton([(0, 1, 'a'),(2,3,'b')], i=0)
            sage: a == b
            False
            sage: a.set_initial_state(0)
            sage: a == b
            True
            sage: a != b
            False
            sage: a < b
            Traceback (most recent call last):
            ...
            NotImplementedError: Comparaison <, >, <= or >= not implemented for DetAutomata.

        """
        from sage.structure.richcmp import (op_EQ, op_NE)
        # (rich_to_bool,
        # op_EQ, op_NE, op_LT, op_LE, op_GT, op_GE)
        cdef int r
        if op != op_EQ and op != op_NE:
            raise NotImplementedError("Comparaison <, >, <= or >= not implemented for DetAutomata.")
        sig_on()
        r = equalsAutomaton(self.a[0], other.a[0], False) and self.A == other.A
        sig_off()
        if op == op_EQ:
            return c_bool(r)
        else:
            return not c_bool(r)

    # give a Sage Automon from the DeterministicAutomaton
    def get_automaton(self):
        r"""
        Give a Sage Automon from the DeterministicAutomaton

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.get_automaton()
            Automaton with 4 states
        """
        return AutomatonToSageAutomaton(self.a[0], self.A)

    # give a Graph from the DeterministicAutomaton
    def get_digraph(self, keep_transitions_labels=True):
        r"""
        Give a DiGraph from the DeterministicAutomaton

        INPUT:

        - ``keep_transition_labels`` -- bool (default: ``True``)
            If false, return the graph without labels on edges.

        EXAMPLES::

            sage: a = dag.AnyLetter(['a', 'b'])
            sage: g = a.get_digraph()
            sage: g
            Looped multi-digraph on 2 vertices
            sage: g.edges()
            [(0, 1, 'a'), (0, 1, 'b')]

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.get_digraph()
            Looped multi-digraph on 4 vertices

        """
        return AutomatonToDiGraph(self.a[0], self.A, keep_edges_labels=keep_transitions_labels)

    def plot(self, int sx=10, int sy=8, vlabels=None,
             html=False, file=None, bool draw=True, verb=False):
        """
        Plot the :class:`DeterministicAutomaton`. Draw using the dot command, if installed on the platform.

        It is strongly recommanded to install the dot command of the Graphviz package in your system in order to get a nice picture.
        Otherwise it will draw using the function plot of :class:`Automaton` of Sage.

        INPUT:

        - ``sx`` - int (default: 10) - width of the picture
        - ``sy`` - int (default: 8) - height of the picture
        - ``vlabels`` - (default: None) - labels of the vertices
        - ``html`` - (default: ``False``) - tell if dot should draw vertices in html mode
        - ``file`` - (default: ``None``) - the address of the .dot file of the drawing (only if dot is installed)
        - ``draw`` - (default: ``True``) - if False, only generate the .dot file (only if dot is installed)
        - ``verb`` - (default: ``False``) - active or not the verbose mode

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.plot()  # random
            <PIL.PngImagePlugin.PngImageFile image mode=RGBA size=189x147 at 0x7F6711E442D0>

            sage: a = dag.Word(['a', 0, None])
            sage: a.plot()  # random

        TESTS::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.plot()  # random
            <PIL.PngImagePlugin.PngImageFile image mode=RGBA size=189x147 at 0x7FD4B6D94390>

        .. PLOT::
           :width: 50%

            a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sphinx_plot(a)
        """
        cdef char *cfile
        cdef char** ll = NULL # labels of edges
        cdef char** vl = NULL # labels of vertices
        cdef int i
        sig_on()
        de = DotExists()
        sig_off()
        if de:
            if file is None:
                from sage.misc.temporary_file import tmp_filename
                file = tmp_filename()+".dot"
            if verb:
                print("alloc %s..."%self.a.na)
            sig_on()
            ll = <char **>malloc(sizeof(char*) * self.a.na)
            sig_off()
            if ll is NULL:
                raise MemoryError("Failed to allocate memory for ll in "
                                  "plot")       
            if vlabels is None:
                if self.S is not None:
                    if verb:
                        print("alloc %s..." % self.a.n)
                    sig_on()
                    vl = <char **>malloc(sizeof(char*) * self.a.n)
                    sig_off()
                    if vl is NULL:
                        raise MemoryError("Failed to allocate memory for vl in "
                                          "plot")
                    strV = []
                    if html:
                        from sage.misc.html import html as htm
                    for i in range(self.a.n):
                        if html:
                            strV.append("<" + str(htm(self.S[i])) + ">")
                        else:
                            strV.append("\"" + str(self.S[i]) + "\"")
                        if verb:
                            print(strV[i])
                        vl[i] = strV[i]
                        if verb:
                            print("i=%s : %s" % (i, vl[i]))
                else:
                    if verb:
                        print("vl = NULL")
                    vl = NULL
            else:
                if verb:
                    print("alloc %s..." % self.a.n)
                vl = <char **>malloc(sizeof(char*) * self.a.n)
                if vl is NULL:
                    raise MemoryError("Failed to allocate memory for vl in "
                                      "plot")
                strV = []
                if verb:
                    print("len %s %s" % (self.a.n, len(vlabels)))
                for i in range(self.a.n):
                    if html:
                        strV.append("<" + vlabels[i] + ">")
                    else:
                        strV.append("\"" + vlabels[i] + "\"")
                    if verb:
                        print(strV[i])
                    vl[i] = strV[i]
                    if verb:
                        print("i=%s : %s" % (i, vl[i]))
            strA = []
            for i in range(self.a.na):
                strA.append(str(self.A[i]))
                ll[i] = strA[i]
            if verb and vl != NULL:
                for i in range(self.a.n):
                    print("i=%s : %s" % (i, vl[i]))
            if verb:
                print("plot...")
            sig_on()
            plotDot(file, self.a[0], ll, "Automaton", sx, sy, vl, html, verb, draw)
            sig_off()
            if verb:
                print("free...plot")
            if ll is not NULL:
                sig_on()
                free(ll)
                sig_off()
            if vl is not NULL: # vlabels is not None and self.S is not None:
                sig_on()
                free(vl)
                sig_off()
            if draw:
                from PIL import Image
                return Image.open(file+'.png')
#other possibility to get the image
#            import ipywidgets as widgets
#            file = open(file+'.png', "rb")
#            image = file.read()
#            w = widgets.Image(
#                value=image,
#                format='png',
##                width=600,
##                height=400,
#            )
#            return w
        else:
            return AutomatonToSageAutomaton(self.a[0], self.A).plot()

    @property
    def alphabet(self):
        """
        To get the :class:`DeterministicAutomaton` attribut alphabet

        OUTPUT:

        Return the alphabet ``A`` of the :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.alphabet
            ['a', 'b']

            sage: a = dag.Word(['a', 1, None])
            sage: a.alphabet
            ['a', 1, None]
        """
        return self.A

    def set_alphabet(self, list A):
        """
        Replace the alphabet by another one.

        INPUT:

        - ``A`` -- list of letters of alphabet

        EXAMPLES::

            sage: a = dag.Word(['a', 'b', 0])
            sage: a.alphabet
            ['a', 0, 'b']
            sage: a.set_alphabet([0, 1, 'c'])
            sage: b = dag.Word([0, 'c', 1])
            sage: a == b
            True

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_alphabet(['a', 'c', 'b'])
            sage: a.alphabet
            ['a', 'c']
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_alphabet(['a','e'])
            sage: a.alphabet
            ['a', 'e']

        TESTS::

            sage: a = dag.Word(['a', 'b', 'c'])
            sage: a.set_alphabet([0, 1])
            Traceback (most recent call last):
            ...
            ValueError: The size 2 of the new alphabet has to be the same (i.e. 3)

        """
        if len(A) < len(self.A):
            raise ValueError("The size %s of the new alphabet has to be the same (i.e. %s)"%(len(A), len(self.A)))
        self.A = A[:self.a.na]

    @property
    def initial_state(self):
        """
        Get the initial state of the :class:`DeterministicAutomaton`

        OUTPUT:

        Return an int, index of the initial state ``i`` of :class:`DeterministicAutomaton`,
        or -1 if there is no initial state.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.initial_state
            -1
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')], i=2)
            sage: a.initial_state
            2
        """
        return self.a.i

    def set_initial_state(self, int i):
        """
        Set the initial state.

        INPUT:

        - ``i`` -- int -- the initial state of the automaton

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_initial_state(2)
            sage: a.initial_state
            2
            sage: a.set_initial_state(6)
            Traceback (most recent call last):
            ...
            ValueError: The initial state must be a state of the automaton or -1: 6 is not in [-1, 3]
        """
        if i < self.a.n and i >= -1:
            self.a.i = i
        else:
            raise ValueError("The initial state must be a state of the automaton or -1: " +
                             "%d is not in [-1, %d]" % (i, self.a.n - 1))

    @property
    def final_states(self):
        """
        Indicate all final states

        OUTPUT:

        Return the list of final states

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.final_states
            [0, 1, 2, 3]
            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
            sage: a.final_states
            [0, 3]
            sage: a = dag.Word(['a', 'b', None])
            sage: a.final_states
            [3]
        """
        cdef int i
        l = []
        for i in range(self.a.n):
            if self.a.e[i].final:
                l.append(i)
        return l

    @property
    def states(self):
        """
        Indicate all states of the automaton

        OUTPUT:

        Return the list of states

        EXAMPLES::

            sage: a = dag.Word(['a', 'b', None])
            sage: a.states
            [0, 1, 2, 3]
        """
        if self.S is None:
            return range(self.a.n)
        else:
            return self.S

    def set_final_states(self, list lf):
        """
        Set the set of final states.

         INPUT:

        - ``lf`` -- list of states to set as final

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_final_states([0,3])
            sage: a.final_states
            [0, 3]
            sage: a.set_final_states([0,4])
            Traceback (most recent call last):
            ...
            ValueError: 4 is not a state!

        """
        cdef int f
        for f in range(self.a.n):
            self.a.e[f].final = 0
        for f in lf:
            if f < 0 or f >= self.a.n:
                raise ValueError("%d is not a state!" % f)
            self.a.e[f].final = 1

    def is_final(self, int e):
        """
        Indicate if the state is final

        INPUT:

        - ``e`` -- int input state to examine as final

        OUTPUT:

        ``True`` if the state ``e`` is final (i.e. ``False`` in the other case)

        EXAMPLES::

            sage: a = dag.Word(['a', 'b'])
            sage: a.is_final(1)
            False

            sage: a = DeterministicAutomaton([(10,11,'a') ,(12,13,'b')])
            sage: a.is_final(1)
            True
            sage: a.is_final(10)
            Traceback (most recent call last):
            ...
            ValueError: 10 is not a state!
        """
        if e >= 0 and e < self.a.n:
            sig_on()
            ans = c_bool(self.a.e[e].final)
            sig_off()
            return ans
        else:
            raise ValueError("%s is not a state!"%e)

    def set_final(self, int e, final=True):
        """
        Set the state as final.

         INPUT:

        - ``e`` -- int state to set as final
        - ``final`` -- (default: ``True``) set the state final or not final

        EXAMPLES::

            sage: a = dag.Word(['a', 'b'])
            sage: a.set_final(1)
            sage: a.final_states
            [1, 2]

            sage: a = DeterministicAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_final(2, False)
            sage: a.final_states
            [0, 1, 3]
            sage: a.set_final(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not the index of a state (i.e. it is not between 0 and 3)!
        """
        if e >= 0 and e < self.a.n:
            self.a.e[e].final = final
        else:
            raise ValueError("%d is not the index of a state (i.e. it is not between 0 and %s)!" % (e, self.a.n-1))

    def successor(self, int i, int j):
        """
        Return the state reached by following the transition with label of index j from state i.
        Return -1 if it doesn't exists.

        INPUT:

        - ``i`` - int - index of the state
        - ``j`` - int - index of the label

        OUTPUT:

        return a int, successor of the state i following edge j

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a'), (2,3,'b')])
            sage: a.successor(0, 1)
            -1
            sage: a.successor(2, 1)
            3
        """
        if i < 0 or i >= self.a.n or j < 0 or j >= self.a.na:
            return -1
        return self.a.e[i].f[j]

    def successors(self, int i):
        """
        Return indices of letters of leaving transitions from state ``i``.

        INPUT:

        - ``i`` -- int - the input state

        OUTPUT:

        return a list of int

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a'), (2,3,'b')])
            sage: a.successors(2)
            [1]
            sage: a.successors(4)
            []

        """
        cdef int j
        if i < 0 or i >= self.a.n:
            return []
        return [j for j in range(self.a.na) if self.a.e[i].f[j] != -1]

    def path(self, list l, i=None):
        """
        Follows the path labeled by a list ``l`` and return the reached state

        INPUT:

        - ``l`` -- list of indices of labels
        - ``i`` -- (default: ``None``) the initial state

        OUTPUT:

        int 
        return the state reached after following the path

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a'), (2,3,'b')], i=2)
            sage: a.path([1])
            3
            sage: a.path([0, 2])
            -1
        """
        if i is None:
            i = self.a.i
        for j in l:
            i = self.successor(i, j)
        return i

    def set_successor(self, int i, int j, int k):
        """
        Set the successor of the state i, following the transition labeled by j.
        In other words, add a transition from state i to state k labeled by A[j].
        If k is -1, this permits to remove the transition from state i labeled by A[j].

        INPUT:

        - ``i`` -- int - the input state
        - ``j`` -- int - the index of a letter
        - ``k`` -- int - the state we reach

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0,1,'a'), (2, 3,'b')], i=2)
            sage: a.set_successor(0, 1, 2)
            sage: a.successors(0)
            [0, 1]
            sage: a.set_successor(0, 4, 2)
            Traceback (most recent call last):
            ...
            ValueError: set_successor(0, 4, 2) : index out of bounds !

        """
        if i < 0 or i >= self.a.n or j < 0 or j >= self.a.na or k < -1 or k >= self.a.n:
            raise ValueError("set_successor(%s, %s, %s) : index out of bounds !" % (i, j, k))
        self.a.e[i].f[j] = k

    def zero_complete_in_place(self, z=None, verb=False):
        """
        Compute an automaton recognizing the language L(l*)^(-1), where L is
        the language of self and l is the letter of index z.
        This language L(l*)^(-1) is the set of words u such that there exists
        an integer n such that ul^n is in the language of self.

        INPUT:

        - ``z`` -- (default: ``None``) - index of the letter l
            If ``None``, take the index of 0.

        - ``verb`` -- (default: ``False``) if True, print debugging informations

        OUTPUT:

        Return a :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = dag.Word([0,1,0])
            sage: a.zero_complete_in_place()
            sage: a.final_states
            [2, 3]

        TESTS::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.zero_complete_in_place()
            Traceback (most recent call last):
            ...
            ValueError: 0 is not in list

        """
        if z is None:
            z = self.A.index(0)
        sig_on()
        ZeroComplete(self.a, z, verb)
        sig_off()

    def concatenate_with_zero_star(self, z=None, sink_state=False, simplify=True, verb=False):
        """
        Compute an automaton recognizing the language L(l*), where L is
        the language of self and l is the letter of index z.

        INPUT:

        - ``z`` -- (default: ``None``) - index of the letter l
            If ``None`` take the index of 0.

        - ``sink_state`` --  (default: ``False``) - give a result with or without sink state

        - ``simplify`` -- (default: ``True``) - prune and minimize the result
            (this removes the sink state)

        - ``verb`` -- (default: ``False``) if True, print debugging informations

        OUTPUT:

        Return a :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = dag.Word([0,1,0])
            sage: a.concatenate_with_zero_star()
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.concatenate_with_zero_star(z=0)
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.concatenate_with_zero_star(z=0, sink_state=True, simplify=False)
            DeterministicAutomaton with 5 states and an alphabet of 2 letters
            sage: b = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')])
            sage: b.concatenate_with_zero_star(True)
            DeterministicAutomaton with 1 state and an alphabet of 2 letters

        TESTS::

            sage: a = dag.Word(['a', 'b'])
            sage: a.concatenate_with_zero_star()
            Traceback (most recent call last):
            ...
            ValueError: 0 is not in list
        """
        cdef Automaton a
        if z is None:
            z = self.A.index(0)
        r = DeterministicAutomaton(None)
        sig_on()
        a = ZeroComplete2(self.a, z, sink_state, verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        if simplify:
            return r.prune().minimize()
        else:
            return r

    def left_concatenate_with_zero_star(self, z=None, simplify=True):
        """
        Compute an automaton recognizing the language (l*)L, where L is
        the language of self and l is the letter of index z.

        INPUT:

        - ``z`` - (default: ``None``) - index of the letter l
            If ``None``, take the index of letter 0.

        - ``simplify`` - (default: ``True``) - if True, prune and minimize the result

        OUTPUT:

        Return a :class:`DeterministicAutomaton` whose language is the set of words of
        the language of a with some zeroes on the beggining
        (we call zero the letter of index z).

        EXAMPLES::

            sage: a = dag.Word([0,1,0])
            sage: a.left_concatenate_with_zero_star()
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.left_concatenate_with_zero_star(0)
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: a.left_concatenate_with_zero_star(1)
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: b = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')])
            sage: b.left_concatenate_with_zero_star(1)
            DeterministicAutomaton with 1 state and an alphabet of 2 letters

        TESTS::

            sage: a = dag.Word(['a', 'b'])
            sage: a.left_concatenate_with_zero_star()
            Traceback (most recent call last):
            ...
            ValueError: 0 is not in list

        """
        cdef Automaton a
        if z is None:
            z = self.A.index(0)
        r = DeterministicAutomaton(None)
        sig_on()
        a = ZeroInv(self.a, z)  # list(self.A).index(self.A[z]))
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        if simplify:
            return r.prune().minimize()
        else:
            return r

    def difference(self, DeterministicAutomaton a, bool det=True, bool simplify=True):
        """
        Compute an automaton whose language is the set of differences of the two languages.

        INPUT:

        - ``a`` - DeterministicAutomaton - the automaton of the language that we substract
        - ``det`` - bool (default: ``True``) - determinize the result
        - ``simplify`` - bool (default: ``True``) - prune and minimize the result

        OUTPUT:

        Return a :class:`DeterministicAutomaton` if ``det`` is True, and return a :class:`Cautomaton` otherwise.

        EXAMPLES::

            sage: a = dag.AnyWord([0, 1])
            sage: a.difference(a)
            DeterministicAutomaton with 1 state and an alphabet of 3 letters

        TESTS::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.difference(a)
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for -: 'str' and 'str'
        """
        cdef DeterministicAutomaton r
        cdef CAutomaton nr
        cdef dict d
        r = self.product(a)
        d = {}
        for i in self.A:
            for j in a.A:
                d[(i, j)] = i-j
        return r.project(d, det=det, simplify=simplify)

    def prune_stop(self, verb=False):
        """
        Prune every stop state : remove all accessible states from which there no infinite way.

        INPUT:

        - ``verb`` -- boolean (default: ``False``) if True, print
          debugging informations

        OUTPUT:

        Return the pruned :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.prune_stop()
            DeterministicAutomaton with 0 state and an alphabet of 2 letters
            sage: b = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')])
            sage: b.prune_stop()
            DeterministicAutomaton with 0 state and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                (10,30,'z')], i=10)
            sage: a.prune_stop()
            DeterministicAutomaton with 3 states and an alphabet of 3 letters

            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                (10,30,'z')], i=10, final_states=[])
            sage: a.prune_stop()
            DeterministicAutomaton with 3 states and an alphabet of 3 letters


        TESTS::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.prune_stop(True)
            induction...
            States counter = 0
            count...
            cpt = 0
            final states...
            DeterministicAutomaton with 0 state and an alphabet of 2 letters

        """
        cdef Automaton a
        r = DeterministicAutomaton(None)
        sig_on()
        a = prune_inf(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        return r

    def prune_non_reachable(self, verb=False):
        """
        Prune the automaton:
        remove all non-reachable states

        INPUT:

        - ``verb`` -- boolean (default: ``False``) if True, print debugging informations

        OUTPUT:

        Return the pruned :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b'), (1, 3, 'b')], i=0)
            sage: a.prune_non_reachable()
            DeterministicAutomaton with 3 states and an alphabet of 2 letters
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b'), (1, 3, 'b')])
            sage: a.prune_non_reachable()
            DeterministicAutomaton with 0 state and an alphabet of 2 letters
            sage: b = DeterministicAutomaton([(0, 1, 'a'), (0, 3, 'b'), (1, 3, 'b')])
            sage: b.prune_non_reachable()
            DeterministicAutomaton with 0 state and an alphabet of 2 letters
        """
        if self.a.i == -1:
            return DeterministicAutomaton([], A=self.alphabet)
        cdef Automaton a
        r = DeterministicAutomaton(None)
        sig_on()
        a = pruneI(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        return r

    def prune(self, verb=False):
        """
        Prune the automaton:
        remove all non-reachable and non-co-reachable states

        INPUT:

        - ``verb`` -- boolean (default: ``False``) if True, print debugging informations

        OUTPUT:

        Return the pruned :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.prune()
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: b = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')])
            sage: b.prune()
            DeterministicAutomaton with 0 state and an alphabet of 2 letters

        """
        cdef Automaton a
        r = DeterministicAutomaton(None)
        sig_on()
        a = prune(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        r.S = None
        return r

    def product(self, DeterministicAutomaton b, dict d=None, bool simplify=True, verb=False):
        """
        Give the product of the :class:`DeterministicAutomaton` and ``a`` an other
        ``DeterministicAutomaton``. Assume the dictionnary ``d`` to be injective.

        INPUT:

        - ``a`` -- :class:`DeterministicAutomaton` to multiply

        - ``d`` -- dict (default: ``None``) - dictionary that associates
            something to couples of letters of each alphabet.
            The alphabet of the result is given by values of this dictionnary.
            If a couple of letters is not a key of this dictionnary,
            then corresponding edges of the product are avoided.

        - ``simplify`` -- bool (default: ``True``) - prune and minimize the result

        - ``verb`` -- boolean (default: ``False``) - if True,
            print debugging informations

        OUTPUT:

        Return the product as a :class:`DeterministicAutomaton`.

        EXAMPLES::

            sage: a = dag.Word(['a','b'])
            sage: b = dag.Word([0,1])
            sage: c = a.product(b)
            sage: c
            DeterministicAutomaton with 3 states and an alphabet of 4 letters
            sage: c.alphabet
            [('a', 0), ('a', 1), ('b', 0), ('b', 1)]

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'c'), (1, 2, 'd')], i=2)
            sage: a.product(b)
            DeterministicAutomaton with 1 state and an alphabet of 4 letters
        
        TESTS::

            sage: a = dag.Word(['a','b'])
            sage: b = dag.Word([0,1])
            sage: c = a.product(b)
            sage: d = dag.Word([('a',0),('b',1)], A = c.alphabet)
            sage: c.has_same_language_as(d)
            True
            
            sage: a = dag.Word(['a', 'b'])
            sage: d = {('a', 'a'):0, ('b', 'b'):0}
            sage: a.product(a, d)
            Traceback (most recent call last):
            ...
            ValueError: The dictionnary is not injective!


            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'c'), (1, 2, 'd')], i=2)
            sage: a.product(b, verb =True)
            {('b', 'c'): ('b', 'c'), ('a', 'd'): ('a', 'd'), ('a', 'c'): ('a', 'c'), ('b', 'd'): ('b', 'd')}
            Av=[('a', 'c'), ('a', 'd'), ('b', 'c'), ('b', 'd')]
            dv={('b', 'c'): 2, ('a', 'd'): 1, ('a', 'c'): 0, ('b', 'd'): 3}
            {'a': 0, 'b': 1}
            {'c': 0, 'd': 1}
            Keys=[('b', 'c'), ('a', 'd'), ('a', 'c'), ('b', 'd')]
            dC=
            [ 0 2 1 3 ]
            Automaton with 4 states, 2 letters.
            0 --0--> 1
            2 --1--> 3
            initial state 0.
            Automaton with 3 states, 2 letters.
            0 --1--> 1
            2 --0--> 1
            initial state 1.
            DeterministicAutomaton with 1 state and an alphabet of 4 letters
            sage: b = DeterministicAutomaton([(3, 2, 'c'), (1, 2, 'd')])
            sage: a.product(b)
            DeterministicAutomaton with 1 state and an alphabet of 4 letters

        """
        if d is None:
            d = {}
            for la in self.A:
                for lb in b.A:
                    d[(la, lb)] = (la, lb)
            if verb:
                print(d)
        else:
            V = d.values()
            if len(V) != len(set(V)):
                raise ValueError("The dictionnary is not injective!")
        cdef Automaton a
        cdef Dict dC
        r = DeterministicAutomaton(None)
        Av = []
        dv = imagProductDict(d, self.A, b.A, Av=Av)
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        dC = getProductDict(d, self.A, b.A, dv=dv, verb=verb)
        if verb:
            print("dC=")
            sig_on()
            printDict(dC)
            sig_off()
        sig_on()
        a = Product(self.a[0], b.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        r.a[0] = a
        r.A = Av
        if simplify:
            return r.prune().minimize()
        else:
            return r

    def intersection(self, DeterministicAutomaton a, simplify=True, verb=False):
        """
        Give a automaton recognizing the intersection of the languages
        of ``self`` and ``a``.

        INPUT:

        - ``a`` -- :class:`DeterministicAutomaton` to intersect

        - ``simplify`` - (default: ``True``) - if True,
            prune and minimize the result

        - ``verb`` -- boolean (default: ``False``) if True,
            print debugging informations

        OUTPUT:

        Return the intersected :class:`DeterministicAutomaton`.

        EXAMPLES::

            #. Intersection of words that contains 'aa' and words that contains 'bb'

                sage: a = dag.AnyWord(['a','b']).concatenate(dag.Word(['a', 'a'])).concatenate(dag.AnyWord(['a','b']))
                sage: b = dag.AnyWord(['a','b']).concatenate(dag.Word(['b', 'b']).concatenate(dag.AnyWord(['a','b'])))
                sage: a.intersection(b)
                DeterministicAutomaton with 8 states and an alphabet of 2 letters

            #. Intersection of two languages whose alphabets are different

                sage: a = dag.AnyWord(['a','b','c'])
                sage: b = dag.AnyWord(['a','b','e'])
                sage: a.intersection(b)
                DeterministicAutomaton with 1 state and an alphabet of 2 letters

            sage: a = dag.AnyWord(['a', 'b'])
            sage: b = dag.Word(['a','b','a'])
            sage: a.intersection(b).has_same_language_as(b)
            True

        TESTS::

            sage: a = dag.AnyWord(['a','b','c'])
            sage: b = dag.AnyWord(['a','b','e'])
            sage: a.intersection(b) == dag.AnyWord(['a', 'b'])
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'c'), (1, 2, 'd')], i=2)
            sage: a.intersection(b)
            DeterministicAutomaton with 1 state and an alphabet of 0 letter
            sage: a.intersection(b, simplify=False)
            DeterministicAutomaton with 12 states and an alphabet of 0 letter

        """
        cdef dict d
        cdef DeterministicAutomaton p
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l, l)] = l
        if verb:
            print("d=%s" % d)
        p = self.product(a, d, simplify=simplify, verb=verb)
        if simplify:
            return p.prune().minimize()
        else:
            return p

    def is_complete(self):
        """
        Determine if the automaton is complete
        (i.e. every state has leaving transitions
        labeled by every letter of the alphabet)

        OUTPUT:

        Return ``True`` if the automaton is complete, ``False`` otherwise.

        EXAMPLES::

            sage: a = dag.AnyWord(['a','b'])
            sage: a.is_complete()
            True
            sage: a = dag.Word(['a','b'])
            sage: a.is_complete()
            False
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.is_complete()
            False
            sage: a = DeterministicAutomaton([(0, 0, 'a')])
            sage: a.is_complete()
            True
        """
        cdef bool res
        sig_on()
        res = IsCompleteAutomaton(self.a[0])
        answ = c_bool(res)
        sig_off()
        return answ

    def complete_in_place(self, label_sink='s'):
        """
        Complete the automaton in-place, by adding a sink state if necessary.

        INPUT:

        - ``label_sink`` -- (default: ``s``) - label of the sink state
          if added and if self has labels of states.

        OUTPUT:

        Return ``True`` if the automaton was not complete
        (a sink state has been added), return ``False`` otherwise.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complete_in_place()
            True
            sage: a
            DeterministicAutomaton with 5 states and an alphabet of 2 letters

            sage: dag.AnyWord(['a','b']).complete_in_place()
            False

        TESTS::

            sage: a = dag.Word(['a','b','a'])
            sage: a.complete_in_place()
            True
            sage: a.has_same_language_as(dag.Word(['a','b','a']))
            True
            sage: a == dag.Word(['a','b','a'])
            False

        """
        sig_on()
        res = CompleteAutomaton(self.a)
        res = c_bool(res)
        sig_off()
        if res:
            if self.S is not None:
                # add a label for the sink state
                self.S.append(label_sink)
        return res

    def prefix(self, w, i=None, int algo=1):
        """
        Give an automaton recognizing the language w(w^(-1)L) where L is the language of self.
        It is the set of words recognized by self and starting with word w.

        INPUT:

        - ``w`` -- a word

        - ``i`` -- int (default: ``None``) - the initial state used

        - ``algo`` -- int (default: ``1``)
          The algorithm used for the computation.

        EXAMPLES::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.prefix(['a', 'a'])
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.prefix(['a', 'a'], algo=2)
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

        TESTS::

            sage: a = dag.Word(['a', 'b', 'a'])
            sage: a.prefix(['c'])
            Traceback (most recent call last):
            ...
            ValueError: 'c' is not in the alphabet ['a', 'b']

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.prefix(['a', 'a']).has_same_language_as(a.prefix(['a', 'a'], algo=2))
            True

        """
        cdef DeterministicAutomaton a
        cdef int* l
        cdef int j
        if algo == 2:
            l = <int*>malloc(sizeof(int)*self.a.n)
            for j in range(len(w)):
                try:
                    l[j] = self.A.index(w[j])
                except:
                    raise ValueError("letter %s of the word %s is not in the alphabet of self"%(w[j], w))
            if i is None:
                i = self.a.i
            a = DeterministicAutomaton(None)
            sig_on()
            a.a[0] = PieceAutomaton(self.a[0], l, len(w), i)
            free(l)
            sig_off()
            a.A = self.A
            return a
        else:
            a = self.copy()
            a.shift_list_in_place(w)
            return dag.Word(w).concatenate(a)

    def prefix_closure(self):
        """
        Give an automaton recognizing the smallest language stable
        by prefix and containing the language of self
        i.e. every states after pruning begins final

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.prefix_closure()
            DeterministicAutomaton with 2 states and an alphabet of 2 letters

            sage: b = dag.Word(['a','b','a'])
            sage: c = b.prefix_closure()
            sage: c
            DeterministicAutomaton with 4 states and an alphabet of 2 letters
            sage: dag.Word(['a', 'b']).included(c)
            True
        """
        cdef int i
        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = prune(self.a[0], False)
        sig_off()
        r.A = self.A
        for i in range(r.a.n):
            r.a.e[i].final = True
        return r

    # Use epsilon-transitions rather than product ?
    def union(self, DeterministicAutomaton a, simplify=True, verb=False):
        """
        Return an automaton recognizing the union of the two languages.
        Warning: there is a side effect, the automata are completed.

        INPUT:

        - ``a`` -- :class:`DeterministicAutomaton`

        - ``simplify`` --  (default: ``True``)
          prune and minimize the result ?

        - ``verb`` -- boolean (default: ``False``) if True,
          print debugging informations

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            #. Union of words that contains 'aa' and words that contains 'bb'

                sage: a = dag.AnyWord(['a','b']).concatenate(dag.Word(['a', 'a'])).concatenate(dag.AnyWord(['a','b']))
                sage: b = dag.AnyWord(['a','b']).concatenate(dag.Word(['b', 'b']).concatenate(dag.AnyWord(['a','b'])))
                sage: a.union(b)
                DeterministicAutomaton with 4 states and an alphabet of 2 letters

            #. Union of two languages whose alphabets are differents

                sage: a = dag.Word(['a','b','c'])
                sage: b = dag.Word(['a','b','a'])
                sage: a.union(b)
                DeterministicAutomaton with 4 states and an alphabet of 3 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.union(b)
            DeterministicAutomaton with 2 states and an alphabet of 3 letters
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')])
            sage: a.union(b)
            DeterministicAutomaton with 2 states and an alphabet of 3 letters

        TESTS::

            sage: a = dag.AnyWord(['a','b']).concatenate(dag.Word(['a', 'a'])).concatenate(dag.AnyWord(['a','b']))
            sage: b = dag.AnyWord(['a','b']).concatenate(dag.Word(['b', 'b']).concatenate(dag.AnyWord(['a','b'])))
            sage: c = DeterministicAutomaton([(0, 0, 'a'), (0, 0, 'b'), (1, 2, 'a'), (1, 3, 'b'), (2, 0, 'a'), (2, 3, 'b'), (3, 2, 'a'), (3, 0, 'b')], i=1, final_states=[0])
            sage: a.union(b) == c
            True

            sage: a = dag.Word(['a','b','c'])
            sage: b = dag.Word(['a','b','a'])
            sage: a.included(a.union(b))
            True
            sage: b.included(a.union(b))
            True

        """
        cdef DeterministicAutomaton a1 = self
        cdef DeterministicAutomaton a2 = a

        # increase the alphabets if necessary
        if set(a1.A) != set(a2.A):
            A = list(set(a1.A+a2.A))
            a1 = a1.bigger_alphabet(A)
            a2 = a2.bigger_alphabet(A)
        else:
            a1 = a1.copy()  # in order to avoid side effect
            a2 = a2.copy()

        # complete the automata
        sig_on()
        CompleteAutomaton(a1.a)
        CompleteAutomaton(a2.a)
        sig_off()

        # make the product
        d = {}
        for l in a1.A:
            d[(l, l)] = l

        cdef Automaton ap
        cdef Dict dC
        r = DeterministicAutomaton(None)
        Av = []
        dv = imagProductDict(d, a1.A, a2.A, Av=Av)
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        dC = getProductDict(d, a1.A, a2.A, dv=dv, verb=verb)
        if verb:
            print("dC=")
            sig_on()
            printDict(dC)
            sig_off()
        sig_on()
        ap = Product(a1.a[0], a2.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        # set final states
        cdef int i, j
        cdef n1 = a1.a.n
        for i in range(n1):
            for j in range(a2.a.n):
                ap.e[i + n1 * j].final = a1.a.e[i].final or a2.a.e[j].final

        r.a[0] = ap
        r.A = Av
        if simplify:
            return r.prune().minimize()
        else:
            return r

    def split(self, DeterministicAutomaton a, bool simplify=True, bool verb=False):
        """
        Split the automaton with respect to a :class:`DeterministicAutomaton` ``a``.
        Return two DeterministicAutomaton recognizing the intersection of the language
        of self with the one of a and with the complementary of the language
        of a.
        We assume that the two automata have the same alphabet, otherwise the
        complementary is taken in the set of words over the intersection of
        the two alphabets.

        INPUT:

        - ``a`` -- :class:`DeterministicAutomaton` - we split ``self`` with respect
          to this automaton.

        - ``simplify`` -- bool (default: True) - if True, prune and
          minimize the result.

        - ``verb`` -- bool (default: False) - if True, display
          informations for debugging.

        OUTPUT:

        Return tuple of two splited automaton, recognizing respectively
        the language intersection of L and La, and the language
        intersection of L and complementary of La,
        where L is the language of self and La is the language of a.

        EXAMPLES::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: b = a.concatenate(dag.Word(['a', 'a'])).concatenate(a)
            sage: a.split(b)
            [DeterministicAutomaton with 3 states and an alphabet of 2 letters,
             DeterministicAutomaton with 2 states and an alphabet of 2 letters]
            sage: a = dag.AnyWord(['a', 'b'])
            sage: b = a.concatenate(dag.Word(['a', 'a'])).concatenate(a)
            sage: c = a.concatenate(dag.Word(['b', 'b'])).concatenate(a)
            sage: c.split(b)
            [DeterministicAutomaton with 8 states and an alphabet of 2 letters,
             DeterministicAutomaton with 5 states and an alphabet of 2 letters]
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.split(b)
            [DeterministicAutomaton with 1 state and an alphabet of 1 letter,
             DeterministicAutomaton with 2 states and an alphabet of 1 letter]
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')])
            sage: a.split(b)
            [DeterministicAutomaton with 1 state and an alphabet of 1 letter,
             DeterministicAutomaton with 2 states and an alphabet of 1 letter]
            sage: print(a)
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

        TESTS::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: b = dag.AnyWord(['c', 'a'])
            sage: a.split(b)
            [DeterministicAutomaton with 1 state and an alphabet of 1 letter,
             DeterministicAutomaton with 1 state and an alphabet of 1 letter]
            sage: print(a)
            DeterministicAutomaton with 1 state and an alphabet of 2 letters
            
        .. WARNING::

            There is a side-effect: the automaton ``a`` is completed.

        """
        cdef Automaton ap, ap2
        cdef Dict dC
        cdef DeterministicAutomaton r, r2
        cdef list Av

        # complete the automaton a
        sig_on()
        CompleteAutomaton(a.a)
        sig_off()
        # make the product
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l, l)] = l

        r = DeterministicAutomaton(None)
        r2 = DeterministicAutomaton(None)
        Av = []
        dv = imagProductDict(d, self.A, a.A, Av=Av)
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        dC = getProductDict(d, self.A, a.A, dv=dv, verb=verb)
        if verb:
            print("dC=")
            sig_on()
            printDict(dC)
            sig_off()
        sig_on()
        ap = Product(self.a[0], a.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        # set final states for the intersection
        cdef int i, j
        cdef n1 = self.a.n
        for i in range(n1):
            for j in range(a.a.n):
                ap.e[i+n1*j].final = self.a.e[i].final and a.a.e[j].final

        # complementary of a in self
        sig_on()
        ap2 = CopyAutomaton(ap, ap.n, ap.na)
        sig_off()
        # set final states
        for i in range(n1):
            for j in range(a.a.n):
                ap2.e[i+n1*j].final = self.a.e[i].final and not a.a.e[j].final

        r.a[0] = ap
        r.A = Av
        r2.a[0] = ap2
        r2.A = Av
        if simplify:
            return [r.prune().minimize(), r2.prune().minimize()]
        else:
            return [r, r2]

    def shift_in_place(self, l, int np=1, verb=False):
        """
        Shift the automaton in-place to recognize the language shifted ``np``
        times by the letter l.
        The new language is the language of words u such that (l^np)u 
        was recognized by self.

        INPUT:

        - ``l`` -- letter to shift
        - ``np`` -- int (default: ``1``) - number of time we shift
        - ``verb`` -- boolean (default: ``False``) - if True, print
          debugging informations

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.shift_in_place('a')
            sage: a.initial_state
            1
            sage: a.shift_in_place('a', 2)
            sage: a.initial_state
            -1

        TESTS::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.shift_in_place('c')
            Traceback (most recent call last):
            ...
            ValueError: 'c' is not a letter of the alphabet ['a', 'b']

        """
        cdef int i
        try:
            l = self.A.index(l)
        except:
            raise ValueError("%r is not a letter of the alphabet %s" % (l, self.A))
        for i in range(np):
            if self.a.i != -1:
                self.a.i = self.a.e[self.a.i].f[l]

    def shift_list_in_place(self, list l):
        """
        Shift the automaton in-place to recognize the language shifted by l (list of letters).
        The new language is the language of words u such that lu 
        was recognized by self.

        INPUT:

        - ``l`` -- list - list of letter to shift

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.shift_list_in_place(['a'])
            sage: a.initial_state
            1

        TESTS::

            sage: a = dag.Word(['a','b','a','b'])
            sage: a.shift_list_in_place(['a','b'])
            sage: a.simplify().has_same_language_as(dag.Word(['a','b']))
            True

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.shift_list_in_place(['a','c'])
            Traceback (most recent call last):
            ...
            ValueError: 'c' is not in the alphabet ['a', 'b']

        """
        for i in l:
            try:
                i = self.A.index(i)
            except:
                raise ValueError("%r is not in the alphabet %s" % (i, self.A))
            if self.a.i != -1:
                self.a.i = self.a.e[self.a.i].f[i]

    def unshift(self, l, int np=1, final=False):
        """
        Unshift the automaton to recognize the language shifted
        ``np`` times by letter ``l``.
        The new language is the languages of words (l^np)u,
        where u is recognized by self and l^np is
        the letter l repeated np times.

        INPUT:

        - ``l`` -- letter to shift

        - ``np``  --  int (default: ``1``)
          Number of times we shift

        - ``final`` -- boolean (default: ``False``)
          if True, the empty word is added to the language

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.unshift('a')
            DeterministicAutomaton with 5 states and an alphabet of 2 letters
            sage: a.unshift('c')
            DeterministicAutomaton with 5 states and an alphabet of 3 letters
            sage: a.unshift('a', 4)
            DeterministicAutomaton with 8 states and an alphabet of 2 letters
        """
        cdef DeterministicAutomaton r
        cdef int i
        cdef int ne

        if np <= 0:
            raise ValueError("np must be positive!")

        r = DeterministicAutomaton(None)
        try:
            l = self.A.index(l)
            r.A = self.A[:]
            sig_on()
            r.a[0] = CopyAutomaton(self.a[0], self.a.n+np, self.a.na)
            sig_off()
        except:
            r.A = self.A+[l]
            l = len(self.A) # new letter
            sig_on()
            r.a[0] = CopyAutomaton(self.a[0], self.a.n+np, self.a.na+1)
            sig_off()
        ne = self.a.n # first new state
        for j in range(np):
            for i in range(r.a.na):
                r.a.e[ne+j].f[i] = -1
            if j > 0:
                r.a.e[ne+j].f[l] = ne+j-1
            else:
                r.a.e[ne+j].f[l] = self.a.i
            r.a.e[ne+j].final = final
            sig_check()
        r.a.i = ne+np-1
        return r

    # this function could be written in a more efficient way
    def unshift_list(self, list l):
        """
        Return a new automaton whose language is the set of words wu,
        where u is recognized by self, and w is the word
        corresponding to the list of indices of letters l.

        INPUT:

        - ``l`` -- list of indices of letters

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.unshift_list(['a', 'b'])
            DeterministicAutomaton with 6 states and an alphabet of 2 letters
            sage: a.unshift_list(['c', 'd'])
            DeterministicAutomaton with 6 states and an alphabet of 4 letters

            sage: a = dag.AnyWord([0,1])
            sage: a.unshift_list([1,0,1])
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

        """
        a = self
        l.reverse()
        for i in l:
            a = a.unshift(i)
        l.reverse()
        return a

    def get_cautomaton(self, verb=False):
        """
        Convert a determinist automaton :class:`DeterministicAutomaton` to
        a non determinist automaton :class:`CAutomaton`

        INPUT:

        - ``verb`` -- boolean (default: ``False``) if True,
          print debugging informations

        OUTPUT:

        Return a :class:`CAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.get_cautomaton()
            CAutomaton with 4 states and an alphabet of 2 letters

            sage: CAutomaton(dag.AnyWord(['a','b']))
            CAutomaton with 1 state and an alphabet of 2 letters

        TESTS::

            sage: a = dag.Word(['a','b'])
            sage: CAutomaton(a)
            CAutomaton with 3 states and an alphabet of 2 letters

            sage: a = dag.Word(['a','b'])
            sage: a.get_cautomaton().determinize() == a
            True

        """
        cdef NAutomaton a
        cdef CAutomaton r

        r = CAutomaton(None)
        sig_on()
        a = CopyN(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A[:]
        return r

    def concatenate(self, DeterministicAutomaton b, det=True, simplify=True, verb=False):
        """
        Return an automaton recognizing the concatenation of the
        languages of self and ``b``.

        INPUT:

        - ``b`` -- :class:`DeterministicAutomaton`  to concatenate to self

        - ``det``  -- boolean (default: ``True``) - if True, determinize

        - ``simplify`` -- boolean (default: ``True``) - if True and if det=True,
          prune and minimize

        - ``verb`` -- boolean (default: ``False``) if True, print
          debugging informations

        OUTPUT:

        Return a :class:`CAutomaton` (if ``det`` is ``False``)
        or  :class:`DeterministicAutomaton` (if ``det`` is ``True``)

        EXAMPLES::

            #. Words starting by 'aa'
                sage: dag.Word(['a','a']).concatenate(dag.AnyWord(['a','b']))
                DeterministicAutomaton with 3 states and an alphabet of 2 letters

            #. words containing 'aa'

                sage: a = dag.AnyWord(['a','b'])
                sage: a.concatenate(dag.Word(['a','a'])).concatenate(a)
                DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.concatenate(b)
            DeterministicAutomaton with 2 states and an alphabet of 3 letters
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.concatenate(b, det=False)
            CAutomaton with 7 states and an alphabet of 3 letters
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')])
            sage: a.concatenate(b)
            DeterministicAutomaton with 1 state and an alphabet of 3 letters

        TESTS::

            sage: a = dag.Word(['a', 'b'])
            sage: b = dag.Word(['a', 'c'])
            sage: a.concatenate(b).has_same_language_as(dag.Word(['a','b','a','c']))
            True

        """
        cdef DeterministicAutomaton a
        cdef NAutomaton na
        cdef CAutomaton r
        cdef DeterministicAutomaton r2

        if self.A != b.A:
            A = list(set(self.A).union(set(b.A)))
            if verb:
                print("Alphabet Changing (%s, %s -> %s)..." % (self.A, b.A, A))
            a = self.bigger_alphabet(A)
            b = b.bigger_alphabet(A)
        else:
            a = self
            A = self.A
        if verb:
            print("a=%s (A=%s)\nb=%s (A=%s)" % (a, a.A, b, b.A))
        r = CAutomaton(None)
        sig_on()
        na = Concat(a.a[0], b.a[0], verb)
        sig_off()
        r.a[0] = na
        r.A = A

        if det:
            if verb:
                print("Determinize...")
            r2 = r.determinize()
            if simplify:
                if verb:
                    print("Prune and minimize...")
                return r2.prune().minimize()
            else:
                return r2
        else:
            return r

    def project(self, dict d, bool det=True, bool simplify=True, bool verb=False):
        """
        Project with respect to the dictionary ``d``.
        Give an automaton where labels are replaced according to ``d``.

        INPUT:

        - ``d`` -- dictionary used for the projection

        - ``det``  --  bool (default: ``true``) - determinize the result or not

        - ``simplify`` -- bool (default: ``True``) - if True and if
          det=True, prune and minimize

        - ``verb`` -- boolean (default: ``False``) - activate or
          desactivate the verbose mode

        OUTPUT:

        Return a :class:`CAutomaton` (when ``det``=``False``)
        or a :class:`DeterministicAutomaton` (when ``det``=``True``)

        EXAMPLES::

            sage: a = dag.Word(['a','b'])
            sage: b = a.product(a)
            sage: d = {('a','a'):'a', ('a','b'):'a', ('b','a'):'b', ('b','b'):'b'}
            sage: b.project(d)
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = {'a' : 0, 'b': 0, 'c':1, 'd':1}
            sage: a.project(d)
            DeterministicAutomaton with 2 states and an alphabet of 1 letter

            sage: a = dag.Word([(0,1), (1,0), (1,1)])
            sage: d = {}
            sage: for i in range(2):
            ....:     for j in range(2):
            ....:         d[(i,j)] = i-j
            sage: a.project(d).has_same_language_as(dag.Word([-1, 1, 0]))
            True

        TESTS::

            sage: a = dag.Word(['a','b'])
            sage: b = a.product(a)
            sage: d = {('a','a'):'a', ('a','b'):'a', ('b','a'):'b', ('b','b'):'b'}
            sage: a.has_same_language_as(b.project(d))
            True

        """
        cdef NAutomaton a
        cdef Dict dC
        cdef CAutomaton r
        cdef DeterministicAutomaton r2

        r = CAutomaton(None)
        A2 = []
        d1 = imagDict(d, self.A, A2=A2)
        if verb:
            print("d1=%s, A2=%s" % (d1, A2))
        dC = getDict(d, self.A, d1=d1)
        sig_on()
        a = Proj(self.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        r.a[0] = a
        r.A = A2
        if det:
            r2 = r.determinize()
            if simplify:
                return r2.prune().minimize()
            else:
                return r2
        else:
            return r

    def project_on_ith(self, int i, bool det=True, bool simplify=True, bool verb=False):
        """
        Assuming that the alphabet of the automaton are iterable, project on
        the ith coordinate.
        Give a new automaton where labels are replaced by the projection on
        the ith coordinate.

        INPUT:

        - ``i`` -- int - coordinate of projection

        - ``det``  --  bool (default: ``true``) - determinize or not the result

        - ``simplify`` -- bool (default: ``True``) - if True and if det=True, prune
          and minimize

        - ``verb`` -- boolean (default: ``False``) - to activate or
          desactivate the verbose mode

        OUTPUT:

        Return a :class:`CAutomaton` (``det``=``False``)
        or :class:`DeterministicAutomaton` (``det``=``True``)


        EXAMPLES::

            sage: a = dag.Word([('a', 'b'), ('b', 'a'), ('c', 'a')])
            sage: a.project_on_ith(1)
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'abc')], i=0)
            sage: b = a.project_on_ith(0)
            sage: b
            DeterministicAutomaton with 2 states and an alphabet of 1 letter
            sage: b == DeterministicAutomaton([(1, 0, 'a')], i=1)
            True

        TESTS::

            sage: a = dag.Word([('a', 'b'), ('b', 'a'), ('c', 'a')])
            sage: a.project_on_ith(1).has_same_language_as(dag.Word(['b', 'a', 'a']))
            True

            sage: a = dag.Word([0,1])
            sage: a.project_on_ith(0)
            Traceback (most recent call last):
            ...
            TypeError: object of type 'sage.rings.integer.Integer' has no len()

        """
        cdef dict d

        if i < 0:
            raise ValueError("index i=%d cannot be negative" % i)
        d = {}
        for l in self.A:
            if i < len(l):
                d[l] = l[i]
            else:
                raise ValueError("index i=%d must be smaller than the dimension of the label %s" % (i, l))
        return self.project(d, det=det, simplify=simplify, verb=verb)

    def duplicate_transitions(self, d, verb=False):
        """
        Replace every transition of self labeled by a letter l
        by a list of transitions labeled by elements of d[l].
        The result is assumed deterministic !!!

        INPUT:

        - ``d``  -- dictionary giving a list a new letters for each letter of the alphabet of self

        - ``verb`` -- boolean (default: ``False``) fix to ``True``
          for activation the verbose mode

        OUTPUT:

        Return a new :class:`DeterministicAutomaton` with new letters

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'ac', 'b': 'cc', 'ca':'c', 'd':'b'}
            sage: b = a.duplicate_transitions(d)
            sage: b.alphabet
            ['a', 'c']

            sage: a = dag.Word(['a', 'b', 'c'])
            sage: d = { 'a' : [0,1], 'b': [0], 'c':[]}
            sage: b = a.duplicate_transitions(d)
            sage: b.alphabet
            [0, 1]
            sage: b.has_empty_language()
            True
        """
        cdef Automaton a
        cdef InvertDict dC
        cdef DeterministicAutomaton r

        r = DeterministicAutomaton(None)
        A2 = []
        d1 = imagDict2(d, self.A, A2=A2)
        if verb:
            print("d1=%s, A2=%s" % (d1, A2))
        dC = getDict2(d, self.A, d1=d1)
        if verb:
            sig_on()
            printInvertDict(dC)
            sig_off()
        sig_on()
        a = Duplicate(self.a[0], dC, len(A2), verb)
        sig_off()
        if verb:
            print("end...")
        sig_on()
        FreeInvertDict(dC)
        sig_off()
        r.a[0] = a
        r.A = A2
        return r

    def relabel_in_place(self, dict d):
        """
        Change letters of the :class:`DeterministicAutomaton` in-place,
        with respect to the dictionnary ``d``.
        The dictionary is assumed to be one-to-one.

        INPUT:

         - ``d``  -- dictionary that gives new letters from the old ones

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'a', 'b': 'c', 'c':'b', 'd':'b'}
            sage: a.relabel_in_place(d)
            sage: a.alphabet
            ['a', 'c']

            sage: a = dag.Word(['a', 'b'])
            sage: d = {'a':0, 'b':1}
            sage: a.relabel_in_place(d)
            sage: a.alphabet
            [0, 1]
            sage: a == dag.Word([0, 1])
            True

        TESTS::

            sage: a = dag.AnyWord([0,1])
            sage: a.relabel_in_place({})
            Traceback (most recent call last):
            ...
            KeyError: 0

        """
        self.A = [d[c] for c in self.A]

    def permut(self, list A, bool verb=False):
        """
        Permutes (and eventually remove) letters of the alphabet,
        and return permuted new :class:`DeterministicAutomaton` with
        the same language restricted to the new alphabet.

        INPUT:

        - ``A``  -- list of letters in the new order
          (with potentially less letters)

        - ``verb`` -- boolean (default: ``False``) set to ``True`` to
          activate the verbose mode

        OUTPUT:

        Return new :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = dag.Word([0,1,2])
            sage: a.alphabet
            [0, 1, 2]
            sage: b = a.permut([2,1,0])
            sage: b.alphabet
            [2, 1, 0]
            sage: a.has_same_language_as(b)
            True

            sage: a = dag.AnyWord(['a', 'b', 'c'])
            sage: a.alphabet
            ['a', 'c', 'b']
            sage: b = a.permut(['a', 'b'])
            sage: b.alphabet
            ['a', 'b']
            sage: b.has_same_language_as(dag.AnyWord(['a', 'b']))
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: l = [ 'b', 'c', 'a']
            sage: b = a.permut(l)
            sage: b.alphabet
            ['b', 'a']

        TESTS::

            sage: a = dag.Word(['a', 'b'])
            sage: b = a.permut([0,1])
            sage: b.alphabet
            []

        """
        cdef Automaton a
        cdef int *l
        cdef int i, na
        cdef DeterministicAutomaton r

        if verb:
            print("A=%s" % A)
        r = DeterministicAutomaton(None)
        l = <int*>malloc(sizeof(int) * len(A))
        if l is NULL:
            raise MemoryError("Failed to allocate memory for l in "
                              "permut")
        for i in range(self.a.na):
            l[i] = -1
        d = {}
        for i, c in enumerate(self.A):
            d[c] = i
        r.A = []
        for i, c in enumerate(A):
            if d.has_key(c):
                l[i] = d[c]  # l gives the old index from the new one
                r.A.append(c)
        if verb:
            str = "l=["
            for i in range(len(A)):
                str += " %s" % l[i]
            str += " ]"
            print(str)
        sig_on()
        a = Permut(self.a[0], l, len(A), verb)
        free(l)
        sig_off()
        r.a[0] = a
        r.a.na = len(A)
        return r

    def permut_in_place(self, list A, bool verb=False):
        """
        Permutes (and eventually remove) letters of the alphabet in-place,
        without changing the language restricted to the new alphabet.

        INPUT:

        - ``A``  -- list of letters in the new order (number can be less to
          the alphabet)
        - ``verb`` -- boolean (default: ``False``) fix to ``True`` for
          activation the verbose mode

        EXAMPLES::

            sage: a = dag.Word([0,1,2])
            sage: a.alphabet
            [0, 1, 2]
            sage: a.permut_in_place([2,1,0])
            sage: a.alphabet
            [2, 1, 0]
            sage: a.has_same_language_as(dag.Word([0,1,2]))
            True

            sage: a = dag.AnyWord(['a', 'b', 'c'])
            sage: a.alphabet
            ['a', 'c', 'b']
            sage: a.permut_in_place(['a', 'b'])
            sage: a.alphabet
            ['a', 'b']
            sage: a.has_same_language_as(dag.AnyWord(['a', 'b']))
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: l = [ 'b', 'c', 'a']
            sage: a.permut_in_place(l)
            sage: a.alphabet
            ['b', 'a']

        TESTS::

            sage: a = dag.Word(['a', 'b'])
            sage: a.permut_in_place([0,1])
            sage: a.alphabet
            []

        """
        cdef int *l
        cdef int i
        cdef int nA = len(A)
        cdef list A2

        if verb:
            print("A=%s" % A)
        sig_on()
        l = <int*>malloc(sizeof(int) * nA)
        sig_off()
        if l is NULL:
            raise MemoryError("Failed to allocate memory for l in "
                              "permut_in_place")
        for i in range(self.a.na):
            l[i] = -1
        d = {}
        for i, c in enumerate(self.A):
            d[c] = i
        A2 = []
        for i, c in enumerate(A):
            if d.has_key(c):
                l[i] = d[c]  # l gives the old index from the new one
                A2.append(c)
                sig_check()
        if verb:
            str = "l=["
            for i in range(len(A)):
                str += " %s" % l[i]
            str += " ]"
            print(str)
        sig_on()
        PermutOP(self.a[0], l, len(A), verb)
        free(l)
        sig_off()
        self.A = A2
        self.a.na = len(A2)

    def mirror_deterministic(self):
        """
        Return a :class:`DeterministicAutomaton`, whose language is the mirror of the
        language of self.
        We assume the result to be deterministic!
        If it is not the case, you should rather use mirror()
        to get a correct result.

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = dag.Word(['a', 'b', 'c', 'b'])
            sage: b = a.mirror_deterministic()
            sage: b
            DeterministicAutomaton with 5 states and an alphabet of 3 letters
            sage: b.has_same_language_as(dag.Word(['b', 'c', 'b', 'a']))
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.mirror_deterministic()
            DeterministicAutomaton with 4 states and an alphabet of 2 letters

        TESTS::

            sage: a = DeterministicAutomaton([('a','a',0), ('a','b',1), ('b','a',0)], i='a')
            sage: b = a.mirror_deterministic()
            sage: b
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: b.has_same_language_as(a.mirror().determinize())
            False
            sage: b.states
            ['a', 'b']
        """
        cdef DeterministicAutomaton r

        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = MirrorDet(self.a[0])
        sig_off()
        r.A = self.A
        r.S = self.S
        return r

    def mirror(self):
        """
        Return a :class:`CAutomaton`, whose language is the mirror of the
        language of self.

        OUTPUT:

        Return a :class:`CAutomaton`

        EXAMPLES::

            sage: a = dag.Word(['a', 'b', 'c', 'b'])
            sage: b = a.mirror()
            sage: b
            CAutomaton with 5 states and an alphabet of 3 letters
            sage: b.determinize().has_same_language_as(dag.Word(['b', 'c', 'b', 'a']))
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.mirror()
            CAutomaton with 4 states and an alphabet of 2 letters

        TESTS::
        
            sage: a = DeterministicAutomaton([('a','a',0), ('a','b',1), ('b','a',0)], i='a')
            sage: a.mirror()
            CAutomaton with 2 states and an alphabet of 2 letters

        """
        cdef CAutomaton r

        r = CAutomaton(None)
        sig_on()
        r.a[0] = Mirror(self.a[0])
        sig_off()
        r.A = self.A
        r.S = self.S
        return r

    def strongly_connected_components(self, bool no_trivials=False):
        r"""
        Determine a partition into strongly connected components.
        A strongly connected component is a minimal subset of the set
        of states such that
        there is no path going outside of the subset, from a state of
        the subset to a state of the subset.
        This computation is done by computing an efficient topological
        ordering of the graph.

        INPUT:

        - ``no_trivials`` -- bool (default: ``False``) If True, do not take into
          account components without any transition from itself to itself
          (such component contains only one element).

        OUTPUT:

        Return the list of strongly connected components, given as a list of
        list of states.

        EXAMPLES::

            sage: a = dag.Word(['a', 'b']).concatenate(dag.AnyWord(['a', 'b']))
            sage: a
            DeterministicAutomaton with 3 states and an alphabet of 2 letters
            sage: a.strongly_connected_components()
            [[0], [2], [1]]
            
            sage: a = dag.Random(20, [0, 1, 2])
            sage: a.strongly_connected_components()     # random
            [[4],
             [12, 15],
             [2, 6, 13, 18],
             [3, 5, 7, 8, 10, 11, 16, 19],
             [0, 14],
             [17],
             [1],
             [9]]

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (1,0,'b'), (2, 3, 'b')], i=0)
            sage: a.strongly_connected_components()
            [[0, 1], [3], [2]]

            sage: a = DeterministicAutomaton([('a', 'a', 0), ('a', 'b', 1), ('b', 'a', 0)], i='a')
            sage: a.strongly_connected_components()
            [[0, 1]]

        """
        cdef int* l
        cdef int ncc
        cdef dict l2
        cdef int i

        sig_on()
        l = <int*>malloc(sizeof(int) * self.a.n)
        sig_off()
        if l is NULL:
            raise MemoryError("Failed to allocate memory for l in "
                              "strongly_connected_component.")
        sig_on()
        ncc = StronglyConnectedComponents(self.a[0], l)
        sig_off()
        # inverse la liste
        l2 = {}
        for i in range(self.a.n):
            if not l2.has_key(l[i]):
                l2[l[i]] = []
            l2[l[i]].append(i)
        if no_trivials:
            for i in l2.keys():
                if len(l2[i]) == 1:
                    trivial = True
                    for j in range(len(self.A)):
                        if self.a.e[l2[i][0]].f[j] == l2[i][0]:
                            trivial = False
                            break
                    if trivial:
                        # on retire cette composante qui est triviale
                        l2.pop(i)
        sig_on()
        free(l)
        sig_off()
        return l2.values()

    def sub_automaton(self, list l, bool keep_states_labels=True, bool verb=False):
        """
        Compute the sub automaton whose states are given by the set ``l``.

        INPUT:

        - ``l``  -- list of states to keep

        - ``verb`` -- boolean (default: ``False``) set to ``True`` to activate
          the verbose mode

        - ``keep_states_labels``  -- boolean (default: ``True``) keep the labels of states

        OUTPUT:

        Return a :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = dag.Word(['a', 'b', 'c'])
            sage: a.sub_automaton([0,1])
            DeterministicAutomaton with 2 states and an alphabet of 3 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = a.sub_automaton([0,1])
            sage: b
            DeterministicAutomaton with 2 states and an alphabet of 2 letters
            sage: b.has_same_language_as(a)
            True

        TESTS::

            sage: a = dag.Word(['a','b','a'])
            sage: a.sub_automaton([3,4])
            Traceback (most recent call last):
            ...
            ValueError: Element 4 of the list given is not a state (i.e. between 0 and 4).

        """
        cdef DeterministicAutomaton r
        cdef int i, n
        # check that the list l is correct
        if len(l) != len(set(l)):
            raise ValueError("A state of the list appears several times")
        n = self.a.n
        for i in l:
            if i < 0 or i >= n:
                raise ValueError("Element %s of the list given is not a state (i.e. between 0 and %s)."% (i, n))
        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = SubAutomaton(self.a[0], list_to_Dict(l), verb)
        sig_off()
        r.A = self.A
        if keep_states_labels and self.S is not None:
            r.S = []
            for i in l:
                r.S.append(self.S[i])
        return r

    def simplify(self):
        """
        Prune and minimize self.
        This gives a canonical automaton for the language of self.

        OUTPUT:

        :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = dag.Word(['a','b','a','b'])
            sage: a.shift_list_in_place(['a','b'])
            sage: a
            DeterministicAutomaton with 5 states and an alphabet of 2 letters
            sage: a.simplify()
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

        """
        return self.prune().minimize()

    def minimize(self, verb=False):
        """
        Compute the minimal automaton
        by Hopcroft's algorithm
        see [Hopcroft]

        INPUT:

         - ``verb`` -- boolean (default: ``False``) set to ``True`` to
           activate the verbose mode


        OUTPUT:

        Return a minimized :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([('a','b',0),('a','c',1),('b','a',0),
            ....:                   ('b','d',1),('c','e',0),('c','f',1),
            ....:                   ('d','e',0),('d','f',1),('e','e',0),
            ....:                   ('e','f',1),('f','f',0),('f','f',1)],
            ....:                   i='a', final_states=['c','d','e'])
            sage: a.minimize()
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0,1,0), (0,5,1), (1,2,1), (1,6,0),
            ....:                   (2,0,0), (2,2,1), (3,2,0), (3,6,1),
            ....:                   (4,7,0), (4,5,1), (5,2,0), (5,6,1),
            ....:                   (6,6,0), (6,4,1), (7,6,0), (7,2,1)],
            ....:                   i=0, final_states=[2,6])
            sage: a.minimize()
            DeterministicAutomaton with 5 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.minimize()
            DeterministicAutomaton with 3 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                    (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                    (10,30,'z')], i=10)
            sage: a.minimize()
            DeterministicAutomaton with 1 state and an alphabet of 3 letters

        TESTS::

            sage: a = DeterministicAutomaton([('a','b',0),('a','c',1),('b','a',0),
            ....:                   ('b','d',1),('c','e',0),('c','f',1),
            ....:                   ('d','e',0),('d','f',1),('e','e',0),
            ....:                   ('e','f',1),('f','f',0),('f','f',1)],
            ....:                   i='a', final_states=['c','d','e'])
            sage: b = DeterministicAutomaton([(0, 0, 0), (0, 2, 1), (1, 1, 0),
            ....:                   (1, 0, 1), (2, 2, 0), (2, 2, 1)],
            ....:                  A=[0, 1], i=1, final_states=[0])
            sage: a.minimize() == b
            True
            
            sage: a = DeterministicAutomaton([(0,1,0), (0,5,1), (1,2,1), (1,6,0),
            ....:                   (2,0,0), (2,2,1), (3,2,0), (3,6,1),
            ....:                   (4,7,0), (4,5,1), (5,2,0), (5,6,1),
            ....:                   (6,6,0), (6,4,1), (7,6,0), (7,2,1)],
            ....:                   i=0, final_states=[2,6])
            sage: b = DeterministicAutomaton([(0, 0, 0), (0, 2, 1), (1, 3, 0),
            ....:                   (1, 0, 1), (2, 4, 0), (2, 1, 1),
            ....:                   (3, 2, 0), (3, 3, 1), (4, 0, 0),
            ....:                   (4, 3, 1)], i=2, final_states=[0, 3])
            sage: a.minimize() == b
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.minimize(True)
            transition i[0][0] = [ ]
            transition i[0][1] = [ ]
            transition i[1][0] = [ 0 ]
            transition i[1][1] = [ ]
            transition i[2][0] = [ ]
            transition i[2][1] = [ ]
            transition i[3][0] = [ ]
            transition i[3][1] = [ 2 ]
            transition i[4][0] = [ 1 2 3 4 ]
            transition i[4][1] = [ 0 1 3 4 ]
            partition = [ 0 1 2 3 4 ]
            partitioni = [ 0 1 2 3 4 ]
            Initial partition :
            class 0 : 0 1 2 3
            class 1 : 4
            split 1 0...
            new visited class : 0 (1 parent of 4)
            re-visited class : 0 (2 parent of 4)
            re-visited class : 0 (3 parent of 4)
            new visited class : 1 (4 parent of 4)
            class 0 : 1 2 3 0
            class 1 : 4
            2 class encountered
            class 0 : l = 0 3 4 = h
            class 1 : l = 4 5 5 = h
            split 1 1...
            new visited class : 2 (0 parent of 4)
            new visited class : 0 (1 parent of 4)
            re-visited class : 0 (3 parent of 4)
            new visited class : 1 (4 parent of 4)
            class 0 : 1 3 2
            class 1 : 4
            class 2 : 0
            3 class encountered
            class 2 : l = 3 4 4 = h
            class 0 : l = 0 2 3 = h
            class 1 : l = 4 5 5 = h
            split 3 0...
            class 0 : 1 3
            class 1 : 4
            class 2 : 0
            class 3 : 2
            0 class encountered
            split 3 1...
            class 0 : 1 3
            class 1 : 4
            class 2 : 0
            class 3 : 2
            0 class encountered
            split 2 0...
            class 0 : 1 3
            class 1 : 4
            class 2 : 0
            class 3 : 2
            0 class encountered
            split 2 1...
            class 0 : 1 3
            class 1 : 4
            class 2 : 0
            class 3 : 2
            0 class encountered
            Final partition :
            class 0 : 1 3
            class 1 : 4
            class 2 : 0
            class 3 : 2
            a.i = 0 class 2
            removes the sink state  1...
            DeterministicAutomaton with 3 states and an alphabet of 2 letters
        """
        cdef DeterministicAutomaton r

        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = Minimise(self.a[0], verb)
        sig_off()
        r.A = self.A
        return r

    def adjacency_matrix(self, sparse=None):
        """
        Compute the adjacency matrix of the :class:`DeterministicAutomaton`

        INPUT:

        - ``sparse`` -- indicate if the return matrix is sparse or not

        OUTPUT:

        Return the corresponding adjacency matrix

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.adjacency_matrix()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 0 1]
            [0 0 0 0]

            sage: a = dag.Word(['a','b','a','c'])
            sage: a.adjacency_matrix()
            [0 1 0 0 0]
            [0 0 1 0 0]
            [0 0 0 1 0]
            [0 0 0 0 1]
            [0 0 0 0 0]

            sage: a = dag.AnyLetter(['a','b','c'])
            sage: a.adjacency_matrix()
            [0 3]
            [0 0]

            sage: a = dag.Random(500, ['a', 'b', 'c'])
            sage: a.adjacency_matrix()
            500 x 500 sparse matrix over Integer Ring (use the '.str()' method to see the entries)

        """
        cdef int i, j, f
        cdef dict d

        if sparse is None:
            if self.a.n <= 128:
                sparse = False
            else:
                sparse = True

        d = {}
        for i in range(self.a.n):
            for j in range(self.a.na):
                f = self.a.e[i].f[j]
                if f != -1:
                    if d.has_key((i, f)):
                        d[(i, f)] += 1
                    else:
                        d[(i, f)] = 1
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import IntegerRing
        return matrix(IntegerRing(), self.a.n, self.a.n, d, sparse=sparse)

    def delete_state(self, int i):
        """
        Gives a copy of the :class:`DeterministicAutomaton` but
        without the state ``i``.

        INPUT:

        - ``i``  - int - the state to remove

        OUTPUT:

        :class:`DeterministicAutomaton`

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.delete_state(2)
            DeterministicAutomaton with 3 states and an alphabet of 2 letters
            sage: a.delete_state(2).delete_state(1)
            DeterministicAutomaton with 2 states and an alphabet of 2 letters

        TESTS::

            sage: a = DeterministicAutomaton([('a', 'a', 0), ('a', 'b', 1), ('b', 'a', 0), ('b', 'c', 1), ('c', 'a', 0)], i='a')
            sage: a = a.delete_state(a.states.index('b'))
            sage: b = DeterministicAutomaton([('a', 'a', 0), ('c', 'a', 0)], A=[0, 1], i='a')
            sage: a == b
            True

            sage: a = dag.AnyWord(['a','b'])
            sage: a.delete_state(1)
            Traceback (most recent call last):
            ...
            ValueError: 1 is not a state (should be between 0 and 0)

        """
        cdef DeterministicAutomaton r
        if i < 0 or i >= self.a.n:
            raise ValueError("%s is not a state (should be between 0 and %s)"%(i, self.a.n-1))
        r = DeterministicAutomaton(None)
        sig_on()
        r.a[0] = DeleteVertex(self.a[0], i)
        sig_off()
        r.A = self.A
        if self.S is not None:
            r.S = [s for j,s in enumerate(self.S) if j != i]
        return r

    def delete_state_in_place(self, int i):
        """
        Delete vertex ``i`` in-place.

        INPUT:

        - ``i``  -- int number of vertex to remove

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.delete_state_in_place(2)
            sage: a
            DeterministicAutomaton with 3 states and an alphabet of 2 letters
            sage: a.delete_state_in_place(1)
            sage: a
            DeterministicAutomaton with 2 states and an alphabet of 2 letters

        TESTS::

            sage: a = dag.AnyWord(['a','b'])
            sage: a.delete_state_in_place(1)
            Traceback (most recent call last):
            ...
            ValueError: 1 is not a state (should be between 0 and 0)

        """
        if i < 0 or i >= self.a.n:
            raise ValueError("%s is not a state (should be between 0 and %s)"%(i, self.a.n-1))
        sig_on()
        DeleteVertexOP(self.a, i)
        sig_off()
        if self.S is not None:
            self.S = [s for j,s in enumerate(self.S) if j != i]

    def spectral_radius(self, bool approx=True, bool couple=False, bool only_non_trivial=False, bool verb=False):
        """
        Return the spectral radius of the underlying graph.

        INPUT:

        - ``approx`` - boolean (default: ``True``) If True gives an approximation,
          otherwise gives the exact value as an algebraic number. 

        - ``only_non_trivial`` - boolean (default: ``False``) - if True,
          don't take into account strongly connected components of
          cardinality one.

        - ``couple``- boolean (default: ``False``) - return a interval
        containing the spectral radius (only when ``approx`` is ``True``)

        - ``verb`` - boolean (default: ``False``) - set to ``True`` to activate
          the verbose mode

        OUTPUT:

        A positive real algebraic number if approx is False,
        and an float number otherwise.

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.spectral_radius()
            0.0000000000...
            
            sage: a = DeterministicAutomaton([(10,10,'a'),(10,20,'b'),(20,10,'a')], i=10)
            sage: a.spectral_radius()
            1.61803398...
            
            sage: a = DeterministicAutomaton([(10,10,'a'),(10,20,'b'),(20,10,'a')], i=10)
            sage: a.spectral_radius(couple=True)
            (1.61803398..., 1.61803398...)
            
            sage: a = DeterministicAutomaton([(10,10,'a'),(10,20,'b'),(20,10,'a')], i=10)
            sage: a.spectral_radius(approx=False).minpoly()
            x^2 - x - 1

        """
        a = self.minimize()
        if verb:
            print("minimal Automata : %s" % a)
        l = a.strongly_connected_components()
        if verb:
            print("%s component strongly connex." % len(l))
        r = 0  # valeur propre maximale trouvée
        spm = (0,0) #encadrement de l'approximation maximale trouvée
        for c in l:
            if not only_non_trivial or len(c) > 1:
                if verb:
                    print("component with %s states..." % len(c))
                b = a.sub_automaton(c)
                if approx:
                    g = b.get_digraph()
                    if verb:
                        print("g=%s" % g)
                    if g.is_aperiodic():
                        if verb:
                            print("aperiodic")
                        sp = g.spectral_radius()
                        if verb:
                            print("sp=%s" % sp)
                    else:
                        if verb:
                            print("non aperiodic")
                        r = max(b.adjacency_matrix().charpoly().real_roots())
                        sp = (r,r)
                    spm = max([sp, spm], key=lambda x:x[0])
                    if verb:
                        print("spm=%s" % spm)
                else:
                    m = b.adjacency_matrix()
                    cp = m.charpoly()
                    fs = cp.factor()
                    if verb:
                        print(fs)
                    for f in fs:
                        if verb:
                            print(f)
                        from sage.functions.other import real_part
                        from sage.rings.qqbar import AlgebraicRealField
                        r = max([ro[0] for ro in f[0].roots(ring=AlgebraicRealField())] + [r])
        if approx:
            if couple:
                return spm
            else:
                return spm[1]
        else:
            return r

    def has_empty_language(self):
        r"""
        Test if the  :class:`DeterministicAutomaton` has an empty language.

        OUTPUT:

        Return ``True`` the :class:`DeterministicAutomaton` has a empty language
        ``False`` if not

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.has_empty_language()
            False

        """
        sig_on()
        res = emptyLanguage(self.a[0])
        answ = c_bool(res)
        sig_off()
        return answ

    def has_same_language_as(self, DeterministicAutomaton a, bool minimized=False,
                        bool pruned=False, bool verb=False):
        """
        Test if the languages of :class:`DeterministicAutomaton` ``self`` and ``a`` are
        equal or not.
        If the alphabets are differents it returns False, even if the automata describe the same languages.

        INPUT:

        - ``a``  -- the :class:`DeterministicAutomaton` to compare

        - ``minimized``  -- (default: ``False``) if minimization is
          required or not

        - ``pruned``  -- (default: ``False``) if emondation is required or not

        - ``verb`` -- boolean (default: ``False``) set to ``True`` to activate
          the verbose mode

        OUTPUT:

        Return ``True`` if the both :class:`DeterministicAutomaton` have
        the same language and same alphabet ``False`` if not

        EXAMPLES::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=3)
            sage: c = DeterministicAutomaton([(3, 2, 'd'), (1, 2, 'c')], i=2)
            sage: a.has_same_language_as(b)
            False
            sage: a.has_same_language_as(c)
            False
            sage: c = DeterministicAutomaton([(3, 2, 'd'), (1, 2, 'c')])
            sage: a.has_same_language_as(c)
            False

            sage: a = DeterministicAutomaton([(0,0,0), (0,1,1), (1,1,0), (1,1,1)], i=0)
            sage: dag.AnyWord([0,1]).has_same_language_as(a)
            True
        """
        cdef Dict d
        cdef int i, j
        if set(self.A) != set(a.A):
            return False
        sig_on()
        d = NewDict(self.a.na)
        sig_off()
        for i in range(self.a.na):
            for j in range(a.a.na):
                if self.A[i] == a.A[j]:
                    d.e[i] = j
                    if verb:
                        print("%d -> %d" % (i, j))
                    break
                sig_check()
        if verb:
            sig_on()
            printDict(d)
            sig_off()
        sig_on()
        res = equalsLanguages(self.a, a.a, d, minimized, pruned, verb)
        answ = c_bool(res)
        FreeDict(&d)
        sig_off()
        sig_on()
        FreeDict(&d)
        sig_off()
        return answ

    def intersect(self, DeterministicAutomaton a, bool verb=False):
        """
        Determine if the languages of the :class:`DeterministicAutomaton` ``self``
        and ``a`` have a non-empty intersection.

        INPUT:

        -  ``a``  -- the :class:`Detautomaton` to intersect

        - ``verb`` -- boolean (default: ``False``) True to activate
          the verbose mode

        OUTPUT:

        Return ``True`` if the intersection of the languages is non-empty,
        return ``False`` otherwise.

        EXAMPLES::

            sage: a = dag.Word("ababb")
            sage: b = dag.Word("aba")
            sage: a.intersect(b)
            False
            sage: a.prefix_closure().intersect(b)
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.intersect(b)
            True
            sage: b = DeterministicAutomaton([(3, 2, 'a'), (1, 2, 'd')])
            sage: a.intersect(b)
            False
        """
        sig_on()
        res = Intersect(self.a[0], a.a[0], verb)
        answ = c_bool(res)
        sig_off()
        return answ

    def find_word(self, bool verb=False):
        """
        Find a word in the language of the automaton

        INPUT:

        - ``verb`` -- (default: ``False``)  the verbose parameter

        OUTPUT:

        return a word of the language as a list
        of letters if it exists, otherwise return None

        EXAMPLES::

            sage: a = dag.Word("gabian")
            sage: a.find_word()
            ['g', 'a', 'b', 'i', 'a', 'n']

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.find_word()
            []
            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')])
            sage: a.find_word()

            sage: a = DeterministicAutomaton([(0, 0, 'x'), (0, 1, 'y')], i=0, final_states=[1])
            sage: a.find_word()
            ['y']
        """
        cdef Dict w
        cdef list r

        sig_on()
        res = findWord(self.a[0], &w, verb)
        sig_off()
        if not res:
            return None
        r = []
        for i in range(w.n):
            r.append(self.A[w.e[i]])
        sig_on()
        FreeDict(&w)
        sig_off()
        return r

    def shortest_word(self, i=None, f=None, bool verb=False):
        """
        Compute a shortest words of the language of self.

        INPUT:

        - ``i`` -- (default: None)  the initial state

        - ``f`` -- (default: None)  the final state

        - ``verb`` -- (default: False)  the verbose parameter

        OUTPUT:

        return a word, as list of letters, or None if the language is empty.

        EXAMPLES::

            sage: a = dag.Word("gabian")
            sage: a.shortest_word()
            ['g', 'a', 'b', 'i', 'a', 'n']

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')])
            sage: a.shortest_word(i=2, f=3)
            ['b']

            sage: a = DeterministicAutomaton([(0, 0, 'x'), (0, 1, 'y')], i=0, final_states=[1])
            sage: a.shortest_word()
            ['y']

            sage: a = dag.Empty(['a', 'b'])
            sage: a.shortest_word()
            

        TESTS::

            sage: a = dag.AnyWord("ab")
            sage: a.shortest_word(i=-2)
            Traceback (most recent call last):
            ...
            ValueError: -2 is not a state (i.e. between 0 and 0)

            sage: a = dag.AnyWord("ab")
            sage: a.shortest_word(f=-2)
            Traceback (most recent call last):
            ...
            ValueError: -2 is not a state (i.e. between 0 and 0)

        """
        cdef Dict w
        if i is None:
            i = self.a.i
        else:
            if i < -1 or i >= self.a.n:
                raise ValueError("%s is not a state (i.e. between 0 and %s)" % (i, self.a.n-1))
        if f is None:
            f = -1
        else:
            if f < -1 or f >= self.a.n:
                raise ValueError("%s is not a state (i.e. between 0 and %s)" % (f, self.a.n-1))
        sig_on()
        res = shortestWord(self.a[0], &w, i, f, verb)
        sig_off()
        if not res:
            return None
        r = []
        for i in range(w.n):
            r.append(self.A[w.e[i]])
        sig_on()
        FreeDict(&w)
        sig_off()
        return r

    def shortest_words(self, i=None, bool verb=False):
        """
        Compute the list of shortest words from the initial state to every state.

        INPUT:

        - ``i`` -- (default: None) - the initial state

        - ``verb`` -- bool (default: False) - the verbose parameter

        OUTPUT:

        return a list of words, as a list of list of letters.

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.shortest_words()
            [[], ['a'], ['a', 'b'], ['a', 'b', 'a'], ['a', 'b', 'a', 'a']]

            sage: a = dag.AnyWord("ab")
            sage: a.shortest_words()
            [[]]

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.shortest_words()
            [[], ['a'], [], []]

        TESTS::

            sage: a = dag.AnyWord("ab")
            sage: a.shortest_words(i=-2)
            Traceback (most recent call last):
            ...
            ValueError: -2 is not a state (i.e. between 0 and 0)
        """
        cdef Dict* w
        sig_on()
        w = <Dict*>malloc(sizeof(Dict) * self.a.n)
        sig_off()
        if w is NULL:
            raise MemoryError("Failed to allocate memory for w in "
                              "shortest_words")
        if i is None:
            i = self.a.i
        else:
            if i < -1 or i >= self.a.n:
                raise ValueError("%s is not a state (i.e. between 0 and %s)" % (i, self.a.n-1))
        sig_on()
        res = shortestWords(self.a[0], w, i, verb)
        sig_off()
        if not res:
            return None
        rt = []
        for j in range(self.a.n):
            r = []
            for i in range(w[j].n):
                r.append(self.A[w[j].e[i]])
            rt.append(r)
            sig_on()
            FreeDict(&w[j])
            sig_off()
        sig_on()
        free(w)
        sig_off()
        return rt

    def recognized_word(self, w):
        """
        Determine if the word ``w`` is recognized or nor not by the automaton

        INPUT:

        - ``w`` -- a list of letters

        OUTPUT:

        return ``True`` if the word is recognized, otherwise ``False``

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.recognized_word("aba")
            False
            sage: a.recognized_word("abaa")
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.recognized_word(['a', 'b', 'b'])
            False
            sage: a.recognized_word("a")
            True

        TESTS::

            sage: a = dag.Word("ab")
            sage: a.recognized_word(['a', 0])
            Traceback (most recent call last):
            ...
            KeyError: 0

        """
        cdef int i, e
        cdef dict d     
        e = self.a.i
        if e == -1:
            return False
        d = {}
        for i, a in enumerate(self.A):
            d[a] = i
        for a in w:
            e = self.a.e[e].f[d[a]]
            if e == -1:
                return False
        return c_bool(self.a.e[e].final)

    def add_state(self, bool final, label=""):
        """
        Add a state in the automaton

        INPUT:

        - ``final`` -- boolean indicate if the added state is final
        
        - ``label`` (default: "") -- label of the new state

        OUTPUT:

        return the new state (which is an integer)

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: s = a.add_state(true)
            sage: a.add_transition(2,'b',s)
            sage: a.recognized_word("abb")
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.add_state(True)
            4
            sage: a.add_state(False)
            5

            sage: a = DeterministicAutomaton([(10,10,'x'),(10,20,'y'),(20,20,'z'),\
                (20,10,'y'),(20,30,'x'),(30,30,'y'),(30,10,'z'),(30,20,'x'),\
                (10,30,'z')], i=10)
            sage: a
            DeterministicAutomaton with 3 states and an alphabet of 3 letters
            sage: a.add_state(True)
            3
            sage: a
            DeterministicAutomaton with 4 states and an alphabet of 3 letters

        """
        sig_on()
        AddState(self.a, final)
        sig_off()
        if self.S is not None:
            self.S.append(label)
        return self.a.n-1

    def add_transition(self, int i, l, int j):
        """
        Add an edge in the automaton

        INPUT:

        - ``i`` - int - the first state

        - ``l`` -- the label of edge

        - ``j`` - int - the second state

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.add_transition(2,'a',4)
            sage: a.recognized_word("aba")
            True

            sage: a = dag.Word("abaa")
            sage: s = a.add_state(true)
            sage: a.add_transition(2,'b',s)
            sage: a.recognized_word("abb")
            True

        TESTS::

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.add_transition(2,'a',1)
            sage: a == DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b'), (2, 1, 'a')], i=0)
            True
            sage: a.add_transition(2,'v',1)
            Traceback (most recent call last):
            ...
            ValueError: The letter v doesn't exist.
            sage: a.add_transition(2,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state 6 doesn't exist.
            sage: a.add_transition(5,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state 5 doesn't exist.

        """
        if i < 0 or i >= self.a.n:
            raise ValueError("The state %s doesn't exist." % i)
        if j< 0 or j >= self.a.n:
            raise ValueError("The state %s doesn't exist." % j)
        try:
            k = self.A.index(l)
        except Exception:
            raise ValueError("The letter %s doesn't exist." % l)
        sig_on()
        self.a.e[i].f[k] = j
        sig_off()

    @property
    def n_states(self):
        """
        return the numbers of states

        OUTPUT:

        int

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.n_states
            5

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.n_states
            4

        """
        return self.a.n

    @property
    def n_letters(self):
        """
        return the numbers of letters of the alphabet.

        OUTPUT:

        int

        EXAMPLES::

            sage: a = dag.Word("abaaca")
            sage: a.n_letters
            3

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.n_letters
            2

        """
        return self.a.na

    def bigger_alphabet(self, A):
        """
        Gives a copy of the :class:`DeterministicAutomaton`, but with the
        bigger alphabet ``A``.

        INPUT:

        - ``A`` --  alphabet of the new automaton. We assume that it
          contains the current alphabet of the automaton.

        OUTPUT:

        return a :class:`DeterministicAutomaton` with a bigger alphabet ``nA``

        EXAMPLES::

            sage: a = dag.AnyLetter(['a', 'b'])
            sage: a.bigger_alphabet(['a','b','c'])
            DeterministicAutomaton with 2 states and an alphabet of 3 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')])
            sage: a.bigger_alphabet(['a','b','c'])
            DeterministicAutomaton with 4 states and an alphabet of 3 letters

        TESTS::

            sage: a = dag.AnyWord(['a', 'b'])
            sage: a.bigger_alphabet(['a','1','c'])
            Traceback (most recent call last):
            ...
            ValueError: Letter 'b' not found in the new alphabet

        """
        cdef Dict d
        cdef int i
        cdef DeterministicAutomaton r
        r = DeterministicAutomaton(None)
        sig_on()
        d = NewDict(self.a.na)
        sig_off()
        for i in range(self.a.na):
            try:
                d.e[i] = A.index(self.A[i])
            except Exception:
                raise ValueError("Letter %r not found in the new alphabet"%self.A[i])
            sig_check()
        sig_on()
        r.a[0] = BiggerAlphabet(self.a[0], d, len(A))
        sig_off()
        r.A = A
        return r

    def complementary_in_place(self):
        """
        Change the language of the automaton to the complementary in-place.
        If the language of self was L and the alphabet is A,
        then the new language of self is A^*\L
        (i.e. words over the alphabet A that are not in L)

        OUTPUT:

        return None
        (the operation is in-place)

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.complementary_in_place()
            sage: a
            DeterministicAutomaton with 6 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complementary_in_place()
            sage: a
            DeterministicAutomaton with 5 states and an alphabet of 2 letters
        """
        cdef i
        self.complete_in_place()
        for i in range(self.a.n):
            self.a.e[i].final = not self.a.e[i].final

    def complementary(self):
        """
        Gives an automaton whose language is the complementary.

        OUTPUT:

        :class:`DeterministicAutomaton`
        return  a new automaton whose language is the complementary
        of the language of ``self``

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.complementary()
            DeterministicAutomaton with 6 states and an alphabet of 2 letters

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complementary()
            DeterministicAutomaton with 5 states and an alphabet of 2 letters
        """
        a = self.copy()
        a.complementary_in_place()
        return a

    def included(self, DeterministicAutomaton a, bool pruned=False, bool verb=False):
        r"""
        Test if the language of self is included in the language of ``a``

        INPUT:

        - ``a`` --  a :class:`DeterministicAutomaton`

        - ``pruned`` -- (default: False) - set to True if the automaton self is
          already pruned, in order to avoid unuseful computation.

        - ``verb`` -- (default: False) - verbose parameter

        OUTPUT:

        Boolean
        return ``True`` if the language of ``self`` is included in the
        language of ``a``, otherwise return ``False`` 

        EXAMPLES::

            sage: a = dag.AnyWord("ab")
            sage: b = dag.Word("aba")
            sage: b.included(a)
            True
            sage: a.included(b)
            False

            sage: a = dag.AnyWord("abc")
            sage: b = dag.AnyWord("ab")
            sage: a.included(b)
            False
            sage: b.included(a)
            True

            sage: a = DeterministicAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.included(a)
            True
            sage: b = DeterministicAutomaton([(0, 1, 'c')], i=0)
            sage: a.included(b)
            False
            sage: b = DeterministicAutomaton([(0, 1, 'a')], i=0)
            sage: b.included(a)
            True
            sage: b = DeterministicAutomaton([(0, 1, 'c')])
            sage: b.included(a)
            True
        """
        cdef DeterministicAutomaton b, a2
        cdef list A
        if self.A != a.A:
            A = list(set(a.A+self.A))
            b = self.bigger_alphabet(A)
            a2 = a.bigger_alphabet(A)
        else:
            b = self
            a2 = a
        sig_on()
        res = Included(b.a[0], a2.a[0], pruned, verb)
        answ = c_bool(res)
        sig_off()
        return answ

#        d = {}
#        for l in self.A:
#            if l in a.A:
#                d[(l,l)] = l
#        if verb:
#            print("d=%s"%d)
#        a.complete_in_place()
#        cdef DeterministicAutomaton p = self.product(a, d, verb=verb)
#
#        #set final states
#        cdef int i,j
#        cdef n1 = self.a.n
#        for i in range(n1):
#            for j in range(a.a.n):
#                p.a.e[i+n1*j].final = self.a.e[i].final and not a.a.e[j].final
#
#        if step == 1:
#            return p;
#
#        return p.has_empty_language()


    def random_word(self, int nmin=-1, int nmax=100):
        r"""
        Return a random word recognized by the automaton, by following a
        random path in the automaton from the initial state. If we don't fall
        into the set of final states before reaching the maximal length
        ``nmax``, then return ``word not found!``.

        INPUT:

        - ``nmin`` -- (default: -1) - minimal length of the word

        - ``nmax`` -- (default: 100) - maximal length of the word

        OUTPUT:

        Return a random word of length between ``nmin`` and ``nmax`` if found.
        Otherwise return ``word not found!``.

        EXAMPLES::

            sage: a = dag.Word("abaa")
            sage: a.random_word()
            ['a', 'b', 'a', 'a']

            sage: a = dag.AnyWord("ab")
            sage: a.random_word()
            []

            sage: a = dag.Word("aba").union(dag.Word("bab"))
            sage: a.random_word()   # random
            ['a', 'b', 'a']

            sage: a = dag.Empty(['a', 'b'])
            sage: a.random_word()
            'word not found!'

        """
        cdef int i, j, l, na
        cdef list w, li
        i = self.a.i
        if i == -1:
            return "word not found!"
        w = []
        na = len(self.A)
        if nmin < 0:
            nmin = 0
        from sage.misc.prandom import random
        for j in range(nmin):
            li = [l for l in range(na) if self.successor(i, l) != -1]
            l = li[(int)(random() * len(li))]
            w.append(self.A[l])
            i = self.successor(i, l)
        # continue the word to get into a final state
        for j in range(nmax-nmin):
            if self.a.e[i].final:
                break
            li = [l for l in range(na) if self.successor(i, l) != -1]
            if li == []:
                return "word not found!"
            l = li[(int)(random() * len(li))]
            w.append(self.A[l])
            i = self.successor(i, l)
        if i < 0 or not self.a.e[i].final:
            return "word not found!"
        return w
