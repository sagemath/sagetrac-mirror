# coding=utf8
"""
Fast automaton of Finite state machines using C
FastAutomaton for determinist automata and NFastAutomaton for non determinist


AUTHORS:

- Paul Mercat (2013)- I2M AMU Aix-Marseille Universite - initial version
- Dominique Benielli (2018) Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in -SageMath

REFERENCES:

.. [Hopcroft] "Around Hopcroft’s Algorithm"  Manuel of BACLET and
    Claire PAGETTI.

"""

#*****************************************************************************
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
#*****************************************************************************
from __future__ import print_function
from libc.stdlib cimport malloc, free

cimport sage.combinat.words.cautomata

from cysignals.signals cimport sig_on, sig_off, sig_check

#ctypedef Automate Automaton

cdef extern from "automataC.h":

    cdef cppclass Dict:
        int* e
        int n
    cdef cppclass InvertDict:
        Dict* d
        int n

    bool DotExists ()
    #    Automaton NewAutomaton (int n, int na)
    #    void FreeAutomaton (Automaton *a)
    int hashAutomaton(Automaton a)
    void FreeNAutomaton(NAutomaton *a)
    Automaton CopyAutomaton(Automaton a, int nalloc, int naalloc)
    Automaton PieceAutomaton(Automaton a, int *w, int n, int e)
    void init(Automaton *a)
    void printAutomaton(Automaton a)
    void plotDot(const char *file, Automaton a, const char **labels, const char *graph_name, double sx, double sy, const char **vlabels, bool html, bool verb, bool run_dot)
    void NplotDot (const char *file, NAutomaton a, const char **labels, const char *graph_name, double sx, double sy, bool run_dot)
    Automaton Product(Automaton a1, Automaton a2, Dict d, bool verb)
    Automaton Determinize(Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
    Automaton DeterminizeN(NAutomaton a, bool puits, int verb)
    NAutomaton Concat(Automaton a, Automaton b, bool verb)
    NAutomaton CopyN(Automaton a, bool verb)
    void AddEdgeN(NAutomaton *a, int e, int f, int l)
    void AddPathN(NAutomaton *a, int e, int f, int *l, int len, bool verb)
    NAutomaton Proj(Automaton a, Dict d, bool verb)
    void ZeroComplete(Automaton *a, int l0, bool verb)
    Automaton ZeroComplete2(Automaton *a, int l0, bool etat_puits, bool verb)
    Automaton ZeroInv(Automaton *a, int l0)
    Automaton prune_inf(Automaton a, bool verb)
    Automaton prune(Automaton a, bool verb)
    Automaton pruneI(Automaton a, bool verb)
    void AccCoAcc(Automaton *a, int *coa)
    void CoAcc(Automaton *a, int *coa)
    bool equalsAutomaton(Automaton a1, Automaton a2)
    Dict NewDict(int n)
    void FreeDict(Dict *d)
    void printDict(Dict d)
    InvertDict NewInvertDict(int n)
    void FreeInvertDict(InvertDict id)
    void printInvertDict(InvertDict id)
    Automaton Duplicate(Automaton a, InvertDict id, int na2, bool verb)
    Automaton TransposeDet(Automaton a)
    NAutomaton Transpose(Automaton a)
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
    void AddEtat(Automaton *a, bool final)
    bool IsCompleteAutomaton(Automaton a)
    bool CompleteAutomaton(Automaton *a)
    Automaton BiggerAlphabet(Automaton a, Dict d, int nna) #copy the automaton with a new bigger alphabet
    bool findWord(Automaton a, Dict *w, bool verb)
    bool shortestWord(Automaton a, Dict *w, int i, int f, bool verb)
    bool shortestWords(Automaton a, Dict *w, int i, bool verb)
    bool rec_word(Automaton a, Dict d)
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
    r = NewDict(len(A))
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
    r = NewInvertDict(len(A))
    cdef int i
    if d1 is None:
        d1 = imagDict2(d, A)
    # print(d1)
    for i in range(r.n):
        if d.has_key(A[i]):
            r.d[i] = NewDict(len(d[A[i]]))
            for j in range(r.d[i].n):
                r.d[i].e[j] = d1[d[A[i]][j]]
        else:
            r.d[i].n = 0
    return r

cdef imagProductDict(dict d, list A1, list A2, list Av=[]):
    """
    Dictionary which is numbering the projeted alphabet
    """
    dv = {}
    i = 0
    for a1 in A1:
        for a2 in A2:
            if d.has_key((a1,a2)):
                if not dv.has_key(d[(a1,a2)]):
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
    r = NewDict(n1*n2)
    Keys = d.keys()
    if verb:
        print("Keys=%s" % Keys)
    for (a1, a2) in Keys:
        if d1.has_key(a1) and d2.has_key(a2):
            r.e[d1[a1]+d2[a2]*n1] = dv[d[(a1, a2)]]
    return r

cdef Automaton getAutomaton (a, initial=None, F=None, A=None):
    d = {}
    da = {}
    if F is None:
        if not hasattr(a, 'F'):
            F = a.vertices()
        else:
            F = a.F
    cdef Automaton r

    if A is None:
        A = list(a.Alphabet)
    V = list(a.vertices())
    cdef int n = len(V)
    cdef int na = len(A)

    sig_on()
    r = NewAutomaton(n, na)
    init(&r)
    sig_off()
    for i in range(na):
        da[A[i]] = i
    for i in range(n):
        r.e[i].final = 0
        d[V[i]] = i
    for v in F:
        if not d.has_key(v):
            sig_on()
            FreeAutomaton(&r)
            r = NewAutomaton(0,0)
            sig_off()
            raise ValueError("Incorrect set of final states.")
        r.e[d[v]].final = 1

    if initial is None:
        if not hasattr(a, 'I'):
            I = []
        else:
            I = list(a.I)
        if len(I) > 1:
            raise ValueError("The automata must be determist (I=%s)" % a.I)
        if len(I) >= 1:
            r.i = d[I[0]]
        else:
            r.i = -1
    else:
        r.i = d[initial]

    for e, f, l in a.edges():
        r.e[d[e]].f[da[l]] = d[f]
    a.dA = da
    a.S = V
    a.dS = d
    return r

cdef AutomatonGet(Automaton a, A):
    """
    Transform an Automaton a with an alphabet A to a DiGraph
    """
    from sage.graphs.digraph import DiGraph
    r = DiGraph(multiedges=True, loops=True)
    cdef int i, j
    r.F = []
    for i in range(a.n):
        for j in range(a.na):
            if a.e[i].f[j] != -1:
                r.add_edge((i, a.e[i].f[j], A[j]))
        if a.e[i].final:
            r.F.append(i)
    r.I = [a.i]
    return r

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

cdef AutomatonToDiGraph(Automaton a, A):
    from sage.graphs.digraph import DiGraph
    L = []
    cdef int i, j
    for i in range(a.n):
        for j in range(a.na):
            if a.e[i].f[j] != -1:
                L.append((i, a.e[i].f[j], A[j]))
    return DiGraph(L, loops=True, multiedges=True)

cdef Bool(int x):
    if x:
        return True
    return False

cdef class NFastAutomaton:
    """
    Class :class:`NFastAutomaton`, this class encapsulates a C structure for Automata and
    implement methods to manipulate non-determinist automata.

    INPUT:

    - ``a`` -- automaton must be a :class:`FastAutomaton`

    OUTPUT:

    Return a instance of :class:`NFastAutomaton`.

    EXAMPLES::

        sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
        sage: b = NFastAutomaton(a)
        sage: b
        NFastAutomaton with 4 states and an alphabet of 2 letters

    """
    def __cinit__(self):
        self.a = <NAutomaton *>malloc(sizeof(NAutomaton))
        self.a.e = NULL
        self.a.n = 0
        self.a.na = 0
        self.A = []

    def _initialise_automaton(self, a):
        """
        Transform a determinist  :class:`FastAutomaton` to a non determinist
        :class:`NFastAutomaton`

        INPUT:

        - ``a`` -- automaton must be a :class:`FastAutomaton`

        OUTPUT:

        Return a instance of :class:`NFastAutomaton` initialized with ``a``

        EXAMPLES::
                sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
                sage: b = NFastAutomaton(a)
                sage: b
                NFastAutomaton with 4 states and an alphabet of 2 letters
        """
        if type(a) == FastAutomaton:
            a.copyn(self)
        else:
            raise NotImplementedError("Cannot construct directly a NFastAutomaton for the moment, except from a deterministic one.")
        return self

    def __init__(self, a): # TO DO i=None, final_states=None, A=None
        """
        TESTS::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=1)
            sage: NFastAutomaton(a)
            NFastAutomaton with 4 states and an alphabet of 2 letters
        """
        if a is None:
            pass
        else:
            self = self._initialise_automaton(a)

    def __dealloc__(self):
        FreeNAutomaton(self.a)
        free(self.a)

    def __repr__(self):
        return "NFastAutomaton with %d states and an alphabet of %d letters" % (self.a.n, self.a.na)

    def _latex_(self):
        r"""
        Return a latex representation of the automaton.

        TESTS::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = NFastAutomaton(a)
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
            raise ModuleNotFoundError("dot2tex must be installed in order to have the LaTeX representation of the NFastAutomaton.\n\
                You can install it by doing './sage -i dot2tex' in a shell in the sage directory, or by doing 'install_package(package='dot2tex')' in the notebook.")
        cdef char** ll
        ll = <char **>malloc(sizeof(char*) * self.a.na)
        cdef int i
        strA = []
        for i in range(self.a.na):
            strA.append(str(self.A[i]))
            ll[i] = strA[i]
        sig_on()
        NplotDot(file, self.a[0], ll, "Automaton", sx, sy, False)
        sig_off()
        dotfile = open(file_name)
        return LatexExpr(dot2tex(dotfile.read()))

    @property
    def n_states(self):
        """
        return the numbers of states

        OUTPUT:

        return the numbers of states

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.n_states
            4
        """
        return self.a.n

    def n_succs(self, int i):
        """
        INPUT:

        - ``i`` -- int successor number

        OUTPUT:

        return the numbers of succussor of state ``i``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.n_succs(0)
            1
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        return self.a.e[i].n

    # give the state at end of the jth edge of the state i
    def succ(self, int i, int j):
        """
        Give the state at end of the ``j`` th edge of the state ``i``

        INPUT:

        - ``i`` int state number
        - ``j`` int edge number

        OUTPUT:

        return the state at end of the ``j``th edge of the state ``i``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.succ(1, 0)
            2
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !"%i)
        if j >= self.a.e[i].n or j < 0:
            raise ValueError("The state %s has no edge number %s !"%(i,j))
        return self.a.e[i].a[j].e

    # give the label of the jth edge of the state i
    def label(self, int i, int j):
        """
        Give the label of the ``j``th edge of the state ``i``

        INPUT:

        -``i`` -- int state number
        -``j`` -- int edge number

        OUTPUT:

        return the label index of the ``j``th edge of the state ``i``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
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
        Return ``True``/``False`` if ``i`` state  is/or not  final

        INPUT:

        -``i`` -- int state number

        OUTPUT:

        Return ``True``/``False`` if ``i`` state  is/or not  final

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.is_final(1)
            True
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        return Bool(self.a.e[i].final)

    def is_initial(self, int i):
        """
        Return `True``/``False`` if ``i`` state  is/or not  initial

        INPUT:

        -``i`` -- int state number

        OUTPUT:

        Return ``True``/``False`` if ``i`` state  is/or not  initial

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'),(1, 2, 'c'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.is_initial(1)
            False
        """
        if i >= self.a.n or i < 0:
            raise ValueError("There is no state %s !" % i)
        return Bool(self.a.e[i].initial)

    @property
    def initial_states(self):
        """
        Get the initial state :class:`NFastAutomaton` attribut

        OUTPUT:

        Return the initial state ``i``  of  :class:`NFastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = NFastAutomaton(a)
            sage: b.initial_states
            []
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], i=2)
            sage: b = NFastAutomaton(a)
            sage: b.initial_states
            [2]
        """
        l = []
        for i in range(self.a.n):
            if self.a.e[i].initial:
                l.append(i)
        return l

    def final_states(self):
        """
        Indicate all final states

        OUTPUT:

        Return the list of final states

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = NFastAutomaton(a)
            sage: b.final_states()
            [0, 1, 2, 3]
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], )
            sage: b = NFastAutomaton(a)
            sage: b.final_states()
            [0, 1, 2, 3]
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
            sage: b = NFastAutomaton(a)
            sage: b.final_states()
            [0, 3]
        """
        l = []
        for i in range(self.a.n):
            if self.a.e[i].final:
                l.append(i)
        return l

    @property
    def Alphabet(self):
        """
        To get the :class:`NFastAutomaton` attribut Alphabet

        OUTPUT:

        Return a the alphabet ``A`` of  :class:`NFastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = NFastAutomaton(a)
            sage: b.Alphabet
            ['a', 'b']

        """
        return self.A

    def set_initial_state(self, int i, bool initial=True):
        """
        Set the initial state.

        INPUT:

        - ``i` -- int the initial state of the automaton
        - ``initial`` -- (default: ``True``) in the case is initial

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = NFastAutomaton(a)
            sage: b.set_initial_state(2)
            sage: b.initial_states
            [2]
            sage: b.set_initial_state(6)
            Traceback (most recent call last):
            ...
            ValueError: initial state must be a current state : 6 not in [-1, 3]

        """
        if i < 0 or i >= self.a.n:
            raise ValueError("initial state must be a current state : " +
                             "%d not in [-1, %d]" % (i, self.a.n - 1))
        self.a.e[i].initial = initial

    def add_edge(self, i, l, f):
        """
        Add a edge in the automaton

        INPUT:

        - ``i`` -- the first state
        - ``l`` -- the label of edge
        - ``f`` -- the second state


        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.add_edge(2,'a',1)
            sage: b.add_edge(2,'v',1)
            Traceback (most recent call last):
            ...
            ValueError: The letter v doesn't exist.
            sage: b.add_edge(2,'v',6)
            Traceback (most recent call last):
            ValueError: The state  6 doesn't exist.
            sage: b.add_edge(5,'v',6)
            Traceback (most recent call last):
            ValueError: The state  5 doesn't exist.

        """
        if i >= self.a.n:
            raise ValueError("The state %s doesn't exist." % i)
        if f >= self.a.n:
            raise ValueError("The state  %s doesn't exist." % f)
        try:
            k = self.A.index(l)
        except:
            # La lettre %s n'existe pas.
            raise ValueError("The letter %s doesn't exist." % l)

        sig_on()
        AddEdgeN(self.a, i, f, k)
        sig_off()

    def add_state(self, bool final):
        """
        Add a state in the automaton

        INPUT:

        - ``final`` -- boolean indicate if the added state is final

        OUTPUT:

        return the numbers of states

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.add_state(True)  # not implemented
            TypeError                                 Traceback (most recent call last)
            ...
            TypeError: 'NotImplementedType' object is not callable
        """
        raise NotImplemented()

    def add_path(self, int e, int f, list li, verb=False):
        """
        Add a path between states ``e`` and ``f``
        of :class:`NFastAutomaton` following ``li``

        INPUT:

        - ``e`` -- int the input state
        - ``f`` -- int the final state 
        - ``li`` -- list of states
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode


        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a'), (2,3,'b')], i=2)
            sage: b = NFastAutomaton(a)
            sage: b.add_path(1, 2, [1])

        """
        cdef int *l = <int *>malloc(sizeof(int)*len(li));
        for i in range(len(li)):
            l[i] = li[i]
        sig_on()
        AddPathN(self.a, e, f, l, len(li), verb)
        sig_off()

    def determinize(self, puits=False, verb=0):
        """
        Determines a non determinist automaton with the same alphabet of ``self``

        INPUT:

        - ``puits``  -- (default: ``False``)
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return a non determinist automaton  :class:`NFastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: b.determinize()
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef Automaton a
        sig_on()
        r = FastAutomaton(None)
        a = DeterminizeN(self.a[0], puits, verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        return r

    def plot(self, int sx=10, int sy=8, verb=False):
        """
        plot a representation of the :class:`NFastAutomaton`.

        INPUT:

        - ``sx`` -- int (default: 10)
        - ``sy`` -- int (default: 8)
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: #g = b.plot()   # random

        .. PLOT::

            a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            b = NFastAutomaton(a)
            sphinx_plot(b)

        """
        cdef char** ll
        cdef int i
        cdef char *file
        if DotExists ():
            ll = <char **>malloc(sizeof(char*) * self.a.na)
            strA = []
            for i in range(self.a.na):
                strA.append(str(self.A[i]))
                ll[i] = strA[i]
            from sage.misc.temporary_file import tmp_filename
            file_name = tmp_filename()
            file = file_name
            if verb:
                print("file=%s" % file_name)
            sig_on()
            NplotDot(file, self.a[0], ll, "Automaton", sx, sy, True)
            free(ll)
            sig_off()
            from PIL import Image
            return Image.open(file_name+'.png')
        else:
            raise NotImplementedError("You cannot plot the NFastAutomaton without dot. Install the dot command of the GraphViz package.")

# cdef set_FastAutomaton (FastAutomaton a, Automaton a2):
#    a.a[0] = a2

cdef class FastAutomaton:
    r"""
    Class :class:`FastAutomaton`, this class encapsulates a C structure for Automata and 
    implement methods to manipulate determinist automata.

    EXAMPLES::

        sage: FastAutomaton([(0,1,'a') ,(2,3,'b')])
        FastAutomaton with 4 states and an alphabet of 2 letters
        sage: d = DiGraph({0: [1,2,3], 1: [0,2]})
        sage: FastAutomaton(d)
        FastAutomaton with 4 states and an alphabet of 1 letters
        sage: g = DiGraph({0:{1:'x',2:'z',3:'a'}, 2:{5:'o'}})
        sage: FastAutomaton(g)
        FastAutomaton with 5 states and an alphabet of 4 letters
        sage: a = FastAutomaton([(0, 1,'a') ,(2, 3,'b')], i = 2)
        sage: a
        FastAutomaton with 4 states and an alphabet of 2 letters
        sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
        sage: a
        FastAutomaton with 4 states and an alphabet of 2 letters
    """

#    cdef Automaton* a
#    cdef list A

    def __cinit__(self):
        """

        """
        # print("cinit")
        self.a = <Automaton *>malloc(sizeof(Automaton))
        # initialise
        self.a.e = NULL
        self.a.n = 0
        self.a.na = 0
        self.a.i = -1
        self.A = []
        self.S = None
        self.dA = None
        self.dS = None

    def __init__(self, a, i=None, final_states=None, A=None, keep_S=True):
        r"""
        INPUT:

        - ``i`` - (default: None) - initial state

        - ``final_states`` - (default: None) - list of final states

        - ``A`` - (default: None) - Alphabet

        - ``keep_S``- (default: `True`) - Keep labels of states or not

        OUTPUT:

        Return a instance of :class:`FastAutomaton`.

        TESTS::

            sage: FastAutomaton([(0,1,'a') ,(2,3,'b')])
            FastAutomaton with 4 states and an alphabet of 2 letters
        """
        # print("init")
        if a is None:
            return
        from sage.graphs.digraph import DiGraph
        if isinstance(a, list):
            a = DiGraph(a, multiedges=True, loops=True)
        if isinstance(a, DiGraph):
            if A is None:
                A = list(set(a.edge_labels()))
            self.A = A
            self.a[0] = getAutomaton(a, initial=i, F=final_states, A=self.A)
            self.dA = a.dA
            if keep_S:
                self.S = a.S
                self.dS = a.dS
        elif isinstance(a, FastAutomaton):
            self = a
        else:
            raise ValueError("Cannot convert the input to FastAutomaton.")

    def __dealloc__(self):
        """
        Desalloc  Automaton  Overwrite built-in function

        TESTS::

        """
        # print("free (%s etats) "%self.a.n)
        sig_on()
        FreeAutomaton(self.a)
        # print("free self.a")
        free(self.a)
        sig_off()

    def __repr__(self):
        """
        Return a representation of automaton,  Overwrite built-in function

        OUTPUT:

        Return a representation of automaton

        TESTS:

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: repr(a)
            'FastAutomaton with 4 states and an alphabet of 2 letters'

        """
        return "FastAutomaton with %d states and an alphabet of %d letters" % (self.a.n, self.a.na)

    def _latex_(self):
        r"""
        Return a latex representation of the automaton.

        OUTPUT:

        string - latex representation of the automaton.

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
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
            print("dot2tex must be installed in order to have the LaTeX representation of the FastAutomaton.")
            print("You can install it by doing './sage -i dot2tex' in a shell in the sage directory, or by doing 'install_package(package='dot2tex')' in the notebook.")
            return None
        cdef char** ll # labels of edges
        cdef char** vl # labels of vertices
        cdef int i
        ll = <char **>malloc(sizeof(char*) * self.a.na)
        if vlabels is None:
            vl = NULL
        else:
            if verb:
                print("alloc %s..."%self.a.n)
            vl = <char **>malloc(sizeof(char*) * self.a.n)
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
        free(ll)
        if vlabels is not None:
            free(vl)
        dotfile = open(file_name)
        return LatexExpr(dot2tex(dotfile.read()))

    def __hash__(self):
        r"""
        Hash automaton,  Overwrite built-in function

        OUTPUT:

        Return the hash code of the :class:`FastAutomaton`

        TESTS::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: hash(a)
            -3816799034168020408
        """
        h = hash(tuple(self.A))
        sig_on()
        h += hashAutomaton(self.a[0])
        sig_off()
        # print "hash=%s"%h
        return h
    
    def _richcmp_(self, FastAutomaton other, int op):
        r"""
        Compare function, Overwrite built-in function

        INPUT:

        - ``other`` -- other :class:`FastAutomaton` to compare

        OUTPUT:

        Return the result of test (``True`` or ``False``)

        TESTS::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: b = FastAutomaton([(0, 1, 'a'),(1,2,'c')], i=0)
            sage: a == b
            False
        """
        from sage.structure.richcmp import (op_EQ, op_NE)
        #(rich_to_bool,
        #op_EQ, op_NE, op_LT, op_LE, op_GT, op_GE)
        cdef int r
        if op != op_EQ and op != op_NE:
            raise NotImplementedError("Comparaison <, >, <= or >= not implemented for FastAutomata.")
        sig_on()
        r = equalsAutomaton(self.a[0], other.a[0]) and self.A == other.A
        sig_off()
        if op == op_EQ:
            return (r == 1)
        else:
            return (r == 0)

#     def Automaton(self):
#         return AutomatonGet(self.a[0], self.A)

#    cdef set_a(self, Automaton a):
#        self.a[0] = a

    #give a Sage Automon from the FastAutomaton
    def get_automaton(self):
        return AutomatonToSageAutomaton(self.a[0], self.A)
    
    #give a Sage Automon from the FastAutomaton
    def get_DiGraph(self):
        return AutomatonToDiGraph(self.a[0], self.A)
    
    # give a FastAutomaton recognizing the full language over A.
    def full(self, list A):
        """
        Give a :class:`FastAutomaton` recognizing the full language over A.

        INPUT:

        - ``A`` -- list of letters of alphabet

        OUTPUT:

        Return a :class:`FastAutomaton` recognizing the full language over A.

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a
            FastAutomaton with 4 states and an alphabet of 2 letters
            sage: a.full(['a'])
            FastAutomaton with 1 states and an alphabet of 1 letters
            sage: a.full(['a','b'])
            FastAutomaton with 1 states and an alphabet of 2 letters
            sage: a.full(['a','b','c'])
            FastAutomaton with 1 states and an alphabet of 3 letters
        """
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = NewAutomaton(1, len(A))
        sig_off()
        for i in range(len(A)):
            a.e[0].f[i] = 0
        a.e[0].final = True
        a.i = 0
        r.a[0] = a
        r.A = A

        return r

    def plot(self, int sx=10, int sy=8, vlabels=None, html=False, file=None, bool draw=True, verb=False):
        """
        Plot the :class:`FastAutomaton`. Draw using the dot command, if installed on the platform.

        It is recommanded to install the dot command of the Graphviz package in your system in order to get a nice picture. Otherwise it will draw using the function plot of :class:`Automaton` of Sage.

        INPUT:

        - ``sx`` - int (default: 10) - width of the picture
        - ``sy`` - int (default: 8) - height of the picture
        - ``vlabels`` - (default: None) - labels of the vertices
        - ``html`` - (default: ``False``) - tell if dot should draw vertices in html mode
        - ``file`` - (default: ``None``) - the adress of the .dot file of the drawing (only if dot is installed)
        - ``draw`` - (default: ``True``) - if False, only generate the .dot file (only if dot is installed)
        - ``verb`` - (default: ``False``) - active or not the verbose mode

        TESTS::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.plot()  # random
            <PIL.PngImagePlugin.PngImageFile image mode=RGBA size=189x147 at 0x7FD4B6D94390>

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: #a.plot()   # random

        .. PLOT::

            a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sphinx_plot(a)
        """
        cdef char *cfile
        cdef char** ll # labels of edges
        cdef char** vl # labels of vertices
        cdef int i
        if DotExists():
            from sage.misc.temporary_file import tmp_filename
            file_name = tmp_filename()+".dot"
            file = file_name
            ll = <char **>malloc(sizeof(char*) * self.a.na)
            if vlabels is None:
                if self.S is not None:
                    if verb:
                        print("alloc %s..." % self.a.n)
                    vl = <char **>malloc(sizeof(char*) * self.a.n)
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
                    vl = NULL
            else:
                if verb:
                    print("alloc %s..." % self.a.n)
                vl = <char **>malloc(sizeof(char*) * self.a.n)
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
            plotDot(file, self.a[0], ll, "Automaton", sx, sy, vl, html, verb, True)
            sig_off()
            if verb:
                print("free...plot")
            free(ll)
            if vlabels is not None:
                free(vl)
            from PIL import Image
            return Image.open(file_name+'.png')
        else:
            return AutomatonToSageAutomaton(self.a[0], self.A).plot()
            #raise NotImplementedError("You cannot plot the FastAutomaton without dot. Install the dot command of the GraphViz package.")

    @property
    def Alphabet(self):
        """
        To get the :class:`FastAutomaton` attribut Alphabet

        OUTPUT:

        Return a the alphabet ``A`` of  :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.Alphabet
            ['a', 'b']
        """
        return self.A

    def setAlphabet(self, list A):
        """
        Set the alphabet

        INPUT:

        - ``A`` -- list of letters of alphabet

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.setAlphabet(['a', 'b', 'c'])
            sage: a.Alphabet
            ['a', 'b', 'c']
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.setAlphabet(['a','e'])
            sage: a.Alphabet
            ['a', 'e']
        """
        self.A = A
        self.a[0].na = len(A)

    @property
    def initial_state(self):
        """
        Get the initial state :class:`FastAutomaton` attribut

        OUTPUT:

        Return the initial state ``i``  of  :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.initial_state
            -1
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], i=2)
            sage: a.initial_state
            2
        """
        return self.a.i

    def set_initial_state(self, int i):
        """
        Set the initial state.

        INPUT:

        - ``i`` -- int the initial state of the automaton

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_initial_state(2)
            sage: a.initial_state
            2
            sage: a.set_initial_state(6)
            Traceback (most recent call last):
            ...
            ValueError: initial state must be a current state : 6 not in [-1, 3]
        """
        if i < self.a.n and i >= -1:
            self.a.i = i
        else:
            raise ValueError("initial state must be a current state : " +
                             "%d not in [-1, %d]" % (i, self.a.n - 1))

    def final_states(self):
        """
        Indicate all final states

        OUTPUT:

        Return the list of final states

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.final_states()
            [0, 1, 2, 3]
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], )
            sage: a.final_states()
            [0, 1, 2, 3]
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
            sage: a.final_states()
            [0, 3]
        """
        l = []
        for i in range(self.a.n):
            if self.a.e[i].final:
                l.append(i)
        return l

    def states(self):
        """
        Indicate all states of the automaton

        OUTPUT:

        Return the list of states

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.states()
            [0, 1, 2, 3]
        """
        if self.S is None:
            return range(self.a.n)
        else:
            return self.S

    def set_final_states(self, lf):
        """
        Set the final states.

         INPUT:

        - ``lf`` -- list of states to set as final

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_final_states([0,3])
            sage: a.final_states()
            [0, 3]
            sage: a.set_final_states([0,4])
            Traceback (most recent call last):
            ...
            ValueError: 4 is not a state !

        """
        cdef int f
        for f in range(self.a.n):
            self.a.e[f].final = 0
        for f in lf:
            if f < 0 or f >= self.a.n:
                raise ValueError("%d is not a state !" % f)
            self.a.e[f].final = 1

    def is_final(self, int e):
        """
        Indicate if the state is final

         INPUT:

        - ``e`` -- int input state to examine as final

        OUTPUT:

        ``True`` if the state ``e`` is final (i.e. ``False`` in the other case)

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.is_final(3)
            True
            sage: a.is_final(4)
            False
            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')], final_states=[0,3])
            sage: a.is_final(2)
            False
        """
        if e >= 0 and e < self.a.n:
            return Bool(self.a.e[e].final)
        else:
            return False

    def set_final_state(self, int e, final=True):
        """
        Set the final state.

         INPUT:

        - ``e`` -- int state to set as final
        - ``final`` -- (default: ``True``) in the case is final

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a') ,(2,3,'b')])
            sage: a.set_final_state(3)
            sage: a.final_states()
            [0, 1, 2, 3]
            sage: a.set_final_state(4)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not a state !
        """
        if e >= 0 and e < self.a.n:
            self.a.e[e].final = final
        else:
            raise ValueError("%d is not a state !" % e)

    def succ(self, int i, int j):
        """
        Return the state reached by following the edge j from state i.

        INPUT:

        - ``i`` - int - number of the input state
        - ``j`` - int - number of the edge

        OUTPUT:

        return successor of state

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a'), (2,3,'b')])
            sage: a.succ(0, 1)
            -1
            sage: a.succ(2,1)
            3
        """
        if i < 0 or i >= self.a.n or j < 0 or j >= self.a.na:
            return -1
        return self.a.e[i].f[j]

    # donne les fils de l'état i
    def succs(self, int i):
        """
        return lines of state ``i``.

        INPUT:

        - ``i`` -- int the input state

        OUTPUT:

        return lines of state ``i``

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a'), (2,3,'b')])
            sage: a.succs(2)
            [1]
            sage: a.succs(4)
            []

        """
#        if i is None:
#            i = self.a.i
#        el
        if i < 0 or i >= self.a.n:
            return []
        return [j for j in range(self.a.na) if self.a.e[i].f[j] != -1]

    # suit le chemin étiqueté par l et rend l'état atteint
    def path(self, list l, i=None):
        """
        Follows the path labeled by ``l`` and return the reached state

        INPUT:

        - ``l`` -- list indicate the  way label
        - ``i`` -- (default: ``None``) the initial state

        OUTPUT:

        return the state reached after the following way

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a'), (2,3,'b')], i=2)
            sage: a.path([1])
            3
            sage: a.path([0, 2])
            -1
        """
        if i is None:
            i = self.a.i
        for j in l:
            i = self.succ(i, j)
        return i

    def set_succ(self, int i, int j, int k):
        """
        Set the successor state

        INPUT:

        - ``i`` -- int the input state
        - ``j`` -- int the output state

        EXAMPLES::

            sage: a = FastAutomaton([(0,1,'a'), (2, 3,'b')], i=2)
            sage: a.set_succ(0, 1, 2)
            sage: a.succs(0)
            [0, 1]
            sage: a.set_succ(0, 4, 2)
            Traceback (most recent call last):
            ...
            ValueError: set_succ(0, 4) : index out of bounds !

        """
        if i < 0 or i >= self.a.n or j < 0 or j >= self.a.na:
            raise ValueError("set_succ(%s, %s) : index out of bounds !" % (i, j))
        self.a.e[i].f[j] = k

    def zero_completeOP(self, verb=False):
        """
        zero-complete automaton

        INPUT:

        - ``verb`` -- (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the zero-complete :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.zero_completeOP(True)
            l0 = 0
            state 0 ..
            state 1 ..
            state 2 ..
        """
        sig_on()
        ZeroComplete(self.a, list(self.A).index(self.A[0]), verb)
        sig_off()

    def zero_complete2(self, etat_puits=False, verb=False):
        """
        zero-complete automaton in the other way

        INPUT:

        - ``etat_puits`` --  (default: ``False``)
        - ``verb`` -- (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the zero-complete :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.zero_complete2(True)
            FastAutomaton with 2 states and an alphabet of 2 letters

        """
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = ZeroComplete2(self.a, list(self.A).index(self.A[0]), etat_puits, verb)
        sig_off()
        r.a[0] = a
        r.A = self.A

        return r.prune().minimize()

    def zero_inv(self, z=0, simplify=True):
        """
        Compute an automaton recognizing the language (A[z]*)L, where L is the language of self.

        INPUT:

        - ``z`` - (default: 0) - index of the letter
        - ``simplify`` - (default: True) - if true, prune and minimize the result

        OUTPUT:

        Return a :class:`FastAutomaton` whose language is the set of words of
        the language of a with some zeroes on the beggining
        (we call zero the letter of index z).

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.zero_inv(0)
            FastAutomaton with 2 states and an alphabet of 2 letters
            sage: a.zero_inv(1)
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = ZeroInv(self.a, list(self.A).index(self.A[z]))
        sig_off()
        r.a[0] = a
        r.A = self.A
        if simplify:
            return r.prune().minimize()
        else:
            return r

    # change the final states of the automaton
    # new final states are the one in a strongly connected component containing a final state, others states are not final
    # this function can be accelerated
    def prune_inf2OP(self, verb=False):
        """
        Compute the emondation of the automaton
        change the final states of the automaton
        new final states are the one in a strongly
        connected component containing a final state,
        others states are not final
        this function can be accelerated

        INPUT:

        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the pruned :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.prune_inf2OP(True)
            sage: a
            FastAutomaton with 3 states and an alphabet of 2 letters
        """
        cc = self.strongly_connected_components()
        f = []
        for c in cc:
            # test que l'on peut boucler dans cette composante
            ok = False
            for i in range(self.a.na):
                if self.a.e[c[0]].f[i] in c:
                    ok = True
                    break
            if not ok:
                continue
            for i in c:
                if self.a.e[i].final:
                    f += c
                    break
        self.set_final_states(f)

    # new final states are the ones in strongly connected components
    def prune_inf(self, verb=False):
        """
        Compute the emondation of the automaton
        remove all states from which there no infinite way

        INPUT:

        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the pruned :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b')], i=0)
            sage: a.prune_inf(True)
            recurrence...
            States counter = 0
            count...
            cpt = 0
            final states...
            FastAutomaton with 0 states and an alphabet of 2 letters

        """
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = prune_inf(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        return r

    def prune_i(self, verb=False):
        """
        Compute the emondation of the automaton
        remove all not accessible states

        INPUT:

        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the pruned :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b'), (0, 3, 'b')], i=0)
            sage: a.prune_i(True)
            deleted States : [ ]
            FastAutomaton with 3 states and an alphabet of 2 letters
            sage: a = FastAutomaton([(0, 1, 'a'), (0, 3, 'b'), (0, 3, 'b')])
            sage: a.prune_i(True)
            FastAutomaton with 0 states and an alphabet of 2 letters
        """
        if self.initial_state == -1:
            empty = FastAutomaton([])
            empty.setAlphabet(self.Alphabet)
            return empty
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = pruneI(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        return r

    def prune(self, verb=False):
        """
        Compute the emondation of the automaton
        remove all not accessible and not co-accessible states

        INPUT:

        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the pruned :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.prune(True)
            4 components : [ 1 0 3 2 ]
            0 : [ 1 ]
            1 : [ 0 ]
            2 : [ 3 ]
            3 : [ 2 ]
            0 co-acc
            1 co-acc
            2 co-acc
            3 co-acc
            rec...
            l : [ 0(7) 1(7) -1(5) -1(5) ]
            create the new automaton 2 2...
            pass 2
            pass 3
            deleted States : [ 2( non-acc ) 3( non-acc ) ]
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = prune(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        return r

#    def equals (self, FastAutomaton b):
#        return Bool(equalsAutomaton(self.a[0], b.a[0]))

    # assume that the dictionnary d is injective !!!
    def product(self, FastAutomaton b, dict d=None, verb=False):
        """
        Give the product of the :class:`FastAutomaton` and ``a`` an other
        ``FastAutomaton``.

        INPUT:

        - ``a`` -- :class:`FastAutomaton` to multiply
        - ``d`` -- dict (default: ``None``) dictionary to translate
          language of automaton
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the product as a :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'c'), (1, 2, 'd')], i=2)
            sage: a.product(b)
            FastAutomaton with 12 states and an alphabet of 4 letters
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
            initial State 0.
            Automaton with 3 states, 2 letters.
            0 --1--> 1
            2 --0--> 1
            initial State 1.
            FastAutomaton with 12 states and an alphabet of 4 letters
        """
        if d is None:
            d = {}
            for la in self.A:
                for lb in b.A:
                    d[(la, lb)] = (la, lb)
            if verb:
                print(d)
        cdef Automaton a
        cdef Dict dC
        r = FastAutomaton(None)
        Av = []
        sig_on()
        dv = imagProductDict(d, self.A, b.A, Av=Av)
        sig_off()
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        sig_on()
        dC = getProductDict(d, self.A, b.A, dv=dv, verb=verb)
        sig_off()
        if verb:
            print("dC=")
            printDict(dC)
        sig_on()
        a = Product(self.a[0], b.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        r.a[0] = a
        r.A = Av

        return r

    def intersection(self, FastAutomaton a, verb=False, simplify=True):
        """
        Give the intersection of the :class:`FastAutomaton` and ``a`` an other
        :class:`FastAutomaton`.

        INPUT:

        - ``a`` -- :class:`FastAutomaton` to intersect
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode
        - ``simplify`` - (default: ``True``) - if True, prune and minimize the result

        OUTPUT:

        Return the intersected :class:`FastAutomaton`.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'c'), (1, 2, 'd')], i=2)
            sage: a.intersection(b)
            FastAutomaton with 1 states and an alphabet of 0 letters
            sage: a.intersection(b, simplify=False)
            FastAutomaton with 12 states and an alphabet of 0 letters
        """
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l, l)] = l
        if verb:
            print("d=%s" % d)
        p = self.product(a, d, verb=verb)
        if simplify:
            return p.prune().minimize()
        else:
            return p

    # determine if the automaton is complete (i.e. with his hole state)
    def is_complete(self):
        """
        Determine if the automaton is complete (i.e. with his hole state)

        OUTPUT:

        Return ``True`` if the automaton is complete (i.e. ``False`` if not)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.is_complete()
            False
            sage: a = FastAutomaton([(0, 0, 'a')])
            sage: a.is_complete()
            True
        """
        sig_on()
        res = IsCompleteAutomaton(self.a[0])
        sig_off()
        return Bool(res)

    # give a complete automaton (i.e. with his hole state)
    def complete(self):
        """
        Give the complete automaton (i.e. with his hole state)

        OUTPUT:

        Return ``True`` if the automaton is complete (i.e. ``False`` if not)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complete()
            True
        """
        sig_on()
        res = CompleteAutomaton(self.a)
        sig_off()
        return Bool(res)

    # give the smallest language stable by prefix containing the language of self
    # i.e. every states begin finals
    def prefix_closure(self):
        """
        give the smallest language stable by prefix containing the language of self
        i.e. every states begin finals

        OUTPUT:

        Return the smallest language :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.prefix_closure()
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef int i
        cdef Automaton a
        r = FastAutomaton(None)
        sig_on()
        a = prune(self.a[0], False)
        sig_off()
        r.a[0] = a
        r.A = self.A
        for i in range(a.n):
            a.e[i].final = True
        return r

    # FastAutomaton
    def union(self, FastAutomaton a, simplify=True, verb=False):
        """
        re-union :class:`FastAutomaton`

        INPUT:

        - ``a`` -- :class:`FastAutomaton` to union
        - ``simplify`` --  (default: ``True``) if simplification
          is required
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return the union :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.union(b, verb=True)
            Av=['a']
            dv={'a': 0}
            {'a': 0, 'b': 1}
            {'a': 0, 'd': 1}
            Keys=[('a', 'a')]
            dC=
            [ 0 -1 -1 -1 ]
            Automaton with 5 states, 2 letters.
            0 --0--> 1
            0 --1--> 4
            1 --0--> 4
            1 --1--> 4
            2 --0--> 4
            2 --1--> 3
            3 --0--> 4
            3 --1--> 4
            4 --0--> 4
            4 --1--> 4
            initial State 0.
            Automaton with 4 states, 2 letters.
            0 --0--> 3
            0 --1--> 1
            1 --0--> 3
            1 --1--> 3
            2 --0--> 1
            2 --1--> 3
            3 --0--> 3
            3 --1--> 3
            initial State 1.
            FastAutomaton with 2 states and an alphabet of 1 letters
        """
        # complete the automata
        sig_on()
        CompleteAutomaton(self.a)
        CompleteAutomaton(a.a)
        sig_off()

        # make the product
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l, l)] = l

        cdef Automaton ap
        cdef Dict dC
        r = FastAutomaton(None)
        Av = []
        sig_on()
        dv = imagProductDict(d, self.A, a.A, Av=Av)
        sig_off()
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        sig_on()
        dC = getProductDict(d, self.A, a.A, dv=dv, verb=verb)
        sig_off()
        sig_on()
        if verb:
            print("dC=")
            printDict(dC)
        ap = Product(self.a[0], a.a[0], dC, verb)
        FreeDict(&dC)
        sig_off()
        # set final states
        cdef int i, j
        cdef n1 = self.a.n
        for i in range(n1):
            for j in range(a.a.n):
                ap.e[i + n1 * j].final = self.a.e[i].final or a.a.e[j].final

        r.a[0] = ap
        r.A = Av
        if simplify:
            return r.prune().minimize()
        else:
            return r

    # split the automaton with respect to a FastAutomaton
    def split (self, FastAutomaton a, simplify=True, verb=False):
        """
        Split the automaton with respect to ``a`` a :class:`FastAutomaton`

        INPUT:

        - ``a`` -- :class:`FastAutomaton` in respect to split
        - ``simplify`` - (default: True) - if True, prune and minimize the result
        - ``verb`` - (default: False) - if True, display informations for debugging

        OUTPUT:

        Return tuple of two splited automaton, recognizing respectively
        the language intersection of L and La, and the language
        intersection of L and complementary of La,
        where L is the language of self and La is the language of a.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.split(b, verb=True)
            Av=['a']
            dv={'a': 0}
            {'a': 0, 'b': 1}
            {'a': 0, 'd': 1}
            Keys=[('a', 'a')]
            dC=
            [ 0 -1 -1 -1 ]
            Automaton with 4 states, 2 letters.
            0 --0--> 1
            2 --1--> 3
            initial State 0.
            Automaton with 4 states, 2 letters.
            0 --0--> 3
            0 --1--> 1
            1 --0--> 3
            1 --1--> 3
            2 --0--> 1
            2 --1--> 3
            3 --0--> 3
            3 --1--> 3
            initial State 1.
            [FastAutomaton with 1 states and an alphabet of 1 letters,
             FastAutomaton with 2 states and an alphabet of 1 letters]
        """
        # complete the automaton a
        sig_on()
        CompleteAutomaton(a.a)
        sig_off()
        # make the product
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l, l)] = l

        cdef Automaton ap
        cdef Dict dC
        r = FastAutomaton(None)
        r2 = FastAutomaton(None)
        Av = []
        sig_on()
        dv = imagProductDict(d, self.A, a.A, Av=Av)
        sig_off()
        if verb:
            print("Av=%s" % Av)
            print("dv=%s" % dv)
        sig_on()
        dC = getProductDict(d, self.A, a.A, dv=dv, verb=verb)

        if verb:
            print("dC=")
            printDict(dC)

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
        cdef Automaton ap2
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

    # modify the automaton to recognize the langage shifted by l (letter given by its index)
    def shift1OP (self, int l, verb=False):
        """
        Shift the automaton ON PLACE to recognize the language shifted
        by ``l`` (letter given by its index).
        The new language is the language of words u such that lu was recognized by self.

        INPUT:

        - ``l`` -- int  index of letter to shift
        - ``verb`` - boolean (default: ``False``) - if True, display informations for debugging

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.shift1OP(0, verb=True)
            sage: a.initial_state
            1
        """
        if self.a.i != -1:
            self.a.i = self.a.e[self.a.i].f[l]

    # modify the automaton to recognize the langage shifted by l^np (where l is a letter given by its index)
    def shiftOP(self, l, int np, verb=False):
        """
        Shift the automaton ON PLACE to recognize the language shifted ``np`` times
        by l (letter given by its index).
        The new language is the language of words u such that (l^np)u was recognized by self.

        INPUT:

        - ``l`` -- int  index of letter to shift
        - ``np`` -- int  number of shift
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` to activate the verbose mode

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.shiftOP(0, 2)
            sage: a.initial_state
            -1
        """
        cdef int i
        for i in range(np):
            if self.a.i != -1:
                self.a.i = self.a.e[self.a.i].f[l]

    def unshift1(self, int l, final=False):
        """
        Unshift the automaton to recognize the language shifted
        by ``l`` (letter given by its index).
        The new language is the languages of words lu, where u is recognized by self.

        INPUT:

        - ``l`` -- int  index of letter to shift
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` to activate the verbose mode

        OUTPUT:

        Return a :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.unshift1(1)
            FastAutomaton with 5 states and an alphabet of 2 letters
        """
        r = FastAutomaton(None)
        cdef Automaton aut
        sig_on()
        aut = CopyAutomaton(self.a[0], self.a.n+1, self.a.na)
        sig_off()
        cdef int i
        cdef int ne = self.a.n
        for i in range(aut.na):
            aut.e[ne].f[i] = -1
        aut.e[ne].f[l] = self.a.i
        aut.e[ne].final = final
        aut.i = ne
        r.a[0] = aut
        r.A = self.A
        return r

    # this function could be written in a more efficient way
    def unshiftl(self, list l):
        """
        Return a new automaton whose language is the set of words wu,
        where u is recognized by self, and w is the word
        corresponding to the list of indices of letters l.  
        
        INPUT:

        - ``l`` -- list of indices of letters

        OUTPUT:

        Return a unshifted  :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.initial_state
            0
            sage: a.shiftOP(0, 2)
            sage: a.unshiftl([0, 1])
            FastAutomaton with 6 states and an alphabet of 2 letters
        """
        a = self
        l.reverse()
        for i in l:
            a = a.unshift1(i)
        l.reverse()
        return a

    def unshift(self, int l, int np, final=None):
        """
        Return a new automaton whose language is the set of words wu,
        where u is recognized by self, and w is the word
        corresponding to np times the letter of index l.

        INPUT:

        - ``l`` -- int  
        - ``np``  --  int  
        - ``final`` -- (default: ``None``) if final or not

        OUTPUT:

        Return a :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.unshift(0, 2)
            FastAutomaton with 6 states and an alphabet of 2 letters
        """
        if np == 0:
            return self
        if final is None:
            if self.a.i == -1:
                r = FastAutomaton(None)
                r.A = self.A
                return r
            final = self.a.e[self.a.i].final
        r = FastAutomaton(None)
        cdef Automaton aut
        sig_on()
        aut = CopyAutomaton(self.a[0], self.a.n+np, self.a.na)
        sig_off()
        cdef int i
        cdef int ne = self.a.n
        for j in range(np):
            for i in range(aut.na):
                aut.e[ne+j].f[i] = -1
            if j > 0:
                aut.e[ne+j].f[l] = ne+j-1
            else:
                aut.e[ne+j].f[l] = self.a.i
            aut.e[ne+j].final = final
        aut.i = ne+np-1
        r.a[0] = aut
        r.A = self.A
        return r

    def copyn(self, NFastAutomaton r, verb=False):
        """
        Convert  a determinist automaton :class:`FastAutomaton` to
        a non determinist automaton :class:`NFastAutomaton`

        INPUT:

        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode
        - ``r`` -- :class:`NFastAutomaton` to replace

        OUTPUT:

        Return the :class:`NFastAutomaton` copy like of the :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = NFastAutomaton(a)
            sage: a.copyn(b)
            NFastAutomaton with 4 states and an alphabet of 2 letters
        """
        cdef NAutomaton a
        sig_on()
        a = CopyN(self.a[0], verb)
        sig_off()
        r.a[0] = a
        r.A = self.A
        return r

    def concat(self, FastAutomaton b, det=True, simplify=True, verb=False):
        """
        Concatenates :class:`FastAutomaton`, ``b``

        INPUT:

        - ``b`` -- :class:`FastAutomaton`  to concatenate
        - ``det``  -- (default: ``True``) - if True, determinize
        - ``simplify`` -- (default: ``True``) - if True and if det=True, prune and minimize
        - ``verb`` -- boolean (default: ``False``) fix
          to ``True`` for activation the verbose mode

        OUTPUT:

        Return a :class:`NFastAutomaton` (if ``det``=``False``)
        or  :class:`FastAutomaton` (if ``det``=``True``)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.concat(b)
            FastAutomaton with 2 states and an alphabet of 3 letters
            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.concat(b, det=False)
            NFastAutomaton with 7 states and an alphabet of 3 letters
        """
        cdef FastAutomaton a
        if self.A != b.A:
            A = list(set(self.A).union(set(b.A)))
            if verb:
                print("Alphabet Changing (%s, %s -> %s)..." %(self.A, b.A, A))
            a = self.bigger_alphabet(A)
            b = b.bigger_alphabet(A)
            # raise ValueError("Error : concatenation of automaton having differents alphabets.")
        else:
            a = self
            A = self.A
        if verb:
            print("a=%s (A=%s)\nb=%s (A=%s)" % (a, a.A, b, b.A))
        cdef NAutomaton na
        r = NFastAutomaton(None)
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

    def proj(self, dict d, det=True, simplify=True, verb=False):
        """
        Project following the dictionary ``d``.
        Give an automaton where labels are replaced according to ``d``

        INPUT:

        - ``d`` -- dictionary to determine projection
        - ``det``  --  (default: ``true``) determinize the result or not
        - ``simplify`` -- (default: ``True``) - if True and if det=True, prune and minimize
        - ``verb`` -- boolean (default: ``False``) activate or desactivate the verbose mode

        OUTPUT:

        Return a new projected :class:`NFastAutomaton` (``det``=``False``)
        or  :class:`FastAutomaton` (``det``=``True``)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'a', 'b': 'c', 'c':'c', 'd':'b'}
            sage: a.proj(d)
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef NAutomaton a
        cdef Dict dC
        r = NFastAutomaton(None)
        A2 = []
        sig_on()
        d1 = imagDict(d, self.A, A2=A2)
        sig_off()
        if verb:
            print("d1=%s, A2=%s" % (d1, A2))
        sig_on()
        dC = getDict(d, self.A, d1=d1)
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

    def proji (self, int i, det=True, simplify=True, verb=False):
        """
        Assuming that the alphabet of the automaton are tuples, project on the ith coordinate.
        Give a new automaton where labels are replaced by the projection on the ith coordinate.

        INPUT:

        - ``i`` -- int index of the label projection
        - ``det``  --  (default: ``true``) determinize or not the result
        - ``simplify`` -- (default: ``True``) - if True and if det=True, prune and minimize
        - ``verb`` -- boolean (default: ``False``) to activate or desactivate the verbose mode

        OUTPUT:

        Return a new projected :class:`NFastAutomaton` (``det``=``False``)
        or  :class:`FastAutomaton` (``det``=``True``)


        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.proji(0)
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        d = {}
        for l in self.A:
            if i < len(l):
                d[l] = l[i]
            else:
                raise ValueError("index i %d must be smaller then label size  %d" % (i, len(l)))
        return self.proj(d, det=det, simplify=simplify, verb=verb)

    def determinize_proj(self, d, noempty=True,
                         onlyfinals=False, nof=False, verb=False):
        """
        Project following the dictonary ``d`` and determinize the result.

        INPUT:

        - ``d`` -- dictionary to determine projection
        - ``noempty``  -- (default: ``True``) If ``True`` won't put the empty set state,
            if ``False`` put it and gives a complete automaton.
        - ``onlyfinals``  -- (default: ``False``) 
        - ``nof``  -- (default: ``False``)
        - ``verb`` -- boolean (default: ``False``) activate or desactivate the verbose mode

        OUTPUT:

        Return a new projected  :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'a', 'b': 'c', 'c':'c', 'd':'b'}
            sage: a.determinize_proj(d)
            FastAutomaton with 2 states and an alphabet of 2 letters
        """
        cdef Automaton a
        cdef Dict dC
        if noempty and not onlyfinals and not nof:
            return self.proj(d=d, verb=verb)
        else:
            r = FastAutomaton(None)
            A2 = []
            sig_on()
            d1 = imagDict(d, self.A, A2=A2)
            if verb:
                print("d1=%s, A2=%s" % (d1, A2))
            dC = getDict(d, self.A, d1=d1)
            a = Determinize(self.a[0], dC, noempty, onlyfinals, nof, verb)
            FreeDict(&dC)
            sig_off()
            # FreeAutomaton(self.a[0])
            r.a[0] = a
            r.A = A2
            return r

    # change lettes according to d, duplicating edges if necessary
    # the result is assumed deterministic !!!
    def duplicate(self, d, verb=False):
        """
        Change letters according to dictionnary ``d``, with duplication of edge
        if necessary, the result is assumed deterministic !!!

        INPUT:

        - ``d``  -- dictionary for relabel
        - ``verb`` -- boolean (default: ``False``) fix to ``True`` for activation
          the verbose mode

        OUTPUT:

        Return a new :class:`FastAutomaton` with new letters

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'a', 'b': 'c', 'c':'c', 'd':'b'}
            sage: b = a.duplicate(d)
            sage: b.Alphabet
            ['a', 'c']
        """
        cdef Automaton a
        cdef InvertDict dC
        r = FastAutomaton(None)
        A2 = []
        sig_on()
        d1 = imagDict2(d, self.A, A2=A2)
        if verb:
            print("d1=%s, A2=%s" % (d1, A2))
        dC = getDict2(d, self.A, d1=d1)
        if verb:
            printInvertDict(dC)
        a = Duplicate(self.a[0], dC, len(A2), verb)
        if verb:
            print("end...")
        FreeInvertDict(dC)
        sig_off()
        r.a[0] = a
        r.A = A2
        return r

    def relabel(self, d):
        """
        Change letters of the :class:`FastAutomaton` ON PLACE,
        the dictionary is assumed to be one-to-one.

        INPUT:

         - ``d``  -- dictionary for relabel

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: d = { 'a' : 'a', 'b': 'c', 'c':'b', 'd':'b'}
            sage: a.relabel(d)
            sage: a.Alphabet
            ['a', 'c']
        """
        self.A = [d[c] for c in self.A]

    def permut(self, list A, verb=False):
        """
        Permutes (and eventually remove) letters and return permuted new :class:`FastAutomaton`

        INPUT:

         - ``A``  -- list of letters in the new order (number can be less to the Alphabet)
         - ``verb`` -- boolean (default: ``False``) fix to ``True`` to activate the verbose mode

        OUTPUT:

        Return permuted new :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: l = [ 'b', 'c', 'a']
            sage: b = a.permut(l, verb=True)
            A=['b', 'c', 'a']
            l=[ 1 -1 0 ]
            l = [ 1 -1 0 ]
            sage: b.Alphabet
            ['b', 'c', 'a']

        """
        if verb:
            print("A=%s" % A)
        cdef Automaton a
        r = FastAutomaton(None)
        cdef int *l = <int*>malloc(sizeof(int) * len(A))
        cdef int i
        for i in range(self.a.na):
            l[i] = -1
        d = {}
        for i, c in enumerate(self.A):
            d[c] = i
        for i, c in enumerate(A):
            if d.has_key(c):
                l[i] = d[c]  # l gives the old index from the new one
        if verb:
            str = "l=["
            for i in range(len(A)):
                str += " %s" % l[i]
            str += " ]"
            print(str)
        sig_on()
        a = Permut(self.a[0], l, len(A), verb)
        sig_off()
        free(l)
        r.a[0] = a
        r.A = A

        return r

    # permute les lettres SUR PLACE
    # A = liste des lettres dans le nouvel ordre (il peut y en avoir moins)
    def permut_op(self, list A, verb=False):
        """
        Permutes (and eventually remove) letters ON PLACE.

        INPUT:

         - ``A``  -- list of letters in the new order (number can be less to the Alphabet)
         - ``verb`` -- boolean (default: ``False``) fix to ``True`` for activation
           the verbose mode

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.Alphabet
            ['a', 'b']
            sage: l = [ 'b', 'c', 'a']
            sage: a.permut_op(l, verb=True)
            A=['b', 'c', 'a']
            l=[ 1 -1 0 ]
            l = [ 1 -1 0 ]
            sage: a.Alphabet
            ['b', 'c', 'a']

        """
        if verb:
            print("A=%s" % A)
        cdef int *l = <int*>malloc(sizeof(int) * len(A))
        cdef int i
        for i in range(self.a.na):
            l[i] = -1
        d = {}
        for i, c in enumerate(self.A):
            d[c] = i
        for i, c in enumerate(A):
            if d.has_key(c):
                l[i] = d[c]  # l gives the old index from the new one
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
        self.A = A

    # Compute the transposition, assuming it is deterministic
    def transpose_det(self):
        """
        Transposes :class:`FastAutomaton`, assuming the result to be deterministic

        OUTPUT:

        Return the transposed :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.transpose_det()
            FastAutomaton with 4 states and an alphabet of 2 letters

        """
        r = FastAutomaton(None)
        sig_on()
        r.a[0] = TransposeDet(self.a[0])
        r.A = self.A
        sig_off()
        return r

    def transpose(self):
        """
        Transpose :class:`FastAutomaton` and return a ``NFastAutomaton``

        OUTPUT:

        Return a :class:`NFastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.transpose()
            NFastAutomaton with 4 states and an alphabet of 2 letters
        """
        r = NFastAutomaton(None)
        sig_on()
        r.a[0] = Transpose(self.a[0])
        r.A = self.A
        sig_off()
        return r

    def strongly_connected_components(self, no_trivials=False):
        r"""
        Determine a partition into strongly connected components.
        A strongly connected component is a minimal subset of the set of states such that
        there is no path going outside of the subset, from a state of the subset to a state of the subset.

        INPUT:

        - ``no_trivials`` -- (default: ``False``) If True, do not take into account components
            without any transition from itself to itself (such component contains only one element).

        OUTPUT:

        Return the list of strongly connected components, given as a list of list of states.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.strongly_connected_components()
            [[1], [0], [3], [2]]
        """
        cdef int* l = <int*>malloc(sizeof(int) * self.a.n)
        sig_on()
        cdef int ncc = StronglyConnectedComponents(self.a[0], l)
        sig_off()
        # inverse la liste
        l2 = {}
        cdef int i
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
        free(l)
        return l2.values()

    def acc_and_coacc(self):
        """
        Determine accessible and coaccessible states

        OUTPUT:

        Return the list of accessible and coaccessible states

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.acc_and_coacc()
            [0, 1]
        """
        sig_on()
        cdef int* l = <int*>malloc(sizeof(int) * self.a.n)
        AccCoAcc(self.a, l)
        sig_off()
        return [i for i in range(self.a.n) if l[i] == 1]

    def coaccessible_states(self):
        """
        Compute the co-accessible states of the :class:`FastAutomaton`

        OUTPUT:

        Return list of  co-accessible states of the :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.coaccessible_states()
            [0, 1, 2, 3]
        """
        sig_on()
        cdef int* l = <int*>malloc(sizeof(int) * self.a.n)
        CoAcc(self.a, l)
        sig_off()
        return [i for i in range(self.a.n) if l[i] == 1]

    def sub_automaton(self, l, keep_vlabels=True, verb=False):
        """
        Compute the sub automaton whose states are given by the list ``l``.

        INPUT:

         - ``l``  -- list of states to keep
         - ``verb`` -- boolean (default: ``False``) fix to ``True`` to activate
          the verbose mode
         - ``keep_vlabels``  -- boolean (default: ``True``)

        OUTPUT:

        Return a  :class:`FastAutomaton` sub automaton

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: l = [0, 1]
            sage: a.sub_automaton(l)
            FastAutomaton with 2 states and an alphabet of 2 letters
            sage: a
            FastAutomaton with 4 states and an alphabet of 2 letters
        """
        sig_on()
        r = FastAutomaton(None)
        r.a[0] = SubAutomaton(self.a[0], list_to_Dict(l), verb)
        r.A = self.A
        if keep_vlabels and self.S is not None:
            r.S = [self.S[i] for i in l]
        sig_off()
        return r

    def minimize (self, verb=False):
        """
        Compute the minimal automaton
        by Hopcroft's algorithm
        see [Hopcroft]

        INPUT:

         - ``verb`` -- boolean (default: ``False``) set to ``True`` to activate the verbose mode


        OUTPUT:

        Return a minimized :class:`FastAutomaton` automaton

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
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
            removes the hole state  1...
            FastAutomaton with 3 states and an alphabet of 2 letters
        """
        sig_on()
        r = FastAutomaton(None)
        r.a[0] = Minimise(self.a[0], verb)
        r.A = self.A
        sig_off()
        return r

    def adjacency_matrix(self, sparse=None):
        """
        Compute the adjacency matrix of the :class:`FastAutomaton`

        INPUT:

        - ``sparse`` -- indicate if the return matrix is sparse or not

        OUTPUT:

        Return the corresponding adjacency matrix

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.adjacency_matrix()
            [0 1 0 0]
            [0 0 0 0]
            [0 0 0 1]
            [0 0 0 0]
        """
        if sparse is None:
            if self.a.n <= 128:
                sparse = False
            else:
                sparse = True

        d = {}
        cdef int i, j, f
        for i in range(self.a.n):
            for j in range(self.a.na):
                f = self.a.e[i].f[j]
                if f != -1:
                    if d.has_key((i,f)):
                        d[(i, f)] += 1
                    else:
                        d[(i, f)] = 1
        from sage.matrix.constructor import matrix
        from sage.rings.integer_ring import IntegerRing
        return matrix(IntegerRing(), self.a.n, self.a.n, d, sparse=sparse)

    def delete_vertex(self, int i):
        """
        Gives a copy of the :class:`FastAutomaton` but without the vertex ``i``.

        INPUT:

        - ``i``  - int - the vertex to remove

        OUTPUT:

        Return a copy of :class:`FastAutomaton` with deleted vertex

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.delete_vertex(2)
            FastAutomaton with 3 states and an alphabet of 2 letters

        """
        sig_on()
        r = FastAutomaton(None)
        r.a[0] = DeleteVertex(self.a[0], i)
        r.A = self.A
        sig_off()
        return r

    def delete_vertex_op(self, int i):
        """
        Delete vertex ``i`` on place.

        INPUT:

        - ``i``  -- int number of vertex to remove

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.delete_vertex_op(2)
            sage: a
            FastAutomaton with 3 states and an alphabet of 2 letters

        """
        sig_on()
        DeleteVertexOP(self.a, i)
        sig_off()

    def spectral_radius(self, only_non_trivial=False, verb=False):
        """
        Return the spectral radius of the underlying graph.

        INPUT:

        - ``only_non_trivial``  - (default: ``False``) - if True, don't take into account strongly connected components of cardinality one.
        - ``verb`` - boolean (default: ``False``) - fix to ``True`` to activate
          the verbose mode

        OUTPUT:

        Return spectral radius of strongly connex component

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.spectral_radius(only_non_trivial=False, verb=True)
            minimal Automata : FastAutomaton with 3 states and an alphabet of 2 letters
            3 component strongly connex.
            component with 1 states...
            x
            (x, 1)
            component with 1 states...
            x
            (x, 1)
            component with 1 states...
            x
            (x, 1)
            0

        """
        a = self.minimize()
        if verb:
            print("minimal Automata : %s" % a)
        l = a.strongly_connected_components()
        if verb:
            print("%s component strongly connex." % len(l))
        r = 0 # valeur propre maximale trouvée
        for c in l:
            if not only_non_trivial or len(c) > 1:
                if verb:
                    print("component with %s states..." % len(c))
                b = a.sub_automaton(c)
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
        return r

    def copy(self):
        """
        Do a copy of the :class:`FastAutomaton`.

        OUTPUT:

        Return a copy of the :class:`FastAutomaton`

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.copy()
            FastAutomaton with 4 states and an alphabet of 2 letters

        """
        r = FastAutomaton(None)
        sig_on()
        r.a[0] = CopyAutomaton(self.a[0], self.a.n, self.a.na)
        sig_off()
        r.A = self.A
        return r

    def has_empty_langage(self):
        r"""
        Test if the  :class:`FastAutomaton` has an empty language.

        OUTPUT:

        Return ``True`` the :class:`FastAutomaton` has a empty language
        ``False`` if not

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.has_empty_langage()
            False

        """
        sig_on()
        res = emptyLanguage(self.a[0])
        sig_off()
        return Bool(res)

    def equals_langages(self, FastAutomaton a2, minimized=False, pruned=False, verb=False):
        """
        Test if the languages of :class:`FastAutomaton` ``self`` and ``a2`` are equal or not.

        INPUT:

        - ``a2``  -- the :class:`FastAutomaton` to compare
        - ``minimized``  -- (default: ``False``) if minimization is
          required or not
        - ``pruned``  -- (default: ``False``) if emondation is required or not
        - ``verb`` -- boolean (default: ``False``) fix to ``True`` to activate
          the verbose mode

        OUTPUT:

        Return ``True`` if the both :class:`FastAutomaton` have
        the same language ``False`` if not

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=3)
            sage: c = FastAutomaton([(3, 2, 'd'), (1, 2, 'c')], i=2)
            sage: a.equals_langages(b)
            True
            sage: a.equals_langages(c)
            False
        """
        sig_on()
        cdef Dict d = NewDict(self.a.na)
        cdef int i, j
        for i in range(self.a.na):
            for j in range(a2.a.na):
                if self.A[i] == a2.A[j]:
                    d.e[i] = j
                    if verb:
                        print("%d -> %d"%(i, j))
                    break
        if verb:
            printDict(d)
        res = equalsLanguages(self.a, a2.a, d, minimized, pruned, verb)
        sig_off()
        return Bool(res)

#    def empty_product (self, FastAutomaton a2, d=None, verb=False):
#        if d is None:
#            return self.has_empty_langage() or a2.has_empty_langage()
#        sig_on()
#        cdef Dict dC
#        Av = []
#        dv = imagProductDict(d, self.A, a2.A, Av=Av)
#        if verb:
#            print("Av=%s"%Av
#            print("dv=%s"%dv
#        dC = getProductDict(d, self.A, a2.A, dv=dv, verb=verb)
#        if verb:
#            print("dC="
#            printDict(dC)
#        res = EmptyProduct(self.a[0], a2.a[0], dC, verb)
#        sig_off()
#        return Bool(res)

    def intersect(self, FastAutomaton a2, bool verb=False):
        """
        Determine if the languages of the :class:`FastAutomaton` ``self`` and ``a2`` have a non-empty intersection.

        INPUT:

        -  ``a2``  -- the :class:`Fastautomaton` to intersect
        - ``verb`` -- boolean (default: ``False``) True to active
          the verbose mode

        OUTPUT:

        Return ``True`` if the intersection of the languages is non-empty, return ``False`` otherwise.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: b = FastAutomaton([(3, 2, 'a'), (1, 2, 'd')], i=2)
            sage: a.intersect(b)
            True
        """
        sig_on()
        res = Intersect(self.a[0], a2.a[0], verb)
        sig_off()
        return Bool(res)

#    def intersect (self, FastAutomaton a2, pruned=False, verb=False):
#        sig_on()
#        cdef Dict d = NewDict(self.a.na)
#        cdef int i,j
#        for i in range(self.a.na):
#            for j in range(a2.a.na):
#                if self.A[i] == a2.A[j]:
#                    d.e[i] = j
#                    if verb: print("%d -> %d"%(i, j))
#                    break
#        if verb:
#            printDict(d)
#        res = intersectLanguage(self.a, a2.a, d, pruned, verb)
#        sig_off()
#        return Bool(res)

    def find_word(self, bool verb=False):
        """
        Find a word in the language of the automaton

        INPUT:

        - ``verb`` -- (default: ``False``)  the verbose parameter

        OUTPUT:

        return a word of the language of the Automaton as list of letters

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.find_word()
            []
        """
        sig_on()
        cdef Dict w
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
        Compute a shortest words of the automaton

        INPUT:

        - ``i`` -- (default: None)  the initial state
        - ``f`` -- (default: None)  the final state
        - ``verb`` -- (default: False)  the verbose parameter

        OUTPUT:

        return a word, as list of letters

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.shortest_word(i=2, f=1)
            []
        """
        cdef Dict w
        if i is None:
            i = self.a.i
        if f is None:
            f = -1
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

    def shortest_words(self, i=None, verb=False):
        """
        Compute the shortest words of the automaton

        INPUT:

        - ``i`` -- (default: None)  the initial state
        - ``verb`` -- (default: False)  the verbose parameter

        OUTPUT:

        return a list of words, as a list of list of letters.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.shortest_words()
            [[], ['a'], [], []]

        """
        cdef Dict* w = <Dict*>malloc(sizeof(Dict) * self.a.n)
        if i is None:
            i = self.a.i
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
        free(w)
        return rt

    # determine if the word is recognized by the automaton or not
    def rec_word2(self, list w):
        """
        Determine if the word ``w`` is recognized or nor not by the automaton

        INPUT:

        - ``w`` -- a list of letters

        OUTPUT:

        return 1 if the word is recognized, otherwise 0

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: w = ['a', 'b', 'b']
            sage: a.rec_word2(w)
            0

        """
        rd = {}
        for i, a in enumerate(self.A):
            rd[a] = i
        sig_on()
        cdef Dict d = NewDict(len(w))
        cdef bool res
        for i, a in enumerate(w):
            d.e[i] = rd[a]
        res = rec_word(self.a[0], d)
        sig_off()
        return res

    # determine if the word is recognized by the automaton or not
    def rec_word(self, list w):
        """
        Determine if the word ``w`` is recognized or nor not by the automaton

        INPUT:

        - ``w`` -- a list of letters

        OUTPUT:

        return ``True`` if the word is recognized, otherwise ``False``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: w = ['a', 'b', 'b']
            sage: a.rec_word(w)
            False


        """
        cdef int e = self.a.i
        if e == -1:
            return False
        d = {}
        for i, a in enumerate(self.A):
            d[a] = i
        for a in w:
            e = self.a.e[e].f[d[a]]
            if e == -1:
                return False
        return Bool(self.a.e[e].final)

    def add_state(self, bool final):
        """
        Add a state in the automaton

        INPUT:

        - ``final`` -- boolean indicate if the added state is final

        OUTPUT:

        return the new state (which is an integer)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.add_state(True)
            4
            sage: a.add_state(False)
            5
        """
        sig_on()
        AddEtat(self.a, final)
        sig_off()
        return self.a.n-1

    def add_edge(self, int i, l, int j):
        """
        Add an edge in the automaton

        INPUT:

        - ``i`` - int - the first state
        - ``l`` -- the label of edge
        - ``j`` - int - the second state

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.add_edge(2,'a',1)

            sage: a.add_edge(2,'v',1)
            Traceback (most recent call last):
            ...
            ValueError: The letter v doesn't exist.
            sage: a.add_edge(2,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state  6 doesn't exist.
            sage: a.add_edge(5,'v',6)
            Traceback (most recent call last):
            ...
            ValueError: The state  5 doesn't exist.
        """
        if i >= self.a.n:
            raise ValueError("The state %s doesn't exist." % i)
        if j >= self.a.n:
            raise ValueError("The state  %s doesn't exist." %j)
        try:
            k = self.A.index(l)
        except:
            # La lettre %s n'existe pas.
            raise ValueError("The letter %s doesn't exist." % l)
        self.a.e[i].f[k] = j

    @property
    def n_states(self):
        """
        return the numbers of states

        OUTPUT:

        return the numbers of states

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.n_states
            4
        """
        return self.a.n

    def bigger_alphabet(self, nA):
        """
        Gives a copy of the :class:`FastAutomaton`, but with the bigger alphabet ``nA``

        INPUT:

        - ``nA`` --  Alphabet of the new automaton. We assume that it contains the current alphabet of the automaton.

        OUTPUT:

        return a :class:`FastAutomaton` with a bigger alphabet ``nA``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.bigger_alphabet(['a','b','c'])
            FastAutomaton with 4 states and an alphabet of 3 letters
        """
        cdef Dict d
        r = FastAutomaton(None)
        sig_on()
        d = NewDict(self.a.na)
        for i in range(self.a.na):
            if self.A[i] in nA:
                d.e[i] = nA.index(self.A[i])
        r.a[0] = BiggerAlphabet(self.a[0], d, len(nA))
        sig_off()
        r.A = nA
        return r

    def complementary_op(self):
        """
        Change the language of the automaton to the complementary ON PLACE.

        OUTPUT:

        return None
        (the operation is on place)

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complementary_op()
            sage: a
            FastAutomaton with 5 states and an alphabet of 2 letters

        """
        self.complete()
        cdef i
        for i in range(self.a.n):
            self.a.e[i].final = not self.a.e[i].final

    def complementary(self):
        """
        Gives an automaton whose language is the complementary.

        OUTPUT:

        return  a new automaton whose language is the complementary of the language of ``self``

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.complementary()
            FastAutomaton with 5 states and an alphabet of 2 letters
        """
        a = self.copy()
        a.complementary_op()
        return a

    def included(self, FastAutomaton a, bool verb=False, pruned=False):
        r"""
        Test if the language of self is included in the language of ``a``

        INPUT:

        - ``a`` --  a :class:`FastAutomaton`
        - ``verb`` -- (default: False) verbose parameter
        - ``pruned`` -- (default: False) set to True if the automaton is already pruned, in order to avoid unuseful computation.

        OUTPUT:

        return  ``True`` if the language of ``self`` is included in the language of ``a``
        return    ``False`` otherwise

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.included(a)
            True
            sage: b = FastAutomaton([(0, 1, 'c')], i=0)
            sage: a.included(b)
            False
            sage: b = FastAutomaton([(0, 1, 'a')], i=0)
            sage: b.included(a)
            True
        """
        cdef FastAutomaton b, a2
        if self.A != a.A:
            A = list(set(a.A+self.A))
            b = self.bigger_alphabet(A)
            a2 = a.bigger_alphabet(A)
        else:
            b = self
            a2 = a
        sig_on()
        res = Included(b.a[0], a2.a[0], pruned, verb)
        sig_off()
        return Bool(res)

#        d = {}
#        for l in self.A:
#            if l in a.A:
#                d[(l,l)] = l
#        if verb:
#            print("d=%s"%d)
#        a.complete()
#        cdef FastAutomaton p = self.product(a, d, verb=verb)
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
#        return p.has_empty_langage()

#    # donne un automate reconnaissant w(w^(-1)L) où L est le langage
#    # de a partant de e
#    def piece(self, w, e=None):
#        """
#        return a automaton recognizing ``w`` as :math:`w (w^{-1}L)` where ``L``
#        is the language of automaton a from e entry state
#
#        INPUT:
#
#        - ``w`` --  word
#        - ``e`` -- (default: None) the entry state
#
#        OUTPUT:
#
#        return  a automaton recognizing ``w``
#
#        EXAMPLES::
#
#            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
#            sage: a.piece([1])
#            FastAutomaton with 0 states and an alphabet of 2 letters
#
#        """
#        cdef int* l = <int*>malloc(sizeof(int)*self.a.n)
#        cdef int i
#        if type(w) != list:
#            w = [int(w)]
#        for i in range(len(w)):
#            l[i] = w[i]
#        if e is None:
#            e = self.a.i
#        r = FastAutomaton(None)
#        sig_on()
#        r.a[0] = PieceAutomaton(self.a[0], l, len(w), e)
#        sig_off()
#        free(l)
#        r.A = self.A
#        return r

    # tell if the language of the automaton is empty
    # (this function is not very efficient)
    def is_empty(self):
        """
        Examines if the language of the automaton is empty

        OUTPUT:

        return ``True`` if the language of the automaton is empty, ``False`` otherwise

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b')], i=0)
            sage: a.is_empty()
            False

        """
        return (self.find_word() is None)

    def random_word(self, nmin=-1, nmax=100):
        r"""
        Return a random word recognized by the automaton, by following a random path in the automaton from the initial state. If we don't fall into the set of final states before reaching the maximal length ``nmax``, then return ``word not found !``.

        INPUT:

        - ``nmin`` -- (default: -1) minimal length of the word
        - ``nmax`` -- (default: 100) maximal length of the word

        OUTPUT:

        Return a random word of length between ``nmin``and ``nmax`` if found.
        Otherwise return ``word not found !``.

        EXAMPLES::

            sage: a = FastAutomaton([(0, 1, 'a'), (2, 3, 'b'), (0, 3, 'c')], i=0)
            sage: a.random_word() # random
            ['a']


        """
        cdef int i = self.a.i
        w = []
        na = len(self.A)
        if nmin < 0:
            nmin = 1
        from sage.misc.prandom import random
        for j in range(nmin):
            li = [l for l in range(na) if self.succ(i, l) != -1]
            l = li[(int)(random() * len(li))]
            w.append(self.A[l])
            i = self.succ(i, l)
        # continue the word to get into a final state
        for j in range(nmax-nmin):
            if self.a.e[i].final:
                break
            li = [l for l in range(na) if self.succ(i, l) != -1]
            l = li[(int)(random() * len(li))]
            w.append(self.A[l])
            i = self.succ(i, l)
        if not self.a.e[i].final:
            print("word not found !")  # "Mot non trouvé !"
        return w
