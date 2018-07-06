# coding=utf8
"""
Wrapper for automatonC  of Finite state machines using C


AUTHORS:

- Paul Mercat (2013)- I2M AMU Aix-Marseille Universite - initial version
- Dominique Benielli (2018) Labex Archimede - I2M -
  AMU Aix-Marseille Universite - Integration in -SageMath


"""

#*****************************************************************************
#       Copyright (C) 2014 Paul Mercat <mercatp@icloud.com>
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

cimport sage.combinat.words.automatonC_wrapper
from cysignals.signals cimport sig_on, sig_off, sig_check
from cpython cimport bool as c_bool
from .cautomata import  FastAutomaton, NFastAutomaton
from .cautomata cimport getAutomaton, list_to_Dict

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
    void FreeNAutomaton(NAutomaton *a)
    Automaton CopyAutomaton(Automaton a, int nalloc, int naalloc)
    Automaton PieceAutomaton(Automaton a, int *w, int n, int e)
    void init(Automaton *a)
    void printAutomaton(Automaton a)
    void plotDot(const char *file, Automaton a, const char **labels, const char *graph_name, double sx, double sy, const char **vlabels, bool html, bool verb, bool run_dot)
    void NplotDot(const char *file, NAutomaton a, const char **labels, const char *graph_name, double sx, double sy, bool run_dot)
    Automaton Product(Automaton a1, Automaton a2, Dict d, bool verb)
    Automaton Determinise(Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
    Automaton DeterminiseN(NAutomaton a, bool puits, int verb)
    NAutomaton Concat(Automaton a, Automaton b, bool verb)
    NAutomaton CopyN(Automaton a, bool verb)
    void AddEdgeN(NAutomaton *a, int e, int f, int l)
    void AddPathN(NAutomaton *a, int e, int f, int *l, int len, bool verb)
    NAutomaton Proj(Automaton a, Dict d, bool verb)
    void ZeroComplete(Automaton *a, int l0, bool verb)
    Automaton ZeroComplete2(Automaton *a, int l0, bool etat_puits, bool verb)
    Automaton ZeroInv(Automaton *a, int l0)
    Automaton emonde_inf(Automaton a, bool verb)
    Automaton emonde(Automaton a, bool verb)
    Automaton emondeI(Automaton a, bool verb)
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
    bool equalsLangages(Automaton *a1, Automaton *a2, Dict a1toa2, bool minimized, bool emonded, bool verb)
    bool Intersect(Automaton a1, Automaton a2, bool verb)
    bool Included(Automaton a1, Automaton a2, bool emonded, bool verb)
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
    # bool intersectLangage (Automaton *a1, Automaton *a2, Dict a1toa2, bool emonded, bool verb)

def _printAutomaton_wrapper(a):
    """
    Test the printAutomaton c Function.

    INPUT:

    - ``a`` -- DiGraph

    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _printAutomaton_wrapper
        sage: li = ([(0, 1, 'a'), (2, 3, 'b')]
        sage: di = DiGraph(li, multiedges=True, loops=True)
        sage: _printAutomaton_wrapper()
        Automaton with 4 states, 2 letters.
        0 --0--> 1
        2 --1--> 3
        initial state -1.


    """
    cdef Automaton  r
    A = list(set(a.edge_labels()))
    F = a.vertices()
    r = getAutomaton(a, initial=None, F=F, A=A)
    sig_on()
    printAutomaton(r)
    sig_off()


def _dotExists_wrapper():
    """
    Test the Doexist c Function.

    OUTPUT:

    Return ``True`` or `False `` if file exist

    TESTS::

        sage: _dotExists_wrapper()
        sage: True


    """
    sig_on()
    if DotExists():
        sig_off()
        return True
    sig_off()
    return False


def _copyAutomaton_wrapper(fa):
    """
    Test the copyAutomaton c Function.


    INPUT:

    - ``fa`` -- DiGraph


    OUTPUT:

    Return a tuple with a dictionnary and the interger initial state filled from 
    a copy of c automaton


    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _copyAutomaton_wrapper
        sage: li = ([(0, 1, 'a'), (2, 3, 'b')]
        sage: di = DiGraph(li, multiedges=True, loops=True)
        sage: _copyAutomaton_wrapper(a)
        ({0: {0: 1, 1: -1, 2: 1},
          1: {0: -1, 1: -1, 2: 1},
          2: {0: -1, 1: 3, 2: 1},
          3: {0: -1, 1: -1, 2: 1}},
         -1)

    """
    cdef Automaton aut, r
    A = list(set(fa.edge_labels()))
    F = fa.vertices()
    a = getAutomaton(fa, initial=None, F=F, A=A)
    sig_on()
    aut = CopyAutomaton(a, a.n, a.na)
    sig_off()
    e = {}
    for i in range(aut.n):
        pf = {}
        for j in range(aut.na):
            pf[j] = aut.e[i].f[j]
        pf[j+1] = aut.e[i].final
        e[i] = pf
    return e, aut.i


def _NewDict(n):
    """
    Test the _NewDict c Function.


    INPUT:

    - ``n`` -- number of elements of the new dict


    OUTPUT:

    Return a new dict of c automaton

    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _NewDict
        sage: r = _NewDict(3)
        sage: r
        {0: -1, 1: -1, 2: -1}

    """
    cdef Dict r
    d = _NewDict(n)
    sig_on()
    r = NewDict(n)
    sig_off()
    for i in range(r.n):
        d[i] = r.e[i]
    sig_on()
    FreeDict(&r)
    sig_off()
    return d

def _NewInvertDict(n):
    """
    Test the _NewInvertDict c Function.


    INPUT:

    - ``n`` -- number of elements of the new dict


    OUTPUT:

    Return a dict of new InvertDict of c automaton

    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _NewInvertDict
        sage: r = _NewInvertDict(3)
        sage: r
        {0: {}, 1: {}, 2: {}}

    """
    cdef InvertDict r
    d = {}
    sig_on()
    r = NewInvertDict(n)
    sig_off()
    for i in range(r.n):
        df = {}
        for j in range(r.d[i].n):
            df[j] = r.d[i].e[j]
        d[i] = df
    sig_on()
    FreeInvertDict(r)
    sig_off()
    return d


# def _InvertDict(li):
#     """
#     Test the _NewInvertDict c Function.
# 
# 
#     INPUT:
# 
#     - ``n`` -- number of elements of the new dict
# 
# 
#     OUTPUT:
# 
#     Return a dict of new InvertDict of c automaton
# 
#     TESTS::
# 
#         sage: from  sage.combinat.words.automatonC_wrapper import _NewInvertDict
#         sage: r = _InvertDict([0, 1, 2, 3])
#         sage: r
#         {0: {}, 1: {}, 2: {}}
# 
#     """
#     cdef InvertDict r
#     cdef Dict d
#     dr = {}
#     d = list_to_Dict(li)
#     sig_on()
#     r = invertDict(d)
#     sig_off()
#     for i in range(r.n):
#         df = {}
#         for j in range(r.d[i].n):
#             df[j] = r.d[i].e[j]
#         dr[i] = df
#     sig_on()
#     FreeInvertDict(r)
#     sig_off()
#     return dr

def _RecWord(fa, li):
    """
    Test the _NewDict c Function.


    INPUT:

    - ``a`` -- DiGraph

    - `li` -- list

    OUTPUT:

    Return ``True`` or `False `` if 


    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _RecWord
        sage: a = [(0, 1, 'a'), (2, 3, 'b')]
        sage: fa = DiGraph(a, multiedges=True, loops=True)
        sage: li = [0, 1, 2]
        sage: b = _RecWord(fa, li)
        sage: b
        False

    """
    cdef bool res
    cdef Automaton a
    cdef Dict d
    A = list(set(fa.edge_labels()))
    F = fa.vertices()
    a = getAutomaton(fa, initial=None, F=F, A=A)
    d = list_to_Dict(li)
    sig_on()
    res = rec_word(a, d)
    sig_off()
    return c_bool(res)

def _emptyLanguage_wrapper(a):
    """
    Test the _IsCompleteAutomaton c Function.


    INPUT:

    - ``a`` -- DiGraph

    OUTPUT:

    Return ``True`` or `False `` if automaton of ``DiGraph``
    has a empty language.


    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _emptyLanguage_wrapper
        sage: a = [(0, 1, 'a'), (2, 3, 'b')]
        sage: di = DiGraph(a, multiedges=True, loops=True)
        sage: _emptyLanguage_wrapper(di)
        True
    """
    cdef bool res
    cdef Automaton r
    A = list(set(a.edge_labels()))
    F = a.vertices()
    r = getAutomaton(a, initial=None, F=F, A=A)
    sig_on()
    res = emptyLanguage(r)
    sig_off()
    return c_bool(res)


def _IsCompleteAutomaton_wrapper(a):
    """
    Test the _IsCompleteAutomaton c Function.


    INPUT:

    - ``a`` -- DiGraph

    OUTPUT:

    Return ``True`` or `False `` if automaton of ``DiGraph``
    is complete or not


    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _IsCompleteAutomaton_wrapper
        sage: a = [(0, 1, 'a'), (2, 3, 'b')]
        sage: di = DiGraph(a, multiedges=True, loops=True)
        sage: _IsCompleteAutomaton_wrapper(di)
        False
    """
    cdef bool res
    cdef Automaton r
    A = list(set(a.edge_labels()))
    F = a.vertices()
    r = getAutomaton(a, initial=None, F=F, A=A)
    sig_on()
    res = IsCompleteAutomaton(r)
    sig_off()
    return c_bool(res)


def _CompleteAutomaton_wrapper(a):
    """
    Test the _IsCompleteAutomaton c Function.


    INPUT:

    - ``a`` -- DiGraph

    OUTPUT:

    Return ``True`` or `False `` if automaton of ``DiGraph``
    is complete or not


    TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _CompleteAutomaton_wrapper
        sage: a = [(0, 1, 'a'), (2, 3, 'b')]
        sage: di = DiGraph(a, multiedges=True, loops=True)
        sage: _CompleteAutomaton_wrapper(di)
        True
    """
    cdef bool res
    cdef Automaton r
    A = list(set(a.edge_labels()))
    F = a.vertices()
    r = getAutomaton(a, initial=None, F=F, A=A)
    sig_on()
    res = CompleteAutomaton(&r)
    sig_off()
    return c_bool(res)

def _equalsAutomaton_wrapper(a1, a2):
    """
    Test the _IsCompleteAutomaton c Function.


    INPUT:

    - ``a1`` -- DiGraph
    - ``a2`` -- DiGraph

    OUTPUT:

    Return ``True`` or `False `` if automaton of input ``DiGraph``
    are equal or not
        TESTS::

        sage: from  sage.combinat.words.automatonC_wrapper import _equalsAutomaton_wrapper
        sage: a1 = [(0, 1, 'a'), (2, 3, 'b')]
        sage: di1 = DiGraph(a1, multiedges=True, loops=True)
        sage: _equalsAutomaton_wrapper(di1, di1)
        True
        sage: a2 = [(0, 1, 'a'), (0, 3, 'b')]
        sage: di2 = DiGraph(a2, multiedges=True, loops=True)
        sage: _equalsAutomaton_wrapper(di1, di2)
        False
    """
    cdef bool res
    cdef Automaton r1, r2
    A1 = list(set(a1.edge_labels()))
    F1 = a1.vertices()
    r1 = getAutomaton(a1, initial=None, F=F1, A=A1)
    A2 = list(set(a2.edge_labels()))
    F2 = a2.vertices()
    r2 = getAutomaton(a2, initial=None, F=F2, A=A2)
    sig_on()
    res = equalsAutomaton(r1, r2)
    sig_off()
    return c_bool(res)

#     """
#     Test the _IsCompleteAutomaton c Function.
# 
# 
#     INPUT:
# 
#     - ``a1`` -- DiGraph
#     - ``a2`` -- DiGraph
# 
#     OUTPUT:
# 
#     Return ``True`` or `False `` if automaton of input ``DiGraph``
#     are equal or not
#         TESTS::
# 
#         sage: from  sage.combinat.words.automatonC_wrapper import _equalsAutomaton_wrapper
#         sage: a1 = [(0, 1, 'a'), (2, 3, 'b')]
#         sage: di1 = DiGraph(a1, multiedges=True, loops=True)
#         sage: _equalsAutomaton_wrapper(di1, di1)
#         True
#         sage: a2 = [(0, 1, 'a'), (0, 3, 'b')]
#         sage: di2 = DiGraph(a2, multiedges=True, loops=True)
#         sage: _equalsAutomaton_wrapper(di1, di2)
#         False
#     """
# bool findWord(Automaton a, Dict *w, bool verb)
