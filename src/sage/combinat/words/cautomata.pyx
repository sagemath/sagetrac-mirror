# coding=utf8
"""
Finite state machines using C

AUTHORS:

- Paul Mercat
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
#*****************************************************************************

from libc.stdlib cimport malloc, free

cimport sage.combinat.words.cautomata

#ctypedef Automate Automaton

cdef extern from "automataC.h":
    ctypedef Automate Automaton
    cdef cppclass Dict:
        int* e
        int n
    cdef cppclass InvertDict:
        Dict* d
        int n
    
    Automaton NewAutomaton (int n, int na)
    void FreeAutomaton (Automaton a)
    void init (Automaton a)
    void printAutomaton (Automaton a)
    void plotTikZ (Automaton a, const char **labels, const char *graph_name, double sx, double sy)
    Automaton Product(Automaton a1, Automaton a2, Dict d)
    Automaton Determinise (Automaton a, Dict d, bool noempty, bool onlyfinals, bool nof, bool verb)
    Automaton emonde_inf (Automaton a)
    Automaton emonde (Automaton a, bool verb)
    Automaton emondeI (Automaton a, bool verb)
    bool equalsAutomaton (Automaton a1, Automaton a2)
    Dict NewDict (int n)
    void FreeDict (Dict d)
    void printDict (Dict d)
    InvertDict NewInvertDict (int n)
    void FreeInvertDict (InvertDict id)
    void printInvertDict (InvertDict id)
    Automaton Duplicate (Automaton a, InvertDict id, int na2, bool verb)
    Automaton Transpose (Automaton a)
    int StronglyConnectedComponents (Automaton a, int *res)
    Automaton SubAutomaton (Automaton a, Dict d, bool verb)

#dictionnaire numérotant l'alphabet projeté
cdef imagDict (dict d, list A, list A2=[]):
    d1 = {}
    i = 0
    for a in A:
        if d.has_key(a):
            if not d1.has_key(d[a]):
                d1[d[a]] = i
                A2.append(d[a])
                i+=1
    return d1

#dictionnaire numérotant le nouvel alphabet
cdef imagDict2 (dict d, list A, list A2=[]):
    #print "d=%s, A=%s"%(d,A)
    d1 = {}
    i = 0
    for a in A:
        if d.has_key(a):
            for v in d[a]:
                if not d1.has_key(v):
                    d1[v] = i
                    A2.append(v)
                    i+=1
    return d1

cdef Dict getDict (dict d, list A, dict d1=None):
    A = list(A)
    cdef Dict r
    r = NewDict(len(A))
    cdef int i
    if d1 is None:
        d1 = imagDict(d, A)
    #print d1
    for i in range(r.n):
        if d.has_key(A[i]):
            r.e[i] = d1[d[A[i]]]
        else:
            r.e[i] = -1
    return r
    
cdef Dict list_to_Dict (list l):
    cdef Dict d = NewDict(len(l))
    cdef int i
    for i in range(len(l)):
        d.e[i] = l[i]
    return d

cdef InvertDict getDict2 (dict d, list A, dict d1=None):
    A = list(A)
    cdef InvertDict r
    r = NewInvertDict(len(A))
    cdef int i
    if d1 is None:
        d1 = imagDict2(d, A)
    #print d1
    for i in range(r.n):
        if d.has_key(A[i]):
            r.d[i] = NewDict(len(d[A[i]]))
            for j in range(r.d[i].n):
                r.d[i].e[j] = d1[d[A[i]][j]]
        else:
            r.d[i].n = 0;
    return r

#dictionnaire numérotant l'alphabet projeté
cdef imagProductDict (dict d, list A1, list A2, list Av=[]):
    dv = {}
    i = 0
    for a1 in A1:
        for a2 in A2:
            if d.has_key((a1,a2)):
                if not dv.has_key(d[(a1,a2)]):
                    dv[d[(a1,a2)]] = i
                    Av.append(d[(a1,a2)])
                    i+=1
    return dv

cdef Dict getProductDict (dict d, list A1, list A2, dict dv=None):
    cdef Dict r
    d1 = {}
    d2 = {}
    cdef int i, n1, n2
    n1 = len(A1)
    for i in range(n1):
        d1[A1[i]] = i
    n2 = len(A2)
    for i in range(n2):
        d2[A2[i]] = i
    if dv is None:
        dv = imagProductDict(d, A1, A2)
    r = NewDict(n1*n2)
    Keys = d.keys()
    for (a1,a2) in Keys:
        r.e[d1[a1]+d2[a2]*n1] = dv[d[(a1,a2)]]
    return r

def TestAutomaton (a):
    cdef Automaton r
    #d = {}
    #da = {}
    r = getAutomaton(a) #, d, da)
    printAutomaton(r)
    #print d, da
    print a.vertices(), list(a.Alphabet())

def TestProduct (a1, a2, di):
    cdef Automaton a,b,c
    a = getAutomaton(a1)
    b = getAutomaton(a2)
    printAutomaton(a)
    print a1.vertices(), a1.Alphabet()
    printAutomaton(b)
    print a2.vertices(), a2.Alphabet()
    cdef Dict d
    d = getProductDict(di, list(a1.Alphabet()), list(a2.Alphabet()))
    print "dictionnaire du produit :"
    printDict(d)
    c = Product(a, b, d)
    print "résultat :"
    printAutomaton(c)

#def TestDeterminise (a, d, noempty=True, verb=True):
#    cdef Dict di = getDict(d, a.Alphabet())
#    cdef Automaton au = getAutomaton(a)
#    if verb:
#        printDict(di)
#    if verb:
#        printAutomaton(au)
#    cdef Automaton r = Determinise(au, di, noempty, verb)
#    printAutomaton(r)

#def TestDeterminiseEmonde (a, d, noempty=True, verb=True):
#    cdef Dict di = getDict(d, a.Alphabet())
#    cdef Automaton au = getAutomaton(a)
#    if verb:
#        printDict(di)
#    if verb:
#        printAutomaton(au)
#    cdef Automaton r = Determinise(au, di, noempty, verb)
#    print "Avant émondation :"
#    printAutomaton(r)
#    cdef Automaton r2 = emonde_inf(r)
#    print "Après émondation :"
#    printAutomaton(r2)
#    if equalsAutomaton(r, r2):
#        print "equals !"
#    else:
#        print "differents !"

def TestEmonde (a, noempty=True, verb=True):
    cdef Automaton au = getAutomaton(a)
    if verb:
        print "Avant émondation :"
        printAutomaton(au)
    cdef Automaton r = emonde_inf(au)
    if verb:
        print "Après émondation :"
        printAutomaton(r)
    if equalsAutomaton(r, au):
        print "equals !"
    else:
        print "differents !"
    return AutomatonGet(r)
    

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
        A = list(a.Alphabet())
    V = list(a.vertices())
    cdef int n = len(V)
    cdef int na = len(A)
    r = NewAutomaton(n, na)
    init(r)
    for i in range(na):
        da[A[i]] = i
    for i in range(n):
        r.e[i].final = 0
        d[V[i]] = i
    for v in F:
        r.e[d[v]].final = 1
    if initial is None:
        if not hasattr(a, 'I'):
            I = []
            #raise ValueError("I must be defined !")
        else:
            I = list(a.I)
        if len(I) > 1:
            print "L'automate doit être déterministe ! (I=%s)"%a.I
        if len(I) >= 1:
            r.i = d[I[0]]
        else:
            r.i = -1
    else:
        r.i = initial
    for e,f,l in a.edges():
        r.e[d[e]].f[da[l]] = d[f]
    return r

cdef AutomatonGet (Automaton a, A=None):
    from sage.combinat.words.automata import Automaton
    r = Automaton(multiedges=True, loops=True)
    cdef int i,j
    r.F = []
    if A is None:
        A = [i for i in range(a.na)]
    for i in range(a.n):
        for j in range(a.na):
            if a.e[i].f[j] != -1:
                r.add_edge((i, a.e[i].f[j], A[j]))
        if a.e[i].final:
            r.F.append(i)
    r.I = [a.i]
    return r

#cdef initFA (Automaton *a):
#    *a = NewAutomaton(1,1)

cdef Bool (int x):
    if x:
        return True
    return False

cdef class FastAutomaton:
    
#    cdef Automaton* a
#    cdef list A

    def __cinit__ (self):
        #print "cinit"
        self.a = <Automaton *>malloc(sizeof(Automaton))
        #initialise
        self.a[0].e = NULL
        self.a[0].n = 0
        self.a[0].na = 0
        self.a[0].i = -1
        self.A = []
    
    def __init__(self, a, i=None, F=None, A=None):
        #print "init"
        if a is None:
            return
        from sage.graphs.digraph import DiGraph
        if isinstance(a, list):
            a = DiGraph(a)
        if isinstance(a, DiGraph):
            if A is None:
                if hasattr(a, 'A'):
                    A = list(a.A)
                else:
                    A = list(set(a.edge_labels()))
            self.A = A
            self.a[0] = getAutomaton(a, initial=i, F=F, A=self.A)
        else:
            raise ValueError("Cannot convert the input to FastAutomaton.")
    
    def __dealloc__ (self):
        #print "free"
        FreeAutomaton(self.a[0])
        free(self.a)
    
    def __repr__ (self):
        return "FastAutomaton with %d states and an alphabet of %d letters"%(self.a.n, self.a.na)
    
    def Automaton (self):
        return AutomatonGet(self.a[0], self.A)
    
    def plot (self, int sx=10, int sy=8):
        cdef char** ll
        ll = <char **>malloc(sizeof(char*)*self.a.na)
        cdef int i
        strA = []
        for i in range(self.a.na):
            strA.append(str(self.A[i]))
            ll[i] = strA[i]
        plotTikZ(self.a[0], ll, "Automaton", sx, sy)
        free(ll);
        
        #self.Automaton().plot2()
    
    def Alphabet (self):
        return self.A
    
    def initial_state (self):
        return self.a.i
    
    def final_states (self):
        l = []
        for i in range(self.a.n):
            if self.a.e[i].final:
                l.append(i)
        return l
    
    def emonde_inf (self):
        #sig_on()
        cdef Automaton a
        r = FastAutomaton(None)
        a = emonde_inf(self.a[0])
        #FreeAutomaton(self.a[0])
        r.a[0] = a
        r.A = self.A
        #sig_off()
        return r
    
    def emonde_i (self, verb=False):
        #sig_on()
        cdef Automaton a
        r = FastAutomaton(None)
        a = emondeI(self.a[0], verb)
        r.a[0] = a
        r.A = self.A
        #sig_off()
        return r
    
    def emonde (self, verb=False):
        #sig_on()
        cdef Automaton a
        r = FastAutomaton(None)
        a = emonde(self.a[0], verb)
        r.a[0] = a
        r.A = self.A
        #sig_off()
        return r
    
    def equals (self, FastAutomaton b):
        return Bool(equalsAutomaton(self.a[0], b.a[0]))
    
    #assume that the dictionnary d is injective !!!
    def product (self, FastAutomaton b, dict d, verb=False):
        #sig_on()
        cdef Automaton a
        cdef Dict dC
        r = FastAutomaton(None)
        Av = []
        dv = imagProductDict(d, self.A, b.A, Av=Av)
        if verb:
            print "Av=%s"%Av
            print "dv=%s"%dv
        dC = getProductDict(d, self.A, b.A, dv=dv)
        if verb:
            print "dC="
            printDict(dC)
        a = Product(self.a[0], b.a[0], dC)
        FreeDict(dC)
        r.a[0] = a
        r.A = Av
        #sig_off()
        return r
    
    def intersection (self, FastAutomaton a, verb=False):
        d = {}
        for l in self.A:
            if l in a.A:
                d[(l,l)] = l
        if verb:
            print "d=%s"%d
        return self.product(a, d, verb=verb)
        
    def determinise_proj (self, dict d, noempty=True, onlyfinals=False, nof=False, verb=False):
        #sig_on()
        cdef Automaton a
        cdef Dict dC
        r = FastAutomaton(None)
        A2 = []
        d1 = imagDict(d, self.A, A2=A2)
        if verb:
            print "d1=%s, A2=%s"%(d1,A2)
        dC = getDict(d, self.A, d1=d1)
        a = Determinise (self.a[0], dC, noempty, onlyfinals, nof, verb)
        FreeDict(dC)
        #FreeAutomaton(self.a[0])
        r.a[0] = a
        r.A = A2
        #sig_off()
        return r
    
    #change les lettres selon d, en dupliquant les arêtes si nécessaire
    #the result is assumed deterministic !!!
    def duplicate (self, dict d, verb=False):
        #sig_on()
        cdef Automaton a
        cdef InvertDict dC
        r = FastAutomaton(None)
        A2 = []
        d1 = imagDict2(d, self.A, A2=A2)
        if verb:
            print "d1=%s, A2=%s"%(d1,A2)
        dC = getDict2(d, self.A, d1=d1)
        if verb:
            printInvertDict(dC)
        a = Duplicate (self.a[0], dC, len(A2), verb)
        if verb:
            print "end..."
        FreeInvertDict(dC)
        r.a[0] = a
        r.A = A2
        #sig_off()
        return r
    
    def transpose (self):
        #sig_on()
        r = FastAutomaton(None)
        r.a[0] = Transpose(self.a[0])
        r.A = self.A
        #sig_off()
        return r
    
    def strongly_connected_components (self):
        #sig_on()
        cdef int* l = <int*>malloc(sizeof(int)*self.a.n)
        cdef int ncc = StronglyConnectedComponents(self.a[0], l)
        #inverse la liste
        l2 = {}
        cdef int i
        for i in range(self.a.n):
            if not l2.has_key(l[i]):
                l2[l[i]] = []
            l2[l[i]].append(i)
        free(l)
        #sig_off()
        return l2.values()
    
    def sub_automaton (self, l, verb=False):
        #sig_on()
        r = FastAutomaton(None)
        r.a[0] = SubAutomaton(self.a[0], list_to_Dict(l), verb)
        r.A = self.A
        #sig_off()
        return r
        