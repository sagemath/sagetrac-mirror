# coding=utf8
"""
Finite state machines

AUTHORS:

- Paul Mercat
"""

#*****************************************************************************
#       Copyright (C) 2013 Paul Mercat <mercatp@icloud.com>
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

from sage.graphs.digraph import DiGraph
from sage.sets.set import Set

cdef extern from "automata_tools.c":
    void TikZ (const char *data, const char *graph_name, double sx, double sy)

def SaveTikZ (data, name, sx, sy):
    #print len(data), data
    TikZ(data, name, sx, sy)

#test
#
#cdef emonde0_simplify_rec2 (aut, e, a, d, int *n):
#    keep = False
#    *n = (*n) + 1
#    for u,v,l in aut.outgoing_edges(e):
#        if emonde0_simplify_rec(aut, e=v, a=a, d=d, n=n):
#            keep = True
    

class Automaton (DiGraph):
#    def __init__(self, I, F, **args):
#        DiGraph.__init__(self, **args)
#        self._I = I
#        self._F = F
#        return 
    def __repr__(self):
        return "Finite automaton with %d states"%self.num_verts()
        
    def _latex_(self):
        return "Finite automaton with $%d$ states"%self.num_verts()
    
    def plot (self):
        r"""
        Plot the automaton.
    
        EXAMPLES::
    
            sage: from sage.combinat.words.morphism import CallableDict
            sage: d = CallableDict({1:'one', 2:'zwei', 3:'trois'})
            sage: d(1), d(2), d(3)
            ('one', 'zwei', 'trois')
        """
        return DiGraph.plot(self, edge_labels=True)
    
    def plot2 (self, sx=10, sy=8, name="Automaton", verb=False):
        r"""
        Plot the automaton.
        """
        I = []
        if hasattr(self, 'I'):
            I = self.I
        F = []
        if hasattr(self, 'F'):
            F = self.F
        if verb:
            print I,F
        txt = '{\n'
        for e in self.vertices():
            txt += '    "'+str(e)+'" [shape='
            if e in F:
                txt += 'doublecircle'
            else:
                txt += 'circle'
            txt += ', style='
            if e in I:
                txt += 'bold'
            else:
                txt += 'solid'
            txt += ', fontsize=20, margin=0]\n'    
        txt += '\n'
        for e,f,l in self.edges():
            txt += '    "'+str(e)+'" -> "'+str(f)+'" [label="'+str(l)+'"];\n'
        txt += '}\n'
        SaveTikZ(txt, name, sx, sy)
    
    def plot3 (self, name="Automaton"):
        r"""
        Plot the automaton.
        """
        SaveTikZ(self.graphviz_string(edge_labels=True), name)
    
    def copy (self):
        a = DiGraph.copy(self)
        if hasattr(self, 'I'):
            a.I = self.I
        if hasattr(self, 'F'):
            a.F = self.F
        if hasattr(self, 'A'):
            a.A = self.A
        return a
    
    def Alphabet (self):
        if hasattr(self, 'A'):
            return self.A
        else:
            return set(self.edge_labels())
    
    def product2_ (self, A, d=None, verb=False):
        L = self.Alphabet()
        LA = A.Alphabet()
        if d is None:
            d = dict([])
            for k in L:
                for ka in LA:
                    d[(k, ka)] = [(k, ka)] #alphabet produit par défaut
#        if dA is None:
#            dA = dict([])
#            for k in L:
#                dA[k] = k
#        if d.keys() != dA.keys():
#            print "Error : incompatible keys : %s and %s"%(d,dA)
        P = Automaton(loops=True, multiedges=True)
        V = self.vertices()
        VA = A.vertices()
        for v in V:
            for va in VA:
                P.add_vertex((v,va))
                for v1, v2, l in self.outgoing_edges(v, labels=True):
                    for va1, va2, la in A.outgoing_edges(va, labels=True):
                        for e in d[(l, la)]:
                            if verb:
                                print "ajout de %s, %s, %s"%((v1, va1), (v2, va2), e)
                            P.add_edge((v1, va1), (v2, va2), e)
        if hasattr(self, 'I') and hasattr(A, 'I'):
            P.I = [(u,v) for u in self.I for v in A.I]
        if hasattr(self, 'F') and hasattr(A, 'F'):
            P.F = [(u,v) for u in self.F for v in A.F]
        P.A = set([a for l in L for la in LA if d.has_key((l,la)) for a in d[(l,la)]])
        #[u for u in set(d.values()) if u is not None]
        return P
    
    def product2 (self, A, d=None, verb=False):
        L = self.Alphabet()
        LA = A.Alphabet()
        if d is None:
            d = dict([])
            for k in L:
                for ka in LA:
                    d[(k, ka)] = (k, ka) #alphabet produit par défaut
        P = Automaton(loops=True, multiedges=True)
        S = set([(i, i2) for i in self.I for i2 in A.I]) #ensemble des sommets à explorer
        E = set([]) #ensemble des sommets déjà explorés
        while len(S) != 0:
            (v,va) = S.pop()
            P.add_vertex((v, va))
            for v1, v2, l in self.outgoing_edges(v, labels=True):
                for va1, va2, la in A.outgoing_edges(va, labels=True):
                    if d.has_key((l, la)):
                        #print "has_key (%s,%s)"%(l,la)
                        for e in d[(l, la)]:
                            if verb:
                                print "ajout de %s, %s, %s"%((v1, va1), (v2, va2), e)
                            P.add_edge((v1, va1), (v2, va2), e)
                        if len(d[(l, la)]) > 0 and (v2, va2) not in E:
                            S.add((v2, va2))
            E.add((v,va))
        if hasattr(self, 'I') and hasattr(A, 'I'):
            P.I = [(u,v) for u in self.I for v in A.I]
        if hasattr(self, 'F') and hasattr(A, 'F'):
            P.F = [(u,v) for u in self.F for v in A.F]
        #print "A=..."
        P.A = set([a for l in L for la in LA if d.has_key((l,la)) for a in d[(l,la)]])
        #print P.A
        #[u for u in set(d.values()) if u is not None]
        return P
    
    def product (self, A, d=None, verb=False):
        L = self.Alphabet()
        LA = A.Alphabet()
        if d is None:
            d = dict([])
            for k in L:
                for ka in LA:
                    d[(k, ka)] = (k, ka) #alphabet produit par défaut
        P = Automaton(loops=True, multiedges=True)
        S = set([(i, ia) for i in self.I for ia in A.I]) #ensemble des sommets à explorer
        E = set([]) #ensemble des sommets déjà explorés
        while len(S) != 0:
            if verb:
                print S, E
            (v,va) = S.pop()
            P.add_vertex((v, va))
            for v1, v2, l in self.outgoing_edges(v, labels=True):
                for va1, va2, la in A.outgoing_edges(va, labels=True):
                    if d[(l, la)] is not None:
                        P.add_edge((v1, va1), (v2, va2), d[(l, la)])
                        if (v2, va2) not in E:
                            S.add((v2, va2))
            E.add((v,va))
        if hasattr(self, 'I') and hasattr(A, 'I'):
            P.I = [(u,v) for u in self.I for v in A.I]
        if hasattr(self, 'F') and hasattr(A, 'F'):
            P.F = [(u,v) for u in self.F for v in A.F]
        P.A = [d[(l,la)] for l in L for la in LA if d[(l,la)] is not None]
        return P
    
    def product_ (self, A, d=None):
        L = self.Alphabet()
        LA = A.Alphabet()
        if d is None:
            d = dict([])
            for k in L:
                for ka in LA:
                    d[(k, ka)] = (k, ka) #alphabet produit par défaut
#        if dA is None:
#            dA = dict([])
#            for k in L:
#                dA[k] = k
#        if d.keys() != dA.keys():
#            print "Error : incompatible keys : %s and %s"%(d,dA)
        P = Automaton(loops=True, multiedges=True)
        V = self.vertices()
        VA = A.vertices()
        for v in V:
            for va in VA:
                P.add_vertex((v,va))
                for v1, v2, l in self.outgoing_edges(v, labels=True):
                    for va1, va2, la in A.outgoing_edges(va, labels=True):
                        if d[(l, la)] is not None:
                            #print "ajout de %s, %s, %s"%((v1, va1), (v2, va2), d[(l, la)])
                            P.add_edge((v1, va1), (v2, va2), d[(l, la)])
        if hasattr(self, 'I') and hasattr(A, 'I'):
            P.I = [(u,v) for u in self.I for v in A.I]
        if hasattr(self, 'F') and hasattr(A, 'F'):
            P.F = [(u,v) for u in self.F for v in A.F]
        P.A = [d[(l,la)] for l in L for la in LA if d[(l,la)] is not None]
        #[u for u in set(d.values()) if u is not None]
        return P
   
    def intersection (self, A, verb=False):
        L = self.Alphabet()
        LA = A.Alphabet()
        d = dict([])
        for k in L:
            for ka in LA:
                if k == ka:
                    d[(k, ka)] = k
                else:
                    d[(k, ka)] = None
        if verb:
            print d
        res = self.product(A=A, d=d, verb=verb)
        #if self.I is not None and A.I is not None:
        #if hasattr(self, 'I') and hasattr(A, 'I'):
        #    res.I = [(u, v) for u in self.I for v in A.I]
        return res
    
    def emonde0_rec (self, e):
        global E
        #print "e=%s, E=%s"%(e,E)
        #from sage.sets.set import set
        E.add(e)
        keep = False
        O = set([b for a,b,l in self.outgoing_edges(e)])
        #print "O=%s"%O
        for b in O:
            if self.has_vertex(b): #au cas où le sommet ait été supprimé entre-temps
                if b not in E:
                    self.emonde0_rec (e=b)
                if b in E:
                    keep = True
        if not keep:
            E.remove(e)
            self.delete_vertex(e)
    
    def emonde0 (self, I=None): #retire les etats puits
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                print "The set I of initial states must be defined !"
                return
        global E
           #from sage.sets.set import set
        E = set()
        for e in I:
            self.emonde0_rec(e=e)
    
    def emonde0_simplify_rec (self, e, a, d):
        keep = False
        og = self.outgoing_edges(e)
        O = set([v for u,v,l in og])
        cdef int ne = a.add_vertex()
        d[e] = ne
        for b in O:
            if self.has_vertex(b): #au cas où le sommet ait été supprimé entre-temps
                if not d.has_key(b):
                    self.emonde0_simplify_rec (e=b, a=a, d=d)
                if d.has_key(b):
                    keep = True
        if not keep:
            a.delete_vertex(d.pop(e))
        else:
            for u,v,l in og:
                if d.has_key(v):
                    if not a.has_edge(ne, d[v], l):
                        a.add_edge(ne, d[v], l)
    
    def emonde0_simplify (self, I=None): #rend l'automate sans états puits et avec des états plus simples
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                raise ValueError("The set I of initial states must be defined !")
        a = Automaton(loops=True, multiedges=True)
        d = dict()
        cdef int n = 0
        for e in I:
            #emonde0_simplify_rec(self, e=e, a=a, d=d, n=&n)
            self.emonde0_simplify_rec(e=e, a=a, d=d)
        a.I = [d[i] for i in I if d.has_key(i)]
        if hasattr(self, 'F'):
            a.F = [d[f] for f in self.F if d.has_key(f)]
        if hasattr(self, 'A'):
            a.A = self.A
        return a
    
    def emonde0_simplify_ (self, I=None, verb=False): #rend l'automate sans états puits et avec des états plus simples
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                raise ValueError("The set I of initial states must be defined !")
        a = Automaton(loops=True, multiedges=True)
        I = list(set(I))
        d = {}
        cdef int n = 0
        cdef int ni
        if verb:
            print "I=%s"%I
        while len(I) != 0:
            i = I.pop()
            if not d.has_key(i):
                ni = n
                d[i] = ni
                a.add_vertex(ni)
                n += 1
            else:
                ni = d[i]
            for u,v,l in self.outgoing_edges(i):
                if not d.has_key(v):
                    d[v] = n
                    a.add_vertex(n)
                    I.append(v)
                    n += 1
                a.add_edge(ni, d[v], l)
        a.I = [d[i] for i in self.I]
        if hasattr(self, 'F'):
            a.F = [d[f] for f in self.F if d.has_key(f)]
        if hasattr(self, 'A'):
            a.A = self.A
        return a
    
    def emondeI (self, I=None):
        #from sage.sets.set import set
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            else:
                raise ValueError("The initial states set I must be defined !")
        
        G = set()
        E = set(I)
        while E:
            #for e in E:
            #    E.remove(e)
            e = E.pop()
            if e not in G:
                G.add(e)
                E = E.union(set([b for (a,b,f) in self.outgoing_edges(e)]))
        for e in self.vertices():
            if e not in G:
                self.delete_vertex(e)
    
    def emondeF (self, F=None):
        #from sage.sets.set import set
        if F is None:
            if hasattr(self, 'F'):
                F = self.F
            else:
                raise ValueError("The final states set F must be defined !")
        
        G = set()
        E = set(F)
        while E:
            #for e in E:
            #    E.remove(e)
            e = E.pop()
            if e not in G:
                G.add(e)
                E = E.union(set([a for (a,b,f) in self.incoming_edges(e)]))
        for e in self.vertices():
            if e not in G:
                self.delete_vertex(e)
    
    def emonde (self, I=None, F=None):
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            else:
                raise ValueError("The initial states set I must be defined !")
        if F is None:
            if hasattr(self, 'F'):
                F = self.F
            else:
                raise ValueError("The final states set F must be defined !")
        
        self.emondeI(I=I)
        self.emondeF(F=F)
    
    def prefix (self, w, i=None, i2=None, verb=False): #retourne un automate dont le langage est w^(-1)L (self est supposé déterministe)
        if i is None:
            if hasattr(self, 'I'):
                i = self.I[0]
            else:
                raise ValueError("The initial state i must be given !")
        a = self.copy()
        e = i #etat initial
        if i2 is None:
            c2 = a.add_vertex()
            i2 = c2
        else:
            a.add_vertex(i2)
            c2 = i2
        for l in w:
            found = False
            if verb: print "e=%s"%e
            for u,v,l2 in self.outgoing_edges(e):
                if l2 == l:
                    found = True
                    e = v
                    c = c2
                    c2 = a.add_vertex()
                    a.add_edge(c, c2, l)
                    break
            if not found:
                print "Erreur : le mot %s n'est pas reconnu par l'automate !"%w
        #connecte
        a.add_edge(c, e, l)
        a.delete_vertex(c2)
        a.I = [i2]
        if hasattr(self, 'A'):
            a.A = self.A
        return a
    
    def succ (self, i, a): #donne un état atteint en partant de i et en lisant a (rend None si n'existe pas)
        if i is None:
            return
        for u,v,l in self.outgoing_edges(i):
            if l == a:
                return v
    
    def transpose (self):
        a  = Automaton(loops=True, multiedges=True)
        for u,v,l in self.edges():
            a.add_edge(v,u,l)
        if hasattr(self, 'F'):
            a.I = self.F
        if hasattr(self, 'I'):
            a.F = self.I
        if hasattr(self, 'A'):
            a.A = self.A
        return a
       
    def save_graphviz (self, file):
        a = self.copy()
        d=dict()
        for u in a.Alphabet():
            d[u]=[str(u)]
        a.graphviz_to_file_named(file, edge_labels=True)
           
    def relabel2 (self, d):
        self.allow_multiple_edges(True)
        for u,v,l in self.edges():
            self.delete_edge(u, v, l)
            if d.has_key(l):
                for l2 in d[l]:
                    self.add_edge(u, v, l2)
            else:
                print "Key %s not in the dictionary !"%l
        self.A = [a for A in d.values() for a in A]
    
    def determinize_rec (self, S, A, a, nof, noempty, verb):
        #from sage.sets.set import set
        #A2 = set([])
        
        if verb: print S
        
        if not a.has_key(S):
            a[S] = dict([])
            o = dict([]) #associe à chaque lettre l'ensemble des états partants de S
            for s in S:
                for f, d, l in self.outgoing_edges(s):
                    #A2 = A2.union(set([l]))
                    #if verb: print "l=%s"%l
                    if not o.has_key(l):
                        o[l] = set()
                    o[l].add(d)
           # for l in A:
           #     if l not in A2:
           #         if not a[S].has_key(set([])):
           #             a[S][set([])] = []
           #         a[S][set([])] += [l]
                        
            #for l in o.keys():
            for l in A:
                if o.has_key(l) or not noempty:
                    if not o.has_key(l):
                        o[l] = set()
                    #if verb: print "S=%s, l=%s, o[l]=%s"%(S, l, o[l])
                    if set(nof).isdisjoint(o[l]):
                        o[l] = Set(o[l])
                        if not a[S].has_key(o[l]):
                            a[S][o[l]] = []
                        a[S][o[l]] += [l] 
            for r in set(a[S].keys()):
                a = self.determinize_rec(S=r, A=A, a=a, nof=nof, noempty=noempty, verb=verb)
        return a
        
    def determinize (self, I=None, A=None, nof=set(), noempty=True, verb=False):
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                raise ValueError("The set I of initial states must be defined !")
        if A is None:
            if hasattr(self, 'A'):
                A = self.A
            if A is None:
                A = set(self.edge_labels())
                #raise ValueError("The alphabet A must be defined !")
        #from sage.sets.set import set
        
        a = dict()
        res = Automaton(self.determinize_rec (S=Set(I), A=A, a=a, nof=nof, noempty=noempty, verb=verb), multiedges=True, loops=True)
        res.I = [Set(I)]
        res.A = A
        if hasattr(self, 'F'):
            res.F = set([v for v in res.vertices() for f in self.F if f in v])
            #res.F = set()
            #for f in self.F:
            #    for v in res.vertices():
            #        if f in v:
            #            res.F.add(v)
        return res
    
    def determinize2 (self, I=None, A=None, nof=set(), noempty=True, verb=False):
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                raise ValueError("The set I of initial states must be defined !")
        if A is None:
            A = self.Alphabet()
                #raise ValueError("The alphabet A must be defined !")
        #from sage.sets.set import set
        
        nof = set(nof)
        
        c = 0
        a = Automaton(loops=True, multiedges=True)
        SS = set([Set(I)])
        a.I = set([Set(I)])
        a.A = A
        a.add_vertex(Set(I))
        while len(SS) != 0:
            S = SS.pop()
            o = dict([])
            for l in A:
                o[l] = set([])
            for s in S:
                for f, d, l in self.outgoing_edges(s):
                    o[l].add(d)
            for l in A:
                if o[l] != set([]) or not noempty:
                    if nof.isdisjoint(o[l]):
                        o[l] = Set(o[l])
                        if not a.has_vertex(o[l]):
                            a.add_vertex(o[l])
                            SS.add(o[l])
                            c+=1
                            if c%1000 == 0 and verb:
                                print c
                        a.add_edge(S, o[l], l)
        if hasattr(self, 'F'):
            a.F = set([v for v in a.vertices() for f in self.F if f in v])
        return a
    
    def complete (self, name=None, verb=False):
        r"""
        Complete the given automaton.
        If the given automaton is not complete, it will add a hole state named name.

        EXAMPLES::

            sage: BetaAdicMonoid(3, {0, 1, 3/2}).relations_automaton().complete()
            Finite automaton with 4 states
        """
        self.allow_multiple_edges(True)
        #from sage.sets.set import set
        L = set(self.Alphabet())
        V = self.vertices()
        ns = True
        for v in V:
            L2 = []
            for f, d, l in self.outgoing_edges(v):
                L2 += [l]
            if set(L2) != set(L):
                if verb: print "L=%s, L2=%s"%(L,L2)
                for l in L.difference(L2):
                    if ns:
                        ns = False
                        if name is None:
                            name = self.add_vertex()
                            V.append(name)
                        elif name in V:
                            raise ValueError("the name of the hole state must be unused")
                    self.add_edge(v, name, l)
        return name
    
    def minimize (self, F=None, verb=False):
        r"""
        Returns the minimal automaton.
        The given automaton is supposed to be deterministic.

        EXAMPLES::

            sage: BetaAdicMonoid(3, {0, 1, 3/2}).relations_automaton().complete().minimize()
            Finite automaton with 4 states
        """
        puit = None
        if F is None:
            if hasattr(self, 'F'):
                F = self.F
            if F is None:
                F = self.vertices()
        puit = self.complete()
        
        if verb: print "self=%s"%self
        
        #from sage.sets.set import set
        F = Set(F)
        Q = Set(self.vertices())
        P = set([F, Q.difference(F)])
        W = set([F])
        if hasattr(self, 'L'):
            L = self.L
        else:
            L = self.Alphabet()
        while W:
            if verb: print "P=%s, W=%s"%(P,W)
            A = W.pop()
            if verb: print "A=%s"%A
            X = dict()
            for s in A:
                for f, d, l in self.incoming_edges(s):
                    #for l in L2:
                    if not X.has_key(l):
                        X[l] = set()
                    X[l].add(f)
            #P2 = set([p for p in P])
            #while P2:
            #    Y = P2.pop()
            stop = False
            for l in L:
                if verb: print "P=%s, W=%s, l=%s"%(P,W,l)
                if X.has_key(l):
                    X[l] = Set(X[l])
                    
                    cont = True
                    while cont:
                        cont=False
                        for Y in P:
                            I = Y.intersection(X[l])
                            D = Y.difference(X[l])
                            if I and D:
                                if verb: print "I=%s, D=%s"%(I,D)
                                P.discard(Y)
                                P = P.union([I,D])
                                if Y in W:
                                    W.discard(Y)
                                    W |= set([I,D])
                                else:
                                    if len(I) <= len(D):
                                        W = W.union(set([I]))
                                    else:
                                        W = W.union(set([D]))
                                #stop = True
                                cont = True
                                break
                #if stop:
                #   break
#                    else:
#                        raise ValueError("the automaton must be complete")
        #construct the new automaton
        cs = dict()
        for S in P:
            for s in S:
                cs[s] = Set(S)
        a = dict()
        for S in P:
            S = Set(S)
            a[S] = dict()
            if not S.is_empty():
                for f, d, l in self.outgoing_edges(S.an_element()):
                    if not a[S].has_key(cs[d]):
                        a[S][cs[d]] = []
                    a[S][cs[d]] += [l]
        res = Automaton(a) #, loops=True, multiedges=True)
        res.F = [Set(S) for S in P if not set(S).isdisjoint(set(F))]
        res.A = L
        if hasattr(self, 'I'):
            if verb: print "I=%s"%I
            res.I = [Set(S) for S in P for i in self.I if i in S]
            if verb: print "res.I=%s"%res.I
            #for i in self.I:
            #    for v in res.vertices():
            #        if i in v:
            #            res.I += [v]
        if verb: print "res=%s"%res
        #delete the puit if it was not
        if puit is not None:
            if verb: print "effacement de l'etat puit %s"%puit
            self.delete_vertex(puit)
            for S in P:
                if puit in S:
                    if verb: print "delete %s"%S
                    res.delete_vertex(Set(S))
        res = res.emonde0_simplify()
        return res
     
    #compute the automaton recognizing the complementary
    def complementary (self):
        self.complete()
        if hasattr(self, 'F'):
            self.F = [f for f in self.vertices() if f not in self.F]
        else:
            self.F = self.vertices()
        
#    def plot (self, edge_labels=True, **options):
#        DiGraph.plot(self, edge_labels=edge_labels, **options)
    
#    def add_state (self, *args): #ajoute l'état s à l'automate (le réinitialise si déjà présent)
#        DiGraph.add_vertex(self, *args)   
