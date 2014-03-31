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
    
    def copy (self):
        a = DiGraph.copy(self)
        if hasattr(self, 'I'):
            a.I = self.I
        if hasattr(self, 'F'):
            a.F = self.F
        if hasattr(self, 'A'):
            a.A = self.A
        return a
    
    def product (self, A, d=None):
        if d is None:
            L = self.edge_labels()
            LA = A.edge_labels()
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
        P.A = [u for u in set(d.values()) if u is not None]
        return P
   
    def intersection(self, A):
        L = self.edge_labels()
        LA = A.edge_labels()
        d = dict([])
        for k in L:
            for ka in LA:
                if k == ka:
                    d[(k, ka)] = k
                else:
                    d[(k, ka)] = None
        res = self.product(A=A, d=d)
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
        #from sage.sets.set import set
        O = set([v for u,v,l in self.outgoing_edges(e)])
        #print "O=%s"%O
        #d[e] = "a%s"%a.num_verts()
        #a.add_vertex(d[e])
        d[e] = a.add_vertex()
        for b in O:
            if self.has_vertex(b): #au cas où le sommet ait été supprimé entre-temps
                if not d.has_key(b):
                    self.emonde0_simplify_rec (e=b, a=a, d=d)
                if d.has_key(b):
                    keep = True
        if not keep:
            a.delete_vertex(d.pop(e))
        else:
            for u,v,l in self.outgoing_edges(e):
                if d.has_key(v):
                    if not a.has_edge(d[u], d[v], l):
                        a.add_edge(d[u], d[v], l)
    
    def emonde0_simplify (self, I=None): #rend l'automate sans états puits et avec des états plus simples
        if I is None:
            if hasattr(self, 'I'):
                I = self.I
            if I is None:
                raise ValueError("The set I of initial states must be defined !")
        a = Automaton(loops=True, multiedges=True)
        d = dict()
        for e in I:
            self.emonde0_simplify_rec(e=e, a=a, d=d)
        a.I = [d[i] for i in I if d.has_key(i)]
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
        for u,v,l in self.outgoing_edges(i):
            if l == a:
                return v
    
    def transpose (self):
        a = self.copy()
        for u,v,l in self.edges():
            a.delete_edge(u,v,l)
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
        for u in a.edge_labels():
            d[u]=[str(u)]
        a.graphviz_to_file_named(file, edge_labels=True)
           
    def relabel2 (self, d):
        self.allow_multiple_edges(True)
        for u,v,l in self.edges():
            self.delete_edge(u, v, l)
            for l2 in d[l]:
                self.add_edge(u, v, l2)
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
        L = set(self.edge_labels())
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
            L = self.edge_labels()
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
        
#    def plot (self, edge_labels=True, **options):
#        DiGraph.plot(self, edge_labels=edge_labels, **options)
    
#    def add_state (self, *args): #ajoute l'état s à l'automate (le réinitialise si déjà présent)
#        DiGraph.add_vertex(self, *args)   
