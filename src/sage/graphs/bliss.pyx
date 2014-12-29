r"""

Interface for graph automorphisms with bliss.

Implemented functions:

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Returns the automorphism group of the given (di)graph
    :meth:`canonical_form` | Computes a canonical certificate for the given (di) graph. 
    :meth:`is_isomorpic` | Tests whether the passed (di) graphs are isomorphic. 


AUTHORS:

    - Jernej Azarija 

"""
include "sage/ext/interrupt.pxi"
include 'sage/ext/stdsage.pxi'
include 'sage/ext/cdefs.pxi'

cdef extern from "graph.hh" namespace "bliss":

    cdef cppclass Stats:
        pass
    
    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int, const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int, const unsigned int*), void*)
        # FIXME: See is_isomorphic
        int cmp(Graph& )

    cdef cppclass Digraph(AbstractGraph):
        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*)(void* , unsigned int, const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);
        const unsigned int* canonical_form(Stats&, void (*)(void*,unsigned int, const unsigned int*), void*)
        
# What we recieve here is a permutation given as a bijection and we create the respective
# cycle notation for it. NOTE it would be nice if this could be made as fast as possible.

# FIXME. One can use Permutation().to_cylce here I am not sure I like it though since we'd have
# to translate labelings back and forth and as I said I want this to be efficent.
# I think at this point it would be best to keep this and modify Permutation to support labeled stuff.
# Maybe add a parameter to Permutation apply_labeling=dict or something.   
cdef void add_gen(void *user_param , unsigned int n, const unsigned int *aut):

    covered = cur = 0
    perm = []       
    done = [False]*n

    gens, d2 = <object> <PyObject *> user_param
    
    while covered < n:
        while cur < n and done[cur]:
            cur+=1

        cycle = [d2[aut[cur]]]
        sec = aut[aut[cur]]
        done[aut[cur]] = done[cur] = True

        covered+=1

        while d2[sec] != cycle[0]:
            cycle.append(d2[sec])
            done[sec] = True
            sec = aut[sec]
            covered+=1
        perm+=[tuple(cycle)]
    gens += [perm]

# The function accepts a graph, a coloring  of its vertices (possibly None)
# and two empty dictionaries. The entries of the dicitionary are later set
# to record the labeling of our graph. They are taken as arguments to avoid
# technicalities of returning Python objects in Cython functions.
cdef Graph *bliss_graph(G, partition, vert2int, int2vert):    

    n = G.order()
    
    cdef Graph *g = new Graph(n)

    # FIXME if g == NULL ....

    for i,v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v 

    for x,y in G.edges(labels=False):
       g.add_edge(vert2int[x],vert2int[y])     

    if partition:
        for i in xrange(1,len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

# This is the same function as bliss_graph with the only exception
# being that it returns a digraph. This code duplication is of course
# ugly and I am open for suggestions (that do not waste time)            
cdef Digraph *bliss_digraph(G, partition, vert2int, int2vert):    

    n = G.order()
    cdef Digraph *g = new Digraph(n) 

    for i,v in enumerate(G):
        vert2int[v] = i
        int2vert[i] = v 

    # FIXME if g == NULL ....

    for x,y in G.edges(labels=False):
        g.add_edge(vert2int[x],vert2int[y])     

    if partition:
        for i in xrange(1,len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

# FIXME unknown sppliting herustic abortion error if a digraph is passed a few times??
"""
sage: D = digraphs.Kautz(10,2,1)
    sage: from sage.graphs.bliss import *
    sage: D.automorphism_group().is_isomorphic(automorphism_group(D))
    True
    sage: D.automorphism_group().is_isomorphic(automorphism_group(D))
    Bliss fatal error: Internal error - unknown splitting heuristics
    Aborting!

"""    
def automorphism_group(G, partition=None):

    cv = 0 
    n = G.order()
    vert2int = {}
    int2vert = {}

    cdef Graph *g = NULL
    cdef Digraph *d = NULL
    cdef Stats s

    gens = [] 
    data = (gens, int2vert)        

    if G.is_directed(): 
        d = bliss_digraph(G, partition, vert2int, int2vert)
        d.find_automorphisms(s, add_gen, <PyObject *> data)
    else:
        g = bliss_graph(G, partition, vert2int, int2vert) 
        g.find_automorphisms(s, add_gen, <PyObject *> data)
    del g 
    del d

    from sage.groups.perm_gps.permgroup import PermutationGroup

    return PermutationGroup(gens)

cdef void empty_hook(void *user_param , unsigned int n, const unsigned int *aut):
    return

# Returns a Python object such that G and H are isomorphic
# if and only if canonical_form(G) == canonical_form(H).
# If return_graph=True the respective object is a Graph

# If return_labeling=True returns the actual labeling that
# gives you the canonical form of G.
def canonical_form(G, partition=None, return_graph=False, return_labeling=False):

    cdef Graph *g = NULL
    cdef Digraph *d = NULL
    cdef const unsigned int *aut
    cdef Stats s

    vert2int = {}
    
    if G.is_directed():
        d = bliss_digraph(G, partition, vert2int, {}) 
        aut = d.canonical_form(s, empty_hook, NULL)
    else:
        g = bliss_graph(G, partition, vert2int, {})
        aut = g.canonical_form(s, empty_hook, NULL)     
     
    # FIXME is copy + relabel perhaps faster?
    edges = [(aut[vert2int[x]], aut[vert2int[y]]) for x,y in G.edges(labels=False)] 

    del g
    del d 

    if return_graph:
        if G.is_directed():
            from sage.graphs.graph import DiGraph
            G = DiGraph(edges)
        else:
            from sage.graphs.graph import Graph
            G = Graph(edges)

        return G, vert2int if return_labeling else G

    if return_labeling:
        return sorted(edges),vert2int

    return sorted(edges)

# Returns true if and only if G and H are isomorphic
# If cert=True return the respective automorphism in case G
# and H are isomorphic and  None otherwise.
# NOTE we assume G is directed <=> H is directed

# FIXME if cert=False then the best way to do this would actually
# be to compute the two canonical forms of G and H and use bliss::Graph::cmp
# which I somehow cannot make it work? Can someone look into that?
def is_isomorphic(G,H, cert=False):
    
    c1,lab1 = canonical_form(G, return_labeling=True)
    c2,lab2 = canonical_form(H, return_labeling=True)
   
    if c1 == c2:
        if cert:
            lab2inv = { lab2[key]:key for key in lab2}
            return { v: lab2inv[lab1[v]] for v in G}
        else:
            return True
    return False  
