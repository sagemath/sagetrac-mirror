r"""
Interface for graph automorphisms with bliss.

Index
-----

**Cython functions**

**Python functions**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`automorphism_group` | Returns the automorphism group of the given graph 

Functions
---------

Author:
- Jernej Azarija (2014)

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
        void find_automorphisms(Stats&, void (*hook)(void* , unsigned int, const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);

    cdef cppclass Digraph(AbstractGraph):
        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*hook)(void* , unsigned int, const unsigned int*), void*)
        void change_color(const unsigned int, const unsigned int);


# What we recieve here is a permutation given as a bijection and we create the respective
# cycle notation for it. NOTE it would be nice if this could be made as fast as possible.

# FIXME. One can use Permutation().to_cylce here I am not sure I like it though since we'd have
# to translate labelings back and forth and as I said I want this to be fast.
# I think at this point it would be best to keep this and modify Permutation to support labeled stuff        
cdef void add_gen(void *user_param , unsigned int n, const unsigned int *aut):

    covered = cur = 0
    perm = []       
    done = [0]*n

    gens, d2 = <object> <PyObject *> user_param
    
    while covered < n:
        while cur < n and done[cur] == -1:
            cur+=1

        cycle = [d2[aut[cur]]]
        sec = aut[aut[cur]]
        done[aut[cur]] = -1
        done[cur] = -1
        covered+=1

        while d2[sec] != cycle[0]:
            cycle.append(d2[sec])
            done[sec] = -1
            sec = aut[sec]
            covered+=1
        perm+=[tuple(cycle)]
    gens += [perm]


cdef Graph *bliss_graph(G,partition, vert2int, int2vert):    

    n = G.order()
    cdef Graph *g = new Graph(n)

    for x,y in G.edges(labels=False):
        if x not in vert2int:
            vert2int[x] = cv
            int2vert[cv] = x
            cv+=1
        if y not in vert2int:
            vert2int[y] = cv
            int2vert[cv] = y
            cv+=1
        g.add_edge(vert2int[x],vert2int[y])     
    if partition:
        for i in xrange(1,len(partition)):
            for v in partition[i]:
                g.change_color(vert2int[v], i)
    return g

def automorphism_group(G,partition=None):

    cv = 0 
    n = G.order()
    d = {}
    d2 = {}

    cdef Graph *g = NULL
    cdef Digraph *h = NULL
    cdef Stats s

    # FIXME can this directed/undirected thing be solved in a more elegant way?
    # Its really really ugly. I guess the solution would be to patch bliss to make add_edge
    # and change_color part of the AbstractGraph class?
    isDir = G.is_directed() 

    if G.is_directed(): 
        h = new Digraph(n)
        for x,y in G.edges(labels=False):
            if x not in d:
                d[x] = cv
                d2[cv] = x
                cv+=1
            if y not in d:
                d[y] = cv
                d2[cv] = y
                cv+=1
            h.add_edge(d[x],d[y])     
        if partition:
            for i in xrange(1,len(partition)):
                for v in partition[i]:
                    h.change_color(d2[v], i)
    else:
        g = new Graph(n)
        for x,y in G.edges(labels=False):
            if x not in d:
                d[x] = cv
                d2[cv] = x
                cv+=1
            if y not in d:
                d[y] = cv
                d2[cv] = y
                cv+=1
            g.add_edge(d[x],d[y])     

        if partition:
            for i in xrange(1,len(partition)):
                for v in partition[i]:
                    g.change_color(d2[v], i)

    gens = [] 
    data = (gens, d2)        

    if G.is_directed():
       h.find_automorphisms(s, add_gen, <PyObject *> data)
    else:
       g.find_automorphisms(s, add_gen, <PyObject *> data)
    del g 
    del h

    from sage.groups.perm_gps.permgroup import PermutationGroup
    return PermutationGroup(gens)
