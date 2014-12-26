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

cdef extern from "/home/azi/bliss-0.72/graph.hh" namespace "bliss":

    cdef cppclass Stats:
        pass
    
    cdef cppclass AbstractGraph:
        pass

    cdef cppclass Graph(AbstractGraph):
        Graph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*hook)(void* , unsigned int, const unsigned int*), void*)

    cdef cppclass Digraph(AbstractGraph):
        Digraph(const unsigned int)
        void add_edge(const unsigned int, const unsigned int)
        void find_automorphisms(Stats&, void (*hook)(void* , unsigned int, const unsigned int*), void*)



# What we recieve here is a permutation given as a bijection and we create the respective
# cycle notation for it. NOTE it would be nice if this could be made as fast as possible.

# FIXME. One can use Permutation().to_cylce here I am not sure I like it though since we'd have
# to translate labelings back and forth and as I said I want this to be fast.
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


def automorphism_group(G,partition=None):

    cv = 0 
    n = G.order()
    d = {}
    d2 = {}

    cdef Graph *g = NULL
    cdef Digraph *h = NULL
    cdef Stats s

    # FIXME can this directed/undirected thing be solved in a more elegant way?
    isDir = G.is_directed() 

    if isDir:
        h = new Digraph(n)
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
        if isDir: 
            h.add_edge(d[x],d[y])     
        else:
            g.add_edge(d[x],d[y])     

    gens = [] 
    data = (gens, d2)        

    if isDir:
       h.find_automorphisms(s, add_gen, <PyObject *> data)
    else:
       g.find_automorphisms(s, add_gen, <PyObject *> data)
    
    del g 
    del h

    from sage.groups.perm_gps.permgroup import PermutationGroup
    return PermutationGroup(gens)

