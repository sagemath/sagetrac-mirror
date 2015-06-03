#*****************************************************************************
#       Copyright (C) 2015 Michele Borassi michele.borassi@imtlucca.it
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


cdef extern from "boost/graph/adjacency_list.hpp" namespace "boost":
    cdef cppclass BoostVecGenGraph "boost::adjacency_list<boost::vecS, boost::vecS> ":
        pass
    cdef cppclass BoostVecGraph "boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> ":
        pass
    cdef cppclass BoostVecDiGraph "boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> ":
        pass
    
cdef extern from "boost/graph/graph_traits.hpp" namespace "boost":
    cdef cppclass graph_traits[T]:
        cppclass vertex_descriptor:
            pass
        cppclass edge_descriptor:
            pass
        cppclass degree_size_type:
            pass

ctypedef graph_traits[BoostVecGraph].vertex_descriptor vertex_descriptor_und
ctypedef graph_traits[BoostVecDiGraph].vertex_descriptor vertex_descriptor_dir

ctypedef graph_traits[BoostVecGraph].edge_descriptor edge_descriptor_und
ctypedef graph_traits[BoostVecDiGraph].edge_descriptor edge_descriptor_dir

ctypedef graph_traits[BoostVecGenGraph].degree_size_type degree_size_type
ctypedef graph_traits[BoostVecDiGraph].degree_size_type degree_size_type_dir

cdef extern from "boost/graph/adjacency_list.hpp" namespace "boost":
    int num_vertices(BoostVecGraph)
    int num_vertices(BoostVecDiGraph)

    int num_edges(BoostVecGraph)
    int num_edges(BoostVecDiGraph)

    vertex_descriptor_und add_vertex(BoostVecGraph)
    vertex_descriptor_dir add_vertex(BoostVecDiGraph)
    
    edge_descriptor_und add_edge(vertex_descriptor_und, 
                                 vertex_descriptor_und, 
                                 BoostVecGraph)
    edge_descriptor_dir add_edge(vertex_descriptor_dir, 
                                 vertex_descriptor_dir, 
                                 BoostVecDiGraph)
    
    cdef vertex_descriptor_und source(edge_descriptor_und, BoostVecGraph)
    cdef vertex_descriptor_und target(edge_descriptor_und, BoostVecGraph)

    cdef vertex_descriptor_dir source(edge_descriptor_dir, BoostVecDiGraph)
    cdef vertex_descriptor_dir target(edge_descriptor_dir, BoostVecDiGraph)
    
cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        T& at(int)        
        void push_back(T&) except +
        T& operator[](int)

ctypedef vector[edge_descriptor_und] edge_container_t_und
ctypedef vector[edge_descriptor_dir] edge_container_t_dir

cdef extern from "<iterator>" namespace "std":
    cdef cppclass back_insert_iterator[T]:
        back_insert_iterator(T)
    back_insert_iterator[T] back_inserter[T](T e)


cdef extern from "boost/graph/edge_connectivity.hpp" namespace "boost":
    graph_traits[T].degree_size_type edge_connectivity[T,TT](T G, TT GG)

cdef class BoostGraph(object):
    cdef BoostVecGraph *graph
    cdef vector[vertex_descriptor_und] vertices
    cpdef add_vertex(self)
    cpdef add_edge(self, int, int)
    cpdef add_edge_unsafe(self, int, int)
    cpdef num_verts(self)
    cpdef num_edges(self)
    cpdef edge_connectivity(self)
    
    
cdef class BoostDiGraph(object):
    cdef BoostVecDiGraph *graph
    cdef vector[vertex_descriptor_dir] vertices
    cpdef add_vertex(self)
    cpdef add_edge(self, int, int)
    cpdef add_edge_unsafe(self, int, int)
    cpdef num_verts(self)
    cpdef num_edges(self)
    cpdef edge_connectivity(self)
