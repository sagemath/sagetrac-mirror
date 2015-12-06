#*****************************************************************************
#       Copyright (C) 2012 Christian Kuper <christian.kuper@t-online.de>
#       Copyright (C) 2015 Jeroen Demeyer <jdemeyer@cage.ugent.be>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.glpk.types cimport glp_graph

ctypedef struct c_v_data:
         double rhs
         double pi
         double es
         double ls
         long cut

ctypedef struct c_a_data:
         double low
         double cap
         double cost
         double x


cdef class GLPKGraphBackend(object):
    cdef glp_graph * graph
    cpdef add_vertex(self, char* name = ?)
    cpdef list add_vertices(self, vertices)
    cpdef __add_vertices_sage(self, g)
    cpdef dict get_vertex(self, char* vertex)
    cpdef dict get_vertices(self, verts)
    cpdef set_vertex_demand(self, char* vertex, param)
    cpdef set_vertices_demand(self, list pairs)
    cpdef list vertices(self)
    cpdef add_edge(self, char* u, char* v, dict params = ?)
    cpdef __add_edges_sage(self, g)
    cpdef list add_edges(self, edges)
    cpdef delete_edge(self, char* u, char* v, dict params = ?)
    cpdef tuple get_edge(self, char* u, char* v)
    cpdef list edges(self)
    cpdef delete_vertex(self, char* vert)
    cpdef delete_vertices(self, list verts)
    cpdef int _find_vertex(self, char *)
    cpdef int write_graph(self, char *fname)
    cpdef int write_ccdata(self, char *fname)
    cpdef int write_mincost(self, char *fname)
    cpdef double mincost_okalg(self) except -1
    cdef int s
    cdef int t
    cpdef int write_maxflow(self, char *fname) except -1
    cpdef double maxflow_ffalg(self, u = ?, v = ?) except -1
    cpdef double cpp(self)
