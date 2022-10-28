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


cdef class GLPKGraphBackend():
    cdef glp_graph * graph
    cdef add_vertex(self, name = ?)
    cdef list add_vertices(self, vertices)
    cdef __add_vertices_sage(self, g)
    cdef dict get_vertex(self, vertex)
    cdef dict get_vertices(self, verts)
    cdef set_vertex_demand(self, vertex, param)
    cdef set_vertices_demand(self, list pairs)
    cdef list vertices(self)
    cdef add_edge(self, u, v, dict params = ?)
    cdef __add_edges_sage(self, g)
    cdef list add_edges(self, edges)
    cdef delete_edge(self, u, v, dict params = ?)
    cdef tuple get_edge(self, u, v)
    cdef list edges(self)
    cdef delete_vertex(self, vert)
    cdef delete_vertices(self, list verts)
    cdef int _find_vertex(self, vert)
    cdef int write_graph(self, fname)
    cdef int write_ccdata(self, fname)
    cdef int write_mincost(self, fname)
    cdef double mincost_okalg(self) except -1
    cdef int s
    cdef int t
    cdef int write_maxflow(self, fname) except -1
    cdef double maxflow_ffalg(self, u = ?, v = ?) except -1
    cdef double cpp(self)
