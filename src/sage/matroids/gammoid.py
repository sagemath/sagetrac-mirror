r"""
Gammoids
"""
from __future__ import absolute_import
#*****************************************************************************
#       Copyright (C) 2013 Rudi Pendavingh <rudi.pendavingh@gmail.com>
#       Copyright (C) 2013 Michael Welsh <michael@welsh.co.nz>
#       Copyright (C) 2013 Stefan van Zwam <stefanvanzwam@gmail.com>
#
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from .matroid import Matroid
from .utilities import sanitize_contractions_deletions, setprint_s

from copy import copy, deepcopy

class Gammoid(Matroid):
    r"""
    Gammoid

    INPUT:

    - ``D`` -- A DiGraph representing the gammoid.
    - ``roots`` -- A subset of the vertices.
    - ``groundset`` -- (optional) A subset of the vertices. If not specified,
      the entire vertex set is used (and the gammoid will be strict).
    """

    def __init__(self, D, roots, groundset = None):
        """
        Documentation will be under class definition.
        """
        self._roots = frozenset(roots)
        vertices = frozenset(D.vertices())
        if not self._roots.issubset(vertices):
            raise ValueError("roots must be a subset of the vertices")

        if groundset is None:
            self._groundset = frozenset(vertices)
        else:
            self._groundset = frozenset(groundset)
            if not self._groundset.issubset(vertices):
                raise ValueError("ground set must be a subset of the vertices")

        self._G = D.copy(immutable=True)
        self._D = copy(D)
        self._rootv = self._D.add_vertex()
        for v in roots:
            self._D.add_edge(v, self._rootv)

    def digraph(self):
        """
        Return the DiGraph associated with the gammoid.
        """
        return self._G.copy(data_structure='sparse')

    def graph_plot(self):
        """
        Plot the graph with color-coded vertices.
        """
        self._inter = frozenset(self._roots.intersection(self._groundset))
        self._starting = frozenset(self._roots.difference(self._inter))
        self._ending = frozenset(self._groundset.difference(self._inter))
        self._therest = vertices.difference(self._roots)
        self._therest = frozenset(self._therest.difference(self._groundset))
        # Vertices just in starting set are red "#D55E00"
        # Vertices just in the ending set are blue "#0072B2"
        # Vertices in both are pink "#CC79A7"
        # Vertices in neither are grey "#999999"
        d = {"#D55E00": list(self._starting), "#CC79A7": list(self._inter),
            "#0072B2": list(self._ending), "#999999": list(self._therest)}
        return self._G.plot(vertex_colors = d)

    def _rank(self, X):
        """
        Return the rank of a set.

        EXAMPLES::

            sage: from sage.matroids.gammoid import Gammoid
            sage: D = digraphs.TransitiveTournament(4)
            sage: M = Gammoid(D, roots=[2,3])
            sage: M.rank([2])
            1
            sage: M.full_rank()
            2
            sage: M.rank(M.groundset())
            2

        """
        source = self._D.add_vertex()
        for x in X:
            self._D.add_edge(source, x)

        rank = len(self._D.vertex_disjoint_paths(source, self._rootv))
        self._D.delete_vertex(source)
        return rank

    def groundset(self):
        return self._groundset
