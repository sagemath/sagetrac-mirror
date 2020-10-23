r"""
View classes

This module implements views for (di)graphs. A view is a read-only iterable
container enabling operations like ``for e in E`` and ``e in E``. It is updated
as the graph is updated. Hence, the graph should not be updated while iterating
through a view. Views can be iterated multiple times.

.. TODO::

    - View of neighborhood to get open/close neighborhood of a vertex/set of
      vertices, and also the vertex boundary

Classes
-------
"""
# ****************************************************************************
#       Copyright (C) 2019 David Coudert <david.coudert@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from __future__ import print_function, absolute_import

from itertools import islice
from sys import maxsize as sys_maxsize
from sage.graphs.generic_graph_pyx cimport GenericGraph_pyx

cdef class EdgesView:
    r"""
    EdgesView class.

    This class implements a read-only iterable container of edges enabling
    operations like ``for e in E`` and ``e in E``. An :class:`EdgesView` can be
    iterated multiple times, and checking membership is done in constant
    time. It avoids the construction of edge lists and so consumes little
    memory. It is updated as the graph is updated. Hence, the graph should not
    be updated while iterating through an :class:`EdgesView`.

    INPUT:

    - ``G`` -- a (di)graph

    - ``vertices`` -- list (default: ``None``); an iterable container of
      vertices or ``None``. When set, consider only edges incident to specified
      vertices.

    - ``labels`` -- boolean (default: ``True``); if ``False``, each edge is
      simply a pair ``(u, v)`` of vertices

    - ``ignore_direction`` -- boolean (default: ``False``); only applies to
      directed graphs. If ``True``, searches across edges in either direction.

    - ``sort`` -- boolean (default: ``None``); whether to sort edges

      - if ``None``, sort edges according to the default ordering and give a
        deprecation warning as sorting will be set to ``False`` by default in
        the future

      - if ``True``, edges are sorted according the ordering specified with
        parameter ``key``

      - if ``False``, edges are not sorted. This is the fastest and less memory
        consuming method for iterating over edges. This will become the default
        behavior in the future.

    - ``key`` -- a function (default: ``None``); a function that takes an edge
      (a pair or a triple, according to the ``labels`` keyword) as its one
      argument and returns a value that can be used for comparisons in the
      sorting algorithm. This parameter is ignored when ``sort = False``.

    .. WARNING::

        Since any object may be a vertex, there is no guarantee that any two
        vertices will be comparable, and thus no guarantee how two edges may
        compare. With default objects for vertices (all integers), or when all
        the vertices are of the same simple type, then there should not be a
        problem with how the vertices will be sorted. However, if you need to
        guarantee a total order for the sorting of the edges, use the ``key``
        argument, as illustrated in the examples below.

    EXAMPLES::

        sage: from sage.graphs.views import EdgesView
        sage: G = Graph([(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')])
        sage: E = EdgesView(G, sort=True); E
        [(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')]
        sage: (1, 2) in E
        False
        sage: (1, 2, 'B') in E
        True
        sage: E = EdgesView(G, labels=False, sort=True); E
        [(0, 1), (0, 2), (1, 2)]
        sage: (1, 2) in E
        True
        sage: (1, 2, 'B') in E
        False
        sage: [e for e in E]
        [(0, 1), (0, 2), (1, 2)]

    An :class:`EdgesView` can be iterated multiple times::

        sage: G = graphs.CycleGraph(3)
        sage: print(E)
        [(0, 1), (0, 2), (1, 2)]
        sage: print(E)
        [(0, 1), (0, 2), (1, 2)]
        sage: for e in E:
        ....:     for ee in E:
        ....:         print((e, ee))
        ((0, 1), (0, 1))
        ((0, 1), (0, 2))
        ((0, 1), (1, 2))
        ((0, 2), (0, 1))
        ((0, 2), (0, 2))
        ((0, 2), (1, 2))
        ((1, 2), (0, 1))
        ((1, 2), (0, 2))
        ((1, 2), (1, 2))

    We can check if a view is empty::

        sage: E = EdgesView(graphs.CycleGraph(3), sort=False)
        sage: if E:
        ....:     print('not empty')
        not empty
        sage: E = EdgesView(Graph(), sort=False)
        sage: if not E:
        ....:     print('empty')
        empty

    When ``sort`` is ``True``, edges are sorted by default in the default
    fashion::

        sage: G = Graph([(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')])
        sage: E = EdgesView(G, sort=True); E
        [(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')]

    This can be overridden by specifying a key function. This first example just
    ignores the labels in the third component of the triple::

        sage: G = Graph([(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')])
        sage: E = EdgesView(G, sort=True, key=lambda x: (x[1], -x[0])); E
        [(0, 1, 'C'), (1, 2, 'B'), (0, 2, 'A')]

    We can also sort according to the labels::

        sage: G = Graph([(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')])
        sage: E = EdgesView(G, sort=True, key=lambda x: x[2]); E
        [(0, 2, 'A'), (1, 2, 'B'), (0, 1, 'C')]

    With a directed graph::

        sage: G = digraphs.DeBruijn(2, 2)
        sage: E = EdgesView(G, labels=False, sort=True); E
        [('00', '00'), ('00', '01'), ('01', '10'), ('01', '11'),
         ('10', '00'), ('10', '01'), ('11', '10'), ('11', '11')]
        sage: E = EdgesView(G, labels=False, sort=True, key=lambda e:(e[1], e[0])); E
        [('00', '00'), ('10', '00'), ('00', '01'), ('10', '01'),
         ('01', '10'), ('11', '10'), ('01', '11'), ('11', '11')]

    We can consider only edges incident to a specified set of vertices::

        sage: G = graphs.CycleGraph(5)
        sage: E = EdgesView(G, vertices=[0, 1], labels=False, sort=True); E
        [(0, 1), (0, 4), (1, 2)]
        sage: E = EdgesView(G, vertices=[0], labels=False, sort=True); E
        [(0, 1), (0, 4)]
        sage: E = EdgesView(G, vertices=None, labels=False, sort=True); E
        [(0, 1), (0, 4), (1, 2), (2, 3), (3, 4)]

        sage: G = digraphs.Circuit(5)
        sage: E = EdgesView(G, vertices=[0, 1], labels=False, sort=True); E
        [(0, 1), (1, 2)]

    We can ignore the direction of the edges of a directed graph, in which case
    we search accross edges in either direction::

        sage: G = digraphs.Circuit(5)
        sage: E = EdgesView(G, vertices=[0, 1], labels=False, sort=True, ignore_direction=False); E
        [(0, 1), (1, 2)]
        sage: (1, 0) in E
        False
        sage: E = EdgesView(G, vertices=[0, 1], labels=False, sort=True, ignore_direction=True); E
        [(0, 1), (0, 1), (1, 2), (4, 0)]
        sage: (1, 0) in E
        True
        sage: G.has_edge(1, 0)
        False

    A view is updated as the graph is updated::

        sage: G = Graph()
        sage: E = EdgesView(G, vertices=[0, 3], labels=False, sort=True); E
        []
        sage: G.add_edges([(0, 1), (1, 2)])
        sage: E
        [(0, 1)]
        sage: G.add_edge(2, 3)
        sage: E
        [(0, 1), (2, 3)]

    Hence, the graph should not be updated while iterating through a view::

        sage: G = Graph([('a', 'b'), ('b', 'c')])
        sage: E = EdgesView(G, labels=False, sort=False); E
        [('a', 'b'), ('b', 'c')]
        sage: for u, v in E:
        ....:     G.add_edge(u + u, v + v)
        Traceback (most recent call last):
        ...
        RuntimeError: dictionary changed size during iteration

    Two :class:`EdgesView` are considered equal if they report either both
    directed, or both undirected edges, they have the same settings for
    ``ignore_direction``, they have the same settings for ``labels``, and they
    report the same edges in the same order::

        sage: G = graphs.HouseGraph()
        sage: EG = EdgesView(G, sort=False)
        sage: H = Graph(EG)
        sage: EH = EdgesView(H, sort=False)
        sage: EG == EH
        True
        sage: G.add_edge(0, 10)
        sage: EG = EdgesView(G, sort=False)
        sage: EG == EH
        False
        sage: H.add_edge(0, 10)
        sage: EH = EdgesView(H, sort=False)
        sage: EG == EH
        True
        sage: H = G.strong_orientation()
        sage: EH = EdgesView(H, sort=False)
        sage: EG == EH
        False

    The sum of two :class:`EdgesView` is a list containing the edges in both
    :class:`EdgesView`::

        sage: E1 = EdgesView(Graph([(0, 1)]), labels=False, sort=False)
        sage: E2 = EdgesView(Graph([(2, 3)]), labels=False, sort=False)
        sage: E1 + E2
        [(0, 1), (2, 3)]
        sage: E2 + E1
        [(2, 3), (0, 1)]

    Recall that a :class:`EdgesView` is read-only and that this method
    returns a list::

        sage: E1 += E2
        sage: type(E1) is list
        True

    It is also possible to get the sum a :class:`EdgesView` with itself `n`
    times::

        sage: E = EdgesView(Graph([(0, 1), (2, 3)]), labels=False, sort=True)
        sage: E * 3
        [(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)]
        sage: 3 * E
        [(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)]

    Recall that a :class:`EdgesView` is read-only and that this method
    returns a list::

        sage: E *= 2
        sage: type(E) is list
        True

    We can ask for the `i`-th edge, or a slice of the edges as a list::

        sage: E = EdgesView(graphs.HouseGraph(), labels=False, sort=True)
        sage: E[0]
        (0, 1)
        sage: E[2]
        (1, 3)
        sage: E[-1]
        (3, 4)
        sage: E[1:-1]
        [(0, 2), (1, 3), (2, 3), (2, 4)]
        sage: E[::-1]
        [(3, 4), (2, 4), (2, 3), (1, 3), (0, 2), (0, 1)]
    """
    cdef readonly GenericGraph_pyx _graph
    cdef list _vertices
    cdef frozenset _vertex_set
    cdef readonly bint _labels
    cdef readonly bint _ignore_direction
    cdef bint _sort_edges
    cdef _sort_edges_key

    def __init__(self, G, vertices=None, labels=True, ignore_direction=False,
                     sort=None, key=None):
        """
        Construction of this :class:`EdgesView`.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = Graph([(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')])
            sage: E = EdgesView(G, sort=True); E
            [(0, 1, 'C'), (0, 2, 'A'), (1, 2, 'B')]

        TESTS::

            sage: EdgesView(Graph(), sort=False, key=lambda x:0)
            Traceback (most recent call last):
            ...
            ValueError: sort keyword is not True, yet a key function is given
        """
        self._graph = <GenericGraph_pyx?>G

        if vertices is None:
            self._vertices = None
            self._vertex_set = None
        else:
            self._vertices = list(vertices)
            self._vertex_set = frozenset(self._vertices)

        # None and 0 are interpreted as False, 1 as True
        self._labels = labels
        self._ignore_direction = ignore_direction
        self._sort_edges = sort
        if not sort and key is not None:
            raise ValueError('sort keyword is not True, yet a key function is given')
        self._sort_edges_key = key

    def __len__(self):
        """
        Return the number of edges in ``self``.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: E = EdgesView(G)
            sage: len(E)
            6
            sage: len(E) == G.size()
            True
            sage: E = EdgesView(G, vertices=[0], labels=False, sort=True); E
            [(0, 1), (0, 2)]
            sage: len(E)
            2
            sage: G = digraphs.Circuit(4)
            sage: len(EdgesView(G, ignore_direction=False, sort=False))
            4
            sage: len(EdgesView(G, ignore_direction=True, sort=False))
            8
        """
        if self._vertices is self._graph:
            if self._graph._directed and self._ignore_direction:
                return 2 * self._graph.size()
            return self._graph.size()
        else:
            return sum(1 for _ in self)

    def __repr__(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: E = EdgesView(G, labels=False)
            sage: repr(E)
            '[(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]'
            sage: E
            [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
        """
        return "[%s]" % ', '.join(map(repr, self))

    def _iter_unsorted(self, vertices):
        """
        Iterator over the unsorted edges in ``self``.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: E = EdgesView(G, labels=False, sort=False)
            sage: list(E)
            [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
        """
        if self._graph._directed:
            yield from self._graph._backend.iterator_out_edges(vertices, self._labels)
            if self._ignore_direction:
                yield from self._graph._backend.iterator_in_edges(vertices, self._labels)
        else:
            yield from self._graph._backend.iterator_edges(vertices, self._labels)

    def __iter__(self):
        """
        Iterator over the edges in ``self``.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: E = EdgesView(G, labels=False, sort=True)
            sage: list(E)
            [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
            sage: sum(1 for e in E for ee in E) == len(E) * len(E)
            True
            sage: E = EdgesView(G, labels=False)
            sage: sum(1 for e in E for ee in E) == G.size() * G.size()
            True
        """
        if self._vertices is None:
            vertices = self._graph
        else:
            vertices = self._vertices
        if self._sort_edges:
            yield from sorted(self._iter_unsorted(vertices), key=self._sort_edges_key)
        else:
            yield from self._iter_unsorted(vertices)

    def __bool__(self):
        """
        Check whether ``self`` is empty.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(Graph())
            sage: bool(E)
            False
            sage: E = EdgesView(graphs.HouseGraph())
            sage: bool(E)
            True
        """
        if self._vertices is None:
            vertices = self._graph
        else:
            vertices = self._vertices
        if self._graph._directed:
            for _ in self._graph._backend.iterator_out_edges(vertices, False):
                return True
            if self._ignore_direction:
                for _ in self._graph._backend.iterator_in_edges(vertices, False):
                    return True
        else:
            for _ in self._graph._backend.iterator_edges(vertices, False):
                return True
        return False

    def __eq__(self, right):
        """
        Check whether ``self`` and ``right`` are equal.

        Do not call this method directly. That is, for ``E1.__eq__(E2)`` write
        ``E1 == E2``.

        Two :class:`EdgesView` are considered equal if the following hold:
        - they report either both directed, or both undirected edges;
        - they have the same settings for ``ignore_direction``;
        - they have the same settings for ``labels``;
        - they report the same edges in the same order.

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: EG = EdgesView(G)
            sage: H = Graph(EG)
            sage: EH = EdgesView(H)
            sage: EG == EH
            True
            sage: G.add_edge(0, 10)
            sage: EG = EdgesView(G)
            sage: EG == EH
            False
            sage: H.add_edge(0, 10)
            sage: EH = EdgesView(H)
            sage: EG == EH
            True
            sage: H = G.strong_orientation()
            sage: EH = EdgesView(H)
            sage: EG == EH
            False

        TESTS::

            sage: from sage.graphs.views import EdgesView
            sage: G = graphs.HouseGraph()
            sage: E = EdgesView(G)
            sage: E == E
            True
            sage: E == G
            False
            sage: G == E
            False

        Check that :trac:`29180` is fixed::

            sage: G = graphs.CycleGraph(4)
            sage: E = graphs.EmptyGraph()
            sage: G.edges() == E.edges()
            False
        """
        if not isinstance(right, EdgesView):
            return NotImplemented
        cdef EdgesView other = <EdgesView>right
        if self is other:
            return True
        # Check parameters
        if (self._graph._directed != other._graph._directed or
            self._ignore_direction != other._ignore_direction or
            self._labels != other._labels):
            return False
        # Check that self and other have the same number of edges
        if len(self) != len(other):
            return False
        # Check that the same edges are reported in the same order
        return all(es == eo for es, eo in zip(self, other))

    def __contains__(self, e):
        """
        Check whether edge ``e`` is part of ``self``.

        INPUT:

        - ``e`` -- tuple; when ``self.labels`` is ``True``, the expected form is
          ``(u, v, label)``, while when ``self.labels`` is ``False``, the
          expected form is ``(u, v)``

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: G = Graph([(0, 1)])
            sage: E = EdgesView(G, labels=False)
            sage: print(E)
            [(0, 1)]
            sage: (0, 1) in E
            True
            sage: (0, 1, None) in E
            False
            sage: E = EdgesView(G, labels=True)
            sage: print(E)
            [(0, 1, None)]
            sage: (0, 1) in E
            False
            sage: (0, 1, None) in E
            True
        """
        if self._labels:
            try:
                u, v, label = e
            except Exception:
                return False
        else:
            try:
                u, v = e
            except Exception:
                return False
            label = None
        if (self._vertex_set is not None
                and u not in self._vertex_set and v not in self._vertex_set):
            return False
        if self._graph._directed and self._ignore_direction:
            return (self._graph._backend.has_edge(u, v, label)
                        or self._graph._backend.has_edge(v, u, label))
        return self._graph._backend.has_edge(u, v, label)

    def __getitem__(self, i):
        r"""
        Return the `i`-th edge in ``self``.

        This method takes time `O(i)`. When several calls to this method are
        done, prefer making ``list`` from ``self`` before querying items.

        INPUT:

        - ``i`` -- integer or slice

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(graphs.HouseGraph(), labels=False)
            sage: E[0]
            (0, 1)
            sage: E[2]
            (1, 3)
            sage: E[1:-1]
            [(0, 2), (1, 3), (2, 3), (2, 4)]
            sage: E[-1]
            (3, 4)

        TESTS::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(graphs.HouseGraph(), labels=False)
            sage: E[10]
            Traceback (most recent call last):
            ...
            IndexError: index out of range
        """
        cdef Py_ssize_t start, stop, step
        if isinstance(i, slice):
            start, stop, step = i.start or 0, i.stop or sys_maxsize, i.step or 1
            if start >= 0 and stop >= 0 and step >= 0:
                return list(islice(self, start, stop, step))
            else:
                return list(self)[i]
        elif i < 0:
            return list(self)[i]
        else:
            i = int(i)  # For Python < 3.7 where islice doesn't support non-int
            try:
                return next(islice(self, i, i + 1, 1))
            except StopIteration:
                raise IndexError('index out of range')

    def __add__(left, right):
        """
        Return a list containing the edges of ``left`` and ``right``.

        The returned list contains the edges of ``left`` with prescribed order
        followed by the edges of ``right`` in prescribed order.

        INPUT:

        - ``left,right`` -- :class:`EdgesView` or list of edges

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: E1 = EdgesView(Graph([(0, 1)]), labels=False)
            sage: E2 = EdgesView(Graph([(2, 3)]), labels=False)
            sage: E1 + E2
            [(0, 1), (2, 3)]
            sage: E2 + E1
            [(2, 3), (0, 1)]
            sage: E1 + E2 + E1
            [(0, 1), (2, 3), (0, 1)]

        Recall that a :class:`EdgesView` is read-only and that this method
        returns a list::

            sage: E1 += E2
            sage: type(E1) is list
            True

        TESTS::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(graphs.HouseGraph())
            sage: E + 'foo'
            Traceback (most recent call last):
            ...
            TypeError: unsupported operand type(s) for +: 'sage.graphs.views.EdgesView' and 'str'
        """
        if not isinstance(left, (list, EdgesView)) or not isinstance(right, (list, EdgesView)):
            return NotImplemented
        cdef list L = list(left)
        L.extend(right)
        return L

    def __mul__(left, right):
        r"""
        Return the sum of ``left`` with itself ``right`` times.

        INPUT:

        - ``right`` -- integer

        EXAMPLES::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(Graph([(0, 1), (2, 3)]), labels=False, sort=True)
            sage: E * 3
            [(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)]
            sage: 3 * E
            [(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)]

        Recall that a :class:`EdgesView` is read-only and that this method
        returns a list::

            sage: E *= 2
            sage: type(E) is list
            True

        TESTS::

            sage: from sage.graphs.views import EdgesView
            sage: E = EdgesView(Graph([(0, 1)]))
            sage: E * (-1)
            []
            sage: E * 1.5
            Traceback (most recent call last):
            ...
            TypeError: can...t multiply sequence by non-int of type 'sage.rings.real_mpfr.RealLiteral'
        """
        if isinstance(left, EdgesView):
            return list(left) * right
        else:
            # Case __rmul__
            return list(right) * left
