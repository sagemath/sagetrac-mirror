r"""
Convex H-Polyhedral Sets
"""

# ****************************************************************************
#       Copyright (C) 2008-2012 Marshall Hampton <hamptonio@gmail.com>
#       Copyright (C) 2011-2015 Volker Braun <vbraun.name@gmail.com>
#       Copyright (C) 2012-2018 Frederic Chapoton
#       Copyright (C) 2013      Andrey Novoseltsev
#       Copyright (C) 2014-2017 Moritz Firsching
#       Copyright (C) 2014-2019 Thierry Monteil
#       Copyright (C) 2015      Nathann Cohen
#       Copyright (C) 2015-2017 Jeroen Demeyer
#       Copyright (C) 2015-2017 Vincent Delecroix
#       Copyright (C) 2015-2018 Dima Pasechnik
#       Copyright (C) 2015-2020 Jean-Philippe Labbe <labbe at math.huji.ac.il>
#       Copyright (C) 2015-2022 Matthias Koeppe
#       Copyright (C) 2016-2019 Daniel Krenn
#       Copyright (C) 2017      Marcelo Forets
#       Copyright (C) 2017-2018 Mark Bell
#       Copyright (C) 2019      Julian Ritter
#       Copyright (C) 2019-2020 Laith Rastanawi
#       Copyright (C) 2019-2020 Sophia Elia
#       Copyright (C) 2019-2021 Jonathan Kliem <jonathan.kliem@fu-berlin.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.geometry.convex_set import ConvexSet_base
from sage.misc.abstract_method import abstract_method


class ConvexHPolyhedralSet_base(ConvexSet_base):
    """
    Abstract base class for convex H-polyhedral sets.

    Such a set is a convex set defined by finitely many linear equations and
    linear inequalities, some of which may be strict.

    A ``ConvexHPolyhedralSet_base`` instances has a distinguished H-representation,
    which is not necessarily minimal.

    Equality of ``ConvexHPolyhedralSet_base`` instances is equality as sets,
    not equality of representations.
    """

    @abstract_method
    def equation_generator(self):
        """
        Return a generator for the linear equations satisfied by the
        polyhedron.

        EXAMPLES::

            sage: p = polytopes.regular_polygon(8,base_ring=RDF)
            sage: p3 = Polyhedron(vertices = [x+[0] for x in p.vertices()], base_ring=RDF)
            sage: next(p3.equation_generator())
            An equation (0.0, 0.0, 1.0) x + 0.0 == 0
        """

    def equations(self):
        """
        Return all linear constraints of the polyhedron.

        OUTPUT:

        A tuple of equations.

        EXAMPLES::

            sage: test_p = Polyhedron(vertices = [[1,2,3,4],[2,1,3,4],[4,3,2,1],[3,4,1,2]])
            sage: test_p.equations()
            (An equation (1, 1, 1, 1) x - 10 == 0,)
        """
        return tuple(self.equation_generator())

    def n_equations(self):
        """
        Return the number of equations.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]])
            sage: p.n_equations()
            1
        """
        return len(self.equations())

    @abstract_method
    def inequality_generator(self):
        """
        Return a generator for the defining non-strict inequalities.

        OUTPUT:

        A generator of the inequality H-representation objects.

        EXAMPLES::

            sage: triangle = Polyhedron(vertices=[[1,0],[0,1],[1,1]])
            sage: for v in triangle.inequality_generator(): print(v)
            An inequality (1, 1) x - 1 >= 0
            An inequality (0, -1) x + 1 >= 0
            An inequality (-1, 0) x + 1 >= 0
            sage: [ v for v in triangle.inequality_generator() ]
            [An inequality (1, 1) x - 1 >= 0,
             An inequality (0, -1) x + 1 >= 0,
             An inequality (-1, 0) x + 1 >= 0]
            sage: [ [v.A(), v.b()] for v in triangle.inequality_generator() ]
            [[(1, 1), -1], [(0, -1), 1], [(-1, 0), 1]]
        """

    def inequalities(self):
        """
        Return all non-strict inequalities.

        OUTPUT:

        A tuple of inequalities.

        EXAMPLES::

            sage: p = Polyhedron(vertices = [[0,0,0],[0,0,1],[0,1,0],[1,0,0],[2,2,2]])
            sage: p.inequalities()[0:3]
            (An inequality (1, 0, 0) x + 0 >= 0,
             An inequality (0, 1, 0) x + 0 >= 0,
             An inequality (0, 0, 1) x + 0 >= 0)

            sage: p3 = Polyhedron(vertices=Permutations([1, 2, 3, 4]))          # optional - sage.combinat
            sage: ieqs = p3.inequalities()                                      # optional - sage.combinat
            sage: ieqs[0]                                                       # optional - sage.combinat
            An inequality (0, 1, 1, 1) x - 6 >= 0
            sage: list(_)                                                       # optional - sage.combinat
            [-6, 0, 1, 1, 1]
        """
        return tuple(self.inequality_generator())

    def n_inequalities(self):
        """
        Return the number of non-strict inequalities.

        EXAMPLES::

            sage: p = Polyhedron(vertices=[[1,0,0], [0,1,0], [0,0,1]])
            sage: p.n_inequalities()
            3

            sage: p = Polyhedron(vertices=[[t,t^2,t^3] for t in range(6)])
            sage: p.n_facets()
            8
        """
        return len(self.inequalities())
