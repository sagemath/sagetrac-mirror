r"""
Points of Topological Manifolds

The class :class:`ManifoldPoint` implements points of a
topological manifold.
A :class:`ManifoldPoint` object can have coordinates in
various charts defined on the manifold. Two points are declared equal if they
have the same coordinates in the same chart.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- [Lee11]_ J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed.,
  Springer (New York) (2011)
- [Lee13]_ J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed.,
  Springer (New York, 2013)

EXAMPLES:

Defining a point in `\RR^3` by its spherical coordinates::

    sage: M = Manifold(3, 'R^3', structure='topological')
    sage: U = M.open_subset('U')  # the complement of the half-plane (y=0, x>=0)
    sage: c_spher.<r,th,ph> = U.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')
    sage: p = U((1, pi/2, pi), name='P') # coordinates in U's default chart (c_spher)
    sage: p
    Point P on the 3-dimensional topological manifold R^3
    sage: latex(p)
    P
    sage: p in U
    True
    sage: p.parent()
    Open subset U of the 3-dimensional topological manifold R^3
    sage: c_spher(p)
    (1, 1/2*pi, pi)
    sage: p.coord(c_spher) # equivalent to above
    (1, 1/2*pi, pi)

Computing the point's coordinates in a new chart::

    sage: c_cart.<x,y,z> = U.chart() # Cartesian coordinates on U
    sage: spher_to_cart = c_spher.transition_map(c_cart,
    ....:                    [r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th)])
    sage: c_cart(p)  # evaluate P's Cartesian coordinates
    (-1, 0, 0)

Points can be compared::

    sage: p1 = U((1, pi/2, pi))
    sage: p == p1
    True
    sage: q = U((1,2,3), chart=c_cart, name='Q') # point defined by its Cartesian coordinates
    sage: p == q
    False

"""

#*****************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element

class ManifoldPoint(Element):
    r"""
    Point of a topological manifold.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.manifold.TopologicalManifold`
    or :class:`~sage.manifolds.subset.ManifoldSubset`

    INPUT:

    - ``parent`` -- the manifold subset to which the point belongs
    - ``coords`` -- (default: ``None``) the point coordinates (as a tuple
      or a list) in the chart ``chart``
    - ``chart`` -- (default: ``None``) chart in which the coordinates are
      given; if  ``None``, the coordinates are assumed to refer to the
      default chart of ``parent``
    - ``name`` -- (default: ``None``) name given to the point
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the point;
      if ``None``, the LaTeX symbol is set to ``name``
    - ``check_coords`` -- (default: ``True``) determines whether ``coords`` are
      valid coordinates for the chart ``chart``; for symbolic coordinates, it
      is recommended to set ``check_coords`` to ``False``

    EXAMPLES:

    A point on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M', structure='topological')
        sage: c_xy.<x,y> = M.chart()
        sage: (a, b) = var('a b') # generic coordinates for the point
        sage: p = M.point((a, b), name='P'); p
        Point P on the 2-dimensional topological manifold M
        sage: p.coord()  # coordinates of P in the subset's default chart
        (a, b)

    Since points are Sage *elements*, the *parent* of which being the
    subset on which they are defined, it is equivalent to write::

        sage: p = M((a, b), name='P'); p
        Point P on the 2-dimensional topological manifold M

    A point is an element of the manifold subset in which it has been defined::

        sage: p in M
        True
        sage: p.parent()
        2-dimensional topological manifold M
        sage: U = M.open_subset('U', coord_def={c_xy: x>0})
        sage: q = U.point((2,1), name='q')
        sage: q.parent()
        Open subset U of the 2-dimensional topological manifold M
        sage: q in U
        True
        sage: q in M
        True

    By default, the LaTeX symbol of the point is deduced from its name::

        sage: latex(p)
        P

    But it can be set to any value::

        sage: p = M.point((a, b), name='P', latex_name=r'\mathcal{P}')
        sage: latex(p)
        \mathcal{P}

    Points can be drawn in 2D or 3D graphics thanks to the method :meth:`plot`.

    """
    def __init__(self, parent, coords=None, chart=None, name=None,
                 latex_name=None, check_coords=True):
        r"""
        Construct a manifold point.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,3), name='p'); p
            Point p on the 2-dimensional topological manifold M
            sage: TestSuite(p).run()
            sage: U = M.open_subset('U', coord_def={X: x<0})
            sage: q = U((-1,2), name='q'); q
            Point q on the 2-dimensional topological manifold M
            sage: TestSuite(q).run()

        """
        Element.__init__(self, parent)
        self._manifold = parent.manifold()  # a useful shortcut
        self._coordinates = {} # dictionary of the point coordinates in various
                               # charts, with the charts as keys
        if coords is not None:
            if len(coords) != parent.manifold().dimension():
                raise ValueError("the number of coordinates must be equal " +
                                 "to the manifold's dimension")
            from sage.manifolds.manifold import TopologicalManifold
            if chart is None:
                chart = parent._def_chart
            elif isinstance(parent, TopologicalManifold):
                if chart not in parent._atlas:
                    raise ValueError("the {} has not been".format(chart) +
                                     "defined on the {}".format(parent))
            if check_coords:
                if not chart.valid_coordinates(*coords):
                    raise ValueError("the coordinates {}".format(coords) +
                                     " are not valid on the {}".format(chart))
            for schart in chart._supercharts:
                self._coordinates[schart] = tuple(coords)
            for schart in chart._subcharts:
                if schart != chart:
                    if schart.valid_coordinates(*coords):
                        self._coordinates[schart] = tuple(coords)
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name

    def _repr_(self):
        r"""
        Return a string representation of the point.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3))
            sage: p._repr_()
            'Point on the 2-dimensional topological manifold M'
            sage: p = M((2,-3), name='p')
            sage: p._repr_()
            'Point p on the 2-dimensional topological manifold M'
            sage: repr(p)  # indirect doctest
            'Point p on the 2-dimensional topological manifold M'

        """
        description = "Point"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._manifold)
        return description

    def _latex_(self):
        r"""
        Return a LaTeX representation of the point.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3))
            sage: p._latex_()
            '\\mbox{Point on the 2-dimensional topological manifold M}'
            sage: p = M((2,-3), name='p')
            sage: p._latex_()
            'p'
            sage: p = M((2,-3), name='p', latex_name=r'\mathcal{P}')
            sage: p._latex_()
            '\\mathcal{P}'
            sage: latex(p)  # indirect doctest
            \mathcal{P}

        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        return self._latex_name


    def coord(self, chart=None, old_chart=None):
        r"""
        Return the point coordinates in the specified chart.

        If these coordinates are not already known, they are computed from
        known ones by means of change-of-chart formulas.

        An equivalent way to get the coordinates of a point is to let the
        chart acting of the point, i.e. if ``X`` is a chart and ``p`` a
        point, one has ``p.coord(chart=X) == X(p)``.

        INPUT:

        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart
        - ``old_chart`` -- (default: ``None``) chart from which the
          coordinates in ``chart`` are to be computed; if ``None``, a chart
          in which the point's coordinates are already known will be picked,
          privileging the subset's default chart

        EXAMPLES:

        Spherical coordinates of a point on `\RR^3`::

            sage: M = Manifold(3, 'M', structure='topological') # the part of R^3 covered by spherical coordinates
            sage: c_spher.<r,th,ph> = M.chart(r'r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi') # spherical coordinates
            sage: p = M.point((1, pi/2, pi))
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, pi)
            sage: p.coord(c_spher) # with the chart c_spher specified (same result as above since this is the default chart)
            (1, 1/2*pi, pi)

        An alternative way to get the coordinates is to let the chart act
        on the point (from the very definition of a chart)::

            sage: c_spher(p)
            (1, 1/2*pi, pi)

        Computing the Cartesian coordinates from the spherical ones::

            sage: c_cart.<x,y,z> = M.chart()  # Cartesian coordinates
            sage: c_spher.transition_map(c_cart, [r*sin(th)*cos(ph),
            ....:                                 r*sin(th)*sin(ph), r*cos(th)])
            Change of coordinates from Chart (M, (r, th, ph)) to Chart (M, (x, y, z))
            sage: p.coord(c_cart)  # the computation is performed by means of the above change of coordinates
            (-1, 0, 0)
            sage: p.coord(c_cart) == c_cart(p)
            True

        Coordinates of a point on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: c_xy.<x,y> = M.chart()
            sage: (a, b) = var('a b') # generic coordinates for the point
            sage: p = M.point((a, b), name='P')
            sage: p.coord()  # coordinates of P in the manifold's default chart
            (a, b)

        Coordinates of P in a new chart::

            sage: c_uv.<u,v> = M.chart()
            sage: ch_xy_uv = c_xy.transition_map(c_uv, [x-y, x+y])
            sage: p.coord(c_uv)
            (a - b, a + b)

        Coordinates of P in a third chart::

            sage: c_wz.<w,z> = M.chart()
            sage: ch_uv_wz = c_uv.transition_map(c_wz, [u^3, v^3])
            sage: p.coord(c_wz, old_chart=c_uv)
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)

        Actually, in the present case, it is not necessary to specify
        old_chart='uv'::

            sage: p.set_coord((a-b, a+b), c_uv) # erases all the coordinates except those in the chart c_uv
            sage: p._coordinates
            {Chart (M, (u, v)): (a - b, a + b)}
            sage: p.coord(c_wz)
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)
            sage: p._coordinates  # random (dictionary output)
            {Chart (M, (u, v)): (a - b, a + b),
             Chart (M, (w, z)): (a^3 - 3*a^2*b + 3*a*b^2 - b^3,
                                 a^3 + 3*a^2*b + 3*a*b^2 + b^3)}

        """
        if chart is None:
            dom = self.parent()
            chart = dom._def_chart
            def_chart = chart
        else:
            dom = chart._domain
            def_chart = dom._def_chart
            if self not in dom:
                raise ValueError("the point does not belong to the domain " +
                                 "of {}".format(chart))
        if chart not in self._coordinates:
            # Check whether chart corresponds to a superchart of a chart
            # in which the coordinates are known:
            for ochart in self._coordinates:
                if chart in ochart._supercharts or chart in ochart._subcharts:
                    self._coordinates[chart] = self._coordinates[ochart]
                    return self._coordinates[chart]
            # If this point is reached, some change of coordinates must be
            # performed
            if old_chart is not None:
                s_old_chart = old_chart
                s_chart = chart
            else:
                # A chart must be found as a starting point of the computation
                # The domain's default chart is privileged:
                if def_chart in self._coordinates \
                        and (def_chart, chart) in dom._coord_changes:
                    old_chart = def_chart
                    s_old_chart = def_chart
                    s_chart = chart
                else:
                    for ochart in self._coordinates:
                        for subchart in ochart._subcharts:
                            if (subchart, chart) in dom._coord_changes:
                                old_chart = ochart
                                s_old_chart = subchart
                                s_chart = chart
                                break
                        if old_chart is not None:
                            break
                if old_chart is None:
                    # Some search involving the subcharts of chart is
                    # performed:
                    for schart in chart._subcharts:
                        for ochart in self._coordinates:
                            for subchart in ochart._subcharts:
                                if (subchart, schart) in dom._coord_changes:
                                    old_chart = ochart
                                    s_old_chart = subchart
                                    s_chart = schart
                                    break
                            if old_chart is not None:
                                break
                        if old_chart is not None:
                            break
            if old_chart is None:
                raise ValueError("the coordinates of {}".format(self) +
                          " in the {}".format(chart) + " cannot be computed " +
                          "by means of known changes of charts.")
            else:
                chcoord = dom._coord_changes[(s_old_chart, s_chart)]
                self._coordinates[chart] = chcoord(*self._coordinates[old_chart])
        return self._coordinates[chart]

    def set_coord(self, coords, chart=None):
        r"""
        Sets the point coordinates in the specified chart.

        Coordinates with respect to other charts are deleted, in order to
        avoid any inconsistency. To keep them, use the method :meth:`add_coord`
        instead.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart

        EXAMPLES:

        Setting coordinates to a point on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point()
            sage: p.set_coord((2,-3)) # coordinates on the manifold's default chart
            sage: p.coord()
            (2, -3)
            sage: X(p)
            (2, -3)

        Let us introduce a second chart on the manifold::

            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])

        If we set the coordinates of p in the chart Y, those in the chart X
        are lost::

            sage: Y(p)
            (-1, 5)
            sage: p.set_coord(Y(p), chart=Y)
            sage: p._coordinates
            {Chart (M, (u, v)): (-1, 5)}

        """
        self._coordinates.clear()
        self.add_coord(coords, chart)

    def add_coord(self, coords, chart=None):
        r"""
        Adds some coordinates in the specified chart.

        The previous coordinates with respect to other charts are kept. To
        clear them, use :meth:`set_coord` instead.

        INPUT:

        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart`` -- (default: ``None``) chart in which the coordinates are
          given; if none is provided, the coordinates are assumed to refer to
          the subset's default chart

        .. WARNING::

           If the point has already coordinates in other charts, it
           is the user's responsibility to make sure that the coordinates
           to be added are consistent with them.

        EXAMPLES:

        Setting coordinates to a point on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point()
            sage: p.add_coord((2,-3)) # coordinates on the manifold's default chart
            sage: p.coord()
            (2, -3)
            sage: X(p)
            (2, -3)

        Let us introduce a second chart on the manifold::

            sage: Y.<u,v> = M.chart()
            sage: X_to_Y = X.transition_map(Y, [x+y, x-y])

        If we add the coordinates of p in the chart Y, those in the chart X
        are kept::

            sage: p.add_coord((-1,5), chart=Y)
            sage: p._coordinates  # random (dictionary output)
            {Chart (M, (u, v)): (-1, 5), Chart (M, (x, y)): (2, -3)}

        On the contrary, with the method :meth:`set_coord`, the coordinates
        in charts different from Y would be lost::

            sage: p.set_coord((-1,5), chart=Y)
            sage: p._coordinates
            {Chart (M, (u, v)): (-1, 5)}

        """
        if len(coords) != self.parent().manifold()._dim:
            raise ValueError("the number of coordinates must be equal to " +
                             "the manifold's dimension.")
        if chart is None:
            chart = self.parent()._def_chart
        else:
            if chart not in self.parent()._atlas:
                raise ValueError("the {}".format(chart) + " has not been " +
                                 "defined on the {}".format(self.parent()))
        self._coordinates[chart] = coords

    def __eq__(self, other):
        r"""
        Compares the current point with another one.

        EXAMPLES:

        Comparison with coordinates in the same chart::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3), chart=X)
            sage: q = M((2,-3), chart=X)
            sage: p == q
            True
            sage: q = M((-2,-3), chart=X)
            sage: p == q
            False

        Comparison with coordinates of other in a subchart::

            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: XU = X.restrict(U)
            sage: q = U((2,-3), chart=XU)
            sage: p == q and q == p
            True
            sage: q = U((1,-3), chart=XU)
            sage: p == q or q == p
            False

        Comparison requiring a change of chart::

            sage: Y.<u,v> = U.chart()
            sage: XU_to_Y = XU.transition_map(Y, (ln(x), x+y))
            sage: XU_to_Y.inverse()(u,v)
            (e^u, v - e^u)
            sage: q = U((ln(2),-1), chart=Y)
            sage: p == q and q == p
            True
            sage: q = U((ln(3),1), chart=Y)
            sage: p == q or q == p
            False

        """
        if other is self:
            return True
        if not isinstance(other, ManifoldPoint):
            return False
        if other.parent().manifold() != self.parent().manifold():
            return False
        # Search for a common chart to compare the coordinates
        common_chart = None
        # the subset's default chart is privileged:
        # FIXME: Make this a better test
        if hasattr(self.parent(), '_def_chart'):  # self.parent() is open
            def_chart = self.parent()._def_chart
        else:
            def_chart = self.parent().manifold()._def_chart
        if def_chart in self._coordinates and def_chart in other._coordinates:
            common_chart = def_chart
        else:
            for chart in self._coordinates:
                if chart in other._coordinates:
                    common_chart = chart
                    break
        if common_chart is None:
            # At this stage, a commont chart is searched via a coordinate
            # transformation:
            for chart in self._coordinates:
                try:
                    other.coord(chart)
                    common_chart = chart
                    break
                except ValueError:
                    pass
            else:
                # Attempt a coordinate transformation in the reverse way:
                for chart in other._coordinates:
                    try:
                        self.coord(chart)
                        common_chart = chart
                        break
                    except ValueError:
                        pass
        if common_chart is None:
            return False
            #!# Another option would be:
            # raise ValueError("no common chart has been found to compare " +
            #                  "{} and {}".format(self, other))
        return self._coordinates[common_chart] == other._coordinates[common_chart]

    def __ne__(self, other):
        r"""
        Non-equality operator.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3), chart=X)
            sage: q = M((0,1), chart=X)
            sage: p != q
            True
            sage: p != M((2,-3), chart=X)
            False

        """
        return not (self == other)

    def __cmp__(self, other):
        r"""
        Old-style (Python 2) comparison operator.

        This is provisory, until migration to Python 3 is achieved.

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3), chart=X)
            sage: q = M((2,-3), chart=X)
            sage: p.__cmp__(q)
            0
            sage: q = M((0,1), chart=X)
            sage: p.__cmp__(q)
            -1

        """
        if self == other:
            return 0
        else:
            return -1

    def __hash__(self):
        r"""
        This hash function is set to constant on a given manifold, to fulfill
        Python's credo:
        p == q  ==>  hash(p) == hash(q)
        This is necessary since p and q may be created in different coordinate
        systems and nevertheless be equal

        TESTS::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M((2,-3), chart=X)
            sage: hash(p)  # random
            8791657334475
            sage: hash(p) == hash(M)
            True

        """
        return hash(self.parent().manifold())


    def plot(self, chart=None, ambient_coords=None, mapping=None, size=10,
             color='black', label=None, label_color=None, fontsize=10,
             label_offset=0.1, parameters=None):
        r"""
        For real manifolds, plot the current point in a Cartesian graph based
        on the coordinates of some ambient chart.

        The point is drawn in terms of two (2D graphics) or three (3D graphics)
        coordinates of a given chart, called hereafter the *ambient chart*.
        The domain of the ambient chart must contain the point, or its image
        by a continuous manifold map `\Phi`.

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, the ambient chart is set the default chart of
          ``self.parent()``
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2 or 3
          coordinates of the ambient chart in terms of which the plot is
          performed; if ``None``, all the coordinates of the ambient chart are
          considered
        - ``mapping`` -- (default: ``None``) continuous manifold map `\Phi`
          (instance of :class:`~sage.manifolds.continuous_map.ContinuousMap`)
          providing the link between the current point `p` and the ambient
          chart ``chart``: the domain of ``chart`` must contain `\Phi(p)`;
          if ``None``, the identity map is assumed
        - ``size`` -- (default: 10) size of the point once drawn as a small
          disk or sphere
        - ``color`` -- (default: 'black') color of the point
        - ``label`` -- (default: ``None``) label printed next to the point;
          if ``None``, the point's name is used.
        - ``label_color`` -- (default: ``None``) color to print the label;
          if ``None``, the value of ``color`` is used
        - ``fontsize`` -- (default: 10) size of the font used to print the
          label
        - ``label_offset`` -- (default: 0.1) determines the separation between
          the point and its label
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the point coordinates

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of the ambient chart) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of the ambient chart)

        EXAMPLES:

        Drawing a point on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M', structure='topological')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,3), name='p')
            sage: g = p.plot(X)
            sage: print g
            Graphics object consisting of 2 graphics primitives
            sage: gX = X.plot(max_range=4) # plot of the coordinate grid
            sage: show(g+gX) # display of the point atop the coordinate grid

        .. PLOT::

            M = Manifold(2, 'M', structure='topological')
            X = M.chart('x y'); x,y = X[:]
            p = M.point((1,3), name='p')
            g = p.plot(X)
            gX = X.plot(max_range=4)
            sphinx_plot(g+gX)

        Actually, since ``X`` is the default chart of the open set in which
        ``p`` has been defined, it can be skipped in the arguments of
        ``plot``::

            sage: g = p.plot()
            sage: show(g+gX)

        Call with some options::

            sage: g = p.plot(chart=X, size=40, color='green', label='$P$',
            ....:            label_color='blue', fontsize=20, label_offset=0.3)
            sage: show(g+gX)

        .. PLOT::

            M = Manifold(2, 'M', structure='topological')
            X = M.chart('x y'); x,y = X[:]
            p = M.point((1,3), name='p')
            g = p.plot(chart=X, size=40, color='green', label='$P$', \
                       label_color='blue', fontsize=20, label_offset=0.3)
            gX = X.plot(max_range=4)
            sphinx_plot(g+gX)

        Use of the ``parameters`` option to set a numerical value of some
        symbolic variable::

            sage: a = var('a')
            sage: q = M.point((a,2*a), name='q')
            sage: gq = q.plot(parameters={a:-2}, label_offset=0.2)
            sage: show(g+gX+gq)

        .. PLOT::

            M = Manifold(2, 'M', structure='topological')
            X = M.chart('x y'); x,y = X[:]
            p = M.point((1,3), name='p')
            g = p.plot(chart=X, size=40, color='green', label='$P$', \
                       label_color='blue', fontsize=20, label_offset=0.3)
            var('a')
            q = M.point((a,2*a), name='q')
            gq = q.plot(parameters={a:-2}, label_offset=0.2)
            gX = X.plot(max_range=4)
            sphinx_plot(g+gX+gq)

        The numerical value is used only for the plot::

            sage: q.coord()
            (a, 2*a)

        Drawing a point on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', structure='topological')
            sage: X.<x,y,z> = M.chart()
            sage: p = M.point((2,1,3), name='p')
            sage: g = p.plot()
            sage: print g
            Graphics3d Object
            sage: gX = X.plot(nb_values=5) # coordinate mesh cube
            sage: show(g+gX) # display of the point atop the coordinate mesh

        Call with some options::

            sage: g = p.plot(chart=X, size=40, color='green', label='P_1',
            ....:            label_color='blue', fontsize=20, label_offset=0.3)
            sage: show(g+gX)

        An example of plot via a mapping: plot of a point on a 2-sphere viewed
        in the 3-dimensional space ``M``::

            sage: S2 = Manifold(2, 'S^2', structure='topological')
            sage: U = S2.open_subset('U') # the open set covered by spherical coord.
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: p = U.point((pi/4, pi/8), name='p')
            sage: F = S2.continuous_map(M, {(XS, X): [sin(th)*cos(ph),
            ....:                           sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display()
            F: S^2 --> M
            on U: (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: g = p.plot(chart=X, mapping=F)
            sage: gS2 = XS.plot(chart=X, mapping=F, nb_values=9)
            sage: show(g+gS2)

        Use of the option ``ambient_coords`` for plots on a 4-dimensional
        manifold::

            sage: M = Manifold(4, 'M', structure='topological')
            sage: X.<t,x,y,z> = M.chart()
            sage: p = M.point((1,2,3,4), name='p')
            sage: g = p.plot(X, ambient_coords=(t,x,y), label_offset=0.4)  # the coordinate z is skipped
            sage: gX = X.plot(X, ambient_coords=(t,x,y), nb_values=5)
            sage: show(g+gX) # 3D plot
            sage: g = p.plot(X, ambient_coords=(t,y,z), label_offset=0.4)  # the coordinate x is skipped
            sage: gX = X.plot(X, ambient_coords=(t,y,z), nb_values=5)
            sage: show(g+gX) # 3D plot
            sage: g = p.plot(X, ambient_coords=(y,z), label_offset=0.4)  # the coordinates t and x are skipped
            sage: gX = X.plot(X, ambient_coords=(y,z))
            sage: show(g+gX) # 2D plot

        .. PLOT::

            M = Manifold(4, 'M', structure='topological')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            p = M.point((1,2,3,4), name='p')
            g = p.plot(X, ambient_coords=(y,z), label_offset=0.4)
            gX = X.plot(X, ambient_coords=(y,z))
            sphinx_plot(g+gX)

        """
        from sage.plot.point import point2d
        from sage.plot.text import text
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, text3d
        from sage.manifolds.chart import Chart
        if self._manifold.base_field_type() != 'real':
            raise NotImplementedError('plot of points on manifolds over ' +
                                  'fields different from R is not implemented')
        # The ambient chart:
        if chart is None:
            chart = self.parent().default_chart()
        elif not isinstance(chart, Chart):
            raise TypeError("the argument 'chart' must be a coordinate chart")
        # The effective point to be plotted:
        if mapping is None:
            eff_point = self
        else:
            eff_point = mapping(self)
        # The coordinates of the ambient chart used for the plot:
        if ambient_coords is None:
            ambient_coords = chart[:]
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca !=3:
            raise TypeError("bad number of ambient coordinates: {}".format(nca))
        # The point coordinates:
        coords = eff_point.coord(chart)
        xx = chart[:]
        xp = [coords[xx.index(c)] for c in ambient_coords]
        if parameters is not None:
            xps = [coord.substitute(parameters) for coord in xp]
            xp = xps
        xlab = [coord + label_offset for coord in xp]
        if label_color is None:
            label_color = color
        resu = Graphics()
        if nca == 2:
            if label is None:
                label = r'$' + self._latex_name + r'$'
            resu += point2d(xp, color=color, size=size) + \
                    text(label, xlab, fontsize=fontsize, color=label_color)
        else:
            if label is None:
                label = self._name
            resu += point3d(xp, color=color, size=size) + \
                    text3d(label, xlab, fontsize=fontsize, color=label_color)
        return resu
