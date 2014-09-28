r"""
Points on a manifold

The class :class:`Point` implements the concept of point on a manifold, in a 
coordinate independent manner: a :class:`Point` object can have coordinates in 
various charts defined on the manifold. Two points are declared equal if they 
have the same coordinates in the same chart. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013) : initial version

EXAMPLES: 

Defining a point on `\RR^3` by its spherical coordinates::

    sage: M = Manifold(3, 'R3', r'\mathcal{M}') 
    sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
    sage: p = M.point((1, pi/2, 0), name='P') # coordinates in the manifold's default chart
    sage: p
    point 'P' on 3-dimensional manifold 'R3'
    sage: latex(p) 
    P

Computing the coordinates of the point in a new chart::

    sage: c_cart.<x,y,z> = M.chart()        
    sage: ch = c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
    sage: p.coord(c_cart) # evaluate P's Cartesian coordinates
    (1, 0, 0)

Points can be compared::

    sage: p1 = M.point((1, pi/2, 0))
    sage: p == p1
    True
    sage: q = M.point((1,2,3), c_cart, name='Q') # point defined by its Cartesian coordinates
    sage: p == q
    False

Listing all the coordinates of a point in different charts::

    sage: p._coordinates # random (dictionary output)
    {chart (R3, (r, th, ph)): (1, 1/2*pi, 0), chart (R3, (x, y, z)): (1, 0, 0)}
    sage: q._coordinates
    {chart (R3, (x, y, z)): (1, 2, 3)}

"""

#*****************************************************************************
#       Copyright (C) 2013 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element import Element   

class Point(Element):
    r"""
    Class for points on a manifold.

    INPUT:
    
    - ``domain`` -- the manifold domain to which the point belongs (can be 
      the entire manifold)
    - ``coords`` -- (default: None) the point coordinates (as a tuple or a list)
    - ``chart`` -- (default: None) chart in which the coordinates are given; 
      if none is provided, the coordinates are assumed 
      to refer to the domain's default chart
    - ``name`` -- (default: None) name given to the point
    - ``latex_name`` -- (default: None) LaTeX symbol to denote the point; if 
      none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:
    
    A point on a 2-dimensional manifold::
    
        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: (a, b) = var('a b') # generic coordinates for the point
        sage: p = M.point((a, b), name='P') ; p
        point 'P' on 2-dimensional manifold 'M'
        sage: p.coord()  # coordinates of P in the domain's default chart
        (a, b)

    Since points are Sage 'Element', the 'Parent' of which being the domain
    on which they are defined, it is equivalent to write::
    
        sage: p = M((a, b), name='P') ; p
        point 'P' on 2-dimensional manifold 'M'
    
    A point is an element of the manifold domain on which it has been defined::
    
        sage: p in M
        True
        sage: p.parent()
        2-dimensional manifold 'M'
        
    By default, the LaTeX symbol of the point is deduced from its name::
    
        sage: latex(p)
        P
        
    But it can be set to any value::
    
        sage: p = M.point((a, b), name='P', latex_name=r'\mathcal{P}')
        sage: latex(p)
        \mathcal{P}
    
    Points can be drawn in 2D or 3D graphics thanks to the method :meth:`plot`.
    
    """
    def __init__(self, domain, coords=None, chart=None, name=None, 
                 latex_name=None): 
        Element.__init__(self, domain)
        self._manifold = domain._manifold
        self._domain = domain
        self._coordinates = {}
        if coords is not None: 
            if len(coords) != self._manifold._dim: 
                raise ValueError("The number of coordinates must be equal" +
                                 " to the manifold dimension.")
            if chart is None: 
                chart = self._domain._def_chart
            else: 
                if chart not in self._domain._atlas: 
                    raise ValueError("The " + str(chart) +
                            " has not been defined on the " + str(self._domain))
            #!# The following check is not performed for it would fail with 
            # symbolic coordinates:
            # if not chart.valid_coordinates(*coords):
            #    raise ValueError("The coordinates " + str(coords) + 
            #                     " are not valid on the " + str(chart))
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
        self._tangent_space = None # tangent vector space at the point (not
                                   # constructed yet)
        self._frame_bases = {} # dictionary of bases of the tangent vector 
                               # derived from vector frames (keys: vector 
                               # frames)

    def _repr_(self):
        r"""
        Special Sage function for the string representation of the object.
        """
        description = "point"
        if self._name is not None:
            description += " '%s'" % self._name
        description += " on " + str(self._manifold)
        return description

    def _latex_(self):
        r"""
        Special Sage function for the LaTeX representation of the object.
        """
        if self._latex_name is None:
            return r'\mbox{' + str(self) + r'}'
        else:
           return self._latex_name

    def coord(self, chart=None, old_chart=None):
        r"""
        Return the point coordinates in the specified chart.

        If these coordinates are not already known, they are computed from 
        known ones by means of change-of-chart formulas. 

        INPUT:
    
        - ``chart`` -- (default: None) chart in which the coordinates are 
          given; if none is provided, the coordinates are assumed to refer to
          the domain's default chart
        - ``old_chart`` -- (default: None) chart from which the coordinates in 
          ``chart`` are to be computed. If None, a chart in which the point's 
          coordinates are already known will be picked, priveleging the 
          domain's default chart.

        EXAMPLES: 

        Spherical coordinates of a point on `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi') # spherical coordinates
            sage: p = M.point((1, pi/2, 0)) 
            sage: p.coord()    # coordinates on the manifold's default chart
            (1, 1/2*pi, 0)
            sage: p.coord(c_spher) # with the chart c_spher specified (same result as above since this is the default chart)
            (1, 1/2*pi, 0)

        Computing the Cartesian coordinates from the spherical ones::

            sage: c_cart.<x,y,z> = M.chart()  # Cartesian coordinates   
            sage: c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            coordinate change from chart (R3, (r, th, ph)) to chart (R3, (x, y, z))
            sage: p.coord(c_cart)  # the computation is performed by means of the above change of coordinates
            (1, 0, 0)

        Coordinates of a point on a 2-dimensional manifold::
    
            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: (a, b) = var('a b') # generic coordinates for the point
            sage: p = M.point((a, b), name='P')
            sage: p.coord()  # coordinates of P in the manifold's default chart
            (a, b)
            
        Coordinates of P in a new chart::
        
            sage: c_uv.<u,v> = M.chart()
            sage: ch_xy_uv = c_xy.coord_change(c_uv, x-y, x+y)
            sage: p.coord(c_uv)
            (a - b, a + b)

        Coordinates of P in a third chart::
        
            sage: c_wz.<w,z> = M.chart()
            sage: ch_uv_wz = c_uv.coord_change(c_wz, u^3, v^3)   
            sage: p.coord(c_wz, old_chart=c_uv)
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)

        Actually, in the present case, it is not necessary to specify 
        old_chart='uv'::
        
            sage: p.set_coord((a-b, a+b), c_uv) # erases all the coordinates except those in the chart c_uv
            sage: p._coordinates  
            {chart (M, (u, v)): (a - b, a + b)}              
            sage: p.coord(c_wz)                
            (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)
            sage: p._coordinates # random (dictionary output)
            {chart (M, (u, v)): (a - b, a + b), chart (M, (w, z)): (a^3 - 3*a^2*b + 3*a*b^2 - b^3, a^3 + 3*a^2*b + 3*a*b^2 + b^3)}

        """
        if chart is None:
            dom = self._domain 
            chart = dom._def_chart
            def_chart = chart
        else:
            dom = chart._domain
            def_chart = dom._def_chart
            if self not in dom:
                raise ValueError("The point does not belong to the domain " + 
                                 "of " + str(chart))
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
                raise ValueError("The coordinates of " + str(self) + \
                    " in the " + str(chart) + " cannot be computed" + \
                    " by means of known changes of charts.")
            else:
                chcoord = dom._coord_changes[(s_old_chart, s_chart)]
                self._coordinates[chart] = \
                                    chcoord(*self._coordinates[old_chart])
        return self._coordinates[chart]
        
    def set_coord(self, coords, chart=None):
        r"""
        Sets the point coordinates in the specified chart.
        
        Coordinates with respect to other charts are deleted, in order to 
        avoid any inconsistency. To keep them, use the method :meth:`add_coord` 
        instead.

        INPUT:
        
        - ``coords`` -- the point coordinates (as a tuple or a list)
        - ``chart`` -- (default: None) chart in which the coordinates are 
          given; if none is provided, the coordinates are assumed to refer to 
          the domain's default chart
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart()
            sage: p = M.point()
            sage: p.set_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: q = M.point()
            sage: q.set_coord((1,2,3), c_spher)
            sage: cart_from_spher = c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
        
        If we set the coordinates of q in the chart c_cart, those in the chart c_spher
        are lost::
        
            sage: q.set_coord( cart_from_spher(*q.coord(c_spher)), c_cart)
            sage: q._coordinates
            {chart (R3, (x, y, z)): (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p._coordinates
            {chart (R3, (x, y, z)): (1, 2, 3)}
            
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
        - ``chart`` -- (default: None) chart in which the coordinates are 
          given; if none is provided, the coordinates are assumed to refer to 
          the domain's default chart
          
        .. WARNING::
        
           If the point has already coordinates in other charts, it 
           is the user's responsability to make sure that the coordinates
           to be added are consistent with them. 
    
        EXAMPLES: 

        Setting coordinates to a point on `\RR^3`::

            sage: Manifold._clear_cache_() # for doctests only
            sage: M = Manifold(3, 'R3', r'\mathcal{M}')
            sage: c_cart.<x,y,z> = M.chart()
            sage: p = M.point()
            sage: p.add_coord((1,2,3))  # coordinates on the manifold's default chart
            sage: p.coord()
            (1, 2, 3)
            
        A point defined in another coordinate system::
        
            sage: c_spher.<r,th,ph> = M.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: q = M.point()
            sage: q.add_coord((1,2,3), c_spher)
            sage: cart_from_spher = c_spher.coord_change(c_cart, r*sin(th)*cos(ph), r*sin(th)*sin(ph), r*cos(th))
            sage: q.add_coord( cart_from_spher(*q.coord(c_spher)), c_cart)
            sage: q._coordinates # random (dictionary output)
            {chart (R3, (r, th, ph)): (1, 2, 3), chart (R3, (x, y, z)): (cos(3)*sin(2), sin(3)*sin(2), cos(2))}
            sage: p._coordinates
            {chart (R3, (x, y, z)): (1, 2, 3)}
            sage: p == q  # p and q should differ because the coordinates (1,2,3) are on different charts
            False     
            
        Contrary to :meth:`set_coord`, the method :meth:`add_coord` does not 
        the coordinates in other charts::
        
            sage: p = M.point((1,2,3), c_spher)
            sage: p._coordinates
            {chart (R3, (r, th, ph)): (1, 2, 3)}
            sage: p.set_coord((4,5,6), c_cart)
            sage: p._coordinates
            {chart (R3, (x, y, z)): (4, 5, 6)}
            sage: p.add_coord((7,8,9), c_spher)
            sage: p._coordinates # random (dictionary output) 
            {chart (R3, (x, y, z)): (4, 5, 6), chart (R3, (r, th, ph)): (7, 8, 9)}
            
        """
        if len(coords) != self._manifold._dim: 
            raise ValueError("The number of coordinates must be equal " + 
                             "to the manifold dimension.")
        if chart is None: 
            chart = self._domain._def_chart
        else: 
            if chart not in self._domain._atlas:
                raise ValueError("The " + str(chart) +
                    " has not been defined on the " + str(self._domain))
        self._coordinates[chart] = coords

    def __eq__(self, other):
        r"""
        Compares the current point with another one.
        """
        if not isinstance(other, Point):
            return False
        # Search for a common chart to compare the coordinates
        common_chart = None
        # the domain's default chart is privileged:
        def_chart = self._domain._def_chart
        if def_chart in self._coordinates and def_chart in other._coordinates:
            common_chart = def_chart
        else:
            for chart in self._coordinates:
                if chart in other._coordinates:
                    common_chart = chart
                    break
        if common_chart is None:
            raise ValueError("No common chart has been found to compare " +
                             str(self) + " and " + str(other))
        return self._coordinates[common_chart] == \
                                              other._coordinates[common_chart]
        
    def tangent_space(self): 
        r"""
        Returns the tangent space at self.
        """
        from tangentspace import TangentSpace
        if self._tangent_space is not None:
            return self._tangent_space
        else: 
            res = TangentSpace(self)
            (self._manifold)._tangent_spaces[self] = res
            return res 

    def plot(self, ambient_chart, ambient_coords=None, size=10, color='black', 
             label=None, label_color=None, fontsize=10, label_offset=0.1, 
             parameters=None):
        r"""
        Plot the point in a Cartesian graph based on the coordinates of
        some chart, called hereafter the *ambient chart*.

        INPUT:
        
        - ``ambient_chart`` -- the ambient chart (see above)
        - ``ambient_coords`` -- (default: None) tuple containing the 2 or 3 
          coordinates of the ambient chart in terms of which the plot is 
          performed; if None, all the coordinates of the ambient chart are 
          considered
        - ``size`` -- (default: 10) size of the point once drawn as a small
          disk or sphere
        - ``color`` -- (default: 'black') color of the point
        - ``label`` -- (default: None) label printed next to the point; if None, 
          the point's name is used.  
        - ``label_color`` -- (default: None) color to print the label; if None,
          the value of ``color`` is used
        - ``fontsize`` -- (default: 10) size of the font used to print the 
          label
        - ``label_offset`` -- (default: 0.1) determines the separation between
          the point and its label
        - ``parameters`` -- (default: None) dictionary giving the numerical
          values of the parameters that may appear in the point coordinates

        OUTPUT:
        
        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of the ambient chart) or an instance of 
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e. 
          based on 3 coordinates of the ambient chart)
        
        EXAMPLES:
        
        Drawing a point on a 2-dimensional manifold::
        
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,3), name='p')
            sage: g = p.plot(X)
            sage: print g
            Graphics object consisting of 2 graphics primitives
            sage: gX = X.plot(X) # plot of the coordinate grid
            sage: show(g+gX) # display of the point atop the coordinate grid
            
        Call with some options:: 
        
            sage: g = p.plot(X, size = 40, color='green', label='$P$', label_color='blue', fontsize=20, label_offset=0.3)
            sage: show(g+gX)
            
        Use of the ``parameters`` option to set a numerical value of some 
        symbolic variable::
        
            sage: a = var('a')
            sage: q = M.point((a,2*a), name='q')
            sage: gq = q.plot(X, parameters={a:-2})
            sage: show(g+gX+gq)
        
        The numerical value is used only for the plot::

            sage: q.coord()
            (a, 2*a)

        Drawing a point on a 3-dimensional manifold::

            sage: M = Manifold(3, 'M')
            sage: X.<x,y,z> = M.chart()
            sage: p = M.point((2,1,3), name='p')
            sage: g = p.plot(X)
            sage: print g
            Graphics3d Object
            sage: gX = X.plot(X, nb_values=5) # coordinate mesh cube
            sage: show(g+gX) # display of the point atop the coordinate mesh
            
        Call with some options:: 
        
            sage: g = p.plot(X, size = 40, color='green', label='P_1', label_color='blue', fontsize=20, label_offset=0.3)
            sage: show(g+gX)

        Use of the option ``ambient_coords`` for plots on a 4-dimensional 
        manifold::
        
            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: p = M.point((1,2,3,4), name='p')
            sage: g = p.plot(X, ambient_coords=(t,x,y))  # the coordinate z is skipped
            sage: gX = X.plot(X, ambient_coords=(t,x,y), nb_values=5)
            sage: show(g+gX) # 3D plot
            sage: g = p.plot(X, ambient_coords=(t,y,z))  # the coordinate x is skipped
            sage: gX = X.plot(X, ambient_coords=(t,y,z), nb_values=5)
            sage: show(g+gX) # 3D plot
            sage: g = p.plot(X, ambient_coords=(y,z))  # the coordinates t and x are skipped
            sage: gX = X.plot(X, ambient_coords=(y,z))
            sage: show(g+gX) # 2D plot

        """
        from sage.plot.point import point2d
        from sage.plot.text import text
        from sage.plot.graphics import Graphics
        from sage.plot.plot3d.shapes2 import point3d, text3d
        if ambient_coords is None:
            ambient_coords = ambient_chart._xx
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca !=3:
            raise TypeError("Bad number of ambient coordinates: " + str(nca))
        coords = self.coord(ambient_chart)
        xx = ambient_chart[:]
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
