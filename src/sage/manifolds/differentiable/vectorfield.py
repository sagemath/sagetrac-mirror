r"""
Vector fields

Given two differentiable manifolds `U` and `M` over the same topological field
`K` and a differentiable map

.. MATH::

    \Phi:\ U \longrightarrow  M

we define a *vector field along* `U` *with values on* `M` to be a
differentiable map

.. MATH::

    v:\ U  \longrightarrow TM

(`TM` being the tangent bundle of `M`) such that

.. MATH::

    \forall p \in U,\ v(p) \in T_{\Phi(p)}M

The standard case of vector fields *on* a differentiable manifold corresponds
to `U=M` and `\Phi = \mathrm{Id}_M`. Other common cases are `\Phi`
being an immersion and `\Phi` being a curve in `M` (`U` is then an open
interval of `\RR`).

Vector fields are implemented via two classes: :class:`VectorFieldParal` and
:class:`VectorField`, depending respectively whether the manifold `M`
is parallelizable or not, i.e. whether the bundle `TM` is trivial or not.


AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013); :doi:`10.1007/978-1-4419-9982-5`
- B O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego) (1983)
- M. Berger & B. Gostiaux : *Differential Geometry: Manifolds, Curves and
  Surfaces*, Springer (New York) (1988); :doi:`10.1007/978-1-4612-1033-7`

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_tensor import FiniteRankFreeModuleElement
from sage.manifolds.differentiable.tensorfield import TensorField
from sage.manifolds.differentiable.tensorfield_paral import TensorFieldParal

class VectorField(TensorField):
    r"""
    Vector field along a differentiable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a differentiable manifold `M`, via a
    differentiable map `U \rightarrow M`. More precisely, given a
    differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U=M` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    If `M` is parallelizable, the class
    :class:`~sage.manifolds.differentiable.vectorfield.VectorFieldParal`
    must be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldModule`.

    INPUT:

    - ``vector_field_module`` -- module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a non-parallelizable 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: U = M.open_subset('U') ; V = M.open_subset('V')
        sage: M.declare_union(U,V)   # M is the union of U and V
        sage: c_xy.<x,y> = U.chart() ; c_tu.<t,u> = V.chart()
        sage: transf = c_xy.transition_map(c_tu, (x+y, x-y), intersection_name='W',
        ....:                              restrictions1= x>0, restrictions2= t+u>0)
        sage: inv = transf.inverse()
        sage: W = U.intersection(V)
        sage: eU = c_xy.frame() ; eV = c_tu.frame()
        sage: c_tuW = c_tu.restrict(W) ; eVW = c_tuW.frame()
        sage: v = M.vector_field('v') ; v
        Vector field v on the 2-dimensional differentiable manifold M
        sage: v.parent()
        Module X(M) of vector fields on the 2-dimensional differentiable
         manifold M

    The vector field is first defined on the domain `U` by means of its
    components w.r.t. the frame eU::

        sage: v[eU,:] = [-y, 1+x]

    The components w.r.t the frame eV are then deduced by continuation of the
    components w.r.t. the frame eVW on the domain `W=U\cap V`, expressed in
    terms on the coordinates covering `V`::

        sage: v[eV,0] = v[eVW,0,c_tuW].expr()
        sage: v[eV,1] = v[eVW,1,c_tuW].expr()

    At this stage, the vector field is fully defined on the whole manifold::

        sage: v.display(eU)
        v = -y d/dx + (x + 1) d/dy
        sage: v.display(eV)
        v = (u + 1) d/dt + (-t - 1) d/du

    The vector field acting on scalar fields::

        sage: f = M.scalar_field({c_xy: (x+y)^2, c_tu: t^2}, name='f')
        sage: s = v(f) ; s
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: s.display()
        v(f): M --> R
        on U: (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
        on V: (t, u) |--> 2*t*u + 2*t

    Some checks::

        sage: v(f) == f.differential()(v)
        True
        sage: v(f) == f.lie_der(v)
        True

    The result is defined on the intersection of the vector field's domain and
    the scalar field's one::

        sage: s = v(f.restrict(U)) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s == v(f).restrict(U)
        True
        sage: s = v(f.restrict(W)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) |--> 2*t*u + 2*t
        sage: s = v.restrict(U)(f) ; s
        Scalar field v(f) on the Open subset U of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): U --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
        on W: (t, u) |--> 2*t*u + 2*t
        sage: s = v.restrict(U)(f.restrict(V)) ; s
        Scalar field v(f) on the Open subset W of the 2-dimensional
         differentiable manifold M
        sage: s.display()
        v(f): W --> R
           (x, y) |--> 2*x^2 - 2*y^2 + 2*x + 2*y
           (t, u) |--> 2*t*u + 2*t

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a non-parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorField``, to fit with the category framework::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: XM = M.vector_field_module()
            sage: a = XM.element_class(XM, name='a'); a
            Vector field a on the 2-dimensional differentiable manifold M
            sage: a[c_xy.frame(),:] = [x, y]
            sage: a[c_uv.frame(),:] = [-u, -v]
            sage: TestSuite(a).run(skip='_test_pickling')

        Construction with ``DifferentiableManifold.vector_field``::

            sage: a1 = M.vector_field(name='a'); a1
            Vector field a on the 2-dimensional differentiable manifold M
            sage: type(a1) == type(a)
            True

        .. TODO::

            fix _test_pickling (in the superclass TensorField)

        """
        TensorField.__init__(self, vector_field_module, (1,0), name=name,
                             latex_name=latex_name)
        # Initialization of derived quantities:
        TensorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self) :
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M

        """
        description = "Vector field "
        if self._name is not None:
            description += self._name + " "
        return self._final_repr(description)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True

        """
        return type(self)(self._vmodule)

    def _init_dependencies(self):
        r"""
        Initialize list of quantities that depend on ``self``

        TEST::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._init_dependencies()

        """
        self._lie_der_along_self = {}

    def _del_dependencies(self):
        r"""
        Clear list of quantities that depend on ``self``

        TEST::

            sage: M = Manifold(2, 'M')
            sage: v = M.vector_field(name='v')
            sage: v._del_dependencies()

        """
        if self._lie_der_along_self != {}:
            for idtens, tens in self._lie_der_along_self.iteritems():
                del tens._lie_derivatives[id(self)]
            self._lie_der_along_self.clear()

    def __call__(self, scalar):
        r"""
        Action on a scalar field (or on a 1-form)

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        TESTS::

            sage: M = Manifold(2, 'M') # the 2-dimensional sphere S^2
            sage: U = M.open_subset('U') # complement of the North pole
            sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            sage: V = M.open_subset('V') # complement of the South pole
            sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            sage: M.declare_union(U,V)   # S^2 is the union of U and V
            sage: a = M.vector_field(name='a')
            sage: a[c_xy.frame(),:] = [x, y]
            sage: a[c_uv.frame(),:] = [-u, -v]
            sage: f = M.scalar_field({c_xy: atan(x^2+y^2), c_uv: pi/2-atan(u^2+v^2)},
            ....:                    name='f')
            sage: s = a.__call__(f); s
            Scalar field a(f) on the 2-dimensional differentiable manifold M
            sage: s.display()
            a(f): M --> R
            on U: (x, y) |--> 2*(x^2 + y^2)/(x^4 + 2*x^2*y^2 + y^4 + 1)
            on V: (u, v) |--> 2*(u^2 + v^2)/(u^4 + 2*u^2*v^2 + v^4 + 1)
            sage: s == f.differential()(a)
            True

        """
        if scalar._tensor_type == (0,1):
            # This is actually the action of the vector field on a 1-form,
            # as a tensor field of type (1,0):
            return scalar(self)
        if scalar._tensor_type != (0,0):
            raise TypeError("the argument must be a scalar field")
        #!# Could it be simply
        # return scalar.differential()(self)
        # ?
        dom_resu = self._domain.intersection(scalar._domain)
        self_r = self.restrict(dom_resu)
        scalar_r = scalar.restrict(dom_resu)
        if scalar_r._is_zero:
            return dom_resu._zero_scalar_field
        if isinstance(self_r, VectorFieldParal):
            return self_r(scalar_r)
        # Creation of the result:
        if self._name is not None and scalar._name is not None:
            resu_name = self._name + "(" + scalar._name + ")"
        else:
            resu_name = None
        if self._latex_name is not None and scalar._latex_name is not None:
            resu_latex = self._latex_name + r"\left(" + scalar._latex_name + \
                        r"\right)"
        else:
            resu_latex = None
        resu = dom_resu.scalar_field(name=resu_name, latex_name=resu_latex)
        for dom, rst in self_r._restrictions.iteritems():
            resu_rst = rst(scalar_r.restrict(dom))
            for chart, funct in resu_rst._express.iteritems():
                resu._express[chart] = funct
        return resu

    def plot(self, chart=None, ambient_coords=None, mapping=None,
             chart_domain=None, fixed_coords=None, ranges=None, max_range=8,
             nb_values=None, steps=None,scale=1, color='blue', parameters=None,
             label_axes=True, **extra_options):
        r"""
        Plot the vector field in a Cartesian graph based on the coordinates
        of some ambient chart.

        The vector field is drawn in terms of two (2D graphics) or three
        (3D graphics) coordinates of a given chart, called hereafter the
        *ambient chart*.
        The vector field's base points `p` (or their images `\Phi(p)` by some
        differentiable mapping `\Phi`) must lie in the ambient chart's domain.

        INPUT:

        - ``chart`` -- (default: ``None``) the ambient chart (see above); if
          ``None``, the default chart of the vector field's domain is used
        - ``ambient_coords`` -- (default: ``None``) tuple containing the 2 or 3
          coordinates of the ambient chart in terms of which the plot is
          performed; if ``None``, all the coordinates of the ambient chart are
          considered
        - ``mapping`` -- (default: ``None``) differentiable map `\Phi`
          (instance of
          :class:`~sage.manifolds.differentiable.diff_map.DiffMap`)
          providing the link between the vector field's domain and
          the ambient chart ``chart``; if ``None``, the identity map is
          assumed
        - ``chart_domain`` -- (default: ``None``) chart on the vector
          field's domain to define the points at which vector arrows are to be
          plotted; if ``None``, the default chart of the vector field's domain
          is used
        - ``fixed_coords`` -- (default: ``None``) dictionary with keys the
          coordinates of ``chart_domain`` that are kept fixed and with values
          the value of these coordinates; if ``None``, all the coordinates of
          ``chart_domain`` are used
        - ``ranges`` -- (default: ``None``) dictionary with keys the
          coordinates of ``chart_domain`` to be used and values
          tuples ``(x_min,x_max)`` specifying the
          coordinate range for the plot; if ``None``, the entire coordinate
          range declared during the construction of ``chart_domain`` is
          considered (with ``-Infinity`` replaced by ``-max_range`` and
          ``+Infinity`` by ``max_range``)
        - ``max_range`` -- (default: 8) numerical value substituted to
          ``+Infinity`` if the latter is the upper bound of the range of a
          coordinate for which the plot is performed over the entire coordinate
          range (i.e. for which no specific plot range has been set in
          ``ranges``); similarly ``-max_range`` is the numerical valued
          substituted for ``-Infinity``
        - ``nb_values`` -- (default: ``None``) either an integer or a dictionary
          with keys the coordinates of ``chart_domain`` to be used and values
          the number of values of the coordinate for sampling
          the part of the vector field's domain involved in the plot ; if
          ``nb_values`` is a single integer, it represents the number of
          values for all coordinates; if ``nb_values`` is ``None``, it is set
          to 9 for a 2D plot and to 5 for a 3D plot
        - ``steps`` -- (default: ``None``) dictionary with keys the coordinates
          of ``chart_domain`` to be used and values the step between each
          constant value of the coordinate; if ``None``, the step is computed
          from the coordinate range (specified in ``ranges``) and ``nb_values``.
          On the contrary, if the step is provided for some coordinate, the
          corresponding number of values is deduced from it and the coordinate
          range.
        - ``scale`` -- (default: 1) value by which the lengths of the arrows
          representing the vectors is multiplied
        - ``color`` -- (default: 'blue') color of the arrows representing the
          vectors
        - ``parameters`` -- (default: ``None``) dictionary giving the numerical
          values of the parameters that may appear in the coordinate expression
          of the vector field (see example below)
        - ``label_axes`` -- (default: ``True``) boolean determining whether the
          labels of the coordinate axes of ``chart`` shall be added to the
          graph; can be set to ``False`` if the graph is 3D and must be
          superposed with another graph.
        - ``**extra_options`` -- extra options for the arrow plot, like
          ``linestyle``, ``width`` or ``arrowsize`` (see
          :func:`~sage.plot.arrow.arrow2d` and
          :func:`~sage.plot.plot3d.shapes.arrow3d` for details)

        OUTPUT:

        - a graphic object, either an instance of
          :class:`~sage.plot.graphics.Graphics` for a 2D plot (i.e. based on
          2 coordinates of ``chart``) or an instance of
          :class:`~sage.plot.plot3d.base.Graphics3d` for a 3D plot (i.e.
          based on 3 coordinates of ``chart``)

        EXAMPLES:

        Plot of a vector field on a 2-dimensional manifold::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: v = M.vector_field(name='v')
            sage: v[:] = -y, x ; v.display()
            v = -y d/dx + x d/dy
            sage: v.plot()
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(name='v'); v[:] = -y, x
            g = v.plot()
            sphinx_plot(g)

        Plot with various options::

            sage: v.plot(scale=0.5, color='green', linestyle='--', width=1,
            ....:        arrowsize=6)
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(name='v'); v[:] = -y, x
            g = v.plot(scale=0.5, color='green', linestyle='--', width=1, arrowsize=6)
            sphinx_plot(g)

        ::

            sage: v.plot(max_range=4, nb_values=5, scale=0.5)
            Graphics object consisting of 24 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(name='v'); v[:] = -y, x
            g = v.plot(max_range=4, nb_values=5, scale=0.5)
            sphinx_plot(g)

        Plots along a line of fixed coordinate::

            sage: v.plot(fixed_coords={x: -2})
            Graphics object consisting of 9 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(name='v'); v[:] = -y, x
            g = v.plot(fixed_coords={x: -2})
            sphinx_plot(g)

        ::

            sage: v.plot(fixed_coords={y: 1})
            Graphics object consisting of 9 graphics primitives

        .. PLOT::

            M = Manifold(2, 'M')
            X = M.chart('x y'); x, y = X[:]
            v = M.vector_field(name='v'); v[:] = -y, x
            g = v.plot(fixed_coords={y: 1})
            sphinx_plot(g)

        Let us now consider a vector field on a 4-dimensional manifold::

            sage: M = Manifold(4, 'M')
            sage: X.<t,x,y,z> = M.chart()
            sage: v = M.vector_field(name='v')
            sage: v[:] = (t/8)^2, -t*y/4, t*x/4, t*z/4 ; v.display()
            v = 1/64*t^2 d/dt - 1/4*t*y d/dx + 1/4*t*x d/dy + 1/4*t*z d/dz

        We cannot make a 4D plot directly::

            sage: v.plot()
            Traceback (most recent call last):
            ...
            ValueError: the number of ambient coordinates must be either 2 or 3, not 4

        Rather, we have to select some coordinates for the plot, via
        the argument ``ambient_coords``. For instance, for a 3D plot::

            sage: g = v.plot(ambient_coords=(x, y, z), fixed_coords={t: 1})
            sage: print(g)
            Graphics3d Object
            sage: g = v.plot(ambient_coords=(x, y, t), fixed_coords={z: 0},
            ....:            ranges={x: (-2,2), y: (-2,2), t: (-1, 4)},
            ....:            nb_values=4)
            sage: print(g)
            Graphics3d Object

        or, for a 2D plot::

            sage: v.plot(ambient_coords=(x, y), fixed_coords={t: 1, z: 0})
            Graphics object consisting of 80 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            v = M.vector_field(name='v')
            v[:] = v[:] = (t/8)**2, -t*y/4, t*x/4, t*z/4
            g = v.plot(ambient_coords=(x, y), fixed_coords={t: 1, z: 0})
            sphinx_plot(g)

        ::

            sage: v.plot(ambient_coords=(x, t), fixed_coords={y: 1, z: 0})
            Graphics object consisting of 72 graphics primitives

        .. PLOT::

            M = Manifold(4, 'M')
            X = M.chart('t x y z'); t,x,y,z = X[:]
            v = M.vector_field(name='v')
            v[:] = v[:] = (t/8)**2, -t*y/4, t*x/4, t*z/4
            g = v.plot(ambient_coords=(x, t), fixed_coords={y: 1, z: 0})
            sphinx_plot(g)

        An example of plot via a differential mapping: plot of a vector field
        tangent to a 2-sphere viewed in `\RR^3`::

            sage: S2 = Manifold(2, 'S^2')
            sage: U = S2.open_subset('U') # the open set covered by spherical coord.
            sage: XS.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            sage: R3 = Manifold(3, 'R^3')
            sage: X3.<x,y,z> = R3.chart()
            sage: F = S2.diff_map(R3, {(XS, X3): [sin(th)*cos(ph),
            ....:                       sin(th)*sin(ph), cos(th)]}, name='F')
            sage: F.display() # the standard embedding of S^2 into R^3
            F: S^2 --> R^3
            on U: (th, ph) |--> (x, y, z) = (cos(ph)*sin(th), sin(ph)*sin(th), cos(th))
            sage: v = XS.frame()[1] ; v
            Vector field d/dph on the Open subset U of the 2-dimensional
             differentiable manifold S^2
            sage: graph_v = v.plot(chart=X3, mapping=F, label_axes=False)
            sage: graph_S2 = XS.plot(chart=X3, mapping=F, nb_values=9)
            sage: show(graph_v + graph_S2)

        """
        from sage.rings.infinity import Infinity
        from sage.misc.functional import numerical_approx
        from sage.misc.latex import latex
        from sage.plot.graphics import Graphics
        from sage.manifolds.chart import RealChart
        from sage.manifolds.utilities import set_axes_labels
        #
        # 1/ Treatment of input parameters
        #    -----------------------------
        if chart is None:
            chart = self._domain.default_chart()
        elif not isinstance(chart, RealChart):
            raise TypeError("{} is not a chart on a real ".format(chart) +
                            "manifold")
        if chart_domain is None:
            chart_domain = self._domain.default_chart()
        elif not isinstance(chart_domain, RealChart):
            raise TypeError("{} is not a chart on a ".format(chart_domain) +
                            "real manifold")
        elif not chart_domain.domain().is_subset(self._domain):
            raise ValueError("the domain of {} is not ".format(chart_domain) +
                             "included in the domain of {}".format(self))
        if fixed_coords is None:
            coords = chart_domain._xx
        else:
            fixed_coord_list = fixed_coords.keys()
            coords = []
            for coord in chart_domain._xx:
                if coord not in fixed_coord_list:
                    coords.append(coord)
            coords = tuple(coords)
        if ambient_coords is None:
            ambient_coords = chart[:]
        elif not isinstance(ambient_coords, tuple):
            ambient_coords = tuple(ambient_coords)
        nca = len(ambient_coords)
        if nca != 2 and nca !=3:
            raise ValueError("the number of ambient coordinates must be " +
                             "either 2 or 3, not {}".format(nca))
        if ranges is None:
            ranges = {}
        ranges0 = {}
        for coord in coords:
            if coord in ranges:
                ranges0[coord] = (numerical_approx(ranges[coord][0]),
                                  numerical_approx(ranges[coord][1]))
            else:
                bounds = chart_domain._bounds[chart_domain[:].index(coord)]
                xmin0 = bounds[0][0]
                xmax0 = bounds[1][0]
                if xmin0 == -Infinity:
                    xmin = numerical_approx(-max_range)
                elif bounds[0][1]:
                    xmin = numerical_approx(xmin0)
                else:
                    xmin = numerical_approx(xmin0 + 1.e-3)
                if xmax0 == Infinity:
                    xmax = numerical_approx(max_range)
                elif bounds[1][1]:
                    xmax = numerical_approx(xmax0)
                else:
                    xmax = numerical_approx(xmax0 - 1.e-3)
                ranges0[coord] = (xmin, xmax)
        ranges = ranges0
        if nb_values is None:
            if nca == 2: # 2D plot
                nb_values = 9
            else:   # 3D plot
                nb_values = 5
        if not isinstance(nb_values, dict):
            nb_values0 = {}
            for coord in coords:
                nb_values0[coord] = nb_values
            nb_values = nb_values0
        if steps is None:
            steps = {}
        for coord in coords:
            if coord not in steps:
                steps[coord] = (ranges[coord][1] - ranges[coord][0])/ \
                               (nb_values[coord]-1)
            else:
                nb_values[coord] = 1 + int(
                           (ranges[coord][1] - ranges[coord][0])/ steps[coord])
        #
        # 2/ Plots
        #    -----
        dom = chart_domain.domain()
        vector = self.restrict(dom)
        if vector.parent().destination_map() is dom.identity_map():
            if mapping is not None:
                vector = mapping.pushforward(vector)
                mapping = None
        nc = len(chart_domain[:])
        ncp = len(coords)
        xx = [0] * nc
        if fixed_coords is not None:
            if len(fixed_coords) != nc - ncp:
                raise ValueError("bad number of fixed coordinates")
            for fc, val in fixed_coords.iteritems():
                xx[chart_domain[:].index(fc)] = val
        index_p = [chart_domain[:].index(cd) for cd in coords]
        resu = Graphics()
        ind = [0] * ncp
        ind_max = [0] * ncp
        ind_max[0] = nb_values[coords[0]]
        xmin = [ranges[cd][0] for cd in coords]
        step_tab = [steps[cd] for cd in coords]
        while ind != ind_max:
            for i in range(ncp):
                xx[index_p[i]] = xmin[i] + ind[i]*step_tab[i]
            if chart_domain.valid_coordinates(*xx, tolerance=1e-13,
                                              parameters=parameters):
                point = dom(xx, chart=chart_domain)
                resu += vector.at(point).plot(chart=chart,
                                              ambient_coords=ambient_coords,
                                              mapping=mapping, scale=scale,
                                              color=color, print_label=False,
                                              parameters=parameters,
                                              **extra_options)
            # Next index:
            ret = 1
            for pos in range(ncp-1,-1,-1):
                imax = nb_values[coords[pos]] - 1
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = 0
                        ret = 1
        if label_axes:
            if nca==2:  # 2D graphic
                # We update the dictionary _extra_kwds (options to be passed
                # to show()), instead of using the method
                # Graphics.axes_labels() since the latter is not robust w.r.t.
                # graph addition
                resu._extra_kwds['axes_labels'] = [r'$'+latex(ac)+r'$'
                                                   for ac in ambient_coords]
            else: # 3D graphic
                labels = [str(ac) for ac in ambient_coords]
                resu = set_axes_labels(resu, *labels)
        return resu

#******************************************************************************

class VectorFieldParal(FiniteRankFreeModuleElement, TensorFieldParal,
                       VectorField):

    r"""
    Vector field along a differentiable manifold, with values on a
    parallelizable manifold.

    An instance of this class is a vector field along a differentiable
    manifold `U` with values on a parallelizable manifold `M`, via a
    differentiable map `U \rightarrow M`. More precisely, given a
    differentiable map

    .. MATH::

        \Phi:\ U \longrightarrow M,

    a *vector field along* `U` *with values on* `M` is a differentiable map

    .. MATH::

        v:\ U  \longrightarrow TM

    (`TM` being the tangent bundle of `M`) such that

    .. MATH::

        \forall p \in U,\ v(p) \in T_{\Phi(p)}M

    The standard case of vector fields *on* a differentiable manifold
    corresponds to `U=M` and `\Phi = \mathrm{Id}_M`. Other common cases are
    `\Phi` being an immersion and `\Phi` being a curve in `M` (`U` is then an
    open interval of `\RR`).

    If `M` is not parallelizable, the class
    :class:`~sage.manifolds.differentiable.vectorfield.VectorField`
    must be used instead.

    This is a Sage *element* class, the corresponding *parent* class being
    :class:`~sage.manifolds.differentiable.vectorfield_module.VectorFieldFreeModule`.

    INPUT:

    - ``vector_field_module`` -- free module `\mathcal{X}(U,\Phi)` of vector
      fields along `U` with values on `M\supset\Phi(U)`
    - ``name`` -- (default: ``None``) name given to the vector field
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the vector
      field; if none is provided, the LaTeX symbol is set to ``name``

    EXAMPLES:

    A vector field on a parallelizable 3-dimensional manifold::

        sage: M = Manifold(3, 'M')
        sage: c_xyz.<x,y,z> = M.chart()
        sage: v = M.vector_field('V') ; v
        Vector field V on the 3-dimensional differentiable manifold M
        sage: latex(v)
        V

    Vector fields are considered as elements of a module over the ring
    (algebra) of scalar fields on `M`::

        sage: v.parent()
        Free module X(M) of vector fields on the 3-dimensional differentiable
         manifold M
        sage: v.parent().base_ring()
        Algebra of differentiable scalar fields on the 3-dimensional
         differentiable manifold M
        sage: v.parent() is M.vector_field_module()
        True

    A vector field is a tensor field of rank 1 and of type (1,0)::

        sage: v.tensor_rank()
        1
        sage: v.tensor_type()
        (1, 0)

    Components of a vector field with respect to a given frame::

        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v[0], v[1], v[2] = (1, 4, 9)  # components on M's default frame (e)
        sage: v.comp()
        1-index components w.r.t. Vector frame (M, (e_0,e_1,e_2))

    The totality of the components are accessed via the operator [:]::

        sage: v[:] = (1, 4, 9)  # equivalent to v[0], v[1], v[2] = (1, 4, 9)
        sage: v[:]
        [1, 4, 9]

    The components are also read on the expansion on the frame 'e', as provided
    by the method
    :meth:`~sage.tensor.modules.free_module_tensor.FreeModuleTensor.display`::

        sage: v.display()  # displays the expansion on the manifold's default frame (e)
        V = e_0 + 4 e_1 + 9 e_2

    A subset of the components can be accessed by means of Python's slice
    notation::

        sage: v[1:] = (-2, -3)
        sage: v[:]
        [1, -2, -3]
        sage: v[:2]
        [1, -2]

    The components are instances of class
    :class:`~sage.tensor.modules.comp.Components`::

        sage: type(v.comp())
        <class 'sage.tensor.modules.comp.Components'>

    Components in another frame::

        sage: f = M.vector_frame('f')
        sage: for i in range(3):
        ....:     v.set_comp(f)[i] = (i+1)**3
        ....:
        sage: v.comp(f)[2]
        27
        sage: v[f, 2]  # equivalent to above
        27
        sage: v.display(f)
        V = f_0 + 8 f_1 + 27 f_2

    The range of the indices depends on the convention set for the manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: c_xyz.<x,y,z> = M.chart()
        sage: e = M.vector_frame('e') ; M.set_default_frame(e)
        sage: v = M.vector_field('V')
        sage: (v[1], v[2], v[3]) = (1, 4, 9)
        sage: v[0]
        Traceback (most recent call last):
        ...
        IndexError: index out of range: 0 not in [1, 3]

    A vector field acts on scalar fields (derivation along the vector field)::

        sage: M = Manifold(2, 'M')
        sage: c_cart.<x,y> = M.chart()
        sage: f = M.scalar_field(x*y^2, name='f')
        sage: v = M.vector_field('v')
        sage: v[:] = (-y, x)
        sage: v.display()
        v = -y d/dx + x d/dy
        sage: v(f)
        Scalar field v(f) on the 2-dimensional differentiable manifold M
        sage: v(f).expr()
        2*x^2*y - y^3
        sage: latex(v(f))
        v\left(f\right)

    Example of a vector field associated with a non-trivial map `\Phi`:
    vector field along a curve in `M`::

        sage: R = Manifold(1, 'R')  # R as a 1-dimensional manifold
        sage: T.<t> = R.chart()  # canonical chart on R
        sage: Phi = R.diff_map(M, [cos(t), sin(t)], name='Phi') ; Phi
        Differentiable map Phi from the 1-dimensional differentiable manifold R
         to the 2-dimensional differentiable manifold M
        sage: Phi.display()
        Phi: R --> M
           t |--> (x, y) = (cos(t), sin(t))
        sage: w = R.vector_field('w', dest_map=Phi) ; w
        Vector field w along the 1-dimensional differentiable manifold R with
         values on the 2-dimensional differentiable manifold M
        sage: w.parent()
        Free module X(R,Phi) of vector fields along the 1-dimensional
         differentiable manifold R mapped into the 2-dimensional differentiable
         manifold M
        sage: w[:] = (-sin(t), cos(t))
        sage: w.display()
        w = -sin(t) d/dx + cos(t) d/dy

    Value at a given point::

        sage: p = R((0,), name='p') ; p
        Point p on the 1-dimensional differentiable manifold R
        sage: w.at(p)
        Tangent vector w at Point Phi(p) on the 2-dimensional differentiable
         manifold M
        sage: w.at(p).display()
        w = d/dy
        sage: w.at(p) == v.at(Phi(p))
        True

    """
    def __init__(self, vector_field_module, name=None, latex_name=None):
        r"""
        Construct a vector field with values on a parallelizable manifold.

        TESTS:

        Construction via ``parent.element_class``, and not via a direct call
        to ``VectorFieldParal``, to fit with the category framework::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: XM = M.vector_field_module()  # the parent
            sage: v = XM.element_class(XM, name='v'); v
            Vector field v on the 2-dimensional differentiable manifold M
            sage: v[:] = (-y, x)
            sage: v.display()
            v = -y d/dx + x d/dy
            sage: TestSuite(v).run()

        Construction via ``DifferentiableManifold.vector_field``::

            sage: u = M.vector_field(name='u'); u
            Vector field u on the 2-dimensional differentiable manifold M
            sage: type(u) == type(v)
            True
            sage: u.parent() is v.parent()
            True
            sage: u[:] = (1+x, 1-y)
            sage: TestSuite(u).run()

        """
        FiniteRankFreeModuleElement.__init__(self, vector_field_module,
                                             name=name, latex_name=latex_name)
        # TensorFieldParal attributes:
        self._domain = vector_field_module._domain
        self._ambient_domain = vector_field_module._ambient_domain
        # VectorField attributes:
        self._vmodule = vector_field_module
        # Initialization of derived quantities:
        TensorFieldParal._init_derived(self)
        VectorField._init_derived(self)
        # Initialization of list of quantities depending on self:
        self._init_dependencies()

    def _repr_(self) :
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._repr_()
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: repr(v)  # indirect doctest
            'Vector field v on the 2-dimensional differentiable manifold M'
            sage: v  # indirect doctest
            Vector field v on the 2-dimensional differentiable manifold M

        """
        return VectorField._repr_(self)

    def _new_instance(self):
        r"""
        Create an instance of the same class as ``self`` on the same module.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: u = v._new_instance(); u
            Vector field on the 2-dimensional differentiable manifold M
            sage: u.parent() is v.parent()
            True

        """
        return type(self)(self._fmodule)

    def _del_derived(self, del_restrictions=True):
        r"""
        Delete the derived quantities

        INPUT:

        - ``del_restrictions`` -- (default: True) determines whether the
          restrictions of ``self`` to subdomains are deleted.

        TEST::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()  # makes M parallelizable
            sage: v = M.vector_field(name='v')
            sage: v._del_derived()

        """
        TensorFieldParal._del_derived(self, del_restrictions=del_restrictions)
        VectorField._del_derived(self)
        self._del_dependencies()

    def __call__(self, scalar):
        r"""
        Action on a scalar field.

        INPUT:

        - ``scalar`` -- scalar field `f`

        OUTPUT:

        - scalar field representing the derivative of `f` along the vector
          field, i.e. `v^i \frac{\partial f}{\partial x^i}`

        EXAMPLES:

        Action of a vector field on a scalar field on a 2-dimensional
        manifold::

            sage: M = Manifold(2, 'M')
            sage: c_cart.<x,y> = M.chart()
            sage: f = M.scalar_field(x*y^2)
            sage: v = M.vector_field()
            sage: v[:] = (-y, x)
            sage: v(f)
            Scalar field on the 2-dimensional differentiable manifold M
            sage: v(f).display()
            M --> R
            (x, y) |--> 2*x^2*y - y^3

        """
        from sage.manifolds.differentiable.vectorframe import CoordFrame
        if scalar._tensor_type == (0,1):
            # This is actually the action of the vector field on a 1-form,
            # as a tensor field of type (1,0):
            return scalar(self)
        if scalar._tensor_type != (0,0):
            raise TypeError("the argument must be a scalar field")
        #!# Could it be simply
        # return scalar.differential()(self)
        # ?
        dom_resu = self._domain.intersection(scalar._domain)
        self_r = self.restrict(dom_resu)
        scalar_r = scalar.restrict(dom_resu)
        if scalar_r._is_zero:
            return dom_resu._zero_scalar_field
        # Creation of the result:
        if self._name is not None and scalar._name is not None:
            resu_name = self._name + "(" + scalar._name + ")"
        else:
            resu_name = None
        if self._latex_name is not None and scalar._latex_name is not None:
            resu_latex = self._latex_name + r"\left(" + scalar._latex_name + \
                        r"\right)"
        else:
            resu_latex = None
        resu = dom_resu.scalar_field(name=resu_name, latex_name=resu_latex)
        # Search for common charts for the computation:
        common_charts = set()
        for chart in scalar_r._express:
            try:
                self_r.comp(chart._frame)
                common_charts.add(chart)
            except ValueError:
                pass
        for frame in self_r._components:
            if isinstance(frame, CoordFrame):
                chart = frame._chart
                try:
                    scalar_r.coord_function(chart)
                    common_charts.add(chart)
                except ValueError:
                    pass
        if not common_charts:
            raise ValueError("no common chart found")
        # The computation:
        manif = scalar._manifold
        for chart in common_charts:
            v = self_r.comp(chart._frame)
            f = scalar_r.coord_function(chart)
            res = 0
            for i in manif.irange():
                res += v[i, chart] * f.diff(i)
            resu._express[chart] = res
        return resu
