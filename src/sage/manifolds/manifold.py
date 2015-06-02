r"""
Topological manifolds

Given a topological field `K` (in most applications, `K = \RR` or
`K = \CC`) and a non-negative integer `n`, a *topological manifold of
dimension* `n` *over K* is a topological space `M` such that

- `M` is a Hausdorff space,
- `M` is second countable,
- every point in `M` has a neighborhood homeomorphic to `K^n`

Topological manifolds are implemented via the class :class:`TopManifold`.
Open subsets of topological manifolds are also implemented via
:class:`TopManifold`, since they are topological manifolds by themselves.

In the current setting, topological manifolds are mostly described by means of
charts (see :class:`~sage.manifolds.chart.Chart`).

:class:`TopManifold` serves as a base class for more specific manifold classes.

.. RUBRIC:: Example 1: the 2-sphere as a topological manifold of dimension
  2 over `\RR`

One starts by declaring `S^2` as a 2-dimensional topological manifold::

    sage: M = TopManifold(2, 'S^2')
    sage: M
    2-dimensional topological manifold S^2

Since the base topological field has not been specified in the argument list
of ``TopManifold``, `\RR` is assumed::

    sage: M.base_field()
    'real'
    sage: dim(M)
    2

Let us consider the complement of a point, the "North pole" say; this is an
open subset of `S^2`, which we call `U`::

    sage: U = M.open_subset('U'); U
    Open subset U of the 2-dimensional topological manifold S^2

A standard chart on `U` is provided by the stereographic projection from the
North pole to the equatorial plane::

    sage: stereoN.<x,y> = U.chart(); stereoN
    Chart (U, (x, y))

Thanks to the operator ``<x,y>`` on the left-hand side, the coordinates
declared in a chart (here `x` and `y`), are accessible by their names; they are
Sage's symbolic variables::

    sage: y
    y
    sage: type(y)
    <type 'sage.symbolic.expression.Expression'>

The South pole is the point of coordinates `(x,y)=(0,0)` in the above
chart::

    sage: S = U.point((0,0), chart=stereoN, name='S'); S
    Point S on the 2-dimensional topological manifold S^2

Let us call `V` the open subset that is the complement of the South pole and
let us introduce on it the chart induced by the stereographic projection from
the South pole to the equatorial plane::

    sage: V = M.open_subset('V'); V
    Open subset V of the 2-dimensional topological manifold S^2
    sage: stereoS.<u,v> = V.chart(); stereoS
    Chart (V, (u, v))

The North pole is the point of coordinates `(u,v)=(0,0)` in this chart::

    sage: N = V.point((0,0), chart=stereoS, name='N'); N
    Point N on the 2-dimensional topological manifold S^2

To fully construct the manifold, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the charts ``stereoN`` = `(U, (x, y))`
and ``stereoS`` = `(V, (u, v))`, denoting by W the intersection of U and V
(W is the subset of U defined by `x^2+y^2\not=0`, as well as the subset of V
defined by`u^2+v^2\not=0`)::

    sage: stereoN_to_S = stereoN.transition_map(stereoS, [x/(x^2+y^2), y/(x^2+y^2)],
    ....:                intersection_name='W', restrictions1= x^2+y^2!=0, restrictions2= u^2+v^2!=0)
    sage: stereoN_to_S
    Change of coordinates from Chart (W, (x, y)) to Chart (W, (u, v))
    sage: stereoN_to_S.display()
    u = x/(x^2 + y^2)
    v = y/(x^2 + y^2)

We give the name ``W`` to the Python variable representing `W=U\cap V`::

    sage: W = U.intersection(V)

The inverse of the transition map is computed by the method ``inverse()``::

    sage: stereoN_to_S.inverse()
    Change of coordinates from Chart (W, (u, v)) to Chart (W, (x, y))
    sage: stereoN_to_S.inverse().display()
    x = u/(u^2 + v^2)
    y = v/(u^2 + v^2)

At this stage, we have four open subsets on `S^2`::

    sage: M.list_of_subsets()
    [2-dimensional topological manifold S^2,
     Open subset U of the 2-dimensional topological manifold S^2,
     Open subset V of the 2-dimensional topological manifold S^2,
     Open subset W of the 2-dimensional topological manifold S^2]

`W` is the open subset that is the complement of the two poles::

    sage: N in W or S in W
    False


The North pole lies in `V` and the South pole in `U`::

    sage: N in V, N in U
    (True, False)
    sage: S in U, S in V
    (True, False)

The manifold's (user) atlas contains four charts, two of them
being restrictions of charts to a smaller domain::

    sage: M.atlas()
    [Chart (U, (x, y)), Chart (V, (u, v)), Chart (W, (x, y)), Chart (W, (u, v))]

Let us consider the point of coordinates (1,2) in the chart ``stereoN``::

    sage: p = M.point((1,2), chart=stereoN, name='p'); p
    Point p on the 2-dimensional topological manifold S^2
    sage: p.parent()
    2-dimensional topological manifold S^2
    sage: p in W
    True

The coordinates of `p` in the chart ``stereoS`` are::

    sage: stereoS(p)
    (1/5, 2/5)

Given the definition of `p`, we have of course::

    sage: stereoN(p)
    (1, 2)

Similarly::

    sage: stereoS(N)
    (0, 0)
    sage: stereoN(S)
    (0, 0)


.. RUBRIC:: Example 2: the Riemann sphere as a topological manifold of
  dimension 1 over `\CC`

We declare the Riemann sphere `\CC^*` as a 1-dimensional topological manifold
over `\CC`::

    sage: M = TopManifold(1, 'C*', field='complex'); M
    Complex 1-dimensional topological manifold C*

We introduce a first open subset, which is actually
`\CC = \CC^*\setminus\{\infty\}` if we interpret `\CC^*` as the Alexandroff
one-point compactification of `\CC`::

    sage: U = M.open_subset('U')

A natural chart on `U` is then nothing but the identity map of `\CC`, hence
we denote the associated coordinate by `z`::

    sage: Z.<z> = U.chart()

The origin of the complex plane is the point of coordinate `z=0`::

    sage: O = U.point((0,), chart=Z, name='O'); O
    Point O on the Complex 1-dimensional topological manifold C*

Another open subset of `\CC^*` is `V = \CC^*\setminus\{O\}`::

    sage: V = M.open_subset('V')

We define a chart on `V` such that the point at infinity is the point of
coordinate 0 in this chart::

    sage: W.<w> = V.chart(); W
    Chart (V, (w,))
    sage: inf = M.point((0,), chart=W, name='inf', latex_name=r'\infty')
    sage: inf
    Point inf on the Complex 1-dimensional topological manifold C*

To fully construct the Riemann sphere, we declare that it is the union of `U`
and `V`::

    sage: M.declare_union(U,V)

and we provide the transition map between the two charts as `w=1/z` on
on `A = U\cap V`::

    sage: Z_to_W = Z.transition_map(W, 1/z, intersection_name='A',
    ....:                           restrictions1= z!=0, restrictions2= w!=0)
    sage: Z_to_W
    Change of coordinates from Chart (A, (z,)) to Chart (A, (w,))
    sage: Z_to_W.display()
    w = 1/z
    sage: Z_to_W.inverse()
    Change of coordinates from Chart (A, (w,)) to Chart (A, (z,))
    sage: Z_to_W.inverse().display()
    z = 1/w

Let consider the complex number `i` as a point of the Riemann sphere::

    sage: i = M((I,), chart=Z, name='i'); i
    Point i on the Complex 1-dimensional topological manifold C*

Its coordinates w.r.t. the charts ``Z`` and ``W`` are::

    sage: Z(i)
    (I,)
    sage: W(i)
    (-I,)

and we have::

    sage: i in U
    True
    sage: i in V
    True

The following subsets and charts have been defined::

    sage: M.list_of_subsets()
    [Open subset A of the Complex 1-dimensional topological manifold C*,
     Complex 1-dimensional topological manifold C*,
     Open subset U of the Complex 1-dimensional topological manifold C*,
     Open subset V of the Complex 1-dimensional topological manifold C*]
    sage: M.atlas()
    [Chart (U, (z,)), Chart (V, (w,)), Chart (A, (z,)), Chart (A, (w,))]


AUTHORS:

- Eric Gourgoulhon (2015): initial version

REFERENCES:

- J.M. Lee : *Introduction to Topological Manifolds*, 2nd ed., Springer (New
  York) (2011)

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

from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.fields import Fields
#*# Before #18175:
from sage.categories.sets_cat import Sets
#*# After #18175, this should become
# from sage.categories.manifolds import Manifolds
from sage.manifolds.subset import TopManifoldSubset

class TopManifold(TopManifoldSubset):
    r"""
    Topological manifold over a topological field `K`.

    Given a topological field `K` (in most applications, `K = \RR` or
    `K = \CC`) and a non-negative integer `n`, a *topological manifold of
    dimension* `n` *over K* is a topological space `M` such that

    - `M` is a Hausdorff space,
    - `M` is second countable,
    - every point in `M` has a neighborhood homeomorphic to `K^n`

    This is a Sage *parent* class, the corresponding *element*
    class being :class:`~sage.manifolds.point.TopManifoldPoint`.

    INPUT:

    - ``n`` -- positive integer; dimension of the manifold
    - ``name`` -- string; name (symbol) given to the manifold
    - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
      manifold; if none is provided, it is set to ``name``
    - ``field`` -- (default: 'real') field `K` on which the manifold is
      defined; allowed values are

        - 'real' for a manifold over `\RR`
        - 'complex' for a manifold over `\CC`
        - any object in the category of fields (see
          :class:`~sage.categories.fields.Fields`) for more general manifolds

    - ``start_index`` -- (default: 0) integer; lower value of the range of
      indices used for "indexed objects" on the manifold, e.g. coordinates
      in a chart.
    - ``category`` -- (default: ``None``) to specify the categeory; the
      default being ``Sets()`` (``Manifolds()`` after :trac:`18175` is
      implemented)
    - ``ambient_manifold`` -- (default: ``None``) if not ``None``, the created
      object is considered as an open subset of the topological manifold
      ``ambient_manifold``

    EXAMPLES:

    A 4-dimensional topological manifold (over `\RR`)::

        sage: M = TopManifold(4, 'M', latex_name=r'\mathcal{M}')
        sage: M
        4-dimensional topological manifold M
        sage: latex(M)
        \mathcal{M}
        sage: M.base_field()
        'real'
        sage: dim(M)
        4

    The input parameter ``start_index`` defines the range of indices on the
    manifold::

        sage: M = TopManifold(4, 'M')
        sage: list(M.irange())
        [0, 1, 2, 3]
        sage: M = TopManifold(4, 'M', start_index=1)
        sage: list(M.irange())
        [1, 2, 3, 4]
        sage: list(TopManifold(4, 'M', start_index=-2).irange())
        [-2, -1, 0, 1]

    A complex manifold::

        sage: N = TopManifold(3, 'N', field='complex'); N
        Complex 3-dimensional topological manifold N

    A manifold over `\QQ`::

        sage: N = TopManifold(6, 'N', field=QQ); N
        6-dimensional topological manifold N over the Rational Field

    A manifold is a Sage *parent* object, in the category of sets::

        sage: isinstance(M, Parent)
        True
        sage: M.category()
        Category of sets
        sage: M in Sets()
        True

    The corresponding Sage *elements* are points::

        sage: X.<t, x, y, z> = M.chart()
        sage: p = M.an_element(); p
        Point on the 4-dimensional topological manifold M
        sage: p.parent()
        4-dimensional topological manifold M
        sage: M.is_parent_of(p)
        True
        sage: p in M
        True

    The manifold's points are instances of class
    :class:`~sage.manifolds.point.TopManifoldPoint`::

        sage: isinstance(p, sage.manifolds.point.TopManifoldPoint)
        True

    Manifolds are unique, as long as they are created with the same arguments::

        sage: M is TopManifold(4, 'M', start_index=1)
        True
        sage: M is TopManifold(4, 'M')
        False
        sage: M is TopManifold(4, 'M', latex_name='M', start_index=1)
        False

    Since an open subset of a topological manifold `M` is itself a topological
    manifold, open subsets of `M` are instances of the class
    :class:`TopManifold`::

        sage: U = M.open_subset('U'); U
        Open subset U of the 4-dimensional topological manifold M
        sage: isinstance(U, TopManifold)
        True
        sage: U.base_field() == M.base_field()
        True
        sage: dim(U) == dim(M)
        True

    The manifold passes all the tests of the test suite relative to the
    category of Sets::

        sage: TestSuite(M).run(verbose=True)
        running ._test_an_element() . . . pass
        running ._test_category() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass


    """
    def __init__(self, n, name, latex_name=None, field='real', start_index=0,
                 category=None, ambient_manifold=None):
        r"""
        Construct a topological manifold.

        TESTS::

            sage: M = TopManifold(3, 'M', latex_name=r'\mathbb{M}', start_index=1)
            sage: M
            3-dimensional topological manifold M
            sage: latex(M)
            \mathbb{M}
            sage: dim(M)
            3
            sage: X.<x,y,z> = M.chart()
            sage: TestSuite(M).run()

        """
        # Initialization of the attributes _dim, _field and _start_index:
        from sage.rings.integer import Integer
        if not isinstance(n, (int, Integer)):
            raise TypeError("the manifold dimension must be an integer")
        if n<1:
            raise ValueError("the manifold dimension must be strictly " +
                             "positive")
        self._dim = n
        if field not in ['real', 'complex']:
            if field not in Fields():
                raise TypeError("the argument 'field' must be a field")
        self._field = field
        if not isinstance(start_index, (int, Integer)):
            raise TypeError("the starting index must be an integer")
        self._sindex = start_index
        if category is None:
            category = Sets()
            #*# After #18175, this should become
            # category = Manifolds()
        if ambient_manifold is None:
            ambient_manifold = self
        elif not isinstance(ambient_manifold, TopManifold):
            raise TypeError("the argument 'ambient_manifold' must be " +
                            " a topological manifold")
        # Initialization as a subset of the ambient manifold (possibly itself):
        TopManifoldSubset.__init__(self, ambient_manifold, name,
                                   latex_name=latex_name, category=category)
        self._is_open = True
        self._atlas = []  # list of charts defined on subsets of self
        self._top_charts = []  # list of charts defined on subsets of self
                        # that are not subcharts of charts on larger subsets
        self._def_chart = None  # default chart
        self._coord_changes = {} # dictionary of transition maps
        # list of charts that individually cover self, i.e. whose
        # domains are self (if non-empty, self is a coordinate domain):
        self._covering_charts = []
        # algebra of scalar fields defined on self:
        #*# self._scalar_field_algebra = ScalarFieldAlgebra(self)
        # the zero scalar field:
        #*# self._zero_scalar_field = self._scalar_field_algebra.zero()
        # the identity map on self:
        #*# self._identity_map = Hom(self, self).one()

    def _repr_(self):
        r"""
        String representation of the manifold.

        TESTS::

            sage: M = TopManifold(3, 'M')
            sage: M._repr_()
            '3-dimensional topological manifold M'
            sage: repr(M)  # indirect doctest
            '3-dimensional topological manifold M'
            sage: M  # indirect doctest
            3-dimensional topological manifold M
            sage: M = TopManifold(3, 'M', field='complex')
            sage: M._repr_()
            'Complex 3-dimensional topological manifold M'
            sage: M = TopManifold(3, 'M', field=QQ)
            sage: M._repr_()
            '3-dimensional topological manifold M over the Rational Field'

        If the manifold is actually an open subset of a larger manifold, the
        string representation is different::

            sage: U = M.open_subset('U')
            sage: U._repr_()
            'Open subset U of the 3-dimensional topological manifold M over the Rational Field'

        """
        if self._manifold is self:
            if self._field == 'real':
                return "{}-dimensional topological manifold {}".format(
                                                         self._dim, self._name)
            elif self._field == 'complex':
                return "Complex {}-dimensional topological manifold {}".format(
                                                         self._dim, self._name)
            return "{}-dimensional topological manifold {} over the {}".format(
                                            self._dim, self._name, self._field)
        else:
            return "Open subset {} of the {}".format(self._name,
                                                     self._manifold)

    def _latex_(self):
        r"""
        LaTeX representation of the manifold.

        TESTS::

            sage: M = TopManifold(3, 'M')
            sage: M._latex_()
            'M'
            sage: latex(M)
            M
            sage: M = TopManifold(3, 'M', latex_name=r'\mathcal{M}')
            sage: M._latex_()
            '\\mathcal{M}'
            sage: latex(M)
            \mathcal{M}

        """
        return self._latex_name

    def _an_element_(self):
        r"""
        Construct some point on the manifold.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p.coord()
            (0, 0)
            sage: U = M.open_subset('U', coord_def={X: y>1}); U
            Open subset U of the 2-dimensional topological manifold M
            sage: p = U._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in U
            True
            sage: p.coord()
            (0, 2)
            sage: V = U.open_subset('V', coord_def={X.restrict(U): x<-pi})
            sage: p = V._an_element_(); p
            Point on the 2-dimensional topological manifold M
            sage: p in V
            True
            sage: p.coord()
            (-pi - 1, 2)

        """
        from sage.rings.infinity import Infinity
        if self._def_chart is None:
            return self.element_class(self)
        # Attempt to construct a point in the domain of the default chart
        chart = self._def_chart
        if self._field == 'real':
            coords = []
            for coord_range in chart._bounds:
                xmin = coord_range[0][0]
                xmax = coord_range[1][0]
                if xmin == -Infinity:
                    if xmax == Infinity:
                        x = 0
                    else:
                        x = xmax - 1
                else:
                    if xmax == Infinity:
                        x = xmin + 1
                    else:
                        x = (xmin + xmax)/2
                coords.append(x)
        else:
            coords = self._dim*[0]
        if not chart.valid_coordinates(*coords):
            # Attempt to construct a point in the domain of other charts
            if self._field == 'real':
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    coords = []
                    for coord_range in ch._bounds:
                        xmin = coord_range[0][0]
                        xmax = coord_range[1][0]
                        if xmin == -Infinity:
                            if xmax == Infinity:
                                x = 0
                            else:
                                x = xmax - 1
                        else:
                            if xmax == Infinity:
                                x = xmin + 1
                            else:
                                x = (xmin + xmax)/2
                        coords.append(x)
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    # A generic element with specific coordinates could not be
                    # automatically generated, due to too complex cooordinate
                    # conditions. An element without any coordinate set is
                    # returned instead:
                    return self.element_class(self)
            else:
                # Case of manifolds over a field different from R
                for ch in self._atlas:
                    if ch is self._def_chart:
                        continue # since this case has already been attempted
                    if ch.valid_coordinates(*coords):
                        chart = ch
                        break
                else:
                    return self.element_class(self)
        # The point is constructed with check_coords=False since the check
        # has just been performed above:
        return self.element_class(self, coords=coords, chart=chart,
                                  check_coords=False)

    def __contains__(self, point):
        r"""
        Check whether a point is contained in the manifold.

        EXAMPLES::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,2), chart=X)
            sage: M.__contains__(p)
            True
            sage: p in M  # indirect doctest
            True
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: U.__contains__(p)
            True
            sage: p in U  # indirect doctest
            True
            sage: V = U.open_subset('V', coord_def={X.restrict(U): y<0})
            sage: V.__contains__(p)
            False
            sage: p in V  # indirect doctest
            False

        """
        # for efficiency, a quick test first:
        if point._subset is self:
            return True
        if point._subset.is_subset(self):
            return True
        for chart in self._atlas:
            if chart in point._coordinates:
                if chart.valid_coordinates( *(point._coordinates[chart]) ):
                    return True
        for chart in point._coordinates:
            for schart in chart._subcharts:
                if schart in self._atlas and schart.valid_coordinates(
                                          *(point._coordinates[chart]) ):
                    return True
        return False

    def dimension(self):
        r"""
        Return the dimension of the manifold over its base field.

        EXAMPLE::

            sage: M = TopManifold(2, 'M')
            sage: M.dimension()
            2

        A shortcut is ``dim()``::

            sage: M.dim()
            2

        The Sage global function ``dim`` can also be used::

            sage: dim(M)
            2

        """
        return self._dim

    dim = dimension

    def base_field(self):
        r"""
        Return the field on which the manifolds is defined.

        OUTPUT:

        - a field or one of the strings 'real' and 'complex', since there
          is no exact representation of the fields `\RR` and `\CC` in Sage

        EXAMPLES::

            sage: M = TopManifold(3, 'M')
            sage: M.base_field()
            'real'
            sage: M = TopManifold(3, 'M', field='complex')
            sage: M.base_field()
            'complex'
            sage: M = TopManifold(3, 'M', field=QQ)
            sage: M.base_field()
            Rational Field

        """
        return self._field

    def start_index(self):
        r"""
        Return the first value of the index range used on the manifold.

        This is the parameter ``start_index`` passed at the construction of
        the manifold.

        OUTPUT:

        - the integer `i_0` such that all indices of indexed objects on the
          manifold range from `i_0` to `i_0 + n - 1`, where `n` is the
          manifold's dimension.

        EXAMPLES::

            sage: M = TopManifold(3, 'M')
            sage: M.start_index()
            0
            sage: M = TopManifold(3, 'M', start_index=1)
            sage: M.start_index()
            1

        """
        return self._sindex

    def irange(self, start=None):
        r"""
        Single index generator.

        INPUT:

        - ``start`` -- (default: ``None``) initial value `i_0` of the index; if
          none is provided, the value returned by :meth:`start_index()` is
          assumed.

        OUTPUT:

        - an iterable index, starting from `i_0` and ending at
          `i_0 + n - 1`, where `n` is the manifold's dimension.

        EXAMPLES:

        Index range on a 4-dimensional manifold::

            sage: M = TopManifold(4, 'M')
            sage: for i in M.irange():
            ...       print i,
            ...
            0 1 2 3
            sage: for i in M.irange(2):
            ...       print i,
            ...
            2 3
            sage: list(M.irange())
            [0, 1, 2, 3]

        Index range on a 4-dimensional manifold with starting index=1::

            sage: M = TopManifold(4, 'M', start_index=1)
            sage: for i in M.irange():
            ...       print i,
            ...
            1 2 3 4
            sage: for i in M.irange(2):
            ...      print i,
            ...
            2 3 4

        In general, one has always::

            sage: M.irange().next() == M.start_index()
            True

        """
        si = self._sindex
        imax = self._dim + si
        if start is None:
            i = si
        else:
            i = start
        while i < imax:
            yield i
            i += 1


    def index_generator(self, nb_indices):
        r"""
        Generator of index series.

        INPUT:

        - ``nb_indices`` -- number of indices in a series

        OUTPUT:

        - an iterable index series for a generic component with the specified
          number of indices

        EXAMPLES:

        Indices on a 2-dimensional manifold::

            sage: M = TopManifold(2, 'M', start_index=1)
            sage: for ind in M.index_generator(2):
            ...       print ind
            ...
            (1, 1)
            (1, 2)
            (2, 1)
            (2, 2)

        Loops can be nested::

            sage: for ind1 in M.index_generator(2):
            ...       print ind1, " : ",
            ...       for ind2 in M.index_generator(2):
            ...           print ind2,
            ...       print ""
            ...
            (1, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (1, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 1)  :  (1, 1) (1, 2) (2, 1) (2, 2)
            (2, 2)  :  (1, 1) (1, 2) (2, 1) (2, 2)

        """
        si = self._sindex
        imax = self._dim - 1 + si
        ind = [si for k in range(nb_indices)]
        ind_end = [si for k in range(nb_indices)]
        ind_end[0] = imax+1
        while ind != ind_end:
            yield tuple(ind)
            ret = 1
            for pos in range(nb_indices-1,-1,-1):
                if ind[pos] != imax:
                    ind[pos] += ret
                    ret = 0
                elif ret == 1:
                    if pos == 0:
                        ind[pos] = imax + 1 # end point reached
                    else:
                        ind[pos] = si
                        ret = 1

    def atlas(self):
        r"""
        Return the list of charts that have been defined on the manifold.

        EXAMPLES:

        Charts on subsets of `\RR^2`::

            sage: M = TopManifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: M.atlas()
            [Chart (R^2, (x, y))]
            sage: U = M.open_subset('U', coord_def={c_cart: (y!=0,x<0)}) # U = R^2 \ half line {y=0,x>=0}
            sage: U.atlas()
            [Chart (U, (x, y))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y))]
            sage: c_spher.<r, ph> = U.chart(r'r:(0,+oo) ph:(0,2*pi):\phi') # spherical (polar) coordinates on U
            sage: U.atlas()
            [Chart (U, (x, y)), Chart (U, (r, ph))]
            sage: M.atlas()
            [Chart (R^2, (x, y)), Chart (U, (x, y)), Chart (U, (r, ph))]

        """
        return self._atlas

    def top_charts(self):
        r"""
        Return the list of charts defined on subsets of the current set
        that are not subcharts of charts on larger subsets.

        OUTPUT:

        - list of charts defined on open subsets of ``self`` but not on
          larger subsets

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_()  # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: Y.<u,v> = U.chart()
            sage: M.top_charts()
            [Chart (M, (x, y)), Chart (U, (u, v))]

        Note that the (user) atlas contains one more chart: (U, (x,y)), which
        is not a "top" chart::

            sage: M.atlas()
            [Chart (M, (x, y)), Chart (U, (x, y)), Chart (U, (u, v))]

        """
        return self._top_charts

    def default_chart(self):
        r"""
        Return the default chart defined on the manifold.

        Unless changed via :meth:`set_default_chart`, the *default chart*
        is the first one defined on a subset of the manifold (possibly itself).

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.Chart`
          representing the default chart.

        EXAMPLES:

        Default chart on a 2-dimensional manifold and on some subsets::

            sage: TopManifold._clear_cache_()  # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: M.chart('x y')
            Chart (M, (x, y))
            sage: M.chart('u v')
            Chart (M, (u, v))
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: A = M.open_subset('A')
            sage: A.chart('t z')
            Chart (A, (t, z))
            sage: A.default_chart()
            Chart (A, (t, z))

        """
        return self._def_chart

    def set_default_chart(self, chart):
        r"""
        Changing the default chart on ``self``.

        INPUT:

        - ``chart`` -- a chart (must be defined on some subset ``self``)

        EXAMPLES:

        Charts on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: M.default_chart()
            Chart (M, (x, y))
            sage: M.set_default_chart(c_uv)
            sage: M.default_chart()
            Chart (M, (u, v))

        """
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError(str(chart) + " is not a chart.")
        if chart._domain is not self:
            if self.is_manifestly_coordinate_domain():
                raise TypeError("The chart domain must coincide with the " +
                                str(self) + ".")
            if chart not in self._atlas:
                raise ValueError("The chart must be defined on the " +
                                 str(self))
        self._def_chart = chart

    def coord_change(self, chart1, chart2):
        r"""
        Return the change of coordinates (transition map) between two charts
        defined on the manifold.

        The change of coordinates must have been defined previously, for
        instance by the method
        :meth:`~sage.manifolds.chart.Chart.transition_map`.

        INPUT:

        - ``chart1`` -- chart 1
        - ``chart2`` -- chart 2

        OUTPUT:

        - instance of :class:`~sage.manifolds.chart.CoordChange`
          representing the transition map from chart 1 to chart 2

        EXAMPLES:

        Change of coordinates on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: c_xy.transition_map(c_uv, (x+y, x-y)) # defines the coordinate change
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))
            sage: M.coord_change(c_xy, c_uv) # returns the coordinate change defined above
            Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))

        """
        if (chart1, chart2) not in self._coord_changes:
            raise TypeError("The change of coordinates from " + str(chart1) +
                            " to " + str(chart2) + " has not been " +
                            "defined on the " + str(self))
        return self._coord_changes[(chart1, chart2)]


    def coord_changes(self):
        r"""
        Return the changes of coordinates (transition maps) defined on
        subsets of the manifold.

        OUTPUT:

        - dictionary of changes of coordinates, with pairs of charts as keys

        EXAMPLES:

        Various changes of coordinates on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: c_xy.<x,y> = M.chart()
            sage: c_uv.<u,v> = M.chart()
            sage: xy_to_uv = c_xy.transition_map(c_uv, [x+y, x-y])
            sage: M.coord_changes()
            {(Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: uv_to_xy = xy_to_uv.inverse()
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: c_rs.<r,s> = M.chart()
            sage: uv_to_rs = c_uv.transition_map(c_rs, [-u+2*v, 3*u-v])
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v))}
            sage: xy_to_rs = uv_to_rs * xy_to_uv
            sage: M.coord_changes()  # random (dictionary output)
            {(Chart (M, (u, v)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (r, s)),
             (Chart (M, (u, v)),
              Chart (M, (x, y))): Change of coordinates from Chart (M, (u, v)) to Chart (M, (x, y)),
             (Chart (M, (x, y)),
              Chart (M, (u, v))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (u, v)),
             (Chart (M, (x, y)),
              Chart (M, (r, s))): Change of coordinates from Chart (M, (x, y)) to Chart (M, (r, s))}

        """
        return self._coord_changes


    def is_manifestly_coordinate_domain(self):
        r"""
        Returns ``True`` if the manifold is known to be the domain of some
        coordinate chart and ``False`` otherwise.

        If ``False`` is returned, either the manifold cannot be the domain of
        some coordinate chart or no such chart has been declared yet.

        EXAMPLES::

            sage: TopManifold._clear_cache_()  # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart()
            sage: U.is_manifestly_coordinate_domain()
            True
            sage: M.is_manifestly_coordinate_domain()
            False
            sage: Y.<u,v> = M.chart()
            sage: M.is_manifestly_coordinate_domain()
            True

        """
        return not self._covering_charts == []

    def open_subset(self, name, latex_name=None, coord_def={}):
        r"""
        Create an open subset of the manifold.

        An open subset is a set that is (i) included in the manifold and (ii)
        open with respect to the manifold's topology.

        INPUT:

        - ``name`` -- name given to the open subset
        - ``latex_name`` --  (default: ``None``) LaTeX symbol to denote the
          subset; if none is provided, it is set to ``name``
        - ``coord_def`` -- (default: {}) definition of the subset in
          terms of coordinates; ``coord_def`` must a be dictionary with keys
          charts on ``self`` and values the symbolic expressions formed by the
          coordinates to define the subset.

        OUTPUT:

        - the open subset, as an instance of :class:`TopManifold`.

        EXAMPLES:

        Creating an open subset of a manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: A = M.open_subset('A'); A
            Open subset A of the 2-dimensional topological manifold M

        As an open subset of a topological manifold, A is itself a topological
        manifold, on the same topological field and of the same dimension a
        M::

            sage: isinstance(A, TopManifold)
            True
            sage: A.base_field() == M.base_field()
            True
            sage: dim(A) == dim(M)
            True

        Creating an open subset of A::

            sage: B = A.open_subset('B'); B
            Open subset B of the 2-dimensional topological manifold M

        We have then::

            sage: A.subsets()  # random (set output)
            {Open subset B of the 2-dimensional topological manifold M,
             Open subset A of the 2-dimensional topological manifold M}
            sage: B.is_subset(A)
            True
            sage: B.is_subset(M)
            True

        Defining an open subset by some coordinate restrictions: the open
        unit disk in `\RR^2`::

            sage: M = TopManifold(2, 'R^2')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: U = M.open_subset('U', coord_def={c_cart: x^2+y^2<1}); U
            Open subset U of the 2-dimensional topological manifold R^2

        Since the argument ``coord_def`` has been set, U is automatically
        provided with a chart, which is the restriction of the Cartesian one
        to U::

            sage: U.atlas()
            [Chart (U, (x, y))]

        Therefore, one can immediately check whether a point belongs to U::

            sage: M.point((0,0)) in U
            True
            sage: M.point((1/2,1/3)) in U
            True
            sage: M.point((1,2)) in U
            False

        """
        resu = TopManifold(self._dim, name, latex_name=latex_name,
                           field=self._field, start_index=self._sindex,
                           category=self.category(),
                           ambient_manifold=self._manifold)
        resu._supersets.update(self._supersets)
        for sd in self._supersets:
            sd._subsets.add(resu)
        self._top_subsets.add(resu)
        for chart, restrictions in coord_def.iteritems():
            if chart not in self._atlas:
                raise ValueError("The " + str(chart) + "does not belong to " +
                    "the atlas of " + str(self))
            chart.restrict(resu, restrictions)
        return resu

    def chart(self, coordinates='', names=None):
        r"""
        Define a chart the domain of which is the manifold.

        A *chart* is a pair `(U,\varphi)`, where `U` is the current manifold
        and `\varphi: U \rightarrow V \subset K^n`
        is a homeomorphism from `U` to an open subset `V` of `K^n`, `K` being
        the field on which the manifold containing the open set is defined.

        The components `(x^1,\ldots,x^n)` of `\varphi`, defined by
        `\varphi(p) = (x^1(p),\ldots,x^n(p))`, are called the *coordinates*
        of the chart `(U,\varphi)`.

        See :class:`~sage.manifolds.chart.Chart` for a complete
        documentation.

        INPUT:

        - ``coordinates`` -- single string defining the coordinate symbols and
          ranges: the coordinates are separated by ' ' (space) and each
          coordinate has at most three fields, separated by ':':

          1. The coordinate symbol (a letter or a few letters)
          2. (optional, only for manifolds over `\RR`) The interval `I`
             defining the coordinate range: if not
             provided, the coordinate is assumed to span all `\RR`; otherwise
             `I` must be provided in the form (a,b) (or equivalently ]a,b[)
             The bounds a and b can be +/-Infinity, Inf, infinity, inf or oo.
             For *singular* coordinates, non-open intervals such as [a,b] and
             (a,b] (or equivalently ]a,b]) are allowed.
             Note that the interval declaration must not contain any space
             character.
          3. (optional) The LaTeX spelling of the coordinate; if not provided the
             coordinate symbol given in the first field will be used.

          The order of the fields 2 and 3 does not matter and each of them can
          be omitted.
          If it contains any LaTeX expression, the string ``coordinates`` must
          be declared with the prefix 'r' (for "raw") to allow for a proper
          treatment of the backslash character (see examples below).
          If no interval range and no LaTeX spelling is to be provided for any
          coordinate, the argument ``coordinates`` can be omitted when the
          shortcut operator <,> is used via Sage preparser (see examples below)
        - ``names`` -- (default: ``None``) unused argument, except if
          ``coordinates`` is not provided; it must then be a tuple containing
          the coordinate symbols (this is guaranted if the shortcut operator <,>
          is used).

        OUTPUT:

        - the created chart, as an instance of
          :class:`~sage.manifolds.chart.Chart` or of the subclass
          :class:`~sage.manifolds.chart.RealChart` for manifolds over `\RR`.

        EXAMPLES:

        Chart on a 2-dimensional manifold::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X = U.chart('x y'); X
            Chart (U, (x, y))
            sage: X[0]
            x
            sage: X[1]
            y
            sage: X[:]
            (x, y)

        The declared coordinates are not known at the global level::

            sage: y
            Traceback (most recent call last):
            ...
            NameError: name 'y' is not defined

        They can be recovered by the operator ``[:]`` applied to the chart::

            sage: (x, y) = X[:]
            sage: y
            y
            sage: type(y)
            <type 'sage.symbolic.expression.Expression'>

        But a shorter way to proceed is to use the operator ``<,>`` in the
        left-hand side of the chart declaration (there is then no need to
        pass the string 'x y' to chart())::

            sage: TopManifold._clear_cache_() # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: U = M.open_subset('U')
            sage: X.<x,y> = U.chart(); X
            Chart (U, (x, y))

        Indeed, the declared coordinates are then known at the global level::

            sage: y
            y
            sage: (x,y) == X[:]
            True

        Actually the instruction ``X.<x,y> = U.chart()`` is
        equivalent to the combination of the two instructions
        ``X = U.chart('x y')`` and ``(x,y) = X[:]``.

        See the documentation of class
        :class:`~sage.manifolds.chart.Chart` for more examples,
        especially regarding the coordinates ranges and restrictions.

        """
        from sage.manifolds.chart import Chart, RealChart
        if self._field == 'real':
            return RealChart(self, coordinates=coordinates, names=names)
        return Chart(self, coordinates=coordinates, names=names)

    #~ def scalar_field_algebra(self):
        #~ r"""
        #~ Returns the algebra of scalar fields defined on ``self``.
#~
        #~ See :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
        #~ for a complete documentation.
#~
        #~ OUTPUT:
#~
        #~ - instance of
          #~ :class:`~sage.manifolds.scalarfield_algebra.ScalarFieldAlgebra`
          #~ representing the algebra `C^\infty(U)` of all scalar fields defined
          #~ on `U` = ``self``.
#~
        #~ EXAMPLE:
#~
        #~ Scalar algebra of a 3-dimensional open subset::
#~
            #~ sage: TopManifold._clear_cache_() # for doctests only
            #~ sage: M = TopManifold(3, 'M')
            #~ sage: U = M.open_subset('U')
            #~ sage: CU = U.scalar_field_algebra(); CU
            #~ algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            #~ sage: CU.category()
            #~ Category of commutative algebras over Symbolic Ring
            #~ sage: CU.zero()
            #~ scalar field 'zero' on the open subset 'U' of the 3-dimensional manifold 'M'
#~
        #~ """
        #~ return self._scalar_field_algebra
#~
    #~ def _Hom_(self, other, category=None):
        #~ r"""
        #~ Construct the set of morphisms (i.e. continuous maps)
        #~ ``self`` --> ``other``.
#~
        #~ INPUT:
#~
        #~ - ``other`` -- an open subset of some manifold
        #~ - ``category`` -- (default: ``None``) not used here (to ensure
          #~ compatibility with generic hook ``_Hom_``)
#~
        #~ OUTPUT:
#~
        #~ - the homset Hom(U,V), where U is ``self`` and V is ``other``
#~
        #~ See class
        #~ :class:`~sage.manifolds.manifold_homset.ManifoldHomset`
        #~ for more documentation.
#~
        #~ """
        #~ from sage.manifolds.manifold_homset import TopManifoldHomset
        #~ return TopManifoldHomset(self, other)
#~
    #~ def scalar_field(self, coord_expression=None, chart=None, name=None,
                     #~ latex_name=None):
        #~ r"""
        #~ Define a scalar field on the open set.
#~
        #~ See :class:`~sage.manifolds.scalarfield.ScalarField` for a
        #~ complete documentation.
#~
        #~ INPUT:
#~
        #~ - ``coord_expression`` -- (default: ``None``) coordinate expression(s)
          #~ of the scalar field; this can be either
#~
          #~ - a single coordinate expression; if the argument ``chart`` is
            #~ ``'all'``, this expression is set to all the charts defined
            #~ on the open set; otherwise, the expression is set in the
            #~ specific chart provided by the argument ``chart``
          #~ - a dictionary of coordinate expressions, with the charts as keys.
#~
          #~ If ``coord_expression`` is ``None`` or does not fully specified the
          #~ scalar field, other coordinate expressions can be added subsequently
          #~ by means of the methods
          #~ :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr`,
          #~ :meth:`~sage.manifolds.scalarfield.ScalarField.add_expr_by_continuation`,
          #~ or :meth:`~sage.manifolds.scalarfield.ScalarField.set_expr`
        #~ - ``chart`` -- (default: ``None``) chart defining the coordinates used
          #~ in ``coord_expression`` when the latter is a single coordinate
          #~ expression; if none is provided (default), the default chart of the
          #~ open set is assumed. If ``chart=='all'``, ``coord_expression`` is
          #~ assumed to be independent of the chart (constant scalar field).
        #~ - ``name`` -- (default: ``None``) name given to the scalar field
        #~ - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          #~ field; if none is provided, the LaTeX symbol is set to ``name``
#~
        #~ OUTPUT:
#~
        #~ - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          #~ representing the defined scalar field.
#~
        #~ EXAMPLES:
#~
        #~ A scalar field defined by its coordinate expression in the open
        #~ set's default chart::
#~
            #~ sage: TopManifold._clear_cache_() # for doctests only
            #~ sage: M = TopManifold(3, 'M')
            #~ sage: U = M.open_subset('U')
            #~ sage: c_xyz.<x,y,z> = U.chart()
            #~ sage: f = U.scalar_field(sin(x)*cos(y) + z, name='F'); f
            #~ scalar field 'F' on the open subset 'U' of the 3-dimensional manifold 'M'
            #~ sage: f.display()
            #~ F: U --> R
               #~ (x, y, z) |--> cos(y)*sin(x) + z
            #~ sage: f.parent()
            #~ algebra of scalar fields on the open subset 'U' of the 3-dimensional manifold 'M'
            #~ sage: f in U.scalar_field_algebra()
            #~ True
#~
        #~ Equivalent definition with the chart specified::
#~
            #~ sage: f = U.scalar_field(sin(x)*cos(y) + z, chart=c_xyz, name='F')
            #~ sage: f.display()
            #~ F: U --> R
               #~ (x, y, z) |--> cos(y)*sin(x) + z
#~
        #~ Equivalent definition with a dictionary of coordinate expression(s)::
#~
            #~ sage: f = U.scalar_field({c_xyz: sin(x)*cos(y) + z}, name='F')
            #~ sage: f.display()
            #~ F: U --> R
               #~ (x, y, z) |--> cos(y)*sin(x) + z
#~
        #~ See the documentation of class
        #~ :class:`~sage.manifolds.scalarfield.ScalarField` for more
        #~ examples.
#~
        #~ .. SEEALSO::
#~
            #~ :meth:`constant_scalar_field`, :meth:`zero_scalar_field`
#~
        #~ """
        #~ if isinstance(coord_expression, dict):
            #~ # check validity of entry
            #~ for chart in coord_expression:
                #~ if not chart._domain.is_subset(self):
                    #~ raise ValueError("The " + str(chart) + " is not defined " +
                                     #~ "on some subset of the " + str(self))
        #~ elif coord_expression is not None and chart != 'all':
            #~ # coord_expression is valid only in a specific chart
            #~ if chart is None:
                #~ chart = self._def_chart
            #~ coord_expression = {chart: coord_expression}
        #~ return self.scalar_field_algebra()._element_constructor_(
                                            #~ coord_expression=coord_expression,
                                            #~ name=name, latex_name=latex_name)
#~
    #~ def constant_scalar_field(self, value, name=None, latex_name=None):
        #~ r"""
        #~ Define a constant scalar field on the open set.
#~
        #~ INPUT:
#~
        #~ - ``value`` -- constant value of the scalar field, either a numerical
          #~ value or a symbolic expression not involving any chart coordinates
        #~ - ``name`` -- (default: ``None``) name given to the scalar field
        #~ - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the scalar
          #~ field; if none is provided, the LaTeX symbol is set to ``name``
#~
        #~ OUTPUT:
#~
        #~ - instance of :class:`~sage.manifolds.scalarfield.ScalarField`
          #~ representing the scalar field whose constant value is ``value``
#~
        #~ EXAMPLES:
#~
        #~ A constant scalar field on the 2-sphere::
#~
            #~ sage: M = TopManifold(2, 'M') # the 2-dimensional sphere S^2
            #~ sage: U = M.open_subset('U') # complement of the North pole
            #~ sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
            #~ sage: V = M.open_subset('V') # complement of the South pole
            #~ sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
            #~ sage: M.declare_union(U,V)   # S^2 is the union of U and V
            #~ sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)),
            #~ ....:                                intersection_name='W',
            #~ ....:                                restrictions1= x^2+y^2!=0,
            #~ ....:                                restrictions2= u^2+v^2!=0)
            #~ sage: uv_to_xy = xy_to_uv.inverse()
            #~ sage: f = M.constant_scalar_field(1); f
            #~ scalar field on the 2-dimensional manifold 'M'
            #~ sage: f.display()
            #~ M --> R
            #~ on U: (x, y) |--> 1
            #~ on V: (u, v) |--> 1
#~
        #~ We have::
#~
            #~ sage: f.restrict(U) == U.constant_scalar_field(1)
            #~ True
            #~ sage: M.constant_scalar_field(0) is M.zero_scalar_field()
            #~ True
#~
        #~ .. SEEALSO::
#~
            #~ :meth:`zero_scalar_field`
        #~ """
        #~ return self.scalar_field_algebra()._element_constructor_(
                                              #~ coord_expression=value,
                                              #~ name=name, latex_name=latex_name)
#~
    #~ def zero_scalar_field(self):
        #~ r"""
        #~ Return the zero scalar field defined on the open set.
#~
        #~ EXAMPLE::
#~
            #~ sage: TopManifold._clear_cache_() # for doctests only
            #~ sage: M = TopManifold(2, 'M')
            #~ sage: X.<x,y> = M.chart()
            #~ sage: f = M.zero_scalar_field(); f
            #~ scalar field 'zero' on the 2-dimensional manifold 'M'
            #~ sage: f.display()
            #~ zero: M --> R
               #~ (x, y) |--> 0
            #~ sage: f.parent()
            #~ algebra of scalar fields on the 2-dimensional manifold 'M'
            #~ sage: f is M.scalar_field_algebra().zero()
            #~ True
#~
        #~ """
        #~ return self._zero_scalar_field
#~
#~
    #~ def curve(self, coord_expression, param, chart=None, name=None,
              #~ latex_name=None):
        #~ r"""
        #~ Define a curve in ``self``.
#~
        #~ See :class:`~sage.manifolds.curve.ManifoldCurve` for details.
#~
        #~ INPUT:
#~
        #~ - ``coord_expression`` -- either
#~
          #~ - (i) a dictionary whose keys are charts on ``self`` and values
            #~ the coordinate expressions (as lists or tuples) of the curve in
            #~ the given chart
          #~ - (ii) a single coordinate expression in a given chart on ``self``,
            #~ the latter being provided by the argument ``chart``
#~
          #~ In both cases, if the dimension of the arrival manifold is 1,
          #~ a single coordinate expression can be passed instead of a tuple with
          #~ a single element
        #~ - ``param`` -- a tuple of the type ``(t, t_min, t_max)``, where ``t``
          #~ is the curve parameter used in ``coord_expression``, ``t_min`` is its
          #~ minimal value and ``t_max`` its maximal value; if ``t_min=-Infinity``
          #~ and ``t_max=+Infinity``, they can be omitted and ``t`` can be passed
          #~ for ``param``, instead of the tuple ``(t, t_min, t_max)``
        #~ - ``chart`` -- (default: ``None``) chart on ``self`` used for case (ii)
          #~ above; if ``None`` the default chart of ``self`` is assumed.
        #~ - ``name`` -- (default: ``None``) string; symbol given to the curve
        #~ - ``latex_name`` -- (default: ``None``) string; LaTeX symbol to denote the
          #~ the curve; if none is provided, ``name`` will be used
#~
        #~ OUTPUT:
#~
        #~ - instance of
          #~ :class:`~sage.manifolds.curve.ManifoldCurve`
#~
        #~ EXAMPLES:
#~
        #~ The lemniscate of Gerono in the 2-dimensional Euclidean plane::
#~
            #~ sage: M = TopManifold(2, 'M')
            #~ sage: X.<x,y> = M.chart()
            #~ sage: R.<t> = RealLine()
            #~ sage: c = M.curve([sin(t), sin(2*t)/2], (t, 0, 2*pi), name='c'); c
            #~ Curve 'c' in the 2-dimensional manifold 'M'
#~
        #~ The same definition with the coordinate expression passed as a
        #~ dictionary::
#~
            #~ sage: c = M.curve({X: [sin(t), sin(2*t)/2]}, (t, 0, 2*pi), name='c'); c
            #~ Curve 'c' in the 2-dimensional manifold 'M'
#~
        #~ An example of definition with ``t_min`` and ``t_max`` omitted: a helix
        #~ in `\RR^3`::
#~
            #~ sage: R3 = TopManifold(3, 'R^3')
            #~ sage: X.<x,y,z> = R3.chart()
            #~ sage: c = R3.curve([cos(t), sin(t), t], t, name='c'); c
            #~ Curve 'c' in the 3-dimensional manifold 'R^3'
            #~ sage: c.domain() # check that t is unbounded
            #~ field R of real numbers
#~
        #~ See the documentation of
        #~ :class:`~sage.manifolds.curve.ManifoldCurve` for more
        #~ examples.
#~
        #~ """
        #~ from sage.rings.infinity import Infinity
        #~ from sage.manifolds.smooth.manifold import RealLine
        #~ if not isinstance(param, (tuple, list)):
            #~ param = (param, -Infinity, Infinity)
        #~ elif len(param) != 3:
            #~ raise TypeError("the argument 'param' must be of the type " +
                            #~ "(t, t_min, t_max)")
        #~ t = param[0]
        #~ t_min = param[1]
        #~ t_max = param[2]
        #~ real_field = RealLine(names=(repr(t),))
        #~ interval = real_field.open_interval(t_min, t_max)
        #~ curve_set = Hom(interval, self)
        #~ if not isinstance(coord_expression, dict):
            #~ # Turn coord_expression into a dictionary:
            #~ if chart is None:
                #~ chart = self._def_chart
            #~ elif chart not in self._atlas:
                #~ raise ValueError("the {} has not been".format(chart) +
                                     #~ " defined on the {}".format(self))
            #~ if isinstance(coord_expression, (tuple, list)):
                #~ coord_expression = {chart: coord_expression}
            #~ else:
                #~ # case self.dim()=1
                #~ coord_expression = {chart: (coord_expression,)}
        #~ return curve_set(coord_expression, name=name, latex_name=latex_name)
#~
    #~ def continuous_mapping(self, codomain, coord_functions=None, chart1=None,
                     #~ chart2=None, name=None, latex_name=None):
        #~ r"""
        #~ Define a continuous mapping between ``self`` and another
        #~ subset (possibly on another manifold).
#~
        #~ See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        #~ complete documentation.
#~
        #~ INPUT:
#~
        #~ - ``codomain`` -- mapping's codomain (the arrival manifold or some
          #~ subset of it)
        #~ - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          #~ either
#~
          #~ - (i) a dictionary of
            #~ the coordinate expressions (as lists (or tuples) of the
            #~ coordinates of the image expressed in terms of the coordinates of
            #~ the considered point) with the pairs of charts (chart1, chart2)
            #~ as keys (chart1 being a chart on ``self`` and chart2 a chart on
            #~ ``codomain``)
          #~ - (ii) a single coordinate expression in a given pair of charts, the
            #~ latter being provided by the arguments ``chart1`` and ``chart2``
#~
          #~ In both cases, if the dimension of the arrival manifold is 1,
          #~ a single coordinate expression can be passed instead of a tuple with
          #~ a single element
        #~ - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          #~ on ``self`` defining the start coordinates involved in
          #~ ``coord_functions`` for case (ii); if none is provided, the
          #~ coordinates are assumed to refer to the default chart of ``self``
        #~ - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          #~ on ``codomain`` defining the arrival coordinates involved in
          #~ ``coord_functions`` for case (ii); if none is provided, the
          #~ coordinates are assumed to refer to the default chart of ``codomain``
        #~ - ``name`` -- (default: ``None``) name given to the differentiable
          #~ mapping
        #~ - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          #~ differentiable mapping; if none is provided, the LaTeX symbol is set to
          #~ ``name``
#~
        #~ OUTPUT:
#~
        #~ - the differentiable mapping, as an instance of
          #~ :class:`~sage.manifolds.diffmapping.DiffMapping`
#~
        #~ EXAMPLES:
#~
        #~ A mapping between an open subset of `S^2` covered by regular spherical
        #~ coordinates and `\RR^3`::
#~
            #~ sage: TopManifold._clear_cache_() # for doctests only
            #~ sage: M = TopManifold(2, 'S^2')
            #~ sage: U = M.open_subset('U')
            #~ sage: c_spher.<th,ph> = U.chart(r'th:(0,pi):\theta ph:(0,2*pi):\phi')
            #~ sage: N = TopManifold(3, 'R^3', r'\RR^3')
            #~ sage: c_cart.<x,y,z> = N.chart()  # Cartesian coord. on R^3
            #~ sage: Phi = U.diff_mapping(N, (sin(th)*cos(ph), sin(th)*sin(ph), cos(th)),
            #~ ....:                      name='Phi', latex_name=r'\Phi'); Phi
            #~ differentiable mapping 'Phi' from the open subset 'U' of the
             #~ 2-dimensional manifold 'S^2' to the 3-dimensional manifold 'R^3'
#~
        #~ The same definition, but with a dictionary with pairs of charts as
        #~ keys (case (i) above)::
#~
            #~ sage: Phi1 = U.diff_mapping(N,
            #~ ....:        {(c_spher, c_cart): (sin(th)*cos(ph), sin(th)*sin(ph), cos(th))},
            #~ ....:        name='Phi', latex_name=r'\Phi')
            #~ sage: Phi1 == Phi
            #~ True
#~
        #~ The differentiable mapping acting on a point::
#~
            #~ sage: p = U.point((pi/2, pi)); p
            #~ point on 2-dimensional manifold 'S^2'
            #~ sage: Phi(p)
            #~ point on 3-dimensional manifold 'R^3'
            #~ sage: Phi(p).coord(c_cart)
            #~ (-1, 0, 0)
            #~ sage: Phi1(p) == Phi(p)
            #~ True
#~
        #~ See the documentation of class
        #~ :class:`~sage.manifolds.diffmapping.DiffMapping` for more
        #~ examples.
#~
        #~ """
        #~ homset = Hom(self, codomain)
        #~ if coord_functions is None:
            #~ coord_functions = {}
        #~ if not isinstance(coord_functions, dict):
            #~ # Turn coord_functions into a dictionary:
            #~ if chart1 is None:
                #~ chart1 = self._def_chart
            #~ elif chart1 not in self._atlas:
                #~ raise ValueError("{} is not a chart ".format(chart1) +
                                 #~ "defined on the {}".format(self))
            #~ if chart2 is None:
                #~ chart2 = codomain._def_chart
            #~ elif chart2 not in codomain._atlas:
                #~ raise ValueError("{} is not a chart ".format(chart2) +
                                 #~ " defined on the {}".format(codomain))
            #~ coord_functions = {(chart1, chart2): coord_functions}
        #~ return homset(coord_functions, name=name, latex_name=latex_name)
#~
    #~ def homeomorphism(self, codomain, coord_functions=None, chart1=None,
                       #~ chart2=None, name=None, latex_name=None):
        #~ r"""
        #~ Define a homeomorphism between ``self`` and another open subset
        #~ (possibly on another manifold).
#~
        #~ See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        #~ complete documentation.
#~
        #~ INPUT:
#~
        #~ - ``codomain`` -- mapping's codomain (the arrival manifold or some
          #~ subset of it)
        #~ - ``coord_functions`` -- (default: ``None``) if not ``None``, must be
          #~ either
#~
          #~ - (i) a dictionary of
            #~ the coordinate expressions (as lists (or tuples) of the
            #~ coordinates of the image expressed in terms of the coordinates of
            #~ the considered point) with the pairs of charts (chart1, chart2)
            #~ as keys (chart1 being a chart on ``self`` and chart2 a chart on
            #~ ``codomain``)
          #~ - (ii) a single coordinate expression in a given pair of charts, the
            #~ latter being provided by the arguments ``chart1`` and ``chart2``
#~
          #~ In both cases, if the dimension of the arrival manifold is 1,
          #~ a single coordinate expression can be passed instead of a tuple with
          #~ a single element
        #~ - ``chart1`` -- (default: ``None``; used only in case (ii) above) chart
          #~ on ``self`` defining the start coordinates involved in
          #~ ``coord_functions`` for case (ii); if none is provided, the
          #~ coordinates are assumed to refer to the default chart of ``self``
        #~ - ``chart2`` -- (default: ``None``; used only in case (ii) above) chart
          #~ on ``codomain`` defining the arrival coordinates involved in
          #~ ``coord_functions`` for case (ii); if none is provided, the
          #~ coordinates are assumed to refer to the default chart of ``codomain``
        #~ - ``name`` -- (default: ``None``) name given to the diffeomorphism
        #~ - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the
          #~ diffeomorphism; if none is provided, the LaTeX symbol is set to
          #~ ``name``
#~
        #~ OUTPUT:
#~
        #~ - the diffeomorphism, as an instance of
          #~ :class:`~sage.manifolds.diffmapping.DiffMapping`
#~
        #~ EXAMPLE:
#~
        #~ A diffeomorphism between two 2-dimensional subsets::
#~
            #~ sage: TopManifold._clear_cache_() # for doctests only
            #~ sage: M = TopManifold(2, 'M', r'{\cal M}')
            #~ sage: U = M.open_subset('U')
            #~ sage: c_xv.<x,y> = U.chart(r'x:(-pi/2,+oo) y:(-pi/2,+oo)')
            #~ sage: N = TopManifold(2, 'N', r'{\cal N}')
            #~ sage: V = N.open_subset('V')
            #~ sage: c_zt.<z,t> = V.chart(r'z t')
            #~ sage: Phi = U.diffeomorphism(V, (arctan(x), arctan(y)), name='Phi',
            #~ ....:                        latex_name=r'\Phi')
#~
        #~ See the documentation of class
        #~ :class:`~sage.manifolds.diffmapping.DiffMapping` for more
        #~ examples.
#~
        #~ """
        #~ homset = Hom(self, codomain)
        #~ if coord_functions is None:
            #~ coord_functions = {}
        #~ if not isinstance(coord_functions, dict):
            #~ # Turn coord_functions into a dictionary:
            #~ if chart1 is None:
                #~ chart1 = self._def_chart
            #~ elif chart1 not in self._atlas:
                #~ raise ValueError("{} is not a chart ".format(chart1) +
                                 #~ "defined on the {}".format(self))
            #~ if chart2 is None:
                #~ chart2 = codomain._def_chart
            #~ elif chart2 not in codomain._atlas:
                #~ raise ValueError("{} is not a chart ".format(chart2) +
                                 #~ " defined on the {}".format(codomain))
            #~ coord_functions = {(chart1, chart2): coord_functions}
        #~ return homset(coord_functions, name=name, latex_name=latex_name,
                      #~ is_diffeomorphism=True)
#~
    #~ def identity_map(self):
        #~ r"""
        #~ Identity map on ``self``.
#~
        #~ See :class:`~sage.manifolds.diffmapping.DiffMapping` for a
        #~ complete documentation.
#~
        #~ OUTPUT:
#~
        #~ - the identity map, as an instance of
          #~ :class:`~sage.manifolds.diffmapping.DiffMapping`
#~
        #~ """
        #~ return self._identity_map

