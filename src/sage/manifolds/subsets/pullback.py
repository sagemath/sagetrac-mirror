r"""
Manifold Subsets Defined as Pullbacks of Subsets under Continuous Maps
"""


# ****************************************************************************
#       Copyright (C) 2021 Matthias Koeppe <mkoeppe@math.ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.sets_cat import Sets
from sage.categories.modules_with_basis import ModulesWithBasis as FreeModules
from sage.categories.metric_spaces import MetricSpaces
from sage.modules.free_module import is_FreeModule
from sage.rings.infinity import infinity, minus_infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.manifolds.subset import ManifoldSubset
from sage.manifolds.continuous_map import ContinuousMap
from sage.manifolds.chart import Chart
from sage.manifolds.scalarfield import ScalarField
from sage.sets.real_set import RealSet
from sage.geometry.polyhedron.base import is_Polyhedron

class ManifoldSubsetPullback(ManifoldSubset):

    """
    Manifold subset defined as a pullback of a subset under a continuous map.

    INPUT:

    - ``map`` - an instance of :class:`ContinuousMap` or
      :class:`ScalarField` or :class:`Chart`

    - ``codomain_subset`` - an instance of :class:`ManifoldSubset`,
      :class:`RealSet`, :class:`Polyhedron_base`,
      :class:`C_Polyhedron`, :class:`NNC_Polyhedron`

    EXAMPLES::

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: M = Manifold(2, 'R^2', structure='topological')
        sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

    Pulling back a real interval under a scalar field::

        sage: r_squared = M.scalar_field(x^2+y^2)
        sage: r_squared.set_immutable()
        sage: cl_I = RealSet((1, 4)); cl_I
        (1, 4)
        sage: cl_O = ManifoldSubsetPullback(r_squared, None, I); cl_O
        Subset f_inv_I of the 2-dimensional topological manifold R^2
        sage: M.point((0, 0)) in cl_O
        True
        sage: M.point((1, 1)) in cl_O
        False

    Pulling back an open real interval gives an open subset::

        sage: I = RealSet((1, 4)); I
        (1, 4)
        sage: O = ManifoldSubsetPullback(r_squared, None, I); O
        Open subset f_inv_(1, 4) of the 2-dimensional topological manifold R^2
        sage: M.point((1, 0)) in O
        False
        sage: M.point((1, 1)) in O
        True

    Pulling back a polytope under a chart::

        sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        sage: S = ManifoldSubsetPullback(c_cart, None, P); S
        Subset x_y_inv_P of the 2-dimensional topological manifold R^2
        sage: M((1, 2)) in S
        True
        sage: M((2, 0)) in S
        False

    Using the embedding map of a submanifold::

        sage: M = Manifold(3, 'M', structure="topological")
        sage: N = Manifold(2, 'N', ambient=M, structure="topological")
        sage: N
        2-dimensional topological submanifold N immersed in the 3-dimensional topological manifold M
        sage: CM.<x,y,z> = M.chart()
        sage: CN.<u,v> = N.chart()
        sage: t = var('t')
        sage: phi = N.continuous_map(M, {(CN,CM): [u,v,t+u^2+v^2]})
        sage: phi_inv = M.continuous_map(N, {(CM,CN): [x,y]})
        sage: phi_inv_t = M.scalar_field({CM: z-x^2-y^2})
        sage: N.set_immersion(phi, inverse=phi_inv, var=t,
        ....:                 t_inverse={t: phi_inv_t})
        sage: N.declare_embedding()

        sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
        sage: S = M.open_subset('S', coord_def={CM: z<1})
        sage: phi_without_t = N.continuous_map(M, {(CN, CM): [expr.subs(t=0) for expr in phi.expr()]}); phi_without_t
        Continuous map
         from the 2-dimensional topological submanifold N
          embedded in the 3-dimensional topological manifold M
         to the 3-dimensional topological manifold M
        sage: phi_without_t.expr()
        (u, v, u^2 + v^2)
        sage: D = ManifoldSubsetPullback(phi_without_t, codomain_subset=S); D
        Subset f_inv_S of the 2-dimensional topological submanifold N embedded in the 3-dimensional topological manifold M
        sage: N.point((2,0)) in D
        False

    """
    @staticmethod
    def __classcall_private__(cls, map, inverse=None, codomain_subset=None,
                              name=None, latex_name=None):
        try:
            is_mutable = map.is_mutable()
        except AttributeError:
            pass
        else:
            if is_mutable:
                map = map.copy()
                map.set_immutable()
        try:
            is_mutable = inverse.is_mutable()
        except AttributeError:
            pass
        else:
            if is_mutable:
                inverse = inverse.copy()
                inverse.set_immutable()

        if codomain_subset is None:
            try:
                codomain_subset = map.codomain()
            except AttributeError:
                if isinstance(codomain_subset, Chart):
                    codomain_subset = FreeModule(self.base_field(), map.domain().dimension())

        if inverse is None:
            if isinstance(map, Chart):
                from sage.misc.latex import latex
                inverse_latex_name = '(' + ','.join(str(latex(x)) + '^{-1}' for x in map) + ')'
                inverse_name = '_'.join(repr(x) for x in map) + '_inv'
            else:
                map_name = map._name or 'f'
                map_latex_name = map._latex_name or map_name
                inverse_name = map_name + '_inv'
                inverse_latex_name = map_latex_name + r'^{-1}'
        else:
            inverse_name = inverse._name
            inverse_latex_name = inverse._latex_name
        try:
            codomain_subset_latex_name = codomain_subset._latex_name
            codomain_subset_name = codomain_subset._name
        except AttributeError:
            from sage.misc.latex import latex
            codomain_subset_latex_name = str(latex(codomain_subset))
            s = repr(codomain_subset)
            if len(s) > 10:
                codomain_subset_name = 'P'
            else:
                codomain_subset_name = s
        if latex_name is None:
            if name is None:
                latex_name = inverse_latex_name + '(' + codomain_subset_latex_name + ')'
            else:
                latex_name = name
        if name is None:
            name = inverse_name + '_' + codomain_subset_name

        if cls._is_open(codomain_subset):

            try:
                coord_def = cls._coord_def(map, codomain_subset)
            except NotImplementedError:
                pass
            else:
                return map.domain().open_subset(name=name, latex_name=latex_name,
                                                coord_def=coord_def)

        self = super().__classcall__(cls, map, inverse, codomain_subset, name, latex_name)

        return self

    @staticmethod
    def _is_open(codomain_subset):
        """
        Return whether ``codomain_subset`` is (known to be) an open subset of its ambient space.

        EXAMPLES:

        Manifolds and subsets::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: R2 = Manifold(2, 'R^2', structure='topological'); R2
            2-dimensional topological manifold R^2
            sage: ManifoldSubsetPullback._is_open(R2)
            True
            sage: A = R2.subset('A'); A
            Subset A of the 2-dimensional topological manifold R^2
            sage: ManifoldSubsetPullback._is_open(A)
            False

        :class:`RealSet` instances::

            sage: I = RealSet.open(1, 2); I
            (1, 2)
            sage: ManifoldSubsetPullback._is_open(I)
            True
            sage: cl_I = RealSet.closed(1, 2); cl_I
            [1, 2]
            sage: ManifoldSubsetPullback._is_open(cl_I)
            False

        Polyhedra::

            sage: Empty = Polyhedron(ambient_dim=2); Empty
            The empty polyhedron in ZZ^2
            sage: ManifoldSubsetPullback._is_open(Empty)
            True
            sage: C = polytopes.cube(); C
            A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 8 vertices
            sage: ManifoldSubsetPullback._is_open(C)
            False

        PPL polyhedra and not-necessarily-closed polyhedra::

            sage: from ppl import Variable, C_Polyhedron, NNC_Polyhedron, Constraint_System
            sage: u = Variable(0)
            sage: v = Variable(1)
            sage: CS = Constraint_System()
            sage: CS.insert(0 < u)
            sage: CS.insert(u < 1)
            sage: CS.insert(0 < v)
            sage: CS.insert(v < 1)
            sage: CS.insert(u + v <= 3)       # redundant inequality
            sage: P = NNC_Polyhedron(CS); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 4 closure_points
            sage: ManifoldSubsetPullback._is_open(P)
            True
            sage: CS.insert(u + v <= 1)
            sage: T = NNC_Polyhedron(CS); T
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 1 point, 3 closure_points
            sage: ManifoldSubsetPullback._is_open(T)
            False

        """

        if isinstance(codomain_subset, ManifoldSubset):
            return codomain_subset.is_open()

        if isinstance(codomain_subset, RealSet):
            return codomain_subset.is_open()

        if is_Polyhedron(codomain_subset):
            return codomain_subset.is_empty() or codomain_subset.is_universe()

        if codomain_subset in Sets().Finite():
            return codomain_subset.cardinality() == 0

        if hasattr(codomain_subset, 'minimized_constraints'):
            try:
                from ppl import NNC_Polyhedron, C_Polyhedron
            except ImportError:
                pass
            else:
                if isinstance(codomain_subset, (NNC_Polyhedron, C_Polyhedron)):
                    cs = codomain_subset.minimized_constraints()
                    if cs.has_equalities():
                        return False
                    if any(constraint.is_nonstrict_inequality()
                           for constraint in cs):
                        return False
                    return True

        return False

    @staticmethod
    def _interval_restriction(expr, interval):
        """
        Return a restriction expressing that ``expr`` lies in ``interval``.

        INPUT:

        - ``expr`` -- a symbolic expression
        - ``interval`` -- an instance of :class:`~sage.sets.real_set.InternalRealInterval`

        OUTPUT:

        - A restriction suitable as input to :meth:`~sage.manifolds.chart.restrict`:
          lists are conjunctions, tuples are disjunctions

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _interval_restriction = ManifoldSubsetPullback._interval_restriction
            sage: var('t')
            t
            sage: assume(t >= -2)
            sage: assume(t <= 5)
            sage: _interval_restriction(t, RealSet(3, 4)[0])
            [t > 3, t < 4]
            sage: _interval_restriction(t, RealSet.unbounded_below_closed(2)[0])
            t <= 2
            sage: _interval_restriction(t, RealSet.closed(-5, 5)[0])
            []
            sage: _interval_restriction(t, RealSet.unbounded_below_closed(-5)[0])
            ()
            sage: _interval_restriction(t, RealSet.unbounded_above_closed(6)[0])
            ()
            sage: _interval_restriction(t^2, RealSet.unbounded_above_closed(0)[0])
            []

        """

        conjunction = []
        if interval.lower() != minus_infinity:
            if interval.lower_closed():
                condition = (expr >= interval.lower())
                negation  = (expr <  interval.lower())
            else:
                condition = (expr >  interval.lower())
                negation  = (expr <= interval.lower())
            if negation:
                # known to be false
                return ()
            if not condition:
                # not known to be true
                conjunction.append(condition)

        if interval.upper() != infinity:
            if interval.upper_closed():
                condition = (expr <= interval.upper())
                negation  = (expr >  interval.upper())
            else:
                condition = (expr <  interval.upper())
                negation  = (expr >= interval.upper())
            if negation:
                # known to be false
                return ()
            if not condition:
                # not known to be true
                conjunction.append(condition)

        if len(conjunction) == 1:
            return conjunction[0]
        else:
            # lists express 'and'
            return conjunction

    @staticmethod
    def _realset_restriction(expr, realset):
        """
        Return a restriction expressing that ``expr`` lies in ``realset``.

        INPUT:

        - ``expr`` -- a symbolic expression
        - ``interval`` -- an instance of :class:`~sage.sets.real_set.RealSet`

        OUTPUT:

        - A restriction suitable as input to :meth:`~sage.manifolds.chart.restrict`:
          lists are conjunctions, tuples are disjunctions

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: _realset_restriction = ManifoldSubsetPullback._realset_restriction
            sage: var('t')
            t
            sage: assume(t >= -2)
            sage: assume(t <= 5)
            sage: _realset_restriction(t, RealSet(-oo, oo))
            []
            sage: _realset_restriction(t, RealSet())
            ()
            sage: _realset_restriction(t, RealSet([-5, -4], (-1, 1), [3, 4], [6, 7]))
            ([t > -1, t < 1], [t >= 3, t <= 4])

        """
        disjunction = []
        for interval in realset:
            condition = ManifoldSubsetPullback._interval_restriction(expr, interval)
            if condition == []:
                return []
            if condition != ():
                disjunction.append(condition)

        if len(disjunction) == 1:
            return disjunction[0]
        else:
            # tuples express 'or'
            return tuple(disjunction)

    @staticmethod
    def _coord_def(map, codomain_subset):

        #if isinstance(map, ContinuousMap) and isinstance(codomain_subset, Manifold):

        if isinstance(map, ScalarField) and isinstance(codomain_subset, RealSet):

            return {chart: ManifoldSubsetPullback._realset_restriction(func.expr(),
                                                                       codomain_subset)
                    for chart, func in map._express.items()}

        raise NotImplementedError


    def __init__(self, map, inverse, codomain_subset, name, latex_name):
        r"""
        Construct a manifold subset that is a pullback.

        """
        self._map = map
        self._inverse = inverse
        self._codomain_subset = codomain_subset
        base_manifold = map.domain()
        ManifoldSubset.__init__(self, base_manifold, name, latex_name=latex_name)

    def __contains__(self, point):
        r"""
        Check whether ``point`` is contained in ``self``.

        """
        if super().__contains__(point):
            return True
        return self._map(point) in self._codomain_subset

    def is_open(self):
        """
        Return if ``self`` is an open set.

        """
        # Because the map is continuous, the pullback is open if the
        # codomain_subset is open.  But because other code assumes
        # that open subsets are instances of Manifold, we do not use this
        # fact here. Instead, the constructor is responsible for creating
        # an instance of the appropriate subclass.
        return super().is_open()

    def is_closed(self):
        """
        Return if ``self`` is (known to be) a closed subset of the manifold.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2

        The pullback of a closed real interval under a scalar field is closed::

            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: cl_I = RealSet([1, 2]); cl_I
            [1, 2]
            sage: cl_O = ManifoldSubsetPullback(r_squared, None, cl_I); cl_O
            Subset f_inv_[1, 2] of the 2-dimensional topological manifold R^2
            sage: cl_O.is_closed()
            True

        The pullback of a (closed convex) polyhedron under a chart is closed::

            sage: P = Polyhedron(vertices=[[0, 0], [1, 2], [3, 4]]); P
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
            sage: McP = ManifoldSubsetPullback(c_cart, None, P, name='McP'); McP
            Subset McP of the 2-dimensional topological manifold R^2
            sage: McP.is_closed()
            True

            sage: from ppl import Variable, NNC_Polyhedron, Constraint_System
            sage: u = Variable(0)
            sage: v = Variable(1)
            sage: CS = Constraint_System()
            sage: CS.insert(0 <= u)
            sage: CS.insert(u <= 1)
            sage: CS.insert(0 <= v)
            sage: CS.insert(v <= 1)
            sage: CS.insert(u + v < 3)
            sage: P = NNC_Polyhedron(CS); P
            A 2-dimensional polyhedron in QQ^2 defined as the convex hull of 4 points
            sage: S = ManifoldSubsetPullback(c_cart, None, P); S
            Traceback (most recent call last):
            ...
            TypeError: unhashable type: 'NNC_Polyhedron'
            sage: S.is_closed()
            Traceback (most recent call last):
            ...
            NameError: name 'S' is not defined

        The pullback of real vector subspaces under a chart is closed::

            sage: V = span([[1, 2]], RR); V
            Vector space of degree 2 and dimension 1 over Real Field with 53 bits of precision
            Basis matrix:
            [1.00000000000000 2.00000000000000]
            sage: McV = ManifoldSubsetPullback(c_cart, None, V, name='McV'); McV
            Subset McV of the 2-dimensional topological manifold R^2
            sage: McV.is_closed()
            True

        The pullback of point lattices under a chart is closed::

            sage: W = span([[1, 0], [3, 5]], ZZ); W
            Free module of degree 2 and rank 2 over Integer Ring
            Echelon basis matrix:
            [1 0]
            [0 5]
            sage: McW = ManifoldSubsetPullback(c_cart, None, W, name='McW'); McW
            Subset McW of the 2-dimensional topological manifold R^2
            sage: McW.is_closed()
            True

        The pullback of finite sets is closed::

            sage: F = Family([vector(QQ, [1, 2], immutable=True), vector(QQ, [2, 3], immutable=True)])
            sage: McF = ManifoldSubsetPullback(c_cart, None, F, name='McF'); McF
            Subset McF of the 2-dimensional topological manifold R^2
            sage: McF.is_closed()
            False

        """
        if self.manifold().dimension() == 0:
            return True
        if isinstance(self._codomain_subset, ManifoldSubset):
            if self._codomain_subset.is_closed():
                # known closed
                return True
        elif isinstance(self._codomain_subset, RealSet):
            # RealSet can decide closedness authoritatively
            return self._codomain_subset.is_closed()
        elif is_Polyhedron(self._codomain_subset):
            # Regardless of their base_ring, we treat polyhedra as closed
            # convex subsets of R^n
            return True
        elif is_FreeModule(self._codomain_subset) and self._codomain_subset.rank() != infinity:
            if self._codomain_subset.base_ring() in MetricSpaces().Complete():
                # Closed topological vector subspace
                return True
            if self._codomain_subset.base_ring() == ZZ:
                if self._codomain_subset.coordinate_ring().is_subring(QQ):
                    # Discrete subgroup of R^n
                    return True
                if self._codomain_subset.rank() == self._codomain_subset.base_extend(RR).dimension():
                    # Discrete subgroup of R^n
                    return True
        else:
            if hasattr(self._codomain_subset, 'is_topologically_closed'):
                try:
                    from ppl import NNC_Polyhedron, C_Polyhedron
                except ImportError:
                    pass
                else:
                    if isinstance(self._codomain_subset, (NNC_Polyhedron, C_Polyhedron)):
                        # ppl polyhedra can decide closedness authoritatively
                        return self._codomain_subset.is_topologically_closed()
        return super().is_closed()

    def closure(self, name=None, latex_name=None):
        """
        Return the topological closure of ``self`` in the manifold.

        Because ``self`` is a pullback of some subset under a continuous map,
        the closure of ``self`` is the pullback of the closure.

        EXAMPLES::

            sage: from sage.manifolds.subsets.pullback import ManifoldSubsetPullback
            sage: M = Manifold(2, 'R^2', structure='topological')
            sage: c_cart.<x,y> = M.chart() # Cartesian coordinates on R^2
            sage: r_squared = M.scalar_field(x^2+y^2)
            sage: r_squared.set_immutable()
            sage: I = RealSet.open_closed(1, 2); I
            (1, 2]
            sage: O = ManifoldSubsetPullback(r_squared, None, I); O
            Subset f_inv_(1, 2] of the 2-dimensional topological manifold R^2
            sage: latex(O)
            f^{-1}((1, 2])
            sage: cl_O = O.closure(); cl_O
            Subset f_inv_[1, 2] of the 2-dimensional topological manifold R^2
            sage: cl_O.is_closed()
            True

        """
        if self.is_closed():
            return self
        try:
            codomain_subset_closure = self._codomain_subset.closure()
        except AttributeError:
            return super().closure()
        closure = ManifoldSubsetPullback(self._map, self._inverse,
                                         codomain_subset_closure,
                                         name=name, latex_name=latex_name)
        closure.declare_superset(self)
        return closure
