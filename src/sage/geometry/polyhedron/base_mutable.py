r"""
Base class for mutable polyhedra.

Just like vectors and matrices they can be set immutable.
The constructor does this by default.
"""

from sage.misc.abstract_method import abstract_method

from .base import Polyhedron_base
from .representation import VERTEX, RAY, LINE, INEQUALITY, EQUATION


class Polyhedron_mutable(Polyhedron_base):
    """
    Base class for polyhedra that allow mutability.

    This should not be used directly.
    """

    def __hash__(self):
        r"""
        TESTS::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: set([p])
            Traceback (most recent call last):
            ...
            TypeError: mutable polyhedra are unhashable
            sage: p.set_immutable()
            sage: set([p])
            {A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex}
        """
        if self._is_mutable:
            raise TypeError("mutable polyhedra are unhashable")
        return Polyhedron_base.__hash__(self)

    def _clear_cache(self):
        r"""
        Clear the Vrepresentation and Hrepresentation data of ``self``.

        TESTS::

            sage: p = polytopes.permutahedron(4)
            sage: P = p.parent()
            sage: q = P._element_constructor_(p, mutable=True)
            sage: TestSuite(q).run()
            sage: q._clear_cache()
            sage: TestSuite(q).run()

        ::

            sage: q.set_immutable()
            sage: q._clear_cache()
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra
        """
        if not self._is_mutable:
            raise TypeError("cannot clear cache of immutable polyhedra")

        # Invalidate object pointing towards this polyhedron (faces etc.).
        for ob in self._dependent_objects:
            ob._polyhedron = None
        backend_object = self.__dict__["_" + self._backend_object_name]
        del self.__dict__
        self.__dict__["_" + self._backend_object_name] = backend_object
        self._is_mutable = True
        self._dependent_objects = []

    def _add_dependent_object(self, ob):
        r"""
        Add an object that has ``self`` has attribute ``_polyhedron``.

        When ``self`` is modified, we delete this attribute to invalidate those objects.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: class foo:
            ....:     def __init__(self, p):
            ....:         self._polyhedron = p
            ....:
            sage: a = foo(p)
            sage: a.__dict__
            {'_polyhedron': A 0-dimensional polyhedron in ZZ^2 defined as the convex hull of 1 vertex}
            sage: p._add_dependent_object(a)
            sage: p._clear_cache()
            sage: a.__dict__
            {'_polyhedron': None}

        TESTS::

            sage: from sage.geometry.newton_polygon import NewtonPolygon
            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: n = NewtonPolygon(p)
            sage: n
            Finite Newton polygon with 1 vertex: (1, 1)
            sage: n = NewtonPolygon(p)
            sage: p._clear_cache()
            sage: n
            <repr(<sage.geometry.newton_polygon.ParentNewtonPolygon_with_category.element_class at ...>) failed: AttributeError: 'NoneType' object has no attribute 'vertices'>

        ::

            sage: f = p.faces(0)[0]; f
            A 0-dimensional face of a Polyhedron in ZZ^2 defined as the convex hull of 1 vertex
            sage: p._clear_cache()
            sage: f
            <repr(<sage.geometry.polyhedron.face.PolyhedronFace at ...>) failed: AttributeError: 'NoneType' object has no attribute 'parent'>

        ::

            sage: v = p.vertices()[0]
            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: v = p.Vrepresentation(0); v
            A vertex at (1, 1)
            sage: h = p.Hrepresentation(0); h
            An equation (0, 1) x - 1 == 0
            sage: p._clear_cache()
            sage: v.polyhedron() is None
            True
            sage: h.polyhedron() is None
            True

        ::

            sage: p = Polyhedron([[1, 0], [0, 1]], mutable=True)
            sage: r = p.relative_interior()
            sage: p._clear_cache()
            sage: r
            Relative interior of None
        """
        if ob._polyhedron is not self:
            raise ValueError
        self._dependent_objects.append(ob)

    def is_mutable(self):
        r"""
        Return True if the polyhedron is mutable, i.e. it can be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_mutable()
            True
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_mutable()
            False
        """
        return self._is_mutable

    def is_immutable(self):
        r"""
        Return True if the polyhedron is immutable, i.e. it cannot be modified in place.

        EXAMPLES::

            sage: p = Polyhedron([[1, 1]], mutable=True)
            sage: p.is_immutable()
            False
            sage: p = Polyhedron([[1, 1]], mutable=False)
            sage: p.is_immutable()
            True
        """
        return not self._is_mutable

    @abstract_method
    def set_immutable(self):
        r"""
        Make this polyhedron immutable. This operation cannot be undone.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.set_immutable(p)
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """

    @abstract_method
    def Vrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Vrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Vrepresentation(p)
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """
        # A derived class must implemented it to recalculate, if necessary.

    @abstract_method
    def Hrepresentation(self):
        r"""
        A derived class must overwrite such that it restores Hrepresentation
        after clearing it.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable.Hrepresentation(p)
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """
        # A derived class must implemented it to recalculate, if necessary.

    def _modification_method(method):
        r"""
        Proper wrapper around a modification method to make it safe and respect immutability.

        Modification methods should all be written to modify the polyhedron in place and
        this wrapper adds a keyword ``inplace`` such that

        - the new method returns a copy if ``inplace`` is ``False`` (default)

        - or the new method returns ``None`` and modifies the polyhedron in place
          if ``inplace`` is ``True``; but checks for mutability first.

        EXAMPLES::

            sage: def add_something(self, vertex, ray=None):
            ....:     '''
            ....:     Great documentation.
            ....:     '''
            ....:     self._modified = True
            ....:     print('added {}'.format(vertex))
            ....:     if ray is not None:
            ....:         print('added {}'.format(ray))
            ....:
            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: Polyhedron_mutable.add_something = Polyhedron_mutable._modification_method(add_something)
            sage: P = polytopes.cube()
            sage: Q = P.add_something('myVertex')
            added myVertex
            sage: Q._modified
            True
            sage: hasattr(P, '_modified')
            False
            sage: Q.is_mutable() == P.is_mutable()
            True
            sage: Q = P.add_something('myVertex', 'myRay')
            added myVertex
            added myRay
            sage: Q = P.add_something('myVertex', ray='myRay')
            added myVertex
            added myRay
            sage: P.add_something('myVertex', ray='myRay', inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra
            sage: parent = P.parent()
            sage: P_copy = parent._element_constructor_(P, mutable=True)
            sage: P_copy.add_something('myVertex', ray='myRay', inplace=True)
            added myVertex
            added myRay
            sage: P_copy._modified
            True
            sage: P_copy.add_something('myVertex', inplace=True)
            added myVertex
            sage: P_copy.add_something('myVertex', inplace=True, ray='myRay')
            added myVertex
            added myRay
            sage: P_copy.add_something('myVertex', 'myRay', inplace=True)
            added myVertex
            added myRay
            sage: from inspect import signature
            sage: signature(P.add_something)
            <Signature (vertex, ray=None, *, inplace=False)>
            sage: P.add_something.__doc__
            '\n    Great documentation.\n    '
            sage: from sage.misc.sageinspect import sage_getsource
            sage: sage_getsource(P.add_something)
            "def add_something(self, vertex, ray=None):\n    '''\n    Great documentation.\n    '''\n    self._modified = True\n    print('added {}'.format(vertex))\n    if ray is not None:\n        print('added {}'.format(ray))\n"
            sage: del Polyhedron_mutable.add_something
        """
        def modification_method(self, *args, inplace=False, **kwds):
            r"""
            Modify the polyhedron according to the method name.

            INPUT:

            - ``args`` -- arguments according to the function name

            - ``inplace`` -- boolean (default: ``False``);
              if ``True`` modify ``self`` and return ``None``;
              if ``False`` return a modified copy
            """
            if not inplace:
                mutable = self.is_mutable()
                P = self.parent()
                self = P._element_constructor_(self, mutable=True)
            self._clear_cache()
            method.__call__(self, *args, **kwds)
            if not inplace:
                if not mutable:
                    self.set_immutable()
                return self

        modification_method.__doc__ = method.__doc__

        # Merge the signatures.
        from inspect import signature
        sig1 = signature(method)
        sig2 = signature(modification_method)
        new_sig = sig1.replace(parameters=tuple(sig1.parameters.values()) + (sig2.parameters['inplace'],))
        modification_method.__signature__ = new_sig

        from sage.misc.sageinspect import sage_getsource
        def get_source():
            return sage_getsource(method)

        modification_method._is_modification_method = True

        modification_method._sage_src_ = get_source

        return modification_method

    def _test_modification_methods(self, tester=None, **option):
        r"""
        Check that all modification methods document the usage of ``inplace``.

        EXAMPLES::

            sage: P = polytopes.cube()
            sage: TestSuite(P).run()
            sage: def add_something(self, vertex, ray=None):
            ....:     self._modified = True
            ....:     print('added {}'.format(vertex))
            ....:     if ray is not None:
            ....:         print('added {}'.format(ray))
            ....:
            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: Polyhedron_mutable.add_something = Polyhedron_mutable._modification_method(add_something)
            sage: TestSuite(P).run()
            Failure in _test_modification_methods:
            ...
            AssertionError: None is not an instance of <class 'str'>
            ------------------------------------------------------------
            The following tests failed: _test_modification_methods
            sage: def add_something(self, vertex, ray=None):
            ....:     '''
            ....:     Incomplete documentation.
            ....:     '''
            ....:     self._modified = True
            ....:     print('added {}'.format(vertex))
            ....:     if ray is not None:
            ....:         print('added {}'.format(ray))
            sage: Polyhedron_mutable.add_something = Polyhedron_mutable._modification_method(add_something)
            sage: TestSuite(P).run()
            Failure in _test_modification_methods:
            ...
            AssertionError: '- ``inplace`` -- boolean (default: ``False``)' not found in '\n    Incomplete documentation.\n    '
            ------------------------------------------------------------
            The following tests failed: _test_modification_methods
            sage: def add_something(self, vertex, ray=None):
            ....:     '''
            ....:     Still incomplete documentation.
            ....:
            ....:     INPUT:
            ....:
            ....:     - ``vertex`` -- vertex
            ....:
            ....:     - ``ray`` -- ray
            ....:
            ....:     - ``inplace`` -- boolean (default: ``False``)
            ....:     '''
            ....:     self._modified = True
            ....:     print('added {}'.format(vertex))
            ....:     if ray is not None:
            ....:         print('added {}'.format(ray))
            sage: Polyhedron_mutable.add_something = Polyhedron_mutable._modification_method(add_something)
            sage: TestSuite(P).run()
            Failure in _test_modification_methods:
            ...
            AssertionError: 'inplace=True' not found in '\n    Still incomplete documentation.\n\n    INPUT:\n\n    - ``vertex`` -- vertex\n\n    - ``ray`` -- ray\n\n    - ``inplace`` -- boolean (default: ``False``)\n    '
            ------------------------------------------------------------
            The following tests failed: _test_modification_methods
            sage: def add_something(self, vertex, ray=None):
            ....:     '''
            ....:     Somewhat complete documentation.
            ....:
            ....:     INPUT:
            ....:
            ....:     - ``vertex`` -- vertex
            ....:
            ....:     - ``ray`` -- ray
            ....:
            ....:     - ``inplace`` -- boolean (default: ``False``)
            ....:
            ....:     EXAMPLES:
            ....:
            ....:         sage: P = polytopes.cube()
            ....:         sage: Q = P.add_something('myVertex', ray='myRay')
            ....:         added myVertex
            ....:         added myRay
            ....:         sage: P.add_something('myVertex', ray='myRay', inplace=True)
            ....:         Traceback (most recent call last):
            ....:         ...
            ....:         TypeError: cannot clear cache of immutable polyhedra
            ....:     '''
            ....:     self._modified = True
            ....:     print('added {}'.format(vertex))
            ....:     if ray is not None:
            ....:         print('added {}'.format(ray))
            sage: del Polyhedron_mutable.add_something
        """
        if tester is None:
            tester = self._tester(**options)

        for name in dir(self):
            try:
                attr = getattr(self, name)
            except (NotImplementedError, AttributeError):
                pass
            else:
                if hasattr(attr, '_is_modification_method'):
                    tester.assertIsInstance(attr.__doc__, str)
                    tester.assertIn('- ``inplace`` -- boolean (default: ``False``)', attr.__doc__)
                    tester.assertIn('inplace=True', attr.__doc__)

    @_modification_method
    def add_vertex(self, vertex):
        r"""
        Add a vertex to the polyhedron.

        INPUT:

        - ``vertex`` -- a point

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: P = polytopes.cube()
            sage: Q = P.add_vertex([4, 5, 6])
            sage: Q.vertices()
            (A vertex at (-1, -1, -1),
             A vertex at (-1, -1, 1),
             A vertex at (-1, 1, -1),
             A vertex at (-1, 1, 1),
             A vertex at (1, -1, -1),
             A vertex at (1, -1, 1),
             A vertex at (1, 1, -1),
             A vertex at (4, 5, 6))

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_vertex([4, 5, 6], inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.vertices()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A vertex at (1, 1, 1))
            sage: P_copy.add_vertex([4, 5, 6], inplace=True)
            sage: P_copy.vertices()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A vertex at (4, 5, 6))
            sage: P_copy.add_vertex([-4, -5, -6], inplace=True)
            sage: P_copy.vertices()
            (A vertex at (-4, -5, -6),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A vertex at (4, 5, 6))
        """
        self._add_generator(vertex, VERTEX)

    @_modification_method
    def add_ray(self, ray):
        r"""
        Add a ray to the polyhedron.

        INPUT:

        - ``ray`` -- a ray

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: P = polytopes.cube()
            sage: Q = P.add_ray([4, 5, 6])
            sage: Q.Vrepresentation()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A ray in the direction (4, 5, 6))

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_ray([4, 5, 6], inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.Vrepresentation()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A vertex at (1, 1, 1))
            sage: P_copy.add_ray([4/2, 5/2, 6/2], inplace=True)
            sage: P_copy.Vrepresentation()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A ray in the direction (4, 5, 6))
            sage: P_copy.add_ray([-4/2, 5/2, 6/2], inplace=True)
            sage: P_copy.Vrepresentation()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A ray in the direction (-4, 5, 6),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A ray in the direction (4, 5, 6))
        """
        self._add_generator(ray, RAY)

    @_modification_method
    def add_line(self, line):
        r"""
        Add a line to the polyhedron.

        INPUT:

        - ``line`` -- a line

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: P = polytopes.cube()
            sage: Q = P.add_line([1, 1, 1])
            sage: Q.Vrepresentation()
            (A line in the direction (1, 1, 1),
            A vertex at (2, 2, 0),
            A vertex at (-2, -2, 0),
            A vertex at (0, 2, 0),
            A vertex at (-2, 0, 0),
            A vertex at (2, 0, 0),
            A vertex at (0, -2, 0))

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_line([1, 1, 1], inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.Vrepresentation()
            (A vertex at (-1, -1, -1),
            A vertex at (-1, -1, 1),
            A vertex at (-1, 1, -1),
            A vertex at (-1, 1, 1),
            A vertex at (1, -1, -1),
            A vertex at (1, -1, 1),
            A vertex at (1, 1, -1),
            A vertex at (1, 1, 1))
            sage: P_copy.add_line([1, 1, 1], inplace=True)
            sage: P_copy.Vrepresentation()
            (A line in the direction (1, 1, 1),
             A vertex at (-2, -2, 0),
             A vertex at (0, 2, 0),
             A vertex at (-2, 0, 0),
             A vertex at (2, 0, 0),
             A vertex at (0, -2, 0),
             A vertex at (2, 2, 0))
        """
        self._add_generator(line, LINE)

    @_modification_method
    def add_Vrepresentatives(self, vertices, rays, lines):
        r"""
        Add Vrepresentatives to the polyhedron.

        INPUT:

        - ``vertices`` -- iterable of vertices or None

        - ``rays`` -- iterable of rays or None

        - ``lines`` -- iterable of lines or None

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: from itertools import product
            sage: P = polytopes.cross_polytope(3)
            sage: Q = P.add_Vrepresentatives([p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None, None)
            sage: Q.Vrepresentation()
            (A vertex at (-1, -1, 0),
            A vertex at (-1, 0, -1),
            A vertex at (-1, 0, 0),
            A vertex at (0, -1, -1),
            A vertex at (0, -1, 0),
            A vertex at (0, 0, -1),
            A vertex at (0, 0, 1),
            A vertex at (0, 1, 0),
            A vertex at (0, 1, 1),
            A vertex at (1, 0, 0),
            A vertex at (1, 0, 1),
            A vertex at (1, 1, 0))

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_Vrepresentatives([p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None, None, inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.f_vector()
            (1, 6, 12, 8, 1)
            sage: P_copy.add_Vrepresentatives([p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None, None, inplace=True)
            sage: P_copy.f_vector()
            (1, 12, 24, 14, 1)
        """
        self._add_generators(vertices, rays, lines)

    @_modification_method
    def add_inequality(self, inequality):
        r"""
        Add an inequality to the polyhedron.

        INPUT:

        - ``inequality`` -- an inequality

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: P = polytopes.cube().base_extend(QQ)
            sage: Q = P.add_inequality([1, -1, 1, 0])
            sage: Q.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (0, 1, 0) x + 1 >= 0,
             An inequality (1, 0, 0) x + 1 >= 0,
             An inequality (-1, 1, 0) x + 1 >= 0)

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_inequality([1, -1, 1, 0], inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.Hrepresentation()
            (An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0)
            sage: P_copy.add_inequality([1, -1, 1, 0], inplace=True)
            sage: P_copy.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (-1, 1, 0) x + 1 >= 0)
            sage: P_copy.add_inequality([1, -1, 1, 1], inplace=True)
            sage: P_copy.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (-1, 1, 0) x + 1 >= 0,
             An inequality (0, -1, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (0, 1, 0) x + 1 >= 0,
             An inequality (1, 0, 0) x + 1 >= 0,
             An inequality (-1, 1, 1) x + 1 >= 0)
        """
        self._add_constraint(inequality, INEQUALITY)

    @_modification_method
    def add_equation(self, equation):
        r"""
        Add an equation to the polyhedron.

        INPUT:

        - ``equation`` -- an equation

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: P = polytopes.cube().base_extend(QQ)
            sage: Q = P.add_equation([1, -1, 1, 0])
            sage: Q.Hrepresentation()
            (An equation (1, -1, 0) x - 1 == 0,
             An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (1, 0, 0) x + 0 >= 0)

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_equation([1, -1, 1, 0], inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.Hrepresentation()
            (An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0)
            sage: P_copy.add_equation([1, -1, 1, 0], inplace=True)
            sage: P_copy.Hrepresentation()
            (An equation (1, -1, 0) x - 1 == 0,
             An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (0, 0, -1) x + 1 >= 0,
             An inequality (0, 0, 1) x + 1 >= 0,
             An inequality (1, 0, 0) x + 0 >= 0)
            sage: P_copy.add_equation([1, -1, 1, 1], inplace=True)
            sage: P_copy.Hrepresentation()
            (An equation (0, 0, 1) x + 0 == 0,
             An equation (1, -1, 0) x - 1 == 0,
             An inequality (-1, 0, 0) x + 1 >= 0,
             An inequality (1, 0, 0) x + 0 >= 0)
        """
        self._add_constraint(equation, EQUATION)

    @_modification_method
    def add_Hrepresentatives(self, ieqs, eqns):
        r"""
        Add Hrepresentatives to the polyhedron.

        INPUT:

        - ``ieqs`` -- iterable of inequalities or None

        - ``eqns`` -- iterable of equations or None

        - ``inplace`` -- boolean (default: ``False``);
          if ``True`` modify ``self`` and return ``None``;
          if ``False`` return a modified copy

        EXAMPLES:

        Construct a modified copy::

            sage: from itertools import product
            sage: P = polytopes.cube().base_extend(QQ)
            sage: Q = P.add_Hrepresentatives([(1,) + p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None)
            sage: Q.Hrepresentation()
            (An inequality (-1, 0, 0) x + 1 >= 0,
            An inequality (0, -1, 0) x + 1 >= 0,
            An inequality (0, 0, -1) x + 1 >= 0,
            An inequality (0, 0, 1) x + 1 >= 0,
            An inequality (0, 1, 0) x + 1 >= 0,
            An inequality (1, 0, 0) x + 1 >= 0,
            An inequality (-1, -1, 0) x + 1 >= 0,
            An inequality (-1, 0, -1) x + 1 >= 0,
            An inequality (0, -1, -1) x + 1 >= 0,
            An inequality (0, 1, 1) x + 1 >= 0,
            An inequality (1, 0, 1) x + 1 >= 0,
            An inequality (1, 1, 0) x + 1 >= 0)

        To modify the polyhedron in place, it must be mutable:

            sage: P.add_Hrepresentatives([(1,) + p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None, inplace=True)
            Traceback (most recent call last):
            ...
            TypeError: cannot clear cache of immutable polyhedra

        Make a mutable copy to modify the polyhedron in place:

            sage: P_copy = P.parent()._element_constructor_(P, mutable=True)
            sage: P_copy.f_vector()
            (1, 8, 12, 6, 1)
            sage: P_copy.add_Hrepresentatives([(1,) + p for p in product([1, -1, 0], repeat=3) if sum(p) in (2, -2)], None, inplace=True)
            sage: P_copy.f_vector()
            (1, 14, 24, 12, 1)
        """
        self._add_constraints(ieqs, eqns)

    @abstract_method
    def _add_generator(self, generator, typ):
        r"""
        Add a generator to the polyhedron in place.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable._add_generator()
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """

    @abstract_method
    def _add_generators(self, vertices, rays, lines):
        r"""
        Add generators to the polyhedron in place.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable._add_generators()
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """

    @abstract_method
    def _add_constraint(self, constraint, typ):
        r"""
        Add a constraint to the polyhedron in place.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable._add_constraint()
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """

    @abstract_method
    def _add_constraints(self, ieqs, eqns):
        r"""
        Add constraints to the polyhedron in place.

        TESTS::

            sage: from sage.geometry.polyhedron.base_mutable import Polyhedron_mutable
            sage: p = polytopes.cube()
            sage: Polyhedron_mutable._add_constraints()
            Traceback (most recent call last):
            ...
            TypeError: 'AbstractMethod' object is not callable
        """
