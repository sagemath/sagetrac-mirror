# -*- coding: utf-8 -*-
r"""
Generic cell complexes

AUTHORS:

- John H. Palmieri (2009-08)

This module defines a class of abstract finite cell complexes.  This
is meant as a base class from which other classes (like
:class:`~sage.homology.simplicial_complex.SimplicialComplex`,
:class:`~sage.homology.cubical_complex.CubicalComplex`, and
:class:`~sage.homology.delta_complex.DeltaComplex`) should derive.  As
such, most of its properties are not implemented.  It is meant for use
by developers producing new classes, not casual users.

.. NOTE::

    Keywords for :meth:`~GenericCellComplex.chain_complex`,
    :meth:`~GenericCellComplex.homology`, etc.: any keywords given to
    the :meth:`~GenericCellComplex.homology` method get passed on to
    the :meth:`~GenericCellComplex.chain_complex` method and also to
    the constructor for chain complexes in
    :class:`sage.homology.chain_complex.ChainComplex_class <ChainComplex>`,
    as well as its associated
    :meth:`~sage.homology.chain_complex.ChainComplex_class.homology` method.
    This means that those keywords should have consistent meaning in
    all of those situations.  It also means that it is easy to
    implement new keywords: for example, if you implement a new
    keyword for the
    :meth:`sage.homology.chain_complex.ChainComplex_class.homology` method,
    then it will be automatically accessible through the
    :meth:`~GenericCellComplex.homology` method for cell complexes --
    just make sure it gets documented.
"""


from sage.structure.sage_object import SageObject
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.combinat.free_module import CombinatorialFreeModule, CombinatorialFreeModuleElement
from sage.misc.abstract_method import abstract_method

class GenericCellComplex(SageObject):
    r"""
    Class of abstract cell complexes.

    This is meant to be used by developers to produce new classes, not
    by casual users.  Classes which derive from this are
    :class:`~sage.homology.simplicial_complex.SimplicialComplex`,
    :class:`~sage.homology.delta_complex.DeltaComplex`, and
    :class:`~sage.homology.cubical_complex.CubicalComplex`.

    Most of the methods here are not implemented, but probably should
    be implemented in a derived class.  Most of the other methods call
    a non-implemented one; their docstrings contain examples from
    derived classes in which the various methods have been defined.
    For example, :meth:`homology` calls :meth:`chain_complex`; the
    class :class:`~sage.homology.delta_complex.DeltaComplex`
    implements
    :meth:`~sage.homology.delta_complex.DeltaComplex.chain_complex`,
    and so the :meth:`homology` method here is illustrated with
    examples involving `\Delta`-complexes.

    EXAMPLES:

    It's hard to give informative examples of the base class, since
    essentially nothing is implemented. ::

        sage: from sage.homology.cell_complex import GenericCellComplex
        sage: A = GenericCellComplex()
    """
    def __eq__(self,right):
        """
        Comparisons of cell complexes are not implemented.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A == B # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def __ne__(self,right):
        """
        Comparisons of cell complexes are not implemented.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A != B # indirect doctest
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    ############################################################
    # self.cells() and related methods
    ############################################################

    @abstract_method
    def cells(self, subcomplex=None):
        """
        The cells of this cell complex, in the form of a dictionary:
        the keys are integers, representing dimension, and the value
        associated to an integer `d` is the set of `d`-cells.  If the
        optional argument ``subcomplex`` is present, then return only
        the faces which are *not* in the subcomplex.

        :param subcomplex: a subcomplex of this cell complex.  Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default None

        This is not implemented in general; it should be implemented
        in any derived class.  When implementing, see the warning in
        the :meth:`dimension` method.

        This method is used by various other methods, such as
        :meth:`n_cells` and :meth:`f_vector`.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.cells()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method cells at ...>
        """

    def dimension(self):
        """
        The dimension of this cell complex: the maximum
        dimension of its cells.

        .. WARNING::

          If the :meth:`cells` method calls :meth:`dimension`,
          then you'll get an infinite loop.  So either don't use
          :meth:`dimension` or override :meth:`dimension`.

        EXAMPLES::

            sage: simplicial_complexes.RandomComplex(d=5, n=8).dimension()
            5
            sage: delta_complexes.Sphere(3).dimension()
            3
            sage: T = cubical_complexes.Torus()
            sage: T.product(T).dimension()
            4
        """
        try:
            return max([x.dimension() for x in self._facets])
        except AttributeError:
            return max(self.cells())

    def n_cells(self, n, subcomplex=None):
        """
        List of cells of dimension ``n`` of this cell complex.
        If the optional argument ``subcomplex`` is present, then
        return the ``n``-dimensional faces which are *not* in the
        subcomplex.

        :param n: the dimension
        :type n: non-negative integer
        :param subcomplex: a subcomplex of this cell complex. Return
           the cells which are not in this subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES::

            sage: delta_complexes.Torus().n_cells(1)
            [(0, 0), (0, 0), (0, 0)]
            sage: cubical_complexes.Cube(1).n_cells(0)
            [[1,1], [0,0]]
        """
        if n in self.cells(subcomplex):
            return list(self.cells(subcomplex)[n])
        else:
            # don't barf if someone asks for n_cells in a dimension where there are none
            return []

    def f_vector(self):
        """
        The `f`-vector of this cell complex: a list whose `n^{th}`
        item is the number of `(n-1)`-cells.  Note that, like all
        lists in Sage, this is indexed starting at 0: the 0th element
        in this list is the number of `(-1)`-cells (which is 1: the
        empty cell is the only `(-1)`-cell).

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle().f_vector()
            [1, 8, 24, 16]
            sage: delta_complexes.KleinBottle().f_vector()
            [1, 1, 3, 2]
            sage: cubical_complexes.KleinBottle().f_vector()
            [1, 42, 84, 42]
        """
        return [self._f_dict()[n] for n in range(-1, self.dimension()+1)]

    def _f_dict(self):
        """
        The `f`-vector of this cell complex as a dictionary: the
        item associated to an integer `n` is the number of the
        `n`-cells.

        EXAMPLES::

            sage: simplicial_complexes.KleinBottle()._f_dict()[1]
            24
            sage: delta_complexes.KleinBottle()._f_dict()[1]
            3
        """
        answer = {}
        answer[-1] = 1
        for n in range(self.dimension() + 1):
            answer[n] = len(self.cells()[n])
        return answer

    def euler_characteristic(self):
        r"""
        The Euler characteristic of this cell complex: the
        alternating sum over `n \geq 0` of the number of
        `n`-cells.

        EXAMPLES::

            sage: simplicial_complexes.Simplex(5).euler_characteristic()
            1
            sage: delta_complexes.Sphere(6).euler_characteristic()
            2
            sage: cubical_complexes.KleinBottle().euler_characteristic()
            0
        """
        return sum([(-1)**n * self.f_vector()[n+1] for n in range(self.dimension() + 1)])

    ############################################################
    # end of methods using self.cells()
    ############################################################

    @abstract_method
    def product(self, right, rename_vertices=True):
        """
        The (Cartesian) product of this cell complex with another one.

        Products are not implemented for general cell complexes.  They
        may be implemented in some derived classes (like simplicial
        complexes).

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.product(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method product at ...>
        """

    @abstract_method
    def disjoint_union(self, right):
        """
        The disjoint union of this cell complex with another one.

        :param right: the other cell complex (the right-hand factor)

        Disjoint unions are not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.disjoint_union(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method disjoint_union at ...>
        """

    @abstract_method
    def wedge(self, right):
        """
        The wedge (one-point union) of this cell complex with
        another one.

        :param right: the other cell complex (the right-hand factor)

        Wedges are not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.wedge(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method wedge at ...>
        """

    ############################################################
    # self.join() and related methods
    ############################################################

    @abstract_method
    def join(self, right, **kwds):
        """
        The join of this cell complex with another one.

        :param right: the other cell complex (the right-hand factor)

        Joins are not implemented for general cell complexes.  They
        may be implemented in some derived classes (like simplicial
        complexes).

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex(); B = GenericCellComplex()
            sage: A.join(B)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method join at ...>
        """

    # for some classes, you may want * to mean join:
    ###
    # __mul__ = join

    # the cone on X is the join of X with a point.  See
    # simplicial_complex.py for one implementation.
    ###
    # def cone(self):
    #     return self.join(POINT)

    # the suspension of X is the join of X with the 0-sphere (two
    # points).  See simplicial_complex.py for one implementation.
    ###
    # def suspension(self, n=1):
    #     """
    #     The suspension of this cell complex.
    #
    #     INPUT:
    #
    #     -  ``n`` - positive integer (optional, default 1): suspend this
    #        many times.
    #     """
    #     raise NotImplementedError

    ############################################################
    # end of methods using self.join()
    ############################################################

    ############################################################
    # chain complexes, homology
    ############################################################

    @abstract_method
    def chain_complex(self, **kwds):
        """
        This is not implemented for general cell complexes.

        Some keywords to possibly implement in a derived class:

        - ``subcomplex`` -- a subcomplex: compute the relative chain complex
        - ``augmented`` -- a bool: whether to return the augmented complex
        - ``verbose`` -- a bool: whether to print informational messages as
          the chain complex is being computed
        - ``check_diffs`` -- a bool: whether to check that the each
          composite of two consecutive differentials is zero
        -  ``dimensions`` -- if ``None``, compute the chain complex in all
           dimensions.  If a list or tuple of integers, compute the
           chain complex in those dimensions, setting the chain groups
           in all other dimensions to zero.

        Definitely implement the following:

        -  ``base_ring`` -- commutative ring (optional, default ZZ)
        -  ``cochain`` -- a bool: whether to return the cochain complex

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.chain_complex()
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method chain_complex at ...>
        """

    def homology(self, dim=None, **kwds):
        r"""
        The (reduced) homology of this cell complex.

        :param dim: If None, then return the homology in every
           dimension.  If ``dim`` is an integer or list, return the
           homology in the given dimensions.  (Actually, if ``dim`` is
           a list, return the homology in the range from ``min(dim)``
           to ``max(dim)``.)
        :type dim: integer or list of integers or None; optional,
           default None
        :param base_ring: commutative ring, must be ZZ or a field.
        :type base_ring: optional, default ZZ
        :param subcomplex: a subcomplex of this simplicial complex.
           Compute homology relative to this subcomplex.
        :type subcomplex: optional, default empty
        :param generators: If ``True``, return generators for the homology
           groups along with the groups.  NOTE: Since :trac:`6100`, the result
           may not be what you expect when not using CHomP since its return
           is in terms of the chain complex.
        :type generators: boolean; optional, default False
        :param cohomology: If True, compute cohomology rather than homology.
        :type cohomology: boolean; optional, default False
        :param algorithm: The options are 'auto', 'dhsw', 'pari' or 'no_chomp'.
           See below for a description of what they mean.
        :type algorithm: string; optional, default 'auto'
        :param verbose: If True, print some messages as the homology is
           computed.
        :type verbose: boolean; optional, default False
        :param reduced: If ``True``, return the reduced homology.
        :type reduced: boolean; optional, default ``True``

        .. note::

            The keyword arguments to this function get passed on to
            :meth:`chain_complex` and its homology.

        ALGORITHM:

        If ``algorithm`` is set to 'auto' (the default), then use
        CHomP if available.  (CHomP is available at the web page
        http://chomp.rutgers.edu/.  It is also an experimental package
        for Sage.)

        CHomP computes homology, not cohomology, and only works over
        the integers or finite prime fields.  Therefore if any of
        these conditions fails, or if CHomP is not present, or if
        ``algorithm`` is set to 'no_chomp', go to plan B: if this complex
        has a ``_homology`` method -- each simplicial complex has
        this, for example -- then call that.  Such a method implements
        specialized algorithms for the particular type of cell
        complex.

        Otherwise, move on to plan C: compute the chain complex of
        this complex and compute its homology groups.  To do this: over a
        field, just compute ranks and nullities, thus obtaining
        dimensions of the homology groups as vector spaces.  Over the
        integers, compute Smith normal form of the boundary matrices
        defining the chain complex according to the value of
        ``algorithm``.  If ``algorithm`` is 'auto' or 'no_chomp', then
        for each relatively small matrix, use the standard Sage
        method, which calls the Pari package.  For any large matrix,
        reduce it using the Dumas, Heckenbach, Saunders, and Welker
        elimination algorithm: see
        :func:`sage.homology.matrix_utils.dhsw_snf` for details.

        Finally, ``algorithm`` may also be 'pari' or 'dhsw', which
        forces the named algorithm to be used regardless of the size
        of the matrices and regardless of whether CHomP is available.

        As of this writing, CHomP is by far the fastest option,
        followed by the 'auto' or 'no_chomp' setting of using the
        Dumas, Heckenbach, Saunders, and Welker elimination algorithm
        for large matrices and Pari for small ones.

        EXAMPLES::

            sage: P = delta_complexes.RealProjectivePlane()
            sage: P.homology()
            {0: 0, 1: C2, 2: 0}
            sage: P.homology(reduced=False)
            {0: Z, 1: C2, 2: 0}
            sage: P.homology(base_ring=GF(2))
            {0: Vector space of dimension 0 over Finite Field of size 2,
             1: Vector space of dimension 1 over Finite Field of size 2,
             2: Vector space of dimension 1 over Finite Field of size 2}
            sage: S7 = delta_complexes.Sphere(7)
            sage: S7.homology(7)
            Z
            sage: cubical_complexes.KleinBottle().homology(1, base_ring=GF(2))
            Vector space of dimension 2 over Finite Field of size 2

        If CHomP is installed, Sage can compute generators of homology
        groups::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S2.homology(dim=2, generators=True, base_ring=GF(2))  # optional - CHomP
            (Vector space of dimension 1 over Finite Field of size 2, [(0, 1, 2) + (0, 1, 3) + (0, 2, 3) + (1, 2, 3)])

        When generators are computed, Sage returns a pair for each
        dimension: the group and the list of generators.  For
        simplicial complexes, each generator is represented as a
        linear combination of simplices, as above, and for cubical
        complexes, each generator is a linear combination of cubes::

            sage: S2_cub = cubical_complexes.Sphere(2)
            sage: S2_cub.homology(dim=2, generators=True)  # optional - CHomP
            (Z, [-[[0,1] x [0,1] x [0,0]] + [[0,1] x [0,1] x [1,1]] - [[0,0] x [0,1] x [0,1]] - [[0,1] x [1,1] x [0,1]] + [[0,1] x [0,0] x [0,1]] + [[1,1] x [0,1] x [0,1]]])
        """
        from sage.interfaces.chomp import have_chomp, homcubes, homsimpl
        from sage.homology.cubical_complex import CubicalComplex
        from sage.homology.simplicial_complex import SimplicialComplex
        from sage.modules.all import VectorSpace
        from sage.homology.homology_group import HomologyGroup

        base_ring = kwds.get('base_ring', ZZ)
        cohomology = kwds.get('cohomology', False)
        subcomplex = kwds.pop('subcomplex', None)
        verbose = kwds.get('verbose', False)
        algorithm = kwds.get('algorithm', 'auto')
        reduced = kwds.get('reduced', True)

        if dim is not None:
            if isinstance(dim, (list, tuple)):
                low = min(dim) - 1
                high = max(dim) + 2
            else:
                low = dim - 1
                high = dim + 2
            dims = range(low, high)
        else:
            dims = None

        # try to use CHomP if computing homology (not cohomology) and
        # working over Z or F_p, p a prime.
        if (algorithm == 'auto' and cohomology is False
            and (base_ring == ZZ or (base_ring.is_prime_field()
                                     and base_ring != QQ))):
            # homcubes, homsimpl seems fastest if all of homology is computed.
            H = None
            if isinstance(self, CubicalComplex):
                if have_chomp('homcubes'):
                    H = homcubes(self, subcomplex, **kwds)
            elif isinstance(self, SimplicialComplex):
                if have_chomp('homsimpl'):
                    H = homsimpl(self, subcomplex, **kwds)

            # now pick off the requested dimensions
            if H:
                answer = {}
                if not dims:
                    dims =range(self.dimension() + 1)
                for d in dims:
                    answer[d] = H.get(d, HomologyGroup(0, base_ring))
                if dim is not None:
                    if not isinstance(dim, (list, tuple)):
                        answer = answer.get(dim, HomologyGroup(0, base_ring))
                return answer

        # Derived classes can implement specialized algorithms using a
        # _homology_ method.  See SimplicialComplex for one example.
        if hasattr(self, '_homology_'):
            return self._homology_(dim, subcomplex=subcomplex, **kwds)

        C = self.chain_complex(cochain=cohomology, augmented=reduced,
                               dimensions=dims, subcomplex=subcomplex, **kwds)
        answer = C.homology(**kwds)
        if dim is None:
            dim = range(self.dimension()+1)
        zero = HomologyGroup(0, base_ring)
        if isinstance(dim, (list, tuple)):
            return dict([d, answer.get(d, zero)] for d in dim)
        return answer.get(dim, zero)

    def cohomology(self, dim=None, **kwds):
        r"""
        The reduced cohomology of this cell complex.

        The arguments are the same as for the :meth:`homology` method,
        except that :meth:`homology` accepts a ``cohomology`` key
        word, while this function does not: ``cohomology`` is
        automatically true here.  Indeed, this function just calls
        :meth:`homology` with ``cohomology`` set to ``True``.

        :param dim:
        :param base_ring:
        :param subcomplex:
        :param algorithm:
        :param verbose:
        :param reduced:

        EXAMPLES::

            sage: circle = SimplicialComplex([[0,1], [1,2], [0, 2]])
            sage: circle.cohomology(0)
            0
            sage: circle.cohomology(1)
            Z
            sage: P2 = SimplicialComplex([[0,1,2], [0,2,3], [0,1,5], [0,4,5], [0,3,4], [1,2,4], [1,3,4], [1,3,5], [2,3,5], [2,4,5]])   # projective plane
            sage: P2.cohomology(2)
            C2
            sage: P2.cohomology(2, base_ring=GF(2))
            Vector space of dimension 1 over Finite Field of size 2
            sage: P2.cohomology(2, base_ring=GF(3))
            Vector space of dimension 0 over Finite Field of size 3

            sage: cubical_complexes.KleinBottle().cohomology(2)
            C2

        Relative cohomology::

            sage: T = SimplicialComplex([[0,1]])
            sage: U = SimplicialComplex([[0], [1]])
            sage: T.cohomology(1, subcomplex=U)
            Z

        A `\Delta`-complex example::

            sage: s5 = delta_complexes.Sphere(5)
            sage: s5.cohomology(base_ring=GF(7))[5]
            Vector space of dimension 1 over Finite Field of size 7
        """
        return self.homology(dim=dim, cohomology=True, **kwds)

    def betti(self, dim=None, subcomplex=None):
        r"""
        The Betti numbers of this simplicial complex as a dictionary
        (or a single Betti number, if only one dimension is given):
        the ith Betti number is the rank of the ith homology group.

        :param dim: If ``None``, then return every Betti number, as
           a dictionary with keys the non-negative integers.  If
           ``dim`` is an integer or list, return the Betti number for
           each given dimension.  (Actually, if ``dim`` is a list,
           return the Betti numbers, as a dictionary, in the range
           from ``min(dim)`` to ``max(dim)``.  If ``dim`` is a number,
           return the Betti number in that dimension.)
        :type dim: integer or list of integers or ``None``; optional,
           default ``None``
        :param subcomplex: a subcomplex of this cell complex.  Compute
           the Betti numbers of the homology relative to this subcomplex.
        :type subcomplex: optional, default ``None``

        EXAMPLES:

        Build the two-sphere as a three-fold join of a
        two-point space with itself::

            sage: S = SimplicialComplex([[0], [1]])
            sage: (S*S*S).betti()
            {0: 1, 1: 0, 2: 1}
            sage: (S*S*S).betti([1,2])
            {1: 0, 2: 1}
            sage: (S*S*S).betti(2)
            1

        Or build the two-sphere as a `\Delta`-complex::

            sage: S2 = delta_complexes.Sphere(2)
            sage: S2.betti([1,2])
            {1: 0, 2: 1}

        Or as a cubical complex::

            sage: S2c = cubical_complexes.Sphere(2)
            sage: S2c.betti(2)
            1
        """
        dict = {}
        H = self.homology(dim, base_ring=QQ, subcomplex=subcomplex)
        try:
            for n in H.keys():
                dict[n] = H[n].dimension()
                if n == 0:
                    dict[n] += 1
            return dict
        except AttributeError:
            return H.dimension()

    def is_acyclic(self, base_ring=ZZ):
        """
        True if the reduced homology with coefficients in ``base_ring`` of
        this cell complex is zero.

        INPUT:

        - ``base_ring`` -- optional, default ``ZZ``. Compute homology
          with coefficients in this ring.

        EXAMPLES::

            sage: RP2 = simplicial_complexes.RealProjectivePlane()
            sage: RP2.is_acyclic()
            False
            sage: RP2.is_acyclic(QQ)
            True

        This first computes the Euler characteristic: if it is not 1,
        the complex cannot be acyclic. So this should return ``False``
        reasonably quickly on complexes with Euler characteristic not
        equal to 1::

            sage: K = cubical_complexes.KleinBottle()
            sage: C = cubical_complexes.Cube(2)
            sage: P = K.product(C)
            sage: P
            Cubical complex with 168 vertices and 1512 cubes
            sage: P.euler_characteristic()
            0
            sage: P.is_acyclic()
            False
        """
        if self.euler_characteristic() != 1:
            return False
        H = self.homology(base_ring=base_ring)
        if base_ring == ZZ:
            return all(len(x.invariants()) == 0 for x in H.values())
        else:
            # base_ring is a field.
            return all(x.dimension() == 0 for x in H.values())

    def n_chains(self, n, base_ring=None, cochains=False):
        r"""
        Return the free module of chains in degree ``n`` over ``base_ring``.

        INPUT:

        - ``n`` -- integer
        - ``base_ring`` -- ring (optional, default `\ZZ`)
        - ``cochains`` -- boolean (optional, default ``False``); if
          ``True``, return cochains instead

        The only difference between chains and cochains is
        notation. In a simplicial complex, for example, a simplex
        ``(0,1,2)`` is written as "(0,1,2)" in the group of chains but
        as "\chi_(0,1,2)" in the group of cochains.

        EXAMPLES::

            sage: S2 = simplicial_complexes.Sphere(2)
            sage: S2.n_chains(1, QQ)
            Free module generated by {(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)} over Rational Field
            sage: list(simplicial_complexes.Sphere(2).n_chains(1, QQ, cochains=False).basis())
            [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
            sage: list(simplicial_complexes.Sphere(2).n_chains(1, QQ, cochains=True).basis())
            [\chi_(0, 1), \chi_(0, 2), \chi_(0, 3), \chi_(1, 2), \chi_(1, 3), \chi_(2, 3)]
        """
        return Chains(tuple(self.n_cells(n)), base_ring, cochains)

    def algebraic_topological_model(self, base_ring=None):
        r"""
        Algebraic topological model for this cell complex with
        coefficients in ``base_ring``.

        The term "algebraic topological model" is defined by Pilarczyk
        and Réal [PR]_.

        This is not implemented for generic cell complexes. For any
        classes deriving from this one, when this method is
        implemented, it should essentially just call either
        :func:`~sage.homology.algebraic_topological_model.algebraic_topological_model`
        or
        :func:`~sage.homology.algebraic_topological_model.algebraic_topological_model_delta_complex`.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.algebraic_topological_model(QQ)
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    def homology_with_basis(self, base_ring=None, cohomology=False):
        r"""
        Return the unreduced homology of this complex with
        coefficients in ``base_ring`` with a chosen basis.

        This is implemented for simplicial, cubical, and
        `\Delta`-complexes, not for arbitrary generic cell complexes.

        INPUT:

        - ``base_ring`` -- coefficient ring (optional, default
          ``QQ``); must be a field
        - ``cohomology`` -- boolean (optional, default ``False``); if
          ``True``, return cohomology instead of homology

        Homology basis elements are named 'h_{dim,i}' where i ranges
        between 0 and `r-1`, if `r` is the rank of the homology
        group. Cohomology basis elements are denoted `h^{dim,i}`
        instead.

        .. SEEALSO::

            If ``cohomology`` is ``True``, this returns the cohomology
            as a graded module. For the ring structure, use
            :meth:`cohomology_ring`.

        EXAMPLES::

            sage: K = simplicial_complexes.KleinBottle()
            sage: H = K.homology_with_basis(QQ); H
            Homology module of Minimal triangulation of the Klein bottle
             over Rational Field
            sage: sorted(H.basis(), key=str)
            [h_{0,0}, h_{1,0}]
            sage: H = K.homology_with_basis(GF(2)); H
            Homology module of Minimal triangulation of the Klein bottle
             over Finite Field of size 2
            sage: sorted(H.basis(), key=str)
            [h_{0,0}, h_{1,0}, h_{1,1}, h_{2,0}]

        The homology is constructed as a graded object, so for
        example, you can ask for the basis in a single degree::

            sage: H.basis(1)
            Finite family {(1, 0): h_{1,0}, (1, 1): h_{1,1}}
            sage: S3 = delta_complexes.Sphere(3)
            sage: H = S3.homology_with_basis(QQ, cohomology=True)
            sage: list(H.basis(3))
            [h^{3,0}]
        """
        from homology_vector_space_with_basis import HomologyVectorSpaceWithBasis
        if base_ring is None:
            base_ring = QQ
        return HomologyVectorSpaceWithBasis(base_ring, self, cohomology)

    def cohomology_ring(self, base_ring=None):
        r"""
        Return the unreduced cohomology with coefficients in
        ``base_ring`` with a chosen basis.

        This is implemented for simplicial, cubical, and
        `\Delta`-complexes, not for arbitrary generic cell complexes.
        The resulting elements are suitable for computing cup
        products. For simplicial complexes, they should be suitable
        for computing cohomology operations; so far, only mod 2
        cohomology operations have been implemented.

        INPUT:

        - ``base_ring`` -- coefficient ring (optional, default
          ``QQ``); must be a field

        The basis elements in dimension ``dim`` are named 'h^{dim,i}'
        where `i` ranges between 0 and `r-1`, if `r` is the rank of
        the cohomology group.

        .. NOTE::

            For all but the smallest complexes, this is likely to be
            slower than :meth:`cohomology` (with field coefficients),
            possibly by several orders of magnitute. This and its
            companion :meth:`homology_with_basis` carry extra
            information which allows computation of cup products, for
            example, but because of speed issues, you may only wish to
            use these if you need that extra information.

        EXAMPLES::

            sage: K = simplicial_complexes.KleinBottle()
            sage: H = K.cohomology_ring(QQ); H
            Cohomology ring of Minimal triangulation of the Klein bottle
             over Rational Field
            sage: sorted(H.basis(), key=str)
            [h^{0,0}, h^{1,0}]
            sage: H = K.cohomology_ring(GF(2)); H
            Cohomology ring of Minimal triangulation of the Klein bottle
             over Finite Field of size 2
            sage: sorted(H.basis(), key=str)
            [h^{0,0}, h^{1,0}, h^{1,1}, h^{2,0}]

            sage: X = delta_complexes.SurfaceOfGenus(2)
            sage: H = X.cohomology_ring(QQ); H
            Cohomology ring of Delta complex with 3 vertices and 29 simplices
             over Rational Field
            sage: sorted(H.basis(1), key=str)
            [h^{1,0}, h^{1,1}, h^{1,2}, h^{1,3}]

            sage: H = simplicial_complexes.Torus().cohomology_ring(QQ); H
            Cohomology ring of Minimal triangulation of the torus
             over Rational Field
            sage: x = H.basis()[1,0]; x
            h^{1,0}
            sage: y = H.basis()[1,1]; y
            h^{1,1}

        You can compute cup products of cohomology classes::

            sage: x.cup_product(y)
            -h^{2,0}
            sage: x * y # alternate notation
            -h^{2,0}
            sage: y.cup_product(x)
            h^{2,0}
            sage: x.cup_product(x)
            0

        Cohomology operations::

            sage: RP2 = simplicial_complexes.RealProjectivePlane()
            sage: K = RP2.suspension()
            sage: K.set_immutable()
            sage: y = K.cohomology_ring(GF(2)).basis()[2,0]; y
            h^{2,0}
            sage: y.Sq(1)
            h^{3,0}

        To compute the cohomology ring, the complex must be
        "immutable". This is only relevant for simplicial complexes,
        and most simplicial complexes are immutable, but certain
        constructions make them mutable. The suspension is one
        example, and this is the reason for calling
        ``K.set_immutable()`` above. Another example::

            sage: S1 = simplicial_complexes.Sphere(1)
            sage: T = S1.product(S1)
            sage: T.is_immutable()
            False
            sage: T.cohomology_ring()
            Traceback (most recent call last):
            ...
            ValueError: This simplicial complex must be immutable. Call set_immutable().
            sage: T.set_immutable()
            sage: T.cohomology_ring()
            Cohomology ring of Simplicial complex with 9 vertices and
            18 facets over Rational Field
        """
        from homology_vector_space_with_basis import CohomologyRing
        if base_ring is None:
            base_ring = QQ
        return CohomologyRing(base_ring, self)

    @abstract_method
    def alexander_whitney(self, cell, dim_left):
        r"""
        The decomposition of ``cell`` in this complex into left and right
        factors, suitable for computing cup products. This should
        provide a cellular approximation for the diagonal map `K \to K
        \times K`.

        This method is not implemented for generic cell complexes, but
        must be implemented for any derived class to make cup products
        work in ``self.cohomology_ring()``.

        INPUT:

        - ``cell`` -- a cell in this complex
        - ``dim_left`` -- the dimension of the left-hand factors in
          the decomposition

        OUTPUT: a list containing triples ``(c, left, right)``.
        ``left`` and ``right`` should be cells in this complex, and
        ``c`` an integer. In the cellular approximation of the
        diagonal map, the chain represented by ``cell`` should get
        sent to the sum of terms `c (left \otimes right)` in the
        tensor product `C(K) \otimes C(K)` of the chain complex for
        this complex with itself.

        This gets used in the method
        :meth:`~sage.homology.homology_vector_space_with_basis.CohomologyRing.product_on_basis`
        for the class of cohomology rings.

        For simplicial and cubical complexes, the decomposition can be
        done at the level of individual cells: see
        :meth:`~sage.homology.simplicial_complex.Simplex.alexander_whitney`
        and
        :meth:`~sage.homology.cubical_complex.Cube.alexander_whitney`. Then
        the method for simplicial complexes just calls the method for
        individual simplices, and similarly for cubical complexes. For
        `\Delta`-complexes, the method is instead defined at the level
        of the cell complex.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.alexander_whitney(None, 2)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method alexander_whitney at ...>
        """

    ############################################################
    # end of chain complexes, homology
    ############################################################

    def face_poset(self):
        r"""
        The face poset of this cell complex, the poset of
        nonempty cells, ordered by inclusion.

        This uses the :meth:`cells` method, and also assumes that for
        each cell ``f``, all of ``f.faces()``, ``tuple(f)``, and
        ``f.dimension()`` make sense.  (If this is not the case in
        some derived class, as happens with `\Delta`-complexes, then
        override this method.)

        EXAMPLES::

            sage: P = SimplicialComplex([[0, 1], [1,2], [2,3]]).face_poset(); P
            Finite poset containing 7 elements
            sage: P.list()
            [(3,), (2,), (2, 3), (1,), (1, 2), (0,), (0, 1)]

            sage: S2 = cubical_complexes.Sphere(2)
            sage: S2.face_poset()
            Finite poset containing 26 elements
        """
        from sage.combinat.posets.posets import Poset
        from sage.misc.flatten import flatten
        covers = {}
        # The code for posets seems to work better if each cell is
        # converted to a tuple.
        all_cells = flatten([list(f) for f in self.cells().values()])

        for C in all_cells:
            if C.dimension() >= 0:
                covers[tuple(C)] = []
        for C in all_cells:
            for face in C.faces():
                if face.dimension() >= 0:
                    covers[tuple(face)].append(tuple(C))
        return Poset(covers)

    def graph(self):
        """
        The 1-skeleton of this cell complex, as a graph.

        This is not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.graph()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError

    @abstract_method
    def n_skeleton(self, n):
        """
        The `n`-skeleton of this cell complex: the cell
        complex obtained by discarding all of the simplices in
        dimensions larger than `n`.

        :param n: non-negative integer

        This is not implemented for general cell complexes.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: A = GenericCellComplex()
            sage: A.n_skeleton(3)
            Traceback (most recent call last):
            ...
            NotImplementedError: <abstract method n_skeleton at ...>
        """

    def _string_constants(self):
        """
        Tuple containing the name of the type of complex, and the
        singular and plural of the name of the cells from which it is
        built.  This is used in constructing the string representation.

        :return: tuple of strings

        This returns ``('Cell', 'cell', 'cells')``, as in "Cell
        complex", "1 cell", and "24 cells", but in other classes it
        could be overridden, as for example with ``('Cubical', 'cube',
        'cubes')`` or ``('Delta', 'simplex', 'simplices')``.  If for a
        derived class, the basic form of the print representation is
        acceptable, you can just modify these strings.

        EXAMPLES::

            sage: from sage.homology.cell_complex import GenericCellComplex
            sage: GenericCellComplex()._string_constants()
            ('Cell', 'cell', 'cells')
            sage: delta_complexes.Sphere(0)._string_constants()
            ('Delta', 'simplex', 'simplices')
            sage: cubical_complexes.Sphere(0)._string_constants()
            ('Cubical', 'cube', 'cubes')
        """
        return ('Cell', 'cell', 'cells')

    def _repr_(self):
        """
        Print representation.

        :return: string

        EXAMPLES::

            sage: delta_complexes.Sphere(7) # indirect doctest
            Delta complex with 8 vertices and 257 simplices
            sage: delta_complexes.Torus()._repr_()
            'Delta complex with 1 vertex and 7 simplices'
        """
        vertices = len(self.n_cells(0))
        Name, cell_name, cells_name = self._string_constants()
        if vertices != 1:
            vertex_string = "with %s vertices" % vertices
        else:
            vertex_string = "with 1 vertex"
        cells = 0
        for dim in self.cells():
            cells += len(self.cells()[dim])
        if cells != 1:
            cells_string = " and %s %s" % (cells, cells_name)
        else:
            cells_string = " and 1 %s" % cell_name
        return Name + " complex " + vertex_string + cells_string


class Chains(CombinatorialFreeModule):
    r"""
    Class for the free module of chains and/or cochains in a given
    degree.

    INPUT:

    - ``n_cells`` -- tuple of `n`-cells, which thus forms a basis for
      this module
    - ``base_ring`` -- optional (default `\ZZ`)
    - ``cochains`` -- boolean (optional, default ``False``); if
      ``True``, return cochains instead

    One difference between chains and cochains is notation. In a
    simplicial complex, for example, a simplex ``(0,1,2)`` is written
    as "(0,1,2)" in the group of chains but as "\chi_(0,1,2)" in the
    group of cochains.

    Also, since the free modules of chains and cochains are dual,
    there is a pairing `\langle c, z \rangle`, sending a cochain `c`
    and a chain `z` to a scalar.

    EXAMPLES::

        sage: S2 = simplicial_complexes.Sphere(2)
        sage: C_2 = S2.n_chains(1)
        sage: C_2_co = S2.n_chains(1, cochains=True)
        sage: x = C_2.basis()[Simplex((0,2))]
        sage: y = C_2.basis()[Simplex((1,3))]
        sage: z = x+2*y
        sage: a = C_2_co.basis()[Simplex((1,3))]
        sage: b = C_2_co.basis()[Simplex((0,3))]
        sage: c = 3*a-2*b
        sage: z
        (0, 2) + 2*(1, 3)
        sage: c
        -2*\chi_(0, 3) + 3*\chi_(1, 3)
        sage: c.eval(z)
        6
    """
    def __init__(self, n_cells, base_ring=None, cochains=False):
        """
        EXAMPLES::

            sage: T = cubical_complexes.Torus()
            sage: T.n_chains(2, QQ)
            Free module generated by {[1,1] x [0,1] x [1,1] x [0,1],
             [0,0] x [0,1] x [0,1] x [1,1], [0,0] x [0,1] x [1,1] x [0,1],
             [0,0] x [0,1] x [0,0] x [0,1], [0,1] x [1,1] x [0,1] x [0,0],
             [0,1] x [0,0] x [0,0] x [0,1], [1,1] x [0,1] x [0,1] x [0,0],
             [0,1] x [1,1] x [0,0] x [0,1], [0,0] x [0,1] x [0,1] x [0,0],
             [0,1] x [0,0] x [0,1] x [0,0], [0,1] x [0,0] x [1,1] x [0,1],
             [0,1] x [1,1] x [1,1] x [0,1], [0,1] x [0,0] x [0,1] x [1,1],
             [1,1] x [0,1] x [0,0] x [0,1], [1,1] x [0,1] x [0,1] x [1,1],
             [0,1] x [1,1] x [0,1] x [1,1]} over Rational Field
            sage: T.n_chains(2).dimension()
            16

        TESTS::

            sage: T.n_chains(2).base_ring()
            Integer Ring
            sage: T.n_chains(8).dimension()
            0
            sage: T.n_chains(-3).dimension()
            0
        """
        if base_ring is None:
            base_ring=ZZ
        self._cochains = cochains
        if cochains:
            CombinatorialFreeModule.__init__(self, base_ring, n_cells,
                                             prefix='\\chi', bracket=['_', ''])
        else:
            CombinatorialFreeModule.__init__(self, base_ring, n_cells,
                                             prefix='', bracket=False)

    class Element(CombinatorialFreeModuleElement):

        def eval(self, other):
            """
            Evaluate this cochain on the chain ``other``.

            INPUT:

            - ``other`` -- a chain for the same cell complex in the
              same dimension with the same base ring

            OUTPUT: scalar

            EXAMPLES::

                sage: S2 = simplicial_complexes.Sphere(2)
                sage: C_2 = S2.n_chains(1)
                sage: C_2_co = S2.n_chains(1, cochains=True)
                sage: x = C_2.basis()[Simplex((0,2))]
                sage: y = C_2.basis()[Simplex((1,3))]
                sage: z = x+2*y
                sage: a = C_2_co.basis()[Simplex((1,3))]
                sage: b = C_2_co.basis()[Simplex((0,3))]
                sage: c = 3*a-2*b
                sage: z
                (0, 2) + 2*(1, 3)
                sage: c
                -2*\chi_(0, 3) + 3*\chi_(1, 3)
                sage: c.eval(z)
                6

            TESTS::

                sage: z.eval(c) # z is not a cochain
                Traceback (most recent call last):
                ...
                ValueError: this element is not a cochain
                sage: c.eval(c) # can't evaluate a cochain on a cochain
                Traceback (most recent call last):
                ...
                ValueError: the elements are not compatible
            """
            if not self.parent()._cochains:
                raise ValueError('this element is not a cochain')
            if not (other.parent().indices() == self.parent().indices()
                    and other.base_ring() == self.base_ring()
                    and not other.parent()._cochains):
                raise ValueError('the elements are not compatible')
            result = sum(coeff * other.coefficient(cell) for cell, coeff in self)
            return result

