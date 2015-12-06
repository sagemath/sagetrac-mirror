r"""
Tangent spaces

The class :class:`TangentSpace` implements tangent vector spaces to a
differentiable manifold.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- Chap. 3 of J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer
  (New York) (2013)

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

from sage.symbolic.ring import SR
from sage.tensor.modules.finite_rank_free_module import FiniteRankFreeModule
from sage.manifolds.differentiable.tangent_vector import TangentVector

class TangentSpace(FiniteRankFreeModule):
    r"""
    Tangent space to a differentiable manifold at a given point.

    Let `M` be a differentiable manifold of dimension `n` over a topological
    field `K`.  Since the tangent space at a given point `p\in M` is a
    `n`-dimensional vector space over `K` without any distinguished basis,
    the class :class:`TangentSpace` inherits from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`,
    which implements free modules of finite rank (hence vector spaces of finite
    dimension) without any distinguished basis.

    This is a Sage *parent* class, the corresponding *element* class being
    :class:`~sage.manifolds.differentiable.tangent_vector.TangentVector`.

    INPUT:

    - ``point`` -- (instance of
      :class:`~sage.manifolds.point.TopologicalManifoldPoint`) point `p` at which the
      tangent space is defined.

    EXAMPLES:

    Tangent space on a 2-dimensional manifold::

        sage: M = Manifold(2, 'M')
        sage: c_xy.<x,y> = M.chart()
        sage: p = M.point((-1,2), name='p')
        sage: Tp = M.tangent_space(p) ; Tp
        Tangent space at Point p on the 2-dimensional differentiable manifold M

    Tangent spaces actually belong to a dynamically generated subclass of
    :class:`TangentSpace`::

        sage: type(Tp)
        <class 'sage.manifolds.differentiable.tangent_space.TangentSpace_with_category'>

    They are free modules of finite rank over Sage's Symbolic Ring (actually
    vector spaces of finite dimension over the manifold base field `K`, with
    `K=\RR` here)::

        sage: isinstance(Tp, FiniteRankFreeModule)
        True
        sage: Tp.base_ring()
        Symbolic Ring
        sage: Tp.category()
        Category of vector spaces over Symbolic Ring
        sage: Tp.rank()
        2
        sage: dim(Tp)
        2

    The tangent space is automatically endowed with bases deduced from the
    vector frames around the point::

        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M]
        sage: M.frames()
        [Coordinate frame (M, (d/dx,d/dy))]

    At this stage, only one basis has been defined in the tangent space, but
    new bases can be added from vector frames on the manifold by means of the
    method :meth:`~sage.manifolds.differentiable.vectorframe.VectorFrame.at`,
    for instance, from the frame associated with some new coordinates::

        sage: c_uv.<u,v> = M.chart()
        sage: c_uv.frame().at(p)
        Basis (d/du,d/dv) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M
        sage: Tp.bases()
        [Basis (d/dx,d/dy) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M,
         Basis (d/du,d/dv) on the Tangent space at Point p on the 2-dimensional
         differentiable manifold M]

    All the bases defined on ``Tp`` are on the same footing. Accordingly the
    tangent space is not in the category of modules with a distinguished
    basis::

        sage: Tp in ModulesWithBasis(SR)
        False

    It is simply in the category of modules, as any instance of
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`::

        sage: Tp in Modules(SR)
        True

    Since the base ring is a field, it is actually in the category of
    vector spaces::

        sage: Tp in VectorSpaces(SR)
        True

    A typical element::

        sage: v = Tp.an_element() ; v
        Tangent vector at Point p on the 2-dimensional differentiable manifold M
        sage: v.display()
        d/dx + 2 d/dy
        sage: v.parent()
        Tangent space at Point p on the 2-dimensional differentiable manifold M

    The zero vector::

        sage: Tp.zero()
        Tangent vector zero at Point p on the 2-dimensional differentiable
         manifold M
        sage: Tp.zero().display()
        zero = 0
        sage: Tp.zero().parent()
        Tangent space at Point p on the 2-dimensional differentiable manifold M

    Tangent spaces are unique::

        sage: M.tangent_space(p) is Tp
        True
        sage: p1 = M.point((-1,2))
        sage: M.tangent_space(p1) is Tp
        True

    even if points are not::

        sage: p1 is p
        False

    Actually ``p1`` and ``p`` share the same tangent space because they
    compare equal::

        sage: p1 == p
        True

    This results from the unique representation behavior of the class
    :class:`TangentSpace`, which is inherited from
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`::

        sage: isinstance(Tp, UniqueRepresentation)
        True

    The tangent-space uniqueness holds even if the points are created in
    different coordinate systems::

        sage: xy_to_uv = c_xy.transition_map(c_uv, (x+y, x-y))
        sage: uv_to_xv = xy_to_uv.inverse()
        sage: p2 = M.point((1, -3), chart=c_uv, name='p_2')
        sage: p2 is p
        False
        sage: M.tangent_space(p2) is Tp
        True
        sage: p2 == p
        True

    See
    :class:`~sage.tensor.modules.finite_rank_free_module.FiniteRankFreeModule`
    for more documentation.

    """

    Element = TangentVector

    def __init__(self, point):
        r"""
        Construct the tangent space at a given point.

        TESTS::

            sage: from sage.manifolds.differentiable.tangent_space import TangentSpace
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = TangentSpace(p) ; Tp
            Tangent space at Point p on the 2-dimensional differentiable
             manifold M
            sage: TestSuite(Tp).run()

        """
        manif = point._manifold
        name = "T_" + str(point._name) + " " + str(manif._name)
        latex_name = r"T_{" + str(point._latex_name) + "}\," + \
                     str(manif._latex_name)
        self._point = point
        self._manif = manif
        FiniteRankFreeModule.__init__(self, SR, manif._dim, name=name,
                                      latex_name=latex_name,
                                      start_index=manif._sindex)
        # Initialization of bases of the tangent space from existing vector
        # frames around the point:
        self._frame_bases = {} # dictionary of bases of the tangent vector
                               # derived from vector frames (keys: vector
                               # frames)
        for frame in point._subset._top_frames:
            if point in frame._domain:
                basis = self.basis(symbol=frame._symbol,
                                   latex_symbol=frame._latex_symbol)
                # Names of basis vectors set to those of the frame vector
                # fields:
                n = manif._dim
                for i in range(n):
                    basis._vec[i]._name = frame._vec[i]._name
                    basis._vec[i]._latex_name = frame._vec[i]._latex_name
                basis._name = "(" + \
                        ",".join([basis._vec[i]._name for i in range(n)]) + ")"
                basis._latex_name = r"\left(" + \
                     ",".join([basis._vec[i]._latex_name for i in range(n)])+ \
                     r"\right)"
                basis._symbol = basis._name
                basis._latex_symbol = basis._latex_name
                # Names of cobasis linear forms set to those of the coframe
                # 1-forms:
                coframe = frame.coframe()
                cobasis = basis.dual_basis()
                for i in range(n):
                    cobasis._form[i]._name = coframe._form[i]._name
                    cobasis._form[i]._latex_name = coframe._form[i]._latex_name
                cobasis._name = "(" + \
                     ",".join([cobasis._form[i]._name for i in range(n)]) + ")"
                cobasis._latex_name = r"\left(" + \
                  ",".join([cobasis._form[i]._latex_name for i in range(n)])+ \
                  r"\right)"
                self._frame_bases[frame] = basis
        # The basis induced by the default frame of the manifold subset
        # in which the point has been created is declared the default
        # basis of self:
        def_frame = point._subset._def_frame
        if def_frame in self._frame_bases:
            self._def_basis = self._frame_bases[def_frame]
        # Initialization of the changes of bases from the existing changes of
        # frames around the point:
        for frame_pair, automorph in point._subset._frame_changes.iteritems():
            if point in automorph.domain():
                frame1 = frame_pair[0] ; frame2 = frame_pair[1]
                fr1, fr2 = None, None
                for frame in self._frame_bases:
                    if frame1 in frame._subframes:
                        fr1 = frame
                        break
                for frame in self._frame_bases:
                    if frame2 in frame._subframes:
                        fr2 = frame
                        break
                if fr1 is not None and fr2 is not None:
                    basis1 = self._frame_bases[fr1]
                    basis2 = self._frame_bases[fr2]
                    auto = self.automorphism()
                    for frame, comp in automorph._components.iteritems():
                        basis = None
                        if frame is frame1:
                            basis = basis1
                        if frame is frame2:
                            basis = basis2
                        if basis is not None:
                            cauto = auto.add_comp(basis)
                            for ind, val in comp._comp.iteritems():
                                cauto._comp[ind] = val(point)
                    self._basis_changes[(basis1, basis2)] = auto

    def _repr_(self):
        r"""
        String representation of the object.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp._repr_()
            'Tangent space at Point p on the 2-dimensional differentiable manifold M'
            sage: repr(Tp)  # indirect doctest
            'Tangent space at Point p on the 2-dimensional differentiable manifold M'

        """
        description = "Tangent space at {}".format(self._point)
        return description

    def _an_element_(self):
        r"""
        Construct some (unnamed) vector in the tangent space

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((3,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp._an_element_()
            Tangent vector at Point p on the 2-dimensional differentiable
             manifold M
            sage: Tp._an_element_().display()
            d/dx + 2 d/dy

        """
        resu = self.element_class(self)
        if self._def_basis is not None:
            resu.set_comp()[:] = range(1, self._rank+1)
        return resu

    def dimension(self):
        r"""
        Return the vector space dimension of the tangent space.

        EXAMPLE::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp.dimension()
            2

        A shortcut is ``dim()``::

            sage: Tp.dim()
            2

        One can also use the global function ``dim``::

            sage: dim(Tp)
            2

        """
        # The dimension is the rank of self as a free module:
        return self._rank

    dim = dimension

    def base_point(self):
        r"""
        Return the manifold point at which the tangent space is defined.

        EXAMPLE::

            sage: from sage.manifolds.differentiable.tangent_space import TangentSpace # for doctests only
            sage: TangentSpace._clear_cache_()  # for doctests only
            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: p = M.point((1,-2), name='p')
            sage: Tp = M.tangent_space(p)
            sage: Tp.base_point()
            Point p on the 2-dimensional differentiable manifold M
            sage: Tp.base_point() is p
            True

        """
        return self._point
