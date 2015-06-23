r"""
Algebra of differentiable scalar fields

The class :class:`DiffScalarFieldAlgebra` implements the commutative algebra
`C^k(U)` of scalar fields on some open subset `U` of a
topological manifold `M` over a topological field `K` (in most applications,
`K = \RR` or `K = \CC`). By *differentiable scalar field*, it
is meant a function of class C^`k` `U\rightarrow K`.
`C^k(U)` is an algebra over `K`, whose ring product is the pointwise
multiplication of `K`-valued functions, which is clearly commutative.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014-2015): initial version

REFERENCES:

- S. Kobayashi & K. Nomizu : *Foundations of Differential Geometry*, vol. 1,
  Interscience Publishers (New York) (1963)
- J.M. Lee : *Introduction to Smooth Manifolds*, 2nd ed., Springer (New York)
  (2013)
- B O'Neill : *Semi-Riemannian Geometry*, Academic Press (San Diego) (1983)

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

from sage.rings.infinity import infinity
from sage.manifolds.scalarfield_algebra import ScalarFieldAlgebra
from sage.manifolds.differentiable.scalarfield import DiffScalarField

class DiffScalarFieldAlgebra(ScalarFieldAlgebra):
    r"""
    Commutative algebra of differentiable scalar fields on some open subset of
    a differentiable manifold.

    If `M` is a topological manifold over a topological field `K` and `U`
    an open subset of `M`, the commutative algebra of scalar fields on `U`
    is the set `C^0(U)` of all continuous maps `U\rightarrow K`.
    `C^0(U)` is an algebra over `K`, whose ring product is the pointwise
    multiplication of `K`-valued functions, which is clearly commutative.

    If `K = \RR` or `K = \CC`, the field `K` over which the
    albegra `C^0(U)` is constructed,
    is represented by Sage's Symbolic Ring SR, since there is no exact
    representation of `\RR` nor `\CC` in Sage.

    The class :class:`ScalarFieldAlgebra` inherits from
    :class:`~sage.structure.parent.Parent`, with the category set to
    :class:`~sage.categories.commutative_algebras.CommutativeAlgebras`.
    The corresponding *element* class is
    :class:`~sage.manifolds.scalarfield.ScalarField`.

    INPUT:

    - ``domain`` -- the manifold open subset `U` on which the scalar fields are
      defined (must be an instance of class
      :class:`~sage.manifolds.manifold.TopManifold`)

    EXAMPLES:

    Algebras of scalar fields on the sphere `S^2` and on some subdomain of it::

        sage: M = TopManifold(2, 'M') # the 2-dimensional sphere S^2
        sage: U = M.open_subset('U') # complement of the North pole
        sage: c_xy.<x,y> = U.chart() # stereographic coordinates from the North pole
        sage: V = M.open_subset('V') # complement of the South pole
        sage: c_uv.<u,v> = V.chart() # stereographic coordinates from the South pole
        sage: M.declare_union(U,V)   # S^2 is the union of U and V
        sage: xy_to_uv = c_xy.transition_map(c_uv, (x/(x^2+y^2), y/(x^2+y^2)), \
                                             intersection_name='W', restrictions1= x^2+y^2!=0, \
                                             restrictions2= u^2+v^2!=0)
        sage: uv_to_xy = xy_to_uv.inverse()
        sage: CM = M.scalar_field_algebra() ; CM
        Algebra of scalar fields on the 2-dimensional topological manifold M
        sage: W = U.intersection(V)  # S^2 minus the two poles
        sage: CW = W.scalar_field_algebra() ; CW
        Algebra of scalar fields on the Open subset W of the 2-dimensional topological manifold M

    `C^0(M)` and `C^0(W)` belong to the category of commutative
    algebras over `\RR` (represented here by Sage's Symbolic Ring)::

        sage: CM.category()
        Category of commutative algebras over Symbolic Ring
        sage: CM.base_ring()
        Symbolic Ring
        sage: CW.category()
        Category of commutative algebras over Symbolic Ring
        sage: CW.base_ring()
        Symbolic Ring

    The elements of `C^0(M)` are scalar fields on `M`::

        sage: CM.an_element()
        Scalar field on the 2-dimensional topological manifold M
        sage: CM.an_element().display()  # this sample element is a constant field
        M --> R
        on U: (x, y) |--> 2
        on V: (u, v) |--> 2

    Those of `C^0(W)` are scalar fields on `W`::

        sage: CW.an_element()
        Scalar field on the Open subset W of the 2-dimensional topological manifold M
        sage: CW.an_element().display()  # this sample element is a constant field
        W --> R
        (x, y) |--> 2
        (u, v) |--> 2

    The zero element::

        sage: CM.zero()
        Scalar field zero on the 2-dimensional topological manifold M
        sage: CM.zero().display()
        zero: M --> R
        on U: (x, y) |--> 0
        on V: (u, v) |--> 0

    ::

        sage: CW.zero()
        Scalar field zero on the Open subset W of the 2-dimensional topological manifold M
        sage: CW.zero().display()
        zero: W --> R
           (x, y) |--> 0
           (u, v) |--> 0

    The unit element::

        sage: CM.one()
        Scalar field 1 on the 2-dimensional topological manifold M
        sage: CM.one().display()
        1: M --> R
        on U: (x, y) |--> 1
        on V: (u, v) |--> 1

    ::

        sage: CW.one()
        Scalar field 1 on the Open subset W of the 2-dimensional topological manifold M
        sage: CW.one().display()
        1: W --> R
        (x, y) |--> 1
        (u, v) |--> 1

    A generic element can be constructed as for any parent in Sage, namely
    by means of the ``__call__`` operator on the parent (here with the dictionary
    of the coordinate expressions defining the scalar field)::

        sage: f = CM({c_xy: atan(x^2+y^2), c_uv: pi/2 - atan(u^2+v^2)}); f
        Scalar field on the 2-dimensional topological manifold M
        sage: f.display()
        M --> R
        on U: (x, y) |--> arctan(x^2 + y^2)
        on V: (u, v) |--> 1/2*pi - arctan(u^2 + v^2)
        sage: f.parent()
        Algebra of scalar fields on the 2-dimensional topological manifold M

    Specific elements can also be constructed in this way::

        sage: CM(0) == CM.zero()
        True
        sage: CM(1) == CM.one()
        True

    Note that the zero scalar field is cached::

        sage: CM(0) is CM.zero()
        True

    Elements can also be constructed by means of the method
    :meth:`~sage.manifolds.manifold.TopManifold.scalar_field` acting on
    the domain (this allows one to set the name of the scalar field at the
    construction)::

        sage: f1 = M.scalar_field({c_xy: atan(x^2+y^2), c_uv: pi/2 - atan(u^2+v^2)}, name='f')
        sage: f1.parent()
        Algebra of scalar fields on the 2-dimensional topological manifold M
        sage: f1 == f
        True
        sage: M.scalar_field(0, chart='all') == CM.zero()
        True

    The algebra `C^0(M)` coerces to `C^0(W)` since `W` is an open
    subset of `M`::

        sage: CW.has_coerce_map_from(CM)
        True

    The reverse is of course false::

        sage: CM.has_coerce_map_from(CW)
        False

    The coercion map is nothing but the restriction to `W` of scalar fields
    on `M`::

        sage: fW = CW(f) ; fW
        Scalar field on the Open subset W of the 2-dimensional topological manifold M
        sage: fW.display()
        W --> R
        (x, y) |--> arctan(x^2 + y^2)
        (u, v) |--> 1/2*pi - arctan(u^2 + v^2)

    ::

        sage: CW(CM.one()) == CW.one()
        True

    The coercion map allows for the addition of elements of `C^0(W)`
    with elements of `C^0(M)`, the result being an element of
    `C^0(W)`::

        sage: s = fW + f
        sage: s.parent()
        Algebra of scalar fields on the Open subset W of the 2-dimensional topological manifold M
        sage: s.display()
        W --> R
        (x, y) |--> 2*arctan(x^2 + y^2)
        (u, v) |--> pi - 2*arctan(u^2 + v^2)

    Other coercions are those from the Symbolic Ring, leading to constant
    scalar fields (if the symbolic expression does not involve any chart
    coordinate)::

        sage: h = CM(pi*sqrt(2)) ; h
        Scalar field on the 2-dimensional topological manifold M
        sage: h.display()
        M --> R
        on U: (x, y) |--> sqrt(2)*pi
        on V: (u, v) |--> sqrt(2)*pi

    TESTS OF THE ALGEBRA LAWS:

    Ring laws::

        sage: s = f + h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> sqrt(2)*pi + arctan(x^2 + y^2)
        on V: (u, v) |--> 1/2*pi*(2*sqrt(2) + 1) - arctan(u^2 + v^2)

    ::

        sage: s = f - h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> -sqrt(2)*pi + arctan(x^2 + y^2)
        on V: (u, v) |--> -1/2*pi*(2*sqrt(2) - 1) - arctan(u^2 + v^2)

    ::

        sage: s = f*h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> sqrt(2)*pi*arctan(x^2 + y^2)
        on V: (u, v) |--> 1/2*sqrt(2)*(pi^2 - 2*pi*arctan(u^2 + v^2))

    ::

        sage: s = f/h ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> 1/2*sqrt(2)*arctan(x^2 + y^2)/pi
        on V: (u, v) |--> 1/4*sqrt(2)*(pi - 2*arctan(u^2 + v^2))/pi

    ::

        sage: f*(h+f) == f*h + f*f
        True

    Ring laws with coercion::

        sage: f - fW == CW.zero()
        True
        sage: f/fW == CW.one()
        True
        sage: s = f*fW ; s
        Scalar field on the Open subset W of the 2-dimensional topological manifold M
        sage: s.display()
        W --> R
        (x, y) |--> arctan(x^2 + y^2)^2
        (u, v) |--> 1/4*pi^2 - pi*arctan(u^2 + v^2) + arctan(u^2 + v^2)^2
        sage: s/f == fW
        True

    Multiplication by a real number::

        sage: s = 2*f ; s
        Scalar field on the 2-dimensional topological manifold M
        sage: s.display()
        M --> R
        on U: (x, y) |--> 2*arctan(x^2 + y^2)
        on V: (u, v) |--> pi - 2*arctan(u^2 + v^2)

    ::

        sage: 0*f == CM.zero()
        True
        sage: 1*f == f
        True
        sage: 2*(f/2) == f
        True
        sage: (f+2*f)/3 == f
        True
        sage: 1/3*(f+2*f) == f
        True

    Sage test suite for algebras is passed::

        sage: TestSuite(CM).run(verbose=True)
        running ._test_additive_associativity() . . . pass
        running ._test_an_element() . . . pass
        running ._test_associativity() . . . pass
        running ._test_category() . . . pass
        running ._test_characteristic() . . . pass
        running ._test_distributivity() . . . pass
        running ._test_elements() . . .
          Running the test suite of self.an_element()
          running ._test_category() . . . pass
          running ._test_eq() . . . pass
          running ._test_nonzero_equal() . . . pass
          running ._test_not_implemented_methods() . . . pass
          running ._test_pickling() . . . pass
          pass
        running ._test_elements_eq_reflexive() . . . pass
        running ._test_elements_eq_symmetric() . . . pass
        running ._test_elements_eq_transitive() . . . pass
        running ._test_elements_neq() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_one() . . . pass
        running ._test_pickling() . . . pass
        running ._test_prod() . . . pass
        running ._test_some_elements() . . . pass
        running ._test_zero() . . . pass

    It is passed also for `C^0(W)`::

        sage: TestSuite(CW).run()

    """

    Element = DiffScalarField

    def __init__(self, domain):
        r"""
        Construct an algebra of scalar fields.

        TESTS::

            sage: TopManifold._clear_cache_()  # for doctests only
            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: from sage.manifolds.scalarfield_algebra import ScalarFieldAlgebra
            sage: CM = ScalarFieldAlgebra(M); CM
            Algebra of scalar fields on the 2-dimensional topological manifold M
            sage: TestSuite(CM).run()

        """
        ScalarFieldAlgebra.__init__(self, domain)

    #### Methods required for any Parent


    def _coerce_map_from_(self, other):
        r"""
        Determine whether coercion to self exists from other parent

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: CM = M.scalar_field_algebra()
            sage: CM._coerce_map_from_(SR)
            True
            sage: U = M.open_subset('U', coord_def={X: x>0})
            sage: CU = U.scalar_field_algebra()
            sage: CM._coerce_map_from_(CU)
            False
            sage: CU._coerce_map_from_(CM)
            True

        """
        if other is SR:
            return True  # coercion from the base ring (multiplication by the
                         # algebra unit, i.e. self.one())
                         # cf. ScalarField._lmul_() for the implementation of
                         # the coercion map
        elif isinstance(other, DiffScalarFieldAlgebra):
            return self._domain.is_subset(other._domain)
        else:
            return False

    #### End of methods required for any Parent

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: CM = M.scalar_field_algebra()
            sage: CM._repr_()
            'Algebra of scalar fields on the 2-dimensional topological manifold M'
            sage: repr(CM)  # indirect doctest
            'Algebra of scalar fields on the 2-dimensional topological manifold M'

        """
        return "Algebra of differentiable scalar fields on the {}".format(self._domain)

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = TopManifold(2, 'M')
            sage: CM = M.scalar_field_algebra()
            sage: CM._latex_()
            'C^0 \\left(M\\right)'
            sage: latex(CM)  # indirect doctest
            C^0 \left(M\right)

        """
        degree = self._domain.diff_degree()
        if degree == infinity:
            latex_degree = r"\infty"  # to skip the "+" in latex(infinity)
        else:
            latex_degree = "{}".format(degree)
        return r"C^{" + latex_degree + r"}\left("  + self._domain._latex_() + \
               r"\right)"
