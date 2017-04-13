"""
Homomorphisms of Lie Algebras

AUTHORS:

- Travis Scrimshaw (07-15-2013): Initial implementation
"""
#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.morphism import Morphism
from sage.categories.morphism import SetMorphism
from sage.categories.homset import Homset
from sage.structure.sequence import Sequence, Sequence_generic

# TODO: Refactor out common functionality with RingHomomorphism_im_gens
class LieAlgebraHomomorphism_im_gens(Morphism):
    r"""
    A homomorphism of Lie algebras.

    Let `\mathfrak{g}` and `\mathfrak{g}^{\prime}` be Lie algebras,
    a linear map `f \colon \mathfrak{g} \to \mathfrak{g}^{\prime}` is a
    homomorphism (of Lie algebras) if `f([x, y]) = [f(x), f(y)]` for all
    `x, y \in \mathfrak{g}`. Thus homomorphisms are completely determined
    by the image of the generators of `\mathfrak{g}`.

    EXAMPLES::

        sage: L = LieAlgebra(QQ, 'x,y,z')
        sage: Lyn = L.Lyndon()
        sage: H = L.Hall()
        sage: phi = Lyn.coerce_map_from(H); phi
        Lie algebra morphism:
          From: Free Lie algebra generated by (x, y, z) over Rational Field in the Hall basis
          To:   Free Lie algebra generated by (x, y, z) over Rational Field in the Lyndon basis
          Defn: x |--> x
                y |--> y
                z |--> z
    """
    def __init__(self, parent, im_gens, check=True):
        """
        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: phi = Lyn.coerce_map_from(H)

        We skip the category test because the Homset's element class
        does not match this class::

            sage: TestSuite(phi).run(skip=['_test_category'])
        """
        Morphism.__init__(self, parent)
        if not isinstance(im_gens, Sequence_generic):
            if not isinstance(im_gens, (tuple, list)):
                im_gens = [im_gens]
            im_gens = Sequence(im_gens, parent.codomain(), immutable=True)
        if check:
            if len(im_gens) != len(parent.domain().lie_algebra_generators()):
                raise ValueError("number of images must equal number of generators")
            # TODO: Implement a (meaningful) _is_valid_homomorphism_()
            #if not parent.domain()._is_valid_homomorphism_(parent.codomain(), im_gens):
            #    raise ValueError("relations do not all (canonically) map to 0 under map determined by images of generators.")
        if not im_gens.is_immutable():
            import copy
            im_gens = copy.copy(im_gens)
            im_gens.set_immutable()
        self.__im_gens = im_gens

    def _repr_type(self):
        """
        TESTS::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: f = Lyn.coerce_map_from(H)
            sage: type(f)
            <class 'sage.algebras.lie_algebras.morphism.LieAlgebraHomomorphism_im_gens'>
            sage: f._repr_type()
            'Lie algebra'
        """
        return "Lie algebra"

    def im_gens(self):
        """
        Return the images of the generators of the domain.

        OUTPUT:

        - ``list`` -- a copy of the list of gens (it is safe to change this)

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: f = Lyn.coerce_map_from(H)
            sage: f.im_gens()
            [x, y, z]
        """
        return list(self.__im_gens)

    def __hash__(self):
        """
        Return the hash of this morphism.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: phi = Lyn.coerce_map_from(H)
            sage: hash(phi) == hash(phi)
            True
        """
        return hash(self.__im_gens)

    def _repr_defn(self):
        """
        Used in constructing string representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: phi = Lyn.coerce_map_from(H)
            sage: print phi._repr_defn()
            x |--> x
            y |--> y
            z |--> z
        """
        D = self.domain()
        ig = self.__im_gens
        return '\n'.join(['%s |--> %s'%(x, ig[i]) for i, x in enumerate(D.gens())])

    def _call_(self, x):
        """
        Evaluate this homomorphism at ``x``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: phi = Lyn.coerce_map_from(H)
            sage: a = H.graded_basis(5)[12]; a
            [z, [y, [x, [x, z]]]]
            sage: phi(a)
            [x, [[x, z], [y, z]]] + [x, [[[x, z], z], y]]
             + [[x, y], [[x, z], z]] + [[x, [y, z]], [x, z]]
        """
        return x._im_gens_(self.codomain(), self.im_gens())

class LieAlgebraHomset(Homset):
    """
    Homset between two Lie algebras.

    .. TODO::

        This is a very minimal implementation which does not
        have coercions of the morphisms.
    """
    def __init__(self, X, Y, category=None, base=None, check=True):
        """
        Initalize ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: HS = Hom(Lyn, H)

        We skip the elements test since homsets are not proper parents::

            sage: TestSuite(HS).run(skip=["_test_elements"])
        """
        if base is None:
            base = X.base_ring()
        Homset.__init__(self, X, Y, category, base, check)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: Hom(Lyn, H)
            Set of Lie algebra morphisms
             from Free Lie algebra generated by (x, y, z) over Rational Field in the Lyndon basis
             to Free Lie algebra generated by (x, y, z) over Rational Field in the Hall basis
        """
        return "Set of Lie algebra morphisms from {} to {}".format(
                                self.domain(), self.codomain())

    def __call__(self, im_gens, check=True):
        """
        Construct an morphism from ``im_gens``.

        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: HS = Hom(Lyn, H)
            sage: x,y,z = Lyn.gens()
            sage: phi = HS([z,x,y])
            sage: phi(x + Lyn[z,y])
            z - [x, y]

            sage: HS = Hom(Lyn, Lyn)
            sage: phi = HS([z,x,y])
            sage: a = Lyn([z, [x, [[y, z], x]]]); a
            [x, [x, [[y, z], z]]] + [x, [[x, z], [y, z]]] - [[x, [y, z]], [x, z]]
            sage: phi(a)
            [[x, [[y, z], z]], y] + 2*[[[x, z], [y, z]], y] + [[[[x, z], z], y], y]
            sage: phi(phi(phi(a))) == a
            True
        """
        if isinstance(im_gens, Morphism):
            return im_gens
        return LieAlgebraHomomorphism_im_gens(self, im_gens)

    @cached_method
    def zero(self):
        """
        Return the zero morphism.
        EXAMPLES::

            sage: L = LieAlgebra(QQ, 'x,y,z')
            sage: Lyn = L.Lyndon()
            sage: H = L.Hall()
            sage: HS = Hom(Lyn, H)
            sage: HS.zero()
            Generic morphism:
              From: Free Lie algebra generated by (x, y, z) over Rational Field in the Lyndon basis
              To:   Free Lie algebra generated by (x, y, z) over Rational Field in the Hall basis
        """
        return SetMorphism(self, lambda x: self.codomain().zero())

    _an_element_ = zero

