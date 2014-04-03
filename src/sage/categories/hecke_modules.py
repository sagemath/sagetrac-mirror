r"""
Hecke modules
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2005      David Kohel <kohel@maths.usyd.edu>
#                          William Stein <wstein@math.ucsd.edu>
#                2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_types import Category_module
from sage.misc.cachefunc import cached_method
from sage.categories.category import HomCategory
from sage.categories.all import ModulesWithBasis

class HeckeModules(Category_module):
    r"""
    The category of Hecke modules.

    A Hecke module is a module `M` over the \emph{anemic} Hecke
    algebra, i.e., the Hecke algebra generated by Hecke operators
    `T_n` with `n` coprime to the level of `M`.  (Every Hecke module
    defines a level function, which is a positive integer.)  The
    reason we require that `M` only be a module over the anemic Hecke
    algebra is that many natural maps, e.g., degeneracy maps,
    Atkin-Lehner operators, etc., are `\Bold{T}`-module homomorphisms; but
    they are homomorphisms over the anemic Hecke algebra.

    EXAMPLES:

    We create the category of Hecke modules over `\QQ`::

        sage: C = HeckeModules(RationalField()); C
        Category of Hecke modules over Rational Field

    TODO: check that this is what we want::

        sage: C.super_categories()
        [Category of modules with basis over Rational Field]

    # [Category of vector spaces over Rational Field]

    Note that the base ring can be an arbitrary commutative ring::

        sage: HeckeModules(IntegerRing())
        Category of Hecke modules over Integer Ring
        sage: HeckeModules(FiniteField(5))
        Category of Hecke modules over Finite Field of size 5

    The base ring doesn't have to be a principal ideal domain::

        sage: HeckeModules(PolynomialRing(IntegerRing(), 'x'))
        Category of Hecke modules over Univariate Polynomial Ring in x over Integer Ring

    TESTS::

        sage: TestSuite(HeckeModules(ZZ)).run()
    """
    def __init__(self, R):
        """
        TESTS::

            sage: TestSuite(HeckeModules(ZZ)).run()

            sage: HeckeModules(Partitions(3)).run()
            Traceback (most recent call last):
            ...
            TypeError: R (=Partitions of the integer 3) must be a commutative ring
        """
        from .commutative_rings import CommutativeRings
        if R not in CommutativeRings():
            raise TypeError("R (=%s) must be a commutative ring"%R)
        Category_module.__init__(self, R, "Hecke modules")

    def super_categories(self):
        """
        EXAMPLES::

            sage: HeckeModules(QQ).super_categories()
            [Category of modules with basis over Rational Field]
        """
        R = self.base_ring()
        return [ModulesWithBasis(R)]


    class ParentMethods:

        def _Hom_(self, Y, category):
            r"""
            Returns the homset from ``self`` to ``Y`` in the category ``category``

            INPUT::

            - ``Y`` -- an Hecke module
            - ``category`` -- a subcategory of :class:`HeckeModules`() or None

            The sole purpose of this method is to construct the homset
            as a :class:`~sage.modular.hecke.homspace.HeckeModuleHomspace`. If
            ``category`` is specified and is not a subcategory of
            :class:`HeckeModules`, a ``TypeError`` is raised instead

            This method is not meant to be called directly. Please use
            :func:`sage.categories.homset.Hom` instead.

            EXAMPLES::

                sage: M = ModularForms(Gamma0(7), 4)
                sage: H = M._Hom_(M, category = HeckeModules(QQ)); H
                Set of Morphisms from Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) of weight 4 over Rational Field to Modular Forms space of dimension 3 for Congruence Subgroup Gamma0(7) of weight 4 over Rational Field in Category of Hecke modules over Rational Field
                sage: H.__class__
                <class 'sage.modular.hecke.homspace.HeckeModuleHomspace_with_category'>
                sage: TestSuite(H).run(skip=["_test_zero", "_test_elements", "_test_an_element", "_test_additive_associativity", "_test_elements_eq", "_test_elements_eq_reflexive", "_test_elements_eq_transitive", "_test_elements_eq_symmetric", "_test_elements_neq", "_test_some_elements"])

            Fixing :meth:`_test_zero` (``__call__`` should accept a
            function as input) and :meth:`_test_elements*` (modular
            form morphisms elements should inherit from categories) is
            :trac:`12879`.

            TESTS::

                sage: H = M._Hom_(M, category = HeckeModules(GF(5))); H
                Traceback (most recent call last):
                ...
                TypeError: Category of Hecke modules over Finite Field of size 5 is not a subcategory of Category of Hecke modules over Rational Field

            """
            # TODO: double check that it's the correct HeckeModules category below:
            if category is not None and not category.is_subcategory(HeckeModules(self.base_ring())):
                raise TypeError("%s is not a subcategory of %s"%(category, HeckeModules(self.base_ring())))
            from sage.modular.hecke.homspace import HeckeModuleHomspace
            return HeckeModuleHomspace(self, Y, category = category)

    class HomCategory(HomCategory):
        def extra_super_categories(self):
            """
            EXAMPLES::

                sage: HeckeModules(ZZ).hom_category().extra_super_categories()
                []
            """
            return [] # FIXME: what category structure is there on Homsets of hecke modules?


        def base_ring(self):
            """
            EXAMPLES::

                sage: HeckeModules(QQ).hom_category().base_ring()
                Rational Field
            """
            return self.base().base_ring()

        class ParentMethods:
            pass
