r"""
Infinity Crystals
"""
#*****************************************************************************
#  Copyright (C) 2014    Ben Salisbury <salis1bt AT cmich DOT edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.categories.category_singleton import Category_singleton
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.highest_weight_crystals import HighestWeightCrystals

from sage.rings.infinity import Infinity


class InfinityCrystals(Category_singleton):
    r"""
    Category of infinity crystals; i.e., all models of crystal `B(\infty)`.

    EXMAPLES::

        sage: C = InfinityCrystals()
        sage: C.super_categories()
        [Category of infinite enumerated sets, Category of highest weight crystals]
        sage: C.example(n=2)
        The infinity crystal of tableaux of type ['A', 2]

    TESTS::

        sage: C = InfinityCrystals()
        sage: TestSuite(C).run()
        sage: B = C.example()
        sage: TestSuite(B).run(verbose=True)
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
        running ._test_enumerated_set_contains() . . . pass
        running ._test_enumerated_set_iter_cardinality() . . . pass
        running ._test_enumerated_set_iter_list() . . . pass
        running ._test_eq() . . . pass
        running ._test_not_implemented_methods() . . . pass
        running ._test_pickling() . . . pass
        running ._test_some_elements() . . . pass
    """

    def super_categories(self):
        r"""
        EXAMPLES::

            sage: InfinityCrystals().super_categories()
            [Category of infinite enumerated sets, Category of highest weight crystals]
        """
        return [InfiniteEnumeratedSets(), HighestWeightCrystals()]

    def example(self, affine=False, n=3):
        r"""
        Return an example of `B(\infty)`.

        EXAMPLES::

            sage: C = InfinityCrystals()
            sage: C.example()
            The infinity crystal of tableaux of type ['A', 3]
            sage: C.example(n=4)
            The infinity crystal of tableaux of type ['A', 4]
            sage: C.example(n=4,affine=True)
            Infinity Crystal of modified Nakajima monomials of type ['C', 3, 1]^*
        """
        from sage.combinat.crystals.infinity_crystals import InfinityCrystalOfTableaux
        from sage.combinat.crystals.monomial_crystals import InfinityCrystalOfNakajimaMonomials
        if affine == False:
            return InfinityCrystalOfTableaux(['A',n])
        else:
            return InfinityCrystalOfNakajimaMonomials(['D',n,2])


    class ParentMethods:

        def cardinality(self):
            r"""
            Returns the cardinality of the crystal: `\infty`.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("A3")
                sage: B.cardinality()
                +Infinity

                sage: M = crystals.infinity.NakajimaMonomials(['E',6,1])
                sage: M.cardinality()
                +Infinity

                sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
                sage: Y.cardinality()
                +Infinity
            """
            return Infinity


    class ElementMethods:

        def f_star(self,i):
            r"""
            Returns `f_i^*(x)` if it exists or ``None`` otherwise.

            EXAMPLES::

                sage: B = crystals.infinity.Tableaux("A3")
                sage: b = B.highest_weight_vector().f_string([1,2,1,1,2,3,3,3,3,2,1,2,3])
                sage: for i in B.index_set():
                ....:     print b.f_star(i)
                ....:
                [3, 2, 1, 4, 2, 1, 4, 2, 1, 4, 2, 1, 2, 1, 1, 2, 3, 3, 4, 4]
                [3, 2, 1, 4, 2, 1, 4, 2, 1, 2, 1, 4, 1, 1, 3, 3, 4, 4]
                [3, 2, 1, 4, 2, 1, 4, 2, 1, 4, 2, 1, 4, 2, 1, 2, 1, 1, 3, 3, 4, 4]
                sage: for i in B.index_set():
                    print b.f_star(i) == b.f(i)
                ....:
                True
                False
                True
            """
            from sage.combinat.crystals.elementary_crystals import ElementaryCrystal
            from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
            Binf = self.parent()
            t0 = Binf.highest_weight_vector()
            Bi = ElementaryCrystal(Binf.cartan_type(),i)
            b0 = Bi(0)
            tens = TensorProductOfCrystals(Bi,Binf)
            gen = tens(b0,t0)
            embedding = Binf.crystal_morphism({t0:gen})
            pullback = tens.crystal_morphism({gen:t0})
            image = embedding(self)
            return pullback(tens(image[0].f(i),image[1]))

        def e_star(self,i):
            r"""
            Returns `e_i^*(x)` if it exists or ``None`` otherwise.

            EXAMPLES::

                sage:
            """
            from sage.combinat.crystals.elementary_crystals import ElementaryCrystal
            from sage.combinat.crystals.tensor_product import TensorProductOfCrystals
            Binf = self.parent()
            t0 = Binf.highest_weight_vector()
            Bi = ElementaryCrystal(Binf.cartan_type(),i)
            b0 = Bi(0)
            tens = TensorProductOfCrystals(Bi,Binf)
            gen = tens(b0,t0)
            embedding = Binf.crystal_morphism({t0:gen})
            pullback = tens.crystal_morphism({gen:t0})
            image = embedding(self)
            if image[0].e(i) is not None:
                return pullback(tens(image[0].e(i),image[1]))
            else:
                return None

        def f_star_string(self, list):
            r"""
            Applies `f_{i_r}^* ... f_{i_1}^*` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage: Y = crystals.infinity.GeneralizedYoungWalls(3)
                sage: Y.highest_weight_vector().f_star_string([1,3,0,2,1,0])
                [[0], [1, 0, 3], [2, 1]]
            """
            b = self
            for i in list:
                b = b.f_star(i)
                if b is None:
                    return None
            return b

        def e_star_string(self, list):
            r"""
            Applies `e_{i_r}^* ... e_{i_1}^*` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage:
            """
            b = self
            for i in list:
                b = b.e_star(i)
                if b is None:
                    return None
            return b

        def epsilon_star(self,i):
            r"""
            Returns `\varepsilon_i^*(x)` where `x` is ``self``.

            EXAMPLES::

                sage:
            """
            e = 0
            b = self
            while b.e_star(i) is not None:
                e += 1
                b = b.e_star(i)
            return e

        def phi_star(self,i):
            r"""
            Returns `\varphi_i^*(x)` where `x` is ``self``.

            EXAMPLES::

                sage:
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return epsilon_star(self,i) + self.weight().scalar(h[i])
