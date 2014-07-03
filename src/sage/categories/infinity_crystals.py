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
from sage.combinat.crystals.elementary_crystals import ElementaryCrystal as Elementary
from sage.combinat.crystals.tensor_product import TensorProductOfCrystals as TensorProduct

from sage.rings.infinity import Infinity


class InfinityCrystals(Category_singleton):
    r"""
    Category of infinity crystals.

    EXMAPLES::

        sage:
    """
    def super_categories(self):
        r"""
        EXAMPLES::

            sage:
        """
        return [InfiniteEnumeratedSets(), HighestWeightCrystals()]


    class ParentMethods:

        def cardinality(self):
            r"""
            Returns the cardinality of the crystal: `\infty`.

            EXAMPLES::
            """
            return Infinity


    class ElementMethods:

        def f_star(self,i):
            r"""
            Returns `f_i^*(x)` if it exists or ``None`` otherwise.

            EXAMPLES::

                sage:
            """
            Binf = self.parent()
            t0 = Binf.highest_weight_vector()
            Bi = crystals.elementary.Elementary(Binf.cartan_type(),i)
            b0 = Bi(0)
            tens = crystals.TensorProduct(Bi,Binf)
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
            Binf = self.parent()
            t0 = Binf.highest_weight_vector()
            Bi = crystals.elementary.Elementary(Binf.cartan_type(),i)
            b0 = Bi(0)
            tens = crystals.TensorProduct(Bi,Binf)
            gen = tens(b0,t0)
            embedding = Binf.crystal_morphism({t0:gen})
            pullback = tens.crystal_morphism({gen:t0})
            image = embedding(self)
            return pullback(tens(image[0].e(i),image[1]))

        def f_star_string(self, list):
            r"""
            Applies `f_{i_r}^* ... f_{i_1}^*` to self for `list = [i_1, ..., i_r]`

            EXAMPLES::

                sage:
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
            Returns `\varepsilon_i(x)` where `x` is ``self``.

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
            Returns `\varphi_i(x)` where `x` is ``self``.

            EXAMPLES::

                sage:
            """
            P = self.parent().weight_lattice_realization()
            h = P.simple_coroots()
            return epsilon_star(self,i) + self.weight().scalar(h[i])
