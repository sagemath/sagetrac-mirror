r"""
Kirillov-Reshetikhin Crystals
"""
#*****************************************************************************
#  Copyright (C) 2015   Travis Scrimshaw <tscrim at ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.category_singleton import Category_singleton
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.tensor import TensorProductsCategory
from sage.functions.other import ceil

class KirillovReshetikhinCrystals(Category_singleton):
    """
    The category of Kirillov-Reshetikhin crystals.

    EXAMPLES::

        sage: C = KirillovReshetikhinCrystals()
        sage: C
        Category of finite crystals
        sage: C.super_categories()
        [Category of crystals, Category of finite enumerated sets]
        sage: C.example()
        Highest weight crystal of type A_3 of highest weight omega_1

    TESTS::

        sage: TestSuite(C).run()
        sage: B = FiniteCrystals().example()
        sage: TestSuite(B).run(verbose = True)
    """
    def super_categories(self):
        r"""
        EXAMPLES::

            sage: WeylGroups().super_categories()
            [Category of coxeter groups]
        """
        return [FiniteCrystals(), RegularCrystals()]

    def additional_structure(self):
        r"""
        Return  ``None``.

        Indeed, the category of Weyl groups defines no additional
        structure: Weyl groups are a special class of Coxeter groups.

        .. SEEALSO:: :meth:`Category.additional_structure`

        .. TODO:: Should this category be a :class:`CategoryWithAxiom`?

        EXAMPLES::

            sage: WeylGroups().additional_structure()
        """
        return None

    def example(self):
        """
        Returns an example of highest weight crystals, as per
        :meth:`Category.example`.

        EXAMPLES::

            sage: B = FiniteCrystals().example(); B
            Highest weight crystal of type A_3 of highest weight omega_1
        """
        from sage.combinat.crystals.kirillov_reshetikhin import KirillovReshetikhinCrystal
        return KirillovReshetikhinCrystal(['A', 2, 1], 2, 2)

    class ParentMethods:

        def weight_lattice_realization(self):
            """
            Returns the weight lattice realization used to express weights.

            This default implementation uses the ambient space of the
            root system for (non relabelled) finite types and the
            weight lattice otherwise. This is a legacy from when
            ambient spaces were partially implemented, and may be
            changed in the future.

            EXAMPLES::

                sage: K = crystals.KirillovReshetikhin(['A',2,1], 1, 1)
                sage: K.weight_lattice_realization()
                Weight lattice of the Root system of type ['A', 2, 1]
            """
            return self.cartan_type().root_system().weight_lattice(extended=False)

    class TensorProducts(TensorProductsCategory):
        """
        The category of finite crystals constructed by tensor
        product of finite crystals.
        """
        class ElementMethods:
            def energy_function(self):
                r"""
                Return the energy function of ``self``.

                In this implementation, it is assumed that ``self`` is an
                element of a tensor product of perfect crystals of the same
                level, see Theorem 7.5 in [SchillingTingley2011]_.

                INPUT:

                - ``self`` -- an element of a tensor product of perfect
                  Kirillov-Reshetkhin crystals of the same level

                OUTPUT: an integer

                REFERENCES:

                .. [SchillingTingley2011] A. Schilling, P. Tingley.
                   Demazure crystals, Kirillov-Reshetikhin crystals, and
                   the energy function. Electronic Journal of Combinatorics.
                   **19(2)**. 2012. :arXiv:`1104.2359`

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:    print b, b.energy_function()
                    [[[1]], [[1]], [[1]]] 0
                    [[[1]], [[2]], [[1]]] 2
                    [[[2]], [[1]], [[1]]] 1
                    [[[3]], [[2]], [[1]]] 3

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,2)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:  # long time (5s on sage.math, 2011)
                    ....:     print b, b.energy_function()
                    [[], []] 4
                    [[], [[1, 1]]] 1
                    [[[1, 1]], []] 3
                    [[[1, 1]], [[1, 1]]] 0
                    [[[1, 2]], [[1, 1]]] 1
                    [[[2, 2]], [[1, 1]]] 2
                    [[[-1, -1]], [[1, 1]]] 2
                    [[[1, -1]], [[1, 1]]] 2
                    [[[2, -1]], [[1, 1]]] 2

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K)
                    sage: t = T.module_generators[0]
                    sage: t.energy_function()
                    Traceback (most recent call last):
                    ...
                    ValueError: All crystals in the tensor product need to be perfect of the same level
                """
                C = self.parent().crystals[0]
                ell = ceil(C.s() / C.cartan_type().c()[C.r()])
                if any(ell != K.s() / K.cartan_type().c()[K.r()] for K in self.parent().crystals):
                    raise ValueError("All crystals in the tensor product need to be perfect of the same level")
                t = self.parent()(*[K.module_generator() for K in self.parent().crystals])
                d = t.affine_grading()
                return d - self.affine_grading()

            def affine_grading(self):
                r"""
                Return the affine grading of `self`.

                The affine grading is calculated by finding a path from
                ``self`` to a ground state path using the helper method
                :meth:`e_string_to_ground_state` and counting the number of
                affine Kashiwara operators `e_0` applied on the way.

                OUTPUT: an integer

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.affine_grading()
                    1

                    sage: K = crystals.KirillovReshetikhin(['A',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:     print b, b.affine_grading()
                    [[[1]], [[1]], [[1]]] 3
                    [[[1]], [[2]], [[1]]] 1
                    [[[2]], [[1]], [[1]]] 2
                    [[[3]], [[2]], [[1]]] 0

                    sage: K = crystals.KirillovReshetikhin(['C',2,1],1,1)
                    sage: T = crystals.TensorProduct(K,K,K)
                    sage: hw = [b for b in T if all(b.epsilon(i)==0 for i in [1,2])]
                    sage: for b in hw:
                    ....:     print b, b.affine_grading()
                    [[[1]], [[1]], [[1]]] 2
                    [[[1]], [[2]], [[1]]] 1
                    [[[1]], [[-1]], [[1]]] 0
                    [[[2]], [[1]], [[1]]] 1
                    [[[-2]], [[2]], [[1]]] 0
                    [[[-1]], [[1]], [[1]]] 1
                """
                return self.e_string_to_ground_state().count(0)

            @cached_method
            def e_string_to_ground_state(self):
                r"""
                Return a string of integers in the index set
                `(i_1, \ldots, i_k)` such that `e_{i_k} \cdots e_{i_1}`
                of ``self`` is the ground state.

                This method calculates a path from ``self`` to a ground
                state path using Demazure arrows as defined in Lemma 7.3
                in [SchillingTingley2011]_.

                OUTPUT: a tuple of integers `(i_1,\ldots,i_k)`

                EXAMPLES::

                    sage: K = crystals.KirillovReshetikhin(['A',2,1], 1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]
                    sage: t.e_string_to_ground_state()
                    (0, 2)

                    sage: K = crystals.KirillovReshetikhin(['C',2,1], 1,1)
                    sage: T = crystals.TensorProduct(K,K)
                    sage: t = T.module_generators[0]; t
                    [[[1]], [[1]]]
                    sage: t.e_string_to_ground_state()
                    (0,)
                    sage: x=t.e(0)
                    sage: x.e_string_to_ground_state()
                    ()
                    sage: y=t.f_string([1,2,1,1,0]); y
                    [[[2]], [[1]]]
                    sage: y.e_string_to_ground_state()
                    ()
                """
                I = self.cartan_type().classical().index_set()
                ell = max(ceil(K.s()/K.cartan_type().c()[K.r()]) for K in self.parent().crystals)
                for i in I:
                    if self.epsilon(i) > 0:
                        return (i,) + (self.e(i)).e_string_to_ground_state()
                if self.epsilon(0) > ell:
                    return (0,) + (self.e(0)).e_string_to_ground_state()
                return ()

