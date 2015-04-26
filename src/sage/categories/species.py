# -*- coding: utf-8 -*-
"""
Species category

References
----------

.. [BBL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle and Pierre Leroux,
  1998, Cambridge University Press

"""
#*******************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>,
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*******************************************************************************
from sage.categories.category import Category
from sage.categories.objects import Objects
from sage.misc.abstract_method import abstract_method
from sage.rings.infinity import Infinity


class Species(Category):
    """
    The category of species (of structures).


    A *species* (of structures) is a rule `F` which

     - produces, for each finite set `U`, a finite set `F[U]`,
     - produces, for each bijection `\sigma : U \to V`, a function

    MATH::

        F[\sigma] : F[U] \to F[V]\,.

    The functions `F[\sigma]` should further satisfy the following functorial properties:

     - for all bijections `\sigma : U \to V` and `\tau : V \to W`,

    MATH::

        F[\tau \circ \sigma] = F[\tau] \circ F[\sigma]\,,

     - for the identity map `Id_U : U \to U`,

    MATH::

        F[Id_U] = Id_{F[U]}\,.

    (Definition 3, section 1.1, _[BBL])
    """

    def super_categories(self):
        return [Objects()]

    def one(self):
        """
        The species `1`, characteristic of the *empty set*, defined by

        MATH::

            1[U] := \begin{dcases*}
                \{U\} & if `U = \emptyset`,\\
                \emptyset & otherwise,
            \end{dcases*}

        for any finite set `U`.
        """
        from sage.combinat.species2.one import OneSpecies
        return OneSpecies()

    def zero(self):
        """
        The species `0` defined by

        MATH::

            0[U] := \emptyset

        for any finite set `U`.
        """
        from sage.combinat.species2.zero import ZeroSpecies
        return ZeroSpecies()

    def singletons(self):
        """
        The species `X`, characteristic of *singletons*

        MATH::

            X[U] = \begin{dcases*}
                \{U\} & if `|U| = 1`,\\
                \emptyset & otherwise,
            \end{dcases*}

        for any finite set `U`.

        """
        from sage.combinat.species2.singletons import SingletonsSpecies
        return SingletonsSpecies()

    def sets(self):
        """
        The species `E`, of *sets*, defined by

        MATH::

            E[U] = \{U\}

        for any finite set `U`.
        """
        from sage.combinat.species2.sets import SetsSpecies
        return SetsSpecies()

    def recursive_species(self, name="F"):
        """
        Return an (not defined) instance of a recursive species

        :param name: a string

        EXAMPLES::

            sage: Sp = Species()
            sage: O, X = Sp.one(), Sp.singletons()
            sage: B = Sp.recursive_species(name="B")
            sage: B.define(O + B*X*B); B

            sage: list(B.type_generating_series().coefficients(10))
        """
        from sage.combinat.species2.operations.recursive_species import RecursiveSpecies
        return RecursiveSpecies(name=name)

    class ParentMethods:

        @abstract_method
        def structures(self, U):
            """
            The *`F`-structures* on the finite set `U`.

            :param U: a finite set

            An element `s \in F[U]` is called an *`F`-structure* on `U`.
            (Definition 3, section 1.1, _[BBL])
            """

        @abstract_method
        def transport(self, sigma):
            """
            The *transport* of `F`-structures along `\sigma`.

            :param sigma: a bijection `U \to V`

            The function `F[\sigma]` is called the *transport* of `F`-structures along `\sigma`.
            (Definition 3, section 1.1, _[BBL])

            We follow the _[BBL] notation: `\sigma \cdot s` to designate `F[\sigma](s)` the transport of `F`-structure
            along `\sigma`.
            """

        def is_isomorphism_of(self, sigma, s1, s2):
            """
            Return if `\sigma` is an isomorphism of `s_1` to `s_2`.
            (Return if `s_1` and `s_2` have same *isomorphism type*)

            :param sigma: a bijection `U \to V`
            :param s1: a `F`-structure on `U`
            :param s2: a `F`-structure on `V`

            The bijection `\sigma` is an *isomorphism* of `s_1` to `s_2` if `s_2 = \sigma \cdot s_1 = F[\sigma](s_1)`.
            (Definition 4, section 1.1, _[BBL])
            """
            return self.transport(sigma)(s1) == s2

        def is_automorphism_of(self, sigma, s):
            """
            Return if `\sigma` is an automorphism of `s`.

            :param sigma: a bijection `U \to V`
            :param s: a `F`-structure on `U`

            An *automorphism* of `s` is an isomorphism of `s` to `s`.
            (Definition 4, section 1.1, _[BBL])
            """
            return self.is_isomorphism_of(sigma, s, s)

        def isomorphism_types(self, n=None):
            """
            The *isomorphism types* of `F`-structures of order `n`.

            :param n: a non-negative integer

            Let `\sim` be the equivalence relation on `F[n]` defined by `s \sim t` *iff* `s` and `t` have same
            isomorphism type.

            An *isomorphism type* of `F`-structures of order `n` is an equivalence class (modulo `\sim`) of
            `F`-structure on `[n]` (with `[n] := \{1, \cdots, n\}`).

            (see section 1.2 Type generating series, _[BBL])

            This method return the quotient `T(F[n]) := F[n]/\sim`.
            """
            from sage.combinat.species2 import ClassOfIsoTypes
            if n != None:
                return ClassOfIsoTypes(self).graded_component(n)
            return ClassOfIsoTypes(self)

        def Fix(self, sigma):
            """
            The set of fixed point of the transport of `F`-structures along `\sigma`.

            :param sigma: a permutation on $U$

            Let `\sigma` be a permutation on `U`.

            MATH::

                \mathtt{Fix}\, F[\sigma] := \left\{ s \in F[U] \,\mid\, \sigma \cdot s = s \right\}\,.

            (Definition 5, section 1.2, _[BBL])
            """
            assert(sigma.domain() == sigma.codomain()), "`sigma` should be a bijection U -> U"
            return (s for s in self.structures(sigma.domain())
                    if s == self.transport(sigma)(s))

        def fix(self, sigma):
            """
            The number of fixed point of the transport of `F`-structures along `\sigma`.

            ..see:: `:meth:sage.categories.species.Fix`

            MATH::

                \mathtt{fix}\, F[\sigma] := |\mathtt{Fix}\, F[\sigma]|\,.

            (Definition 5, section 1.2, _[BBL])
            """
            return len(self.Fix(sigma))

        ###################################################################
        #####                   Associated series                    ######
        ###################################################################

        def generating_series(self):
            """
            The *(exponential) generating series* of `F`.

            The generating series of a species `F` is a formal power series

            MATH::

                F(x) = \sum_{n \geqslant 0} f_n \frac{x^n}{n!}

            where `f_n` is the number of `F`-structures on a set of cardinality `n`.
            (Definition 1, section 1.2, _[BBL])

            Theorem 8 _[BBL]:

            MATH::

                \tilde{F}(x) = Z_F(x, 0, 0, \cdots)\,.

            """
            return self.cycle_index_series().generating_series()

        def exponential_generating_series(self): return self.generating_series()
        def egs(self): return self.generating_series()

        def type_generating_series(self):
            """
            The *(isomorphism) type generating series* of `F`.

            It is the formal power series

            MATH::

                \tilde{F}(x) = \sum_{n \geqslant 0} \tilde{f_n} x^n\,,

            where `\tilde{f_n}` is the number of *isomorphism type* of `F`-structures of order `n`.
            (Definition 3, section 1.2, _[BBL])

            Theorem 8 _[BBL]:

            MATH::

                \tilde{F}(x) = Z_F(x, x^2, x^3, \cdots)\,.

            """
            return self.cycle_index_series().isomorphism_type_generating_series()

        def ordinary_generating_series(self): return self.type_generating_series()
        def ogs(self): return self.type_generating_series()
        def isomorphism_type_generating_series(self): return self.type_generating_series()

        def cycle_index_series(self):
            """
            The *cycle index series* of `F`.

            The *cycle index series* of `F` is a formal power series (in symmetric functions)

            MATH::

                Z_F(p_1, p_2, \cdots) = \sum_{n \geqslant 0} \frac{1}{n!} \left(
                    \sum_{\sigma \in S_n} \mathtt{fix}\, F[\sigma]\, p_1^{\sigma_1} p_2^{\sigma_2} \cdots
                \right) \,.

            where `p_i` denotes the power sum symmetric functions, `S_n` denotes the group of permutations of `[n]` and
            `\mathtt{fix}\, F[\sigma]` the number of `F`-structures on `[n]` fixed by `F[\sigma]`.
            (Definition 6 and Remark 10, section 1.2, _[BBL])

            TESTS::

                sage: P = Permutations()
                sage: P.cis = Species().ParentMethods.cycle_index_series.f(P)
                sage: for n in range(5):
                ....:     assert(P.cis().Frobenius_characteristic(n) == P.cycle_index_series().Frobenius_characteristic(n))

            """
            from sage.combinat.species2.cycle_index_series.misc import genericCIS
            return genericCIS(self)

        def cis(self): return self.cycle_index_series()

        ###################################################################
        #####                   Operation                            ######
        ###################################################################

        def add(self, G):
            """
            The *sum* of `F (:= self)` and `G`.

            :param G: a species

            The *sum* of `F` and `G`, noted `F + G`, is defined as follows: an `(F+G)`-structure on `U` is an
            `F`-structures or (exclusive) a `G`-structures on `U`:

            MATH::

                (F + G)[U] = F[U] \sqcup G[U]\,.

            The transport of structures along a bijection `\sigma : U \to V` is carried out by setting

            MATH::

                (F + G)[\sigma](s) = \begin{dcases*}
                    F[\sigma](s) & if `s \in F[U]`,\\
                    G[\sigma](s) & if `s \in G[U]`,
                \end{dcases*}

            for each `(F + G)`-structure `s`.
            (section 1.3, _[BBL])
            """
            from sage.combinat.species2.operations.add import Add
            return Add(self, G)

        __add__ = add

        # TODO use the method sum to define some tricky parametred species (like generalized parking functions)

        def restricted(self, min=0, max=Infinity):
            """
            The restriction of the species `F` to finite set of cardinal `n` with `min \leqslant n \leqstant max`.

            :param min, max: both are non-decreasing integer (or `\max` could be *None*)
            If `\max` is *None*, we assume it is `+ \infty`.
            :assert: min <= max

            That `F_n` is the species `F` restricted to finite set of cardinal `F`:

            MATH::

                F_n[U] = \begin{dcases*}
                    F[U] & if `|U| = n`,\\
                    \emptyset & otherwise.
                \end{dcases*}

            (section 1.3, _[BBL])

            This method defines

            MATH::

                F_{[\min, \max]} := F_\min + F_{\min+1} + \cdots + F_{\max}\,.
            """
            from sage.combinat.species2.operations.restriction import Restriction
            return Restriction(self, min=min, max=max)

        def product(self, G):
            """
            The *product* of `F (:= self)` and `G`.

            :param G: a species

            The species `F \cdot G` is the *product* of `F` and `G` defined as follows: an `(F\cdot G)`-structure on `U`
            is an ordered pair `s = (f,g)` where

             - `f` is an `F`-structure on a subset `U_1 \subset U`,
             - `g` is a `G`-structure on `U_2` the complementary of `U_1`.

            In other words,

            MATH::

                (F \cdot G)[U] = \sum_{U_1 \sqcup U_2 = U} F[U_1] \times G[U_2]\,.

            The transport along a bijection `sigma : U \to V` is carried out by setting

            MATH::

                (F \cdot G)[\sigma](s) = (F[\sigma_1](f), G[\sigma_2](g))

            where `sigma_i = \sigma_{|U_i}` is the restriction of `\sigma` on `U_i`, for each `(F + G)`-structure `s`.
            (section 1.3, _[BBL])
            """
            from sage.combinat.species2.operations.product import Prod
            return Prod(self, G)

        _mul_ = product

        def __pow__(self, k):
            """
            The *exponentiation* of `F` to the power `k`: `F^k`.
            :param k: a non-negative integer

            MATH::

                F^k = F \cdot F^{k-1}

            with `F^1 = F` and `F^0 = 1`.

            """
            assert(k >= 0), "`k` must be a non-negative integer."

            if k == 0:
                return Species().one()

            from sage.combinat.species2.operations.product import Prod
            return Prod(*([self]*k))

        def composite(self, G):
            """
            (Partitional) Composite of species.

            Let `F` and `G` be two species such that `G[\emptyset] = \emptyset`.
            The species `F \circ G`, also denoted `F(G)`, the (partitional) composite of `G` in `F`, is defined as follows:
            An `(F \circ G)`-structure on `U` is a triplet `s = (\pi, \varphi, \gamma)`, where

             - `\pi` is a partition of `U`,
             - `\varphi` is an `F`-structure on the set of classes of `\pi`,
             - `\gamma = (\gamma_p)_{p \in \pi}`, where for each class `p` of `\pi`, `\gamma_p` is a `G`-structure on `p`.

            In other words, for any finite set `U`, one has

            MATH::

                (F \circ G)[U] = \sum_{\pi \text{ partition of } U} F[\pi] \times \prod_{p \in \pi} G[p]\,,

            the (disjoint) sum being taken over the set of partitions `\pi` of `U` (*i.e.*, `\pi \in Par[U]`).

            The transport along a bijection `\sigma : U \to V` is carried out by setting, for any `(F \circ G)`-structure
            `s = (\pi, \varphi, (\gamma_p)_{p \in \pi})` on `U`,

            MATH::

                (F \circ G)[\sigma](s) = (\bar\pi, \bar\varphi, (\bar\gamma_{\bar{p}})_{\bar{p} \in \bar\pi})\,,

            where

             - `\bar\pi` is the partition of `V` obtained by the transport of `\pi` along `\sigma`,
             - for each `\bar{p} = \sigma(p) \in \bar\pi`, the structure `\bar\gamma_{\bar{p}}` is obtained from the structure
             `gamma_p` by `G`-transport along `\sigma_{|p}`,
             - the structure `\bar\varphi` is obtained from the transport `\varphi` by `F`-transport along the bijection
             `\bar\sigma` induced on `\pi` by `\sigma`.

            (section 1.4, _[BBL])
            """
            from sage.combinat.species2.operations.composite import Composite
            return Composite(self, G)

        # FIXME: Could we use the _call_ method such that we can to use F(G).
        ## the coercion system seems to forbid overload...

        def derivative(self):
            """
            Derivative of species

            The *derivative* of species `F`, noted `F'` is defined as follows:
            An `F'`-structure on `U` is an `F`-structure on `U^+ = U \cup \{\ast\}` where `\ast` is an element chosen
            outside of `U`. In other words, one sets `F'[U] = F[U^+]`.
            The transport along a bijection `\sigma : U \to V` is carried out by setting, by setting

            MATH::

                F'[\sigma](s) = F[\sigma^+](s)

            where `\sigma^+ : U \sqcup \{\ast\} \to V \sqcup \{\ast\}` is the canonical extension of `\sigma` obtained
            by setting `\sigma^+(u) = \sigma(u)` if `u \in U` and `\sigma^+(\ast) = \ast`.

            (section 1.4, _[BBL])
            """
            from sage.combinat.species2.operations.derivative import Derivative
            return Derivative(self)

        def pointing(self):
            """
            The point of species `F^\bullet` is defined as follows: An `F^\bullet`-structure on `U` is a pair
            `s = (f, u)`, where

             - `f` is an `F`-structure,
             - `u \in U` (a *distinguished element).

            The operations of pointing and derivative are related by the combinatorial equation

            MATH::

                F^\bullet = X \cdot F'\,.

            MATH::

                F^\bullet = F \times (X \cdot E)

            (section 2.1, _[BBL])
            """
            # TODO implement the pointing operator.
            return Species().singletons() * self.derivative()

        def is_pointing(self):
            """
            Test if `F` (*self*) is a pointing of species.

            This method
            """
            # Default implementation
            return False

        def cartesian_product(F, G):
            """
            The *cartesian product* of species

            The species `F \times G`, called *Cartesian product* of `F` and `G`, is defined as follows:
            An `(F \times G)`-structure on a finite set `U` is a pair `s = (f, g)`, where

             - `f` is an `F`-structure on `U`,
             - `g` is a `G`-structure on `U`.

            In other words, for all finite sets `U`, one has

            MATH::

                (F \times G)[U] = F[U] \times G[U]\,.

            The transport along a bijection `\sigma : U \to V` is carried out by setting

            MATH::

                (F \times G)[\sigma](s) = (F[\sigma](f), G[\sigma](g))\,,

            for any `(F \times G)`-structure `s = (f, g)` on `U`.

            (section 2.1, _[BBL])
            """
            from sage.combinat.species2.operations.cartesian_product import CartesianProduct
            return CartesianProduct(F, G)

        def functorial_composite(F, G):
            """
            The *functorial composite* of species

            The species `F \Box G` (also denoted `F[G]`) is the *functorial composite* of `F` and `G`. It is defined as
            follows: An `(F \Box G)`-structures on `U` is an `F`-structure placed on the set `G[U]` of all the
            `G`-structures on `U`.

            In other words, for any finite set `U`,

            MATH::

                (F \Box G)[U] = F[G[U]]\,.

            The transport along a bijection `\sigma : U \to V` is carried out by setting

            MATH::

                (F \Box G)[\sigma] = F[G[\sigma]]\,.

            (section 2.2, _[BBL])
            """
            from sage.combinat.species2.operations.functorial_composite import FunctorialComposite
            return FunctorialComposite(F, G)

        def _valuation_(self):
            """
            The valuation is the first degree `n` such that `F[n] \neq \emptyset`.
            """
            return self.generating_series()._valuation_()