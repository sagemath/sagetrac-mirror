# -*- coding: utf-8 -*-
"""
Species category

References
----------

 _[BBL] Combinatorial species and tree-like structures,
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
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
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
        return []

    @cached_method
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

    @cached_method
    def zero(self):
        """
        The species `0` defined by

        MATH::

            0[U] := \emptyset

        for any finite set `U`.
        """
        from sage.combinat.species2.zero import ZeroSpecies
        return ZeroSpecies()

    @cached_method
    def singleton(self):
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

        def isomorphism_types(self, n):
            """
            The *isomorphism types* of `F`-structures of order `n`.

            :param U: a finite set

            Let `\sim` be the equivalence relation on `F[n]` defined by `s \sim t` *iff* `s` and `t` have same
            isomorphism type.

            An *isomorphism type* of `F`-structures of order `n` is an equivalence class (modulo `\sim`) of
            `F`-structure on `[n]` (with `[n] := \{1, \cdots, n\}`).

            (see section 1.2 Type generating series, _[BBL])

            This method return the quotient `T(F[n]) := F[n]/\sim`.
            """
            # TODO give a generic implementation

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

        def exponential_generating_series(self):
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
            return self.cycle_index_series().exponential_generating_series()

        egs = generating_series = exponential_generating_series

        def isomorphism_type_generating_series(self):
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

        tgs = type_generating_series = isomorphism_type_generating_series

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
            """
            # TODO give a generic implementation

        cis = cycle_index_series

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
            """
            from sage.combinat.species2.operations.product import Prod
            return Prod(self, G)

        _mul_ = product

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