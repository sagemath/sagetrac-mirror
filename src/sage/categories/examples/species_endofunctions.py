# -*- coding: utf-8 -*-
"""
Category of species of structures

Definition (**Species** _[BLL])
-------------------------------

A *species of structures* is a rule `F` which produces

    i) for each finite set `U`, a finite set `F[U]`,
   ii) for each bijection `\sigma: U \to V`, a function `F[\sigma]: F[U] \to
       F[V]`.

The functions `F[\sigma]` should further satisfy the following functorial
properties:

    a) for all bijections `\sigma : U \to V` and `tau :V \to W`,

    MATH::

        F[\tau \circ \sigma] = F[\tau] \circ F[\sigma]\,,

    b) for the identity map `Id_U : U \to U`,

    MATH::

        F[Id_U] = Id_F[U].

An element `s \in F[U]` is called an *`F`-structure on `U`* (or even a
*structure of species `F` on `U`*).
The function `F [\sigma]` is called the *transport of `F`-structures along
`\sigma`*.

References:
-----------

.. [BLL] Combinatorial species and tree-like structures,
  François Bergeron, Gilbert Labelle, and Pierre Leroux

"""
from sage.misc.bijection import bijection
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.sets.set import Set
from sage.combinat.species import Species
from sage.combinat.structures import Structure


class Endofunction(Structure):

    def __init__(self, parent, dct, check=True):
        Structure.__init__(self, parent)
        self._dct_ = dct

        if check: self.check()

    def check(self):
        assert(Set(self._dct_.values()).issubset(Set(self._dct_.keys())))


    def underlying_set(self):
        """
        TESTS::

        """
        return Set(self._dct_.keys())

    def __call__(self, obj):
        return self._dct_[obj]

    def __eq__(self, other):
        if isinstance(other, Endofunction):
            return other._dct_ == self._dct_
        return False

    def _repr_(self):
        return repr(self._dct_)

    @lazy_class_attribute
    def _auto_parent_(cls):
        """
        TESTS::

            sage: Permutation(['a', 'b', 'c']).parent()
            Species of permutations

        """
        return Endofunctions()


class Endofunctions(Species):
    """
    Species of endofunctions

    TESTS::

    """
    @bijection
    def transport(self, bij):
        """
        TESTS::

            sage: from sage.categories.examples.species_permutations import \
                    Permutations, Permutation
            sage: P = Permutations()
            sage: sig = Permutation([(1,2),(3,)])
            sage: P.transport([1,3,2])(sig)
            [(1, 3), (2,)]

            sage: sig = Permutation([('a', 'c', 'b'), ('d',)])
            sage: P.transport({'a':1, 'b':2, 'c':3, 'd':4})(sig)
            [(1, 3, 2), (4,)]

        """
        def _transport_(f):
            assert(isinstance(f, Endofunction))
            return self._element_constructor_({bij(k):bij(f(k))
                                               for k in f._dct_.keys()})
        return _transport_

    def _repr_(self):
        """
        TESTS::

            sage: from sage.categories.examples.species_permutations import \
                    Permutations
            sage: Permutations()
            Species of permutations

        """
        return "Species of endofunctions"

    class Structures(Species.Structures):

        def __iter__(self):
            """
            TESTS::

            """
            def _rec_(U, V):
                if U.cardinality() == 0:
                    yield ()
                    return
                u = U[0]
                for end in _rec_(U.difference((u,)), V):
                    for v in V:
                        yield end + ((u,v),)

            U = self.underlying_set()
            for end in _rec_(U, U):
                yield self._element_constructor_({u: v for (u,v) in end})

    Element = Endofunction