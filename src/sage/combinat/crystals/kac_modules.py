"""
Crystals for Kac modules of the general-linear Lie superalgebra

REFERENCES:

.. [Kwon2012] Jae-Hoon Kwon. *Crystal bases of* `q`-*deformed Kac Modules
              over the Quantum Superalgebra* `U_q(\mathfrak{gl}(m|n))`.
              International Mathematics Research Notices. Vol. 2014, No. 2,
              pp. 512-550 (2012)
"""

#*****************************************************************************
#       Copyright (C) 2017 Travis Scrimshaw <tcscrims at gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.parent import Parent
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.all import ZZ

from sage.categories.regular_supercrystals import RegularSuperCrystals
from sage.combinat.crystals.tensor_product import CrystalOfTableaux
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.partition import _Partitions

from sage.combinat.crystals.letters import CrystalOfBKKLetters
from sage.combinat.crystals.tensor_product_element import CrystalOfBKKTableauxElement

class CrystalOfOddNegativeRoots(UniqueRepresentation, Parent):
    """
    Crystal of the set of odd negative roots.

    This is the crystal structure on the set of negative roots
    as given by Kwon. This is half of the Kac module crystals.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: from sage.combinat.crystals.kac_modules import CrystalOfOddNegativeRoots
            sage: S1 = crystals.OddNegativeRoots(['A', [2,1]])
            sage: S2 = crystals.OddNegativeRoots(CartanType(['A', [2,1]]))
            sage: S1 is S2
            True
        """
        return super(CrystalOfOddNegativeRoots, cls).__classcall__(cls, CartanType(cartan_type))

    def __init__(self, cartan_type):
        """
        Initialize ``self``.

        TESTS::

            sage: S = crystals.OddNegativeRoots(['A', [2,1]])
            sage: TestSuite(S).run()
        """
        self._cartan_type = cartan_type
        Parent.__init__(self, category=RegularSuperCrystals())
        self.module_generators = (self.element_class(self, frozenset()),)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: crystals.OddNegativeRoots(['A', [2,1]])
            Crystal of odd negative roots of type ['A', [2, 1]]
        """
        return "Crystal of odd negative roots of type {}".format(self._cartan_type)

    class Element(ElementWrapper):
        """
        An element of the crystal of odd negative roots.

        TESTS:

        Check that `e_i` and `f_i` are psuedo-inverses::

            sage: S = crystals.OddNegativeRoots(['A', [2,1]])
            sage: for x in S:
            ....:     for i in S.index_set():
            ....:         y = x.f(i)
            ....:         assert y is None or y.e(i) == x

        Check that we obtain the entire powerset of negative odd roots::

            sage: S = crystals.OddNegativeRoots(['A', [2,3]])
            sage: S.cardinality()
            4096
            sage: 2^len(S.weight_lattice_realization().positive_odd_roots())
            4096
        """
        def _repr_(self):
            return ('{'
                    + ", ".join("-e[{}]+e[{}]".format(*i)
                                for i in sorted(self.value))
                    + '}')

        def _latex_(self):
            return ('\{'
                    + ", ".join("-e_{{{}}}+e_{{{}}}".format(*i)
                                for i in sorted(self.value))
                    + '\}')

        def e(self, i):
            if i == 0:
                if (-1,1) not in self.value:
                    return None
                return type(self)(self.parent(), self.value.difference([(-1,1)]))

            count = 0
            act_val = None
            if i < 0:
                lst = sorted(self.value, key=lambda x: (x[1], -x[0]))
                for val in lst:
                    # We don't have to check val[1] because this is an odd root
                    if val[0] == i - 1:
                        if count == 0:
                            act_val = val
                        else:
                            count -= 1
                    elif val[0] == i:
                        count += 1
                if act_val is None:
                    return None
                ret = self.value.difference([act_val]).union([(i, act_val[1])])
                return type(self)(self.parent(), ret)

            # else i > 0
            lst = sorted(self.value, key=lambda x: (-x[0], -x[1]))
            for val in reversed(lst):
                # We don't have to check val[0] because this is an odd root
                if val[1] == i + 1:
                    if count == 0:
                        act_val = val
                    else:
                        count -= 1
                elif val[1] == i:
                    count += 1
            if act_val is None:
                return None
            ret = self.value.difference([act_val]).union([(act_val[0], i)])
            return type(self)(self.parent(), ret)

        def f(self, i):
            if i == 0:
                if (-1,1) in self.value:
                    return None
                return type(self)(self.parent(), self.value.union([(-1,1)]))

            count = 0
            act_val = None
            if i < 0:
                lst = sorted(self.value, key=lambda x: (x[1], -x[0]))
                for val in reversed(lst):
                    # We don't have to check val[1] because this is an odd root
                    if val[0] == i:
                        if count == 0:
                            act_val = val
                        else:
                            count -= 1
                    elif val[0] == i - 1:
                        count += 1
                if act_val is None:
                    return None
                ret = self.value.difference([act_val]).union([(i-1, act_val[1])])
                return type(self)(self.parent(), ret)

            # else i > 0
            lst = sorted(self.value, key=lambda x: (-x[0], -x[1]))
            for val in lst:
                # We don't have to check val[0] because this is an odd root
                if val[1] == i:
                    if count == 0:
                        act_val = val
                    else:
                        count -= 1
                elif val[1] == i + 1:
                    count += 1
            if act_val is None:
                return None
            ret = self.value.difference([act_val]).union([(act_val[0], i+1)])
            return type(self)(self.parent(), ret)

        def epsilon(self, i):
            r"""
            Return `\varepsilon_i` of ``self``.

            TESTS::

                sage: S = crystals.OddNegativeRoots(['A', [2,1]])
                sage: def count_e(x, i):
                ....:     ret = -1
                ....:     while x is not None:
                ....:         x = x.e(i)
                ....:         ret += 1
                ....:     return ret
                sage: for x in S:
                ....:     for i in S.index_set():
                ....:         assert x.epsilon(i) == count_e(x, i)
            """
            if i == 0:
                return ZZ.one() if (-1,1) in self.value else ZZ.zero()

            count = 0
            ret = 0
            if i < 0:
                lst = sorted(self.value, key=lambda x: (x[1], -x[0]))
                for val in lst:
                    # We don't have to check val[1] because this is an odd root
                    if val[0] == i - 1:
                        if count == 0:
                            ret += 1
                        else:
                            count -= 1
                    elif val[0] == i:
                        count += 1

            else: # i > 0
                lst = sorted(self.value, key=lambda x: (-x[0], -x[1]))
                for val in reversed(lst):
                    # We don't have to check val[0] because this is an odd root
                    if val[1] == i + 1:
                        if count == 0:
                            ret += 1
                        else:
                            count -= 1
                    elif val[1] == i:
                        count += 1
            return ret

        def phi(self, i):
            r"""
            Return `\varphi_i` of ``self``.

            TESTS::

                sage: S = crystals.OddNegativeRoots(['A', [2,1]])
                sage: def count_f(x, i):
                ....:     ret = -1
                ....:     while x is not None:
                ....:         x = x.f(i)
                ....:         ret += 1
                ....:     return ret
                sage: for x in S:
                ....:     for i in S.index_set():
                ....:         assert x.phi(i) == count_f(x, i)
            """
            if i == 0:
                return ZZ.zero() if (-1,1) in self.value else ZZ.one()

            count = 0
            ret = 0
            if i < 0:
                lst = sorted(self.value, key=lambda x: (x[1], -x[0]))
                for val in reversed(lst):
                    # We don't have to check val[1] because this is an odd root
                    if val[0] == i:
                        if count == 0:
                            ret += 1
                        else:
                            count -= 1
                    elif val[0] == i - 1:
                        count += 1

            else: # i > 0
                lst = sorted(self.value, key=lambda x: (-x[0], -x[1]))
                for val in lst:
                    # We don't have to check val[0] because this is an odd root
                    if val[1] == i:
                        if count == 0:
                            ret += 1
                        else:
                            count -= 1
                    elif val[1] == i + 1:
                        count += 1
            return ret

        def weight(self):
            """
            Return the weight of ``self``.

            TESTS::

                sage: S = crystals.OddNegativeRoots(['A', [2,1]])
                sage: al = S.weight_lattice_realization().simple_roots()
                sage: for x in S:
                ....:     for i in S.index_set():
                ....:         y = x.f(i)
                ....:         assert y is None or x.weight() - al[i] == y.weight()
            """
            WLR = self.parent().weight_lattice_realization()
            e = WLR.basis()
            return WLR.sum(-e[i]+e[j] for (i,j) in self.value)

class CrystalOfKacModule(UniqueRepresentation, Parent):
    """
    Crystal of a Kac module.

    .. NOTE::

        Our notation differs slightly from Kwon: our last tableau is
        transposed from that of Kwon.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, la, mu):
        """
        Normalize input to ensure a unique representation.

        TESTS::

            sage: K1 = crystals.KacModule(['A', [2,1]], [2,1], [1])
            sage: K2 = crystals.KacModule(CartanType(['A', [2,1]]), (2,1), (1,))
            sage: K1 is K2
            True
        """
        cartan_type = CartanType(cartan_type)
        la = _Partitions(la)
        mu = _Partitions(mu)
        return super(CrystalOfKacModule, cls).__classcall__(cls, cartan_type, la, mu)

    def __init__(self, cartan_type, la, mu):
        """
        Initialize ``self``.

        TESTS::

            sage: K = crystals.KacModule(['A', [2,1]], [2,1], [1])
            sage: TestSuite(K).run()
        """
        self._cartan_type = cartan_type
        self._la = la
        self._mu = mu
        Parent.__init__(self, category=RegularSuperCrystals())
        self._S = CrystalOfOddNegativeRoots(self._cartan_type)
        self._dual = CrystalOfTableaux(['A', self._cartan_type.m], shape=la)
        self._reg = CrystalOfTableaux(['A', self._cartan_type.n], shape=mu)
        data = (self._S.module_generators[0],
                self._dual.module_generators[0],
                self._reg.module_generators[0])
        self.module_generators = (self.element_class(self, data),)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Crystal of Kac module K({}, {}) of type {}".format(
                    self._cartan_type, self._la, self._mu)

    class Element(ElementWrapper):
        """
        An element of a Kac module crystal.

        TESTS:

        Check that `e_i` and `f_i` are psuedo-inverses::

            sage: K = crystals.KacModule(['A', [2,1]], [2,1], [1])
            sage: for x in K:
            ....:     for i in K.index_set():
            ....:         y = x.f(i)
            ....:         assert y is None or y.e(i) == x
        """
        def _repr_(self):
            return repr((self.value[0], to_dual_tableau(self.value[1]), self.value[2]))

        def _latex_(self):
            from sage.misc.latex import latex
            return " \otimes ".join([latex(self.value[0]),
                                     latex_dual(self.value[1]),
                                     latex(self.value[2])])

        def e(self, i):
            if i == 0:
                x = self.value[0].e(i)
                if x is None:
                    return None
                return type(self)(self.parent(), (x, self.value[1], self.value[2]))
            if i > 0:
                if self.value[0].epsilon(i) > self.value[2].phi(i):
                    x = self.value[0].e(i)
                    if x is None:
                        return None
                    return type(self)(self.parent(), (x, self.value[1], self.value[2]))
                else:
                    x = self.value[2].e(i)
                    if x is None:
                        return None
                    return type(self)(self.parent(), (self.value[0], self.value[1], x))
            # else i < 0
            M = self.parent()._cartan_type.m + 1
            if self.value[0].phi(i) < self.value[1].epsilon(M+i):
                x = self.value[1].e(M+i)
                if x is None:
                    return None
                return type(self)(self.parent(), (self.value[0], x, self.value[2]))
            else:
                x = self.value[0].e(i)
                if x is None:
                    return None
                return type(self)(self.parent(), (x, self.value[1], self.value[2]))

        def f(self, i):
            if i == 0:
                x = self.value[0].f(i)
                if x is None:
                    return None
                return type(self)(self.parent(), (x, self.value[1], self.value[2]))
            if i > 0:
                if self.value[0].epsilon(i) < self.value[2].phi(i):
                    x = self.value[2].f(i)
                    if x is None:
                        return None
                    return type(self)(self.parent(), (self.value[0], self.value[1], x))
                else:
                    x = self.value[0].f(i)
                    if x is None:
                        return None
                    return type(self)(self.parent(), (x, self.value[1], self.value[2]))
            # else i < 0
            M = self.parent()._cartan_type.m + 1
            if self.value[0].phi(i) > self.value[1].epsilon(M+i):
                x = self.value[0].f(i)
                if x is None:
                    return None
                return type(self)(self.parent(), (x, self.value[1], self.value[2]))
            else:
                x = self.value[1].f(M+i)
                if x is None:
                    return None
                return type(self)(self.parent(), (self.value[0], x, self.value[2]))

        def weight(self):
            e = self.parent().weight_lattice_realization().basis()
            M = self.parent()._cartan_type.m + 1
            wt = self.value[0].weight()
            wt += sum(c*e[i-M] for i,c in self.value[1].weight())
            wt += sum(c*e[i+1] for i,c in self.value[2].weight())
            return wt

#####################################################################
## Helper functions

def to_dual_tableau(elt):
    M = elt.parent().cartan_type().rank() + 2
    if not elt:
        return Tableau([])
    tab = [ [elt[0].value-M] ]
    for i in range(1, len(elt)):
        if elt[i-1] < elt[i] or (elt[i-1].value != 0 and elt[i-1] == elt[i]):
            tab.append([elt[i].value-M])
        else:
            tab[len(tab)-1].append(elt[i].value-M)
    for x in tab:
        x.reverse()
    from sage.combinat.tableau import Tableau
    return Tableau(tab).conjugate()

def latex_dual(elt):
    M = elt.parent().cartan_type().rank() + 2
    from sage.combinat.output import tex_from_array
    # Modified version of to_tableau() to have the entries be letters
    #   rather than their values
    if not elt:
        return "{\\emptyset}"

    tab = [ ["\\overline{{{}}}".format(M-elt[0].value)] ]
    for i in range(1, len(elt)):
        if elt[i-1] < elt[i] or (elt[i-1].value != 0 and elt[i-1] == elt[i]):
            tab.append(["\\overline{{{}}}".format(M-elt[i].value)])
        else:
            l = len(tab)-1
            tab[l].append("\\overline{{{}}}".format(M-elt[i].value))
    for x in tab:
        x.reverse()
    from sage.combinat.tableau import Tableau
    T = Tableau(tab).conjugate()
    from sage.combinat.output import tex_from_array
    return tex_from_array([list(row) for row in T])

