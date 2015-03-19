# -*- coding: utf-8 -*-
"""
The restriction of species.

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
from sage.categories.species import Species
from sage.combinat.species2 import SpeciesDesign
from sage.combinat.species2.singletons import SingletonsSpecies
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer


class Restriction(SpeciesDesign):
    """
    The species `F_I` is the restriction of `F` to sets of cardinality `n` with `n \in I`.
    More precisely,

    MATH::

        F_I [U] = \begin{dcases*}
            F[U] & if `|U| \in I`,\\
            \emptyset & otherwise,
        \end{dcases*}

    and set

    MATH::

        F_I[\sigma] = F[\sigma]\,,

    for any finite sets `U` and `V` and any bijection `\sigma : U \to V`.

    (section 1.2, _[BBL])

    EXAMPLES::

        sage: from sage.combinat.species2.sets import SetsSpecies
        sage: E = SetsSpecies()
        sage: E3 = E.restricted(min=3)
        sage: [E3.graded_component(n).list() for n in range(10)]
        [[],
         [],
         [],
         [{1, 2, 3}],
         [{1, 2, 3, 4}],
         [{1, 2, 3, 4, 5}],
         [{1, 2, 3, 4, 5, 6}],
         [{1, 2, 3, 4, 5, 6, 7}],
         [{1, 2, 3, 4, 5, 6, 7, 8}],
         [{1, 2, 3, 4, 5, 6, 7, 8, 9}]]

    TESTS::

        sage: from sage.combinat.species2.sets import SetsSpecies
        sage: E = SetsSpecies()
        sage: E3 = E.restricted(min=3)
        sage: TestSuite(E3).run()

    """

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        F = args[0]
        m = opts["min"]
        M = opts["max"]
        # TODO: to do or not to do, that is the question
        if m == M:
            if m == 0 and F.graded_component(m).cardinality() == 1:
                return Species().one()
            elif m == 1 and F.graded_component(m).cardinality() == 1:
                return SingletonsSpecies()
        elif m == 0 and M == Infinity:
            return F
        else:
            return super(Restriction, cls).__classcall__(cls, F, **opts)

    def __init__(self, F, min=Integer(0), max=Infinity):
        SpeciesDesign.__init__(self)
        self._F_ = F
        self._min_ = min
        self._max_ = max
        self.Element = F.Element

    def transport(self, sigma):

        def Fsigma(s):
            t = self._F_.transport(sigma)(s)
            t._set_parent(self)
            return t

        return Fsigma

    def grading(self, s):
        s._set_parent(self._F_)
        n = self._F_.grading(s)
        s._set_parent(self)
        return n

    def _repr_(self):
        if self._max_ == Infinity:
            if self._min_ == 1:
                s = "+"
            else:
                s = "≥" + repr(self._min_)
        elif self._min_ == 0:
            s = "≤" + repr(self._max_)
        else:
            s = "[%d, %d]"%(self._min_, self._max_)
        return repr(self._F_) + "_{" + s + "}"


    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            min = self.ambient()._min_
            max = self.ambient()._max_
            if self.grading() < min or self.grading() > max:
                return 0
            return self.ambient()._F_.graded_component(self.grading()).cardinality()

        def __iter__(self):
            min = self.ambient()._min_
            max = self.ambient()._max_
            if self.grading() < min or self.grading() > max:
                return

            F = self.ambient()._F_
            for s in F.structures(self.finite_set()):
                s._set_parent(self.ambient())
                yield s
