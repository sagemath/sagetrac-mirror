# -*- coding: utf-8 -*-
"""
Sum of species.

References
----------

 _[BBL] Combinatorial species and tree-like structures,
 François Bergeron, Gilbert Labelle and Pierre Leroux,
 1998, Cambridge University Press

AUTHOR:

- Jean-Baptiste Priez (2015)
"""
#*****************************************************************************
#       Copyright (C) 2015 Jean-Baptiste Priez <jbp@kerios.fr>.
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.species import Species
from sage.combinat.species2 import SpeciesDesign
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.misc import uniq
from sage.structure.element_wrapper import ElementWrapper


class WrapperAddStructure(ElementWrapper):

    def __init__(self, parent, s, i):
        ElementWrapper.__init__(self, parent, s)
        self._indice_ = i

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self._indice_ == other._indice_ and self.value == other.value


class Add(SpeciesDesign):
    """
    Sum of species `+`.

    Let `F` and `G` be two species, we denote `H = F + G` species which is the sum of
    `F` and `G` defined by:

    MATH::

        (F + G)[U] := F[U] + G[U]

    the disjoint union of `F[U]` and `G[U]`.
    The transport along bijection `\sigma : U \to V` is carried out by setting, for any `(F + G)`-structures `s` on
    `[U]`,

    MATH::

        (F+G)[\sigma](s) := \begin{dcases*}
            F[\sigma](s) & if `s \in F[U]`,\\
            G[\sigma](s) & if `s \in G[U]`.
        \end{dcases*}

    The species `0` (..see :mod:`sage.combinat.species2.zero`) is the neutral element for addition:

    MATH::

        F + 0 = 0 + F = F\,,

    for any species `F`.

    The operation of addition is associative and commutative, up to isomorphism.
    (section 1.3, _[BLL])

    TESTS::

        sage: from sage.combinat.species2.zero import ZeroSpecies
        sage: Z = ZeroSpecies()
        sage: Z + Z
        `0`

        sage: from sage.combinat.species2.elements import ElementsSpecies
        sage: e = ElementsSpecies()
        sage: e + e
        2·`e`

        sage: from sage.combinat.species2.sets import SetsSpecies
        sage: E = SetsSpecies()
        sage: e + E
        `e` + `E`
        sage: e + E == E + e
        True
        sage: e + E == E + e + e
        False
        sage: e + E == E + e + Z
        True

        sage: from sage.combinat.species2.singletons import SingletonsSpecies
        sage: X = SingletonsSpecies()
        sage: (E + e) + X == E + (e + X)
        True

        sage: f = X + e
        sage: TestSuite(f).run()

    """
    # FIXME: tests about associativity and commutativity should be somewhere else.

    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **options):
        nargs = []
        ## Filter the species 0
        args = filter(lambda F: F != Species().zero(), args)
        if len(args) == 0:
            return Species.zero()
        elif len(args) == 1:
            return args[0]

        ## Associativity: `F + (G + H) --> F + G + H`
        for F in args:
            if isinstance(F, cls): nargs.extend(F._species_)
            else:                  nargs.append(F)

        ## Commutativity: we sort species such that `F + G = G + F`
        args = tuple(sorted(nargs))

        return super(Add, cls).__classcall__(cls, *args, **options)

    def __init__(self, *species):
        """

        :param F: a sequence of species

        """
        SpeciesDesign.__init__(self)
        self._species_ =  species

    def _repr_(self):
        def repr_species(F):
            if self._species_.count(F) == 1:
                return repr(F)
            return repr(self._species_.count(F)) + "·" + repr(F)

        return " + ".join(map(repr_species, uniq(self._species_)))

    def transport(self, sigma):

        def Fsigma(s):
            assert(isinstance(s, WrapperAddStructure))
            return self._element_constructor_(s.value.parent().transport(sigma)(s),
                                              s._indice_)

        return Fsigma

    def grading(self, s):
        assert(isinstance(s, WrapperAddStructure))
        return s.value.parent().grading(s.value)

    def _element_constructor_(self, *args, **options):
        return self.element_class(self, *args, **options)

    class Structures(SpeciesDesign.Structures):

        def cardinality(self):
            return reduce(lambda acc, F: acc + F.graded_component(self.grading()).cardinality(),
                          self.ambient()._species_, 0)

        def __iter__(self):
            for i, F in enumerate(self.ambient()._species_):
                for s in F.structures(self.finite_set()):
                    yield self._element_constructor_(s, i)

    Element = WrapperAddStructure


