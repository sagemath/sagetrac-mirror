"""
Tensor Products of Crystals

Main entry points:

- :class:`~sage.combinat.crystals.tensor_product.TensorProductOfCrystals`
- :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableaux`

AUTHORS:

- Anne Schilling, Nicolas Thiery (2007): Initial version
- Ben Salisbury, Travis Scrimshaw (2013): Refactored tensor products to handle
  non-regular crystals and created new subclass to take advantage of
  the regularity
"""
#*****************************************************************************
#       Copyright (C) 2007 Anne Schilling <anne at math.ucdavis.edu>
#                          Nicolas Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

import operator
from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.global_options import GlobalOptions
from sage.categories.category import Category
from sage.categories.cartesian_product import cartesian_product
from sage.categories.classical_crystals import ClassicalCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.sets_cat import Sets
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.partition import Partition
from .letters import CrystalOfLetters
from .spins import CrystalOfSpins, CrystalOfSpinsMinus, CrystalOfSpinsPlus
from sage.combinat.crystals.tensor_product_element import (TensorProductOfCrystalsElement,
        TensorProductOfRegularCrystalsElement, CrystalOfTableauxElement)
from sage.misc.flatten import flatten
from sage.structure.element import get_coercion_model

##############################################################################
# Until trunc gets implemented in sage.function.other

from sage.functions.other import floor, ceil
def trunc(i):
    """
    Truncates to the integer closer to zero

    EXAMPLES::

        sage: from sage.combinat.crystals.tensor_product import trunc
        sage: trunc(-3/2), trunc(-1), trunc(-1/2), trunc(0), trunc(1/2), trunc(1), trunc(3/2)
        (-1, -1, 0, 0, 0, 1, 1)
        sage: isinstance(trunc(3/2), Integer)
        True
    """
    if i>= 0:
        return floor(i)
    else:
        return ceil(i)

##############################################################################
# Support classes
##############################################################################

class CrystalOfWords(UniqueRepresentation, Parent):
    """
    Auxiliary class to provide a call method to create tensor product elements.
    This class is shared with several tensor product classes and is also used
    in :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableaux`
    to allow tableaux of different tensor product structures in
    column-reading (and hence different shapes) to be considered elements
    in the same crystal.
    """
    def _element_constructor_(self, *crystalElements):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(1,1)
            [1, 1]
            sage: _.parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: T(C(2), C(1), C(1))
            [2, 1, 1]
        """
        return self.element_class(self, list(crystalElements))

    class Element(TensorProductOfCrystalsElement):
        pass

class TensorProductOfCrystals(CrystalOfWords):
    r"""
    Tensor product of crystals.

    Given two crystals `B` and `B'` of the same Cartan type,
    one can form the tensor product `B \otimes B^{\prime}`. As a set
    `B \otimes B^{\prime}` is the Cartesian product
    `B \times B^{\prime}`. The crystal operators `f_i` and
    `e_i` act on `b \otimes b^{\prime} \in B \otimes B^{\prime}` as
    follows:

    .. MATH::

        f_i(b \otimes b^{\prime}) = \begin{cases}
        f_i(b) \otimes b^{\prime} & \text{if } \varepsilon_i(b) \geq
        \varphi_i(b^{\prime}) \\
        b \otimes f_i(b^{\prime}) & \text{otherwise}
        \end{cases}

    and

    .. MATH::

        e_i(b \otimes b') = \begin{cases}
        e_i(b) \otimes b' & \text{if } \varepsilon_i(b) >
        \varphi_i(b') \\ b \otimes e_i(b') & \text{otherwise.}
        \end{cases}

    We also define:

    .. MATH::

        \begin{aligned}
        \varphi_i(b \otimes b') & = \max\left( \varphi_i(b),
        \varphi_i(b') + \langle \alpha_i^{\vee}, \mathrm{wt}(b) \rangle
        \right),
        \\ \varepsilon_i(b \otimes b') & = \max\left( \varepsilon_i(b'),
        \varepsilon_i(b) - \langle \alpha_i^{\vee}, \mathrm{wt}(b') \rangle
        \right).
        \end{aligned}

    .. NOTE::

        This is the opposite of Kashiwara's convention for tensor
        products of crystals.

    Since tensor products are associative `(\mathcal{B} \otimes \mathcal{C})
    \otimes \mathcal{D} \cong \mathcal{B} \otimes (\mathcal{C} \otimes
    \mathcal{D})` via the natural isomorphism `(b \otimes c) \otimes d \mapsto
    b \otimes (c \otimes d)`, we can generalizing this to arbitrary tensor
    products. Thus consider `B_N \otimes \cdots \otimes B_1`, where each
    `B_k` is an abstract crystal. The underlying set of the tensor product is
    `B_N \times \cdots \times B_1`, while the crystal structure is given
    as follows. Let `I` be the index set, and fix some `i \in I` and `b_N
    \otimes \cdots \otimes b_1 \in B_N \otimes \cdots \otimes B_1`. Define

    .. MATH::

        a_i(k) := \varepsilon_i(b_k) - \sum_{j=1}^{k-1} \langle
        \alpha_i^{\vee}, \mathrm{wt}(b_j) \rangle.

    Then

    .. MATH::

        \begin{aligned}
        \mathrm{wt}(b_N \otimes \cdots \otimes b_1) &=
        \mathrm{wt}(b_N) + \cdots + \mathrm{wt}(b_1),
        \\ \varepsilon_i(b_N \otimes \cdots \otimes b_1) &= \max_{1 \leq k
        \leq n}\left( \sum_{j=1}^k \varepsilon_i(b_j) - \sum_{j=1}^{k-1}
        \varphi_i(b_j) \right)
        \\ & = \max_{1 \leq k \leq N}\bigl( a_i(k) \bigr),
        \\ \varphi_i(b_N \otimes \cdots \otimes b_1) &= \max_{1 \leq k \leq N}
        \left( \varphi_i(b_N) + \sum_{j=k}^{N-1} \big( \varphi_i(b_j) -
        \varepsilon_i(b_{j+1}) \big) \right)
        \\ & = \max_{1 \leq k \leq N}\bigl( \lambda_i + a_i(k) \bigr)
        \end{aligned}

    where `\lambda_i = \langle \alpha_i^{\vee}, \mathrm{wt}(b_N \otimes \cdots
    \otimes b_1) \rangle`. Then for `k = 1, \ldots, N` the action of the
    Kashiwara operators is determined as follows.

    - If `a_i(k) > a_i(j)` for `1 \leq j < k` and `a_i(k) \geq a_i(j)`
      for `k < j \leq N`:

      .. MATH::

          e_i(b_N \otimes \cdots \otimes b_1) = b_N \otimes \cdots \otimes
          e_i b_k \otimes \cdots \otimes b_1.

    - If `a_i(k) \geq a_i(j)` for `1 \leq j < k` and `a_i(k) > a_i(j)`
      for `k < j \leq N`:

      .. MATH::

          f_i(b_N \otimes \cdots \otimes b_1) = b_N \otimes \cdots \otimes
          f_i b_k \otimes \cdots \otimes b_1.

    Note that this is just recursively applying the definition of the tensor
    product on two crystals. Recall that `\langle \alpha_i^{\vee},
    \mathrm{wt}(b_j) \rangle = \varphi_i(b_j) - \varepsilon_i(b_j)` by the
    definition of a crystal.

    .. RUBRIC:: Regular crystals

    Now if all crystals `B_k` are regular crystals, all `\varepsilon_i` and
    `\varphi_i` are non-negative and we can
    define tensor product by the *signature rule*. We start by writing a word
    in `+` and `-` as follows:

    .. MATH::

        \underbrace{- \cdots -}_{\varphi_i(b_N) \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i(b_N) \text{ times}}
        \quad \cdots \quad
        \underbrace{- \cdots -}_{\varphi_i(b_1) \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i(b_1) \text{ times}},

    and then canceling ordered pairs of `+-` until the word is in the reduced
    form:

    .. MATH::

        \underbrace{- \cdots -}_{\varphi_i \text{ times}} \quad
        \underbrace{+ \cdots +}_{\varepsilon_i \text{ times}}.

    Here `e_i` acts on the factor corresponding to the leftmost `+` and `f_i`
    on the factor corresponding to the rightmost `-`. If there is no `+` or
    `-` respectively, then the result is `0` (``None``).

    EXAMPLES:

    We construct the type `A_2`-crystal generated by `2 \otimes 1 \otimes 1`::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])

    It has `8` elements::

        sage: T.list()
        [[2, 1, 1], [2, 1, 2], [2, 1, 3], [3, 1, 3],
         [3, 2, 3], [3, 1, 1], [3, 1, 2], [3, 2, 2]]

    One can also check the Cartan type of the crystal::

        sage: T.cartan_type()
        ['A', 2]

    Other examples include crystals of tableaux (which internally are
    represented as tensor products obtained by reading the tableaux
    columnwise)::

        sage: C = crystals.Tableaux(['A',3], shape=[1,1,0])
        sage: D = crystals.Tableaux(['A',3], shape=[1,0,0])
        sage: T = crystals.TensorProduct(C,D, generators=[[C(rows=[[1], [2]]), D(rows=[[1]])], [C(rows=[[2], [3]]), D(rows=[[1]])]])
        sage: T.cardinality()
        24
        sage: TestSuite(T).run()
        sage: T.module_generators
        ([[[1], [2]], [[1]]], [[[2], [3]], [[1]]])
        sage: [x.weight() for x in T.module_generators]
        [(2, 1, 0, 0), (1, 1, 1, 0)]

    If no module generators are specified, we obtain the full tensor
    product::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C)
        sage: T.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
        sage: T.cardinality()
        9

    For a tensor product of crystals without module generators, the
    default implementation of ``module_generators`` contains all elements
    in the tensor product of the crystals. If there is a subset of
    elements in the tensor product that still generates the crystal,
    this needs to be implemented for the specific crystal separately::

        sage: T.module_generators.list()
        [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    For classical highest weight crystals, it is also possible to list
    all highest weight elements::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)],[C(1),C(2),C(1)]])
        sage: T.highest_weight_vectors()
        ([2, 1, 1], [1, 2, 1])

    Examples with non-regular and infinite crystals (these did not work
    before :trac:`14402`)::

        sage: B = crystals.infinity.Tableaux(['D',10])
        sage: T = crystals.TensorProduct(B,B)
        sage: T
        Full tensor product of the crystals
        [The infinity crystal of tableaux of type ['D', 10],
         The infinity crystal of tableaux of type ['D', 10]]

        sage: B = crystals.infinity.GeneralizedYoungWalls(15)
        sage: T = crystals.TensorProduct(B,B,B)
        sage: T
        Full tensor product of the crystals
        [Crystal of generalized Young walls of type ['A', 15, 1],
         Crystal of generalized Young walls of type ['A', 15, 1],
         Crystal of generalized Young walls of type ['A', 15, 1]]

        sage: La = RootSystem(['A',2,1]).weight_lattice(extended=True).fundamental_weights()
        sage: B = crystals.GeneralizedYoungWalls(2,La[0]+La[1])
        sage: C = crystals.GeneralizedYoungWalls(2,2*La[2])
        sage: D = crystals.GeneralizedYoungWalls(2,3*La[0]+La[2])
        sage: T = crystals.TensorProduct(B,C,D)
        sage: T
        Full tensor product of the crystals
        [Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight Lambda[0] + Lambda[1],
         Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight 2*Lambda[2],
         Highest weight crystal of generalized Young walls of Cartan type ['A', 2, 1] and highest weight 3*Lambda[0] + Lambda[2]]

    There is also a global option for setting the convention (by default Sage
    uses anti-Kashiwara)::

        sage: C = crystals.Letters(['A',2])
        sage: T = crystals.TensorProduct(C,C)
        sage: elt = T(C(1), C(2)); elt
        [1, 2]
        sage: crystals.TensorProduct.options.convention = "Kashiwara"
        sage: elt
        [2, 1]
        sage: crystals.TensorProduct.options._reset()
    """
    @staticmethod
    def __classcall_private__(cls, *crystals, **options):
        """
        Create the correct parent object.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C, C)
            sage: T2 = crystals.TensorProduct(C, C, cartan_type=['A',2])
            sage: T is T2
            True
            sage: T.category()
            Category of tensor products of classical crystals

            sage: T3 = crystals.TensorProduct(C, C, C)
            sage: T3p = crystals.TensorProduct(T, C)
            sage: T3 is T3p
            True
            sage: B1 = crystals.TensorProduct(T, C)
            sage: B2 = crystals.TensorProduct(C, T)
            sage: B3 = crystals.TensorProduct(C, C, C)
            sage: B1 is B2 and B2 is B3
            True

            sage: B = crystals.infinity.Tableaux(['A',2])
            sage: T = crystals.TensorProduct(B, B)
            sage: T.category()
            Category of infinite tensor products of highest weight crystals

        TESTS:

        Check that mismatched Cartan types raise an error::

            sage: A2 = crystals.Letters(['A', 2])
            sage: A3 = crystals.Letters(['A', 3])
            sage: crystals.TensorProduct(A2, A3)
            Traceback (most recent call last):
            ...
            ValueError: all crystals must be of the same Cartan type
        """
        crystals = tuple(crystals)
        if "cartan_type" in options:
            cartan_type = CartanType(options.pop("cartan_type"))
        else:
            if not crystals:
                raise ValueError("you need to specify the Cartan type if the tensor product list is empty")
            else:
                cartan_type = crystals[0].cartan_type()

        if any(c.cartan_type() != cartan_type for c in crystals):
            raise ValueError("all crystals must be of the same Cartan type")

        if "generators" in options:
            generators = tuple(tuple(x) if isinstance(x, list) else x for x in options["generators"])

            if all(c in RegularCrystals() for c in crystals):
                return TensorProductOfRegularCrystalsWithGenerators(crystals, generators, cartan_type)
            return TensorProductOfCrystalsWithGenerators(crystals, generators, cartan_type)

        # Flatten out tensor products
        tp = sum([B.crystals if isinstance(B, FullTensorProductOfCrystals) else (B,)
                  for B in crystals], ())

        if all(c in RegularCrystals() for c in crystals):
            return FullTensorProductOfRegularCrystals(tp, cartan_type=cartan_type)
        return FullTensorProductOfCrystals(tp, cartan_type=cartan_type)

    # add options to class
    options=GlobalOptions('TensorProductOfCrystals', 
        module='sage.combinat.crystals',
        doc=r"""
        Sets the global options for tensor products of crystals. The default is to
        use the anti-Kashiwara convention.

        There are two conventions for how `e_i` and `f_i` act on tensor products,
        and the difference between the two is the order of the tensor factors
        are reversed. This affects both the input and output. See the example
        below.
        """,
        end_doc=r"""

        .. NOTE::

            Changing the ``convention`` also changes how the input is handled.

        .. WARNING::

            Internally, the crystals are always stored using the anti-Kashiwara
            convention.

        If no parameters are set, then the function returns a copy of the
        options dictionary.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: elt = T(C(1), C(2)); elt
            [1, 2]
            sage: crystals.TensorProduct.options.convention = "Kashiwara"
            sage: elt
            [2, 1]
            sage: T(C(1), C(2)) == elt
            False
            sage: T(C(2), C(1)) == elt
            True
            sage: crystals.TensorProduct.options._reset()
        """,
        convention=dict(default="antiKashiwara",
                        description='Sets the convention used for displaying/inputting tensor product of crystals',
                        values=dict(antiKashiwara='use the anti-Kashiwara convention',
                                    Kashiwara='use the Kashiwara convention'),
                            alias=dict(anti="antiKashiwara", opposite="antiKashiwara"),
                            case_sensitive=False)
    )

    def _element_constructor_(self, *crystalElements):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T(1,1)
            [1, 1]
            sage: _.parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: T(C(2), C(1), C(1))
            [2, 1, 1]
        """
        if self.options.convention == "Kashiwara":
            crystalElements = reversed(crystalElements)
        return self.element_class(self, list(crystalElements))

class TensorProductOfCrystalsWithGenerators(TensorProductOfCrystals):
    """
    Tensor product of crystals with a generating set.

    .. TODO::

        Deprecate this class in favor of using
        :meth:`~sage.categories.crystals.Crystals.ParentMethods.subcrystal`.
    """
    def __init__(self, crystals, generators, cartan_type):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C,C,generators=[[C(2),C(1),C(1)]])
            sage: TestSuite(T).run()
        """
        assert isinstance(crystals, tuple)
        assert isinstance(generators, tuple)
        category = Category.meet([crystal.category() for crystal in crystals])
        Parent.__init__(self, category = category)
        self.crystals = crystals
        self._cartan_type = cartan_type
        self.module_generators = tuple([self(*x) for x in generators])

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: crystals.TensorProduct(C,C,generators=[[C(2),C(1)]])
            The tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        if self.options.convention == "Kashiwara":
            st = repr(list(reversed(self.crystals)))
        else:
            st = repr(list(self.crystals))
        return "The tensor product of the crystals {}".format(st)

class FullTensorProductOfCrystals(TensorProductOfCrystals):
    """
    Full tensor product of crystals.

    .. TODO::

        Merge this into :class:`TensorProductOfCrystals`.
    """
    def __init__(self, crystals, **options):
        """
        TESTS::

            sage: from sage.combinat.crystals.tensor_product import FullTensorProductOfCrystals
            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: isinstance(T, FullTensorProductOfCrystals)
            True
            sage: TestSuite(T).run()
        """
        category = Category.meet([crystal.category() for crystal in crystals])
        category = category.TensorProducts()
        if any(c in Sets().Infinite() for c in crystals):
            category = category.Infinite()
        Parent.__init__(self, category=category)
        self.crystals = crystals
        if 'cartan_type' in options:
            self._cartan_type = CartanType(options['cartan_type'])
        else:
            if not crystals:
                raise ValueError("you need to specify the Cartan type if the tensor product list is empty")
            else:
                self._cartan_type = crystals[0].cartan_type()
        self.cartesian_product = cartesian_product(self.crystals)
        self.module_generators = self

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: crystals.TensorProduct(C,C)
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        if self.options.convention == "Kashiwara":
            st = repr(list(reversed(self.crystals)))
        else:
            st = repr(list(self.crystals))
        return "Full tensor product of the crystals {}".format(st)

    # TODO: __iter__ and cardinality should be inherited from EnumeratedSets().CartesianProducts()
    def __iter__(self):
        """
        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: list(T)
            [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
            sage: _[0].parent()
            Full tensor product of the crystals [The crystal of letters for type ['A', 2], The crystal of letters for type ['A', 2]]
        """
        for x in self.cartesian_product:
            yield self(*x)

#    list = CombinatorialClass._CombinatorialClass__list_from_iterator

    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::

            sage: C = crystals.Letters(['A',2])
            sage: T = crystals.TensorProduct(C,C)
            sage: T.cardinality()
            9
        """
        return self.cartesian_product.cardinality()

    @cached_method
    def weight_lattice_realization(self):
        r"""
        Return the weight lattice realization used to express weights.

        The weight lattice realization is the common parent which all
        weight lattice realizations of the crystals of ``self`` coerce
        into.

        EXAMPLES::

            sage: B = crystals.elementary.B(['A',4], 2)
            sage: B.weight_lattice_realization()
            Root lattice of the Root system of type ['A', 4]
            sage: T = crystals.infinity.Tableaux(['A',4])
            sage: T.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 4]
            sage: TP = crystals.TensorProduct(B, T)
            sage: TP.weight_lattice_realization()
            Ambient space of the Root system of type ['A', 4]
        """
        cm = get_coercion_model()
        return cm.common_parent(*[crystal.weight_lattice_realization()
                                  for crystal in self.crystals])

class FullTensorProductOfRegularCrystals(FullTensorProductOfCrystals):
    """
    Full tensor product of regular crystals.
    """
    class Element(TensorProductOfRegularCrystalsElement):
        pass

class TensorProductOfRegularCrystalsWithGenerators(TensorProductOfCrystalsWithGenerators):
    """
    Tensor product of regular crystals with a generating set.
    """
    class Element(TensorProductOfRegularCrystalsElement):
        pass

#########################################################
## Crystal of tableaux

class CrystalOfTableaux(CrystalOfWords):
    r"""
    A class for crystals of tableaux with integer valued shapes

    INPUT:

    - ``cartan_type`` -- a Cartan type
    - ``shape`` -- a partition of length at most ``cartan_type.rank()``
    - ``shapes`` -- a list of such partitions

    This constructs a classical crystal with the given Cartan type and
    highest weight(s) corresponding to the given shape(s).

    If the type is `D_r`, the shape is permitted to have a negative
    value in the `r`-th position. Thus if the shape equals `[s_1,\ldots,s_r]`,
    then `s_r` may be negative but in any case `s_1 \geq \cdots \geq s_{r-1}
    \geq |s_r|`. This crystal is related to that of shape
    `[s_1,\ldots,|s_r|]` by the outer automorphism of `SO(2r)`.

    If the type is `D_r` or `B_r`, the shape is permitted to be of
    length `r` with all parts of half integer value. This corresponds
    to having one spin column at the beginning of the tableau. If
    several shapes are provided, they currently should all or none
    have this property.

    Crystals of tableaux are constructed using an embedding into
    tensor products following Kashiwara and Nakashima [KN94]_. Sage's tensor
    product rule for crystals differs from that of Kashiwara and Nakashima
    by reversing the order of the tensor factors. Sage produces the same
    crystals of tableaux as Kashiwara and Nakashima. With Sage's convention,
    the tensor product of crystals is the same as the monoid operation on
    tableaux and hence the plactic monoid.

    .. SEEALSO::

        :mod:`sage.combinat.crystals.crystals` for general help on
        crystals, and in particular plotting and `\LaTeX` output.

    EXAMPLES:

    We create the crystal of tableaux for type `A_2`, with
    highest weight given by the partition `[2,1,1]`::

        sage: T = crystals.Tableaux(['A',3], shape = [2,1,1])

    Here is the list of its elements::

        sage: T.list()
        [[[1, 1], [2], [3]], [[1, 2], [2], [3]], [[1, 3], [2], [3]],
         [[1, 4], [2], [3]], [[1, 4], [2], [4]], [[1, 4], [3], [4]],
         [[2, 4], [3], [4]], [[1, 1], [2], [4]], [[1, 2], [2], [4]],
         [[1, 3], [2], [4]], [[1, 3], [3], [4]], [[2, 3], [3], [4]],
         [[1, 1], [3], [4]], [[1, 2], [3], [4]], [[2, 2], [3], [4]]]

    Internally, a tableau of a given Cartan type is represented as a
    tensor product of letters of the same type. The order in which the
    tensor factors appear is by reading the columns of the tableaux
    left to right, top to bottom (in French notation). As an example::

        sage: T = crystals.Tableaux(['A',2], shape = [3,2])
        sage: T.module_generators[0]
        [[1, 1, 1], [2, 2]]
        sage: list(T.module_generators[0])
        [2, 1, 2, 1, 1]

    To create a tableau, one can use::

        sage: Tab = crystals.Tableaux(['A',3], shape = [2,2])
        sage: Tab(rows=[[1,2],[3,4]])
        [[1, 2], [3, 4]]
        sage: Tab(columns=[[3,1],[4,2]])
        [[1, 2], [3, 4]]

    .. TODO::

        FIXME:

        - Do we want to specify the columns increasingly or
          decreasingly? That is, should this be
          ``Tab(columns = [[1,3],[2,4]])``?
        - Make this fully consistent with
          :func:`~sage.combinat.tableau.Tableau`!

    We illustrate the use of a shape with a negative last entry in
    type `D`::

        sage: T = crystals.Tableaux(['D',4],shape=[1,1,1,-1])
        sage: T.cardinality()
        35
        sage: TestSuite(T).run()

    We illustrate the construction of crystals of spin tableaux when
    the partitions have half integer values in type `B` and `D`::

        sage: T = crystals.Tableaux(['B',3],shape=[3/2,1/2,1/2]); T
        The crystal of tableaux of type ['B', 3] and shape(s) [[3/2, 1/2, 1/2]]
        sage: T.cardinality()
        48
        sage: T.module_generators
        ([+++, [[1]]],)
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['D',3],shape=[3/2,1/2,-1/2]); T
        The crystal of tableaux of type ['D', 3] and shape(s) [[3/2, 1/2, -1/2]]
        sage: T.cardinality()
        20
        sage: T.module_generators
        ([++-, [[1]]],)
        sage: TestSuite(T).run()

    TESTS:

    Base cases::

        sage: T = crystals.Tableaux(['A',2], shape = [])
        sage: T.list()
        [[]]
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['C',2], shape = [1])
        sage: T.list()
        [[[1]], [[2]], [[-2]], [[-1]]]
        sage: TestSuite(T).run()

        sage: T = crystals.Tableaux(['A',2], shapes = [[],[1],[2]])
        sage: T.list()
        [[], [[1]], [[2]], [[3]], [[1, 1]], [[1, 2]], [[2, 2]], [[1, 3]], [[2, 3]], [[3, 3]]]
        sage: T.module_generators
        ([], [[1]], [[1, 1]])

        sage: T = crystals.Tableaux(['B',2], shape=[3])
        sage: T(rows=[[1,1,0]])
        [[1, 1, 0]]

    Input tests::

        sage: T = crystals.Tableaux(['A',3], shape = [2,2])
        sage: C = T.letters
        sage: list(Tab(rows    = [[1,2],[3,4]])) == [C(3),C(1),C(4),C(2)]
        True
        sage: list(Tab(columns = [[3,1],[4,2]])) == [C(3),C(1),C(4),C(2)]
        True

    For compatibility with
    :func:`~sage.combinat.crystals.tensor_product.TensorProductOfCrystals` we
    need to accept as input the internal list or sequence of elements::

        sage: list(Tab(list    = [3,1,4,2]))     == [C(3),C(1),C(4),C(2)]
        True
        sage: list(Tab(3,1,4,2))                 == [C(3),C(1),C(4),C(2)]
        True

    The next example checks whether a given tableau is in fact a valid
    type `C` tableau or not::

        sage: T = crystals.Tableaux(['C',3], shape = [2,2,2])
        sage: Tab = T(rows=[[1,3],[2,-3],[3,-1]])
        sage: Tab in T.list()
        True
        sage: Tab = T(rows=[[2,3],[3,-3],[-3,-2]])
        sage: Tab in T.list()
        False
    """

    @staticmethod
    def __classcall_private__(cls, cartan_type, shapes = None, shape = None):
        """
        Normalizes the input arguments to ensure unique representation,
        and to delegate the construction of spin tableaux.

        EXAMPLES::

            sage: T1 = crystals.Tableaux(CartanType(['A',3]), shape  = [2,2])
            sage: T2 = crystals.Tableaux(['A',3],             shape  = (2,2))
            sage: T3 = crystals.Tableaux(['A',3],             shapes = ([2,2],))
            sage: T2 is T1, T3 is T1
            (True, True)
        """
        cartan_type = CartanType(cartan_type)
        n = cartan_type.rank()
        # standardize shape/shapes input into a tuple of tuples
        assert operator.xor(shape is not None, shapes is not None)
        if shape is not None:
            shapes = (shape,)
        spin_shapes = tuple( tuple(shape) for shape in shapes )
        try:
            shapes = tuple( tuple(trunc(i) for i in shape) for shape in spin_shapes )
        except Exception:
            raise ValueError("shapes should all be partitions or half-integer partitions")
        if spin_shapes == shapes:
            return super(CrystalOfTableaux, cls).__classcall__(cls, cartan_type, shapes)

        # Handle the construction of a crystals of spin tableaux
        # Caveat: this currently only supports all shapes being half
        # integer partitions of length the rank for type B and D. In
        # particular, for type D, the spins all have to be plus or all
        # minus spins
        if any(len(sh) != n for sh in shapes):
            raise ValueError("the length of all half-integer partition shapes should be the rank")
        if any(2*i % 2 != 1 for shape in spin_shapes for i in shape):
            raise ValueError("shapes should be either all partitions or all half-integer partitions")
        if cartan_type.type() == 'D':
            if all( i >= 0 for shape in spin_shapes for i in shape):
                S = CrystalOfSpinsPlus(cartan_type)
            elif all(shape[-1]<0 for shape in spin_shapes):
                S = CrystalOfSpinsMinus(cartan_type)
            else:
                raise ValueError("in type D spins should all be positive or negative")
        else:
            if any( i < 0 for shape in spin_shapes for i in shape):
                raise ValueError("shapes should all be partitions")
            S = CrystalOfSpins(cartan_type)
        B = CrystalOfTableaux(cartan_type, shapes=shapes)
        T = TensorProductOfCrystals(S, B, generators=[[S.module_generators[0],x] for x in B.module_generators])
        T.rename("The crystal of tableaux of type %s and shape(s) %s"%(cartan_type, list(list(shape) for shape in spin_shapes)))
        T.shapes = spin_shapes
        return T


    def __init__(self, cartan_type, shapes):
        """
        Construct the crystal of all tableaux of the given shapes.

        INPUT:

        - ``cartan_type`` -- (data coercible into) a Cartan type
        - ``shapes``      -- a list (or iterable) of shapes
        - ``shape``       -- a shape

        Shapes themselves are lists (or iterable) of integers.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: TestSuite(T).run()
        """
#        super(CrystalOfTableaux, self).__init__(category = FiniteEnumeratedSets())
        Parent.__init__(self, category = ClassicalCrystals())
        self.letters = CrystalOfLetters(cartan_type)
        self.shapes = shapes
        self.module_generators = tuple(self.module_generator(la) for la in shapes)
        self.rename("The crystal of tableaux of type %s and shape(s) %s"%(cartan_type, list(list(shape) for shape in shapes)))

    def cartan_type(self):
        """
        Returns the Cartan type of the associated crystal

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: T.cartan_type()
            ['A', 3]
        """
        return self.letters.cartan_type()

    def module_generator(self, shape):
        """
        This yields the module generator (or highest weight element) of a classical
        crystal of given shape. The module generator is the unique tableau with equal
        shape and content.

        EXAMPLES::

            sage: T = crystals.Tableaux(['D',3], shape = [1,1])
            sage: T.module_generator([1,1])
            [[1], [2]]

            sage: T = crystals.Tableaux(['D',4],shape=[2,2,2,-2])
            sage: T.module_generator(tuple([2,2,2,-2]))
            [[1, 1], [2, 2], [3, 3], [-4, -4]]
            sage: T.cardinality()
            294
            sage: T = crystals.Tableaux(['D',4],shape=[2,2,2,2])
            sage: T.module_generator(tuple([2,2,2,2]))
            [[1, 1], [2, 2], [3, 3], [4, 4]]
            sage: T.cardinality()
            294
        """
        type = self.cartan_type()
        if type[0] == 'D' and len(shape) == type[1] and shape[type[1]-1] < 0:
            invert = True
            shape = shape[:-1] + (-shape[type[1]-1],)
        else:
            invert = False
        p = Partition(shape).conjugate()
        # The column canonical tableau, read by columns
        module_generator = flatten([[val-i for i in range(val)] for val in p])
        if invert:
            module_generator = [(-x if x == type[1] else x) for x in module_generator]
        return self(list=[self.letters(x) for x in module_generator])

    def _element_constructor_(self, *args, **options):
        """
        Return a
        :class:`~sage.combinat.crystals.tensor_product.CrystalOfTableauxElement`.

        EXAMPLES::

            sage: T = crystals.Tableaux(['A',3], shape = [2,2])
            sage: T(rows=[[1,2],[3,4]])
            [[1, 2], [3, 4]]
            sage: T(columns=[[3,1],[4,2]])
            [[1, 2], [3, 4]]
        """
        return self.element_class(self, *args, **options)

    class Element(CrystalOfTableauxElement):
        pass

# deprecations from trac:18555
from sage.misc.superseded import deprecated_function_alias
TensorProductOfCrystals.global_options=deprecated_function_alias(18555, TensorProductOfCrystals.options)
TensorProductOfCrystalsOptions=deprecated_function_alias(18555, TensorProductOfCrystals.options)

