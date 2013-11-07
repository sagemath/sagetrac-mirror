r"""
Hall-Littlewood Young Tableaux

This is an implementation of Inka Klostermann's algorithm to compute the
coefficient of a Hall-Littlewood polynomial in the monomial basis. This is
a combinatorialization of the Gaussent-Littelmann one-skeleton gallery model
used to compute the coefficients. Currently this only works for types `A_n`,
`B_n`, and `C_n` for all weights and `D_n` when both `\omega_{n-1}` and
`\omega_n` do not appear in the weight.

AUTHORS:

- Travis Scrimshaw (2012-02-07): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
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
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets

from sage.combinat.combinat import CombinatorialObject
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.tableau import SemistandardTableaux
from sage.combinat.partition import Partition
from sage.rings.all import QQ, ZZ, infinity, factorial, gcd
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.sets.set import Set

def two_column_LS_tableaux_iter(ct, len_col1, len_col2):
    """
    Method for iterating over 2 column tableaux.
    """
    if ct.type() == 'A':
        for x in SemistandardTableaux([len_col1, len_col2], max_entry=ct.n+1):
            yield x.conjugate()
        return

    S = Set(ct.index_set())
    for s in S.subsets(len_col1):
        for barred in s.subsets():
            col1 = list(s.symmetric_difference(barred))
            col1.sort()
            barred = list(barred)
            barred.sort()
            for val in reversed(barred):
                col1.append(-val)
            
            for second_set in s.subsets(len_col2):
                for second_barred in second_set.subsets():
                    # Check that there is an even number of barred letters
                    if ct.type() == 'C' and (len(barred) + len(second_barred)) % 2 == 1:
                        continue

                    second_unbarred = list(second_set.symmetric_difference(second_barred))
                    second_unbarred.sort()
                    is_SSYT = True
                    for i, val in enumerate(second_unbarred):
                        if col1[i] > val or col1[i] < 0:
                            is_SSYT = False
                            break

                    if is_SSYT:
                        second_barred = list(second_barred)
                        second_barred.sort()
                        num_unbarred = len(second_unbarred)
                        col2 = second_unbarred[:]
                        for i, val in enumerate(reversed(second_barred)):
                            if val < col1[i+num_unbarred] and col1[i+num_unbarred] < 0:
                                is_SSYT = False
                                break
                            # If
                            col2.append(-val)
                        # For

                        if is_SSYT:
                            yield [col1, col2]

class LSTableau(CombinatorialObject, Element):
    """
    An LS-tableau.

    .. WARNING::

        The tableau is currently inputted in as a list of columns.
    """
    def __init__(self, parent, tab):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: t = L([[1,2,3],[1,3],[2,4]])
            sage: TestSuite(t).run()
        """
        Element.__init__(self, parent)
        CombinatorialObject.__init__(self, tab)

    def get_partial_product(self, i):
        """
        Get the contribution from the `i`-th and `(i+1)`-th columns to the
        product.

        This is equivalent to calling
        :meth:`~sage.combinat.hall_littlewood_yt.HallLittlewoodTree.get_sum()`
        on the Hall-Littlewood tree generated from the `i`-th and `(i+1)`-th
        columns.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: t = L([[1,2,3],[1,3],[2,4]])
            sage: t.get_partial_product(1)
            q^5 - 2*q^4 + q^3
        """
        if i == len(self) - 1:
            return HallLittlewoodTree(self.parent()._cartan_type, self[i], [],
                                      self.parent()._use_weyl_group).get_sum()
        return HallLittlewoodTree(self.parent()._cartan_type, self[i], self[i+1],
                                  self.parent()._use_weyl_group).get_sum()

    @cached_method
    def get_product(self):
        """
        Return the product over ``self``.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: t = L([[1,2,3],[1,3],[2,4]])
            sage: t.get_product()
            q^14 - 2*q^13 + q^12
        """
        prod = 1
        for i in range(len(self)):
            prod *= self.get_partial_product(i)
        return prod

    def get_partial_factorization(self, i):
        """
        Return the factorization of the partial product of the `i`-th and
        `(i+1)`-th columns of ``self``. See :meth:`get_partial_product()`.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: t = L([[1,2,3],[1,3],[2,4]])
            sage: t.get_partial_factorization(1)
            (q - 1)^2 * q^3
            sage: t.get_partial_product().factor()
            (q - 1)^2 * q^3
        """
        product = self.get_partial_product(i)
        if product == 1:
            return product
        return product.factor()

    @cached_method
    def get_factorization(self):
        """
        Return the factorization of the product of ``self``. See
        :meth:`get_product()`.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: t = L([[1,2,3],[1,3],[2,4]])
            sage: t.get_factorization()
            (q - 1)^2 * q^12
        """
        product = self.get_product()
        if product == 1:
            return product
        return product.factor()

class LSTableaux(Parent, UniqueRepresentation):
    """
    All LS tableaux.

    .. TODO::

        Make this work for types `B_n`, `C_n`, and `D_n`.
    """
    @staticmethod
    def __classcall_private__(cls, cartan_type, shape, weight, use_weyl_group=True):
        """
        Normalize inputs to ensure a unique representation.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: L2 = LSTableaux(CartanType(['A',4]), (3,3,1), (2,2,2,1))
            sage: L is L2
            True
        """
        cartan_type = CartanType(cartan_type)
        if not cartan_type.is_finite():
            raise ValueError("The Cartan type needs to be a finite type")
        return super(LSTableaux, cls).__classcall__(cls, cartan_type, Partition(shape),
                                                    Partition(weight), use_weyl_group)

    def __init__(self, cartan_type, shape, weight, use_weyl_group):
        """
        Initialize ``self``.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: TestSuite(L).run()
        """
        self._cartan_type = cartan_type
        self._shape = shape
        self._weight = weight
        self._use_weyl_group = use_weyl_group
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::

            sage: LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            Hall-Littlewood trees of Cartan type ['A', 4], shape [3, 3, 1], and weight [2, 2, 2, 1]
        """
        return "Hall-Littlewood trees of Cartan type %s, shape %s, and weight %s"%(
            repr(self._cartan_type), self._shape, self._weight)

    def _element_constructor_(self, tab):
        """
        Construct an element of ``self``.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: L([[1,2,3],[1,3],[2,4]])
            [[1, 2, 3], [1, 3], [2, 4]]
        """
        #if len(first_col) != self._len_first or len(second_col) != self._len_second:
        #    raise ValueError("Column hieght incorrect")
        return self.element_class(self, tab)

    Element = LSTableau

    def __iter__(self):
        """
        Iterate through ``self``.

        For type `A_n`, this just uses the semistandard tableaux iterator and
        then converts it to columns.

        This is only implemented in types `B_n` and `C_n` only if the shape has
        2 columns.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: for x in L: x
            [[1, 2, 4], [1, 3], [2, 3]]
            [[1, 2, 3], [1, 3], [2, 4]]
            [[1, 2, 3], [1, 2], [3, 4]]
        """
        if self._cartan_type.type() == 'A':
            for x in SemistandardTableaux(shape=self._shape, eval=self._weight):
                yield self.element_class(self, x.conjugate())
            return
        conj = self._shape.conjugate()
        if len(conj) != 2:
            raise NotImplementedError("Only implemented for partitions of width 2")
        for x in two_column_LS_tableaux_iter(self._cartan_type, conj[0], conj[1]):
            yield self.element_class(self, x)

    def factor_iter(self):
        """
        Convienence method for iterating over the factorizations of elements
        in ``self``.

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: for x in L.factor_iter(): x
            (q - 1)^2 * q^12
            (q - 1)^2 * q^12
            (q + 1) * (q - 1)^2 * q^11
        """
        for x in self:
            yield x.get_factorization()

    @cached_method
    def coefficient_polynomial(self):
        r"""
        Return the coefficient polynomial of ``self``.

        Given a shape `\lambda` and weight `\mu`, the coefficient
        polynomial is defined as:

        .. MATH::

            L_{\lambda, \mu} = \sum_{T \in SSYT(\lambda, \mu)} c(C_r, \emptyset)
            \prod_{i=0}^{r-1} c(C_i, C_{i+1})

        where `C_i` is the `i`-th column of `T` and `c(A, B)` is the result
        obtained from
        :meth:`sage.combinat.hall_littlewood_yt.HallLittlewoodTree.get_sum`
        using the two column tableau `(A, B)`. Note that the set `SSYT` depends
        upon the Cartan type.

        These arrise from considering a Hall-Littlewood polynomial `P_{\lambda}`
        and expaning it in terms of the monomial symmetric functions:

        .. MATH::

            P_{\lambda}(x; q) = \sum_{\mu \in X_+^{\check}}
            q^{-\langle \lambda+\mu, \rho \rangle} L_{\lambda, \mu}(q)
            m_{\mu}(x)

        .. NOTE::

            One must substitute `q = t^{-1}` for the standard notation of
            Hall-Littlewood polynomials (for example, to agree with the
            Macdonald formula).

        EXAMPLES::

            sage: L = LSTableaux(['A',4], [3,3,1], [2,2,2,1])
            sage: L.coefficient_polynomial()
            3*q^14 - 5*q^13 + q^12 + q^11
        """
        cur_sum = 0
        for x in self:
            cur_sum += x.get_product()
        return cur_sum

class HallLittlewoodTree(SageObject):
    """
    A Hall-Littlewood tree.

    This is the underlying tree structure that is used to compute the
    contributions of the pairs of columns in Klostermann's algorithm.
    """
    def __init__(self, cartan_type, first_column, second_column, use_weyl_group):
        """
        Initialize ``self``.
        """
        if len(first_column) < len(second_column):
            raise ValueError("The first column must be as least as tall as the second")
        self._cartan_type = cartan_type
        ct_type = cartan_type.type()
        self.root_tableau = [first_column, second_column]
        if ct_type == 'A':
            self.root = HallLittlewoodTreeNodeTypeA(self, self.root_tableau, -1, use_weyl_group)
        elif ct_type == 'B':
            self.root = HallLittlewoodTreeNodeTypeB(self, self.root_tableau, -1, use_weyl_group)
        elif ct_type == 'C':
            self.root = HallLittlewoodTreeNodeTypeC(self, self.root_tableau, -1, use_weyl_group)
        elif ct_type == 'D':
            self.root = HallLittlewoodTreeNodeTypeD(self, self.root_tableau, -1, use_weyl_group)
        else:
            raise NotImplementedError("Only implemented for types A, B, C, and D")
    
    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Hall-Littlewood tree of %s"%repr(self.root_tableau)
    
    def __iter__(self):
        """
        Iterate through ``self``.
        """
        # Do depth first search
        path = [self.root]
        cur_child = [-1] # For adding 1 initially
        yield self.root
        while len(cur_child) > 0:
            cur_child[-1] += 1
            if cur_child[-1] < len(path[-1].children):
                path.append(path[-1].children[cur_child[-1]])
                cur_child.append(-1)
                yield path[-1]
            else:
                path.pop()
                cur_child.pop()

    def get_sum(self):
        """
        Return the sum of the products of the edge labels.
        """
        # Do depth first search to build out all paths
        # If there are no children off the root
        if len(self.root.children) == 0:
            return 1

        ZZq = PolynomialRing(ZZ, ['q'])
        q = ZZq.gens()[0]
        sum = 0

        path = [self.root]
        cur_child = [-1] # For adding 1 initially
        num_type = [0, 0]
        while len(cur_child) > 0:
            cur_child[-1] += 1
            if cur_child[-1] < len(path[-1].children):
                path.append(path[-1].children[cur_child[-1]])
                cur_child.append(-1)
                if path[-1].up_edge_type >= 0:
                    num_type[path[-1].up_edge_type] += 1
            else:
                # If it has no children, then its a leaf
                if len(path[-1].children) == 0:
                    sum += q**num_type[0] * (q - 1)**num_type[1]
                if path[-1].up_edge_type >= 0:
                    num_type[path[-1].up_edge_type] -= 1
                path.pop()
                cur_child.pop()
        return sum

class HallLittlewoodTreeNodeAbstract(CombinatorialObject):
    r"""
    Abstract root class for all of the so-called Hall-Littlewood tree nodes.
    """
    
    def __init__(self, tree, tableau, type, use_weyl_group):
        r"""
        type = -1 undefined or `s_i^-`
        type = 0 is `s_i^+`
        type = 1 is `id_i^+`
        
        use_weyl_group is for types `B_n` and `C_n`
        """
        CombinatorialObject.__init__(self, tableau)
        self.tree = tree
        self.up_edge_type = type
        
        # If we are the empty tableau, there will be no children
        if len(tableau[0]) == 0:
            self.children = []
        else:
            self._build_children(use_weyl_group)
    
    @abstract_method
    def _build_children(self, use_weyl_group):
        r"""
        Build and assign as our children.
        """

class HallLittlewoodTreeNodeTypeA(HallLittlewoodTreeNodeAbstract):
    r"""
    Node of the so-called Hall-Littlewood tree of type `A_n`.
    """
    
    def _build_children(self, use_weyl_group):
        r"""
        Build and assign as our children.
        """
        height = len(self._list[0])
        swap_index = None
        for i in range(height-1):
            if self._list[0][i] + 1 != self._list[0][i+1]:
                swap_index = i
                break
        
        # There is not an increasing in the first height-1 entries
        if swap_index is None:
            # Check to see if we can increase the last entry
            if self._list[0][-1] != self.tree._cartan_type.n + 1:
                swap_index = height - 1
            else:
                self.children = []
                return
        
        child_tableau = [self._list[0][:], self._list[1][:]] # Make a deep copy
        value = child_tableau[0][swap_index]
        
        # Perform the swap on the first column (this will always be returned)
        child_tableau[0][swap_index] += 1
        
        # Check the second column
        for i, val in enumerate(self._list[1]):
            if val == value:
                # If we find the same value, we must be non-decreasing
                # Check to see if the swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != value + 1:
                    child_tableau[1][i] += 1
                self.children = [self.__class__(self.tree, child_tableau,
                                                             0, use_weyl_group)]
                return
            elif val == value + 1:
                # We must be swapping down since columns are strictly increasing
                # Create the branch for `id_i^+` (it will always be created)
                self.children = [self.__class__(self.tree, child_tableau,
                                                             1, use_weyl_group)]
                
                # If we get a SSYT
                if child_tableau[0][i] <= val-1:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    self.children.append(HallLittlewoodTreeNodeTypeA(self.tree,
                                         second_child_tableau, -1, use_weyl_group))
                return
            elif val > value + 1:
                self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]
                return
        self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]

class HallLittlewoodTreeNodeTypeB(HallLittlewoodTreeNodeAbstract):
    r"""
    Node of the so-called Hall-Littlewood tree of type `B_n`.
    """

    def _build_children(self, use_weyl_group):
        r"""
        Build the children nodes.
        
        `s_i` acts here by interchanging `i \leftrightarrow i+1` and 
        `\overline{i} \leftrightarrow \overline{i+1}` for all `i` except the 
        highest index if we are not using the full Weyl group, otherwise as
        long as `i \neq n`. For the exceptional case `k`, we have
        `k \leftrightarrow \overline{k}`.
        """
        height = len(self._list[0])
        # Throughout here we will be using the fact that both i and -i cannot
        # appear in the same column.
        n = self.tree._cartan_type.n
        
        swap_index = None
        if use_weyl_group:
            if self._list[0][0] == -n:
                for i in range(1, height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
            else:
                for i in range(height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
        else:
            if self._list[0][0] < 0:
                # If the 1st entry is barred, then it will be the max index
                #   and will become unbarred (hence decrease) so do not try
                for i in range(1, height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
            else:
                # Otherwise go through normally
                for i in range(height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break

        # There is not an increasing swap in the first height-1 entries
        if swap_index is None:
            # Check to see if we can increase the last entry
            if self._list[0][-1] != -1 and (use_weyl_group or height != 1
                                            or self._list[0][0] > 0):
                swap_index = height - 1
            else:
                self.children = []
                return
        
        child_tableau = [self._list[0][:], self._list[1][:]] # Make a deep copy
        value = child_tableau[0][swap_index]
        
        # Perform the swap on the first column (this will always be returned)
        # Note that the value must be positive to be a bar swap because it
        #   would be decreasing otherwise
        if use_weyl_group:
            if value == n:
                child_tableau[0][swap_index] = -n
                self._bar_swap(child_tableau, n, use_weyl_group)
                return
        else:
            if swap_index < height - 1:
                pos_next_value = child_tableau[0][swap_index+1]
                if pos_next_value < 0:
                    pos_next_value = -pos_next_value

                if pos_next_value < value:
                    child_tableau[0][swap_index] = -value
                    self._bar_swap(child_tableau, value, use_weyl_group)
                    return
            elif value > 0:
                # If it is the last (largest) entry, it must be a bar swap
                # Note that if it is the last entry and positive, all others
                #   must be positive.
                child_tableau[0][swap_index] = -value
                self._bar_swap(child_tableau, value, use_weyl_group)
                return

        # We did not perform a bar swap
        child_tableau[0][swap_index] += 1
        if value > 0:
            # Look for a -(value+1) to swap to -value
            for i in range(swap_index+1, height):
                if child_tableau[0][i] == -value-1:
                    child_tableau[0][i] += 1
                    break
        else:
            value = -value-1

        # Check the second column
        for i, val in enumerate(self._list[1]):
            if val == value:
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != value + 1:
                    child_tableau[1][i] += 1
                
                # Check for a barred swap
                for j in range(i+1, len(self._list[1])):
                    if child_tableau[1][j] == -value - 1:
                        child_tableau[1][j] += 1
                        break
                break
            elif val == value + 1:
                self.children = [self.__class__(self.tree, child_tableau,
                                                             1, use_weyl_group)]
                
                # If we get a SSYT
                if child_tableau[0][i] <= value:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    # Check for barred letters
                    for j in range(i, len(self._list[1])):
                        if child_tableau[1][j] == -value:
                            # Check if the barred letters give us a SSYT
                            if 0 > child_tableau[0][j] and child_tableau[0][j] > -value-1:
                                return
                            second_child_tableau[1][j] -= 1
                            break
                    self.children.append(self.__class__(
                                         self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif val == -(value+1):
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != -value:
                    child_tableau[1][i] += 1
                break
            elif val == -value:
                self.children = [self.__class__(self.tree, child_tableau,
                                                              1, use_weyl_group)]
                # Check if SSYT
                if 0 < child_tableau[0][i] or child_tableau[0][i] <= -value-1:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    self.children.append(self.__class__(
                                         self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif 0 > val and val > -value:
                break
        self.children = [self.__class__(self.tree, child_tableau,
                                                     0, use_weyl_group)]
    
    def _bar_swap(self, child_tableau, value, use_weyl_group):
        r"""
        Perform a bar swap.
        """
        # Recall value > 0
        for i, val in enumerate(self._list[1]):
            # Only 1 operation to do in each column if a bar swap
            if val == value:
                child_tableau[1][i] = -value
                break
            elif val == -value:
                self.children = [self.__class__(self.tree, child_tableau,
                                                              1, use_weyl_group)]
                # Check that it is a SSYT
                if 0 <= child_tableau[0][i] and child_tableau[0][i] <= value:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] = value
                    self.children.append(self.__class__(self.tree,
                                         second_child_tableau, -1, use_weyl_group))
                # If
                return
            elif 0 > val and val > -value:
                break
        self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]

class HallLittlewoodTreeNodeTypeC(HallLittlewoodTreeNodeTypeB):
    r"""
    Node of the so-called Hall-Littlewood tree of type `C_n`.
    
    For the tableau:

    - Double every column except the height 1 column
    - If `i \in C_j` then `\overline{i} \notin C_j`
    - Number of bars in `C_{2k}` to make `C_{2k+1}` must be even
    """

    def _build_children(self, use_weyl_group):
        r"""
        Build the children in the tree.
        
        `s_i` acts as in type `B_n` except when not in the full Weyl group,
        we have `s_m` acts by `m \leftrightarrow \overline{m-1}` and
        `m-1 \leftrightarrow \overline{m}`.
        """
        height = len(self._list[0])
        # Throughout here we will be using the fact that both i and -i cannot
        # appear in the same column.
        n = self.tree._cartan_type.n
        
        swap_index = None
        if use_weyl_group:
            if self._list[0][0] == -n:
                for i in range(1, height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
            else:
                for i in range(height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
        else:
            # We are in the partial Weyl group
            if self._list[0][0] < 0:
                # If the 1st entry is barred, then it will be the max index
                #   and will become unbarred (hence decrease) so do not try
                for i in range(1, height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1]: 
                        swap_index = i
                        break
            else:
                # Otherwise go through normally
                for i in range(height-1):
                    if self._list[0][i] + 1 != self._list[0][i+1] \
                      and -(self._list[0][i] - 1) != self._list[0][i+1]:
                        # Check for special case of m -(m-1) since the swap does nothing
                        swap_index = i
                        break
        
        # There is not an increasing swap in the first height-1 entries
        if swap_index is None:
            # Check to see if we can increase the last entry
            if self._list[0][-1] != -1 and (use_weyl_group or height != 1):
                swap_index = height - 1
            else:
                self.children = []
                return
        
        child_tableau = [self._list[0][:], self._list[1][:]] # Make a deep copy
        value = child_tableau[0][swap_index]
        
        # Perform the swap on the first column (this will always be returned)
        # Note that the value must be positive to be a bar swap because it
        #   would be decreasing otherwise
        if use_weyl_group:
            if value == n:
                child_tableau[0][swap_index] = -n
                self._bar_swap(child_tableau, n, use_weyl_group)
                return
        else:
            if swap_index < height - 1:
                pos_next_value = child_tableau[0][swap_index+1]
                if pos_next_value < 0:
                    pos_next_value = -pos_next_value

                if pos_next_value < value:
                    child_tableau[0][swap_index] = -value+1
                    if swap_index > 0 and child_tableau[0][swap_index-1] == value-1:
                        child_tableau[0][swap_index-1] = -value
                    self._bar_swap(child_tableau, value, use_weyl_group)
                    return
                    # Never need to swap -value+1 back to value (would imply a
                    #   skipped pos index).
            elif value > 0:
                # If it is the last (largest) entry, it must be a bar swap
                # Note that if it is the last entry and positive, all others
                #   must be positive.
                child_tableau[0][swap_index] = -value+1
                self._bar_swap(child_tableau, value, use_weyl_group)
                return

        # We did not perform a bar swap
        child_tableau[0][swap_index] += 1
        if value > 0:
            # Look for a -(value+1) to swap to -value
            for i in range(swap_index+1, height):
                if child_tableau[0][i] == -value-1:
                    child_tableau[0][i] += 1
                    break
        else:
            # Make sure we are taking the correct value
            value = -value-1

        # Check the second column
        for i, val in enumerate(self._list[1]):
            if val == value:
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != value + 1:
                    child_tableau[1][i] += 1
                
                # Check for a barred swap
                for j in range(i+1, len(self._list[1])):
                    if child_tableau[1][j] == -value - 1:
                        child_tableau[1][j] += 1
                        break
                
                break
            elif val == value + 1:
                self.children = [self.__class__(self.tree, child_tableau, 1, use_weyl_group)]
                
                # If we get a SSYT
                if child_tableau[0][i] <= value:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    # Check for barred letters
                    for j in range(i, len(self._list[1])):
                        if child_tableau[1][j] == -value:
                            # Check if the barred letters give us a SSYT
                            if 0 > child_tableau[0][j] and child_tableau[0][j] > -value-1:
                                return
                            second_child_tableau[1][j] -= 1
                            break
                    self.children.append(self.__class__(self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif val == -value-1:
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != -value:
                    child_tableau[1][i] += 1
                break
            elif val == -value:
                self.children = [self.__class__(self.tree, child_tableau, 1, use_weyl_group)]
                # Check if SSYT
                if 0 < child_tableau[0][i] or child_tableau[0][i] <= -value-1:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    self.children.append(self.__class__(self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif 0 > val and val > -value:
                break
        self.children = [self.__class__(self.tree, child_tableau,
                                                     0, use_weyl_group)]
    # _build_children()
    
    def _bar_swap(self, child_tableau, value, use_weyl_group):
        r"""
        Perform a bar swap.
        """
        # Recall value > 0
        for i, val in enumerate(self._list[1]):
            # Only 1 operation to do in each column if a bar swap
            if use_weyl_group:
                if val == value:
                    child_tableau[1][i] = -value
                    break
                elif val == -value:
                    self.children = [self.__class__(self.tree, child_tableau,
                                                                  1, use_weyl_group)]
                    # Check that it is a SSYT
                    if 0 <= child_tableau[0][i] and child_tableau[0][i] <= value:
                        second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                        second_child_tableau[1][i] = value
                        self.children.append(self.__class__(self.tree,
                                             second_child_tableau, -1, use_weyl_group))
                    return
                elif 0 > val and val > -value:
                    break
            else:
                # It is the partial Weyl group
                if val == value-1:
                    # If the next value is not -value, then it is not the max index
                    if i != len(self._list[1])-1:
                        if child_tableau[1][i+1] != -value-1:
                            child_tableau[i] = -value
                            if child_tableau[1][i+1] == value+1:
                                child_tableau[i+1] = -value+1
                    else:
                        child_tableau[i] = -value
                    break
                elif val == value:
                    if i != len(self._list[1])-1 and child_tableau[1][i+1] == -value+1:
                        break
                    
                    child_tableau[1][i] = -value+1
                    break
                elif val == -value:
                    self.children = [self.__class__(self.tree, child_tableau,
                                                                  1, use_weyl_group)]
                    # Check that it is a SSYT
                    if 0 <= child_tableau[0][i] and child_tableau[0][i] <= value-1:
                        second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                        second_child_tableau[1][i] = value-1
                        # Check to see if we need to swap the next entry (because we will
                        #   need to also check if that row is semistandard too).
                        if i != len(self._list[1])-1 and child_tableau[1][i+1] == -(value-1):
                            # Note that this will become the max pos value
                            if 0 <= child_tableau[0][i+1]:
                                second_child_tableau[1][i+1] = value
                                self.children.append(self.__class__(self.tree,
                                                     second_child_tableau, -1, use_weyl_group))
                        else:
                            self.children.append(self.__class__(self.tree,
                                                 second_child_tableau, -1, use_weyl_group))
                    return
                elif val == -value+1:
                    self.children = [self.__class__(self.tree, child_tableau,
                                                                  1, use_weyl_group)]
                    # Check that it is a SSYT
                    if 0 <= child_tableau[0][i]:
                        # It will become the maximum pos value in the left column
                        # So as long as the left column entry is pos, the it
                        #   will be a SSYT.
                        second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                        second_child_tableau[1][i] = value
                        self.children.append(self.__class__(self.tree,
                                             second_child_tableau, -1, use_weyl_group))
                    return
        self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]

class HallLittlewoodTreeNodeTypeD(HallLittlewoodTreeNodeTypeA):
    r"""
    Node of the so-called Hall-Littlewood tree of type `D_n`.
    
    For the tableau:

    - Height of each column at most `n-2` (for now)
    """
    def _build_children(self, use_weyl_group):
        r"""
        Build the children in the tree.
        
        `s_i` acts as in type `B_n` except `s_n` acts by
        `n \leftrightarrow \overline{n-1}` and `n-1 \leftrightarrow \overline{n}`.
        """
        height = len(self._list[0])
        # Throughout here we will be using the fact that both i and -i cannot
        # appear in the same column.
        n = self.tree._cartan_type.n
        
        swap_index = None
        for i in range(height-1):
            if self._list[0][i] + 1 != self._list[0][i+1] \
              and -(self._list[0][i] - 1) != self._list[0][i+1]:
                # Check for special case of m -(m-1) since the swap does nothing
                swap_index = i
                break

        # There is not an increasing swap in the first height-1 entries
        if swap_index is None:
            # Check to see if we can increase the last entry
            if self._list[0][-1] != -1:
                swap_index = height - 1
            else:
                self.children = []
                return
        
        child_tableau = [self._list[0][:], self._list[1][:]] # Make a deep copy
        value = child_tableau[0][swap_index]
        
        # Perform the swap on the first column (this will always be returned)
        # ----------
        
        # Note that the value must be positive to be a bar swap because it
        #   would be decreasing otherwise
        # Perform the bar swap
        if value == n:
            child_tableau[0][swap_index] = -(n-1)
            
            for i, val in enumerate(self._list[1]):
                # Only 1 operation to do in second column
                if val == n-1:
                    # We only increase
                    child_tableau[i] = -n
                    if i != len(self._list[1])-1:
                        if child_tableau[1][i+1] == n+1:
                            child_tableau[i+1] = -(n-1)
                    break
                elif val == n:
                    # Check if we make a net change
                    if i != len(self._list[1])-1 and child_tableau[1][i+1] == -(n-1):
                        break
                    
                    child_tableau[1][i] = -(n-1)
                    break
                elif val == -n:
                    self.children = [self.__class__(self.tree, child_tableau,
                                                                  1, use_weyl_group)]
                    # Check that it is a SSYT
                    if 0 <= child_tableau[0][i] and child_tableau[0][i] <= n-1:
                        second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                        second_child_tableau[1][i] = n-1
                        # Check to see if we need to swap the next entry (because we will
                        #   need to also check if that row is semistandard too).
                        if i != len(self._list[1])-1 and child_tableau[1][i+1] == -(n-1):
                            # Note that this will become the max pos value
                            if 0 <= child_tableau[0][i+1]:
                                second_child_tableau[1][i+1] = n
                                self.children.append(self.__class__(self.tree,
                                                     second_child_tableau, -1, use_weyl_group))
                        else:
                            self.children.append(self.__class__(self.tree,
                                                 second_child_tableau, -1, use_weyl_group))
                    return
                elif val == -(n-1):
                    self.children = [self.__class__(self.tree, child_tableau,
                                                                  1, use_weyl_group)]
                    # Check that it is a SSYT
                    if 0 <= child_tableau[0][i]:
                        # It will become the maximum pos value in the left column
                        # So as long as the left column entry is pos, the it
                        #   will be a SSYT.
                        second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                        second_child_tableau[1][i] = n
                        self.children.append(self.__class__(self.tree,
                                             second_child_tableau, -1, use_weyl_group))
                    return
            self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]
            return

        # We did not perform a bar swap
        # ----------
        
        child_tableau[0][swap_index] += 1
        if value > 0:
            # Look for a -(value+1) to swap to -value
            for i in range(swap_index+1, height):
                if child_tableau[0][i] == -value-1:
                    child_tableau[0][i] += 1
                    break
        else:
            # Make sure we are taking the correct value
            value = -value - 1

        # Check the second column
        for i, val in enumerate(self._list[1]):
            if val == value:
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != value + 1:
                    child_tableau[1][i] += 1
                
                # Check for a barred swap
                for j in range(i+1, len(self._list[1])):
                    if child_tableau[1][j] == -value - 1:
                        child_tableau[1][j] += 1
                        break
                break
            elif val == value + 1:
                self.children = [self.__class__(self.tree, child_tableau, 1, use_weyl_group)]
                
                # If we get a SSYT
                if child_tableau[0][i] <= value:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    # Check for barred letters
                    for j in range(i, len(self._list[1])):
                        if child_tableau[1][j] == -value:
                            # Check if the barred letters give us a SSYT
                            if 0 > child_tableau[0][j] and child_tableau[0][j] > -value-1:
                                return
                            second_child_tableau[1][j] -= 1
                            break
                    self.children.append(self.__class__(self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif val == -value - 1:
                # The swap is non-trivial
                if i == len(self._list[1]) - 1 or child_tableau[1][i+1] != -value:
                    child_tableau[1][i] += 1
                break
            elif val == -value:
                self.children = [self.__class__(self.tree, child_tableau,
                                                              1, use_weyl_group)]
                # Check if SSYT
                if 0 < child_tableau[0][i] or child_tableau[0][i] <= -value-1:
                    second_child_tableau = [child_tableau[0][:], child_tableau[1][:]]
                    second_child_tableau[1][i] -= 1
                    self.children.append(self.__class__(
                                         self.tree, second_child_tableau, -1, use_weyl_group))
                return
            elif 0 > val and val > -value:
                break
        self.children = [self.__class__(self.tree, child_tableau, 0, use_weyl_group)]

