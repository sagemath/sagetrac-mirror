r"""
Crystals of RSK

Crystals which arrise from the RSK bijection:

- :class:`CrystalOfMatrices`
- :class:`CrystalOfBiwords`
- :class:`CrystalOfBitableaux`

AUTHORS:

- Travis Scrimshaw: Initial version

.. WARNING::

    Does not work with :func:`TensorProductOfCrystals`.

REFERENCES:

.. [KHKwon12] J.-H. Kwon.
   RSK correspondence and classically irreducible Kirillov-Reshetikhin crystals.
   :arxiv:`1110.2629v2`.
"""

#*****************************************************************************
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
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

from copy import copy
from sage.misc.abstract_method import abstract_method
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.categories.finite_crystals import FiniteCrystals
from sage.categories.regular_crystals import RegularCrystals
from sage.categories.highest_weight_crystals import HighestWeightCrystals
from sage.matrix.constructor import matrix
from sage.combinat.combinat import CombinatorialObject
from sage.combinat.tableau import Tableau
from sage.combinat.root_system.cartan_type import CartanType
from sage.combinat.root_system.root_system import RootSystem
from sage.combinat.rsk import RSK, RSK_inverse, to_matrix

# This is the ABC
class CrystalFromRSK(Parent, UniqueRepresentation):
    r"""
    Root abstract class for crystals arrising from RSK.

    INPUT:

    - ``n`` -- Type `A_{n-1}^{(1)}`
    - ``r`` -- The number of rows

    EXAMPLES::
    """
    def __init__(self, n, r, s):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: CM = CrystalOfMatrices(5, 2)
            sage: TestSuite(CM).run()
        """
        self._n = n
        self._r = r
        self._s = s
        self._cartan_type = CartanType(['A',n-1,1])
        if s is not None:
            Parent.__init__(self, category=(RegularCrystals(), FiniteCrystals()))
        else:            
            Parent.__init__(self, category=HighestWeightCrystals())

    @abstract_method
    def _check_bounded(self, obj):
        """
        Check if this is bounded by `s`.
        """

####################
# Matrices crystal #
####################

class CrystalOfMatricesElement(Element):
    r"""
    A crystal of matrices element.

    For more information, see :class:`CrystalOfMatrices`.

    EXAMPLES::
    """
    def __init__(self, parent, M):
        r"""
        EXAMPLES::

            sage: Y = CrystalOfMatrices(5, 2)
            sage: mg = Y.module_generators[0]
            sage: TestSuite(mg).run()
        """
        self._M = M
        Element.__init__(self, parent)

    def _repr_(self):
        r"""
        EXAMPLES::
        """
        return repr(self._M)

    def _latex_(self):
        r"""
        Generate LaTeX code for ``self``.  Requires TikZ.

        EXAMPLES::

            sage: x = CrystalOfMatrices(5, 2)
            sage: latex(x)
        """
        return self._M._latex_()

    def __eq__(self, rhs):
        """
        Check equality
        """
        if isinstance(rhs, CrystalOfMatricesElement):
            return self._M == rhs._M
        return self._M == rhs

    def __ne__(self, rhs):
        """
        Check not equals.
        """
        return not self.__eq__(rhs)

    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        EXAMPLES::
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r
        n = self.parent()._n
        m = copy(self._M)
        if i == 0:
            m[r-1,n-r-1] += 1
        elif i == r:
            if m[0,0] == 0:
                return None
            m[0,0] -= 1
        elif i < r:
            i -= 1 # for indexing
            found = False
            for j in reversed(range(n-r)):
                if m[i+1,j] > 0:
                    m[i,j] += 1
                    m[i+1,j] -= 1
                    found = True
                    break
            if not found:
                return None
        else:
            i -= r+1 # for indexing
            found = False
            for j in reversed(range(r)):
                if m[j,i+1] > 0:
                    m[j,i] += 1
                    m[j,i+1] -= 1
                    found = True
                    break
            if not found:
                return None
        if not self.parent()._check_bounded(m):
            return None
        m.set_immutable()
        return self.__class__(self.parent(), m)

    def f(self,i):
        r"""
        Return the action of `f_i` on ``self``.

        EXAMPLES::
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r
        n = self.parent()._n
        m = copy(self._M)
        if i == 0:
            if m[r-1,n-r-1] == 0:
                return None
            m[r-1,n-r-1] -= 1
        elif i == r:
            m[0,0] += 1
        elif i < r:
            i -= 1 # for indexing
            found = False
            for j in reversed(range(n-r)):
                if m[i,j] > 0:
                    m[i,j] -= 1
                    m[i+1,j] += 1
                    found = True
                    break
            if not found:
                return None
        else:
            i -= r+1 # for indexing
            found = False
            for j in reversed(range(r)):
                if m[j,i] > 0:
                    m[j,i] -= 1
                    m[j,i+1] += 1
                    found = True
                    break
            if not found:
                return None
        if not self.parent()._check_bounded(m):
            return None
        m.set_immutable()
        return self.__class__(self.parent(), m)

#    def weight(self):
#        r"""
#        Returns the weight of ``self`` as an element of the root lattice
#        `\bigoplus_{i=0}^n \ZZ \alpha_i`.
#
#        EXAMPLES::
#
#            sage: C = CrystalOfMatrices(5, 2)
#            sage: x = C.module_generators[0]
#            sage: x.weight()
#        """
#        pass

    @cached_method
    def to_biword(self):
        """
        Return ``self`` as a biword.
        """
        raise NotImplementedError

    @cached_method
    def to_bitableaux(self):
        """
        Return ``self`` as a bitableaux.
        """
        return RSK(self._M)

class CrystalOfMatrices(CrystalFromRSK):
    r"""
    The crystal of `r \times (n-r)` matrices defined in [JHKwon12]_.

    INPUT:

    - ``n`` -- Type `A_{n-1}^{(1)}`
    - ``r`` -- The number of rows
    - ``s`` -- (Default: ``None``) If this is not ``None``, then this is
      the maximal row/column word length
    """
    def __init__(self, n, r, s=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: CM = CrystalOfMatrices(5, 2)
            sage: TestSuite(CM).run()
        """
        CrystalFromRSK.__init__(self, n, r, s)
        self.module_generators = [self.element_class( self, matrix.zero(r, n-r) )]

    Element = CrystalOfMatricesElement

    def _element_constructor_(self, data):
        r"""
        Construct an element of ``self`` from ``data``.

        INPUT:

        - ``data`` -- a multilist

        EXAMPLES::

            sage: CM = CrystalOfMatrices(5, 2)
        """
        if isinstance(data, CrystalOfMatricesElement):
            if data.parent() is not self:
                raise ValueError("Unable to change parents")
            return data
        return self.element_class(self, matrix(data))

    def _repr_(self):
        r"""
        EXAMPLES::

            sage: CrystalOfMatrices(5, 2)
            Crystal of 2 x 3 matrices
            sage: CrystalOfMatrices(5, 2, 2)
            Crystal of 2 x 3 matrices bounded by 2
        """
        ret = "Crystal of %s x %s matrices"%(self._r, self._n-self._r)
        if self._s is not None:
            ret += " bounded by %s"%self._s
        return ret

    def _check_bounded(self, m):
        """
        Check if the matrix ``m`` is bounded by `s`.
        """
        if self._s is not None:
            for r in m.rows():
                if sum(r) > self._s:
                    return False
            for r in m.columns():
                if sum(r) > self._s:
                    return False
        return True

###################
# Biwords crystal #
###################

class CrystalOfBiwordsElement(CombinatorialObject, Element):
    """
    An element in a crystal of biwords.
    """
    def __init__(self, parent, top, bottom):
        """
        Initialize ``self``.
        """
        CombinatorialObject.__init__(self, [top, bottom])
        Element.__init__(self, parent)

    def _latex_(self):
        """
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: C = CrystalOfBiwords(8, 4)
            sage: latex(C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3]))
        """
        top = " & ".join(repr(x) for x in self._list[0])
        bot = " & ".join(repr(x) for x in self._list[1])
        return "\\begin{pmatrix} %s \\\\ %s \\end{pmatrix}"%(top, bot)

    def _repr_array(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: C = CrystalOfBiwords(8, 4)
            sage: print C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])._repr_array()
        """
        return ''.join('%3s'%x for x in self._list[0]) \
            + '\n' + ''.join('%3s'%x for x in self._list[1])

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: C = CrystalOfBiwords(8, 4)
            sage: C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3]).pp()
        """
        print self._repr_array()

    def e(self, i):
        """
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = CrystalOfBiwords(8, 4)
            sage: elt = C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])
            sage: elt.e(1)
            sage: elt.e(5)
            [[1, 1, 1, 1, 2, 2, 2], [1, 2, 2, 4, 1, 2, 3]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r

        if i == 0:
            top = self._list[0][:]
            bot = self._list[1][:]
            top.append(r)
            bot.append(self.parent()._n-r)
        elif i == r:
            top = self._list[0][:]
            bot = self._list[1][:]
            if len(top) == 0 or top[0] != 1 or bot[0] != 1:
                return None
            top.pop()
            bot.pop()
        elif i < r:
            top = self._list[0][:]
            bot = self._list[1][:]
            k = len(bot)
            if k == 0:
                return None
            pairring_marker = None # Position where we're doing a pairring comparison
            pos = None # Position of unpairred biletter
            step_pos = None # Where the step from i+1 to i occurs
            for j in reversed(range(k)):
                if top[j] == i+1:
                    if pairring_marker is None:
                        pairring_marker = j
                elif top[j] == i:
                    if step_pos is None:
                        if pairring_marker is None: # No (i+1)'s found
                            return None
                        step_pos = j
                    while pairring_marker > step_pos and bot[pairring_marker] >= bot[j]:
                        pos = pairring_marker
                        pairring_marker -= 1
                    if pairring_marker == step_pos:
                        break
                    else:
                        pairring_marker -= 1
                elif top[j] < i:
                    break
            if pos is None:
                return None
            if pairring_marker > step_pos:
                pos = step_pos + 1
            top[step_pos+1] -= 1
            # Shuffle the entries along
            val = bot[pos]
            pos -= 1
            while pos >= 0 and (pos > step_pos or (bot[pos] > val and top[pos] == i)):
                bot[pos+1] = bot[pos]
                pos -= 1
            bot[pos+1] = val
        else:
            i -= r
            bot = self._list[1][:]
            if len(bot) == 0:
                return None
            pos = None
            count = 0
            for j in reversed(range(len(bot))):
                if bot[j] == i+1:
                    if count == 0:
                        pos = j
                    else:
                        count -= 1
                elif bot[j] == i:
                    count += 1
            if pos is None:
                return None
            bot[pos] -= 1
            top = self._list[0][:]
        if self.parent()._check_bounded([top, bot]):
            return self.__class__(self.parent(), top, bot)
        return None

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = CrystalOfBiwords(8, 4)
            sage: elt = C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])
            sage: elt.f(1)
            sage: elt.f(5)
            [[1, 1, 2, 2, 2, 2, 2], [2, 4, 1, 1, 2, 2, 3]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r
        if i == 0:
            top = self._list[0][:]
            bot = self._list[1][:]
            if len(top) == 0 or top[-1] != r or bot[-1] != self.parent()._n-r:
                return None
            top.pop()
            bot.pop()
        elif i == r:
            top = self._list[0][:]
            bot = self._list[1][:]
            top.insert(0, 1)
            bot.insert(0, 1)
        elif i < r:
            top = self._list[0][:]
            bot = self._list[1][:]
            k = len(bot)
            if k == 0:
                return None
            pairring_marker = None # Position where we're doing a pairring comparison
            pos = None # Position of unpairred biletter
            step_pos = None # Where the step from i to i+1 occurs
            for j in reversed(range(k)):
                if top[j] == i+1:
                    if pairring_marker is None:
                        pairring_marker = j
                elif top[j] == i:
                    if step_pos is None:
                        step_pos = j
                        if pairring_marker is None:
                            pos = j
                            break
                    while pairring_marker > step_pos and bot[pairring_marker] >= bot[j]:
                        pairring_marker -= 1
                    if pairring_marker == step_pos:
                        pos = j
                        break
                    else:
                        pairring_marker -= 1
                elif top[j] < i: # We've found no unpairred braces to flip
                    return None
            if pos is None:
                return None
            top[step_pos] += 1
            # Shuffle the entries along
            val = bot[pos]
            pos += 1
            while pos < k and (pos <= step_pos or (bot[pos] < val and top[pos] == i+1)):
                bot[pos-1] = bot[pos]
                pos += 1
            bot[pos-1] = val
        else:
            i -= r
            bot = self._list[1][:]
            if len(bot) == 0:
                return None
            pos = None
            count = 0
            for j in range(len(bot)):
                if bot[j] == i:
                    if count == 0:
                        pos = j
                    else:
                        count -= 1
                elif bot[j] == i+1:
                    count += 1
            if pos is None:
                return None
            bot[pos] += 1
            top = self._list[0][:]
        if self.parent()._check_bounded([top, bot]):
            return self.__class__(self.parent(), top, bot)
        return None

    @cached_method
    def to_matrix(self):
        """
        Return ``self`` as a matrix.
        """
        return to_matrix(*self._list)

    @cached_method
    def to_bitableaux(self):
        """
        Return ``self`` as a bitableaux.
        """
        return RSK(*self._list)

class CrystalOfBiwords(CrystalFromRSK):
    """
    Crystal of biwords.

    INPUT:

    - ``n`` -- Type `A_{n-1}^{(1)}`
    - ``r`` -- The number of rows
    - ``s`` -- (Default: ``None``) If this is not ``None``, then this is
      the maximal row/column word length
    """
    def __init__(self, n, r, s=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: CB = CrystalOfBiwords(5, 2)
            sage: TestSuite(CB).run()
        """
        CrystalFromRSK.__init__(self, n, r, s)
        self.module_generators = [self.element_class(self, [], [])]

    def _element_constructor_(self, top, bot):
        """
        Construct an element of ``self`` from ``top`` and ``bot``.
        """
        return self.element_class(self, top, bot)

    Element = CrystalOfBiwordsElement

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CrystalOfBiwords(5, 2)
            Crystal of biwords in the alphabets [2] and [3].
        """
        return "Crystal of biwords in the alphabets [%s] and [%s]"%(self._r, self._n-self._r)

    def _check_bounded(self, biword):
        """
        Check if the longest weakly decreasing subsequence of ``biword``
        is bounded by `s`.
        """
        if self._s is not None:
            r = []
            for x in biword[1]:
                if max(r+[0]) > x:
                    y = min(filter(lambda z: z > x, r))
                    r[r.index(y)] = x
                else:
                    r.append(x)
            if len(r) > self._s:
                return False
        return True

######################
# Bitableaux crystal #
######################

class CrystalOfBitableauxElement(CombinatorialObject, Element):
    """
    An element in a crystal of biwords.
    """
    def __init__(self, parent, left, right):
        """
        Initialize ``self``.
        """
        CombinatorialObject.__init__(self, [left, right])
        Element.__init__(self, parent)

    def _latex_(self):
        """
        Return a `\LaTeX` representation of ``self``.

        EXAMPLES::

            sage: C = CrystalOfBitableaux(8, 4)
            sage: latex(C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3]))
        """
        return self.to_tableaux()._latex_()

    def _repr_diagram(self):
        """
        Return a string representation of ``self`` as an array.

        EXAMPLES::

            sage: C = CrystalOfBitableaux(8, 4)
            sage: print C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])._repr_diagram()
        """
        return self.to_tableaux()._repr_diagram()

    def pp(self):
        """
        Pretty print ``self``.

        EXAMPLES::

            sage: C = CrystalOfBitableaux(8, 4)
            sage: C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3]).pp()
        """
        print self._repr_diagram()

    def e(self, i):
        """
        Return the action of `e_i` on ``self``.

        EXAMPLES::

            sage: C = CrystalOfBitableaux(8, 4)
            sage: elt = C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])
            sage: elt.e(1)
            [[1, 1, 1, 1, 2, 2, 2], [1, 2, 2, 4, 1, 2, 3]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r

        if i == 0:
            if len(self._list[0]) == 0:
                left = [self.parent()._n-r]
                right = [r]
            else:
                nr = self.parent()._n-r
                cols = [self.parent()._to_column_list(l) for l in self._list]
                num_plus = 0
                pos = None
                for i in reversed(range(len(cols[0]))):
                    if cols[0][i][0] == nr:
                        if cols[1][i][0] == r:
                            num_plus += 1
                    elif cols[1][i][0] < r:
                        if num_plus == 0:
                            pos = i
                        else:
                            num_plus -= 1
                if pos is None: # We did not find an unmatched minus in the tableaux
                    cols[0].append([nr])
                    cols[1].append([r])
                else:
                    cols[0][pos].insert(0, nr)
                    cols[1][pos].insert(0, r)
                left = flatten(cols[0])
                right = flatten(cols[1])
        elif i == r:
            if len(self._list[0]) == 0:
                return None
            return None # Placeholder
        else:
            if len(self._list[0]) == 0:
                return None
            if i < r:
                edit = self._list[1][:]
                i_less_r = True
            else:
                edit = self._list[0][:]
                i_less_r = False
                i -= r
            pos = None
            count = 0
            for j in reversed(range(len(edit))):
                if edit[j] == i+1:
                    if count == 0:
                        pos = j
                    else:
                        count -= 1
                elif edit[j] == i:
                    count += 1
            if pos is None:
                return None
            edit[pos] -= 1
            if i_less_r:
                left = self._list[0][:]
                right = edit
            else:
                left = edit
                right = self._list[1][:]
        if self.parent()._check_bounded([left, right]):
            return self.__class__(self.parent(), left, right)
        return None

    def f(self, i):
        """
        Return the action of `f_i` on ``self``.

        EXAMPLES::

            sage: C = CrystalOfBitableaux(8, 4)
            sage: elt = C([1,1,1,2,2,2,2],[1,2,4,1,2,2,3])
            sage: elt.f(1)
            [[1, 1, 2, 2, 2, 2, 2], [2, 4, 1, 1, 2, 2, 3]]
        """
        if i not in self.index_set():
            raise ValueError("i must be in the index set")
        r = self.parent()._r
        if i == 0:
            if len(self._list[0]) == 0:
                return None
            nr = self.parent()._n - r
            cols = [self.parent()._to_column_list(l) for l in self._list]
            pos = None
            num_minus = 0
            for i in range(len(cols[0])):
                if cols[0][i][0] == nr:
                    if cols[1][i][0] == r:
                        if num_minus == 0:
                            pos = i
                        else:
                            num_minus -= 1
                elif cols[1][i][0] < r:
                    num_minus += 1
            if pos is None: # We did not find an unmatched plus in the tableaux
                return None
            cols[0][pos].pop(0)
            cols[1][pos].pop(0)
            left = flatten(cols[0])
            right = flatten(cols[1])
        elif i == r:
            return None # Placeholder
        else:
            if len(self[0]) == 0:
                return None
            if i < r:
                edit = self[1][:]
                i_less_r = True
            else:
                edit = self[0][:]
                i_less_r = False
                i -= r
            pos = None
            count = 0
            for j in range(len(edit)):
                if edit[j] == i:
                    if count == 0:
                        pos = j
                    else:
                        count -= 1
                elif edit[j] == i+1:
                    count += 1
            if pos is None:
                return None
            edit[pos] += 1
            if i_less_r:
                left = self._list[0][:]
                right = edit
            else:
                left = edit
                right = self._list[1][:]
        if self.parent()._check_bounded([left, right]):
            return self.__class__(self.parent(), left, right)
        return None

    @cached_method
    def to_biword(self):
        """
        Return ``self`` as a biword.
        """
        t = self.to_tableaux()
        return RSK_inverse(t[0], t[1])

    @cached_method
    def to_matrix(self):
        """
        Return ``self`` as a matrix.
        """
        t = self.to_tableaux()
        return RSK_inverse(t[0], t[1], output='matrix')

    @cached_method
    def to_tableaux(self):
        """
        Return ``self`` as a tableau tuple.
        """
        from sage.combinat.tableau_tuple import TableauTuple
        ret = []
        for l in self._list:
            if len(l) == 0:
                ret.append(Tableau([]))
                continue
            tab = self.parent()._to_column_list(l)
            for x in tab:
                x.reverse()
            ret.append(Tableau(tab).conjugate())
        return TableauTuple(ret)

class CrystalOfBitableaux(CrystalFromRSK):
    """
    Crystal of bitableaux.

    INPUT:

    - ``n`` -- Type `A_{n-1}^{(1)}`
    - ``r`` -- The number of rows
    - ``s`` -- (Default: ``None``) If this is not ``None``, then this is
      the maximal row/column word length
    """
    def __init__(self, n, r, s=None):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: CB = CrystalOfBitableaux(5, 2)
            sage: TestSuite(CB).run()
        """
        CrystalFromRSK.__init__(self, n, r, s)
        self.module_generators = [self.element_class(self, [], [])]

    def _element_constructor_(self, left, right):
        """
        Construct an element of ``self`` from ``left`` and ``right``.
        """
        return self.element_class(self, left, right)

    Element = CrystalOfBitableauxElement

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: CrystalOfBitableaux(5, 2)
            Crystal of bitableaux in the alphabets [3] and [2].
        """
        return "Crystal of bitableaux in the alphabets [%s] and [%s]"%(self._n-self._r, self._r)

    def _to_column_list(self, l):
        """
        Helper function that returns ``l`` as a list of weakly-increasing
        columns.
        """
        tab = [ [l[0]] ]
        for i in range(1, len(l)):
            if l[i-1] < l[i] or (l[i-1] != 0 and l[i-1] == l[i]):
                tab.append([l[i]])
            else:
                tab[-1].append(l[i])
        return tab

    def _check_bounded(self, bitab):
        """
        Check if the number of columns is bounded by `s`.
        """
        if self._s is not None:
            if len(bitab[0]) == 0:
                return True
            return len(self._to_column_list(bitab[0])) <= self._s
        return True

