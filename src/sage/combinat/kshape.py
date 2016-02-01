r"""
k-shapes

The k-interior of a partition is the subpartition whose cells have
hook_length sizes greater than k.  The k-boundary of a partition is the skew
shape given by the partition minus its k-interior.  The row
(resp. column) shape of a partition is the composition whose i-th part
is the number of cells in the i-th row (resp. column) of the
k-boundary. A partition is a k-shape if both its row shape and column
shape are partitions.

EXAMPLES:

We construct the 3-shape from the partition [4,3,2,2,1,1]::

    sage: x = KShape([4,3,2,2,1,1], 3)

    sage: x.outer_shape
    [4, 3, 2, 2, 1, 1]

The k-interior is called inner_shape::

    sage: x.inner_shape
    [2, 2, 1, 1]

    sage: x.row_shape
    [2, 1, 1, 1, 1, 1]

    sage: x.column_shape
    [2, 2, 2, 1]

AUTHORS:

 - Mark Shimozono (2008)    - initial revision
 - Florent Hivert (2010-02) - Cleanup + generation
 - Olivier Mallet (2010-02) - methods for adding or removing a k-rectangle, irreducibility test
"""

#*****************************************************************************
#       Copyright (C) 2008 Mark Shimozono <mshimo@vt.edu>
#                     2008 Mike Hansen <mhansen@gmail.com>
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

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.structure.element import Element
import sage.combinat.skew_partition
from sage.rings.integer import Integer
from sage.combinat.partition import Partition, Partitions
from copy import copy

from sage.combinat.tools import transitive_ideal
from sage.misc.all import cached_function, cached_method

from sage.misc.latex import latex

def partition_part(p, i):
    if i >= len(p):
        return 0
    else:
        return p[i]



def diagonal_index(x):
    """
        Given a cell x=[r,c], return its diagonal index c-r.
    """
    return x[1] - x[0]

def set_2D_print(b):
    """
    EXAMPLES::

       sage: from sage.combinat.kshape import set_2D_print
       sage: set_2D_print(True)
       sage: KShape([4,3,1,1], 4)
       .###
       .##
       #
       #
       sage: set_2D_print(False)
   """
    Partition._has_2D_print = b
    sage.combinat.skew_partition.SkewPartition._has_2D_print = b
    KShape._has_2D_print = b
    sage.combinat.dyck_word.DyckWord._has_2D_print = b

@cached_function
def MatrixSpace_cached(r, c):
    """
    """
    from sage.matrix.all import MatrixSpace
    from sage.rings.all import ZZ
    return MatrixSpace(ZZ, r, c, sparse = True)

class SignedNum(int):
    """
    EXAMPLES::

        sage: from sage.combinat.kshape import SignedNum
        sage: SignedNum(4)
        +4
        sage: SignedNum(-2)
        -2
        sage: SignedNum(0)
        0
        sage: SignedNum(0) < 2
        True
        sage: SignedNum(2) + 3
        5
    """
    def __repr__(self):
        if self > 0: return "+"+"%i"%self
        else: return "%i"%self




################################################
# k-shape class
################################################
class KShape(SageObject):
    def __init__(self, outer_shape, k, inner_shape=None, row_shape=None, column_shape=None):
        """
        Construct a k-shape from its various data.

        INPUT:

         - ``outer_shape`` - a partition (either as a list of integers or a Sage
                        :class:`~sage.combinat.partition.Partition`)
         - ``k`` - an integer
         - ``inner_shape`` - a partition (list or Partition) (optional)
         - ``row_shape`` - a list of integers (optional)
         - ``column_shape`` - a list of integers (optional)

        .. warning::

            There is no check that the various data are consistent. You can
            check it using :meth:`check_k_shape`.


        EXAMPLES::

            sage: k1=KShape([5,4,3,2,2,1,1], 3); k1
            [[5, 4, 3, 2, 2, 1, 1], [3, 2, 2, 1, 1]] with k=3
            sage: k2=KShape([5,4,3,2,2,1,1], 3, inner_shape=[3,2,2,1,1]); k1==k2
            True
            sage: k3=KShape([5,4,3,2,2,1,1], 3, inner_shape=[3,2,2,1,1],
            ...             row_shape=[2,2,1,1,1,1,1], column_shape=[2,2,2,2,1]); k1==k3
            True
        """
        # the partition defining the k-shape. It is a Partition object.
        self.outer_shape = Partition(outer_shape)
        # the k in k-shape:
        self.k = k
        # the k-interior, also a Partition object
        if inner_shape is None or row_shape is None or column_shape is None:
            boundary = self.outer_shape.k_boundary(k)
        if inner_shape is None:
            inner_shape = boundary.inner()
        if row_shape is None:
            row_shape = boundary.row_lengths()
        if column_shape is None:
            column_shape = boundary.column_lengths()
        self.inner_shape = Partition(inner_shape)
        # the row shape, as a list of integers
        self.row_shape = row_shape
        # the column shape, as a list of integers
        self.column_shape = column_shape
    _has_2D_print = False

    def check_k_shape(self):
        """
        Check the consistency of the attributes of ``self``.

        EXAMPLES::

            sage: KShape([5,4,3,2,2,1,1], 3, inner_shape=[3,2,2,1]).check_k_shape()
            Traceback (most recent call last):
            ...
            AssertionError: Wrong inner shape
            sage: KShape([3,3,3,3,3,3,3,3,1,1], 4).check_k_shape()
            Traceback (most recent call last):
            ...
            AssertionError: Row shape is not a partition
        """
        boundary = self.outer_shape.k_boundary(self.k)
        assert self.inner_shape == boundary.inner(), "Wrong inner shape"
        assert self.row_shape == boundary.row_lengths(), "Wrong row shape"
        assert self.row_shape in Partitions(), "Row shape is not a partition"
        assert self.column_shape == boundary.column_lengths(), "Wrong column shape"
        assert self.column_shape in Partitions(), "Column shape is not a partition"

    def _repr_(self):
        """
        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3); x
            [[5, 4, 3, 2, 2, 1, 1], [3, 2, 2, 1, 1]] with k=3

        """
        if self._has_2D_print:
            return self.to_skew_partition().ferrers_diagram("#",".")
        else:
            return "[%r, %r] with k=%r"%(self.outer_shape, self.inner_shape, self.k)

    def __iter__(self):
        """
        EXAMPLES::

            sage: list(KShape([5,4,3,2,2,1,1], 3))  # indirect doctest
            [5, 4, 3, 2, 2, 1, 1]
        """
        return iter(self.outer_shape)

    def __eq__(self,other):
        """
        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3); x==x
            True
            sage: x == Partition([5,4,3,2,2,1,1]).k_shape(3)
            True
            sage: x == Partition([4,3,2,2,1,1]).k_shape(3)
            False
            sage: x = Partition([1]); x.k_shape(2)==x.k_shape(3)
            False
        """
        return (isinstance(other, type(self)) and self.k==other.k
                and self.outer_shape==other.outer_shape)

    def to_skew_partition(self):
        """
        Returns this k-shape as a skew partition.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3); x
            [[5, 4, 3, 2, 2, 1, 1], [3, 2, 2, 1, 1]] with k=3
            sage: x.to_skew_partition()
            [[5, 4, 3, 2, 2, 1, 1], [3, 2, 2, 1, 1]]
        """
        # bypass the tests for speed
        return sage.combinat.skew_partition.SkewPartition(
            [self.outer_shape, self.inner_shape])
    
    def _latex_pspicture(self):
        """
        EXAMPLES::

            sage: print KShape([5,4,2,2,1],4)._latex_pspicture()
            \begin{pspicture}(0,0)(5,5)
               \multirput(2,0)(1,0){3}{\psframe(0,0)(1,1)}
               \multirput(2,1)(1,0){2}{\psframe(0,0)(1,1)}
               \multirput(0,2)(1,0){2}{\psframe(0,0)(1,1)}
               \multirput(0,3)(1,0){2}{\psframe(0,0)(1,1)}
               \multirput(0,4)(1,0){1}{\psframe(0,0)(1,1)}
            \end{pspicture}
        """
        out,inn = self.to_skew_partition()
        if not out:
            return "\\emptyset"
        res = "\\begin{pspicture}(0,0)(%s,%s)\n"%(out[0],len(out))
        for i in range(len(inn)):
            res += "   \multirput(%s,%s)(1,0){%s}{\\psframe(0,0)(1,1)}\n"%(
                inn[i], i, out[i]-inn[i])
        for i in range(len(inn), len(out)):
            res += "   \multirput(%s,%s)(1,0){%s}{\\psframe(0,0)(1,1)}\n"%(
                0, i, out[i])
        res += "\\end{pspicture}"
        return res

    color_dict = {1:"white", 2:"gray", 3:"green",
                  8:"blue", 9:"red"}

    tikz_scale = "0.2"

    def _latex_tikz(self, color=True):
        """
        EXAMPLES::

            \vcenter{\hbox{$\begin{tikzpicture}[scale=0.2]
              \draw(0,2) rectangle (1,3);
              \draw(0,3) rectangle (1,4);
              \draw(0,4) rectangle (1,5);
              \draw(1,2) rectangle (2,3);
              \draw(1,3) rectangle (2,4);
              \draw(2,0) rectangle (3,1);
              \draw(2,1) rectangle (3,2);
              \draw(3,0) rectangle (4,1);
              \draw(3,1) rectangle (4,2);
              \draw(4,0) rectangle (5,1);
            \end{tikzpicture}$}}
        """
        out,inn = self.to_skew_partition()
        if not out:
            return "\\emptyset"
        mat = self.matrix()
        res="\\vcenter{\\hbox{$\\begin{tikzpicture}[scale="+self.tikz_scale+"]\n"
        for (i,j),col in mat.dict().items():
            res += "  \\draw(%s,%s)%s rectangle (%s,%s);\n"%(
                j-1,i-1,
                "[fill = %s]"%(self.color_dict[col]) if color else "",
                j,i)
        res += "\\end{tikzpicture}$}}"
        return res

    latex.add_package_to_preamble_if_available("tikz")
    _latex_ = _latex_tikz


    @cached_method
    def matrix(sh):
        """
        Transform ``self`` into a matrix

        EXAMPLES::

            sage: from sage.combinat.kshape import *
            sage: v = ((0, 1, 0, 0), (1, 0, 1))
            sage: l = IrreducibleKShapes_from_kkm1_vectors(v).list()
            sage: [sh.matrix() for sh in l]
            [
            [0 0 0 0 0 0 0 0 0 0]                       [0 0 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 8 2 2 2 0]  [0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 8 2 2 2 0]
            [0 0 0 0 1 2 1 0 0 0]  [0 0 0 0 8 2 2 2 0]  [0 0 0 0 0 2 1 1 0 0]
            [0 0 9 2 2 0 0 0 0 0]  [0 0 9 2 3 0 0 0 0]  [0 0 9 2 2 0 0 0 0 0]
            [0 0 2 1 0 0 0 0 0 0]  [0 0 2 1 0 0 0 0 0]  [0 0 2 1 1 0 0 0 0 0]
            [0 0 2 1 0 0 0 0 0 0]  [0 0 2 1 0 0 0 0 0]  [0 0 2 1 1 0 0 0 0 0]
            [0 8 3 0 0 0 0 0 0 0]  [0 8 3 0 0 0 0 0 0]  [0 8 3 0 0 0 0 0 0 0]
            [0 2 0 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0 0]
            [0 2 0 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0 0]
            [0 2 0 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0], [0 0 0 0 0 0 0 0 0], [0 0 0 0 0 0 0 0 0 0],
            <BLANKLINE>
            [0 0 0 0 0 0 0 0 0]  [0 0 0 0 0 0 0 0 0 0]                   
            [0 0 0 0 8 2 2 2 0]  [0 0 0 0 0 8 2 2 2 0]  [0 0 0 0 0 0 0 0]
            [0 0 9 2 3 0 0 0 0]  [0 0 0 1 1 2 0 0 0 0]  [0 0 0 8 2 2 2 0]
            [0 0 2 1 0 0 0 0 0]  [0 0 9 2 2 0 0 0 0 0]  [0 0 9 3 2 0 0 0]
            [0 8 3 0 0 0 0 0 0]  [0 8 3 0 0 0 0 0 0 0]  [0 8 3 0 0 0 0 0]
            [0 2 2 0 0 0 0 0 0]  [0 2 2 0 0 0 0 0 0 0]  [0 2 2 0 0 0 0 0]
            [0 2 0 0 0 0 0 0 0]  [0 2 2 0 0 0 0 0 0 0]  [0 2 2 0 0 0 0 0]
            [0 2 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0 0 0]  [0 2 0 0 0 0 0 0]
            [0 0 0 0 0 0 0 0 0], [0 0 0 0 0 0 0 0 0 0], [0 0 0 0 0 0 0 0]
            ]
        """
        k = sh.k
        rs = len(sh.row_shape)+2
        cs = len(sh.column_shape)+2
        empty  = 0
        mat = MatrixSpace_cached(rs, cs)([empty]*rs*cs)
        skp = sh.to_skew_partition()
        for (r, c) in skp.cells():
            mat[r+1, c+1] = 1             # free cell
        for (r, c) in skp.inner_corners():
            arm = skp.outer().arm_length(r, c)
            leg = skp.outer().leg_length(r, c)
            if 1+arm+leg == k-1 or 1+arm+leg == k:
                if 1+arm+leg == k-1:
                    mat[r+1, c+1] = 8     # k-1 corner
                else:
                    mat[r+1, c+1] = 9     # k   corner
                i = c+2
                while mat[r+1, i] != empty:
                    mat[r+1, i] += 1       # hook cell
                    i+=1
                i = r+2
                while mat[i, c+1] != empty:
                    mat[i, c+1] += 1       # hook cell
                    i+=1
        mat.set_immutable()
        return mat

    def diagram(self):
        """
        Returns a string which represents the Ferrers diagram
        of this k-shape.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: print x.diagram()
               ##
              ##
              #
             #
             #
            #
            #

        """
        return self.to_skew_partition().diagram()

    def conjugate(self):
        r"""
        Return a new conjugate k-shape object.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3); x
            [[5, 4, 3, 2, 2, 1, 1], [3, 2, 2, 1, 1]] with k=3
            sage: x.conjugate()
            [[7, 5, 3, 2, 1], [5, 3, 1]] with k=3

        """
        return KShape(self.outer_shape.conjugate(), self.k,
                            inner_shape=self.inner_shape.conjugate(),
                            row_shape=self.column_shape, column_shape=self.row_shape)

    from sage.misc.superseded import deprecated_function_alias
    transpose = deprecated_function_alias(10000, conjugate) # kshape-om.patch, 2010-11-29; insert proper ticket number when available

    def bottom_boundary(self, c):
        r"""
        Give the row index of the bottom cell in column c of the
        k-boundary of a k-shape.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: print x.diagram()
               ##
              ##
              #
             #
             #
            #
            #
            sage: x.bottom_boundary(3)
            0
            sage: x.bottom_boundary(7)
            0
            sage: x.bottom_boundary(0)
            5
        """
        if (self.inner_shape.length() == 0 or
            c >= self.inner_shape[0]):
            return 0
        return self.inner_shape.conjugate()[c]

    def left_boundary(self, r):
        r"""
        Give the column index of the leftmost cell in row r of the
        k-boundary.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: print x.diagram()
               ##
              ##
              #
             #
             #
            #
            #
            sage: x.left_boundary(0)
            3
            sage: x.left_boundary(2)
            2
            sage: x.left_boundary(8)
            0
        """
        return self.inner_shape[r] if r < self.inner_shape.length() else 0

    def is_row_string_starter(self, box):
        r"""
        Check whether the length of the hook of the leftmost cell in
        row r of the k-boundary, equals k.
        """
        r, c = box
        if r >= self.outer_shape.length():
            return False
        c = self.inner_shape[r] if r < self.inner_shape.length() else 0
        return self.outer_shape.hook_length(r, c) == self.k

    def is_row_string_ender(self, box):
        r"""
        Check whether the length of the hook of the (French)
        bottommost cell in column c of the k-boundary, is strictly
        less than k.
        """
        r, c = box
        if r == 0:
            return True
        if self.inner_shape.length() == 0 or c >= self.inner_shape[0]:
            r = 0
        else:
            r = self.inner_shape.conjugate()[c]
        return self.outer_shape.hook_length(r, c) < self.k

    def k_strings_max(self):
        """
        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: x.k_strings_max()
            [[(7, 0), (5, 1), (3, 2), (1, 4)], [(2, 3), (0, 5)]]
        """
        corners = self.outer_shape.outside_corners()
        string_list = []
        for r,c in reversed(corners):
            diagonal_corner = c - r
            for string in string_list:
                rs, cs = string[-1]
                diagonal_s = cs - rs
                if diagonal_corner - self.k - diagonal_s in range(2):
                    string.append((r,c))
                    break
            else:
                string_list.append([(r,c)])
        return string_list


#    def is_row_move_degenerate(self, move):
#    r"""
#    Is a row move on a k-shape degenerate or not?
#    Does the column of the rightmost removed cell, lie just to the left of
#        the column of the leftmost added cell?
#    """


    def row_strings(self):
        r"""
        Returns the list of addable row-type strings.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: x.row_strings()
            [[(1, 4)], [(0, 5)]]

        """
        string_list = self.k_strings_max()
        row_strings = []
        for string in string_list:
            for j, start_box in enumerate(string):
                if not self.is_row_string_starter(start_box):
                    continue
                for m, end_box in enumerate(string[j:]):
                    if not self.is_row_string_ender(end_box):
                        continue
                    row_strings.append(string[j:j+m+1])
        return row_strings

    def is_k_shape(self):
        """
        Checks if the row shape and column shape of ``self`` are partitions.

        EXAMPLES::

            sage: KShape([5,4,2,2,1],4).is_k_shape()
            True
            sage: KShape([5,4,3,3],4).is_k_shape()
            False
        """
        k_boundary = self.to_skew_partition()
        return (k_boundary.row_lengths() in Partitions() and
                k_boundary.column_lengths() in Partitions())


    def is_row_string(self, string):
        r"""
        Returns True if ``string`` is an addable row-type string.  The row
        indices are assumed to be distinct.

        EXAMPLES::

            sage: x = Partition([5,4,3,2,2,1,1]).k_shape(3)
            sage: [x.is_row_string(a) for a in x.row_strings()]
            [True, True]
        """
        for r, c in string:
            if r != 0 and not (r > 0 and self.outer_shape[r] < self.outer_shape[r-1]):
                return False

        return (self.is_row_string_starter(string[0])
                and self.is_row_string_ender(string[-1]))

    def add_string(self, string):
        return self.adjoin_cells(string)

    def adjoin_cells(self, cells):
        """
        cells is a list of pairs of nonnegative integers.  p is a list
        of integers.  For each [r,c] in cells, increment the part r of
        p.  If r >= the current length of the list, adjoin 1 to the
        list.
        """
        q = self.outer_shape[:]
        for r,c in cells:
            if r >= len(q):
                q.append(1)
            else:
                q[r] += 1
        return KShape(q, self.k)

    def row_moves_with_first_string(self, first_string):
        r"""
        Assuming ``first_string`` is an addable row-type string, generate
        all the addable row moves with it as first string.
        """
        moves = []
        curr_shape = self.add_string(first_string)
        if curr_shape.is_k_shape():
            moves.append(copy(curr_shape))
        # for each rank, try to make an addable move of that rank
        if len(first_string) > 1:
            rank_bound = first_string[1][1] - first_string[0][1]
        else:
            rank_bound = self.k - 1
        outer_conj = list(self.outer_shape.conjugate())
        first_cell = first_string[0]
        for rk in range(1,rank_bound):
            translation_vector = [partition_part(outer_conj,first_cell[1]+rk) - first_cell[0],rk]
            if diagonal_index(translation_vector) - diagonal_index(first_cell) >= self.k:
                break
            string = [[translation_vector[0]+x[0],translation_vector[1]+x[1]] for x in first_string]
            if not curr_shape.is_row_string(string):
                break
            curr_shape = curr_shape.add_string(string)
            if curr_shape.is_k_shape():
                moves.append(copy(curr_shape))
        return moves

    def row_moves(self):
        r"""
        Generate all addable row moves.

        TESTS::

            sage: x = KShape([4,3,2,2,1,1], 3)
            sage: x.row_moves()
            [[[5, 3, 2, 2, 1, 1], [3, 2, 1, 1]] with k=3]
        """
        return [i for s in self.row_strings()
                for i in self.row_moves_with_first_string(s)]

    def column_moves(self):
        r"""
        Generate all addable column moves.
        """
        kshtr = self.conjugate()
        tmoves = kshtr.row_moves()
        return [x.conjugate() for x in tmoves]

    def all_moves(self):
        return self.row_moves() + self.column_moves()

    def transitive_moves(self):
        r"""
        Generate all k-shapes that can be gotten from the given k-shape
        by a sequence of moves.
        """
        return transitive_ideal(lambda x: x.all_moves(), self)

    def mixed_equivalences(self):
        r"""
        Generate the list of triples (rm, cm, fm) of mixed equivalences
        starting with ``self``, where ``self`` ---> rm is a row move,
        ``self`` ---> cm is a column move, and fm is the k-shape
        that completes the commutative diamond.
        """
        rms = self.row_moves()
        cms = self.column_moves()
        # for each row move and column move starting at self
        # compute the mixed equivalences they generate, if any
        for rm in rms:
            for cm in cms:
                pass

    def row_equivalences(self):
        r"""
        Give a list of triples (rm, Rm, fm) giving mixed equivalences
        starting with ``self``, where there are row moves
        ``self`` ---> rm, ``self`` ---> Rm and fm completes the commutative diamond.
        """
        pass


    def row_plus_1(self, inverse=False):
        """
        Returns a (k+1)-shape having one more cell in each row than ``self`` 
        (or a (k-1)-shape having one less cell in each row than ``self``, 
        if ``inverse`` is set to ``True``).

        EXAMPLES::

            sage: KShape([8,6,4,3,3,1,1,1], 5).row_plus_1()
            [[11, 9, 6, 5, 5, 2, 2, 2], [6, 5, 2, 2, 2]] with k=6
        """
        outer = self.outer_shape[:]
        outer_orig = outer[:]
        inner = self.inner_shape[:]
        inner = inner+[0]*(len(outer)-len(inner))
        cut = 0
        band = 0
        for i in range(-1, -len(outer)-1, -1):
            if inner[i] > cut: # new band
                band += 1
                cut = outer_orig[i+1]
            if inverse:
                outer[i] -= band+1
                inner[i] -= band
            else:
                outer[i] += band+1
                inner[i] += band
        return KShape(outer, self.k + (-1 if inverse else 1))

    def add_k_rectangle(self,width,height):
        r"""
        Adds a rectangle to a k-shape.
        The largest hook length of the rectangle (``width`` + ``height`` - 1) must be k or k-1.

        INPUT:

         - ``width`` - an integer
         - ``height`` - an integer

        EXAMPLES::

            sage: KShape([8,6,4,3,3,1,1,1], 5).add_k_rectangle(3, 2)
            [[11, 9, 7, 6, 4, 3, 3, 1, 1, 1], [7, 6, 4, 3, 1, 1, 1]] with k=5
            sage: KShape([2,1,1], 4).add_k_rectangle(3, 1)
            [[5, 2, 1, 1], [2]] with k=4
            sage: KShape([5,2,1,1], 4).add_k_rectangle(1,4)
            [[6, 3, 2, 2, 1, 1, 1, 1], [3, 1, 1, 1]] with k=4
            sage: KShape([3,3],4).add_k_rectangle(2,3)
            [[5, 5, 2, 2, 2], [2, 2]] with k=4
            sage: KShape([5,4,3,3], 3).add_k_rectangle(3, 2)
            Traceback (most recent call last):
            ...
            AssertionError: Invalid rectangle size
        """
        assert width+height-1==self.k or width+height-1==self.k-1,"Invalid rectangle size"
        # construct the row shape and column shape of the resulting k-shape
        rows=list(self.row_shape)
        columns=list(self.column_shape)
        index_rows=0
        while index_rows<len(rows) and rows[index_rows]>=width:
            index_rows+=1
        index_columns=0
        while index_columns<len(columns) and columns[index_columns]>=height:
            index_columns+=1
        rows[index_rows:index_rows]=height*[width]
        columns[index_columns:index_columns]=width*[height]
        # construct the resulting k-shape from its row shape and column shape
        outer=[rows[-1]]
        y=0
        for x in rows[-2::-1]:
            while y<outer[0] and Partition(outer).conjugate()[y]+1>columns[y]:
                y+=1
            outer.insert(0,x+y)
        return KShape(outer,self.k)

    def rectangles(self, sorted=False):
        """
        Returns the list of rectangles that it is possible to add to ``self``, 
        i.e. the list of all k- or (k-1)-rectangles.

        EXAMPLES::

            sage: k = KShape([4,3,1,1], 4)
            sage: k.rectangles()
            [(1, 4), (2, 3), (3, 2), (4, 1), (1, 3), (2, 2), (3, 1)]
            sage: KShape([], 0).rectangles()
            []
            sage: KShape([], 1).rectangles()
            [(1, 1)]
            sage: KShape([], 2).rectangles()
            [(1, 2), (2, 1), (1, 1)]
            sage: KShape([], 3).rectangles()
            [(1, 3), (2, 2), (3, 1), (1, 2), (2, 1)]

            sage: k.rectangles(True)
            [(1, 4), (1, 3), (2, 3), (2, 2), (3, 2), (3, 1), (4, 1)]
        """
        if sorted:
            res = []
            for w in range(1,self.k):
                res.append((w,self.k+1-w))
                res.append((w,self.k-w))
            res.append((self.k, 1))
            return res
        else:
            return [(w,self.k+1-w) for w in range(1,self.k+1)] + \
                   [(w,self.k-w) for w in range(1,self.k)]

    def remove_k_rectangle(self,width,height):
        r"""
        Removes a rectangle from a k-shape.
        The largest hook length of the rectangle (``width`` + ``height`` - 1) must be k or k-1.

        INPUT:

         - ``width`` - an integer
         - ``height`` - an integer

        EXAMPLES::

            sage: KShape([11, 9, 7, 6, 4, 3, 3, 1, 1, 1], 5).remove_k_rectangle(3, 2)
            [[8, 6, 4, 3, 3, 1, 1, 1], [4, 3, 1, 1, 1]] with k=5
        """
        assert width+height-1==self.k or width+height-1==self.k-1,"Invalid rectangle size"
        # construct the row shape and column shape of the resulting k-shape
        rows=list(self.row_shape)
        columns=list(self.column_shape)
        for i in range(height):
            try:
                rows.remove(width)
            except ValueError:
                raise ValueError, "Unable to remove rectangle"
        for j in range(width):
            try:
                columns.remove(height)
            except ValueError:
                raise ValueError, "Unable to remove rectangle"
        # construct the resulting k-shape from its row shape and column shape
        if not rows:
            return KShape([],self.k)
        outer=[rows[-1]]
        y=0
        for x in rows[-2::-1]:
            while y<outer[0] and Partition(outer).conjugate()[y]+1>columns[y]:
                y+=1
            outer.insert(0,x+y)
        ks=KShape(outer,self.k)
        if ks.row_shape==rows and ks.column_shape==columns:
            return ks
        else:
            raise ValueError, "Result is not a k-shape"

    def inner_frontier(self):
        """
        Returns the cells on the inner frontier of ``self``.

        EXAMPLES::

            sage: ks = KShape([4,3,1,1], 4)
            sage: ks.inner_frontier()
            [(4, 0), (3, 0), (2, 0), (2, 1), (1, 1), (0, 1), (0, 2), (0, 3), (0, 4)]
        """
        res = self.inner_shape.outer_rim()
        conj = self.conjugate()
        for i in range(partition_part(conj.inner_shape, 0)+1,
                       partition_part(conj.outer_shape, 0)+1):
            res.insert(0, (i, 0))
        for j in range(partition_part(self.inner_shape, 0)+1,
                       partition_part(self.outer_shape, 0)+1):
            res.append((0, j))
        return res

    def rc_sorted_frontier(self):
        """
        Compute the row-column sorted frontier of ``self``.

        EXAMPLES::

            sage: ks = KShape([4,3,1,1,1], 4)
            sage: rc = ks.rc_sorted_frontier()
            sage: for rect in ks.rectangles(True):
            ...       print rect, rc[rect]
            (1, 4) [(5, 0), (4, 0), (3, 0), (2, 0)]
            (1, 3) [(5, 0), (4, 0), (3, 0), (2, 0), (2, 1)]
            (2, 3) [(2, 0), (2, 1), (1, 1)]
            (2, 2) [(2, 1), (1, 1)]
            (3, 2) [(1, 1), (0, 1), (0, 2), (0, 3)]
            (3, 1) [(0, 3), (0, 4)]
            (4, 1) [(0, 3), (0, 4)]
        """
        res = dict()
        rects = self.rectangles()
        for i,j in rects:
            res[i,j] = []
        for i,j in self.inner_frontier():
            if i == 0:
                maxr = self.k
            else:
                maxr = partition_part(self.row_shape, i-1)
            minr = partition_part(self.row_shape, i)
            if j == 0:
                maxc = self.k
            else:
                maxc = partition_part(self.column_shape, j-1)
            minc = partition_part(self.column_shape, j)
            for ir in range(minr, maxr+1):
                for jr in range(minc, maxc+1):
                    if (ir, jr) in rects:
                        res[ir,jr].append((i,j))
        return res

    def rc_code(self):
        """
        """
        rcf = self.rc_sorted_frontier()
        return [len(rcf[rect])-1 for rect in self.rectangles(True)]

    def restrict(self, u, v):
        """
        EXAMPLES::

            sage: ks = KShape([10,7,4,2,2,2,1,1,1,1],4)
            sage: ks.restrict(2,3)
            [[3, 1, 1, 1], [1]] with k=4
            sage: ks.restrict(3,1)
            [[8, 5, 2], [5, 2]] with k=4
        """
        if u+v-1 != self.k and u+v-1 != self.k-1:
            raise ValueError, "Input is not a k- or (k-1)-rectangle"
        outer = list(self.outer_shape)
        if outer == []:
            return KShape([],self.k)
        i = 0
        while i < len(self.row_shape) and self.row_shape[i] > u:
            outer.pop(0)
            i = i+1
        j = 0
        while j < len(self.column_shape) and self.column_shape[j] > v:
            for m in range(len(outer)-1,-1,-1):
                if outer[m] == 1:
                    outer.pop(m)
                else:
                    outer[m] = outer[m] - 1
            j = j+1
        return KShape(outer, self.k)

    def restrict_inv(self,u,v):
        """
        Returns the list of irreducible k-shapes s such that s.restrict(u,v) == self.
        """
        if u+v-1 != self.k and u+v-1 != self.k-1:
            raise ValueError, "Input is not a k- or (k-1)-rectangle"
        if self.outer_shape != [] and (self.row_shape[0] > u or self.column_shape[0] > v):
            raise ValueError, "Impossible to obtain this k-shape via restriction to this rectangle"
        return [s for s in IrreducibleKShapes(self.k) if s.restrict(u,v) == self]

    def is_irreducible_def(self):
        r"""
        Checks if this k-shape is irreducible,
        i.e. if it is possible to remove a k-rectangle or a (k-1)-rectangle from it.

        The algorithm is straightforward from the definition.

        EXAMPLES::

            sage: KShape([4,3,1,1,1], 4).is_irreducible_def()
            False
            sage: KShape([5,3,2,1], 4).is_irreducible_def()
            False
            sage: KShape([4,2,1], 4).is_irreducible_def()
            True
            sage: KShape([5,4,2,2,1], 4).is_irreducible_def()
            True
        """
        for r in self.rectangles():
            try:
                self.remove_k_rectangle(*r)
                return False
            except ValueError:
                pass
        return True

    def is_irreducible(self):
        """
        Fast irreducibility test

        EXAMPLE::

            sage: KShape([4,3,1,1,1], 4).is_irreducible()
            False
            sage: KShape([5,3,2,1], 4).is_irreducible()
            False
            sage: KShape([4,2,1], 4).is_irreducible()
            True
            sage: KShape([5,4,2,2,1], 4).is_irreducible()
            True
        """
        rcf = self.rc_sorted_frontier()
        for rect in self.rectangles():
            if len(rcf[rect]) >= sum(rect):
                return False
        return True

    def add_row(self,row=None):
        """
        Returns the list of k-shapes obtained by adding a row to ``self``.

        INPUT:

            ``row`` - an integer

        OUTPUT:

            The list of k-shapes obtained by adding a row of length ``row`` to ``self``, or the list of k-shapes obtained by adding a row of any admissible length to ``self``, if no value is given for ``row``.

        EXAMPLES::

            sage: bla = KShape([], 4)
            sage: bla.add_row()
            [[[1], []] with k=4, [[2], []] with k=4, [[3], []] with k=4, [[4], []] with k=4]
        """
        if row is None:
            res=[]
            min_row = 1 if not self.row_shape else self.row_shape[0]
            for i in range(min_row,self.k+1):
                res += self.add_row(i)
            return res
        if row>self.k or row<partition_part(self.row_shape,0):
            raise ValueError, "Impossible to add a row of this length"
        inner0 = self.inner_shape[0] if self.inner_shape else 0
        res=[]
        for (outer, inner, k, foot, hmax, parts) in children(
            (self.outer_shape[:],
             self.inner_shape[:],
             self.k,
             self.column_shape[inner0:],
             ([self.k]+self.column_shape)[inner0],
             [row])):
            res.append(KShape(outer,k,inner))
        return res

    def foot(self):
        """
        Returns the foot of ``self``.

        The foot is defined as the k-shape obtained by removing the columns
        which don't touch the ground.

        EXAMPLES::

            sage: KShape([3,2,1,1], 4).foot()
            [2, 1]
            sage: KShape([], 4).foot()
            []
            sage: KShape([1, 1, 1], 4).foot()
            [1, 1, 1]
            sage: KShape([3,3,1,1,1], 5).foot()
            [2, 2]
        """
        if self.inner_shape:
            res = Partition(self.outer_shape).conjugate()[self.inner_shape[0]:]
            return Partition(res).conjugate()
        else:
            return Partition(self.outer_shape)

    def h_max(self):
        """
        Returns the `h_{max}` of ``self``.

        The `h_{max}` is the height of the rigthmost column which doesn't touch the
        ground, or k if there isn't any such column.

        EXAMPLES::

            sage: KShape([4,3,1,1], 4).h_max()
            2
            sage: KShape([], 4).h_max()
            4
            sage: KShape([1, 1, 1], 4).h_max()
            4
            sage: KShape([5, 2, 1], 5).h_max()
            2
            sage: KShape([4,2,1,1], 5).h_max()
            3
            sage: KShape([3,3,1,1,1], 5).h_max()
            3
        """
        if self.inner_shape:
            return (self.column_shape[:])[self.inner_shape[0]-1]
        else:
            return self.k

    def left_right_seq(self, row_length=False):
        """
        sage: KShape([3,3,1,1,1], 5).left_right_seq()
        (0, 1, 0, 0, 0)
        sage: KShape([3,3,1,1,1], 5).left_right_seq(True)
        (2, +2, 1, 1, 1)

        sage: all(k.left_right_seq() == k.row_plus_1().left_right_seq()
        ...       for k in IrreducibleKShapes(3))
        True
        sage: all(k.left_right_seq() == k.row_plus_1().left_right_seq()
        ...       for k in IrreducibleKShapes(4))
        True
        sage: all(k.left_right_seq() == k.row_plus_1().left_right_seq()
        ...       for k in IrreducibleKShapes(5))
        True

        """
        height = len(self.outer_shape)
        restr = [KShape(self.outer_shape[i:], self.k)
                 for i in range(height, -1, -1)]
        res = []
        for i in range(height):
            new_part = self.row_shape[height-i-1]
            nextsh = restr[i].add_row(new_part)
            if row_length:
                if len(nextsh) == 1:
                    res.append(new_part)
                else:
                    if restr[i+1] == nextsh[0]:
                        res.append(SignedNum(-new_part))
                    else:
                        res.append(SignedNum(new_part))
            else:
                if len(nextsh) == 1:
                    res.append(0)
                else:
                    if restr[i+1] == nextsh[0]:
                        res.append(-1)
                    else:
                        res.append(+1)
        return tuple(res[::-1])

    def corner_hooks(self, leng=None):
        """
        Returns the list of corners of hook length ``leng``, 
        or the list of all corners if ``leng`` is not specified.

        EXAMPLES::

            sage: k = 4
            sage: sum(x^(k-len(filter(lambda (x,y): x+y-1 == k,
            ...                set(v.corner_hooks()))))
            ...       for v in IrreducibleKShapes(k))
            6*x^4 + 8*x^3 + 3*x^2
        """
        res   = self.to_skew_partition().inner_corners()
        res   = [(self.outer_shape.arm_length(*c)+1,
                  self.outer_shape.leg_length(*c)+1) for c in res]
        if leng is not None:
            res = filter(lambda (x,y): x+y-1 == leng, res)
        return tuple(res)

    def corner_hook_multiplicities(self, leng=None):
        """
        Returns a dictionary which associate to any hook its multiplicity as a
        corner hook of self

        EXAMPLES::

            sage: KShape([4,3,1,1], 4).corner_hook_multiplicities()
            defaultdict(<type 'int'>, {(1, 2): 1, (3, 2): 1})
            sage: KShape([4,3,1,1], 4).corner_hook_multiplicities(4)
            defaultdict(<type 'int'>, {(3, 2): 1})
            sage: KShape([5,2,2,1,1], 5).corner_hook_multiplicities()
            defaultdict(<type 'int'>, {(3, 1): 1, (2, 4): 1})
        """
        from collections import defaultdict
        res  = defaultdict(int)
        for (x,y) in self.corner_hooks(leng):
            res[x,y] +=1
        return res

    def corner_hook_vector(self):
        """
        Returns the corner hook vector

        For a `k`-shape `l` `V(l) := [v_1,...v_k]` with `v_i = 1` if
        there is a corner in the band `(i,k-i+1)` and `0` otherwise.

        EXAMPLES::

            sage: KShape([4,3,1,1], 4).corner_hook_vector()
            (0, 1)
            sage: KShape([5,2,2,1,1], 5).corner_hook_vector()
            (1, 0, 0)
        """
        k = self.k
        s = set(self.corner_hooks(k))
        return tuple(1 if (i,k-i+1) in s else 0 for i in range(2, k))

    def k_plus_one_corner_hook_vector(self):
        """
        """
        k = self.k
        res = [(self.row_shape[i],self.column_shape[j]) for (i,j) in self.inner_shape.corners()]
        s = set(filter(lambda (x,y):x+y+1==k+1,res))
        return tuple(1 if (i,k-i) in s else 0 for i in range(2, k-1))

    def k_plus_one_corner_hook_multiplicities(self):
        """
        EXAMPLES ::

            sage: s=KShape([4,3,1,1], 4)
            sage: s.k_plus_one_corner_hook_multiplicities()
            (1,)
        """
        k = self.k
        res = [(self.row_shape[i],self.column_shape[j]) for (i,j) in self.inner_shape.corners()]
        s = filter(lambda (x,y):x+y+1==k+1,res)
        return tuple(s.count((i,k-i)) for i in range(2, k-1))

    def kk1_mult(self):
        """
        EXAMPLES::

            sage: s=KShape([14,10,9,8,6,6,5,4,3,3,2,2,1,1,1],6)
            sage: s.kk1_mult()
            ((1, 2, 1, 1), (1, 2, 0))
        """
        return (self.k_corner_hook_multiplicities(),
                self.k_plus_one_corner_hook_multiplicities())

    def kk1_mult_from_left_right_seq(self):
        """
        EXAMPLES::

            sage: s=KShape([14,10,9,8,6,6,5,4,3,3,2,2,1,1,1],6)
            sage: s.kk1_mult()
            ((1, 2, 1, 1), (1, 2, 0))
            sage: s.kk1_mult_from_left_right_seq()
            ((1, 2, 1, 1), (1, 2, 0))
        """
        lrs = self.left_right_seq(True)
        lrs2 = list(lrs)
        for i in range(len(lrs2)-2, -1, -1):
            if isinstance(lrs2[i+1], SignedNum) and lrs2[i] == lrs2[i+1]:
                lrs2[i] = lrs2[i+1]
        res  = [int(i) for i in lrs  if isinstance(i, SignedNum)]
        res2 = [int(i) for i in lrs2 if isinstance(i, SignedNum)]
        return (tuple(res.count(-i) for i in range(2, self.k)),
                tuple(res2.count(i) for i in range(2, self.k-1)))

    def kkm1_mult(self):
        """
        EXAMPLES::

            sage: s=KShape([14,10,9,8,6,6,5,4,3,3,2,2,1,1,1],6)
            sage: s.kkm1_mult()
            ((1, 2, 1, 1), (1, 1, 0))
        """
        k = self.k
        resk   = [0]*(k-2)
        reskm1 = [0]*(k-3)
        for (r,c) in self.to_skew_partition().inner_corners():
            arml = self.outer_shape.arm_length(r,c)
            legl = self.outer_shape.leg_length(r,c)
            if arml+legl+1 == k:
                resk[arml-1] += 1
            if arml+legl+2 == k:
                reskm1[arml-1] += 1
        return (tuple(resk), tuple(reskm1))

    def k_corner_hook_multiplicities(self):
        """
        EXAMPLES ::

            sage: KShape([4,3,1,1], 4).k_corner_hook_multiplicities()
            (0, 1)
            sage: KShape([5,2,2,1,1], 5).k_corner_hook_multiplicities()
            (1, 0, 0)
            sage: s=KShape([14,10,9,8,6,6,5,4,3,3,2,2,1,1,1],6)
            sage: s.k_corner_hook_multiplicities()
            (1, 2, 1, 1)
        """
        k = self.k
        s=self.corner_hooks(k)
        return tuple(s.count((i,k+1-i)) for i in range(2, k))

    def k_minus_one_corner_hook_multiplicities(self):
        """
        EXAMPLES ::

            sage: KShape([4,3,1,1], 4).k_minus_one_corner_hook_multiplicities()
            (0,)
            sage: KShape([5,2,2,1,1], 5).k_minus_one_corner_hook_multiplicities()
            (0, 0)
            sage: s=KShape([14,10,9,8,6,6,5,4,3,3,2,2,1,1,1],6)
            sage: s.k_minus_one_corner_hook_multiplicities()
            (1, 1, 0)
        """
        k = self.k
        s=self.corner_hooks(k-1)
        return tuple(s.count((i,k-i)) for i in range(2, k-1))

    def Gandhi_stat(self):
        """
        Statistic de Gandhi

        EXAMPLES::

            sage: k = 4
            sage: sum(x^v.Gandhi_stat() for v in IrreducibleKShapes(k))
            6*x^4 + 8*x^3 + 3*x^2
        """
        return self.k-len(set(self.corner_hooks(self.k)))

    def kk1_stat(self):
        """
        EXAMPLES::

            sage: sum(x^sh.kk1_stat() for sh in IrreducibleKShapes(4))
            x^3 + 4*x^2 + 7*x + 5
        """
        k = self.k
        hks = self.outer_shape.hooks()
        return hks.count(k+1)+hks.count(k)

    def kkm1_stat(self):
        """
        EXAMPLES::

            sage: sum(x^sh.kkm1_stat() for sh in IrreducibleKShapes(4))
            x^3 + 4*x^2 + 7*x + 5
        """
        k = self.k
        return sum(sum(self.kkm1_mult(), ()))

    def kk1_hooks(self):
        """
        EXAMPLES::

            sage: sh = KShape([8, 6, 4, 3, 2, 1, 1, 1], 5)
            sage: sh.kk1_stat()
            4
            sage: sh.kk1_hooks()
            ((4, 2), (4, 3), (3, 3), (2, 4))

        TESTS::

            sage: all(len(sh.kk1_hooks()) == sh.kk1_stat()
            ...       for sh in IrreducibleKShapes(4))
            True
        """
        k = self.k
        outer = self.outer_shape
        cells = outer.cells()
        return tuple(filter(lambda (i,j): k+1 <= i+j <= k+2,
                            ((outer.arm_length(*cell)+1,
                              outer.leg_length(*cell)+1)
                             for cell in cells)))

    def kk1_vectors(self):
        """
        EXAMPLES::

            sage: sh = KShape([8, 6, 4, 3, 2, 1, 1, 1], 5)
            sage: sh.kk1_vectors()
            ((1, 1, 1), (0, 1))

        TESTS::

            sage: all(sum(sum(sh.kk1_vectors(), ())) == sh.kk1_stat()
            ...       for sh in IrreducibleKShapes(4))
            True
        """
        k = self.k
        res = self.kk1_hooks()
        return (tuple(1 if (i, k+1-i) in res else 0 for i in range(2, k)),
                tuple(1 if (i, k+2-i) in res else 0 for i in range(3, k)))

    def cuts(self):
        """
        Returns the cuts of ``self``, i.e. the outside corners of the outer shape which are adjacent to the inner shape.

        EXAMPLES::

            sage: sum(x^(len(p.cuts())) for p in IrreducibleKShapes(4))
            4*x + 13
        """
        return filter(lambda (x,y):(x-1,y-1) in self.inner_shape.corners(),self.outer_shape.outside_corners())

    def remove_cut(self, i, j):
        """
        Returns the k-shape obtained by removing the cut (``i``,``j``) of ``self``.

        EXAMPLES::

            sage: KShape([11,10,6,6,4,2,2,2,1,1],6).remove_cut(2,6)
            [[10, 6, 5, 4, 2, 2, 2, 1, 1], [5, 2, 2, 1]] with k=6
            sage: KShape([7,4,4,2,2,2,1],5).remove_cut(1,4)
            Traceback (most recent call last):
            ...
            ValueError: Result is not a k-shape
            sage: KShape([7,4,4,2,2,2,1],5).remove_cut(2,2)
            Traceback (most recent call last):
            ...
            ValueError: No such cut in the given k-shape
        """
        if (i,j) not in self.cuts():
            raise ValueError, "No such cut in the given k-shape"
        outer = list(self.outer_shape)
        inner = list(self.inner_shape)
        outer[0:i-1] = [l-1 for l in outer[0:i-1]]
        inner[0:i-1] = [l-1 for l in inner[0:i-1]]
        outer.pop(i-1)
        inner.pop(i-1)
        r = self.column_shape[j-1]
        outer[i:i+r-1] = [l-1 for l in outer[i:i+r-1]]
        if Partition(outer).k_boundary(self.k).inner() == inner:
            return KShape(outer, self.k, inner)
        else:
            raise ValueError, "Result is not a k-shape"

    def contains_rectangle(self, u, v):
        """
        EXAMPLES::

            sage: ks=KShape([7,4,4,2,2,2,1],5)
            sage: ks.contains_rectangle(3,1)
            True
            sage: ks.contains_rectangle(2,2)
            True
            sage: ks.contains_rectangle(3,2)
            False
        """
        l = [0]+[m for m,n in self.cuts()]+[len(self.outer_shape)]
        for a in range(len(l)-1):
            if (self.row_shape)[l[a]:l[a+1]] == [u]*v:
                return True
        return False

    def row_code(self):
        """
        EXAMPLES::

            sage: ks=KShape([7,4,4,2,2,2,1],5)
            sage: ks.row_shape
            [3, 2, 2, 2, 2, 2, 1]
            sage: ks.row_code()
            [1, 5, 1, 0]
        """
        return Partition(self.row_shape).to_exp(self.k-1)

def possible_shift(part, k, foot, hmax):
    foot = foot+[0]
    start = 0
    while foot[start] >= hmax:
        start += 1
    while foot[start] != 0 and foot[start]+part > k:
        # find next interface between two rectangles
        start += 1
        while foot[start] != 0 and foot[start] == foot[start-1]:
            start += 1
    if foot[start]+part < k or foot[start]==0:
        return [start]
    start2 = start+1
    while foot[start2] != 0 and foot[start2] == foot[start2-1]:
            start2 += 1
    return [start, start2]

def children(inp):
    (outer, inner, k, foot, hmax, parts)=inp
    if not parts:
        return []
    else:
        res=[]
        for shift in possible_shift(parts[0],k,foot,hmax):
            inner0 = inner[0] if inner else 0
            res.append(([inner0+shift+parts[0]]+outer,
                        [inner0+shift]+inner,
                        k,
                        [i+1 for i in foot[shift:]] + [1]*(parts[0]-len(foot)+shift),
                        ([hmax]+foot)[shift],
                        parts[1:]))
        return res

from sage.combinat.backtrack import SearchForest
def from_k_and_rows(k, rows):
    res = SearchForest([([],[],k,[],k,rows[::-1])], children)
    for (outer, inner, k, foot, hmax, parts) in res:
        if not parts:
            yield KShape(outer,k,inner,rows)


def from_k_and_weight(k, weight):
    """
    EXAMPLES::

        sage: from sage.combinat.kshape import from_k_and_weight
        sage: list(from_k_and_weight(4, 6))
        [[[6, 2], [2]] with k=4, [[5, 1, 1], [1]] with k=4, [[3, 3], []] with k=4, [[6, 3], [3]] with k=4, [[4, 2, 1], [1]] with k=4, [[5, 2, 1], [2]] with k=4, [[4, 1, 1, 1], [1]] with k=4, [[2, 2, 2], []] with k=4, [[4, 2, 2], [2]] with k=4, [[3, 2, 1, 1], [1]] with k=4, [[3, 3, 1, 1], [1, 1]] with k=4, [[3, 1, 1, 1, 1], [1]] with k=4, [[3, 2, 1, 1, 1], [1, 1]] with k=4, [[2, 2, 1, 1, 1, 1], [1, 1]] with k=4, [[2, 2, 2, 1, 1, 1], [1, 1, 1]] with k=4]
    """
    for rows in Partitions(weight,max_part=k):
        for x in from_k_and_rows(k,rows):
            yield x


def list_irreducible_old(k, s):
    r"""
    """
#    ks=filter(lambda skp :skp.outer().is_k_shape(k) and skp.outer().k_boundary(k)==skp,sage.combinat.skew_partition.SkewPartitions(s).list())
#    ks=map(lambda skp:KShape(skp.outer(),k),ks)
    ks = from_k_and_weight(k, s)
    return filter(lambda ksh:ksh.is_irreducible(),ks)


from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
class IrreducibleKShapes(UniqueRepresentation, SearchForest):
    """
    The enumerated set of irreducible k-shapes.
    """
    def __init__(self, k):
        """
        TESTS::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.cardinality()
            17
            sage: I4.list()
            [[[], []] with k=4, [[1], []] with k=4, [[1, 1], []] with k=4, [[2, 1, 1], []] with k=4, [[3, 2, 1, 1], [1]] with k=4, [[5, 3, 2, 1, 1], [2, 1]] with k=4, [[4, 2, 1, 1], [1]] with k=4, [[3, 1, 1], [1]] with k=4, [[4, 3, 1, 1], [1, 1]] with k=4, [[2, 1], []] with k=4, [[2, 2, 1], []] with k=4, [[4, 2, 2, 1], [2]] with k=4, [[5, 4, 2, 2, 1], [2, 2]] with k=4, [[4, 2, 1], [1]] with k=4, [[3, 1], []] with k=4, [[2], []] with k=4, [[3, 2], []] with k=4]

            sage: TestSuite(I4).run()
        """
        self.k = Integer(k)
        SearchForest.__init__(self, category = FiniteEnumeratedSets())

    def roots(self):
        r"""
        The root of enumerations as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: list(IrreducibleKShapes(4).roots())
            [([[], []] with k=4, [], 4)]
        """
        return [(KShape([],self.k), [], self.k)]

    def _repr_(self):
        """
        TEST::

            sage: IrreducibleKShapes(4)
            Irreducible 4-shapes
        """
        return "Irreducible %s-shapes"%(self.k)

    def children(self, inp):
        """
        The children function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        TEST::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.children((KShape([],4), [], 4))
            [([[1], []] with k=4, [1], 4), ([[2], []] with k=4, [1, 1], 4)]
        """
        (ksh, foot, hmax)=inp
        inner = ksh.inner_shape[:]
        outer = ksh.outer_shape[:]
        res=[]
        for new_parts in range(max(len(foot),1), self.k):
            for shift in possible_shift(new_parts,self.k,foot,hmax):
                inner0 = inner[0] if inner else 0
                newShape = KShape([inner0+shift+new_parts]+outer,
                                  self.k,
                                  [inner0+shift]+inner)
                if newShape.is_irreducible():
                    res.append((newShape,
                                ([i+1 for i in foot[shift:]] +
                                 [1]*(new_parts-len(foot)+shift)),
                                ([hmax]+foot)[shift]))
        return res

    def post_process(self, elt):
        """
        The post_process function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.post_process((KShape([],4), [], 4))
            [[], []] with k=4
        """
        return elt[0]

    def __contains__(self, x):
        """
        TESTS::

            sage: I4 = IrreducibleKShapes(4)
            sage: KShape([2,1], 4) in I4
            True
            sage: KShape([2,1], 3) in I4
            False
            sage: KShape([4], 4) in I4
            False
        """
        return isinstance(x, KShape) and self.k == x.k and x.is_irreducible()

class UnambiguousKShapes(IrreducibleKShapes):

    def __init__(self, k, irred=True):
        """
        TESTS::

            sage: I4 = UnambiguousKShapes(4)
            sage: I4.cardinality()
            5
            sage: TestSuite(I4).run()
            sage: I4bis = UnambiguousKShapes(4, irred=False)
            sage: I4bis.cardinality()
            14
            sage: TestSuite(I4bis).run()
        """
        self.k = Integer(k)
        self._irred = irred
        SearchForest.__init__(self, category = FiniteEnumeratedSets())

    def children(self, inp):
        """
        The children function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        TEST::

            sage: I4 = UnambiguousKShapes(4)
            sage: I4.children((KShape([],4), [], 4))
            [([[1], []] with k=4, [1], 4), ([[2], []] with k=4, [1, 1], 4)]
        """
        (ksh, foot, hmax)=inp
        inner = ksh.inner_shape[:]
        outer = ksh.outer_shape[:]
        res=[]
        for new_parts in range(max(len(foot),1), self.k):
            pshift = possible_shift(new_parts,self.k,foot,hmax)
            if len(pshift) == 1:
                shift = pshift[0]
                inner0 = inner[0] if inner else 0
                newShape = KShape([inner0+shift+new_parts]+outer,
                                  self.k,
                                  [inner0+shift]+inner)
                if not self._irred or newShape.is_irreducible():
                    res.append((newShape,
                                ([i+1 for i in foot[shift:]] +
                                 [1]*(new_parts-len(foot)+shift)),
                                ([hmax]+foot)[shift]))
        return res

    def _repr_(self):
        """
        TEST::

            sage: UnambiguousKShapes(4)
            Irreducible unambiguous 4-shapes
        """
        if self._irred:
            return "Irreducible unambiguous %s-shapes"%(self.k)
        else:
            return "Unambiguous %s-shapes"%(self.k)

    def __contains__(self, x):
        """
        TESTS::

            sage: I4 = UnambiguousKShapes(4)
            sage: KShape([2,1], 4) in I4
            True
            sage: KShape([2,1], 3) in I4
            False
            sage: KShape([4], 4) in I4
            False
        """
        return (isinstance(x, KShape) and
                self.k == x.k and
                (not self._irred or x.is_irreducible())
                and sum(x.left_right_seq()) == 0)

class KShapes_from_k_and_rows(UniqueRepresentation, SearchForest):
    """
    The enumerated set of k-shapes with given row shape.
    """

    @staticmethod
    def __classcall__(cls, k, rows):
        return super(KShapes_from_k_and_rows, cls).__classcall__(cls, k, Partition(rows))   

    def __init__(self, k, rows):
        self.k = k
        self.rows = rows
        SearchForest.__init__(self, category = FiniteEnumeratedSets())

    def roots(self):
        r"""
        The root of enumerations as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: list(KShapes_from_k_and_rows(4,[2,2,1]).roots())
            [([], [], 4, [], 4, [1, 2, 2])]
        """
        return [([],[],self.k,[],self.k,self.rows[::-1])]

    def __repr__(self):
        """
        TEST::

            sage: KShapes_from_k_and_rows(4,[2,2,1])
            4-shapes with row shape [2, 2, 1]
        """
        return "%s-shapes with row shape %s"%(self.k,self.rows)

    def children(self, inp):
        """
        The children function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        TEST::

            sage: S = KShapes_from_k_and_rows(4,[2,2,1])
            sage: S.children(([2, 1], [], 4, [2, 1], 4, [2]))
            [([2, 2, 1], [], 4, [3, 2], 4, []), ([3, 2, 1], [1], 4, [2, 1], 2, [])]
        """
        (outer, inner, k, foot, hmax, parts)=inp
        if not parts:
            return []
        else:
            res=[]
            for shift in possible_shift(parts[0],k,foot,hmax):
                inner0 = inner[0] if inner else 0
                res.append(([inner0+shift+parts[0]]+outer,
                            [inner0+shift]+inner if inner0+shift>0 else inner,
                            k,
                            [i+1 for i in foot[shift:]] + [1]*(parts[0]-len(foot)+shift),
                            ([hmax]+foot)[shift],
                            parts[1:]))
            return res

    def post_process(self, elt):
        """
        EXAMPLES::

            sage: S = KShapes_from_k_and_rows(4,[2,2,1])
            sage: S.post_process(([3, 2, 1], [1], 4, [2, 1], 2, []))
            [[3, 2, 1], [1]] with k=4
        """
        (outer, inner, k, foot, hmax, parts)=elt
        if not parts:
            return KShape(outer,k)
        else:
            return None
        return elt[0]


    def __contains__(self, x):
        """
        TESTS::

            sage: S = KShapes_from_k_and_rows(4,[2,2,1])
            sage: KShape([2,2,1], 4) in S
            True

        """
        return isinstance(x, KShape) and self.k == x.k and self.rows == x.row_shape


from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets, Family
class KShapes_from_k_and_weight(DisjointUnionEnumeratedSets):
    """
    The set of k-shapes with given weight.
    """
    def __init__(self,k,weight):
        """
        TESTS::

            sage: S = KShapes_from_k_and_weight(4,5)
            sage: S.cardinality()
            9
            sage: S.list()
            [[[5, 1], [1]] with k=4, [[3, 2], []] with k=4, [[5, 2], [2]] with k=4, [[4, 1, 1], [1]] with k=4, [[2, 2, 1], []] with k=4, [[3, 2, 1], [1]] with k=4, [[3, 1, 1, 1], [1]] with k=4, [[2, 1, 1, 1, 1], [1]] with k=4, [[2, 2, 1, 1, 1], [1, 1]] with k=4]

        """
        DisjointUnionEnumeratedSets.__init__(self,Family(Partitions(weight,max_part=k), lambda r:KShapes_from_k_and_rows(k,r)))
        self.k = k
        self.weight = weight

    def __repr__(self):
        """
        TEST::

            sage: KShapes_from_k_and_weight(4,10)
            4-shapes of weight 10
            
        """
        return "%s-shapes of weight %s"%(self.k,self.weight)

    def __contains__(self, x):
        """
        TESTS::

            sage: S = KShapes_from_k_and_weight(4,10)
            sage: KShape([9,5,1,1], 4) in S
            True
            sage: KShape([2,2,1], 4) in S
            False
            sage: KShape([8,5,2,2], 3) in S
            False

        """
        return (isinstance(x, KShape) 
            and self.k == x.k 
            and self.weight == x.outer_shape.size()-x.inner_shape.size())
            # self.weight == x.size()


class IrreducibleKShapes_from_kkm1_vectors(UniqueRepresentation, SearchForest):
    """
    The enumerated set of irreducible k-shapes.

    EXAMPLES::

        sage: from sage.combinat.kshape import *
        sage: S = IrreducibleKShapes_from_kkm1_vectors(((1,0,0,1),(1,0,0)))
        sage: def dict_vector_meth(it,meth):
        ...       return dict_vector_fun(it, attrcall(meth))

        sage: def dict_vector_fun(it,fun):
        ...       from collections import defaultdict
        ...       res  = defaultdict(list)
        ...       for s in it:
        ...           res[fun(s)].append(s)
        ...       return res

        sage: k = 5
        sage: dsh  = dict_vector_meth(IrreducibleKShapes(k).list(), "kkm1_mult")
        sage: all(set(dsh[v]) == set(IrreducibleKShapes_from_kkm1_vectors(v))
        ...       for v in dsh)
        True
    """
    def __init__(self, v):
        """
        TESTS::

        """
        self.vk, self.vkm1 = v
        assert len(self.vkm1) + 1 == len(self.vk)
        self.k = 2+len(self.vk)
        SearchForest.__init__(self, category = FiniteEnumeratedSets())

    def roots(self):
        r"""
        The root of enumerations as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: from sage.combinat.kshape import *
            sage: list(IrreducibleKShapes_from_kkm1_vectors(((1,0,0,1),(1,0,0))).roots())
            [([[], []] with k=6, [], 6)]
        """
        return [(KShape([],self.k), [], self.k)]

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.kshape import *
            sage: IrreducibleKShapes_from_kkm1_vectors(((1,0,0,1),(1,0,0)))
            Irreducible 6-shapes of k,k-1 vectors ((1, 0, 0, 1),(1, 0, 0))
        """
        return "Irreducible %s-shapes of k,k-1 vectors (%s,%s)"%(self.k, self.vk, self.vkm1)


    def is_prefix(self, ksh):
        """
        Return if ``ksh`` is the prefix of a k-shape of correct vector

            sage: from sage.combinat.kshape import *
            sage: S=IrreducibleKShapes_from_kkm1_vectors(((1,0,0,1),(1,0,0)))
            sage: S.is_prefix(KShape([6,3,3,3,2,1,1,1,1], 6))
            True
            sage: S.is_prefix(KShape([7,3,3,3,2,1,1,1,1], 6))
            True
            sage: S.is_prefix(KShape([5,2,1,1,1,1], 6))
            False

        .. warning:: we assume that ksh is irreducible.
        """
        vk, vkm1 = ksh.kkm1_mult()
        lfoot = max(partition_part(ksh.row_shape, 0), 1)
        lfoot2 = max(lfoot-2, 0)
        if not(vk[:lfoot2] == self.vk[:lfoot2] and
               vkm1[:lfoot2] == self.vkm1[:lfoot2] and
               vk[lfoot2] <= self.vk[lfoot2] + (
                1 if ksh.foot().hook_length(0,0) == self.k else 0)
               ):
            res = False
        elif lfoot > self.k-2:
            res = True
        else:
            res = (vkm1[lfoot2] <= self.vkm1[lfoot2] + (
                1 if ksh.foot().hook_length(0,0) == self.k-1 else 0))
        #from sage.misc.displayhook import print_obj
        #import sys
        #print_obj(sys.stdout, (ksh, (vk, vkm1), res))
        return res

    def children(self, inp):
        """
        The children function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        TEST::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.children((KShape([],4), [], 4))
            [([[1], []] with k=4, [1], 4), ([[2], []] with k=4, [1, 1], 4)]
        """
        (ksh, foot, hmax)=inp
        inner = ksh.inner_shape[:]
        outer = ksh.outer_shape[:]
        res=[]
        for new_parts in range(max(len(foot),1), self.k):
            for shift in possible_shift(new_parts,self.k,foot,hmax):
                inner0 = inner[0] if inner else 0
                newShape = KShape([inner0+shift+new_parts]+outer,
                                  self.k,
                                  [inner0+shift]+inner)
                if newShape.is_irreducible():
                    res.append((newShape,
                                ([i+1 for i in foot[shift:]] +
                                 [1]*(new_parts-len(foot)+shift)),
                                ([hmax]+foot)[shift]))
        return filter(lambda (ksh, foot, hmax): self.is_prefix(ksh), res)

    def post_process(self, elt):
        """
        The post_process function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.post_process((KShape([],4), [], 4))
            [[], []] with k=4
        """
        if elt[0] in self:
            return elt[0]

    def __contains__(self, x):
        """
        TESTS::

            sage: I4 = IrreducibleKShapes(4)
            sage: KShape([2,1], 4) in I4
            True
            sage: KShape([2,1], 3) in I4
            False
            sage: KShape([4], 4) in I4
            False
        """
        return (isinstance(x, KShape) and
                self.k == x.k and
                x.is_irreducible()
                and (self.vk, self.vkm1) == x.kkm1_mult())

class IrreducibleKShapes_from_kk1_vectors(UniqueRepresentation, SearchForest):
    """
    The enumerated set of irreducible k-shapes.

    EXAMPLES::

        sage: from sage.combinat.kshape import *
        sage: S = IrreducibleKShapes_from_kk1_vectors(((1,0,0,1),(1,0,0)))
        sage: def dict_vector_meth(it,meth):
        ...       return dict_vector_fun(it, attrcall(meth))

        sage: def dict_vector_fun(it,fun):
        ...       from collections import defaultdict
        ...       res  = defaultdict(list)
        ...       for s in it:
        ...           res[fun(s)].append(s)
        ...       return res

        sage: k = 5
        sage: dsh  = dict_vector_meth(IrreducibleKShapes(k).list(), "kk1_mult")
        sage: all(set(dsh[v]) == set(IrreducibleKShapes_from_kk1_vectors(v))
        ...       for v in dsh)
        True
    """
    def __init__(self, v):
        """
        TESTS::

        """
        self.vk, self.vk1 = v
        assert len(self.vk1) + 1 == len(self.vk)
        self.k = 2 + len(self.vk)
        SearchForest.__init__(self, category = FiniteEnumeratedSets())

    def roots(self):
        r"""
        The root of enumerations as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: from sage.combinat.kshape import *
            sage: list(IrreducibleKShapes_from_kk1_vectors(((1,0,0,1),(1,0,0))).roots())
            [([[], []] with k=6, [], 6)]
        """
        return [(KShape([],self.k), [], self.k)]

    def _repr_(self):
        """
        TEST::

            sage: from sage.combinat.kshape import *
            sage: IrreducibleKShapes_from_kk1_vectors(((1,0,0,1),(1,0,0)))
            Irreducible 6-shapes of k,k+1 vectors ((1, 0, 0, 1),(1, 0, 0))
        """
        return "Irreducible %s-shapes of k,k+1 vectors (%s,%s)"%(self.k, self.vk, self.vk1)


    def is_prefix(self, ksh):
        """
        Return if ``ksh`` is the prefix of a k-shape of correct vector

            sage: from sage.combinat.kshape import *
            sage: S=IrreducibleKShapes_from_kk1_vectors(((1,0,0,1),(1,0,0)))
            sage: S.is_prefix(KShape([2,2,1,1,1], 6))
            True
            sage: S.is_prefix(KShape([3,3,2,2,2,1,1], 6))
            True
            sage: S.is_prefix(KShape([6,3,3,3,2,1,1,1,1], 6))
            False

        .. warning:: we assume that ksh is irreducible.
        """
        vk, vk1 = ksh.kk1_mult()
        lfoot = max(partition_part(ksh.row_shape, 0), 1)
        lfoot2 = max(lfoot-2, 0)
        if not(vk[:lfoot2] == self.vk[:lfoot2] and
               vk1[:lfoot2] == self.vk1[:lfoot2] and
               vk[lfoot2] <= self.vk[lfoot2] + (
                1 if ksh.foot().hook_length(0,0) == self.k else 0)
               ):
            res = False
        elif lfoot > self.k-2:
            res = True
        else:
            res = (vk1[lfoot2] <= self.vk1[lfoot2] + (
                1 if ksh.foot().hook_length(0,0) == self.k+1 else 0))
        #from sage.misc.displayhook import print_obj
        #import sys
        #print_obj(sys.stdout, (ksh, (vk, vkm1), res))
        return res

    def children(self, inp):
        """
        The children function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        TEST::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.children((KShape([],4), [], 4))
            [([[1], []] with k=4, [1], 4), ([[2], []] with k=4, [1, 1], 4)]
        """
        (ksh, foot, hmax)=inp
        inner = ksh.inner_shape[:]
        outer = ksh.outer_shape[:]
        res=[]
        for new_parts in range(max(len(foot),1), self.k):
            for shift in possible_shift(new_parts,self.k,foot,hmax):
                inner0 = inner[0] if inner else 0
                newShape = KShape([inner0+shift+new_parts]+outer,
                                  self.k,
                                  [inner0+shift]+inner)
                if newShape.is_irreducible():
                    res.append((newShape,
                                ([i+1 for i in foot[shift:]] +
                                 [1]*(new_parts-len(foot)+shift)),
                                ([hmax]+foot)[shift]))
        return filter(lambda (ksh, foot, hmax): self.is_prefix(ksh), res)

    def post_process(self, elt):
        """
        The post_process function used by the search forest enumeration as per
        :class:`~sage.combinat.backtrack.SearchForest`

        EXAMPLES::

            sage: I4 = IrreducibleKShapes(4)
            sage: I4.post_process((KShape([],4), [], 4))
            [[], []] with k=4
        """
        if elt[0] in self:
            return elt[0]

    def __contains__(self, x):
        """
        TESTS::

            sage: I4 = IrreducibleKShapes(4)
            sage: KShape([2,1], 4) in I4
            True
            sage: KShape([2,1], 3) in I4
            False
            sage: KShape([4], 4) in I4
            False
        """
        return (isinstance(x, KShape) and
                self.k == x.k and
                x.is_irreducible()
                and (self.vk, self.vk1) == x.kk1_mult())

# do two moves intersect?
# contiguity?
