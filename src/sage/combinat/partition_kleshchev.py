r"""
Kleshchev partitions
"""
from __future__ import print_function, absolute_import

from .partition import Partition, Partitions
from .partition_tuple import PartitionTuple, PartitionTuples

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.all import NN, ZZ, IntegerModRing
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation


#--------------------------------------------------
# Kleshchev partition - element classes
#--------------------------------------------------
class KleshchevPartition(Partition):

    def conormal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of ``self`` which are conormal.
        If no residue ``i`` is specified then a list of length ``e``
        is returned which gives the conormal cells for ``0 <= i < e``.

        Following [Kleshchev09]_, the *conormal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and addable
        cells (in that order) from this list. The addable `i`-cells that
        remain at the end of the this process are the conormal `i`-cells.

        When computing conormal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention="regular")
            sage: KP([5,4,4,3,2]).conormal_cells()
            {0: [(1, 4)], 1: [(5, 0), (4, 2)]}
            sage: KP([5,4,4,3,2]).conormal_cells(0)
            [(1, 4)]
            sage: KP([5,4,4,3,2]).conormal_cells(1)
            [(5, 0), (4, 2)]
            sage: KP = KleshchevPartitions(3, convention="restricted")
            sage: KP([5,4,4,3,2]).conormal_cells()
            {0: [(1, 4), (3, 3)], 2: [(0, 5)]}
        """

        from collections import defaultdict
        # We use a dictionary for the conormal nodes as the indexing set is Z when e=0
        conormals = defaultdict(list)   # the conormal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells) - #(addable cells)

        # determine if we read up or down the partition
        KP = self.parent()
        rows = list(range(len(self)+1))
        if KP._convention[1]=='G':
            rows.reverse()

        # work through the rows
        for row in rows:
            if row == len(self): # addable cell at bottom of partition
                res = KP._multicharge[0] - row
                if carry[res] == 0:
                    conormals[res].append((row, 0))
                else:
                    carry[res] += 1
            else:
                res = KP._multicharge[0] + self[row] - row - 1
                if row == len(self)-1 or self[row] > self[row+1]: # removable cell
                    carry[res] -= 1
                if row == 0 or self[row-1] > self[row]:               #addable cell
                    if carry[res+1] >= 0:
                        conormals[res+1].append((row, self[row]))
                    else:
                        carry[res+1] += 1

        # finally return the result
        return dict(conormals) if i is None else conormals[i]

    def cogood_cells(self, i=None):
        r"""
        Return a list of the cells of ``self`` that are cogood.
        If no residue ``i`` is specified then the cogood cells of each
        residue are returned (if they exist).

        The cogood `i`-cell is the 'last' conormal `i`-cell. As with the
        conormal cells we can choose to read either up or down the partition.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention="regular")
            sage: KP([5,4,4,3,2]).cogood_cells()
            {0: (1, 4), 1: (4, 2)}
            sage: KP([5,4,4,3,2]).cogood_cells(0)
            (1, 4)
            sage: KP([5,4,4,3,2]).cogood_cells(1)
            (4, 2)
            sage: KP = KleshchevPartitions(4, convention='restricted')
            sage: KP([5,4,4,3,2]).cogood_cells()
            {1: (0, 5), 2: (4, 2), 3: (1, 4)}
            sage: KP([5,4,4,3,2]).cogood_cells(0)
            sage: KP([5,4,4,3,2]).cogood_cells(2)
            (4, 2)
        """
        conormal_cells = self.conormal_cells(i)
        if i is None:
            return {i: conormal_cells[i][-1] for i in conormal_cells}
        elif not conormal_cells:
            return None
        else:
            return conormal_cells[-1]

    def normal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of the partition which are normal.
        If no residue ``i`` is specified then a list of length ``e``
        is returned which gives the normal cells for ``0 <= i < e``.

        Following [Kleshchev09]_, the *normal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and
        addable cells (in that order) from this list. The removable
        `i`-cells that remain at the end of the this process are the
        normal `i`-cells.

        When computing normal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).normal_cells()
            {1: [(2, 3), (0, 4)]}
            sage: KP([5,4,4,3,2]).normal_cells(1)
            [(2, 3), (0, 4)]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).normal_cells()
            {0: [(4, 1)], 2: [(3, 2)]}
            sage: KP([5,4,4,3,2]).normal_cells(2)
            [(3, 2)]
        """
        from collections import defaultdict
        # We use a dictionary for the normal nodes as the indexing set is Z when e=0
        normals = defaultdict(list)     # the normal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        # determine if we read up or down the partition
        KP = self.parent()
        rows = list(range(len(self)+1))
        if KP._convention[1] == 'S':
            rows.reverse()

        # work through the rows
        for row in rows:
            if row == len(self): # addable cell at bottom of partition
                carry[KP._multicharge[0]-row] += 1
            else:
                res = KP._multicharge[0] + self[row] - row - 1
                if row == len(self) - 1 or self[row] > self[row+1]: # removable cell
                    if carry[res] == 0:
                        normals[res].insert(0, (row, self[row]-1))
                    else:
                        carry[res] -= 1
                if row == 0 or self[row-1] > self[row]:              # addable cell
                  carry[res+1] += 1

        # finally return the result
        return dict(normals) if i is None else normals[i]

    def good_cells(self, i=None):
        """
        Return a list of the cells of ``self`` that are good.
        If no residue ``i`` is specified then the good cells of each
        residue are returned (if they exist).

        The good `i`-cell is the 'first' normal `i`-cell. As with the normal
        cells we can choose to read either up or down the partition.

        EXAMPLES::

            sage: KP3 = KleshchevPartitions(3, convention='regular')
            sage: KP3([5,4,4,3,2]).good_cells()
            {1: (2, 3)}
            sage: KP3([5,4,4,3,2]).good_cells(1)
            (2, 3)
            sage: KP4 = KleshchevPartitions(4, convention='restricted')
            sage: KP4([5,4,4,3,2]).good_cells()
            {1: (2, 3)}
            sage: KP4([5,4,4,3,2]).good_cells(0)
            sage: KP4([5,4,4,3,2]).good_cells(1)
            (2, 3)
        """
        normal_cells = self.normal_cells(i)
        if i is None:
            return {j: normal_cells[j][0] for j in normal_cells}
        elif not normal_cells:
            return None
        else:
            return normal_cells[0]

    def good_residue_sequence(self):
        """
        Return a sequence of good nodes from the empty partition
        to ``self``, or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).good_residue_sequence()
            [0, 2, 1, 1, 0, 2, 0, 2, 1, 1, 0, 2, 0, 2, 2, 0, 1, 1]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).good_residue_sequence()
            [0, 1, 2, 2, 0, 1, 0, 2, 1, 2, 0, 1, 0, 2, 1, 2, 1, 0]
        """
        if not self:
            return []

        good_cells = self.good_cells()
        if not good_cells:
            return None

        res = sorted(good_cells.keys())[0]
        r,c = good_cells[res]
        good_seq = type(self)(self.parent(), self.remove_cell(r,c)).good_residue_sequence()
        if good_seq is None:
            return None
        good_seq.append( self.parent()._index_set(res) )
        return good_seq

    def good_cell_sequence(self):
        """
        Return a sequence of good nodes from the empty partition
        to ``self``, or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([5,4,4,3,2]).good_cell_sequence()
            [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2),
             (3, 0), (2, 1), (1, 2), (3, 1), (0, 3), (1, 3),
             (2, 2), (3, 2), (4, 0), (4, 1), (0, 4), (2, 3)]
            sage: KP = KleshchevPartitions(3, convention='restricted')
            sage: KP([5,4,4,3,2]).good_cell_sequence()
            [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0),
             (0, 3), (2, 1), (1, 2), (1, 3), (3, 0), (3, 1),
             (2, 2), (4, 0), (2, 3), (3, 2), (0, 4), (4, 1)]
        """
        if not self:
            return []

        good_cells = self.good_cells()
        if not good_cells:
            return None

        cell = sorted(good_cells.values())[0]
        good_seq = type(self)(self.parent(), self.remove_cell(*cell)).good_cell_sequence()
        if good_seq is None:
            return None
        good_seq.append(cell)
        return good_seq

    def mullineux_conjugate(self):
        """
        Return the partition tuple which is the Mullineux conjugate of this
        partition tuple, or ``None`` if no such partition tuple exists.

        EXAMPLES::

            sage: KleshchevPartitions(3, convention='regular')([5,4,4,3,2]).mullineux_conjugate()
            [9, 7, 1, 1]
            sage: KleshchevPartitions(3, convention='restricted')([5,4,4,3,2]).mullineux_conjugate()
            [3, 2, 2, 2, 2, 2, 2, 1, 1, 1]
        """
        if not self:
            return self
        good_cells = self.good_cells()
        assert good_cells

        r,c = sorted(good_cells.values())[0]
        mu = type(self)(self.parent(), self.remove_cell(r, c)).mullineux_conjugate()
        # add back on a cogood cell of residue -residue(k,r,c)
        return type(self)(self.parent(), mu.add_cell(*mu.cogood_cells( r-c-self.parent()._multicharge[0]) ))

    def is_regular(self):
        """
        Return ``True`` if ``self`` is a `e`-regular partition tuple.

        A partition tuple is `e`-regular if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the down direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2)
            sage: KP([2,1,1]).is_regular()
            False
            sage: KP = KleshchevPartitions(3)
            sage: KP([2,1,1]).is_regular()
            True
            sage: KP([]).is_regular()
            True
        """
        if self.size() == 0 or self.parent()._e == 0:
            return True
        KP = self.parent()
        return super(KleshchevPartition, self).is_regular(KP._e, KP._multicharge)

    def is_restricted(self):
        """
        Return ``True`` if ``self`` is an `e`-restricted partition tuple.

        A partition tuple is `e`-restricted if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the up direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, convention='regular')
            sage: KP([3,1]).is_restricted()
            False
            sage: KP = KleshchevPartitions(3, convention='regular')
            sage: KP([3,1]).is_restricted()
            True
            sage: KP([]).is_restricted()
            True
        """
        if self.size() == 0 or self.parent()._e == 0:
            return True
        KP = self.parent()
        return super(KleshchevPartition, self).is_restricted(KP._e, KP._multicharge)

class KleshchevPartitionTuple(PartitionTuple):
    def conormal_cells(self, i=None):
        r"""
        Return a dictionary of the cells of the partition that are conormal.
        If no residue ``i`` is specified then a list of length ``e``
        is returned which gives the conormal cells for ``0 <= i < e``.

        Following [Kleshchev09]_, the *conormal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and addable
        cells (in that order) from this list. The addable `i`-cells that
        remain at the end of the this process are the conormal `i`-cells.

        When computing conormal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        EXAMPLES::

            sage: KP = KleshchevPartitions(3, [0,1], convention='left regular')
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells()
            {0: [(1, 2, 1), (1, 1, 3), (1, 0, 5)],
            1: [(1, 3, 0), (0, 2, 0), (0, 1, 2), (0, 0, 4)]}
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(1)
            [(1, 3, 0), (0, 2, 0), (0, 1, 2), (0, 0, 4)]
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(2)
            []
            sage: KP = KleshchevPartitions(3, [0,1], convention='right restricted')
            sage: KP([[4, 2], [5, 3, 1]]).conormal_cells(0)
            [(1, 0, 5), (1, 1, 3), (1, 2, 1)]
        """
        from collections import defaultdict
        # We use a dictionary for the conormal nodes as the indexing set is Z when e=0
        conormals = defaultdict(list)   # the conormal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        # the indices for the rows ending in addable nodes
        KP = self.parent()
        if KP._convention[0] == 'L':
            rows = [(k,r) for k,part in enumerate(self) for r in range(len(part)+1)]
        else:
            rows = [(k,r) for k,part in reversed(list(enumerate(self))) for r in range(len(part)+1)]
        if KP._convention[1] == 'G':
            rows.reverse()

        for row in rows:
            k,r = row
            if r == len(self[k]): # addable cell at bottom of a component
                res = KP._multicharge[k] - r
                if carry[res] == 0:
                    conormals[res].append((k, r, 0))
                else:
                    carry[res] += 1
            else:
                res = KP._multicharge[k] + self[k][r] - r - 1
                if r == len(self[k]) - 1 or self[k][r] > self[k][r+1]: # removable cell
                    carry[res] -= 1
                if r == 0 or self[k][r-1] > self[k][r]:                # addable cell
                    if carry[res+1] == 0:
                        conormals[res+1].append((k, r, self[k][r]))
                    else:
                        carry[res+1] += 1

        # finally return the result
        if i is None:
            return dict(conormals)
        return conormals[i]

    def cogood_cells(self, i=None):
        r"""
        Return a list of the cells of the partition that are cogood.
        If no residue ``i`` is specified then the cogood cells of each
        residue are returned (if they exist).

        The cogood `i`-cell is the 'last' conormal `i`-cell. As with the
        conormal cells we can choose to read either up or down the partition.

        EXAMPLES::

            sage: pt = KleshchevPartitions(3, [0,1])([[4, 2], [5, 3, 1]])
            sage: pt.cogood_cells()
            {0: (1, 2, 1), 1: (1, 3, 0)}
            sage: pt.cogood_cells(0)
            (1, 2, 1)
            sage: pt = KleshchevPartitions(4, [0,1], convention='left regular')([[4, 2], [5, 3, 1]])
            sage: pt.cogood_cells()
            {1: (0, 0, 5), 2: (1, 3, 0), 3: (0, 1, 4)}
            sage: pt.cogood_cells(0)
            sage: pt.cogood_cells(1) is None
            False
        """
        conormal_cells = self.conormal_cells(i)
        if i is None:
            return {j: conormal_cells[j][-1] for j in conormal_cells}
        elif not conormal_cells:
            return None
        else:
            return conormal_cells[-1]

    def normal_cells(self, i=None):
        r"""
        Return a dictionary of the removable cells of the partition that
        are normal. If no residue ``i`` is specified then a list of length
        ``e`` is returned which gives the normal cells for ``0 <= i < e``.

        Following [Kleshchev09]_, the *normal* cells are computed by
        reading up (or down) the rows of the partition and marking all
        of the addable and removable cells of `e`-residue `i` and then
        recursively removing all adjacent pairs of removable and
        addable cells (in that order) from this list. The removable
        `i`-cells that remain at the end of the this process are the
        normal `i`-cells.

        When computing normal cells you can either read the cells in order
        from top to bottom (this corresponds to labeling the simple modules
        of the symmetric group by regular partitions) or from bottom to top
        (corresponding to labeling the simples by restricted partitions).
        By default we read down the partition but this can be changed by
        setting ``convention = 'RS'``.

        EXAMPLES::

            sage: KP=KleshchevPartitions(3, [0,1])
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {0: [(1, 2, 1)], 2: [(1, 1, 2)]}
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells(1)
            []
            sage: KP=KleshchevPartitions(3, [0,1], convention='left regular')
            sage: KP([[4, 2], [5, 3, 1]]).normal_cells()
            {1: [(0, 0, 4)]}
        """
        from collections import defaultdict
        # We use a dictionary for the normal nodes as the indexing set is Z when e=0
        normals = defaultdict(list)     # the normal cells of each residue
        carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)

        KP = self.parent()
        if KP._convention[0] == 'R':
            rows = [(k,r) for k,part in enumerate(self) for r in range(len(part)+1)]
        else:
            rows = [(k,r) for k,part in reversed(list(enumerate(self))) for r in range(len(part)+1)]
        if KP._convention[1] == 'S':
            rows.reverse()

        for row in rows:
            k,r = row
            if r == len(self[k]): # addable cell at bottom of a component
                carry[KP._multicharge[k]-r] += 1
            else:
                res = KP._multicharge[k] + self[k][r] - r - 1
                if r == len(self[k])-1 or self[k][r] > self[k][r+1]: # removable cell
                    if carry[res] == 0:
                        normals[res].insert(0, (k, r, self[k][r]-1))
                    else:
                        carry[res] -= 1
                if r == 0 or self[k][r-1] > self[k][r]:               #addable cell
                    carry[res+1] += 1

        # finally return the result
        if i is None:
            return dict(normals)    # change the defaultdict into a dict
        else:
            return normals[i]

    def good_cells(self, i=None):
        r"""
        Return a list of the cells of the partition tuple which are good.
        If no residue ``i`` is specified then the good cells of each
        residue are returned (if they exist).

        The good `i`-cell is the 'first' normal `i`-cell. As with the
        normal cells we can choose to read either up or down the partition.

        EXAMPLES::

            sage: KP=KleshchevPartitions(3,[0,1])
            sage: pt = KP([[4, 2], [5, 3, 1]])
            sage: pt.good_cells()
            {0: (1, 2, 1), 2: (1, 1, 2)}
            sage: pt.good_cells(0)
            (1, 2, 1)
            sage: KP=KleshchevPartitions(4,[0,1], convention='left regular')
            sage: pt = KP([[4, 2], [5, 3, 1]])
            sage: pt.good_cells()
            {0: (1, 2, 1), 2: (0, 1, 3)}
            sage: pt.good_cells(1) is None
            True
        """
        normal_cells = self.normal_cells(i)
        if i is None:
            return {j: normal_cells[j][0] for j in normal_cells}
        elif not normal_cells:
            return None
        else:
            return normal_cells[0]

    def good_residue_sequence(self):
        """
        Return a sequence of good nodes from the empty partition to ``self``,
        or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP=KleshchevPartitions(3,[0,1])
            sage: KP([[4, 2], [5, 3, 1]]).good_residue_sequence()
            [1, 2, 0, 1, 0, 2, 0, 1, 2, 1, 2, 2, 2, 0, 0]
        """
        if self.size() == 0:
            return []
        good_cells = self.good_cells()
        if not good_cells:
            return None

        res = sorted(good_cells.keys())[0]
        k,r,c = good_cells[res]
        good_seq = type(self)(self.parent(), self.remove_cell(k,r,c)).good_residue_sequence()
        if good_seq is None:
            return None
        good_seq.append( self.parent()._index_set(res) )
        return good_seq

    def good_cell_sequence(self):
        """
        Return a sequence of good nodes from the empty partition to this
        partition, or ``None`` if no such sequence exists.

        EXAMPLES::

            sage: KP=KleshchevPartitions(3,[0,1])
            sage: KP([[4, 2], [5, 3, 1]]).good_cell_sequence()
            [(1, 0, 0),
             (1, 0, 1),
             (0, 0, 0),
             (1, 1, 0),
             (0, 0, 1),
             (1, 1, 1),
             (1, 0, 2),
             (1, 0, 3),
             (0, 1, 0),
             (0, 0, 2),
             (1, 2, 0),
             (1, 1, 2),
             (1, 0, 4),
             (0, 1, 1),
             (0, 0, 3)]
        """
        if self.size() == 0:
            return []
        good_cells = self.good_cells()
        if not good_cells:
            return None

        cell = sorted(good_cells.values())[0]
        good_seq = type(self)(self.parent(), self.remove_cell(*cell)).good_cell_sequence()
        if good_seq is None:
            return None
        good_seq.append(cell)
        return good_seq

    def mullineux_conjugate(self):
        """
        Return the partition tuple which is the Mullineux conjugate of
        ``self``, or ``None`` if no such partition tuple exists.

        EXAMPLES::

            sage: KP=KleshchevPartitions(3,[0,1])
            sage: KP([[4, 2], [5, 3, 1]]).mullineux_conjugate()

        """
        if self.size() == 0:
            return self

        good_cells = self.good_cells()
        if not good_cells:
            return None

        k,r,c = sorted(good_cells.values())[0]
        mu = type(self)(self.parent(), self.remove_cell(k,r,c)).mullineux_conjugate()
        # add back on a cogood cell of residue -residue(k,r,c)
        return type(self)(self.parent(), mu.add_cell(*mu.cogood_cells( r-c-self.parent()._multicharge[k])))

    def is_regular(self):
        """
        Return ``True`` if ``self`` is a `e`-regular partition tuple.

        A partition tuple is `e`-regular if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the down direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, [0,2], convention='right restricted')
            sage: KP([[2,1,1], [3,2]]).is_regular()
            False
            sage: KP = KleshchevPartitions(3, [0,2], convention='right restricted')
            sage: KP([[3,1,1], [3,2]]).is_regular()
            True
            sage: KP([[], []]).is_regular()
            True
        """
        if self.size() == 0:
            return True
        KP = self.parent()
        return is_regular(self.to_list(), KP._multicharge, KP._convention)

    def is_restricted(self):
        """
        Return ``True`` if ``self`` is an `e`-restricted partition tuple.

        A partition tuple is `e`-restricted if we can get to the
        empty partition tuple by successively removing a sequence
        of good cells in the up direction.

        EXAMPLES::

            sage: KP = KleshchevPartitions(2, [0,2], convention='left regular')
            sage: KP([[3,2,1], [3,1,1]]).is_restricted()
            False
            sage: KP = KleshchevPartitions(3, [0,2], convention='left regular')
            sage: KP([[3,2,1], [3,1,1]]).is_restricted()
            True
            sage: KP([[], []]).is_restricted()
            True
        """
        if self.size() == 0:
            return True
        KP = self.parent()
        return is_restricted(self.to_list(), KP._multicharge, KP._convention)

class KleshchevCrystalMixin(object):
    """
    Mixin class for the crystal structure of a Kleshchev partition.
    """
    def epsilon(self, i):
        r"""
        Return the Kashiwara crystal operator `\varepsilon_i` applied to ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: [x.epsilon(i) for i in C.index_set()]
            [0, 3, 0]
        """
        P = self.parent()
        return len(self.normal_cells(i))

    def phi(self, i):
        r""", convention='left regular'
        Return the Kashiwara crystal operator `\varphi_i` applied to ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: [x.phi(i) for i in C.index_set()]
            [3, 2, 0]
        """
        P = self.parent()
        return len(self.conormal_cells(i))

    def Epsilon(self):
        r"""
        Return `\varepsilon` of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.Epsilon()
            3*Lambda[1]
        """
        P = self.parent()
        WLR = P.weight_lattice_realization()
        La = WLR.fundamental_weights()
        n = self.normal_cells()
        return WLR.sum(len(n[i])*La[i] for i in P.index_set() if i in n)

    def Phi(self):
        r"""
        Return `\phi` of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.Phi()
            3*Lambda[0] + 2*Lambda[1]
        """
        P = self.parent()
        WLR = P.weight_lattice_realization()
        La = WLR.fundamental_weights()
        c = self.conormal_cells()
        return WLR.sum(len(c[i])*La[i] for i in P.index_set() if i in c)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.weight()
            3*Lambda[0] - Lambda[1] - 5*delta
            sage: x.Phi() - x.Epsilon()
            3*Lambda[0] - Lambda[1]
        """
        WLR = self.parent().weight_lattice_realization()
        alpha = WLR.simple_roots()
        La = WLR.fundamental_weights()
        r = self.parent()._multicharge
        wt = WLR.sum(La[ZZ(x)] for x in r)
        return wt - WLR.sum(alpha[self.content(*c, multicharge=r)]
                            for c in self.cells())

class KleshchevPartitionCrystal(KleshchevPartition, KleshchevCrystalMixin):
    """
    Kleshchev partition with the crystal structure.
    """
    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, convention='left regular')
            sage: x = C([5,4,1])
            sage: x.e(0)
            sage: x.e(1)
            [5, 4]
        """
        P = self.parent()
        cell = self.good_cells(i)
        if cell is None:
            return None
        r,c = cell
        mu = list(self)
        mu[r] -= 1
        return type(self)(P, mu)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, convention='left regular')
            sage: x = C([5,4,1])
            sage: x.f(0)
            [5, 4, 1, 1]
            sage: x.f(1)
            sage: x.f(2)
            [6, 4, 1]
        """
        P = self.parent()
        cell = self.cogood_cells(i)
        if cell is None:
            return None
        r,c = cell
        mu = list(self)
        if c == 0:
            mu.append(1)
        else:
            mu[r] += 1
        return type(self)(P, mu)

class KleshchevPartitionTupleCrystal(KleshchevPartitionTuple, KleshchevCrystalMixin):
    """
    Kleshchev partition tuple with the crystal structure.
    """
    def e(self, i):
        r"""
        Return the action of `e_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.e(0)
            sage: x.e(1)
            ([5, 4, 1], [2, 2, 1, 1])
        """
        P = self.parent()
        cell = self.good_cells(i)
        if cell is None:
            return None
        k,r,c = cell
        mu = self.to_list()
        mu[k][r] -= 1
        return type(self)(P, mu)

    def f(self, i):
        r"""
        Return the action of `f_i` on ``self``.

        INPUT:

        - ``i`` -- an element of the index set

        EXAMPLES::

            sage: C = crystals.KleshchevPartitions(3, [0,2], convention='left regular')
            sage: x = C([[5,4,1],[3,2,1,1]])
            sage: x.f(0)
            ([5, 5, 1], [3, 2, 1, 1])
            sage: x.f(1)
            ([5, 4, 1], [3, 2, 2, 1])
            sage: x.f(2)
        """
        P = self.parent()
        cell = self.cogood_cells(i)
        if cell is None:
            return None
        k,r,c = cell
        mu = self.to_list()
        if c == 0:
            mu[k].append(1)
        else:
            mu[k][r] += 1
        return type(self)(P, mu)

#--------------------------------------------------
# Kleshchev partitions - parent classes
#--------------------------------------------------
class KleshchevPartitions(PartitionTuples):
    """
    Kleshchev partitions

    A partition (tuple) `\mu` is restricted if it can be recursively
    obtained by adding a sequence of good nodes to the empty
    :class:`PartitionTuple` of the same :meth:`~PartitionTuple.level`.

    EXAMPLES::

        sage: KleshchevPartitions(5,[3,2,1],1,convention='RS')[:]
        [([], [], [1]), ([], [1], []), ([1], [], [])]
        sage: KleshchevPartitions(5,[3,2,1],1,convention='LS')[:]
        [([], [], [1]), ([], [1], []), ([1], [], [])]
        sage: KleshchevPartitions(5,[3,2,1],3)[:]
        [([], [], [2, 1]),
        ([1], [], [1, 1]),
        ([], [], [1, 1, 1]),
        ([], [1], [2]),
        ([], [], [3]),
        ([1], [], [2]),
        ([2], [], [1]),
        ([], [1], [1, 1]),
        ([], [1, 1], [1]),
        ([], [2], [1]),
        ([1], [2], []),
        ([], [3], []),
        ([1], [1], [1]),
        ([1, 1], [1], []),
        ([2], [1], []),
        ([3], [], [])]
        sage: KleshchevPartitions(5,[3,2,1],3,convention="LS")[:]
        [([], [1], [1, 1]),
        ([1], [], [1, 1]),
        ([], [], [1, 1, 1]),
        ([], [1, 1], [1]),
        ([], [1], [2]),
        ([1], [1], [1]),
        ([1, 1], [], [1]),
        ([2], [], [1]),
        ([], [1, 1, 1], []),
        ([], [1, 1], [1]),
        ([1], [1, 1], []),
        ([1, 1], [1], []),
        ([1], [2], []),
        ([2], [1], []),
        ([1, 1, 1], [], []),
        ([1, 1], [1], []),
        ([2, 1], [], []),
        ([3], [], [])]

    REFERENCES:

    .. [Kleshchev09] Alexander Kleshchev.
       *Representation theory of symmetric groups and related Hecke algebras*.
       :arxiv:`0909.4844`.
    """
    @staticmethod
    def __classcall_private__(cls, e, multicharge=(0,), size=None,
            convention='left restricted'):
        r"""
        This is a factory class which returns the appropriate parent based on
        the values of `level` and `size`.
        
        EXAMPLES::

            sage: KleshchevPartitions(5, [3,2,1], 1, convention='RS')[:]
            [([], [], [1]), ([], [1], []), ([1], [], [])]
            sage: KleshchevPartitions(5, [3,2,1], 1, convention='LS')[:]
            [([], [], [1]), ([], [1], []), ([1], [], [])]
        """
        if size is None and multicharge in ZZ:
            size = ZZ(multicharge)
            multicharge = (0,)

        I = IntegerModRing(e)
        multicharge = tuple([I(x) for x in multicharge])

        convention = convention.upper()
        if 'S' in convention:
            convention = convention[0]+'S'
        elif 'G' in convention:
            convention = convention[0]+'G'
        if convention not in ['RG','LG', 'RS', 'LS']:
            raise ValueError('invalid convention')

        if size is None:
            return KleshchevPartitions_all(e, multicharge, convention)
        else:
            return KleshchevPartitions_size(e, multicharge, size, convention)

class KleshchevPartitions_all(KleshchevPartitions):
    r"""
    Class of all Kleshchev partitions.

    .. RUBRIC:: Crystal structure

    We consider type `A_{e-1}^{(1)}` crystals, and let `r = (r_i \mid
    r_i \in \ZZ / e \ZZ)` be a finite sequence of length `k`, which
    is the *level*, and `\lambda = \sum_i \Lambda_{r_i}`. We will
    model the highest weight `U_q(\mathfrak{g})`-crystal `B(\lambda)`
    by a particular subset of partition tuples of level `k`.

    Consider a partition tuple `\mu` with multicharge `r`.
    We define `e_i(\mu)` as the partition tuple obtained after the
    deletion of the `i`-:meth:`good cell <PartitionTuple.good_cell>` 
    to `\mu` and `0` if there is no `i`-good cell. We define `f_i(\mu)`
    as the partition tuple obtained by the addition of the
    `i`-:meth:`cogood cell <PartitionTuple.cogood_cell>` to `\mu`
    and `0` if there is no `i`-good cell.

    The crystal `B(\lambda)` is the crystal generated by the empty
    partition tuple. We can compute the weight of an element `\mu` by taking
    `\lambda - \sum_{i=0}^n c_i \alpha_i` where `c_i` is the number of cells
    of `n`-residue `i` in `\mu`. Partition tuples in the crystal are known
    as *Kleshchev partitions*.

    .. NOTE::

        We can describe normal (not restricted) Kleshchev partition tuples
        in `B(\lambda)` as partition tuples `\mu` such that
        `\mu^{(t)}_{r_t - r_{t+1} + x} < \mu^{(t+1)}_x`
        for all `x \geq 1` and `1 \leq t \leq k - 1`.

    INPUT:

    - ``e`` -- for type `A_{e-1}^{(1)}` or `0`
    - ``multicharge`` -- the multicharge sequence `r`
    - ``convention`` -- (default: ``'LS'``) the reading convention

    EXAMPLES:

    We first do an example of a level 1 crystal::

        sage: C = crystals.KleshchevPartitions(3, [0], convention='left restricted')
        sage: C
        Kleshchev partitions with e=3
        sage: mg = C.highest_weight_vector()
        sage: mg
        []
        sage: mg.f(0)
        [1]
        sage: mg.f(1)
        sage: mg.f(2)
        sage: mg.f_string([0,2,1,0])
        [1, 1, 1, 1]
        sage: mg.f_string([0,1,2,0])
        [2, 2]
        sage: S = C.subcrystal(max_depth=5)
        sage: G = C.digraph(subset=S)
        sage: B = crystals.LSPaths(['A',2,1], [1,0,0])
        sage: SB = B.subcrystal(max_depth=5)
        sage: GB = B.digraph(subset=SB)
        sage: G.is_isomorphic(GB, edge_labels=True)
        True

    Now a higher level crystal::

        sage: C = crystals.KleshchevPartitions(3, [0,2], convention='RS')
        sage: mg = C.highest_weight_vector()
        sage: mg
        ([], [])
        sage: mg.f(0)
        ([1], [])
        sage: mg.f(2)
        ([], [1])
        sage: mg.f_string([0,1,2,0])
        ([2, 2], [])
        sage: mg.f_string([0,2,1,0])
        ([1, 1, 1, 1], [])
        sage: mg.f_string([2,0,1,0])
        ([2], [2])
        sage: S = C.subcrystal(max_depth=3)
        sage: G = C.digraph(subset=S)
        sage: B = crystals.LSPaths(['A',2,1], [1,0,1])
        sage: SB = B.subcrystal(max_depth=3)
        sage: GB = B.digraph(subset=SB)
        sage: G.is_isomorphic(GB, edge_labels=True)
        True

    The ordering of the residues gives a different representation of the
    higher level crystals (but it is still isomorphic)::

        sage: C2 = crystals.KleshchevPartitions(3, [2,0], convention='RS')
        sage: mg2 = C2.highest_weight_vector()
        sage: mg2.f_string([0,1,2,0])
        ([2], [2])
        sage: mg2.f_string([0,2,1,0])
        ([1, 1, 1], [1])
        sage: mg2.f_string([2,0,1,0])
        ([2, 1], [1])
        sage: S2 = C2.subcrystal(max_depth=5)
        sage: G2 = C2.digraph(subset=S)
        sage: G.is_isomorphic(G2, edge_labels=True)
        True

    REFERENCES:

    .. [Ariki2001] Susumu Ariki. On the classification of simple modules for
       cyclotomic Hecke algebras of type `G(m,1,n)` and Kleshchev
       multipartitions. Osaka J. Math. **38** (2001). :arxiv:`9908004v2`.

    .. [Vazirani2002] Monica Vazirani. *Parameterizing Hecek algebra modules:
       Bernstein-Zelevinsky multisegments, Kleshchev multipartitions, and
       crystal graphs*. Transform. Groups **7** (2002). pp. 267-303.
       :arxiv:`0107052v1`, :doi:`10.1007/s00031-002-0014-1`.

    .. [TingleyLN] Peter Tingley. Explicit `\widehat{\mathfrak{sl}}_n` crystal
       maps between cylindric plane partitions, multi-partitions, and
       multi-segments. Lecture notes.
       http://webpages.math.luc.edu/~ptingley/lecturenotes/explicit_bijections.pdf

    .. [Tingley2007] Peter Tingley. Three combinatorial models for
       `\widehat{\mathfrak{sl}}_n` crystals, with applications to cylindric
       plane partitions. International Mathematics Research Notices. (2007).
       :arxiv:`0702062v3`.
    """

    def __init__(self, e, multicharge, convention):
        r"""
        Initializes classes of PartitionTuples.

        EXAMPLES::

            sage: K = KleshchevPartitions(4, [2])
            sage: TestSuite(K).run()  # long time
            sage: K = KleshchevPartitions(4, [0,2,1])
            sage: TestSuite(K).run()  # long time
        """
        if e not in NN or e == 1:
            raise ValueError('e must belong to {0,2,3,4,5,6,...}')
        if e > 0:
            from sage.combinat.root_system.cartan_type import CartanType
            from sage.categories.highest_weight_crystals import HighestWeightCrystals
            from sage.categories.regular_crystals import RegularCrystals
            self._cartan_type = CartanType(['A', e-1, 1])
            cat = (HighestWeightCrystals(), RegularCrystals().Infinite())
        else:
            cat = InfiniteEnumeratedSets()

        self._level = len(multicharge)
        if self._level == 1:
            self.Element = KleshchevPartitionCrystal
            self._element_constructor_ = Partitions._element_constructor_.__func__
        else:
            self.Element = KleshchevPartitionTupleCrystal

        super(KleshchevPartitions_all, self).__init__(category=cat)
        self._e = e   # for printing
        self._index_set = IntegerModRing(e)
        self._multicharge = multicharge
        self._convention = convention
        if e > 0:
            if self._level == 1:
                self.module_generators = (self.element_class(self, []),)
            else:
                self.module_generators = (self.element_class(self, [[]]*self._level),)

    def _repr_(self):
        """
        EXAMPLES::

            sage: KleshchevPartitions(4, [2])
            Kleshchev partitions with e=4
            sage: KleshchevPartitions(3,[0,0,0])
            Kleshchev partitions with e=3 and multicharge=(0,0,0)
            sage: KleshchevPartitions(3,[0,0,1])
            Kleshchev partitions with e=3 and multicharge=(0,0,1)
        """
        if self._level == 1:
            return 'Kleshchev partitions with e=%s' % (self._e)
        else:
            return 'Kleshchev partitions with e=%s and multicharge=(%s)' % (
                        self._e,','.join('%s'%m for m in self._multicharge))

    def __contains__(self, mu):
        """
        Containment test for Kleshchev partitions.

        EXAMPLES::

            sage: PartitionTuple([[3,2],[2]]) in KleshchevPartitions(2, [0,0], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,0,1], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,1,1], 7)
            False
            sage: PartitionTuple([[],[2,1],[3,2]]) in KleshchevPartitions(5, [0,1,1], 8)
            True
            sage: all(mu in PartitionTuples(3,8) for mu in KleshchevPartitions(2, [0,0,0], 8))
            True
        """
        if isinstance(mu, (KleshchevPartition, KleshchevPartitionTuple)):
            if mu.level() != self._level:
                return False
            mu = self.element_class(self, list(mu))
            if self._convention[1] == 'G':
                return mu.is_regular()
            else:
                return mu.is_restricted()

        try:
            mu = self.element_class(self, mu)
        except ValueError:
            return False
        return mu in self

    def __iter__(self):
        r"""
        Iterate over ``self``.

        EXAMPLES::

            sage: KleshchevPartitions(2,[0,1],size=0)[:]
            [([], [])]
            sage: KleshchevPartitions(2,[0,1],size=1)[:]
            [([1], []), ([], [1])]
            sage: KleshchevPartitions(2,[0,1],size=2)[:]
            [([1], [1]), ([], [1, 1])]
            sage: KleshchevPartitions(3,[0,1,2],size=2)[:]
            [([1], [1], []), ([1], [], [1]), ([], [1, 1], []), ([], [1], [1]), ([], [], [2]), ([], [], [1, 1])]
        """
        size = 0
        while True:
            for mu in KleshchevPartitions_size(self._e, self._multicharge,
                                               size, self._convention):
                yield self(mu)
            size += 1

    def _an_element_(self):
        """
        Return a generic element.

        FIXME: This is a test of the wrong class!

        EXAMPLES::

            sage: KleshchevPartitions(3, [0,0,0,0], size=4).an_element()
            ([1], [1], [1], [1])
        """
        return self[12]

class KleshchevPartitions_size(KleshchevPartitions):
    """
    Kleshchev partitions of a fixed size.
    """
    def __init__(self, e, multicharge=(0,), size=0, convention='RS'):
        r"""
        Initialize ``self``.

        EXAMPLES::

            sage: K = KleshchevPartitions(4, 2)
            sage: TestSuite(K).run()
            sage: K = KleshchevPartitions(4, [0,2,1], 4)
            sage: TestSuite(K).run()
        """
        self._level = len(multicharge)
        if self._level == 1:
            self.Element = KleshchevPartition
            self._element_constructor_ = Partitions._element_constructor_.__func__
        else:
            self.Element = KleshchevPartitionTuple
        super(KleshchevPartitions_size, self).__init__(category=FiniteEnumeratedSets())
        self._size = size
        # As lists do not take negative indices the case e=0 needs to be handled
        # differently. Rather than doing this we set e equal to a "really big"
        # number. Mathematically, this is equivalent and it means that we don't
        # have an exception to cater for.
        self._e = e
        self._I = IntegerModRing(e)
        self._multicharge = tuple(self._I(m) for m in multicharge)
        self._convention = convention


    def _repr_(self):
        """
        EXAMPLES::

            sage: KleshchevPartitions(4, [0,0], 3)
            Kleshchev partitions with e=4 and multicharge=(0,0) and size 3
        """
        if self._level == 1:
            return 'Kleshchev partitions with e=%s and size %s' % (self._e, self._size)
        else:
            return 'Kleshchev partitions with e=%s and multicharge=(%s) and size %s' % (
                        self._e,','.join('%s'%m for m in self._multicharge), self._size) 

    def __contains__(self, mu):
        """
        Check if ``mu`` is in ``self``.

        TESTS::

            sage: PartitionTuple([[3,2],[2]]) in KleshchevPartitions(2,[0,0],7)
            False
            sage: PartitionTuple([[3,2],[],[],[],[2]]) in KleshchevPartitions(5,[0,0,0,0,0],7)
            False
            sage: PartitionTuple([[2,1],[],[1,1],[],[2]]) in KleshchevPartitions(5,[0,0,0,0,0],7)
            False
            sage: PartitionTuple([[2,1],[],[1,1],[],[3]]) in KleshchevPartitions(2,[0,0,0,0,0],9)
            False
            sage: all(mu in PartitionTuples(3,8) for mu in KleshchevPartitions(0,[0,0,0],8))
            True
        """
        if isinstance(mu, (KleshchevPartition, KleshchevPartitionTuple)):
            if not (mu.level() == self._level and mu.size() == self._size):
                return False
            mu = self.element_class(self, list(mu))
            if self._convention[1] == 'G':
                return mu.is_regular()
            else:
                return mu.is_restricted()

        try:
            mu = self.element_class(self, mu)
        except ValueError:
            return False
        return mu in self

    def __iter__level_one(self):
        r"""
        Iterate over all Kleshchev partitions of level one and a fixed size.

        EXAMPLES::

            sage: KleshchevPartitions(2,0)[:] #indirect doctest
            [[]]
            sage: KleshchevPartitions(2,1)[:] #indirect doctest
            [[1]]
            sage: KleshchevPartitions(2,2)[:] #indirect doctest
            [[1, 1]]
            sage: KleshchevPartitions(3,2)[:] #indirect doctest
            [[2], [1, 1]]
        """
        if self._size == 0:
            yield self.element_class(self, [])
        else:
            if self._e == 0:
                P = Partitions(self._size)
            elif self._convention[1] == 'G':
                P = Partitions(self._size, regular=self._e)
            else:
                P = Partitions(self._size, restricted=self._e)

            for mu in P:
                yield self.element_class(self, list(mu))

    def __iter__higher_levels(self):
        r"""
        Iterate over all KleshchevPartitions of a fixed level greater than 1
        and a fixed size.

        EXAMPLES::

            sage: KleshchevPartitions(2,[0,0],1)[:] #indirect doctest
            [([], [1])]
            sage: KleshchevPartitions(2,[0,0],2)[:] #indirect doctest
            [([1], [1]), ([], [1, 1])]
            sage: KleshchevPartitions(3,[0,0],2)[:] #indirect doctest
            [([1], [1]), ([], [2]), ([], [1, 1])]
        """
        if self._size == 0:
            yield self.element_class(self, [[]]*len(self._multicharge))
            return

        # For higher levels we have to recursively construct the restricted partitions
        # by adding on co-good nodes to smaller restricted partition. To avoid over 
        # counting we return a new restricted partition only if we added on its lowest
        # good node.
        for mu in KleshchevPartitions_size(self._e, self._multicharge,
                                           size=self._size-1,
                                           convention=self._convention):
            for cell in mu.cogood_cells().values():
                if cell is not None:
                    nu = self.element_class(self, mu.add_cell(*cell))
                    good_cells = nu.good_cells().values()
                    if self._convention[1]=="S":
                        if all(cell >= c for c in good_cells if c is not None):
                            yield nu
                    else:
                        if all(cell <= c for c in good_cells if c is not None):
                            yield nu

    @lazy_attribute
    def __iter__(self):
        """
        Wrapper to return the correct iterator which is different for
        :class:`Partitions` (level 1) and for :class:PartitionTuples`
        (higher levels).

        EXAMPLES::

            sage: KleshchevPartitions(3, 3, convention='RS')[:]
            [[2, 1], [1, 1, 1]]
            sage: KleshchevPartitions(3, 3, convention='RG')[:]
            [[3], [2, 1]]
            sage: KleshchevPartitions(3, [0], 3)[:]
            [[2, 1], [1, 1, 1]]
            sage: KleshchevPartitions(3, [0,0], 3)[:]
            [([1], [2]), ([], [2, 1]), ([1], [1, 1]), ([], [1, 1, 1])]
        """
        if self.level() == 1:
            return self.__iter__level_one
        else:
            return self.__iter__higher_levels

    def _an_element_(self):
        """
        Return a generic element.

        EXAMPLES::

            sage: KleshchevPartitions(4, [0,0,0,0], 4).an_element()
            ([1], [1], [1], [1])
            sage: KleshchevPartitions(4, [2], 4).an_element()
            [3, 1]
        """
        return self[0]

    Element = KleshchevPartitionTuple

#--------------------------------------------------
# helper functions
#--------------------------------------------------

def a_good_cell(kpt, multicharge, convention):
    from collections import defaultdict
    # We use a dictionary for the normal nodes as the indexing set is Z when e=0
    carry = defaultdict(int)        # a tally of #(removable cells)-#(addable cells)
    ret = None

    if convention[0] == 'L':
        rows = [(k,r) for k,part in enumerate(kpt) for r in range(len(part)+1)]
    else:
        rows = [(k,r) for k,part in reversed(list(enumerate(kpt))) for r in range(len(part)+1)]
    if convention[1] == 'S':
        rows.reverse()

    for row in rows:
        k,r = row
        if r == len(kpt[k]): # addable cell at bottom of a component
            carry[multicharge[k]-r] += 1
        else:
            res = multicharge[k] + kpt[k][r] - r - 1
            if r == len(kpt[k])-1 or kpt[k][r] > kpt[k][r+1]: # removable cell
                if carry[res] == 0:
                    ret = (k, r, kpt[k][r]-1)
                else:
                    carry[res] -= 1
            if r == 0 or kpt[k][r-1] > kpt[k][r]:             # addable cell
                carry[res+1] += 1

    # finally return the result
    return ret

def is_regular(kpt, multicharge, convention):
    """
    Return ``True`` if ``kpt`` is a ``multicharge``-regular partition tuple.

    A partition tuple is `e`-regular if we can get to the
    empty partition tuple by successively removing a sequence
    of good cells in the down direction.

    EXAMPLES::

        sage: from sage.combinat.partition_kleshchev import is_regular
        sage: I2 = IntegerModRing(2)
        sage: is_regular([[3,1,1], [3,2]], [I2(0),I2(2)], 'LS')
        False
        sage: I3 = IntegerModRing(3)
        sage: is_regular([[3,1,1], [3,2]], [I3(0),I3(2)], 'LS')
        True
        sage: is_regular([[], []], [I3(0),I3(2)], 'LS')
        True
    """
    if all(part == [] for part in kpt):
        return True
    convention = convention[0] + 'G'
    cell = a_good_cell(kpt, multicharge, convention)
    while cell is not None:
        k,r,c = cell
        if kpt[k][r] == 1:
            kpt[k].pop()
        else:
            kpt[k][r] -= 1
        cell = a_good_cell(kpt, multicharge, convention)
    return all(part == [] for part in kpt)

def is_restricted(kpt, multicharge, convention):
    """
    Return ``True`` if ``kpt`` is an ``multicharge``-restricted partition tuple.

    A partition tuple is `e`-restricted if we can get to the
    empty partition tuple by successively removing a sequence
    of good cells in the up direction.

    EXAMPLES::

        sage: from sage.combinat.partition_kleshchev import is_restricted
        sage: I2 = IntegerModRing(2)
        sage: is_restricted([[3,2,1], [3,1,1]], [I2(0),I2(2)], 'RS')
        False
        sage: I3 = IntegerModRing(3)
        sage: is_restricted([[3,2,1], [3,1,1]], [I3(0),I3(2)], 'RS')
        True
        sage: is_restricted([[], []], [I3(0),I3(2)], 'RS')
        True
    """
    if all(part == [] for part in kpt):
        return True
    convention = convention[0] + 'S'
    cell = a_good_cell(kpt, multicharge, convention)
    while cell is not None:
        k,r,c = cell
        if kpt[k][r] == 1:
            kpt[k].pop()
        else:
            kpt[k][r] -= 1
        cell = a_good_cell(kpt, multicharge, convention)
    return all(part == [] for part in kpt)
