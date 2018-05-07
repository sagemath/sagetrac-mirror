from sage.all import *
from sage.combinat.partition import _Partitions
from sage.combinat.skew_partition import SkewPartition, SkewPartitions
from sage.structure.unique_representation import CachedRepresentation

def row_col_to_skew_partition(rs, cs):
    """
    matt
    """
    outer = []
    inner = []
    current_cs = [0] * len(cs)
    row_index = 0
    for col_coindex, col_length in enumerate(list(reversed(cs))):
        current_col_length = list(reversed(current_cs))[col_coindex]
        num_rows_to_slide = col_length - current_col_length
        if num_rows_to_slide < 0:
            raise ValueError('The inputted (row-shape, col-shape) pair has no possible corresponding skew-shape.')
        # 'col_num' is 1-based index of cols
        col_num = len(cs) - col_coindex
        while num_rows_to_slide > 0:
            if row_index > len(rs) - 1:
                raise ValueError('error more')
            # slide a row
            outer.append(col_num)
            inner.append(col_num - rs[row_index])
            # update params/info
            for c in range(col_num - rs[row_index], col_num):
                current_cs[c] += 1
            row_index += 1
            num_rows_to_slide -= 1
    return SkewPartition([outer, inner])

class kBoundary (SkewPartition, CachedRepresentation):
    """
    matt
    Given a partition l and a positive integer k, the __k-boundary__ of l is the skew-shape obtained from the shape of l by removing all cells of hook-length greater than k.
    """
    @staticmethod
    def __classcall_private__(cls, l, k):
        """
        matt
        Normalize input to ensure unique representation.
        """
        l = Partition(l)
        k = NN(k)
        # this BELOW should fail because SkewPartition(l, k) would fail.
        return super(kBoundary, cls).__classcall__(cls, l, k)

    def __init__(self, l, k):
        # NOTE: THIS FUNCTION IS REDUNDANT WITH Permutation.k_boundary AND THIS SHOULD BE ADDRESSED
        """
        matt
        l: the partition
        k: the largest allowed hook length
        """
        # to make a SkewPartition, we must calculate the inner and outer partitions
        outer_partition = l.to_list()
        # we could make more efficient, but it's simple to get the tableau of all hook lengths and pick out the inner partition from it
        inner_partition = []
        for row_hook_lengths in l.hook_lengths():
            inner_row_length = len([x for x in row_hook_lengths if x > k])
            inner_partition.append(inner_row_length)
        # finally, make the SkewPartition
        SkewPartition.__init__(self, [outer_partition, inner_partition])

    def partition(self):
        """matt
        Return the partition whose k-boundary is self.
        """
        return k_boundary_to_partition(self, strict=False)

    # def __eq__(self):

"""
Given a skew-linked diagram, is it a k-boundary?  (That is, does there exist some partition which - when cells of hook-length > k are removed - becomes the skew-linked diagram.)
"""

def k_boundary_to_partition(skew_shape, strict=True):
    """
    matt
    skew_shape: The skew-shape (a k-boundary) to find the corresponding partition for.
    strict: If true, assert that the skew-shape really is a k-boundary.
    """
    if strict:
        assert is_k_boundary(skew_shape)
    return skew_shape.minimal_containing_partition()
