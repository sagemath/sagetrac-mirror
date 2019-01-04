r"""
Partition/Diagram Algebras
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
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
from __future__ import absolute_import

#####
#A_k#
#####

def SetPartitionsAk(k):
    r"""
    ``SetPartitionsAk`` is deprecated; use :class:`PartitionDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `A_k`.

    EXAMPLES::

        sage: A3 = SetPartitionsAk(3); A3
        doctest:...: DeprecationWarning: SetPartitionsAk is deprecated; use PartitionDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Partition diagrams of order 3

        sage: A3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: A3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: A3.random_element()  #random
        {{1, 3, -3, -1}, {2, -2}}

        sage: A3.cardinality()
        203

        sage: A2p5 = SetPartitionsAk(2.5); A2p5
        Partition diagrams of order 5/2
        sage: A2p5.cardinality()
        52

        sage: A2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: A2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: A2p5.random_element() #random
        {{-1}, {-2}, {3, -3}, {1, 2}}

    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsAk is deprecated; use PartitionDiagrams instead")
    from sage.combinat.diagram_algebras import PartitionDiagrams
    return PartitionDiagrams(k)

#####
#S_k#
#####

def SetPartitionsSk(k):
    r"""
    ``SetPartitionsSk`` is deprecated; use
    :class:`PermutationDiagrams` instead.

    Return the combinatorial class of set partitions of type `S_k`.

    There is a bijection between these set partitions and the
    permutations of `1, \ldots, k`.

    EXAMPLES::

        sage: S3 = SetPartitionsSk(3); S3
        doctest:...: DeprecationWarning: SetPartitionsSk is deprecated; use PermutationDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Permutation diagrams of order 3
        sage: S3.cardinality()
        6

        sage: S3.list()  #random
        [{{2, -2}, {3, -3}, {1, -1}},
         {{1, -1}, {2, -3}, {3, -2}},
         {{2, -1}, {3, -3}, {1, -2}},
         {{1, -2}, {2, -3}, {3, -1}},
         {{1, -3}, {2, -1}, {3, -2}},
         {{1, -3}, {2, -2}, {3, -1}}]
        sage: S3.first() #random
        {{2, -2}, {3, -3}, {1, -1}}
        sage: S3.last() #random
        {{1, -3}, {2, -2}, {3, -1}}
        sage: S3.random_element() #random
        {{1, -3}, {2, -1}, {3, -2}}

        sage: S3p5 = SetPartitionsSk(3.5); S3p5
        Permutation diagrams of order 7/2
        sage: S3p5.cardinality()
        6

        sage: S3p5.list() #random
        [{{2, -2}, {3, -3}, {1, -1}, {4, -4}},
         {{2, -3}, {1, -1}, {4, -4}, {3, -2}},
         {{2, -1}, {3, -3}, {1, -2}, {4, -4}},
         {{2, -3}, {1, -2}, {4, -4}, {3, -1}},
         {{1, -3}, {2, -1}, {4, -4}, {3, -2}},
         {{1, -3}, {2, -2}, {4, -4}, {3, -1}}]
        sage: S3p5.first() #random
        {{2, -2}, {3, -3}, {1, -1}, {4, -4}}
        sage: S3p5.last() #random
        {{1, -3}, {2, -2}, {4, -4}, {3, -1}}
        sage: S3p5.random_element() #random
        {{1, -3}, {2, -2}, {4, -4}, {3, -1}}

    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsSk is deprecated; use PermutationDiagrams instead")
    from sage.combinat.diagram_algebras import PermutationDiagrams
    return PermutationDiagrams(k)

#####
#I_k#
#####

def SetPartitionsIk(k):
    r"""
    ``SetPartitionsIk`` is deprecated; use :class:`IdealDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `I_k`.

    These are set partitions with a propagating number of less than `k`.
    Note that the identity set partition `\{\{1, -1\}, \ldots, \{k, -k\}\}`
    is not in `I_k`.

    EXAMPLES::

        sage: I3 = SetPartitionsIk(3); I3
        doctest:...: DeprecationWarning: SetPartitionsIk is deprecated; use IdealDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Ideal diagrams of order 3
        sage: I3.cardinality()
        197

        sage: I3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: I3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: I3.random_element() #random
        {{-1}, {-3, -2}, {2, 3}, {1}}

        sage: I2p5 = SetPartitionsIk(2.5); I2p5
        Ideal diagrams of order 5/2
        sage: I2p5.cardinality()
        50

        sage: I2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: I2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: I2p5.random_element() #random
        {{-1}, {-2}, {1, 3, -3}, {2}}
    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsIk is deprecated; use IdealDiagrams instead")
    from sage.combinat.diagram_algebras import IdealDiagrams
    return IdealDiagrams(k)

#####
#B_k#
#####

def SetPartitionsBk(k):
    r"""
    ``SetPartitionsBk`` is deprecated; use :class:`BrauerDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `B_k`.

    These are the set partitions where every block has size 2.

    EXAMPLES::

        sage: B3 = SetPartitionsBk(3); B3
        doctest:...: DeprecationWarning: SetPartitionsBk is deprecated; use BrauerDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Brauer diagrams of order 3

        sage: B3.first() #random
        {{2, -2}, {1, -3}, {3, -1}}
        sage: B3.last() #random
        {{1, 2}, {3, -2}, {-3, -1}}
        sage: B3.random_element() #random
        {{2, -1}, {1, -3}, {3, -2}}

        sage: B3.cardinality()
        15

        sage: B2p5 = SetPartitionsBk(2.5); B2p5
        Brauer diagrams of order 5/2

        sage: B2p5.first() #random
        {{2, -1}, {3, -3}, {1, -2}}
        sage: B2p5.last() #random
        {{1, 2}, {3, -3}, {-1, -2}}
        sage: B2p5.random_element() #random
        {{2, -2}, {3, -3}, {1, -1}}

        sage: B2p5.cardinality()
        3
    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsBk is deprecated; use BrauerDiagrams instead")
    from sage.combinat.diagram_algebras import BrauerDiagrams
    return BrauerDiagrams(k)

#####
#P_k#
#####

def SetPartitionsPk(k):
    r"""
    ``SetPartitionsPk`` is deprecated; use :class:`PlanarDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `P_k`.

    These are the planar set partitions.

    EXAMPLES::

        sage: P3 = SetPartitionsPk(3); P3
        doctest:...: DeprecationWarning: SetPartitionsPk is deprecated; use PlanarDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Planar diagrams of order 3
        sage: P3.cardinality()
        132

        sage: P3.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: P3.last() #random
        {{-1}, {-2}, {3}, {1}, {-3}, {2}}
        sage: P3.random_element() #random
        {{1, 2, -1}, {-3}, {3, -2}}

        sage: P2p5 = SetPartitionsPk(2.5); P2p5
        Planar diagrams of order 5/2
        sage: P2p5.cardinality()
        42

        sage: P2p5.first() #random
        {{1, 2, 3, -1, -3, -2}}
        sage: P2p5.last() #random
        {{-1}, {-2}, {2}, {3, -3}, {1}}
        sage: P2p5.random_element() #random
        {{1, 2, 3, -3}, {-1, -2}}

    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsPk is deprecated; use PlanarDiagrams instead")
    from sage.combinat.diagram_algebras import PlanarDiagrams
    return PlanarDiagrams(k)

#####
#T_k#
#####

def SetPartitionsTk(k):
    r"""
    ``SetPartitionsTk`` is deprecated; use :class:`TemperleyLiebDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `T_k`.

    These are planar set partitions where every block is of size 2.

    EXAMPLES::

        sage: T3 = SetPartitionsTk(3); T3
        doctest:...: DeprecationWarning: SetPartitionsTk is deprecated; use TemperleyLiebDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Temperley Lieb diagrams of order 3
        sage: T3.cardinality()
        5

        sage: T3.first() #random
        {{1, -3}, {2, 3}, {-1, -2}}
        sage: T3.last() #random
        {{1, 2}, {3, -1}, {-3, -2}}
        sage: T3.random_element() #random
        {{1, -3}, {2, 3}, {-1, -2}}

        sage: T2p5 = SetPartitionsTk(2.5); T2p5
        Temperley Lieb diagrams of order 5/2
        sage: T2p5.cardinality()
        2

        sage: T2p5.first() #random
        {{2, -2}, {3, -3}, {1, -1}}
        sage: T2p5.last() #random
        {{1, 2}, {3, -3}, {-1, -2}}
    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsTk is deprecated; use TemperleyLiebDiagrams instead")
    from sage.combinat.diagram_algebras import TemperleyLiebDiagrams
    return TemperleyLiebDiagrams(k)

def SetPartitionsRk(k):
    r"""
    ``SetPartitionsRk`` is deprecated; use :class:`RookDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `R_k`.

    EXAMPLES::

        sage: SetPartitionsRk(3)
        doctest:...: DeprecationWarning: SetPartitionsRk is deprecated; use RookDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Rook diagrams of order 3
    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsRk is deprecated; use RookDiagrams instead")
    from sage.combinat.diagram_algebras import RookDiagrams
    return RookDiagrams(k)

def SetPartitionsPRk(k):
    r"""
    ``SetPartitionsPRk`` is deprecated; use :class:`PlanarRookDiagrams`
    instead.

    Return the combinatorial class of set partitions of type `PR_k`.

    EXAMPLES::

        sage: SetPartitionsPRk(3)
        doctest:...: DeprecationWarning: SetPartitionsPRk is deprecated; use PlanarRookDiagrams instead
        See https://trac.sagemath.org/25637 for details.
        Planar rook diagrams of order 3
    """
    from sage.misc.superseded import deprecation
    deprecation(25637, "SetPartitionsPRk is deprecated; use PlanarRookDiagrams instead")
    from sage.combinat.diagram_algebras import PlanarRookDiagrams
    return PlanarRookDiagrams(k)
