"""
Deprecated low-level permutations
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


def PermutationsNK(n, k):
    r"""
    This is deprecated in :trac:`16472`. Use :class:`~permutation.Permutations` instead
    (or ``itertools.permutations`` for iteration).

    EXAMPLES::

        sage: from sage.combinat.permutation_nk import PermutationsNK
        sage: P = PermutationsNK(10,4)
        doctest:...: DeprecationWarning: PermutationsNK is deprecated. Please
        use Permutations instead (or itertools.permutations for iteration).
        See http://trac.sagemath.org/16472 for details.
        sage: [ p for p in PermutationsNK(3,2)]
        [[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]]
        sage: PermutationsNK(3,2).cardinality()
        6
        sage: PermutationsNK(5,4).cardinality()
        120
    """
    from sage.misc.superseded import deprecation
    deprecation(16472, "PermutationsNK is deprecated. Please use Permutations instead (or itertools.permutations for iteration).")
    from permutation import Permutations
    return Permutations(range(n),k)

from permutation import Permutations_nk
class PermutationsNK_backward_compatibility(Permutations_nk):
    r"""
    .. WARNING::

        Not to be used! (backward compatibility for pickling)
    """
    def __getattr__(self, name):
        r"""
        If the attribute

        EXAMPLES::

            sage: from sage.combinat.permutation_nk import PermutationsNK_backward_compatibility
            sage: P = PermutationsNK_backward_compatibility(10,3)
            sage: P
            Permutations of {1,...,10} of length 3
        """
        if name == 'n' or name == 'k':
            n = self._n
            k = self._k
            Permutations_nk.__init__(self, n, k)
            return getattr(self, name)
        raise AttributeError

from sage.structure.sage_object import register_unpickle_override
register_unpickle_override(
    'sage.combinat.permutation_nk',
    'PermutationsNK',
    PermutationsNK_backward_compatibility,
    call_name=('sage.combinat.permutation', 'Permutations'))
