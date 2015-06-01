r"""
Common Sequences


AUTHORS:

- Daniel Krenn (2015-05-22): initial version

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.

Functions and methods
---------------------

"""
#*****************************************************************************
#       Copyright (C) 2015 Daniel Krenn <dev@danielkrenn.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************

import itertools
import sage


def subsequence_by_indices(sequence, indices):
    r"""
    Returns the subsequence of the given sequence determined by the
    specified indices.

    INPUT:

    - ``sequence`` -- an iterable.

    - ``indices`` -- an iterable of increasing non-negative integers.

    OUTPUT:

    An iterator.

    EXAMPLES::

        sage: from sage.combinat.sequences import subsequence_by_indices
        sage: it = xsrange(10, 20)
        sage: tuple(subsequence_by_indices(it, xsrange(0, 10, 2)))
        (10, 12, 14, 16, 18)

    TESTS::

        sage: it = xsrange(10, 12)
        sage: tuple(subsequence_by_indices(it, xsrange(0, 10, 2)))
        (10,)
    """
    it = enumerate(sequence)
    for i in indices:
        for n, s in it:
            if n == i:
                yield s
                break
        else:
            break

