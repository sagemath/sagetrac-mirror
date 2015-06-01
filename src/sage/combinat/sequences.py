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


class subsequence_options(object):
    def __init__(self):
        """
        This decorator of a function (which returns an iterator) adds
        options for extracting subsequences.

        INPUT:

        Nothing.

        EXAMPLES::

            sage: from sage.combinat.sequences import subsequence_options
            sage: from itertools import count
            sage: @subsequence_options()
            ....: def cnt():
            ....:     return count()
            sage: tuple(cnt(stop=3))
            (0, 1, 2)

        ::

            sage: @cached_function
            ....: def f(n):
            ....:     if n == 0 or n == 1:
            ....:         return n
            ....:     return f(n-1) + f(n-2)
            sage: @subsequence_options()
            ....: def seq():
            ....:     return iter(f(n) for n in count())
            sage: tuple(seq(stop=10))
            (0, 1, 1, 2, 3, 5, 8, 13, 21, 34)
        """
        pass


    def __call__(self, func):
        """
        Returns a wrapper around ``func``

        INPUT:

        - ``func`` -- a function.

        OUTPUT:

        A function.

        TESTS::

            sage: from sage.combinat.sequences import subsequence_options
            sage: @subsequence_options()
            ....: def s():
            ....:     return xsrange(7)
            sage: tuple(s())
            (0, 1, 2, 3, 4, 5, 6)
            sage: tuple(s(start=2))
            (2, 3, 4, 5, 6)
            sage: tuple(s(stop=4))
            (0, 1, 2, 3)
            sage: tuple(s(start=3, stop=5))
            (3, 4)
            sage: tuple(s(start=1, stop=6, step=2))
            (1, 3, 5)
            sage: tuple(s(take=is_even))
            (0, 2, 4, 6)
            sage: tuple(s(drop_until=lambda n: n >= 2))
            (2, 3, 4, 5, 6)
            sage: tuple(s(take_while=lambda n: n < 3))
            (0, 1, 2)
        """
        from sage.misc.misc import xsrange
        from sage.rings.integer_ring import ZZ

        @sage.misc.decorators.sage_wraps(func)
        def subsequence(*args, **kwds):

            indices = kwds.pop('indices', None)

            if indices is not None and \
                    any(kwds.has_key(k) for k in ('start', 'stop', 'step')):
                raise ValueError('You cannot specify the parameters '
                                 'indices and '
                                 'start, stop, step at the same time.')

            if indices is None:
                start = kwds.pop('start', ZZ(0))
                stop = kwds.pop('stop', None)
                step = kwds.pop('step', ZZ(1))

                if stop is None:
                    indices = itertools.count(start=start, step=step)
                else:
                    indices = xsrange(start, stop, step)

            take = kwds.pop('take', None)
            drop_until = kwds.pop('drop_until', None)
            take_while = kwds.pop('take_while', None)

            result = subsequence_by_indices(func(*args, **kwds), indices)

            if take is not None:
                result = iter(s for s in result if take(s))

            if drop_until is not None:
                drop_while = lambda e: not drop_until(e)
                result = itertools.dropwhile(drop_while, result)

            if take_while is not None:
                result = itertools.takewhile(take_while, result)

            return result

        return subsequence


class Sequences(object):
    r"""
    A class consisting of constructors for several common sequences.

    A list of all sequences in this database is available via tab
    completion. Type "``sequences.``" and then hit tab to see which
    are there.


# Easy access to the sequence generators from the command line:
sequences = Sequences()
