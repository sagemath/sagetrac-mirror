r"""
Common Sequences

Sequences in SageMath can be built through the ``sequences``
object, which contains common sequences. For example,

::

    sage: sequences.catalan()
    catalan sequence 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, ...

generates the sequence of catalan numbers; the first `10` are shown as
a preview.

**Sequences**

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~Sequences.catalan` | the sequence of catalan numbers
    :meth:`~Sequences.fibonacci` | the sequence of fibonacci numbers


Various
=======

AUTHORS:

- Daniel Krenn (2015, 2016)


ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the Austrian Science Fund (FWF): P 24644-N26.


Classes and Methods
===================
"""
# *****************************************************************************
# Copyright (C) 2015, 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_function
from sage.structure.homogenous_sequence import Sequence


class Sequences(object):
    r"""
    A class consisting of constructors for several common sequences.

    A list of all sequences in this database is available via tab
    completion. Type "``sequences.``" and then hit tab to see which
    are there.

    The sequences currently in this class include:

    - :meth:`~catalan`
    - :meth:`~fibonacci`
    """

    @staticmethod
    @cached_function
    def catalan():
        r"""
        The sequence of Catalan numbers.

        INPUT:

        Nothing.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: sequences.catalan()
            catalan sequence 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, ...
        """
        from sage.combinat.combinat import catalan_number
        from sage.rings.integer_ring import ZZ

        return Sequence(
            lambda n: catalan_number(n),
            universe=ZZ, name='catalan sequence')


    @staticmethod
    @cached_function
    def fibonacci(algorithm=None):
        r"""
        The sequence of Fibonacci numbers.

        INPUT:

        -  ``algorithm`` -- (default: ``None``) passed on to the
           :func:`Fibonacci function <sage.combinat.combinat.fibonacci>`.
           If ``None``, the default algorithm is used.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: sequences.fibonacci()
            fibonacci sequence 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, ...

        ::

            sage: tuple(sequences.fibonacci()[10:20])
            (55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181)

        ::

            sage: sum(sequences.fibonacci()[100:110])
            69919376923075308730013

        TESTS::

            sage: sequences.fibonacci() is sequences.fibonacci()
            True
        """
        from sage.combinat.combinat import fibonacci
        from sage.rings.integer_ring import ZZ
        if algorithm is None:
            algorithm = 'pari'  # TODO: do it in sage.combinat.combinat
        return Sequence(
            lambda n: fibonacci(n, algorithm=algorithm),
            universe=ZZ, name='fibonacci sequence')


# Easy access to the sequence generators:
sequences = Sequences()
r"""

EXAMPLES::

    sage: sequences.fibonacci()
    fibonacci sequence 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, ...
"""
