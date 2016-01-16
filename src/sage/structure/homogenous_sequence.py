r"""
Homogenous Sequences

An immutable sequence with elements of a common guaranteed
universe. The elements are evaluated lazily.

Various
=======

AUTHORS:

- Daniel Krenn (2016)


Classes and Methods
===================
"""

# *****************************************************************************
# Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# http://www.gnu.org/licenses/
# *****************************************************************************

from sage.misc.cachefunc import cached_function
from sage.misc.lazy_list import lazy_list, lazy_list_generic
from sage.structure.sage_object import SageObject


sequence_default_kwds = {
    'name': 'sequence',
    'opening_delimiter': '',
    'closing_delimiter': '',
    'preview': 10}

subsequence_default_prefix = 'subsequence of '


def Sequence(data=None, universe=None, convert=True,
             name=None, **kwds):
    r"""
    EXAMPLES::

        sage: from sage.structure.homogenous_sequence import Sequence
        sage: P = Sequence(Primes(), name='sequence of primes'); P
        sequence of primes 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, ...
        sage: P[1]
        3
        sage: P[5:]
        sequence 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, ...
    """
    from copy import copy

    cls_kwds = copy(sequence_default_kwds)
    if name is not None:
        cls_kwds['name'] = name
    cls_kwds['universe'] = universe
    cls_kwds['convert'] = convert

    return lazy_list(data,
                     cls=HomogenousSequence,
                     cls_kwds=cls_kwds,
                     **kwds)


class HomogenousSequence(SageObject, lazy_list_generic):

    def __init__(self, universe=None, convert=True,
                 name=None, cls_kwds=None,
                 **kwds):
        r"""

        TESTS::

            sage: from sage.structure.homogenous_sequence import Sequence
            sage: Sequence(Primes(), name='sequence of primes')._properties_()       
            {'closing_delimiter': '',
             'cls': <class 'sage.structure.homogenous_sequence.HomogenousSequence'>,
             'cls_kwds': {'closing_delimiter': '',
              'convert': True,
              'name': 'sequence',
              'opening_delimiter': '',
              'preview': 10,
              'universe': None},
             'more': '...',
             'name': 'sequence of primes',
             'opening_delimiter': '',
             'preview': 10,
             'separator': ', ',
             'start': 0,
             'step': 1,
             'stop': 9223372036854775807}
        """
        if cls_kwds is None:
            cls_kwds = {}
        lazy_list_generic.__init__(self, name=name, cls_kwds=cls_kwds, **kwds)

        cls_kwds['name'] = sequence_default_kwds['name']  # subsequence gets the default name

        self.convert = False  # prevent conversion until universe is set for sure
        if universe is None:
            universe = self[0].parent()
        self.universe = universe
        self.convert = convert


    _repr_ = lazy_list_generic.__repr__


    def update_cache_up_to(self, i):
        # The following works correctly if the cache is extended.
        # If some update_cache_up_to (e.g. the one from
        # lazy_list_from_update_function) changes the cache, then
        # a common universe is not guaranteed.
        cache = self._get_cache_()
        length = len(cache)
        result = lazy_list_generic.update_cache_up_to(self, i)
        if self.convert:
            for l in range(length, len(cache)):
                cache[l] = self.universe(cache[l])
        return result

