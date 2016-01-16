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

