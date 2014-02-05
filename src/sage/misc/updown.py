"""
A simple bi-directional directory.

AUTHOR: 

- Martin Albrecht (2010-11, initial version)

.. note::

    The implementation is not very efficient, especially with respect to memory.
"""

##############################################################################
#       Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

class UpDown:
    """
    A simple bi-directional dictionary.

    EXAMPLE::

        sage: from sage.misc.updown import UpDown
        sage: d = UpDown()
        sage: d['a'] = 'b'
        sage: d.down('a')
        'b'
        sage: d.up('b')
        'a'

    .. note::

        Ambigious methods such as :meth:`keys` are biased towards the
        down mapping. For example, :meth:`keys` returns the keys of
        the downward mapping.
    """
    def __init__(self):
        """
        Create a new bi-directional dictionary.

        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d
            down: {}
              up: {}
        """
        self._down = {}
        self._up = {}

    def __setitem__(self, k, v):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.down(0)
            1
            sage: d.up(1)
            0
        """
        self._down[k] = v 
        self._up[v] = k

    def __call__(self, k, v):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d(0,1)
            sage: d.down(0)
            1
            sage: d.up(1)
            0
        """
        self[k] = v

    def __len__(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: len(d)
            1
        """
        return len(self._down)

    def down(self, k):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.down(0)
            1
        """
        return self._down[k]

    def up(self, v):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.up(1)
            0
        """
        return self._up[v]

    def __repr__(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d #indirect doctest
            down: {0: 1}
              up: {1: 0}
        """
        return "down: %s\n  up: %s"%(str(self._down), str(self._up))
    
    def __contains__(self, x):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: 0 in d
            True
            sage: 1 in d
            True
        """
        return x in self._down or x in self._up

    def keys(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.keys()
            [0]
        """
        return self._down.keys()

    def iterkeys(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.iterkeys()
            <dictionary-keyiterator object at 0x...>
        """
        return self._down.iterkeys()

    def values(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.values()
            [1]
        """
        return self._down.values()

    def itervalues(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.itervalues()
            <dictionary-valueiterator object at 0x...>
        """
        return self._down.itervalues()

    def iteritems(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.iteritems()
            <dictionary-itemiterator object at 0x...>
        """
        return self._down.iteritems()
    
    def items(self):
        """
        EXAMPLE::

            sage: from sage.misc.updown import UpDown
            sage: d = UpDown()
            sage: d[0] = 1
            sage: d.items()
            [(0, 1)]
        """
        return self._down.items()
 
