r"""
Examples of set operads
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import SetOperads
from sage.combinat.words.words import Words
from sage.structure.parent import Parent

class AssociativeOperad(Parent):
    r"""
    An example of a set operad: the Associative operad

    This class illustrates a minimal implementation of a set operad.
    """

    def __init__(self):
        """
        EXAMPLES::

            sage: A = SetOperads().example(); A
            An example of a set operad: the Associative operad
            sage: TestSuite(A).run()
        """
        self.element_class = Words()
        Parent.__init__(self, category = SetOperads())

    def _repr_(self):
        """
        EXAMPLES::

            sage: SetOperads().example() # indirect doctest
            An example of a set operad: the Associative operad
        """
        return "An example of a set operad: the Associative operad"

    @cached_method
    def one(self,letter):
        """
        Returns the word of length one, which index the one of this operad,
        as per :meth:`SetOperads.ParentMethods.one`.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.one("a")
            word: a
        """
        return self.element_class([letter])

    def _an_element_(self):
        """
        Returns a word

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A._an_element_()
            word: abcd
        """
        return self.element_class("abcd")

    def composition(self,x,y,i):
        """
        Composition of words, as per
        :meth:`SetOperads.ParentMethods.composition`.

        inserts a word y at position i in a word x and returns a word

        This is the composition of the set-theoretic Associative operad.

        EXAMPLES::

            sage: A = SetOperads().example()
            sage: A.composition(Word("acb"), Word("de"),"c")
            word: adeb

        """
        if x[0]==i:
            return self.element_class(y+x[1:])
        else:
            return self.element_class(x[:1]+self.composition(x[1:],y,i))

Example = AssociativeOperad
