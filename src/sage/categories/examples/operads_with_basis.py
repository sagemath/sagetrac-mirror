r"""
Examples of operads with basis
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Nicolas M. Thiery <nthiery at users.sf.net>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.sets.family import Family
from sage.categories.all import OperadsWithBasis
from sage.combinat.free_module import CombinatorialFreeModule
from sage.combinat.words.words import Words

class AssociativeOperad(CombinatorialFreeModule):
    r"""
    An example of an operad with basis: the Associative operad

    This class illustrates a minimal implementation of an operad with basis.
    """

    def __init__(self, R):
        """
        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example(); A
            An example of an operad with basis: the Associative operad over Rational Field
            sage: TestSuite(A).run()

        """
        CombinatorialFreeModule.__init__(self, R, Words(), category = OperadsWithBasis(R))

    def _repr_(self):
        """
        EXAMPLES::

            sage: OperadsWithBasis(QQ).example() # indirect doctest
            An example of an operad with basis: the Associative operad over Rational Field
        """
        return "An example of an operad with basis: the Associative operad over %s"%(self.base_ring())

    @cached_method
    def one_basis(self,letter):
        """
        Returns the word of length one, which index the one of this operad,
        as per :meth:`OperadsWithBasis.ParentMethods.one_basis`.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: A.one_basis("a")
            word: a
        """
        return self.basis().keys()([letter])

    def map_labels(self,t,f):
        """
        Maps the function `f` on the word `t`.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: m = Words([4,3,2,1])
            sage: A.map_labels(m,lambda u:u)
            word: 4321
        """
        return self.basis().keys()([f(u) for u in t])

    def grafts(self,x,y,i):
        """
        Auxiliary procedure: inserts a word y at position i in a word x
        and returns a word

        This is the composition of the set-theoretic Associative operad.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: A.grafts(Words("acb"), Words("de"),"c")
            word: adeb

        """
        if x[0]==i:
            return y+x[1:]
        else:
            return x[:1]+self.grafts(x[1:],y,i)

    def composition_on_basis(self,x,y,i):
        """
        Composition of basis elements, as per :meth:`OperadsWithBasis.ParentMethods.composition_on_basis`.

        EXAMPLES::

            sage: A = OperadsWithBasis(QQ).example()
            sage: Words = A.basis().keys()
            sage: A.composition_on_basis(Words("acb"), Words("de"),"c")
            B[word: adeb]
        """
        if not(i in x):
            return "The composition index is not present in the first argument."
        else:
            return self.basis()[self.grafts(x,y,i)]


Example = AssociativeOperad
