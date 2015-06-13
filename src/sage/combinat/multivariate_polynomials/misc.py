"""
Misc classes used for multivariate polynomials
"""
#*****************************************************************************
#  Copyright (C) 2008-2009 Viviane Pons < pons at univ-mlv.fr  >
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.free_module import CombinatorialFreeModule


from sage.rings.rational_field import QQ

class MonomialKeyWrapper(Parent, UniqueRepresentation):

    def __init__(self, root_system = None, length = None):
        self._length = length
        self._root_system = root_system

        if length is None and root_system is None:
            raise ValueError("Either length or root_system must be initialize")

        if not root_system is None:
            self._is_typed = True
            self.Element = self.RootSystemWrapper
        else:
            self._is_typed = False
            self.Element = self.SimpleWrapper

    def is_typed(self):
        return self._is_typed

    def length(self):
        return self._length

    def root_sytem(self):
        return


    class SimpleWrapper(Element):

        def __init__(self, parent, wrapped):
            Element.__init__(parent = parent)
            self._wrapped = wrapped

        def __iter__(self):

            for i in self._wrapped:
                yield i

        def __getitem__(self, key):
            return self._wrapped[key]



    class RootSystemWrapper(Element):
        pass

# NT, VP : Any suggestion for the place where this code should go ?
class RelativeIntegerVectors(CombinatorialFreeModule):
    r"""
    """
    def __init__(self, length):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        self._length = length
        CombinatorialFreeModule.__init__(self, QQ, xrange(length))

    def __getitem__(self, x):
        r"""
        TESTS::

            sage: # Fix a nice test
        """
        if( hasattr(x,'__iter__') ):
            return self(sum( [ x[i]*self.monomial(i)  for i in range( len( list(x)) ) ] ))
        else:
            return x * self.monomial(0)

    def __call__(self,x):
        if( type(x) is list or type(x) is tuple):
            return self.__getitem__(x)
        return super(RelativeIntegerVectors,self).__call__(x)

    def length(self):
        r"""
        Returns the length of ``self``.

        EXAMPLES::

            sage: # Fix a nice example
        """
        return self._length

    class Element(CombinatorialFreeModule.Element):
        def coeffs_to_list(self):
            r"""
            TESTS::

                sage: # Fix a nice test
            """
            return [ self[i] for i in xrange(self.parent().length()) ]

        def _repr_(self):
            r"""
            TESTS::

                sage: # Fix a nice test
            """
            return str(self.coeffs_to_list())


