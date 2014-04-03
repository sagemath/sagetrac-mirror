# -*- coding: utf-8 -*-
r"""

AUTHORS:

- Daniel Krenn (2014-04-01): initial version

"""

#*****************************************************************************
#       Copyright (C) 2014 Daniel Krenn <devel@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************






#*****************************************************************************


class CombinatorialStructure(SageObject):
    """
    Abstact base class.

    """
    def __init__(self):
        """

        """
        self.structure = None

    def __iter__(self):
        """
        Returns an iterator of this combinatorial structure.
        """


# be aware: there is sage.combinat.binary_tree.BinaryTrees
class PlaneBinaryTrees(CombinatorialStructure):
    """

    
    EXAMPLES::

        sage: len([T for T in PlaneBinaryTrees(4)])
        

    """
    def __init__(self):
        """

        """
        # T = z + z T^2
        T = CSConstruction()
        T.assign(CSDisjointUnion(CSAtom(), CSCartesianProduct(CSAtom(), T, T)))

        self.structure = T

class NonAdjacentForms(CombinatorialStructure):
    def __init__(self):
        # NAF = (0 + P0 + M0)* (P + M + e)  with (P = 1, M = -1)
        zero = CSAtom('0')
        pone = CSAtom('P')
        mone = CSAtom('M')
        empty = CSEmpty('')
        NAF = CSConstruction((pone + mone + empty) \
                                 CSSequence(zero + zero*pone + zero*mone))

        self.structure = NAF

#*****************************************************************************
# Data Structure
#*****************************************************************************

# TODO: find suitable name instead of CSBase
class CSBase(SageObject):
    """
    Abstact base class.

    """
    def __init__(self, *operands):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self.assign(*operands)

    def assign(self, *operands):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self.operands = operands
        
    def add_operand(self, operand):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self.operands.append(operand)
    
    def __iter__(self):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
 
    def __add__(self, other):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        if is_CSBase(other):
            #if is_CSDisjointUnion(self):
            #    # make copy here and add
                
            return CSDisjointUnion(self, other)
        else:
            raise TypeError, "Operation not supported."

    def __mul__(self, other):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        if is_CSBase(other):
            return CSCartesianProduct(self, other)
        else:
            raise TypeError, "Operation not supported."
            
#*****************************************************************************

class CSConstruction(CSBase)

#*****************************************************************************

# not sure if needed
class CSFiniteSet(CSBase)

#*****************************************************************************

class CSSingleton(CSFiniteSet):
    """

    """
    def __init__(self, singleton, size):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self._set_singleton_(singleton, size)
        
    def _set_singleton_(self, singleton, size):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self.singleton = singleton
        self.size = size

#*****************************************************************************

class CSEmpty(CSSingleton):
    def __init__(self, empty=[]):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self._set_singleton_(singleton=empty, size=0)

#*****************************************************************************

class CSAtom(CSSingleton)
    def __init__(self, atom=[[]]):
        """
        TODO

        EXAMPLES::

            sage: TODO
        """
        self._set_singleton_(singleton=atom, size=1)

#*****************************************************************************

class CSDisjointUnion(CSBase)

#*****************************************************************************

class CSCartesianProduct(CSBase)

#*****************************************************************************

class CSSequence(CSBase)

#*****************************************************************************

class CSMultiSet(CSBase)


#*****************************************************************************

