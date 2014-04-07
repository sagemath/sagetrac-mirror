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

from sage.structure.sage_object import SageObject





#*****************************************************************************


class CombinatorialStructure(SageObject):
    """
    Abstact base class.

    """
    def __init__(self, size=None):
        """

        """
        self.structure = None
        self.size = size

    def __iter__(self):
        """
        Returns an iterator of this combinatorial structure.
        """
        raise NotImplementedError

    def random_element(self):
        """
        Returns a random element of this combinatorial structure.
        """
        raise NotImplementedError


# be aware: there is sage.combinat.binary_tree.BinaryTrees
class PlaneBinaryTrees(CombinatorialStructure):
    """

    
    EXAMPLES::

        sage: from sage.combinat.combinatorial_structures import PlaneBinaryTrees
        sage: len([T for T in PlaneBinaryTrees(4)])  # not tested
        

    """
    def __init__(self, size=None):
        """

        """
        super(PlaneBinaryTrees, self).__init__(size=size)

        # T = z + z T^2
        T = CSConstructionUnlabelled()
        Z = CSAtomUnlabelled()
        T.assign(disjoint_union(Z, cartesian_product(Z, T, T)))

        self.structure = T

class NonAdjacentForms(CombinatorialStructure):
    def __init__(self):
        super(PlaneBinaryTrees, self).__init__(size=size)

        # NAF = (0 + P0 + M0)* (P + M + e)  with (P = 1, M = -1)
        zero = CSAtomUnlabelled('0')
        pone = CSAtomUnlabelled('P')
        mone = CSAtomUnlabelled('M')
        empty = CSEmptyUnlabelled('')
        NAF = CSConstructionUnlabelled((pone + mone + empty) \
                                       * sequence(zero + zero*pone + zero*mone))

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

            sage: TODO  # not tested
        """
        self.assign(*operands)

    def assign(self, *operands):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self.operands = operands
        
    def add_operand(self, operand):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self.operands.append(operand)
    
    def __iter__(self):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
 
    def __add__(self, other):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
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

            sage: TODO  # not tested
        """
        if is_CSBase(other):
            return CSCartesianProduct(self, other)
        else:
            raise TypeError, "Operation not supported."

class CSFlavor(SageObject):
    pass

class CSUnlabelled(CSFlavor):
    pass

class CSLabelled(CSFlavor):
    pass

#*****************************************************************************

def construction(*args, **kwargs):
    return CSConstructionUnlabelled(*args, **kwargs)  # TODO

class CSConstruction(CSBase):
    pass

class CSConstructionUnlabelled(CSBase, CSUnlabelled):
    pass

class CSConstructionLabelled(CSBase, CSLabelled):
    pass

#*****************************************************************************

# not sure if needed
class CSFiniteSetBase(CSBase):
    pass

class CSFiniteSetUnlabelled(CSBase, CSUnlabelled):
    pass

class CSFiniteSetLabelled(CSBase, CSLabelled):
    pass

#*****************************************************************************

def singleton(*args, **kwargs):
    return CSSingletonUnlabelled(*args, **kwargs)  # TODO

class CSSingletonBase(CSFiniteSetBase):
    """

    """
    def __init__(self, singleton, size):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton, size)
        
    def _set_singleton_(self, singleton, size):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self.singleton = singleton
        self.size = size

class CSSingletonLabelled(CSSingletonBase, CSFiniteSetLabelled):
    pass

class CSSingletonUnlabelled(CSSingletonBase, CSFiniteSetUnlabelled):
    pass

#*****************************************************************************

class CSEmptyBase(CSSingletonBase):
    def __init__(self, empty=[]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=empty, size=0)

class CSEmptyUnlabelled(CSEmptyBase, CSSingletonUnlabelled):
    pass

class CSEmptyLabelled(CSEmptyBase, CSSingletonLabelled):
    pass

#*****************************************************************************

class CSAtomBase(CSSingletonBase):
    def __init__(self, atom=[[]]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=atom, size=1)

class CSAtomUnlabelled(CSAtomBase, CSSingletonUnlabelled):
    pass

class CSAtomLabelled(CSAtomBase, CSSingletonLabelled):
    pass

#*****************************************************************************

def disjoint_union(*args, **kwargs):
    return CSDisjointUnionUnlabelled(*args, **kwargs)  # TODO

class CSDisjointUnionBase(CSBase):
    pass

class CSDisjointUnionUnlabelled(CSDisjointUnionBase, CSUnlabelled):
    pass

class CSDisjointUnionLabelled(CSDisjointUnionBase, CSLabelled):
    pass

#*****************************************************************************

def cartesian_product(*args, **kwargs):
    return CSCartesianProductUnlabelled(*args, **kwargs)  # TODO

class CSCartesianProductBase(CSBase):
    pass

class CSCartesianProductUnlabelled(CSCartesianProductBase, CSUnlabelled):
    pass

class CSCartesianProductLabelled(CSCartesianProductBase, CSLabelled):
    pass

#*****************************************************************************

def sequence(*args, **kwargs):
    return CSSequenceUnlabelled(*args, **kwargs)  # TODO

class CSSequenceBase(CSBase):
    pass

class CSSequenceUnlabelled(CSSequenceBase, CSUnlabelled):
    pass

class CSSequenceLabelled(CSSequenceBase, CSLabelled):
    pass

#*****************************************************************************

def multi_set(*args, **kwargs):
    return CSMultiSetUnlabelled(*args, **kwargs)  # TODO

class CSMultiSetBase(CSBase):
    pass

class CSMultiSetUnlabelled(CSMultiSetBase, CSUnlabelled):
    pass

class CSMultiSetLabelled(CSMultiSetBase, CSLabelled):
    pass


#*****************************************************************************

