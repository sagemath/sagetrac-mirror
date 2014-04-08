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





# ----------------------------------------------------------------------------


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
# Data Structures -- Base
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
                
            return disjoint_union(self, other)
        else:
            raise TypeError, "Operation not supported."

    def __mul__(self, other):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        if is_CSBase(other):
            return cartesian_product(self, other)
        else:
            raise TypeError, "Operation not supported."

#*****************************************************************************
# Data Structures -- Flavor
#*****************************************************************************

class CSFlavor(SageObject):
    pass

# ----------------------------------------------------------------------------

def has_unlabelled_flavor(CS):
    """
    Tests whether ``CS`` inherits from :class:`CSUnlabelled` or not.
    """
    return isinstance(CS, CSUnlabelled)

class CSUnlabelled(CSFlavor):
    pass

# ----------------------------------------------------------------------------

def has_labelled_flavor(CS):
    """
    Tests whether ``CS`` inherits from :class:`CSLabelled` or not.
    """
    return isinstance(CS, CSLabelled)

class CSLabelled(CSFlavor):
    pass

# ----------------------------------------------------------------------------

def _process_flavor_(kwargs):
    """
    Returns the flavor (``Labelled`` or ``Unlabelled``) encoded ``kwargs``.

    INPUT:

    - kwargs -- a dictionary

    OUTPUT:

    A string ``Labelled`` or ``Unlabelled``.

    TESTS::

        sage: from sage.combinat.combinatorial_structures import _process_flavor_
        sage: _process_flavor_({})
        'Unlabelled'
        sage: _process_flavor_({'labelled': True})
        'Labelled'
        sage: _process_flavor_({'labelled': False})
        'Unlabelled'
        sage: _process_flavor_({'unlabelled': True})
        'Unlabelled'
        sage: _process_flavor_({'unlabelled': False})
        'Labelled'
        sage: _process_flavor_({'labelled': True, 'unlabelled': True})
        Traceback (most recent call last):
        ...
        ValueError: Arguments incompatible.
        sage: _process_flavor_({'labelled': False, 'unlabelled': True})
        'Unlabelled'
        sage: _process_flavor_({'labelled': True, 'unlabelled': False})
        'Labelled'
        sage: _process_flavor_({'labelled': False, 'unlabelled': False})
        Traceback (most recent call last):
        ...
        ValueError: Arguments incompatible.
    """
    labelled = 'Labelled'
    unlabelled = 'Unlabelled'
    flavor = None

    def assign_flavor(flavor, flavor_new):
        if flavor is not None and flavor_new != flavor:
            raise ValueError, "Arguments incompatible."
        return flavor_new

    if kwargs.has_key('unlabelled'):
        flavor = assign_flavor(flavor, unlabelled if kwargs['unlabelled'] else labelled)
    if kwargs.has_key('labelled'):
        flavor = assign_flavor(flavor, labelled if kwargs['labelled'] else unlabelled)

    try:
        del kwargs['labelled']
        del kwargs['unlabelled']
    except KeyError:
        pass

    if flavor is None:
        flavor = unlabelled  # default
    return flavor


#*****************************************************************************
# Data Structures -- Construction
#*****************************************************************************

# TODO maybe change the word "consruction" to something else

def construction(*args, **kwargs):
    return CSConstructionUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CSConstruction(CSBase):
    pass

# ----------------------------------------------------------------------------

class CSConstructionUnlabelled(CSBase, CSUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSConstructionLabelled(CSBase, CSLabelled):
    pass


#*****************************************************************************
# Data Structures -- FiniteSet
#*****************************************************************************

# not sure if needed
class CSFiniteSetBase(CSBase):
    pass

# ----------------------------------------------------------------------------

class CSFiniteSetUnlabelled(CSBase, CSUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSFiniteSetLabelled(CSBase, CSLabelled):
    pass


#*****************************************************************************
# Data Structures -- Singleton
#*****************************************************************************

def singleton(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) singleton of specified size.

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_structures import singleton
        sage: su = singleton('Z', size=2)
        sage: type(su)
        <class 'sage.combinat.combinatorial_structures.CSSingletonUnlabelled'>
        sage: sl = singleton('Z', size=2, labelled=True)
        sage: type(sl)
        <class 'sage.combinat.combinatorial_structures.CSSingletonLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CSSingleton' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

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

# ----------------------------------------------------------------------------

class CSSingletonLabelled(CSSingletonBase, CSFiniteSetLabelled):
    pass

# ----------------------------------------------------------------------------

class CSSingletonUnlabelled(CSSingletonBase, CSFiniteSetUnlabelled):
    pass


#*****************************************************************************
# Data Structures -- Empty
#*****************************************************************************

def empty(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) empty (i.e. singleton of size `0`).

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_structures import empty
        sage: eu = empty('E')
        sage: type(eu)
        <class 'sage.combinat.combinatorial_structures.CSEmptyUnlabelled'>
        sage: el = empty('E', labelled=True)
        sage: type(el)
        <class 'sage.combinat.combinatorial_structures.CSEmptyLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CSEmpty' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

class CSEmptyBase(CSSingletonBase):
    def __init__(self, empty=[]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=empty, size=0)

# ----------------------------------------------------------------------------

class CSEmptyUnlabelled(CSEmptyBase, CSSingletonUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSEmptyLabelled(CSEmptyBase, CSSingletonLabelled):
    pass


#*****************************************************************************
# Data Structures -- Atom
#*****************************************************************************

def atom(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) atom (i.e. singleton of size `1`).

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_structures import atom
        sage: au = atom('A')
        sage: type(au)
        <class 'sage.combinat.combinatorial_structures.CSAtomUnlabelled'>
        sage: al = atom('A', labelled=True)
        sage: type(al)
        <class 'sage.combinat.combinatorial_structures.CSAtomLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CSAtom' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

class CSAtomBase(CSSingletonBase):
    def __init__(self, atom=[[]]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=atom, size=1)

# ----------------------------------------------------------------------------

class CSAtomUnlabelled(CSAtomBase, CSSingletonUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSAtomLabelled(CSAtomBase, CSSingletonLabelled):
    pass


#*****************************************************************************
# Data Structures -- DisjointUnion
#*****************************************************************************

def disjoint_union(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CSDisjointUnionUnlabelled(*args, **kwargs)  # TODO

class CSDisjointUnionBase(CSBase):
    pass

class CSDisjointUnionUnlabelled(CSDisjointUnionBase, CSUnlabelled):
    pass

class CSDisjointUnionLabelled(CSDisjointUnionBase, CSLabelled):
    pass


#*****************************************************************************
# Data Structures -- CartesianProduct
#*****************************************************************************

def cartesian_product(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CSCartesianProductUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CSCartesianProductBase(CSBase):
    pass

# ----------------------------------------------------------------------------

class CSCartesianProductUnlabelled(CSCartesianProductBase, CSUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSCartesianProductLabelled(CSCartesianProductBase, CSLabelled):
    pass


#*****************************************************************************
# Data Structures -- Sequence
#*****************************************************************************

def sequence(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CSSequenceUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CSSequenceBase(CSBase):
    pass

# ----------------------------------------------------------------------------

class CSSequenceUnlabelled(CSSequenceBase, CSUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSSequenceLabelled(CSSequenceBase, CSLabelled):
    pass


#*****************************************************************************
# Data Structures -- MultiSet
#*****************************************************************************

def multi_set(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CSMultiSetUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CSMultiSetBase(CSBase):
    pass

# ----------------------------------------------------------------------------

class CSMultiSetUnlabelled(CSMultiSetBase, CSUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CSMultiSetLabelled(CSMultiSetBase, CSLabelled):
    pass


#*****************************************************************************

