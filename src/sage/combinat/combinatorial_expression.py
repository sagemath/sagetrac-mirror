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


import sage
from sage.rings.integer import Integer
from sage.misc.lazy_attribute import lazy_attribute
from sage.structure.sage_object import SageObject  # TODO löschen
from itertools import izip

#*****************************************************************************
# Flavor
#*****************************************************************************

class _GenericFlavor_(sage.structure.sage_object.SageObject):

q    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` if labeled, ``False`` if unlabeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: F = _GenericFlavor_()
            sage: F.is_labeled() is None
            True
        """
        return None


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` if unlabeled, ``False`` if labeled, ``None`` otherwise.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _GenericFlavor_)
            sage: F = _GenericFlavor_()
            sage: F.is_unlabeled() is None
            True
        """
        return None


    @staticmethod
    def _class_with_prefix_(prefix, classname):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return globals()[prefix + classname]


    @classmethod
    def class_generic(cls, classname):
        """
        TODO

        INPUT:

        - ```` --

        OUTPUT:

        EXAMPLES::

            sage: TODO  # not tested
        """
        return cls._class_with_prefix_('Generic', classname)


    @classmethod
    def class_unlabeled(cls, classname):
        return cls._class_with_prefix_('Unlabeled', classname)


    @classmethod
    def class_labeled(cls, classname):
        return cls._class_with_prefix_('Labeled', classname)


# ----------------------------------------------------------------------------


class _UnlabeledFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``False`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_labeled()
            False
        """
        return False


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _UnlabeledFlavor_)
            sage: F = _UnlabeledFlavor_()
            sage: F.is_unlabeled()
            True
        """
        return True


# ----------------------------------------------------------------------------


class _LabeledFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_labeled()
            True
        """
        return True


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``False`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _LabeledFlavor_)
            sage: F = _LabeledFlavor_()
            sage: F.is_unlabeled()
            False
        """
        return False


# ----------------------------------------------------------------------------


class _EmptyFlavor_(_GenericFlavor_):

    def is_labeled(self):
        """
        Returns whether combinatorial expression is labeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is labeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _EmptyFlavor_)
            sage: F = _EmptyFlavor_()
            sage: F.is_labeled()
            True
        """
        return True


    def is_unlabeled(self):
        """
        Returns whether combinatorial expression is unlabeled or not.

        INPUT:

        Nothing

        OUTPUT:

        ``True`` since this instance is unlabeled.

        TESTS::

            sage: from sage.combinat.combinatorial_expression import (
            ....:     _EmptyFlavor_)
            sage: F = _EmptyFlavor_()
            sage: F.is_unlabeled()
            True
        """
        return True


#*****************************************************************************
# Data Structures -- Singleton
#*****************************************************************************

def singleton(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) singleton of specified size.

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_expression import singleton
        sage: su = singleton('Z', size=2)
        sage: type(su)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionSingletonUnlabelled'>
        sage: sl = singleton('Z', size=2, labelled=True)
        sage: type(sl)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionSingletonLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CombinatorialExpressionSingleton' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

class CombinatorialExpressionSingletonBase(CombinatorialExpressionFiniteSetBase):
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

class CombinatorialExpressionSingletonLabelled(CombinatorialExpressionSingletonBase, CombinatorialExpressionFiniteSetLabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionSingletonUnlabelled(CombinatorialExpressionSingletonBase, CombinatorialExpressionFiniteSetUnlabelled):
    pass


#*****************************************************************************
# Data Structures -- Empty
#*****************************************************************************

def empty(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) empty (i.e. singleton of size `0`).

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_expression import empty
        sage: eu = empty('E')
        sage: type(eu)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionEmptyUnlabelled'>
        sage: el = empty('E', labelled=True)
        sage: type(el)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionEmptyLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CombinatorialExpressionEmpty' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

class CombinatorialExpressionEmptyBase(CombinatorialExpressionSingletonBase):
    def __init__(self, empty=[]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=empty, size=0)

# ----------------------------------------------------------------------------

class CombinatorialExpressionEmptyUnlabelled(CombinatorialExpressionEmptyBase, CombinatorialExpressionSingletonUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionEmptyLabelled(CombinatorialExpressionEmptyBase, CombinatorialExpressionSingletonLabelled):
    pass


#*****************************************************************************
# Data Structures -- Atom
#*****************************************************************************

def atom(*args, **kwargs):
    """
    Returns the encapsulated (in TODO) atom (i.e. singleton of size `1`).

    TODO

    TESTS::

        sage: from sage.combinat.combinatorial_expression import atom
        sage: au = atom('A')
        sage: type(au)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionAtomUnlabelled'>
        sage: al = atom('A', labelled=True)
        sage: type(al)
        <class 'sage.combinat.combinatorial_expression.CombinatorialExpressionAtomLabelled'>
    """
    flavor = _process_flavor_(kwargs)
    return globals()['CombinatorialExpressionAtom' + flavor](*args, **kwargs)

# ----------------------------------------------------------------------------

class CombinatorialExpressionAtomBase(CombinatorialExpressionSingletonBase):
    def __init__(self, atom=[[]]):
        """
        TODO

        EXAMPLES::

            sage: TODO  # not tested
        """
        self._set_singleton_(singleton=atom, size=1)

# ----------------------------------------------------------------------------

class CombinatorialExpressionAtomUnlabelled(CombinatorialExpressionAtomBase, CombinatorialExpressionSingletonUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionAtomLabelled(CombinatorialExpressionAtomBase, CombinatorialExpressionSingletonLabelled):
    pass


#*****************************************************************************
# Data Structures -- DisjointUnion
#*****************************************************************************

def disjoint_union(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CombinatorialExpressionDisjointUnionUnlabelled(*args, **kwargs)  # TODO

class CombinatorialExpressionDisjointUnionBase(CombinatorialExpressionBase):
    pass

class CombinatorialExpressionDisjointUnionUnlabelled(CombinatorialExpressionDisjointUnionBase, CombinatorialExpressionUnlabelled):
    pass

class CombinatorialExpressionDisjointUnionLabelled(CombinatorialExpressionDisjointUnionBase, CombinatorialExpressionLabelled):
    pass


#*****************************************************************************
# Data Structures -- CartesianProduct
#*****************************************************************************

def cartesian_product(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CombinatorialExpressionCartesianProductUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CombinatorialExpressionCartesianProductBase(CombinatorialExpressionBase):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionCartesianProductUnlabelled(CombinatorialExpressionCartesianProductBase, CombinatorialExpressionUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionCartesianProductLabelled(CombinatorialExpressionCartesianProductBase, CombinatorialExpressionLabelled):
    pass


#*****************************************************************************
# Data Structures -- Sequence
#*****************************************************************************

def sequence(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CombinatorialExpressionSequenceUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CombinatorialExpressionSequenceBase(CombinatorialExpressionBase):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionSequenceUnlabelled(CombinatorialExpressionSequenceBase, CombinatorialExpressionUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionSequenceLabelled(CombinatorialExpressionSequenceBase, CombinatorialExpressionLabelled):
    pass


#*****************************************************************************
# Data Structures -- MultiSet
#*****************************************************************************

def multi_set(*args, **kwargs):
    # TODO: make here the decision if labelled or unlabelled is used
    # (depending on args)
    return CombinatorialExpressionMultiSetUnlabelled(*args, **kwargs)  # TODO

# ----------------------------------------------------------------------------

class CombinatorialExpressionMultiSetBase(CombinatorialExpressionBase):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionMultiSetUnlabelled(CombinatorialExpressionMultiSetBase, CombinatorialExpressionUnlabelled):
    pass

# ----------------------------------------------------------------------------

class CombinatorialExpressionMultiSetLabelled(CombinatorialExpressionMultiSetBase, CombinatorialExpressionLabelled):
    pass


#*****************************************************************************

