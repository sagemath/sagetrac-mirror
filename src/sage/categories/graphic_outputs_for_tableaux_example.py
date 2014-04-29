r"""
"""
#*****************************************************************************
#  Copyright (C) 2014 Adrien Boussicault (boussica@labri.fr), 
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.element_wrapper import ElementWrapper
from sage.structure.list_clone import ClonableList
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.set_factories import (
    SetFactory, ParentWithSetFactory, TopMostParentPolicy
)
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.global_options import GlobalOptions
from sage.sets.set import Set
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_attribute import lazy_attribute
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.rings.integer import Integer
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers 
from sage.misc.cachefunc import cached_method
from sage.misc.latex import latex
from copy import deepcopy
from sage.matrix.constructor import matrix
from sage.combinat.combinat import catalan_number
from sage.combinat.combinatorial_map import combinatorial_map
from sage.categories.graphic_outputs_for_tableaux import GraphicOutputsFor2DTableaux

class TableauExample( ClonableList ):
    r"""
    """
    __metaclass__ = ClasscallMetaclass

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        """
        return cls._auto_parent._element_constructor_( *args, **opts )

    @lazy_class_attribute
    def _auto_parent(cls):
        r"""
        """
        return TableauExamples()

    def __copy__( self ):
        r"""
        """
        return TableauExample( [ self.lower_path(), self.upper_path() ] )

    def check( self ):
        r"""
        """
        pass

    def array(self):
        return self

    def __init__(self, parent, value, check=True):
        r"""
        """
        ClonableList.__init__(self, parent, value )
        if check:
            if not isinstance( value, (list, tuple) ) :
                raise ValueError, "Value %s must be a list or a tuple."%(value)
            self.check()
        self._options = None

    def size( self ):
        return len(self)

    def _repr_( self ):
        return self._repr_tableau()

class TableauExamplesFactory(SetFactory):
    r"""
    """
    def __call__(self, size=None, policy=None):
        r"""
        """
        if policy is None:
            policy = self._default_policy

        if isinstance(size, (Integer, int)):
            return TableauExamples_size(size, policy)
        if size is None:
            return TableauExamples_all(policy)
        raise ValueError, "Invalide argument for Tre-like tableaux Factory."

    def add_constraints(self, cons, (args, opts)):
        r"""
        """
        return cons+args

    @lazy_attribute
    def _default_policy(self):
        return TopMostParentPolicy(self, (), TableauExample)

    def _repr_(self):
        """
        """
        return "Factory for tableau examples"

TableauExamples = TableauExamplesFactory()
TableauExamples.__doc__ = TableauExamplesFactory.__call__.__doc__

class TableauExamples_size(ParentWithSetFactory, UniqueRepresentation):
    r"""
    The tableau examples of size `n`.
    """
    def __init__(self, size, policy):
        r"""
        Construct a set of Parallelogram Polyominoes of a given size.
        """
        self._size = size
        ParentWithSetFactory.__init__(
            self, (size,), policy, category = [ 
                FiniteEnumeratedSets(), GraphicOutputsFor2DTableaux()
            ]
        )

    def _repr_(self):
        r"""
        Return the string representation of the set of tableau examples

        EXAMPLES::
        
            sage: TableauExamples( 3 )
            tableau examples of size 3
        """
        return "tableau examples of size %s"%(self._size)

    def an_element(self):
        r"""
        """
        return next( self.__iter__() )

    def check_element(self, el, check):
        r"""
        """
        if el.size() != self.size():
            raise ValueError, "The parallelogram polyomino have a Wrong size : %s"%(el.size())

    def cardinality( self ):
        r"""
        """
        return catalan_number( self.size() )

    def __iter__(self):
        r"""
        """
        for i in range( self.size() ):
             yield TableauExample( [ i for j in range(self.size()) ] )

    def size( self ):
        r"""
        """
        return self._size

class TableauExamples_all( ParentWithSetFactory, DisjointUnionEnumeratedSets ):
    r"""
    This class enumerate all the tableau examples.
    """
    def __init__(self, policy):
        r"""
        """
        ParentWithSetFactory.__init__(
            self, (), policy, category = [ 
                FiniteEnumeratedSets(), GraphicOutputsFor2DTableaux()
            ]
        )
        DisjointUnionEnumeratedSets.__init__(
            self, Family(
                NonNegativeIntegers(), self._tableau_exemples_size
            ),
            facade=True, keepkey = False,
            category = self.category()
        )

    def _tableau_exemples_size( self, n ):
        return TableauExamples_size( n, policy=self.facade_policy() )

    def _repr_(self):
        r"""
        """
        return "tableau examples"

    def check_element(self, el, check):
        r"""
        """
        pass
