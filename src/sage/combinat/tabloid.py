# -*- coding: utf-8 -*-
r"""
Tabloid

AUTHORS:

- Jackson Criswell (2016): initial version



This file consists of the following major classes:

Element classes:

* :class:`YoungTableau`
* :class:`Tabloid`
* :class:`StandardTabloid`

Factory classes:

* :class:`YoungTableaux`
* :class:`Tabloids`
* :class:`StandardTabloids`

Parent classes:

* :class:`YoungTableaux_all` (facade class)
* :class:`YoungTableaux_size`
* :class:`YoungTableaux_shape`
* :class:`Tabloids_all` (facade class)
* :class:`Tabloids_size`
* :class:`Tabloids_shape`
* :class:`StandardTabloids_all` (facade class)
* :class:`StandardTabloids_size`
* :class:`StandardTabloids_shape`

For display options, see :meth:`Tableaux.global_options`.

.. TODO:

    
"""

#*****************************************************************************
#       Copyright (C) 2016 Jackson Criswell <crisw1ja@cmich.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.global_options import GlobalOptions
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.rings.infinity import PlusInfinity
from sage.arith.all import factorial, binomial
from sage.rings.integer import Integer
from sage.combinat.composition import Composition, Compositions
from integer_vector import IntegerVectors
import sage.libs.symmetrica.all as symmetrica
import sage.misc.prandom as random
import permutation
import itertools
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.misc.all import uniq, prod
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat import Sets
from sage.combinat.combinatorial_map import combinatorial_map
##JACKSONSTART
from sage.sets.set import Set
from sage.combinat.permutation import Permutation
from sage.combinat.permutation import Permutations 
from sage.combinat.combinatorial_algebra import CombinatorialFreeModule
from sage.rings.rational_field import QQ

##JACKSONEND

from sage.combinat.tableau import *
##############
class Tabloid(Tableau):
    """
    Class to model a tabloid of specified shape, as described in Sagan (2000).

   

    INPUT:

    - ``t`` -- a tableau, a list of iterables, or an empty list

    - ``shape`` -- a partition, or list representation of a partition

    OUTPUT:

    - A Tabloid object constructed from ``t``.

    Tabloids can be consider as 1 column generalized tableau, whose entries
    are disjoint sets of integers in [n], in weakly decreasing order by size

    EXAMPLES::



    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableaux`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`

    TESTS::


    """
    @staticmethod
    def __classcall_private__(self, t):
        r"""
        This ensures that a SemistandardTableau is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

            sage: t = SemistandardTableau([[1,1],[2]])
            sage: TestSuite(t).run()

            sage: t.parent()
            Semistandard tableaux
            sage: t.category()
            Category of elements of Semistandard tableaux
            sage: type(t)
            <class 'sage.combinat.tableau.SemistandardTableaux_all_with_category.element_class'>
        """
        if isinstance(t, Tabloid):
            return t
        elif t in Tabloids():
            return Tabloids_all().element_class(Tabloids_all(), t)
        tabloid=[]
        for row in t:
            tabloid+=[[Set(row)]]
#        return Tableau(tabloid)
        return Tabloids_all().element_class(Tabloids_all(), t)
        # t is not a semistandard tableau so we give an appropriate error message
        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c,(int,Integer)) and c>0 for row in t for c in row):
            raise ValueError("entries must be positive integers"%t)


    def __init__(self, parent, t):
        r"""
        Initialize a tabloid.

        TESTS::

            sage: t = Tableaux()([[1,1],[2]])
            sage: s = SemistandardTableaux(3)([[1,1],[2]])
            sage: s==t
            True
            sage: s.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: r = SemistandardTableaux(3)(t); r.parent()
            Semistandard tableaux of size 3 and maximum entry 3
            sage: isinstance(r, Tableau)
            True
            sage: s2 = SemistandardTableaux(3)([(1,1),(2,)])
            sage: s2 == s
            True
            sage: s2.parent()
            Semistandard tableaux of size 3 and maximum entry 3
        """
        tabloid=[]
        for row in t:
            tabloid+=[[Set(row)]]
        super(Tabloid, self).__init__(parent, tabloid)

        # Tableau() has checked that t is tableau, so it remains to check that
        # the entries of t are positive integers which are weakly increasing
        # along rows
#        from sage.sets.positive_integers import PositiveIntegers
#        PI = PositiveIntegers()

#        for row in t:
#            if any(c not in PI for c in row):
#                raise ValueError("the entries of a tabloid must be non-negative integers")
 


    def size(self):
        """
        Return the size of the tableaux contained in the tabloid.



        """""

        return sum(self.shape())

    def shape(self):
        """
        Return the shape of the tableaux contained in the tabloid.



        """
        from sage.combinat.partition import Partition
        shape=[]
        for row in self:
            shape.append(len(row[0]))
        return Partition(shape)
        
        
    def to_tableau(self):
        tableau=[]
        for row in self:
            tableau.append(row[0].list())
        return Tableau(tableau)
        
    def symmetric_group_action_on_entries(self, pi):
        tableau=self.to_tableau()
        tableau=tableau.symmetric_group_action_on_entries(pi)
        
        return Tabloid(tableau)
        
    def column_stabilizer(self):
        return self.to_tableau().column_stabilizer()
        
        
    def row_stabilizer(self):
        return self.to_tableau().row_stabilizer()
    
        
    def symmetric_group_action_on_values(self, pi):
        tableau=self.to_tableau()
        tableau=tableau.symmetric_group_action_on_values(pi)
        
        return Tabloid(tableau)
    
    def is_standard(self):
        return self.to_tableau().is_standard()
        
    def is_semistandard(self):
        return self.to_tableau().is_semistandard()
    
    def entry(self, cell):
        return self.to_tableau().entry(cell)
        
    def entries(self):
        return self.to_tableau().entries()
        
    def reading_word_permutation(self):
        return self.to_tableau().reading_word_permutation()
        
    def to_word(self):
        return self.to_tableau().to_word()
        
    def pp(self):
        return self.to_tableau().pp()
        
    def conjugate(self):
        return self.to_tableau().conjugate()
  

##########################
# Tabloids #
##########################
class Tabloids(Tableaux):
    """
    A factory class for tabloids.

    INPUT:

    Keyword arguments:

    - ``size`` -- The size of the tableaux
    - ``shape`` -- The shape of the tableaux

    Positional arguments:

    - The first argument is interpreted as either ``size`` or ``shape``
      according to whether it is an integer or a partition
    
    OUTPUT:

    - The appropriate class, after checking basic consistency tests. (For
      example, specifying ``eval`` implies a value for `max_entry`).

    A tabloid is a tableau of shape (1,1,1,...) whose entries are sets of weakly decreasing size.  The sets contain the integers [n].
    The n is the size of the tableau contained in the tabloid, and the shape is the partition formed by the sizes of the sets.
   
    Note that Sage uses the English convention for partitions and tableaux;
    the longer rows are displayed on top.


    .. SEEALSO:

        - :class:`Tableaux`
        - :class:`Tableau`
        - :class:`SemistandardTableau`
        - :class:`StandardTableaux`
        - :class:`StandardTableau`
    """
    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`Tabloids`
        for more information.

        TESTS::

  
        """
        from sage.combinat.partition import Partition, _Partitions
        # Process the keyword arguments -- allow for original syntax where
        #   n == size,  p== shape and mu == eval
        n = kwargs.get('n', None)
        size = kwargs.get('size', n)

        p = kwargs.get('p', None)
        shape = kwargs.get('shape', p)

        # Process the positional arguments
        if args:
            # The first arg could be either a size or a shape
            if isinstance(args[0], (int, Integer)):
                if size is not None:
                    raise ValueError( "size was specified more than once" )
                else:
                    size = args[0]
            else:
                if shape is not None:
                    raise ValueError( "the shape was specified more than once" )
                shape = args[0] # we check it's a partition later


        # Consistency checks
        if size is not None:
            if not isinstance(size, (int, Integer)):
                raise ValueError( "size must be an integer" )
            elif size < 0:
                raise ValueError( "size must be non-negative" )

        if shape is not None:
            from sage.combinat.skew_partition import SkewPartitions
            # use in (and not isinstance) below so that lists can be used as
            # shorthand
            if shape in _Partitions:
                shape = Partition(shape)         
            else:
                raise ValueError( "shape must be a partition" )

        if (size is not None) and (shape is not None):
            if sum(shape) != size:
                # This could return an empty class instead of an error
                raise ValueError( "size and shape are different sizes" )

        # Dispatch appropriately

        if (shape is not None):
            return Tabloids_shape(shape)

        if (size is not None):
            return Tabloids_size(size)

        return Tabloids_all()

    Element = Tabloid

    def __init__(self, **kwds):
        """
        Initialize ``self``.

        EXAMPLES::

        """
        Tableaux.__init__(self, **kwds)

    def __getitem__(self, r):
        r"""
        The default implementation of ``__getitem__`` for enumerated sets
        does not allow slices so we override it.

        EXAMPLES::


        TESTS::

     
        """
        if isinstance(r,(int,Integer)):
            return self.unrank(r)
        elif isinstance(r,slice):
            start=0 if r.start is None else r.start
            stop=r.stop
            if stop is None and not self.is_finite():
                raise ValueError( 'infinite set' )
        else:
            raise ValueError( 'r must be an integer or a slice' )
        count=0
        tabs=[]
        for t in self:
            if count==stop:
                break
            if count>=start:
                tabs.append(t)
            count+=1

        # this is to cope with empty slices endpoints like [:6] or [:}
        if count==stop or stop is None:
            return tabs
        raise IndexError('value out of range')

    def __contains__(self, t):
        """
        Return ``True`` if ``t`` can be interpreted as a
        :class:`Tabloid`.

        TESTS::

        """
        if isinstance(t, Tabloid):
            return True
        elif Tableaux.__contains__(self, t):
            entries=[]
            t=Tableau(t)
            tsize=t.size()
            lastlen=0
            for row in t:
                if len(row)>1:
                    return False
                row=row[0]
                if lastlen!=0 and lastlen<len(row):
                    return False                    
                lastlen=len(row)
                entries.append(row)
            entries=Set(entries)
            numbers=Set(range(1,tsize+1))
            if entries!=numbers:
                return False   
            return True
        else:
            return False
            
            

class Tabloids_all(Tabloids, DisjointUnionEnumeratedSets):
    """
    All tabloids.

    .. WARNING::

        Input is not checked; please use :class:`Tabloids` to
        ensure the options are properly parsed.
    """
    def __init__(self):
        r"""
        Initializes the class of all tabloids.

        TESTS::

 
        """
        tabloids_n = lambda n: Tabloids_size(n)
        DisjointUnionEnumeratedSets.__init__( self,
            Family(NonNegativeIntegers(), tabloids_n),
            facade=True, keepkey = False)

    def _repr_(self):
        """
        TESTS::

        """
#        if self.size is not None:
#            return "Tabloids for Young tableaux of size %s"%str(self.size)
        return "Tabloids"


    def list(self):
        """
        TESTS::

            sage: SemistandardTableaux().list()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        raise NotImplementedError


class Tabloids_size(Tabloids):
    """
    Tabloids of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`SemistandardTableaux`
        to ensure the options are properly parsed.
    """
    def __init__(self, n):
        r"""
        Initializes the class of all tabloids whose tableau have size ``n``.

        TESTS::

            sage: TestSuite( StandardTableaux(0) ).run()
            sage: TestSuite( StandardTableaux(3) ).run()
        """
        super(Tabloids_size, self).__init__(
              category = FiniteEnumeratedSets())
        self.size = Integer(n)

    def _repr_(self):
        """
        TESTS::

       
        """
        return "Tabloids whose tableaux have size %s"%str(self.size)

    def __contains__(self, x):
        """
        EXAMPLES::

        """
        return Tabloids.__contains__(self, x) and sum(map(len, x)) == self.size

    def random_element(self):
        """
        
        """
        return StandardTableaux(self.shape.size()).random_element().to_tabloid_representative()

        
        
        
    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::
        """
        from sage.combinat.partition import Partitions
        c = 0
        for part in Partitions(self.size):
            c += Tabloids_shape(part).cardinality()
        return c


    def __iter__(self):
        """
        EXAMPLES::

        """
        from sage.combinat.partition import Partitions
        for part in Partitions(self.size):
            for tabloid in Tabloids_shape(part):
                yield self.element_class(self, tabloid)

class Tabloids_shape(Tabloids):
    """
    Tabloids whose tableaux have fixed shape `p`, and contain the integers [n]
    where `n` is the size of `p`.

    Tabloid contains all positive integers i, i at most n.
    
    INPUT:

    - ``p`` -- A partition

    .. WARNING::

        Input is not checked; please use :class:`Tabloids` to
        ensure the options are properly parsed.
    """
    def __init__(self, p):
        r"""
        Initializes the class of tabloids of shape ``p``.

        TESTS::

        """
        from sage.combinat.partition import Partition
        super(Tabloids_shape, self).__init__(category = FiniteEnumeratedSets())
        self.shape = Partition(p)


    def __iter__(self):
        """
        An iterator for the tabloids of the specified shape.

        EXAMPLES::
        """
        raise NotImplementedError
        
    def list(self):
        raise NotImplementedError
        
        
    def __contains__(self, x):
        """
        EXAMPLES::

        """
        if not(x.shape()==self.shape):
            return false
        else:
            return Tabloids.__contains__(self, x)# and Tabloid(x).shape() == self.shape()

    def _repr_(self):
        """
        TESTS::
        """
        return "Tabloids of shape %s." %str(self.shape)

    def random_element(self):
        """
        
        """
        return StandardTableaux(self.shape).random_element().to_tabloid_representative()

    def cardinality(self):
        raise NotImplementedError


