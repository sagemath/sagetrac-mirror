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

    take over the world
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
from sage.combinat.permutation_module import PermutationModule, SpechtModule
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
        This ensures that a Tabloid is only ever constructed as an
        element_class call of an appropriate parent.

        TESTS::

        """
        if isinstance(t, Tabloid):
            return t
        elif t in Tabloids():
            return Tabloids_all().element_class(Tabloids_all(), t)

        if t not in Tableaux():
            raise ValueError('%s is not a tableau' % t)

        if not all(isinstance(c,(int,Integer)) and c>0 for row in t for c in row):
            raise ValueError("entries must be positive integers"%t)

        tabloid=[]
        for row in t:
            tabloid+=[Set(row)]
#        return Tableau(tabloid)
        return Tabloids_all().element_class(Tabloids_all(), tabloid)
 

    def __init__(self, parent, t):
        r"""
        Initialize a tabloid.

        TESTS::

        """
        
        if isinstance(t, Tabloid):
            super(Tabloid, self).__init__(parent, t)
        else:
            tabloid=[]
            for row in t:
                tabloid+=[[Set(row)]]
            super(Tabloid, self).__init__(parent, tabloid)

 


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
        
    def to_tableaux(self):
        tableaux=Set([])
        tableau=self.to_tableau()
        RS=self.row_stabilizer()
        for pi in RS:
            tableaux+=Set([tableau.permute(Permutation(pi))])
            
        return tableaux.list()
            
        
        
        
    def symmetric_group_action_on_entries(self, pi):
        tableau=self.to_tableau()
        tableau=tableau.symmetric_group_action_on_entries(pi)
        
        return Tabloid(tableau)
        
    def column_stabilizer(self):
        tableaux=self.to_tableaux()

        CS=Set([])
        for T in tableaux:
            CS+=Set([T.column_stabilizer()])

        return CS
        
    def row_stabilizer(self):
        return self.to_tableau().row_stabilizer()
    
        
    def symmetric_group_action_on_values(self, pi):
        tableau=self.to_tableau()
        tableau=tableau.symmetric_group_action_on_values(pi)
        
        return Tabloid(tableau)
    
    def is_standard(self):
        tableaux=self.to_tableaux()
        for T in tableaux:
            if T.is_standard()==True:
                return True
        return False
        
    def is_semistandard(self):
        tableaux=self.to_tableaux()
        for T in tableaux:
            if T.is_semistandard()==True:
                return True
        return False
    
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
                if type(row) == int:
                    length=1
                else:
                    length=len(row)
#                print length
                if length>1:
                    return False
                row=row[0]
                if lastlen!=0 and lastlen<length:
                    return False                    
                
                lastlen=length
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

        """
        raise NotImplementedError


class Tabloids_size(Tabloids):
    """
    Tabloids of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`Tabloids`
        to ensure the options are properly parsed.
    """
    def __init__(self, n):
        r"""
        Initializes the class of all tabloids whose tableau have size ``n``.

        TESTS::

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

        
    def list(self):
        from sage.combinat.partition import Partitions
        tabloids = []
        for part in Partitions(self.size):
            tabloids += Tabloids_shape(part).list()
        return tabloids

        
        
    def cardinality(self):
        """
        Return the cardinality of ``self``.

        EXAMPLES::
        """
        return len(self.list())


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
        self.size=self.shape.size()


    def __iter__(self):
        """
        An iterator for the tabloids of the specified shape.

        EXAMPLES::
        """
        raise NotImplementedError
        
    def list(self):
        SYT=StandardTableaux(self.shape)
        tabloids=Set([])
        for T in SYT:
            for pi in Permutations(T.size()):
                tabloids+= Set([T.permute(pi).to_tabloid()])                
        return tabloids.list()


#        raise NotImplementedError
        
        
    def __contains__(self, x):
        """
        EXAMPLES::

        """
        if not(Tabloid(x).shape()==self.shape):
            return False
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
        return len(self.list())
#        raise NotImplementedError

##############
class StandardTabloid(Tabloid):
    """
    Class to model a Standard tabloid of specified shape, as described in Sagan (2000).

   

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
        This ensures that a :class:`StandardTabloid` is only ever constructed
        as an ``element_class`` call of an appropriate parent.
        """
        if isinstance(t, StandardTabloid):
            return t

        return StandardTabloids_all().element_class(StandardTabloids_all(), t)

    def __init__(self, parent, t):
        r"""
        Initializes a standard tabloid.

        """
        super(StandardTabloid, self).__init__(parent, t)

        # t is semistandard so we only need to check
        # that its entries are in bijection with {1, 2, ..., n}
        flattened_list = [i for row in self.to_tableau() for i in row]
        if sorted(flattened_list) != range(1, len(flattened_list)+1):
            raise ValueError("the entries in a standard tableau must be in bijection with 1,2,...,n")

        if self.is_standard()==False:
            raise ValueError("Tabloid is not standard")


##########################
# Tabloids #
##########################
class StandardTabloids(Tabloids):
    """
    A factory class for Standard tabloids.

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
        arguments.  See the documentation for :class:`StandardTabloids` for
        more information.

        TESTS::

        """
        from sage.combinat.partition import _Partitions, Partition

        if args:
            n = args[0]
        elif 'n' in kwargs:
            n = kwargs[n]
        else:
            n = None

        if n is None:
            return StandardTabloids_all()

        elif n in _Partitions:
            return StandardTabloids_shape(Partition(n))

        if not isinstance(n, (int, Integer)) or n < 0:
            raise ValueError( "The argument must be a non-negative integer or a partition." )

        return StandardTabloids_size(n)

    Element = StandardTabloid

    def __contains__(self, x):
        """
  
        """
        if isinstance(x, StandardTabloid):
            return x.is_standard()
        elif Tabloids.__contains__(self, x):
            return x.is_standard()
        elif Tableaux.__contains__(self, x):
            return x.is_standard()

class StandardTabloids_all(StandardTabloids, DisjointUnionEnumeratedSets):
    """
    All standard tabloids.

    .. WARNING::

        Input is not checked; please use :class:`Tabloids` to
        ensure the options are properly parsed.
    """
    """
    All standard tabloids.
    """
    def __init__(self):
        r"""
        Initializes the class of all standard tabloids.

            sage: TestSuite(ST).run()
        """
        DisjointUnionEnumeratedSets.__init__( self,
                Family(NonNegativeIntegers(), StandardTabloids_size),
                facade=True, keepkey = False)

    def _repr_(self):
        """

        """
        return "Standard tabloids"



    def list(self):
        """
        TESTS::

        """
        raise NotImplementedError


class StandardTabloids_size(Tabloids):
    """
    Standard Tabloids of fixed size `n`.

    .. WARNING::

        Input is not checked; please use :class:`Tabloids`
        to ensure the options are properly parsed.
    """
    def __init__(self, n):
        r"""
        Initializes the class of all standard tabloids of size ``n``.


        """
        super(StandardTabloids_size, self).__init__(
              category = FiniteEnumeratedSets())
        self.size = Integer(n)


    def _repr_(self):
        """
        
        """
        return "Standard tabloids of size %s"%self.size

    def __contains__(self, x):
        """
 
        """
        return StandardTabloids.__contains__(self, x) and x.size()==self.size

    def __iter__(self):
        """
        """
        from sage.combinat.partition import Partitions
        for p in Partitions(self.size):
            for st in StandardTabloids(p):
                yield self.element_class(self, st)

    def cardinality(self):
        r"""
        
        """
        return StandardTableaux(self.size).cardinality()
        
        
    def random_element(self):
    
        r"""
        
        
        """
        return StandardTableaux(self.size).random_element()


class StandardTabloids_shape(Tabloids):
    """
    Standard Tabloids whose tableaux have fixed shape `p`, and contain the integers [n]
    where `n` is the size of `p`.

    Tabloid contains all positive integers i, i at most n.
    
    Standard Tabloids contain at least one standard tableau.
    INPUT:

    - ``p`` -- A partition

    .. WARNING::

        Input is not checked; please use :class:`Tabloids` to
        ensure the options are properly parsed.
    """
    def __init__(self, p):
        r"""
        Initializes the class of all standard tabloids of a given shape.

        """
        from sage.combinat.partition import Partition
        super(StandardTabloids_shape, self).__init__(category = FiniteEnumeratedSets())
        self.shape = Partition(p)


    def __contains__(self, x):
        """
        EXAMPLES::

        """
        return StandardTabloids.__contains__(self, x) and [len(_) for _ in x] == self.shape

    def _repr_(self):
        """
      
        """
        return "Standard tabloids of shape %s"%str(self.shape)

    def cardinality(self):
        return StandardTableaux(self.shape).cardinality()

    def __iter__(self):
        r"""
        An iterator for the standard Young tabloids associated to the
        shape `p` of ``self``.

        EXAMPLES::

            sage: [t for t in StandardTableaux([2,2])]
            [[[1, 3], [2, 4]], [[1, 2], [3, 4]]]
            sage: [t for t in StandardTableaux([3,2])]
            [[[1, 3, 5], [2, 4]],
             [[1, 2, 5], [3, 4]],
             [[1, 3, 4], [2, 5]],
             [[1, 2, 4], [3, 5]],
             [[1, 2, 3], [4, 5]]]
            sage: st = StandardTableaux([2,1])
            sage: st[0].parent() is st
            True
        """

        pi = self.shape
        #Set the initial tableau by filling it in going down the columns
        tableau = [[None]*n for n in pi]
        size = sum(pi)
        row = 0
        col = 0
        for i in range(size):
            tableau[row][col] = i+1

            #If we can move down, then do it;
            #otherwise, move to the next column over
            if ( row + 1 < len(pi) and col < pi[row+1]):
                row += 1
            else:
                row = 0
                col += 1

        yield self.element_class(self, Tableau(tableau).to_tabloid())

        # iterate until we reach the last tableau which is
        # filled with the row indices.
        last_tableau = sum([[row]*l for (row,l) in enumerate(pi)], [])

        #Convert the tableau to "vector format"
        #tableau_vector[i] is the row that number i
        #is in
        tableau_vector = [None]*size
        for row in range(len(pi)):
            for col in range(pi[row]):
                tableau_vector[tableau[row][col]-1] = row

        while tableau_vector!=last_tableau:
            #Locate the smallest integer j such that j is not
            #in the lowest corner of the subtableau T_j formed by
            #1,...,j.  This happens to be first j such that
            #tableau_vector[j]<tableau_vector[j-1].
            #l will correspond to the shape of T_j
            l = [0]*size
            l[0] = 1
            j = 0
            for i in range(1,size):
                l[tableau_vector[i]] += 1
                if ( tableau_vector[i] < tableau_vector[i-1] ):
                    j = i
                    break

            #Find the last nonzero row of l and store it in k
            i = size - 1
            while ( l[i] == 0 ):
                i -= 1
            k = i

            #Find a new row for the letter j (next lowest corner)
            t = l[ 1 + tableau_vector[j] ]
            i = k
            while ( l[i] != t ):
                i -= 1

            #Move the letter j to row i
            tableau_vector[j] = i
            l[i] -= 1

            #Fill in the columns of T_j using 1,...,j-1 in increasing order
            m = 0
            while ( m < j ):
                r = 0
                while ( l[r] != 0 ):
                    tableau_vector[m] = r
                    l[r] -= 1
                    m += 1
                    r += 1

            #Convert the tableau vector back to the regular tableau
            #format
            row_count= [0]*len(pi)
            tableau = [[None]*n for n in pi]

            for i in range(size):
                tableau[tableau_vector[i]][row_count[tableau_vector[i]]] = i+1
                row_count[tableau_vector[i]] += 1

            yield self.element_class(self, Tableau(tableau).to_tabloid())

        return


    def list(self):
        r"""
        Return a list of the standard Young tabloids of the specified shape.


        """
        S=StandardTabloids(self.shape)
        SYT=StandardTableaux(self.shape).list()
        
        return [T.to_tabloid() for T in SYT]


    def random_element(self):
        """
    
        """
        S=StandardTabloids(self.shape)

        return StandardTableaux(self.shape).random_element().to_tabloid()
        
        
        
        
"""

MODULES

- moved into this file to resolve library loading issues
- will figure out more permanent fix
xxxx
-moved back out did not seem to help
"""


        
        
        
