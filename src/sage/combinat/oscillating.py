#Created on Sat May 19 07:15:40 2018
"""
This is an implementation of the category PathTableaux.

The purpose is to define promotion, evacuation and the action of
the cactus group on oscillating tableaux.

AUTHORS:
- Bruce Westbury (2018): initial version
"""

#*****************************************************************************
#       Copyright (C) 2018 Bruce Westbury <bruce.westbury@gmail.com>,
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from six import add_metaclass

from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableList
from sage.structure.parent import Parent
from sage.combinat.perfect_matching import PerfectMatching
from sage.categories.pathtableaux import PathTableaux
from sage.combinat.tableau import Tableau
from sage.combinat.partition import Partition
from sage.modules.free_module_element import vector

@add_metaclass(InheritComparisonClasscallMetaclass)
class OscillatingTableau(ClonableList):

    @staticmethod
    def __classcall_private__(self, ot):

        if isinstance(ot, OscillatingTableau):
            return ot
        
        if isinstance(ot,PerfectMatching):
            return OscillatingTableaux().from_perfectmatching(ot)
        
        if isinstance(ot,(tuple,list)):
            try:
                ot = tuple([ Partition(a) for a in ot ])
                return OscillatingTableaux()(ot)
            except TypeError:
                pass
            
            ot =tuple([Integer(a) for a in ot])
            return OscillatingTableaux().from_word(ot)
            
        raise ValueError("Sorry, I don't know what to do with %s." % str(ot) )

    def _hash_(self):
        return hash(tuple(map(tuple, self)))

    def check(self):
        n = len(self)
        for i in range(n-1):
            h = self[i]
            t = self[i+1]

            if abs(t.size()-h.size()) != 1:
                    raise ValueError("Adjacent partitions differ by more than one box.")
            if t.size() == h.size()+1:
                if not t.contains(h):
                    raise ValueError("Next partition is not obtained by adding a cell.")
            if h.size() == t.size()+1:
                if not h.contains(t):
                    raise ValueError("Next partition is not obtained by removing a cell.")

    @staticmethod
    def _rule(x):
        y = map(list,x)
        m = max([ len(u) for u in y ])
        z = map( lambda u: vector(u + [0]*(m-len(u)) ), y )
        result = z[0]-z[1]+z[2]
        result = map(abs,result)
        result.sort(reverse=True)
        return Partition(result)

    def is_skew(self):
        """
        Returns True if Tableau is skew and False if not.

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.is_skew()
        False
        """
        return self[0] != Partition([])

    def crossing_number(self):
        """
        Returns the crossing number.
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.crossing_number()
        1
        """
        return max( a.length() for a in self )

    def nesting_number(self):
        """
        Returns the nesting number.

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.nesting_number()
        2
        """

        if all( a.length() == 0 for a in self ):
            return 0
        return max( a[0] for a in self if a.length() > 0)

    def to_word(self):
        """
        Converts an oscillating tableau to a word in the alphabet
        ...,-2,-1,1,2,...

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.to_word()
        [1, 1, -1, -1]

        """
        n = len(self)
        m = self.crossing_number()

        result = [0]*(n-1)
        for i in range(n-1):
            u = list(self[i])
            u = u + [0]*(m-len(u))
            v = list(self[i+1])
            v = v + [0]*(m-len(v))
            for j in range(m):
                if u[j] != v[j]:
                    result[i] = v[j]-u[j]
        return result

    def descents(self):
        """
        Returns the descent set. This is defined in
        https://arxiv.org/abs/1303.5850

        EXAMPLE:
        sage: T = OscillatingTableau([[],[1],[2],[1],[]])
        sage: T.descents()
        {2}

        """
        result = set()
        w = self.to_word()

        for i in range(1,len(w)):
            if w[i-1] > 0 and w[i] < 0:
                result.add(i)
            if 0 < w[i-1] < w[i]:
                result.add(i)
            if w[i-1] < w[i] < 0:
                result.add(i)

        return result
    
    def sundaram(self):
        """
        This implements the bijection due to S. Sundaram between
        oscillating tableaux with empty initial shape and pairs
        (S,M) where S is a partial standard tableaux whose shape
        is the final shape of the oscillating tableau and a 
        perfect matching on the complement of the set of entries
        of S.
        
        INPUT: A straight oscillating tableau.
        
        OUTPUT: A pair (S,M); S is a Tableau and M is a PerfectMatching
        
        sage: t = OscillatingTableau([[],[1],[1,1],[1],[]])
        sage: t.sundaram()
        ([], [(1, 3), (2, 4)])
        
        sage: t = OscillatingTableau([[],[1],[1,1],[2,1],[2],[1],[2],[2,1],[2,1,1],[2,1]])
        sage: t.sundaram()
        ([[2, 7], [8]], [(1, 4), (3, 5), (6, 9)])

        sage: s = OscillatingTableau([[],[1],[2],[2,1],[1,1],[1],[1,1],[2,1],[3,1],[2,1]])
        sage: s.sundaram()
        ([[3, 7], [6]], [(1, 5), (2, 4), (8, 9)])

        """
        if self.is_skew():
            raise ValueError("This has only been implemented for straight oscillating tableaux.")
        tb = Tableau([])
        pm = set([])
        
        for i in range(1,len(self)):
            lb = self[i]
            mu = self[i-1]
            if lb.contains(mu):
                cell = [ c for c in lb.corners() if c not in mu.corners() ][0]
                tb = tb.add_entry(cell,i)
            else:
                cell = [ c for c in mu.corners() if c not in lb.corners() ][0]
                tb, x = tb.reverse_bump(cell)
                pm.add((i,x))
                
            
        return (tb,PerfectMatching(pm))

###############################################################################

class OscillatingTableaux(UniqueRepresentation,Parent):

    @staticmethod
    def __classcall_private__(cls):
        return super(OscillatingTableaux, cls).__classcall__(cls)

    def __init__(self):

        Parent.__init__(self, category=PathTableaux())

    def __contains__(self, ot):

        return isinstance(ot, (list, tuple, OscillatingTableau))

    def _element_constructor_(self, ot, check=True):

        if isinstance(ot, OscillatingTableaux) and ot.parent() == self:
            return ot

        return self.element_class(self, list(ot))

    def from_word(self,w):
        """
        This constructs an oscillating tableaux from a list of
        integers. The list may contain positive or negative entries
        but not zero.
        
        This function attempts to construct a straight oscillating
        tableaux. This may fail.
        
        sage: OscillatingTableau([1,-1])
        [[], [1], []]

        """
        if any([a==0 for a in w]):
            raise ValueError("List may not contain zero.")
            
        ot = [Partition([])]*(len(w)+1)
        
        for i,a in enumerate(w):
            if a > 0:
                ot[i+1] = ot[i].add_cell(a-1)
            else:
                pt = list(ot[i])  
                pt[abs(a)-1] -= 1
                ot[i+1] = Partition(pt)
            
        return self.element_class(self,ot)
    
    def from_perfectmatching(self,pm):
        """
        Constructs an oscillating tableau from a perfect matching.
        The intial and final shapes are the empty partition.
        
        sage: pm = PerfectMatching([[1,2]])
        sage: OscillatingTableau(pm)
        [[], [1], []]

        sage: pm = PerfectMatching([[1,3],[2,4]])
        sage: OscillatingTableau(pm)
        [[], [1], [1, 1], [1], []]

        """
        tb = Tableau([])
        n = pm.size()
        ot = [Partition([])]*(n+1)
        
        for i in range(n,1,-1):
            c = tb.cells_containing(i)
            if len(c) > 1:
                raise RuntimeError("Tableau should be standard.")
            if len(c) == 1:
                tc = list(tb)
                tc[c[0][0]] = list(tc[c[0][0]]).remove(i)
                tc = [ a for a in tc if a != None ]
                tb = Tableau(tc)
                ot[i-1] = tb.shape()
            else:
                x = pm.partner(i)
                tb = tb.bump(x)
                ot[i-1] = tb.shape()
                
        return self.element_class(self,ot)
    
    Element = OscillatingTableau

