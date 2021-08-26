r"""


Structure Descriptions of groups with order less than or equal to n


EXAMPLES:

sage: from sage.groups.order_n import order_n

sage: from sage.groups.order_n import _directprod

sage: from sage.groups.order_n import _semidirectprod

sage: order_n(12)                                                               
['1',
 'C2',
 'C3',
 'C4',
 'C2 x C2',
 'C5',
 'S3',
 'C6',
 'C7',
 'C8',
 'C4 x C2',
 'D8',
 'Q8',
 'C2 x C2 x C2',
 'C9',
 'C3 x C3',
 'D10',
 'C10',
 'C11',
 'C3 : C4',
 'C12',
 'A4',
 'D12',
 'C6 x C2']

sage: order_n(20, start = 14)                                                   
['D14',
 'C14',
 'C15',
 'C16',
 'C4 x C4',
 '(C4 x C2) : C2',
 'C4 : C4',
 'C8 x C2',
 'C8 : C2',
 'D16',
 'QD16',
 'Q16',
 'C4 x C2 x C2',
 'C2 x D8',
 'C2 x Q8',
 '(C4 x C2) : C2',
 'C2 x C2 x C2 x C2',
 'C17',
 'D18',
 'C18',
 'C3 x S3',
 '(C3 x C3) : C2',
 'C6 x C3',
 'C19',
 'C5 : C4',
 'C20',
 'C5 : C4',
 'D20',
 'C10 x C2']

sage: order_n(10, "IsCyclic", start = 4)                                        
['C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10']

sage: order_n(20,semidirectproduct=True)                                           
['C3 : C4',
 '(C4 x C2) : C2',
 'C4 : C4',
 'C8 : C2',
 '(C4 x C2) : C2',
 '(C3 x C3) : C2',
 'C5 : C4',
 'C5 : C4']

sage: order_n(20,semidirectproduct=True,directproduct=True)                                                                          
['(C4 x C2) : C2', '(C4 x C2) : C2', '(C3 x C3) : C2']

sage: list1 = order_n(20)                                                       
sage: _semidirectprod(list1)                                                    
['C3 : C4',
 '(C4 x C2) : C2',
 'C4 : C4',
 'C8 : C2',
 '(C4 x C2) : C2',
 '(C3 x C3) : C2',
 'C5 : C4',
 'C5 : C4']

sage: _directprod(list1)                                                        
['C2 x C2',
 'C4 x C2',
 'C2 x C2 x C2',
 'C3 x C3',
 'C6 x C2',
 'C4 x C4',
 '(C4 x C2) : C2',
 'C8 x C2',
 'C4 x C2 x C2',
 'C2 x D8',
 'C2 x Q8',
 '(C4 x C2) : C2',
 'C2 x C2 x C2 x C2',
 'C3 x S3',
 '(C3 x C3) : C2',
 'C6 x C3',
 'C10 x C2']


AUTHORS:

- Alexis Newton (2021-07-21): initial version


"""

# ****************************************************************************
#       Copyright (C) 2021 Alexis Newton alexis.newton@emory.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from sage.interfaces.gap import gap
from sage.libs.gap.libgap import libgap
from sage.interfaces.gap import get_gap_memory_pool_size, set_gap_memory_pool_size
from sage.arith.all import factor

def _directprod(order_n):
    """
    This method takes in a python list of string and gives back only the 
    strings with ``x`` in them (for structure descriptions, this 
    corresponds to direct products)
    - ``order_n`` - list of python strings
    AUTHORS:
    - Alexis Newton (2021-07-21)
    """
    new_order_n = []
    for i in range(0,len(order_n)):
        word=order_n[i].split(' ')
        if 'x' in word:
            new_order_n.append(order_n[i])
    return new_order_n


def _semidirectprod(order_n):
    """
    This method takes in a python list of string and gives back only the 
    strings with ``:`` in them (for structure descriptions, this corresponds 
    to semi direct products)
    - ``order_n`` - list of python strings
    AUTHORS:
    - Alexis Newton (2021-07-21)
    """
    new_order_n = []
    for i in range(0,len(order_n)):
        word=order_n[i].split(' ')
        if ':' in word:
            new_order_n.append(order_n[i])
    return new_order_n


def order_n(n, *parameter, start=1, timeout=60, 
            memory=get_gap_memory_pool_size(), 
            directproduct=False, semidirectproduct=False):
    """
    This method produces a list of structure descriptions of all groups less
     than or equal to some order n. It does this by calling the gap functions 
     and then converting them to python strings.
    
    INPUT:
    - ``n`` - the value of the largest order we want a structure description 
        for
    - ``start`` - the lower bound on the order for the produced structure 
        descriptions
    - ``timeout`` - max time in seconds (default: `60`) that the program will 
        spend on each order before moving on (this is to deal with orders that 
        have many small prime fact)
    - ``*parameter`` -- list of properties we can narrow the output by. These 
        include "IsAbelian","IsCyclic","IsSolvable","IsNilpotent","IsSimple",
        "IsDihedralGroup","IsSymmetricGroup","IsAlternatingGroup",
        "IsPerfectGroup","IsPolycyclicGroup"
    - ``memory`` - size in bytes of desired memory pool 
        (default: `get_gap_memory_pool_size()`)
    -``directproduct`` - (default: `0`) set to 1 if want to give only 
        direct product options
    -``semidirectproduct`` - (default: `0`) set to 1 if want to give only 
        semi direct product options

    OUTPUT: the list of structure descriptions for all groups found that 
    meet the criteria 
    
    AUTHORS:
    - Alexis Newton (2021-07-21)
     """
    set_gap_memory_pool_size(memory)
    below_order_n = gap('EmptyPlist(1)') 
    import time
    import warnings
    timeoutafter = time.time() + timeout
    for i in range(start,n+1):
        factors=list(factor(i))
        numoffactors=sum([factors[j][1] for j in range(0,len(factors))])
        if numoffactors > 6:
            warnings.warn('Expensive command when order has many small prime factors. May need more memory.')
        #documentation on AllSmallGroups command and Structure description
        #is in the GAP System reference guide 
        g=gap(i).AllSmallGroups(*parameter)
        for x in range(1,len(g)+1):  
            r=gap(g[x]).StructureDescription()
            gap.Add(below_order_n,r)
            if time.time() > timeoutafter:
                break
    #convert gap object to python string
    below_order_n = libgap(below_order_n).sage()
    #only print out direct or semi direct prod if asked
    if directproduct == True:
        below_order_n = _directprod(below_order_n)
    if semidirectproduct == True:
        below_order_n = _semidirectprod(below_order_n)
   
    return below_order_n








