r"""
    Functions for computing the profiles of sub-Hopf algebras of the
    mod p Steenrod algebra.

    .. NOTE::
        Palmieri pads all his profiles by adding a [0] to the beginning.
        Furthermore, when scanning over elements defined over an odd
        primary Steenrod Algebra, the Q part comes first, and then the P part.
        When defining a sub-algebra over an odd Steenrod algebra, the P part
        must come first in the profile.

        These functions return the P part first, and then the Q part.

AUTHORS:

    - Robert R. Bruner, Michael J. Catanzaro (2012): initial version

"""

#*****************************************************************************
#       Copyright (C) 2011 Robert R. Bruner <rrb@math.wayne.edu> and
#                          Michael J. Catanzaro <mike@math.wayne.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************



#----------------------------Profile-Functions-----------------------------------

from sage.algebras.steenrod.steenrod_algebra import SteenrodAlgebra
from copy import copy


#import sage.modules.fp_over_steenrod_algebra.utility as Utility
def mod_p_log(n,p):
    r"""
    Input an integer $n$ and a prime $p$
    Output the $k$ so that $p^{k-1} < n <= p^k$

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: mod_p_log(1,4)
        1
        sage: mod_p_log(8,3)
        2
        sage: mod_p_log(9,3)
        3

    """
    k=0
    pow=1
    while n >= pow:
        k += 1
        pow *= p
    return k


def profile_ele(alist,char=2):
    """
    Finds the smallest sub-Hopf algebra containing the element passed.
    `alist' is assumed to be an element of the Steenrod Algebra, (the only other
    possible cases are dealt with immediately) so it is treated as a list.

    INPUT:

    -  ``alist`` -  An element of the Steenrod Algebra (or a sub-Hopf algebra
                 of it). Treated as a list.

    OUTPUT: The profile function corresponding to the smallest sub-Hopf algebra
    containing the element passed.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: A2 = SteenrodAlgebra(2)
        sage: profile_ele(A2.Sq(2))
        (2, 1)
        sage: profile_ele(A2.Sq(4,8))
        (3, 4, 3, 2, 1)

    """
    alist = SteenrodAlgebra(char)(alist)       # Coerce scalars if necessary
    if char == 2:
        alist2 = [e[0] for e in alist]
        maxlength = max([0]+[len(e) for e in alist2])
        alist2 = [list(e) + (maxlength-len(e))*[0] for e in alist2]
        minprofile = [max([alist2[i][j] for i in range(len(alist2))]) \
                                                for j in range(maxlength)]
        minprofile = tuple(map(lambda xx: mod_p_log(xx,char),minprofile))
        return find_min_profile(minprofile,char)
    if char != 2:
        alistQ = [e[0][0] for e in alist]
        alistP = [e[0][1] for e in alist]
        maxlengthQ = max([0]+[len(e) for e in alistQ])
        maxlengthP = max([0]+[len(e) for e in alistP])
        alistQ = [list(e) + (maxlengthQ-len(e))*[0] for e in alistQ]
        alistP = [list(e) + (maxlengthP-len(e))*[0] for e in alistP]
        minprofileQ = [max([alistQ[i][j] for i in range(len(alistQ))]) \
                                            for j in range(maxlengthQ)]
        minprofileP = [max([alistP[i][j] for i in range(len(alistP))]) \
                                            for j in range(maxlengthP)]
        minprofileP = tuple(map(lambda xx: mod_p_log(xx,char),minprofileP))
        if not minprofileQ:
            minpQ=[]
        else:
            minpQ = [1]*(max(minprofileQ)+1)
            for j in minprofileQ:
                minpQ[j] = 2
        return find_min_profile((minprofileP,minpQ),char=char)


def enveloping_profile_elements(alist,char=2):
    """
    Finds the profile function for the smallest sub-Hopf algebra containing
    the list of elements passed. Entries of `alist' are elements of the Steenrod
    Algebra. Hence, each alist[i] is treated as a list. Accepts either a list of
    lists or tuples.

    INPUT:

    -  ``alist``  - A list of Steenrod Algebra elements.

    OUTPUT: The profile function for the minimum sub-algebra containing all the
    elements of `alist`.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: enveloping_profile_elements([Sq(2),Sq(4)])
        (3, 2, 1)
        sage: enveloping_profile_elements([Sq(2,1,2),Sq(7)])
        (3, 2, 2, 1)

    """
    if char == 2:
        alist2 = list(map(profile_ele,[x for x in alist if  x != 0]))
        if not alist2:
             return (0,)
        maxlength = max([len(e) for e in alist2])
        alist2 = [list(e) + (maxlength-len(e))*[0] for e in alist2]
        minprofile = tuple(max([alist2[i][j] for i in range(len(alist2))]) \
                                                for j in range(maxlength))
        return find_min_profile(minprofile)
    else:
        masterlist = [profile_ele(x,char) for x in alist if x != 0]
        alistP = [x[0] for x in masterlist]
        alistQ = [x[1] for x in masterlist]
        if not alistP and not alistQ:
            return ((0,),(0,))
        maxlengthQ = max([len(e) for e in alistQ])
        maxlengthP = max([len(e) for e in alistP])
        alistQ = [list(e) + (maxlengthQ-len(e))*[0] for e in alistQ]
        alistP = [list(e) + (maxlengthP-len(e))*[0] for e in alistP]
        minprofileQ = tuple(max([alistQ[i][j] for i in range(len(alistQ))]) \
                                                for j in range(maxlengthQ))
        minprofileP = tuple(max([alistP[i][j] for i in range(len(alistP))])\
                                                for j in range(maxlengthP))
        return find_min_profile((minprofileP,minprofileQ),char=char)

def enveloping_profile_profiles(alist,char=2):
    """
    Finds the profile function for the smallest sub-Hopf algebra containing
    the sub-algebras corresponding the list of profile functions passed. Accepts
    either a list of lists or tuples.

    INPUT:

    -  ``alist``  - A list of profile functions.

    OUTPUT: The profile function for the minimum sub-algebra containing
    the profile functions in `alist`.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: enveloping_profile_profiles([[1,2,3],[2,4,1,1]])
        (2, 4, 3, 2, 1)
        sage: enveloping_profile_profiles([[4],[1,2,1],[3,2,3]])
        (4, 3, 3, 2, 1)

    """
    if char == 2:
        alist2 = list(copy(alist))
        maxlength = max([len(e) for e in alist2])
        alist2 = [list(e) + (maxlength-len(e))*[0] for e in alist2]
        minprofile = tuple(max([alist2[i][j] for i in range(len(alist2))]) \
                                                for j in range(maxlength))
        return find_min_profile(minprofile)
    else:
        alistP = [copy(alist[i][0]) for i in range(len(alist))]
        alistQ = [copy(alist[i][1]) for i in range(len(alist))]
        maxlengthQ = max([len(e) for e in alistQ])
        maxlengthP = max([len(e) for e in alistP])
        alistQ = [list(e) + (maxlengthQ-len(e))*[0] for e in alistQ]
        alistP = [list(e) + (maxlengthP-len(e))*[0] for e in alistP]
        minprofileQ = tuple(max([alistQ[i][j] for i in range(len(alistQ))]) \
                                                for j in range(maxlengthQ))
        minprofileP = tuple(max([alistP[i][j] for i in range(len(alistP))]) \
                                                for j in range(maxlengthP))
        return find_min_profile((minprofileP,minprofileQ),char=char)


def valid(LL,char=2):
    """
    Determines if the list passed is a valid profile.
    ## When checking at odd primes, the `P`-part must be the 0th entry in LL,
    and the `Q`-part must be the 1st entry in LL.

    INPUT:

    -  ``LL``  - A list of non-negative integers.

    OUTPUT: True or False, depending on whether the list passed is a valid profile.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: valid([3,2,1])
        True
        sage: valid([1,2,3])
        False

    """
    if char == 2:
        L = [0] + list(LL) + [0]*(len(LL)) # Add 0 in beginning to keep rels correct
        value = True                       # so r - i works when i = 0
        for r in range(2,len(L)):
            for i in range(1,r):
                value = value and ((L[r] >= L[r-i] - i) or (L[r] >= L[i]))
        return value
    else:
        (alistP,alistQ) = (LL[0],LL[1])
        M = [0] + list(alistP) + [0]*len(alistP)
        L = list(alistQ) + [1]*(len(alistQ)+1)
        M = M + [0]*abs(len(M)-len(L))  # Pad so they're the same length
        L = L + [1]*abs(len(M)-len(L))
        value = valid(alistP,char=2)  # P part must satisfy same conditions, regardless of prime.
        for r in range(len(L)): # \tau's indexed at 0, unlike \xi 's
            if (L[r] == 1) and value:
                for i in range(r+1):
                    value = value and ((M[i] <= r -i) or (L[r-i] == 1))
        return value

def nextprof(p,n,char=2):
    """
    Takes a possible profile `p' and a base profile `n'. Returns the next
    profile in lexicographic order. After all valid profiles `p' of
    length = len(n) have been checked, n is increased. Intended for internal
    use only. The odd primary piece only alters the Q part of the profile. To
    increment the P part of a profile for odd primes, call nextprof with char =2.
    This works since the P part of the profile is identical to the profile
    function when char == 2.

    INPUT:

    -  ``p``  - A pseudo-profile function which is incremented lexicographically
    and checked for validity.

    -  ``n``  - The base pseudo-profile function.

    OUTPUT: The next lexicographic profile.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: nextprof([1,2],[1,2])
        [2, 2]
        sage: nextprof([2,2],[1,2])
        [1, 2, 1]

        sage: nextprof([2,2,3],[1,2,3])
        [3, 2, 3]
        sage: nextprof([3,2,3],[1,2,3])
        [1, 3, 3]

    """
    if char == 2:
        for i in range(len(p)):
            if p[i] < len(p):     # we increment here because we can without altering length
                p[i] += 1
                return p
            else:
                if i > len(n)-1:               # Past end of n, so reset to 0
                    p[i] = 0
                else:                          # inside n still, so reset to n
                    p[i] = n[i]
        return n + [0]*(len(p)-len(n)) + [1]    # fell off the end
    else:  # odd primes
        pQ = list(p[1])
        nQ = list(n[1])
        for i in range(len(pQ)):
            if pQ[i] < 2:
                pQ[i] += 1
                return pQ
            else:
                if i > len(nQ) -1:
                    pQ[i] = 1
                else:
                    pQ[i] = nQ[i]
        return nQ + [1]*(len(pQ)-len(nQ)) +[1]


def find_min_profile(prof,char=2):
    """
    Given a tuple of integers (a pseudo-profile), this function will
    output the smallest legal profile function containing it. This
    function combines the above functions, and is the only one intended
    for external use.

    INPUT:

    -  ``prof``  - A list or tuple of nonnegative integers.

    OUTPUT:

    -  ``p`` - A valid profile containing ``p``.

    EXAMPLES::

        sage: from sage.modules.fp_over_steenrod_algebra.profile import *
        sage: find_min_profile([1,2])
        (1, 2, 1)
        sage: find_min_profile([2,1])
        (2, 1)
        sage: find_min_profile([1,2,3])
        (1, 2, 3, 1, 1)

    """
    if char == 2:
        prof2 = list(prof)
        if not prof2:
            return (0,)
        r = 0
        for i in range(len(prof2)):
            if prof2[i] != 0:
                r = i
        n = [prof2[i] for i in range(r+1)]
        p = copy(list(n))
        while not valid(p,char):
            p = nextprof(p,n,char)
        return tuple(p)
    else:
        pP,pQ = list(prof[0]), list(prof[1])
        P = find_min_profile(pP,char=2)
        Q = copy(pQ)
        while not valid([P,Q],char):
            Q = nextprof([P,Q],[P,pQ],char)
        return (P,Q)

