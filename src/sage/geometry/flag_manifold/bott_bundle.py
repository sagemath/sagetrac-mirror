r"""
Cohomology calculator for some homogeneous vector bundles on (complex) flag manifolds.

We call Bott bundles to those homogeneous vector bundles on flag manifolds that can be decomposed as a direct sum of homogeneous vector bundles associated to irreducible representations. We can use Bott's theorem to compute the cohomology of such bundles. This code allows the user to create, sum and twist Bott bundles on flag manifolds. Data such as rank and the full cohomology can be computed.

AUTHORS:

- Jorge Caravantes and Alicia Tocino (2016). Initial version.

EXAMPLES:

Generate a homogeneous space::

    sage: from sage.geometry.flag_manifold.bott_bundle import *
    sage: G = GrassForBottBundles(1,4)
    sage: G.explain()
    This is the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space
    Consider the uiversal short exact sequence:
    0 --> S --> O^5 --> Q --> 0
    where rk(S)=2 and rk(Q)=3.
    These bundles can be defined with methods universal_subbundle() and universal_quotient_bundle().
    Their duals can be obtained with methods universal_subbundle_dual() and universal_quotient_bundle_dual()
    Method hyperplane_bundle_multiple(k) gives the k power of the hyperplane section.
    Method cotangent_bundle() gives the cotangent bundle.
    Method tangent_bundle() gives the tangent bundle.

Create the universal quotient bundle::

    sage: Q = G.universal_quotient_bundle()
    sage: print Q
    Homogeneous rank-3 vector bundle on the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space.
    It is the sum of the following irreducible homogeneous vector bundles:
    <BLANKLINE>
    1 time(s) the tensor product of: 
        Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th quotient of tautological subbundles,
        Schur functor of partition [1, 1, 0] of the dual of the last quotient of tautological subbundles,
    <BLANKLINE>

Create the universal subbundle::

    sage: S = G.universal_subbundle()
    sage: print S
    Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space.
    It is the sum of the following irreducible homogeneous vector bundles:
    <BLANKLINE>
    1 time(s) the tensor product of: 
        Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th quotient of tautological subbundles,
        Schur functor of partition [1, 1, 1] of the dual of the last quotient of tautological subbundles,
    <BLANKLINE>

Twist and sum of Bott bundles::

    sage: E = S*Q + Q
    sage: E
    Homogeneous rank-9 vector bundle on the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space.
    It is the sum of the following irreducible homogeneous vector bundles:
    <BLANKLINE>
    1 time(s) the tensor product of: 
        Schur functor of partition [2, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
        Schur functor of partition [2, 2, 1] of the dual of the last quotient of tautological subbundles,
    <BLANKLINE>
    1 time(s) the tensor product of: 
        Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
        Schur functor of partition [1, 1, 0] of the dual of the last quotient of tautological subbundles,
    <BLANKLINE>


Computing the cohomology of a Bott bundle::

    sage: E. cohomology()
    h^ 0  =  5
    All remaining cohomology is zero
    [5]

REFERENCES:

.. [FH] William Fulton, Joe Harris: "Representation Theory - A First Course", GTM, Springer 
    (ISBN 978-1-4612-0979-9), 2004    
.. [M] I.G. Macdonald: "Symmetric Functions and Hall Polynomials", Oxford Mathematical Monographs, 
    (ISBN: 0198504500)
.. [W]  Jerzy Weyman: Cohomology of Vector Bundles and Syzygies, Cambridge Tracts in 
    Mathematics(ISBN: 9780521621977), 2003
    


"""

#*****************************************************************************
#       Copyright (C) 2015 Jorge Caravantes <jorgecaravan@gmail.com>
#       Copyright (C) 2015 Alicia Tocino <aliciatocinosanchez@mat.ucm.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from copy import deepcopy
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.combinat.partition import Partitions


#####################################################################
# First of all, we include some algorithms that will be just useful for
# the main classes and methods
#####################################################################

def bott(alpha):
    r"""
    Bott's algorithm for the general linear group as described in [W]. 
    
    INPUT:
    
    - ``alpha`` -- a sequence of nonnegative integers representing the concatenated 
        partitions of the Schur functors of the quotients of adjacent tautological subbundles.
    
    OUTPUT:
    
    - A pair [i,nu] where:

        - i is possitive iff there is a nonvanishing cohomology group. Such group is H^i 
            by Bott's Theorem, and it is the Schur Functor associated to the partition nu of 
            the vector space of dimension equal to the length of alpha.
        - i=-1 if all cohomology groups vanish. In this case, nu has no relevant 
            information.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: bott([3,1,0,0,0,5,0])
        [-1, [3, 1, 2, 1, 1, 1, 0]]
        sage: bott([3,1,0,0,0,7,0])
        [4, [3, 3, 2, 1, 1, 1, 0]]
        sage: bott([2,1,0,0,0,7,0])
        [-1, [2, 3, 2, 1, 1, 1, 0]]
        sage: bott([2,1,0,0,0,8,0])
        [5, [3, 3, 2, 1, 1, 1, 0]]

    
    
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    
    """
    n = len(alpha) - 1
    nu = deepcopy(alpha)
    i = 0
    goon = True
    while goon:
        j = 0
        while j<=n-1 and nu[j] >= nu[j+1]:
            j += 1
        if j == n:
            return [i, nu]
        if nu[j]+1 == nu[j+1]:
            return [-1, nu]
        nu[j], nu[j+1] = nu[j+1]-1, nu[j]+1
        i += 1


def schur_dimension(l):
    r"""
    Given a partition l, returns the dimension of the Schur/Weyl
    functor of a vector space of dimension len(l)
    
    INPUT:
    
    - ``l`` -- a partition.
    
    OUTPUT:
    
    - The dimension of the Schur functor associated to 'l' of an 'l'-dimensional vector space
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: schur_dimension([4,0,0])
        15
        sage: binomial(6,2)
        15
        sage: schur_dimension([1,1,1,1])
        1
        sage: schur_dimension([3,3,3,3])
        1
        sage: schur_dimension([4,3,2,1])
        64
        sage: schur_dimension([4,3,2,1,0])
        1024

    
    
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    
    """
    n = len(l)
    d = 1
    for i in range(1, n):
        for j in range(i+1, n+1):
            d *= (j-i+l[i-1]-l[j-1]) / (j-i)
    return d




def simplify_decomposition(dec):
    r"""
    Simplifies a redundant decomposition of a Bott bundle

    INPUT:
    
    - ``dec`` -- A list whose elements are pairs of a positive integer (multiplicity) and a list of nonnegative integers.
    
    OUTPUT:

    - A list of pairs as in the INPUT where the second entries of the pairs are pairwise different. The function just sums up the multiplicities of the same list to simplify the data.

    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: simplify_decomposition([[1,[2,1,0,1,0]],[3,[2,1,0,1,0]]])
        [[4, [2, 1, 0, 1, 0]]]


    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    """
    i = 0
    while i < len(dec) - 1:
        j = i + 1
        while j < len(dec):
#            print 'd[',i,']=',dec[i],';  d[',j,']=',dec[j]
            if dec[i][1] == dec[j][1]:
                dec[i][0] = dec[i][0] + dec[j][0]
                dec.pop(j)
            else:
                j += 1
        i += 1
    return dec
    
                


def multiply_decomposition(n, dec):
    r"""
    Multiplies a decomposition (of a Bott bundle) times an integer
    
    INPUT:
    
    - ``n`` -- A nonnegative integer
    - ``dec`` -- A list whose elements are pairs of a positive integer (multiplicity) and a list of nonnegative integers.
    
    OUTPUT:

    - A list of pairs as in the INPUT where the all the multiplicities are multiplied by 'n'.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: multiply_decomposition(3,[[4,[2,1,0,1,0]],[2,[1,0,0,1,0]]])
        [[12, [2, 1, 0, 1, 0]], [6, [1, 0, 0, 1, 0]]]


    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    """
    return [[n*i[0], i[1]] for i in dec]



def twist_of_decompositions(dec1, dec2, a, n):
    r"""
    Exactly what is in __mul__ for bundles (i.e. tensor product), but taking and returning just the decompositions. 

    INPUT:
    
    The data corresponding to a pair of Bott bundles on a flag manifold:
    - ``dec1``, ``dec2`` -- Lists whose elements are pairs of a positive integer (multiplicity) 
    and a list of nonnegative integers.
    - ``a`` -- A list of nonnegative integers (the dimensions of the projective subspaces
    conforming the flags).
    - ``n`` -- A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    
    OUTPUT:

    - A list of pairs as in the INPUT obtained by Littlewood--Richardson rule.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: d1=[[1,[2,0,0,0]],[2,[0,0,1,0]]]
        sage: d2=[[1,[1,1,1,0]],[1,[1,0,1,1]]]
        sage: twist_of_decompositions(d1,d2,[1],3)
        [[1, [3, 1, 1, 0]],
         [1, [3, 0, 1, 1]],
         [1, [2, 1, 1, 1]],
         [2, [1, 1, 2, 0]],
         [2, [1, 1, 1, 1]],
         [2, [1, 0, 2, 1]]]


    .. TODO::

        Eventually, when it is checked, this function will be called in the method __mul__ instead of being written twice.

    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    """
    #We are twisting decomposed bundles. Using distributivity, we twist each
    #pair of summands. 
    #Given a summand if a bundle, it is a twist of Schur functors of 
    #quotient of adjacent tautological subbundles. Using commutativity of the 
    #twist, we can perform Littlewood--Richardson rule on each quotient of adjacent 
    #tautological subbundles and then use distributivity to glue all.
    import sage.libs.lrcalc.lrcalc as lrcalc
    r = len(a)
    result=[]
    for i in dec1:
        for j in dec2:
    #Given a pair of summands (i,j) we have to twist and decompose any pair 
    #of factors associated to the same range of the 
    #list of partitions. Then we put it all together.
            LR_Result = lrcalc.mult(i[1][0:a[0]+1], j[1][0:a[0]+1], a[0]+1)
            SimpleMult = []
            for l in LR_Result.items():
                SimpleMult.append([l[1]*i[0]*j[0], l[0][:]])
    #We have to append zeroes till we reach the length we need
            for s in SimpleMult:
#                print 's=',s,'; long=',len(s[1]),'; a[0]=',a[0]
                for index in range(len(s[1]), a[0]+1):
                    s[1].append(0)
    #We have just calculated the product of the first factors of i and j. 
    #Now the intermediate
            for k in range(1, r):
                LR_Result = lrcalc.mult(i[1][a[k-1]+1:a[k]+1], j[1][a[k-1]+1:a[k]+1], a[k]-a[k-1])
                NewList = []
                for l in LR_Result.items():
                    for s in SimpleMult:
                        aux = deepcopy(s[1])
                        aux.extend(l[0][:])
                        NewList.append([l[1]*s[0], aux])
                SimpleMult = NewList
                for s in SimpleMult:
                    for index in range(len(s[1]), a[k]+1):
                        s[1].append(0)
    #Finally, we compute the product of the last factors of i and j
            NewList = []
            LR_Result = lrcalc.mult(i[1][a[r-1]+1:n+1], j[1][a[r-1]+1:n+1], n-a[r-1])
            for l in LR_Result.items():
                for s in SimpleMult:
                    aux = deepcopy(s[1])
                    aux.extend(l[0][:])
                    NewList.append([l[1]*s[0], aux])
            SimpleMult = NewList
            for s in SimpleMult:
                for index in range(len(s[1]), n+1):
                    s[1].append(0)
            result.extend(SimpleMult)
            result = simplify_decomposition(result)
    return result





def plethysm_of_irred(lmbd, alpha, a, n, s):
    r"""
    Applies Schur functor to an irreducible Bott bundle. Mainly based on [M, Formula (8.9)].

    INPUT:
    
    - ``lmbd`` -- a partition
    - ``alpha`` -- A list of integers: the concatenated Schur functors on the quotient of adjacent tautological subbundles of the flag manifold. The resulting bundle is the twist of all this  Schur functors of bundles
    - ``a`` -- A list of nonnegative integers (the dimensions of the projective subspaces conforming the flags).
    - ``n`` -- A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    - ``universal_subbundle`` -- It should be SymmetricFunctions(QQ).schur()

    OUTPUT:

    - A list  whose elements are pairs of a positive integer (multiplicity)  and a list of nonnegative integers.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: s=SymmetricFunctions(QQ).schur()
        sage: plethysm_of_irred([2,1], [1,0,0,0], [1], 3, s)
        [[1, [2, 1, 0, 0]]]
        sage: plethysm_of_irred([2,1], [0,0,1,0], [1], 3, s)
        [[1, [0, 0, 2, 1]]]
        sage: plethysm_of_irred([2,1], [1,0,1,0],[1],3,s)
        [[1, [2, 1, 3, 0]], [1, [3, 0, 2, 1]], [1, [2, 1, 2, 1]]]
        sage: plethysm_of_irred([4,0], [1,0,1,0], [1], 3, s)
        [[1, [4, 0, 4, 0]], [1, [3, 1, 3, 1]], [1, [2, 2, 2, 2]]]
        sage: plethysm_of_irred([4,0], [2,0,1,1], [1], 3, s)
        [[1, [6, 2, 4, 4]], [1, [4, 4, 4, 4]], [1, [8, 0, 4, 4]]]
    
    
    
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    """
    # The code is based on formula (8.9) of Macdonald's book:
    # [M] "Symmetric Functions and Hall Polynomials".
    if len(a) < 1:
        # In this case, the bundle is te image of an universal bundle by
        # alpha's Schur functor.
        # Therefore, we perform the plethysm with just one partition
        p = s(lmbd).plethysm(s(alpha))
        # Once we have obtained the image of the pletysm, we extract the coefficients 
        # and the partitions
        p = [p.coefficients(), [el[:] for el in p.support()]]
        # Since some of these partitions are longer than the dimension of the original bundle,
        # we have to remove them (for such cases, we get the zero bundle).
        # With partitions that are too short, we have to complete with zeros until we 
        # reach the appropriate length 
        l = 0
        while l < len(p[1]):
            if len(p[1][l]) > len(alpha):
                p[1].pop(l)
                p[0].pop(l)
            else:
                p[1][l].extend([0 for i in range(len(alpha)-len(p[1][l]))])
                l += 1
        return [[p[0][l], p[1][l]] for l in range(len(p[1]))]
    # If our bundle is a twist of several of the first case, we use formula (8.9) of [M]
    Prt = [el[:] for el in Partitions(sum(lmbd))]
    result = []
    for mu in Prt:
        IntPr = s(lmbd).internal_product(s(mu))
        PletInic = IntPr.plethysm(s(alpha[:a[0]+1]))
        PletInic = [PletInic.coefficients(), [el[:] for el in PletInic.support()]]
        # Again we remove too long partitions and we complete with zeros too short ones
        l = 0
        while l < len(PletInic[1]):
            if len(PletInic[1][l]) > a[0]+1:
                PletInic[1].pop(l)
                PletInic[0].pop(l)
            else:
                PletInic[1][l].extend([0 for i in range(a[0]+1-len(PletInic[1][l]))])
                l += 1
        # We have now the first factor of the summand (which is also a sum). 
        # Now, we use recursion to compute the second one
        PletSegund = plethysm_of_irred(mu, alpha[a[0]+1:], [a[i]-a[0] for i in range(1, len(a))], n-a[0]-1, s)
        for l in range(len(PletInic[1])):
            for m in range(len(PletSegund)):
                result.append([PletInic[0][l]*PletSegund[m][0], PletInic[1][l][:]+PletSegund[m][1][:]])
    return simplify_decomposition(result)
            


def plethysmization(lmbd, dec, a, n, s):
    r"""
    Tool to compute the decomposition of the Schur functor of a Bundle. Based on [M, Formula (8.8)] and function plethism_of_irred.
    
    INPUT:
    
    - ``lmbd`` -- a partition
    - ``dec`` -- A list whose elements are pairs of a positive integer (multiplicity) and a list of nonnegative integers.
    - ``a`` -- A list of nonnegative integers (the dimensions of the projective subspaces conforming the flags).
    - ``n`` -- A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    - ``universal_subbundle`` -- It should be SymmetricFunctions(QQ).schur()

    OUTPUT:

    - A list  whose elements are pairs of a positive integer (multiplicity) and a list of nonnegative integers.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: s=SymmetricFunctions(QQ).schur()
        sage: plethysmization([2,1],[[1, [4, 0, 4, 0]], [1, [3, 1, 3, 1]], [1, [2, 2, 2, 2]]],[1],3,s)
        [[37, [8, 4, 8, 4]],
         [11, [8, 4, 6, 6]],
         [12, [8, 4, 10, 2]],
         [1, [8, 4, 12, 0]],
         [25, [8, 4, 9, 3]],
         [12, [10, 2, 8, 4]],
         [5, [10, 2, 6, 6]],
         [8, [10, 2, 10, 2]],
         [1, [10, 2, 12, 0]],
         [12, [10, 2, 9, 3]],
         [7, [11, 1, 8, 4]],
         [3, [11, 1, 6, 6]],
         [6, [11, 1, 10, 2]],
         [2, [11, 1, 12, 0]],
         [9, [11, 1, 9, 3]],
         [42, [7, 5, 8, 4]],
         [13, [7, 5, 6, 6]],
         [12, [7, 5, 10, 2]],
         [1, [7, 5, 12, 0]],
         [30, [7, 5, 9, 3]],
         [36, [9, 3, 8, 4]],
         [7, [9, 3, 6, 6]],
         [11, [9, 3, 10, 2]],
         [1, [9, 3, 12, 0]],
         [19, [9, 3, 9, 3]],
         [15, [6, 6, 8, 4]],
         [5, [6, 6, 10, 2]],
         [3, [6, 6, 11, 1]],
         [23, [6, 6, 7, 5]],
         [7, [6, 6, 9, 3]],
         [5, [8, 4, 11, 1]],
         [53, [8, 4, 7, 5]],
         [4, [10, 2, 11, 1]],
         [11, [10, 2, 7, 5]],
         [5, [11, 1, 11, 1]],
         [7, [11, 1, 7, 5]],
         [7, [7, 5, 11, 1]],
         [76, [7, 5, 7, 5]],
         [3, [12, 0, 8, 4]],
         [3, [12, 0, 10, 2]],
         [6, [12, 0, 11, 1]],
         [3, [12, 0, 7, 5]],
         [3, [12, 0, 9, 3]],
         [5, [9, 3, 11, 1]],
         [58, [9, 3, 7, 5]],
         [4, [6, 6, 6, 6]]]

        
    
    .. WARNING::
    
       INTERNAL FUNCTION! DO NOT USE DIRECTLY!
    
       It uses ``sage.combinat.partition.Partitions``, which in turn uses ``sage.combinat.integer_list.IntegerListsLex``. This last resource seemed bugged when this function was coded.
    
            
    """
    if len(dec) < 2:
        return multiply_decomposition(dec[0][0], plethysm_of_irred(lmbd, dec[0][1], a, n, s))
    P = []
    # First of all, we get all partitions "less than" ``lmbd``, to be included in the sum
    for i in range(sum(lmbd)+1):
        P = P+[el[:] for el in Partitions(i,outer=lmbd).list()]
    # We now get skew Young diagrams associated to suci partitions.
    sups = []
    for p in P:
        sk = s(lmbd).skew_by(s(p))
        sups.append([p,[el[:] for el in sk.support()], sk.coefficients()])
    # Now we use  formula (8.8) of Macdonald's book:
    #"Symmetric Functions and Hall Polynomials".
    P = []
    for p in sups:
        auxdec = []
        # First we apply plethysm to the first summand  
        for i in range(len(p[1])):
            PlI = plethysm_of_irred(p[1][i], dec[0][1], a, n, s)
            auxdec.extend(multiply_decomposition(p[2][i], PlI))
        auxdec = simplify_decomposition(auxdec)
        # We have now the plethysm done to the first summand. 
        # It will be twisted with all other summands.
        # Since the coefficient of this summand is not necessarily 1, it must also be twisted
        # a number of times according to this. Moreover, the copies must be twisted among 
        # themselves
        if dec[0][0] == 1:
            remdec = dec[1:]
        else:
            remdec = [[dec[0,0]-1, dec[0,1]]] + dec[1:]
        # Now the recursive command, remdec is the remaining part of dec
        # after removing the first summand:
        P.extend(twist_of_decompositions(auxdec, plethysmization(p[0], remdec, a, n, s), a, n))
    P=simplify_decomposition(P)
    return P



###################################################################
#Now the class we use for reducible bundles. 
#It intends to be defined for arbitrary flag manifolds.
###################################################################


class BottBundle(SageObject):
    r"""
    This class represents homogeneous vector bundles on a given homogeneous space (flag manifold) that can be decomposed as a sum of irreducible homogeneous vector bundles whose cohomology can be computed by means of Bott's Theorem. 
    
    Such irreducible bundles are twists of Schur functors on the (duals of) quotients of tautological subbundles that are adjacent in the sequence (see [W] and [FH]):

    $0\subset R_{a_0+1}\subset...\subset R_{a_l+1}\subset\mathbb{C}^{n+1}\times\mathrm{Flag}(a_0,...,a_l;n)$
    
    Attributes:

    - ``space`` -- The ``FlagMfldForBottBundles`` on which the bundle is defined

    - ``decomposition`` -- A list whose elements are pairs of a positive integer (multiplicity) and 
        a list of nonnegative integers (irreducible ``BottBundle``).

    - ``aspect`` -- A simpler way to represent the bundle that keeps the way it was constructed.

    Methods:
    
    - Addition is equivalent to direct sum of vector bundles.

    - Multiplication is tensor product.

    - ``cohomology`` -- Computes the whole sheaf cohomology of ``self``.

    - ``rank`` -- returns the rank of ``self``.

    - ``schur`` -- performs a given partition's schur functor of ``self``. 

    - ``sym`` -- a given integer symmetric power of ``self``.

    - ``wedge`` -- a given number exterior power of ``self``.

    INPUT:

    - ``space`` -- a ``FlagMfldForBottBundles``.

    OUTPUT:

    - The zero bundle on ``space``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: F=FlagMfldForBottBundles([1,3,5],7)
        sage: BottBundle(F)
        Homogeneous rank-0 vector bundle on the space of flags of subspaces of dimensions [1, 3, 5] in the 
        7-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
    
    
    WARNING:
    
    The methods ``schur``, ``sym`` and ``wedge`` depend on ``sage.combinat.partition.Partitions``, which in turn uses ``sage.combinat.integer_list.IntegerListsLex``. This last resource seemed bugged when this function was coded. 
    
    .. TODO::
    
        Method _repr_ should be different when self.space is a grassmannian (it should speak about universal bundles) or a projective space (about tangent bundle and hyperplane section).

        A function to compute a concrete cohomology group.

        The attribute ``aspect`` has yet to be adapted to a non ``GrassForBottBundles`` ``FlagMfldForBottBundles``.

        A function that computes the dual of a ``BottBundle``.  
    
    

    """

    def __init__(self, space):
        r"""
        Initialize ``self``
        
        INPUT:
        - ``space`` -- a value of type FlagMfldForBottBundles.
        
        EXAMPLES::
        
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: F=FlagMfldForBottBundles([1,3,5],7)
            sage: BottBundle(F)
            Homogeneous rank-0 vector bundle on the space of flags of subspaces of dimensions [1, 3, 5] in the 
            7-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
        
        """
        self.decomposition = []
        self.space = space
        self.aspect = '0'





    def _repr_(self):
        r"""
        Represent ``self``

        TEST::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: F=FlagMfldForBottBundles([1,3,5],7)
            sage: BottBundle(F)
            Homogeneous rank-0 vector bundle on the space of flags of subspaces of dimensions [1, 3, 5] in the 
            7-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:

        .. TODO::
        
            Special names when ``self.space`` is a grassmannian or a projective space.
        """
        cadena = 'Homogeneous rank-' + str(self.rank()) + ' vector bundle on ' + self.space._repr_()+'.\n' + 'It is the sum of the following irreducible homogeneous vector bundles:' + '\n'
        for i in self.decomposition:
            particiones = ''
            ant = 0
            contador = 1
            for d in self.space.a:
                particiones = particiones + '    Schur functor of partition ' + str(i[1][ant:d+1]) + ' of the dual of the ' + str(contador) + 'st/nd/rd/th  quotient of tautological subbundles,\n'
                contador += 1
                ant = d + 1
            particiones = particiones + '    Schur functor of partition ' + str(i[1][ant:self.space.n+1]) + ' of the dual of the last quotient of tautological subbundles,\n'
            cadena = cadena + '\n' + str(i[0]) + ' time(s) the tensor product of: \n' + particiones
        return cadena 

    def __add__(self, other):
        r"""
        Direct sum of vector bundles ``self`` and ``other``.
        
        INPUT: 

        - ``self``, ``other`` -- of BottBundle class, with the same attribute ``space``.

        OUTPUT:
        
        - The direct sum of ``self`` and ``other``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: F=G.universal_subbundle()
            sage: E+F
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        if self.space <> other.space:
            raise ValueError("Cannot sum bundles on different varieties")
            return
        result = BottBundle(self.space)
        result.aspect = self.aspect + '(+)' + other.aspect
        result.decomposition = simplify_decomposition(deepcopy(self.decomposition) + deepcopy(other.decomposition))
        return result

        


    def __mul__(self, other):
        r"""
        Twists ``self`` times ``other``.

        INPUT:

        - ``self``, ``other`` -- of BottBundle class, with the same attribute ``space``.

        OUTPUT:

        - The tensor product of ``self`` and ``other``.
  
        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: F=G.universal_subbundle()
            sage: E*F
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [2, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        if self.space <> other.space:
            raise ValueError("Cannot twist bundles on different varieties")
            return
        import sage.libs.lrcalc.lrcalc as lrcalc
        n = self.space.n
        a = deepcopy(self.space.a)
        r = len(a)
        result = BottBundle(self.space)
        result.aspect = '(' + self.aspect + ')(x)(' + other.aspect + ')'
    #EVENTUALLY, THE FOLLOWING CODE SHOULD BE REMOVED TO PLACE A CALL TO 
    #THE PROCEDURE twist_of_decompositions
    #Here we take the decompositions of both factors and we 
    #will twist any pair of summands
        for i in self.decomposition:
            for j in other.decomposition:
    #Given a pair of summands (i,j) we have to twist and decompose 
    #any pair of factors associated to the same range of the 
    #list of partitions. Then we put it all together.
                LR_Result = lrcalc.mult(i[1][0:a[0]+1], j[1][0:a[0]+1], a[0]+1)
                SimpleMult = []
                for l in LR_Result.items():
                    SimpleMult.append([l[1]*i[0]*j[0], l[0][:]])
    #We have to append zeroes till we reach the length we need
                for s in SimpleMult:
#                    print 's=',s,'; long=',len(s[1]),'; a[0]=',a[0]
                    for index in range(len(s[1]), a[0]+1):
                        s[1].append(0)
    #We have just calculated the product of the first factors of 
    #i and j. Now the intermediate
                for k in range(1, r):
                    LR_Result = lrcalc.mult(i[1][a[k-1]+1:a[k]+1], j[1][a[k-1]+1:a[k]+1], a[k]-a[k-1])
                    NewList = []
                    for l in LR_Result.items():
                        for s in SimpleMult:
                            aux = deepcopy(s[1])
                            aux.extend(l[0][:])
                            NewList.append([l[1]*s[0],aux])
                    SimpleMult = NewList
                    for s in SimpleMult:
                        for index in range(len(s[1]), a[k]+1):
                            s[1].append(0)
    #Finally, we compute the product of the last factors of i and j
                NewList = []
                LR_Result = lrcalc.mult(i[1][a[r-1]+1:n+1], j[1][a[r-1]+1:n+1], n-a[r-1])
                for l in LR_Result.items():
                    for s in SimpleMult:
                        aux = deepcopy(s[1])
                        aux.extend(l[0][:])
                        NewList.append([l[1]*s[0], aux])
                SimpleMult = NewList
                for s in SimpleMult:
                    for index in range(len(s[1]), n+1):
                        s[1].append(0)
                result.decomposition.extend(SimpleMult)
                result.decomposition = simplify_decomposition(result.decomposition)
        return result    
        


    def schur(self, lmbd):
        r"""
        A given partition's schur functor of ``self``.

        INPUT:

        - ``lmbd`` -- a partition (decreasin list of possitive integers)

        OUTPUT:
    
        - the Schur functor associated to partition ``lmbd`` applied to ``self``. 
        

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: (E*E).schur([2,1])
            Homogeneous rank-28 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>            
            4 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            3 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [5, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        .. WARNING::
        
        
            It uses ``sage.combinat.partition.Partitions``, which in turn uses ``sage.combinat.integer_list.IntegerListsLex``. This last resource seemed bugged when this function was coded.

            Plethism computation is slow, so difficult Schur functors of complex bundles can take forever to compute.
        """
        result = BottBundle(self.space)
        result.aspect = 'Schur_' + lmbd.__str__() + '('+self.aspect+')'
        s = SymmetricFunctions(QQ).schur()
        result.decomposition=plethysmization(lmbd, self.decomposition, self.space.a, self.space.n, s)
        return result
            



    def wedge(self, k):
        r"""
        Computes the ``k``-th exterior power of ``self``
        
        INPUT:

        - ``k`` A positive integer

        OUTPUT:

        - the ``k``-th exterior power of ``self``

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: (E*E).wedge(3)
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>


        .. WARNING::
        
        
            It uses ``sage.combinat.partition.Partitions``, which in turn uses ``sage.combinat.integer_list.IntegerListsLex``. This last resource seemed bugged when this function was coded.

            Plethism computation is slow, so difficult Schur functors of complex bundles can take forever to compute.
        """
        return self.schur([1 for j in range(k)])


    def sym(self, k):            
        r"""
        Computes the ``k``-th symmetric power of ``self``

        INPUT:

        - ``k`` A possitive integer

        OUTPUT:

        - the ``k``-th symmetric power of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: (E*E).sym(3)
            Homogeneous rank-20 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            2 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [6, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            2 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [5, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        .. WARNING::
        
            It uses ``sage.combinat.partition.Partitions``, which in turn uses ``sage.combinat.integer_list.IntegerListsLex``. This last resource seemed bugged when this function was coded.

            Plethism computation is slow, so difficult Schur functors of complex bundles can take forever to compute.
        """
        return self.schur([k])


    def rank(self):
        r"""
        Computes the rank of the vector bundle.

        EXAMPLES::
    
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: E.rank()
            2
            sage: (E+E).rank()
            4
            sage: (E*E*E).rank()
            8
        """
        TotalRank = 0
        for d in self.decomposition:
            PartialRank = 1
            previousa = -1
            for a in self.space.a:
                PartialRank *= schur_dimension(d[1][previousa+1:a+1])
                previousa = a
            PartialRank *= schur_dimension(d[1][previousa+1:self.space.n+1])
            TotalRank += d[0] * PartialRank
        return TotalRank
                
        


    def cohomology(self):
        r"""
        Computes the complete cohomology of the vector bundle.

        OUTPUT:

        A list ``l`` of integers such that ``l[i]`` is the diension of the ``i``-th cohomology space of the bundle. The length of ``l`` is equal to the index of the last nonvanishing cohomology space (+1). Therefore, if the length is zero, all cohomology vanishes.
        
        EXAMPLES::
    
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: E=G.universal_quotient_bundle()
            sage: E.cohomology()
            h^ 0  =  4
            All remaining cohomology is zero
            [4]
            sage: F=G.hyperplane_bundle_multiple(-7)
            sage: (E*F).cohomology()
            h^ 4  =  60
            All remaining cohomology is zero
            [0, 0, 0, 0, 60]
            sage: (E*F+E).cohomology()
            h^ 0  =  4
            h^ 4  =  60
            All remaining cohomology is zero
            [4, 0, 0, 0, 60]

        .. TODO::

            A possible second argument (or a separate function to return a concrete cohomology space).

        """
        n = len(self.decomposition[0][1])
        B = []
        for el in self.decomposition:
            lista = deepcopy(el[1])
            aux = bott(lista)
            aux.append(el[0])            
            B.append(aux)
        B.sort()
        m = B[len(B)-1][0]
        Cohom = []
        for ind in range(0, m+1):
            Cohom.append(0)
        for el in B:
            if el[0] <> -1:
                Cohom[el[0]] += el[2] * schur_dimension(el[1])
        for i in range(0, m+1):
            if Cohom[i] <> 0:
                print "h^", i, " = ", Cohom[i]
        if len(Cohom) == 0:
            print 'All cohomology vanishes'
        else:
            print 'All remaining cohomology is zero'
        return Cohom




class UniversalSubbundleDual(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the dual to the universal subbundle of an instance of ``FlagMfldForBottBundles``

    INPUT:

    - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

    OUTPUT:

    - The dual to the universal subbundle of ``fm``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G=GrassForBottBundles(1,3)
        sage: UniversalSubbundleDual(G)
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: UniversalSubbundleDual(G)
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
        """
        self.decomposition = [[Integer(1),[Integer(1)]+[Integer(0)]*(fm.n)]]
        self.space = fm
        self.aspect = 'S_1^*'


class UniversalSubbundle(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the universal subbundle of an instance of ``GrassForBottBundles``

    INPUT:

    - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

    OUTPUT:

    - The universal subbundle of ``fm``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G=GrassForBottBundles(1,3)
        sage: UniversalSubbundle(G)                                      
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: UniversalSubbundle(G)                                      
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
        """
        self.decomposition = [[Integer(1), [Integer(1)]*(fm.a[0])+[Integer(0)]+[Integer(1)]*(fm.n-fm.a[0])]]
        self.space = fm
        self.aspect = 'S_1'

class UniversalQuotientDual(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the dual to the universal quotient bundle of an instance of ``GrassForBottBundles``

    INPUT:

    - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

    OUTPUT:

    - The dual to the universal quotient bundle of ``fm``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G=GrassForBottBundles(1,3)
        sage: UniversalQuotientDual(G)
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURGAGED
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: UniversalQuotientDual(G)
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
        """
        self.decomposition = [[Integer(1),[Integer(0)]*(fm.a[0]+1)+[Integer(1)]+[Integer(0)]*(fm.n-fm.a[0]-1)]]
        self.space = fm
        self.aspect = 'Q_1^*'

class UniversalQuotient(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the universal quotient bundle of an instance of ``FlagMfldForBottBundles``

    INPUT:

    - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

    OUTPUT:

    - The universal quotient bundle of ``fm``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G=GrassForBottBundles(1,3)
        sage: UniversalQuotient(G)    
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>


    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` -- an instance of ``FlagMfldForBottBundles``, it should be a ``GrassForBottBundles``

        OUTPUT:

        - The universal quotient bundle of ``fm``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: UniversalQuotient(G)    
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
        """
        self.decomposition = [[Integer(1), [Integer(1)]*(fm.n)+[Integer(0)]]]
        self.space = fm
        self.aspect = 'Q_1'

class LineBundleGrass(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize line bundles on grassmannians.

    INPUT:

    - ``grass`` -- an instance of ``FlagMfldForBottBundles`` it should be a ``GrassForBottBundles``

    - ``k`` -- an integer.

    OUTPUT:

    - the ``k``-th multiple of the hyperplane section of ``grass``

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G=GrassForBottBundles(1,3)
        sage: LineBundleGrass(G,3)
        Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>
        sage: LineBundleGrass(G,-2)
        Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    WILL NOT WORK PROPERLY WHEN ``grass`` IS NOT A GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
    """
    def __init__(self, grass, k):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` -- an instance of ``FlagMfldForBottBundles``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: sage.geometry.flag_manifold.bott_bundle.LineBundleGrass(G,3)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            sage: LineBundleGrass(G,-2)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS DISCOURAGED!
        """
        if k > 0:
            self.decomposition = [[Integer(1), [Integer(k)]*(grass.a[0]+1)+[Integer(0)]*(grass.n-grass.a[0])]]
        elif k < 0:
            self.decomposition = [[Integer(1), [Integer(0)]*(grass.a[0]+1)+[Integer(-k)]*(grass.n-grass.a[0])]]
        else:
            self.decomposition = [[Integer(1), [Integer(0)]*(grass.n+1)]]
        self.space = grass
        self.aspect = 'O(' + str(k) + ')'






####################################################################
#Now we introduce the classes referring the homogeneous spaces.
#Since we just use the special linear group, all are flag manifolds
####################################################################



class FlagMfldForBottBundles(SageObject):
    r"""
    This is the class containing all quotients of the general linear group by a parabolic subgroup. 
        
    Attributes:

    - ``a`` -- a list of increasing integers (dimensions of the subspaces of the flag)

    - ``n`` -- possitive integer (the dimension of the ambient projective space of the flag)
    
    
    INPUT:
    
    - ``a`` -- a list of increasing integers (dimensions of the subspaces of the flag)
    
    - ``n`` -- possitive integer (the dimension of the ambient projective space of the flag)
    
    OUTPUT:
    
    - The space of flags of subspaces of dimensions given by ``a`` in the projective space of dimension ``n``.
    
    EXAMPLES::
    
        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: FlagMfldForBottBundles([1,3,5],7)
        the space of flags of subspaces of dimensions [1, 3, 5] in the 7-dimensional projective space    
    
    WARNING:
    
    It only contains properties related with the class ``BottBundle``, so no more geometry is allowed with this class.
    
    .. TODO::
    
        Methods to generate special instances of the class ``BottBundle``, mainly irreducible ones and their duals
        
    """
    def __init__(self, a, n):
        r"""
        Initializes ``self``.

        INPUT:
    
        - ``a`` -- a list of increasing integers (dimensions of the subspaces of the flag)
    
        - ``n`` -- possitive integer (the dimension of the ambient projective space of the flag)
    
        OUTPUT:
    
        - The space of flags of subspaces of dimensions given by ``a`` in the projective space of 
        dimension ``n``.
    
        EXAMPLES::
        
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: FlagMfldForBottBundles([1,3,5],7)
            the space of flags of subspaces of dimensions [1, 3, 5] in the 7-dimensional projective space
    
        """
        if type(a) <> list:
            raise ValueError("First argument must be a nonempty sorted list of nonnegative integers")
        NonInt=False
        for el in a:
            if type(el)<>Integer:
                print type(el), Integer
                NonInt=True
        if NonInt:
            raise ValueError("First argument must be a nonempty sorted list of nonnegative integers")
        if a <> sorted(a):
            raise ValueError("First argument must be a nonempty sorted list of nonnegative integers")
        if len(a) == 0:
            raise ValueError("First argument must be a nonempty sorted list of nonnegative integers")
        if a[0] < 0:
            raise ValueError("First argument must be a nonempty sorted list of nonnegative integers")
        if type(n)<>Integer:
            raise ValueError("Second argument must be a nonnegative integer")
        if a[len(a)-1] > n:
            raise ValueError("Last element of first argument cannot be greater than second argument")
        self.a = a
        self.n = n

    def _repr_(self):
        r"""
        Represents ``self``.
    
        EXAMPLES::
        
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: FlagMfldForBottBundles([1,3,5],7)
            the space of flags of subspaces of dimensions [1, 3, 5] in the 7-dimensional projective space
    
        """
        return 'the space of flags of subspaces of dimensions ' + str(self.a) + ' in the ' + str(self.n) + '-dimensional projective space'

    def is_grassmannian(self):
        r"""
        Determines wether ``self`` is a grassmannian.
    
        EXAMPLES::
        
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: FlagMfldForBottBundles([1,3,5],7).is_grassmannian()                          
            False
            sage: FlagMfldForBottBundles([5],7).is_grassmannian()    
            True
            sage: FlagMfldForBottBundles([0],7).is_grassmannian()
            True
    
        """
        if len(self.a) == 1:
            return True
        else:
            return False

    def is_projective_space(self):
        r"""
        Determines wether ``self`` is a projective space.
    
        EXAMPLES::
        
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: FlagMfldForBottBundles([1,3,5],7).is_projective_space()
            False
            sage: FlagMfldForBottBundles([5],7).is_projective_space()    
            False
            sage: FlagMfldForBottBundles([0],7).is_projective_space()
            True

    
        """
        if self.is_grassmannian():
            if self.a[0] == 0:
                return True
            elif self.a[0] >= self.n-1:
                return True
            else:
                return False
        else:
            return False



class GrassForBottBundles(FlagMfldForBottBundles):
    r"""
    A child class of ``FlagMfldForBottBundles`` to easily initialize grassmannians.

    Methods:

    - ``hyperplane_bundle_multiple`` -- creates a multiple of the hyperplane section

    - ``universal_quotient_bundle`` -- creates the universal quotient bundle

    - ``universal_quotient_bundle_dual`` -- creates the dual of the universal quotient bundle

    - ``universal_subbundle`` -- creates the universal subbundle

    - ``universal_subbundle_dual`` -- creates the dual of the universal subbundle

    - ``tangent_bundle`` -- Creates the tangent bundle

    - ``cotangent_bundle`` -- Creates the cotangent bundle


    INPUT:

    - ``k`` -- a list of increasing integers (dimension of the projective subspace)

    - ``n`` -- possitive integer (the dimension of the ambient projective space of the flag)

    OUTPUT:

    - The grassmannian of ``k``-dimensional projective subspaces in the projective space of dimension ``n``.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: G1=GrassForBottBundles(4,9)                                               
        sage: G2=GrassForBottBundles(2,10)
        sage: G1,G2
        (the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space,
         the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space)

    WARNING:

    It only contains properties related with the class ``BottBundle``, so no more geometry is allowed with this class.

    """
    def __init__(self,k,n):
        r"""
        Initializes ``self``

        INPUT:
    
        - ``k`` -- a list of increasing integers (dimension of the projective subspace)
    
        - ``n`` -- possitive integer (the dimension of the ambient projective space of the flag)
    
        OUTPUT:
    
        - The grassmannian of ``k``-dimensional projective subspaces in the projective space of 
        dimension ``n``.
    
        EXAMPLES::
    
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G1=GrassForBottBundles(4,9)                                               
            sage: G2=GrassForBottBundles(2,10)
            sage: G1,G2
            (the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space,
             the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space)
    
        """
        if type(k)<>Integer or k < 0:
            raise ValueError("First argument must be a nonnegative integer")
        if type(n)<>Integer:
            raise ValueError("Second argument must be a nonnegative integer")
        if k > n:
            raise ValueError("First argument cannot be greater than second argument")
        self.a = [k]
        self.n = n

    def explain(self):
        r"""
        Explains what ``self`` is. Just considering

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(4,9)                                               
            sage: G.explain()
            This is the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^10 --> Q --> 0
            where rk(S)=5 and rk(Q)=5.
            These bundles can be defined with methods universal_subbundle() and universal_quotient_bundle().
            Their duals can be obtained with methods universal_subbundle_dual() and universal_quotient_bundle_dual()
            Method hyperplane_bundle_multiple(k) gives the k power of the hyperplane section.
            Method cotangent_bundle() gives the cotangent bundle.
            Method tangent_bundle() gives the tangent bundle.        

        """
        print 'This is the grassmannian of ' + str(self.a[0]) + '-dimensional subspaces of the ' + str(self.n) + '-dimensional projective space'
        print 'Consider the uiversal short exact sequence:'
        print '0 --> S --> O^' + str(self.n+1) + ' --> Q --> 0'
        print 'where rk(S)=' + str(self.a[0]+1) + ' and rk(Q)=' + str(self.n-self.a[0]) + '.'
        print 'These bundles can be defined with methods universal_subbundle() and universal_quotient_bundle().'
        print 'Their duals can be obtained with methods universal_subbundle_dual() and universal_quotient_bundle_dual()'
        print 'Method hyperplane_bundle_multiple(k) gives the k power of the hyperplane section.'
        print 'Method cotangent_bundle() gives the cotangent bundle.'
        print 'Method tangent_bundle() gives the tangent bundle.'


    def _repr_(self):
        r"""
        Represents ``self``
    
        EXAMPLES::
    
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G1=GrassForBottBundles(4,9)                                               
            sage: G2=GrassForBottBundles(2,10)
            sage: G1,G2
            (the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space,
             the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space)
    
        """
        return 'the grassmannian of ' + str(self.a[0]) + '-dimensional subspaces of the ' + str(self.n) + '-dimensional projective space'

    def hyperplane_bundle_multiple(self, k=Integer(0)):
        r"""
        k-th power of the hyperplane section of the grassmannian

        INPUT:

        - ``k`` -- an integer.

        OUTPUT:

        - the ``k``-th multiple of the hyperplane section of ``self``

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.hyperplane_bundle_multiple(3)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>    
            1 time(s) the tensor product of: 
                Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            sage: G.hyperplane_bundle_multiple(-2)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
            sage: G.hyperplane_bundle_multiple()
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        if type(k)<>Integer:
            raise ValueError("Argument must be an integer")
        return LineBundleGrass(self, k)

    def universal_quotient_bundle(self):
        r"""
        Universal quotient bundle of the grassmannian

        OUTPUT:

        - The universal quotient bundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.universal_quotient_bundle()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return UniversalQuotient(self)

    def universal_subbundle(self):
        r"""
        Universal subbundle of the grassmannian

        OUTPUT:

        - The universal subbundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.universal_subbundle()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        """
        return UniversalSubbundle(self)

    def universal_quotient_bundle_dual(self):
        r"""
        Dual of the universal quotient bundle of the grassmannian

        OUTPUT:

        - The dual to the universal quotient bundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.universal_quotient_bundle_dual()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return UniversalQuotientDual(self)

    def universal_subbundle_dual(self):
        r"""
        Dual of the universal subbundle of the grassmannian

        OUTPUT:

        - The dual to the universal subbundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.universal_subbundle_dual()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return UniversalSubbundleDual(self)

    def cotangent_bundle(self):
        r"""
        Cotangent bundle of the grassmannian

        OUTPUT:

        - The cotangent bundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.cotangent_bundle()
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return self.universal_subbundle()*self.universal_quotient_bundle_dual()

    def tangent_bundle(self):
        r"""
        Tangent bundle of the grassmannian

        OUTPUT:

        - The tangent bundle of ``self``.

        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: G=GrassForBottBundles(1,3)
            sage: G.tangent_bundle()
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [2, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return self.universal_subbundle_dual()*self.universal_quotient_bundle()





class ProjForBottBundles(GrassForBottBundles):
    r"""
    A child class of ``GrassForBottBundles`` to initialize easily projective spaces.

    INPUT:

    - ``n`` -- an integer.

    OUTPUT:

    - The ``n``-dimensional projective space.

    EXAMPLES::

        sage: from sage.geometry.flag_manifold.bott_bundle import *
        sage: ProjForBottBundles(3)
        the projective space of dimension 3
        sage: ProjForBottBundles(7)
        the projective space of dimension 7

    WARNING:

    It only contains properties related with the class ``BottBundle``, so no more geometry is allowed with this class.

    """
    def __init__(self, n):
        r"""
        Initializes ``self``

        EXAMPLES::
    
            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: sage.geometry.flag_manifold.bott_bundle.ProjForBottBundles(3)
            the projective space of dimension 3
            
        """
        if type(n)<>Integer or n < 0:
            raise ValueError("Argument must be a nonnegative integer")
        self.a = [0]
        self.n = n


    def explain_bott_bundles(self):
        r"""


        EXAMPLES::

            sage: from sage.geometry.flag_manifold.bott_bundle import *
            sage: P = sage.geometry.flag_manifold.bott_bundle.ProjForBottBundles(3)
            sage: P.explain_bott_bundles()
            This is the projective space of dimension 3.
            Method hyperplane_bundle_multiple(k) gives the k power of the hyperplane section.
            Method cotangent_bundle() gives the cotangent bundle.
            Method tangent_bundle() gives the tangent bundle.
        """
        print 'This is the projective space of dimension ' + str(self.n) + '.'
        print 'Method hyperplane_bundle_multiple(k) gives the k power of the hyperplane section.'
        print 'Method cotangent_bundle() gives the cotangent bundle.'
        print 'Method tangent_bundle() gives the tangent bundle.'


    def _repr_(self):
        r"""
        Represents ``self``
        
        TEST::
        
            sage: P=sage.geometry.flag_manifold.bott_bundle.ProjForBottBundles(3)
            sage: P
            the projective space of dimension 3
        """
        return 'the projective space of dimension ' + str(self.n)







