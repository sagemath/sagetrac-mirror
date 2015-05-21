r"""
Cohomology calculator for some homogeneous vector bundles on (complex) flag manifolds.

We call Bott bundles to those homogeneous vector bundles on flag manifolds that can be decomposed as a direct sum of homogeneous vector bundles associated to irreducible representations. We can use Bott's theorem to compute the cohomology of such bundles. This code allows the user to create, sum and twist Bott bundles on flag manifolds. Data such as rank and the full cohomology can be computed.

AUTHORS:

- Jorge Caravantes and Alicia Tocino (2015): initial version


EXAMPLES:

Generate a homogeneous space::

    sage: G = Grassmannian(1,4)
    You have defined the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space
    Consider the uiversal short exact sequence:
    0 --> S --> O^5 --> Q --> 0
    where rk(S)=2 and rk(Q)=3.
    These bundles can be defined with methods S() and Q().
    Their duals can be obtained with methods SDual() and Qdual()
    Method O(k) gives the k power of the hyperplane section.
    Method Om() gives the cotangent bundle.
    Method T() gives the tangent bundle.

Create the universal quotient bundle::

    sage: Q = G.Q()
    sage: print Q
    Homogeneous rank-3 vector bundle on the grassmannian of 1-dimensional subspaces of the 4-dimensional projective space.
    It is the sum of the following irreducible homogeneous vector bundles:
    <BLANKLINE>
    1 time(s) the tensor product of: 
        Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th quotient of tautological subbundles,
        Schur functor of partition [1, 1, 0] of the dual of the last quotient of tautological subbundles,
    <BLANKLINE>

Create the universal subbundle::

    sage: S = G.S()
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


from sage.rings.integer import Integer
from copy import deepcopy
from sage.combinat.sf.sf import SymmetricFunctions
from sage.rings.rational_field import QQ
from sage.combinat.partition import Partitions


#####################################################################
#First of all, we include some algorithms that will be just useful for
# the main classes and methods
#####################################################################

def _bott(alpha):
    r"""
    Bott's algorithm for the general linear group as described in [W]. 
    
    INPUT:
    
    - ``alpha`` - a sequence of nonnegative integers representing the concatenated 
    partitions of the Schur functors of the quotients of adjacent tautological subbundles.
    
    OUTPUT:
    
    - A pair [i,nu] where:
        - i is possitive iff there is a nonvanishing cohomology group. Such group is H^i 
        by Bott's Theorem, and it is the Schur Functor associated to the partition nu of 
        the vector space of dimension equal to the length of alpha.
        - i=-1 if all cohomology groups vanish. In this case, nu has no relevant 
        information.
    
    EXAMPLES::
    
        sage: _bott([3,1,0,0,0,5,0])
        [-1, [3, 1, 2, 1, 1, 1, 0]]
        sage: _bott([3,1,0,0,0,7,0])
        [4, [3, 3, 2, 1, 1, 1, 0]]
        sage: _bott([2,1,0,0,0,7,0])
        [-1, [2, 3, 2, 1, 1, 1, 0]]
        sage: _bott([2,1,0,0,0,8,0])
        [5, [3, 3, 2, 1, 1, 1, 0]]

    
    
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    
    REFERENCES:
    
    .. [W]  Jerzy Weyman: Cohomology of Vector Bundles and Syzygies, Cambridge Tracts in 
    Mathematics(ISBN: 9780521621977), 2003
    
    """
    n = len(alpha) - 1
    nu = alpha
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
#        aux=nu[j]
#        nu[j]=nu[j+1]-1
#        nu[j+1]=aux+1
        i += 1


def _schur_dimension(l):
    r"""
    Given a partition l, returns the dimension of the Schur/Weyl
    functor of a vector space of dimension len(l)
    
    INPUT:
    
    - ``l`` - a partition.
    
    OUTPUT:
    
    - The dimension of the Schur functor associated to 'l' of an 'l'-dimensional vector space
    
    EXAMPLES::
    
        sage: _schur_dimension([4,0,0])
        15
        sage: binomial(6,2)
        15
        sage: _schur_dimension([1,1,1,1])
        1
        sage: _schur_dimension([3,3,3,3])
        1
        sage: _schur_dimension([4,3,2,1])
        64
        sage: _schur_dimension([4,3,2,1,0])
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




def _simplify_decomposition(dec):
    r"""
    Simplifies a redundant decomposition of a Bott bundle

    INPUT:
    
    - ``dec`` - A list whose elements are pairs of a positive integer (multiplicity) and a list of
    nonnegative integers.
    
    OUTPUT:

    - A list of pairs as in the INPUT where the second entries of the pairs are pairwise different.
    The function just sums up the multiplicities of the same list to simplify the data.

    EXAMPLES::
    
        sage: _simplify_decomposition([[1,[2,1,0,1,0]],[3,[2,1,0,1,0]]])
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
    
                


def _multiply_decomposition(n, dec):
    r"""
    Multiplies a decomposition (of a Bott bundle) times an integer
    
    INPUT:
    
    - ``n`` - A nonnegative integer
    - ``dec`` - A list whose elements are pairs of a positive integer (multiplicity) and a list of
    nonnegative integers.
    
    OUTPUT:

    - A list of pairs as in the INPUT where the all the multiplicities are multiplied by 'n'.
    
    EXAMPLES::
    
        sage: _multiply_decomposition(3,[[4,[2,1,0,1,0]],[2,[1,0,0,1,0]]])
        [[12, [2, 1, 0, 1, 0]], [6, [1, 0, 0, 1, 0]]]


    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    """
    return [[n*i[0], i[1]] for i in dec]



def _twist_of_decompositions(dec1, dec2, a, n):
    r"""
    Exactly what is in __mul__ for bundles (i.e. tensor product), but taking and returning just
    the decompositions. 

    INPUT:
    
    The data corresponding to a pair of Bott bundles on a flag manifold:
    - ``dec1``, ``dec2`` - Lists whose elements are pairs of a positive integer (multiplicity) 
    and a list of nonnegative integers.
    - ``a`` - A list of nonnegative integers (the dimensions of the projective subspaces
    conforming the flags).
    - ``n`` - A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    
    OUTPUT:

    - A list of pairs as in the INPUT obtained by Littlewood--Richardson rule.
    
    EXAMPLES::
    
        sage: d1=[[1,[2,0,0,0]],[2,[0,0,1,0]]]
        sage: d2=[[1,[1,1,1,0]],[1,[1,0,1,1]]]
        sage: _twist_of_decompositions(d1,d2,[1],3)
        [[1, [3, 1, 1, 0]],
         [1, [3, 0, 1, 1]],
         [1, [2, 1, 1, 1]],
         [2, [1, 1, 1, 1]],
         [2, [1, 1, 2, 0]],
         [2, [1, 0, 2, 1]]]


    TODO:

    Eventually, when it is checked, this function will be called in the method __mul__
    instead of being written twice.

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
            result = _simplify_decomposition(result)
    return result





def _plethism_of_irred(lmbd, alpha, a, n, s):
    r"""
    Applies Schur functor to an irreducible Bott bundle. Mainly based on [M, Formula (8.9)].

    INPUT:
    
    - ``lmbd`` - a partition
    - ``alpha`` - A list of integers: the concatenated Schur functors on the quotient of adjacent 
    tautological subbundles of the flag manifold. The resulting bundle is the twist of all this 
    Schur functors of bundles
    - ``a`` - A list of nonnegative integers (the dimensions of the projective subspaces
    conforming the flags).
    - ``n`` - A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    - ``s`` - It should be SymmetricFunctions(QQ).schur()

    OUTPUT:

    - A list  whose elements are pairs of a positive integer (multiplicity) 
    and a list of nonnegative integers.
    
    EXAMPLES::
    
        sage: s=SymmetricFunctions(QQ).schur()
        sage: _plethism_of_irred([2,1], [1,0,0,0], [1], 3, s)
        [[1, [2, 1, 0, 0]]]
        sage: _plethism_of_irred([2,1], [0,0,1,0], [1], 3, s)
        [[1, [0, 0, 2, 1]]]
        sage: _plethism_of_irred([2,1], [1,0,1,0],[1],3,s)
        [[1, [2, 1, 3, 0]], [1, [2, 1, 2, 1]], [1, [3, 0, 2, 1]]]
        sage: _plethism_of_irred([4,0], [1,0,1,0], [1], 3, s)
        [[1, [4, 0, 4, 0]], [1, [3, 1, 3, 1]], [1, [2, 2, 2, 2]]]
        sage: _plethism_of_irred([4,0], [2,0,1,1], [1], 3, s)
        [[1, [4, 4, 4, 4]], [1, [6, 2, 4, 4]], [1, [8, 0, 4, 4]]]
    
    
    
    .. WARNING::

       INTERNAL FUNCTION! DO NOT USE DIRECTLY!

    REFERENCES:
    
    .. [M] I.G. Macdonald: "Symmetric Functions and Hall Polynomials", Oxford Mathematical Monographs, 
    (ISBN: 0198504500)
    """
    #s=SymmetricFunctions(QQ).schur()
    #Esto se basa en la f'ormula (8.9) del libro de Macdonald:
    #"Symmetric Functions and Hall Polynomials".
    if len(a) < 1:
        #Aqu'i tenemos el caso en que nuestro fibrado el el functor de 
        #Shur asociado a alpha de uno de los fibrados universales.
        #En este caso, con un pletismo tenemos suficiente.
        #print(s(lmbd),s(alpha))
        p = s(lmbd).plethysm(s(alpha))
        #print(p)
        #Una vez obtenido el pletismo, nos quedamos con los coeficientes
        #y las particiones:
        p = [p.coefficients(), [el[:] for el in p.support()]]
        #print(p)
        #El problema, es que algunas de estas particiones son m'as largas
        #que la dimensi'on del fibrado universal y otras sonm'as cortas
        #Quitamos las particiones largas (pues el functor de Schur 
        #saldr'ia cero) y completamos con ceros las cortas.
        l = 0
        while l < len(p[1]):
            #print(l,p)
            if len(p[1][l]) > len(alpha):
                p[1].pop(l)
                p[0].pop(l)
            else:
                p[1][l].extend([0 for i in range(len(alpha)-len(p[1][l]))])
                l += 1
        return [[p[0][l], p[1][l]] for l in range(len(p[1]))]
    #Si nuestro fibrado es un twist de dos o mas functores de Schur de fibrados
    #hay que hacerle pletismo al primero y recursi'on con lo dem'as
    Prt = [el[:] for el in Partitions(sum(lmbd))]
    result = []
    for mu in Prt:
        IntPr = s(lmbd).internal_product(s(mu))
        PletInic = IntPr.plethysm(s(alpha[:a[0]+1]))
        PletInic = [PletInic.coefficients(), [el[:] for el in PletInic.support()]]
        #Quitamos los largos y completamos con ceros
        l = 0
        while l < len(PletInic[1]):
            if len(PletInic[1][l]) > a[0]+1:
                PletInic[1].pop(l)
                PletInic[0].pop(l)
            else:
                PletInic[1][l].extend([0 for i in range(a[0]+1-len(PletInic[1][l]))])
                l += 1
        #Ya tenemos el primer factor del sumando (que tambi'en es una suma). 
        #Ahora toca el segundo, que es donde se usa la recursi'on.
#        print(mu,alpha[a[0]+1:],[a[i]-a[0] for i in range(1,len(a))],n-a[0]-1,PletInic)
        PletSegund = _plethism_of_irred(mu, alpha[a[0]+1:], [a[i]-a[0] for i in range(1, len(a))], n-a[0]-1, s)
        for l in range(len(PletInic[1])):
            for m in range(len(PletSegund)):
                result.append([PletInic[0][l]*PletSegund[m][0], PletInic[1][l][:]+PletSegund[m][1][:]])
    return _simplify_decomposition(result)
            


def _plethysmization(lmbd, dec, a, n, s):
    r"""
    Tool to compute the decomposition of the Schur functor of a Bundle. Based on [M, Formula (8.8)] 
    and function plethism_of_irred.
    
    INPUT:
    
    - ``lmbd`` - a partition
    - ``dec`` - A list whose elements are pairs of a positive integer (multiplicity) and a list of
    nonnegative integers.
    - ``a`` - A list of nonnegative integers (the dimensions of the projective subspaces
    conforming the flags).
    - ``n`` - A nonnegative integer (the dimension of the ambient projectiive space of the flags).
    - ``s`` - It should be SymmetricFunctions(QQ).schur()

    OUTPUT:

    - A list  whose elements are pairs of a positive integer (multiplicity) 
    and a list of nonnegative integers.
    
    EXAMPLES::
    
        sage: s=SymmetricFunctions(QQ).schur()
        sage: _plethysmization([2,1],[[1, [4, 0, 4, 0]], [1, [3, 1, 3, 1]], [1, [2, 2, 2, 2]]],[1],3,s)
        [[9, [7, 5, 6, 6]],
         [26, [7, 5, 8, 4]],
         [21, [7, 5, 9, 3]],
         [12, [7, 5, 10, 2]],
         [1, [7, 5, 12, 0]],
         [11, [8, 4, 6, 6]],
         [32, [8, 4, 8, 4]],
         [26, [8, 4, 9, 3]],
         [15, [8, 4, 10, 2]],
         [2, [8, 4, 12, 0]],
         [7, [9, 3, 6, 6]],
         [26, [9, 3, 8, 4]],
         [21, [9, 3, 9, 3]],
         [13, [9, 3, 10, 2]],
         [1, [9, 3, 12, 0]],
         [5, [10, 2, 6, 6]],
         [15, [10, 2, 8, 4]],
         [13, [10, 2, 9, 3]],
         [9, [10, 2, 10, 2]],
         [1, [10, 2, 12, 0]],
         [2, [11, 1, 6, 6]],
         [6, [11, 1, 8, 4]],
         [6, [11, 1, 9, 3]],
         [4, [11, 1, 10, 2]],
         [1, [11, 1, 12, 0]],
         [9, [6, 6, 7, 5]],
         [11, [6, 6, 8, 4]],
         [7, [6, 6, 9, 3]],
         [5, [6, 6, 10, 2]],
         [2, [6, 6, 11, 1]],
         [24, [7, 5, 7, 5]],
         [5, [7, 5, 11, 1]],
         [26, [8, 4, 7, 5]],
         [6, [8, 4, 11, 1]],
         [21, [9, 3, 7, 5]],
         [6, [9, 3, 11, 1]],
         [12, [10, 2, 7, 5]],
         [4, [10, 2, 11, 1]],
         [5, [11, 1, 7, 5]],
         [2, [11, 1, 11, 1]],
         [1, [12, 0, 7, 5]],
         [2, [12, 0, 8, 4]],
         [1, [12, 0, 9, 3]],
         [1, [12, 0, 10, 2]],
         [1, [12, 0, 11, 1]],
         [4, [6, 6, 6, 6]]]
        
    
    .. WARNING::
    
       INTERNAL FUNCTION! DO NOT USE DIRECTLY!
    
        It uses the class Partitions in 

        /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/partition.py

        Such class uses the class IntegerListsLex in

        /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/integer_list.py

        Which seems bugged
    
    
    REFERENCES:
    
    .. [M] I.G. Macdonald: "Symmetric Functions and Hall Polynomials", Oxford Mathematical Monographs, 
    (ISBN: 0198504500)
        
    """
    #print(dec)
    if len(dec) < 2:
        return _multiply_decomposition(dec[0][0], _plethism_of_irred(lmbd, dec[0][1], a, n, s))
    #s=SymmetricFunctions(QQ).schur()
    P = []
    #Primero cogemos todas las posibles particiones que son "menores" que
    #lambda, que son las que deberemos meter en el sumatorio.
    for i in range(sum(lmbd)+1):
        P = P+[el[:] for el in Partitions(i,outer=lmbd).list()]
    #A continuaci'on creamos los diagramas de Young "skew" esos y
    # los polinomios asociados a ellos. Los ponemos en forma de lista para que 
    # nos permitan hacer ciertas cosas
    sups = []
    for p in P:
        sk = s(lmbd).skew_by(s(p))
#        print [p,[el[:] for el in sk.support()],sk.coefficients()]
        sups.append([p,[el[:] for el in sk.support()], sk.coefficients()])
    #Ahora tenemos lo necesario para la f'ormula (8.8) del libro de Macdonald:
    #"Symmetric Functions and Hall Polynomials".
    #Reutilizamos las variables P y p para ahorrar memoria.
    #Aqu'i empieza la fiesta de la recursi'on.
    P = []
    #Me da la sensaci'on de que donde tenemos que buscar es en sups (estaba sk),
    #pues sk es solo el ultimo que hemos cogido
#    print 'sups=', sups
    for p in sups:
        auxdec = []
#        print 'p=', p
        #Primero hacemos el pletismo del primer sumando:  
        for i in range(len(p[1])):
            PlI = _plethism_of_irred(p[1][i], dec[0][1], a, n, s)
            auxdec.extend(_multiply_decomposition(p[2][i], PlI))
        auxdec = _simplify_decomposition(auxdec)
        #Ya tenemo el pletismo sobre el primer sumando. Como este sumando
        #puede estar repetido, quitarlo no es inmediato. Hay que
        #quitarlo una sola vez, porque esto habr'a que twistarlo con todo,
        #incluidas sus propias copias:
        if dec[0][0] == 1:
            remdec = dec[1:]
        else:
            remdec = [[dec[0,0]-1, dec[0,1]]] + dec[1:]
        #Now the recursive command, remdec is the remaining part of dec
        #after removing the first summand:
        P.extend(_twist_of_decompositions(auxdec, _plethysmization(p[0], remdec, a, n, s), a, n))
    P=_simplify_decomposition(P)
    return P



###################################################################
#Now the class we use for reducible bundles. 
#It intends to be defined for arbitrary flag manifolds.
###################################################################


class BottBundle():
    r"""
    This class represents homogeneous vector bundles on a given homogeneous space (flag manifold) 
    that can be decomposed as a sum of irreducible homogeneous vector bundles whose cohomology 
    can be computed by means of Bott's Theorem. 
    
    Such irreducible bundles are twists of Schur functors on the (duals of) quotients of 
    tautological subbundles that are adjacent in the sequence (see [W] and [FH]):

    $0\subset R_{a_0+1}\subset...\subset R_{a_l+1}\subset 
    \mathbb{C}^{n+1}\times\mathrm{Flag}(a_0,...,a_l;n)$
    
    Attributes:

    - ``space`` - The ``FlagManifold`` on which the bundle is defined

    - ``decomposition`` - A list whose elements are pairs of a positive integer (multiplicity) and a list of nonnegative integers (irreducible ``BottBundle``).

    - ``aspect`` - A simpler way to represent the bundle that keeps the way it was constructed.

    Methods:
    
    - Addition is equivalent to direct sum of vector bundles.

    - Multiplication is tensor product.

    - ``cohomology`` - Computes the whole sheaf cohomology of ``self``.

    - ``rank`` - returns the rank of ``self``.

    - ``schur`` - performs a given partition's schur functor of ``self``. 

    - ``sym`` - a given integer symmetric power of ``self``.

    - ``wedge`` - a given number exterior power of ``self``.
    
    INPUT:

    - ``space`` - a value of type FlagManifold.

    OUTPUT:

    - The zero bundle on ``space``.

    EXAMPLES::
    
        sage: F=FlagManifold([1,3,5],7)
        sage: BottBundle(F)
        Homogeneous rank-0 vector bundle on the space of flags of subspaces of dimensions [1, 3, 5] in the 
        7-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
    
    
    WARNING:
    
    The methods ``schur``, ``sym`` and ``wedge`` depend on a bugged package. 
    
    TODO:
    
    - Method __repr__ should be different when self.space is a grassmannian (it should speak about universal bundles) or a projective space (about tangent bundle and hyperplane section).

    - A function to compute a concrete cohomology group.

    - The attribute ``aspect`` has yet to be adequated to ``FlagManifold``s that are not ``Grassmannian``.

    - A function that computes the dual of a ``BottBundle``.  
    
    REFERENCES:
    
    .. [FH] William Fulton, Joe Harris: "Representation Theory - A First Course", GTM, Springer (ISBN 978-1-4612-0979-9), 2004

    .. [M] I.G. Macdonald: "Symmetric Functions and Hall Polynomials", Oxford Mathematical Monographs, (ISBN: 0198504500), 2nd ed., 1999

    .. [W]  Jerzy Weyman: "Cohomology of Vector Bundles and Syzygies", Cambridge Tracts in Mathematics (ISBN: 9780521621977), 2003
        

    """

    def __init__(self, space):
        r"""
        Initialize ``self``
        
        INPUT:
        - ``space`` - a value of type FlagManifold.
        
        EXAMPLES::
        
            sage: F=FlagManifold([1,3,5],7)
            sage: BottBundle(F)
            Homogeneous rank-0 vector bundle on the space of flags of subspaces of dimensions [1, 3, 5] in the 
            7-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
        
        """
        self.decomposition = []
        self.space = space
        self.aspect = '0'





    def __repr__(self):
        r"""
        Represent ``self``

        TODO::
        
            Special names when ``self.space`` is a grassmannian or a projective space.
        """
        cadena = 'Homogeneous rank-' + str(self.rank()) + ' vector bundle on ' + self.space.__repr__()+'.\n' + 'It is the sum of the following irreducible homogeneous vector bundles:' + '\n'
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

        - ``self``, ``other`` - of BottBundle class, with the same attribute ``space``.

        OUTPUT:
        
        - The direct sum of ``self`` and ``other``.

        EXAMPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: F=G.S()
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
            print 'Error: Cannot sum bundles on different varieties'
            return
        result = BottBundle(self.space)
        result.aspect = self.aspect + '(+)' + other.aspect
        result.decomposition = _simplify_decomposition(deepcopy(self.decomposition) + deepcopy(other.decomposition))
#        for a in other.decomposition:
#            notyet=true
#            for b in result.decomposition:
#                if b[1]==a[1]:
#                    a[0]=a[0]+b[0]
#                    notyet=false
#                    break
#            if notyet:
#                result.decomposition.append(a)
        return result

        


    def __mul__(self, other):
        r"""
        Twists ``self`` times ``other``.

        INPUT:

        - ``self``, ``other`` - of BottBundle class, with the same attribute ``space``.

        OUTPUT:

        - The tensor product of ``self`` and ``other``.
  
        EXAMPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: F=G.S()
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
            print 'Error: Cannot twist bundles on different varieties'
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
                result.decomposition = _simplify_decomposition(result.decomposition)
        return result    
        


    def schur(self, lmbd):
        r"""
        A given partition's schur functor of ``self``.

        INPUT:

        - ``lmbd`` - a partition (decreasin list of possitive integers)

        OUTPUT:
    
        - the Schur functor associated to partition ``lmbd`` applied to ``self``. 
        

        EXAMPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: (E*E).schur([2,1])
            Homogeneous rank-20 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>            
            3 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            2 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [5, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        .. WARNING::
        
        
            It uses the class Partitions in 
    
            /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/partition.py
    
            Such class uses the class IntegerListsLex in
    
            /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/integer_list.py
    
            Which seems bugged.

            Plethism computation is slow, so difficult Schur functors of complex bundles can take forever to compute.
        """
        result = BottBundle(self.space)
        result.aspect = 'Schur_' + lmbd.__str__() + '('+self.aspect+')'
        s = SymmetricFunctions(QQ).schur()
        # Dado que los algoritmos son recursivos, es mejor pasar s como
        # par'ametro para no crear tantas varables iguales...
        result.decomposition=_plethysmization(lmbd, self.decomposition, self.space.a, self.space.n, s)
        return result
            



    def wedge(self, k):
        r"""
        Computes the ``k``-th exterior power of ``self``
        
        INPUT:

        - ``k`` A possitive integer

        OUTPUT:

        - the ``k``-th exterior power of ``self``

        EXAMPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: (E*E).wedge(3)
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>


        .. WARNING::
        
        
            It uses the class Partitions in 
    
            sage.combinat.partition
    
            Such class uses the class IntegerListsLex in
    
            sage.combinat.integer_list
    
            Which seems bugged.

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

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: (E*E).sym(3)
            Homogeneous rank-20 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            2 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [3, 3] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            2 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [4, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [5, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [6, 6] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [6, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        .. WARNING::
        
        
            It uses the class Partitions in 
    
            /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/partition.py
    
            Such class uses the class IntegerListsLex in
    
            /sage-6.5-x86_64-Linux/local/lib/python2.7/site-packages/sage/combinat/integer_list.py
    
            This last class seems bugged.

            Plethism computation is slow, so difficult Schur functors of complex bundles can take forever to compute.
        """
        return self.schur([k])


    def rank(self):
        r"""
        Computes the rank of the vector bundle.

        EXAMPLES:
    
            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
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
                PartialRank *= _schur_dimension(d[1][previousa+1:a+1])
                previousa = a
            PartialRank *= _schur_dimension(d[1][previousa+1:self.space.n+1])
            TotalRank += d[0] * PartialRank
        return TotalRank
                
        


    def cohomology(self):
        r"""
        Computes the complete cohomology of the vector bundle.

        OUTPUT:

        A list ``l`` of integers such that ``l[i]`` is the diension of the ``i``-th cohomology space
        of the bundle. The length of ``l`` is equal to the index of the last nonvanishing cohomology
        space (+1). Therefore, if the length is zero, all cohomology vanishes.
        
        EXAMPLES::
    
            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: E=G.Q()
            sage: E.cohomology()
            h^ 0  =  4
            All remaining cohomology is zero
            [4]
            sage: F=G.O(-7)
            sage: (E*F).cohomology()
            h^ 4  =  60
            All remaining cohomology is zero
            [0, 0, 0, 0, 60]
            sage: (E*F+E).cohomology()
            h^ 0  =  4
            h^ 4  =  60
            All remaining cohomology is zero
            [4, 0, 0, 0, 60]

        TODO:

        - a possible second argument (or a separate function to return a concrete cohomology space).

        """
        n = len(self.decomposition[0][1])
        B = []
        for el in self.decomposition:
            lista = deepcopy(el[1])
            aux = _bott(lista)
            aux.append(el[0])            
            B.append(aux)
        B.sort()
        m = B[len(B)-1][0]
        Cohom = []
        for ind in range(0, m+1):
            Cohom.append(0)
        for el in B:
            if el[0] <> -1:
                Cohom[el[0]] += el[2] * _schur_dimension(el[1])
        for i in range(0, m+1):
            if Cohom[i] <> 0:
                print "h^", i, " = ", Cohom[i]
        if len(Cohom) == 0:
            print 'All cohomology vanishes'
        else:
            print 'All remaining cohomology is zero'
        return Cohom




class _UniversalSubbundleDual(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the dual to the universal subbundle of an instance of
    ``FlagManifold``

    INPUT:

    - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

    OUTPUT:

    - The dual to the universal subbundle of ``fm``.

    EXAPLES::

        sage: G=Grassmannian(1,3)
        You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^4 --> Q --> 0
        where rk(S)=2 and rk(Q)=2.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: _UniversalSubbundleDual(G)
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: _UniversalSubbundleDual(G)
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
        """
        self.decomposition = [[1,[1]+[0]*(fm.n)]]
        self.space = fm
        self.aspect = 'S_1^*'


class _UniversalSubbundle(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the universal subbundle of an instance of ``Grassmannian``

    INPUT:

    - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

    OUTPUT:

    - The universal subbundle of ``fm``.

    EXAPLES::

        sage: G=Grassmannian(1,3)
        You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^4 --> Q --> 0
        where rk(S)=2 and rk(Q)=2.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: _UniversalSubbundle(G)                                      
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: _UniversalSubbundle(G)                                      
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
        """
        self.decomposition = [[1, [1]*(fm.a[0])+[0]+[1]*(fm.n-fm.a[0])]]
        self.space = fm
        self.aspect = 'S_1'

class _UniversalQuotientDual(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the dual to the universal quotient bundle of an instance of ``Grassmannian``

    INPUT:

    - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

    OUTPUT:

    - The dual to the universal quotient bundle of ``fm``.

    EXAPLES::

        sage: G=Grassmannian(1,3)
        You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^4 --> Q --> 0
        where rk(S)=2 and rk(Q)=2.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: _UniversalQuotientDual(G)
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: _UniversalQuotientDual(G)
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        THOUGHT FOR GRASSMANNIANS!

        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
        """
        self.decomposition = [[1,[0]*(fm.a[0]+1)+[1]+[0]*(fm.n-fm.a[0]-1)]]
        self.space = fm
        self.aspect = 'Q_1^*'

class _UniversalQuotient(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize the universal quotient bundle of an instance of
    ``FlagManifold``

    INPUT:

    - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

    OUTPUT:

    - The universal quotient bundle of ``fm``.

    EXAPLES::

        sage: G=Grassmannian(1,3)
        You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^4 --> Q --> 0
        where rk(S)=2 and rk(Q)=2.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: _UniversalQuotient(G)    
        Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>


    WARNING:

    THOUGHT FOR GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
    """
    def __init__(self, fm):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` - an instance of ``FlagManifold``, it should be a ``Grassmannian``

        OUTPUT:

        - The universal quotient bundle of ``fm``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: _UniversalQuotient(G)    
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
        """
        self.decomposition = [[1, [1]*(fm.n)+[0]]]
        self.space = fm
        self.aspect = 'Q_1'

class _LineBundleGrass(BottBundle):
    r"""
    A child class of ``BottBundle`` to easily initialize line bundles on grassmannians

    INPUT:

    - ``grass`` - an instance of ``FlagManifold`` it should be a ``Grassmannian``

    - ``k`` - an integer.

    OUTPUT:

    - the ``k``-th multiple of the hyperplane section of ``grass``

    EXAPLES::

        sage: G=Grassmannian(1,3)
        You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^4 --> Q --> 0
        where rk(S)=2 and rk(Q)=2.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: _LineBundleGrass(G,3)
        Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>
        sage: _LineBundleGrass(G,-2)
        Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
        It is the sum of the following irreducible homogeneous vector bundles:
        <BLANKLINE>
        1 time(s) the tensor product of: 
            Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
            Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
        <BLANKLINE>

    WARNING:

    WILL NOT WORK PROPERLY WHEN ``grass`` IS NOT A GRASSMANNIANS!

    BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
    """
    def __init__(self, grass, k):
        r"""
        Initialize ``self``

        INPUT:

        - ``fm`` - an instance of ``FlagManifold``

        OUTPUT:

        - The dual to the universal subbundle of ``fm``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: _LineBundleGrass(G,3)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            sage: _LineBundleGrass(G,-2)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        WARNING:
    
        BETTER CONSIDERED AN INTERNAL FUNCTION ONLY. DIRECT USAGE IS INADVISABLE
        """
        if k > 0:
            self.decomposition = [[1, [k]*(grass.a[0]+1)+[0]*(grass.n-grass.a[0])]]
        elif k < 0:
            self.decomposition = [[1, [0]*(grass.a[0]+1)+[-k]*(grass.n-grass.a[0])]]
        else:
            self.decomposition = [[1, [0]*(grass.n+1)]]
        self.space = grass
        self.aspect = 'O(' + str(k) + ')'






####################################################################
#Now we introduce the classes referring the homogeneous spaces.
#Since we just use the special linear group, all are flag manifolds
####################################################################

class FlagManifold():
    r"""
    This is the class containing all quotients of the general linear group by
    a parabolic subgroup. 
        
    Attributes:

    - ``a`` - a list of increasing integers (dimensions of the subspaces of the flag)

    - ``n`` - possitive integer (the dimension of the ambient projective space of the flag)
    
    
    INPUT:
    
    - ``a`` - a list of increasing integers (dimensions of the subspaces of the flag)
    
    - ``n`` - possitive integer (the dimension of the ambient projective space of the flag)
    
    OUTPUT:
    
    - The space of flags of subspaces of dimensions given by ``a`` in the projective space of dimension ``n``.
    
    EXAMPLES::
    
        sage: FlagManifold([1,3,5],7)
        the space of flags of subspaces of dimensions [1, 3, 5] in the 7-dimensional projective space    
    
    WARNING:
    
    It only contains properties related with the class ``BottBundle``, so no more geometry is 
    allowed with this class.
    
    TODO:
    
    - Methods to generate special instances of the class ``BottBundle``, mainly irreducible ones and their duals
        
    """
    def __init__(self, a, n):
        r"""
        Initializes ``self``.

        INPUT:
    
        - ``a`` - a list of increasing integers (dimensions of the subspaces of the flag)
    
        - ``n`` - possitive integer (the dimension of the ambient projective space of the flag)
    
        OUTPUT:
    
        - The space of flags of subspaces of dimensions given by ``a`` in the projective space of 
        dimension ``n``.
    
        EXAMPLES::
        
            sage: FlagManifold([1,3,5],7)
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

    def __repr__(self):
        r"""
        Represents ``self``.
    
        EXAMPLES::
        
            sage: FlagManifold([1,3,5],7)
            the space of flags of subspaces of dimensions [1, 3, 5] in the 7-dimensional projective space
    
        """
        return 'the space of flags of subspaces of dimensions ' + str(self.a) + ' in the ' + str(self.n) + '-dimensional projective space'

    def is_grassmannian(self):
        r"""
        Determines wether ``self`` is a grassmannian.
    
        EXAMPLES::
        
            sage: FlagManifold([1,3,5],7).is_grassmannian()                          
            False
            sage: FlagManifold([5],7).is_grassmannian()    
            True
            sage: FlagManifold([0],7).is_grassmannian()
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
        
            sage: FlagManifold([1,3,5],7).is_projective_space()
            False
            sage: FlagManifold([5],7).is_projective_space()    
            False
            sage: FlagManifold([0],7).is_projective_space()
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


class Grassmannian(FlagManifold):
    r"""
    A child class of ``FlagManifold`` to easily initialize grassmannians

    Methods:

    - ``O`` - creates a multiple of the hyperplane section

    - ``Q`` - creates the universal quotient bundle

    - ``Q_dual`` - creates the dual of the universal quotient bundle

    - ``S`` - creates the universal subbundle

    - ``S_dual`` - creates the dual of the universal subbundle

    - ``T`` - Creates the tangent bundle

    - ``Om`` - Creates the cotangent bundle


    INPUT:

    - ``k`` - a list of increasing integers (dimension of the projective subspace)

    - ``n`` - possitive integer (the dimension of the ambient projective space of the flag)

    OUTPUT:

    - The grassmannian of ``k``-dimensional projective subspaces in the projective space of dimension ``n``.

    EXAMPLES::

        sage: G1=Grassmannian(4,9)                                               
        You have defined the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^10 --> Q --> 0
        where rk(S)=5 and rk(Q)=5.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.        
        sage: G2=Grassmannian(2,10)
        You have defined the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space
        Consider the uiversal short exact sequence:
        0 --> S --> O^11 --> Q --> 0
        where rk(S)=3 and rk(Q)=8.
        These bundles can be defined with methods S() and Q().
        Their duals can be obtained with methods SDual() and Qdual()
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        sage: G1,G2
        (the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space,
         the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space)

    WARNING:

    It only contains properties related with the class ``BottBundle``, so no more geometry is 
    allowed with this class.

    """
    def __init__(self,k,n):
        r"""
        Initializes ``self``

        INPUT:
    
        - ``k`` - a list of increasing integers (dimension of the projective subspace)
    
        - ``n`` - possitive integer (the dimension of the ambient projective space of the flag)
    
        OUTPUT:
    
        - The grassmannian of ``k``-dimensional projective subspaces in the projective space of 
        dimension ``n``.
    
        EXAMPLES::
    
            sage: G1=Grassmannian(4,9)                                               
            You have defined the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^10 --> Q --> 0
            where rk(S)=5 and rk(Q)=5.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.        
            sage: G2=Grassmannian(2,10)
            You have defined the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^11 --> Q --> 0
            where rk(S)=3 and rk(Q)=8.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
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
        print 'You have defined the grassmannian of ' + str(self.a[0]) + '-dimensional subspaces of the ' + str(self.n) + '-dimensional projective space'
        print 'Consider the uiversal short exact sequence:'
        print '0 --> S --> O^' + str(self.n+1) + ' --> Q --> 0'
        print 'where rk(S)=' + str(self.a[0]+1) + ' and rk(Q)=' + str(self.n-self.a[0]) + '.'
        print 'These bundles can be defined with methods S() and Q().'
        print 'Their duals can be obtained with methods SDual() and Qdual()'
        print 'Method O(k) gives the k power of the hyperplane section.'
        print 'Method Om() gives the cotangent bundle.'
        print 'Method T() gives the tangent bundle.'


    def __repr__(self):
        r"""
        Represents ``self``
    
        EXAMPLES::
    
            sage: G1=Grassmannian(4,9)                                               
            You have defined the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^10 --> Q --> 0
            where rk(S)=5 and rk(Q)=5.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.        
            sage: G2=Grassmannian(2,10)
            You have defined the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^11 --> Q --> 0
            where rk(S)=3 and rk(Q)=8.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G1,G2
            (the grassmannian of 4-dimensional subspaces of the 9-dimensional projective space,
             the grassmannian of 2-dimensional subspaces of the 10-dimensional projective space)
    
        """
        return 'the grassmannian of ' + str(self.a[0]) + '-dimensional subspaces of the ' + str(self.n) + '-dimensional projective space'

    def O(self, k=0):
        r"""
        k-th power of the hyperplane section of the grassmannian

        INPUT:

        - ``k`` - an integer.

        OUTPUT:

        - the ``k``-th multiple of the hyperplane section of ``self``

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.O(3)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>    
            1 time(s) the tensor product of: 
                Schur functor of partition [3, 3] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
            sage: G.O(-2)
            Homogeneous rank-1 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 2] of the dual of the last quotient of tautological subbundles,
            sage: G.O()
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
        return _LineBundleGrass(self, k)

    def Q(self):
        r"""
        Universal quotient bundle of the grassmannian

        OUTPUT:

        - The universal quotient bundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.Q()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return _UniversalQuotient(self)

    def S(self):
        r"""
        Universal subbundle of the grassmannian

        OUTPUT:

        - The universal subbundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.S()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>

        """
        return _UniversalSubbundle(self)

    def Q_dual(self):
        r"""
        Dual of the universal quotient bundle of the grassmannian

        OUTPUT:

        - The dual to the universal quotient bundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.Q_dual()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [0, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return _UniversalQuotientDual(self)

    def S_dual(self):
        r"""
        Dual of the universal subbundle of the grassmannian

        OUTPUT:

        - The dual to the universal subbundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.S_dual()
            Homogeneous rank-2 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [0, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return _UniversalSubbundleDual(self)

    def Om(self):
        r"""
        Cotangent bundle of the grassmannian

        OUTPUT:

        - The cotangent bundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.Om()
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [1, 0] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [2, 1] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return self.S()*self.Q_dual()

    def T(self):
        r"""
        Tangent bundle of the grassmannian

        OUTPUT:

        - The tangent bundle of ``self``.

        EXAPLES::

            sage: G=Grassmannian(1,3)
            You have defined the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space
            Consider the uiversal short exact sequence:
            0 --> S --> O^4 --> Q --> 0
            where rk(S)=2 and rk(Q)=2.
            These bundles can be defined with methods S() and Q().
            Their duals can be obtained with methods SDual() and Qdual()
            Method O(k) gives the k power of the hyperplane section.
            Method Om() gives the cotangent bundle.
            Method T() gives the tangent bundle.
            sage: G.T()
            Homogeneous rank-4 vector bundle on the grassmannian of 1-dimensional subspaces of the 3-dimensional projective space.
            It is the sum of the following irreducible homogeneous vector bundles:
            <BLANKLINE>
            1 time(s) the tensor product of: 
                Schur functor of partition [2, 1] of the dual of the 1st/nd/rd/th  quotient of tautological subbundles,
                Schur functor of partition [1, 0] of the dual of the last quotient of tautological subbundles,
            <BLANKLINE>
        """
        return self.S_dual()*self.Q()





class Proj(Grassmannian):
    r"""
    A child class of ``Grassmannian`` to initialize easily projective spaces.

    INPUT:

    - ``n`` - an integer.

    OUTPUT:

    - The ``n``-dimensional projective space.

    EXAMPLES::

        sage: Proj(3)
        You have defined the projective space of dimension 3.
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        the projective space of dimension 3
        sage: Proj(7)
        You have defined the projective space of dimension 7.
        Method O(k) gives the k power of the hyperplane section.
        Method Om() gives the cotangent bundle.
        Method T() gives the tangent bundle.
        the projective space of dimension 7

    WARNING:

    It only contains properties related with the class ``BottBundle``, so no more geometry is 
    allowed with this class.

    """
    def __init__(self, n):
        if type(n)<>Integer or n < 0:
            raise ValueError("Argument must be a nonnegative integer")
        self.a = [0]
        self.n = n
        print 'You have defined the projective space of dimension ' + str(self.n) + '.'
        print 'Method O(k) gives the k power of the hyperplane section.'
        print 'Method Om() gives the cotangent bundle.'
        print 'Method T() gives the tangent bundle.'

    def __repr__(self):
        return 'the projective space of dimension ' + str(self.n)




#      A          QQQQ     U     U  I
#     A A        Q    Q    U     U  I
#    A   A      Q      Q   U     U  I
#   A     A     Q    Q Q   U     U  I
#  A A A A A     Q    QQ    U   U   I
# A         A     QQQQ  Q    UUU    I
#A

#The following line were intended to be a tutorial for those loading the package.
#Since they appear when starting sage, probably they are better off.
#
#print 'You have loaded the package of homogeneous bundles cohomology.'
#print 'Assign to a variable either:'
#print '    - Proj(n) to create the n-dimensional projective space, or'
#print '    - Grassmannian(k,n) to create the grassmannian of k-planes in P^n, or'
#print '    - FlagManifold([a_0,...,a_n],n) to create a flag manifold.'
#print 'There is a lot of work to do yet on general flag manifolds...'
