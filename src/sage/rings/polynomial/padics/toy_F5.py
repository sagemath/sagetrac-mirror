"""
This module implements Groebner bases computation via the F5 algorithm.

AUTHOR:

- Tristan Vaccon (2016--2017)
"""


#*****************************************************************************
#  Copyright (C) 2017 Tristan Vaccon <tristan.vaccon@unilim.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from copy import copy
from sage.structure.sequence import Sequence
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix, Matrix

def list_monom(P,d):
    r"""

    Return a decreasingly ordered list of the monomials of P of degree d
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_F5 import list_monom
        sage: list_monom(S,2)
        [x^2, x*y, y^2, x*z, y*z, z^2]
    """
    if d == 0:
        return [P(1)]
    l1 = [ a for a in P.gens()]
    if d == 1:
        return l1
    ld1 = list_monom(P,d-1)
    ld = [[x*u for x in l1] for u in ld1]
    ld = sum(ld,[])
    s = []
    for i in ld:
        if i not in s:
            s.append(i)
    s.sort()
    s.reverse()
    return s

def list_monom_iteratif(P,d,ld1):
    r"""

    Return a decreasingly ordered list of the monomials of P of degree d, built from the corresponding list of degree d-1
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer
    - ``ld1`` -- the list of the monomials of degree d-1

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_F5 import list_monom_iteratif
        sage: list_monom_iteratif(S,2,[x,y,z])
        [x^2, x*y, y^2, x*z, y*z, z^2]
    """
    if d == 0:
        return [P(1)]
    l1 = [ a for a in P.gens()]
    if d == 1:
        return l1
    ld = [[x*u for x in l1] for u in ld1]
    ld = sum(ld,[])
    s = list(set(ld))
    s.sort()
    s.reverse()
    return s

def index_last_variable(P,mon):
    r"""

    Return the index of the last variable dividing mon
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``mon`` -- a monomial

    OUTPUT:

    Return the index of the last variable dividing mon, or -1 for constants

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_F5 import index_last_variable
        sage: index_last_variable(S,x*y)
        1
        sage: index_last_variable(S,S(1))
        -1
    """    
    list_expo = P(mon).exponents()[0]
    l = len(list_expo)
    i = l-1
    while (i >=0):
        if list_expo[i]>0:
            return i
        i=i-1
    return i


def list_monom_iteratif_no_repeat(P,d,ld1):
    r"""

    Return a decreasingly ordered list of the monomials of P of degree d, built from the corresponding list in degree d-1
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer
    - ``ld1`` -- the decreasingly ordered list of the monomials of degree d-1

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_F5 import list_monom_iteratif_no_repeat
        sage: list_monom_iteratif_no_repeat(S,2,[x,y,z])
        [x^2, x*y, y^2, x*z, y*z, z^2]
    """
    if d == 0:
        return [P(1)]
    l1 = [ a for a in P.gens()]
    n = len(l1)
    if d == 1:
        return l1
    ld = []
    for mon in ld1:
        i = index_last_variable(P,mon)
        for j in range(i,n):
            ld.append(mon*l1[j])
    ld.sort()
    ld.reverse()
    return ld

def Make_Macaulay_Matrix(list_poly, monom):
    r"""

    Build the Macaulay matrix for the polynomials of degree d in list_poly. It is represented as a couple of a matrix and the list of the signatures of its rows, assuming compatibility.
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix
        sage: list_poly = [[0,1,x**2+y**2], [1,1,x*y], [2,1,y*z], [3,1,z**2]]
        sage: Mac = Make_Macaulay_Matrix(list_poly,[x**2, x*y, y**2, x*z, y*z, z**2])
        sage: Mac[0]
        [1 0 1 0 0 0]
        [0 1 0 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: Mac[1]
        [[0, 1], [1, 1], [2, 1], [3, 1]]
        
    """

    R = list_poly[0][2].parent()
    d = list_poly[0][2].degree()
    l = len(monom)

    listsign = []
    Mac = Matrix(R.base_ring(), 0, l, [])
    Mac = [Mac, []]
    for fuple in list_poly:
        listsign.append([fuple[0],fuple[1]])
        f = fuple[2]
        list = f.monomials()
        flist = [ R(0) ] * l
        for mon in list:
            if mon in monom:
                i = monom.index(mon)
                flist[i] = f.monomial_coefficient(mon)
        fmat = Matrix(R.base_ring(), 1, l, flist)
        Mac[0] = Mac[0].stack(fmat)
    Mac[1] = listsign
    return Mac

def Make_Macaulay_Matrix_sparse_dict(list_poly, monom):
    r"""

    Build the Macaulay matrix for the polynomials of degree d in list_poly. It is represented as a couple of a matrix and the list of the signatures of its rows, assuming compatibility.
    The matrix is sparse.
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix_sparse_dict
        sage: list_poly = [[0,1,x**2+y**2], [1,1,x*y], [2,1,y*z], [3,1,z**2]]
        sage: Mac = Make_Macaulay_Matrix_sparse_dict(list_poly,[x**2, x*y, y**2, x*z, y*z, z**2])
        sage: Mac[0]
        [1 0 1 0 0 0]
        [0 1 0 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: Mac[1]
        [[0, 1], [1, 1], [2, 1], [3, 1]]
        
    """

    R = list_poly[0][2].parent()
    d = list_poly[0][2].degree()
    l = len(monom)

    dict_monom = {monom[aa]:aa for aa in range(l)}


    nrows = len(list_poly)

    listsign = []
    Mac = MatrixSpace(R.base_ring(),nrows,l,sparse=True)(0)
    for u in range(nrows):
        fuple=list_poly[u]
        listsign.append([fuple[0],fuple[1]])
        f = fuple[2]
        list = f.monomials()
        for mon in list:
            if dict_monom.has_key(mon):
                j = dict_monom[mon]
                Mac[u,j] = f.monomial_coefficient(mon)
    Mac = [Mac, listsign]
    return Mac

def F5_criterion(sign,G):
    r"""
    
    Inside the F5 algorithm
    Check the F5 criterion for the signature sign = (i, x**alpha) and G
    an S-GB up to degree deg(f_i)+deg(x**alpha) and signature strictly less than (i,1)
 
    INPUT:

    - ``sign`` -- a (potential) signature of a polynomial
    - ``G`` -- a list of polynomials with signature
    
    OUTPUT:

    a boolean, testing whether the F5 criterion for sign and G is passed

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import F5_criterion
        sage: sign = [1,x**3]
        sage: G =  [[0,1,x**2], [1,1,x*y], [2,1,y*z], [3,1,z**2]]
        sage: F5_criterion(sign,G)
        True
        sage: sign2 = [0,x**2]
        sage: sign3 = [3,y**2]
        sage: F5_criterion(sign2,G)
        False
        sage: F5_criterion(sign3,G)
        False
        
    """
    for g in G:
        ei = g[0]
        glm = g[2].lm()
        if ei < sign[0]:
            if glm.divides(sign[1]):
                return True
    return False  

def Reducible_resp_sign(f,g):
    r"""
    
    Inside the F5 algorithm, f and g are of the form [i,x**alpha, poly].
    Tests whether g can reduce f while respecting their signatures with strict inequality
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``g`` -- a polynomial with signature
    
    OUTPUT:

    a boolean, testing whether g can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Reducible_resp_sign
        sage: g = [0,x,x**2]
        sage: f = [1,x**2,x**3*y]
        sage: f2 = [0,x**2,x**3]
        sage: f3 = [0,1,x**3]
        sage: Reducible_resp_sign(f,g)
        True
        sage: Reducible_resp_sign(f2,g)
        False
        sage: Reducible_resp_sign(f3,g)
        False
    """
    R = f[2].parent()
    flm = f[2].lm()
    glm = g[2].lm()
    if glm.divides(flm):
        return (f[0]> g[0] or (f[0] == g[0] and R.monomial_quotient(flm,glm)*g[1]<f[1]))
    return False


def Reducible_large_resp_sign(f,g):
    r"""
    
    Inside the F5 algorithm, f and g are of the form [i,x**alpha, poly].
    Tests whether g can reduce f while respecting their signatures with non-strict inequality
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``g`` -- a polynomial with signature
    
    OUTPUT:

    a boolean, testing whether g can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Reducible_large_resp_sign
        sage: g = [0,x,x**2]
        sage: f = [1,x**2,x**3*y]
        sage: f2 = [0,x**2,x**3]
        sage: f3 = [0,1,x**3]
        sage: Reducible_large_resp_sign(f,g)
        True
        sage: Reducible_large_resp_sign(f2,g)
        True
        sage: Reducible_large_resp_sign(f3,g)
        False
    """
    R = f[2].parent()
    flm = f[2].lm()
    glm = g[2].lm()
    if glm.divides(flm):
        return (f[0]> g[0] or (f[0] == g[0] and R.monomial_quotient(flm,glm)*g[1]<=f[1]))
    return False

def Reducible_large_resp_sign_family(f,G):
    r"""
    
    Inside the F5 algorithm, f and the elements of the list G are of the form [i,x**alpha, poly].
    Tests whether G can reduce f while respecting their signatures
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``G`` -- a list of polynomials with signature
    
    OUTPUT:

    a boolean, testing whether G can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: G = [[0,1,x**2], [1,1,x*y], [2,1,y*z], [3,1,z**2]]
        sage: f = [1,x**2,x**3*y]
        sage: f2 = [0,x**2,x**3]
        sage: f3 = [0,1,x**3]
        sage: Reducible_large_resp_sign_family(f,G)
        True
        sage: Reducible_large_resp_sign_family(f2,G)
        True
        sage: Reducible_large_resp_sign_family(f3,G)
        False
    """
    for g in G:
        if Reducible_large_resp_sign(f,g):
            return True
    return False


def Rewritten(R,G,g):
    r"""
    
    Find the latest element of G such it has a multiple with same signature as g.

    
    INPUT:

    - ``R`` -- a polynomial ring
    - ``G`` -- a list of polynomials with signature, ordered by signature
    - ``g`` -- a polynomial with signature
    - ``paires`` -- a list of pairs of polynomials with signature

    
    OUTPUT:

    A multiple of an element of G with same signature and leading monomial
    as g. It is obtained with the latest element of G possible

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: G = [[0,y,x**2+2*z**2], [1,x,x*y+y**2], [2,1,y*z+z**2], [3,1,z**2]]
        sage: g = [2,x,2*x*y*z+x*z**2+y**2*z]
        sage: Rewritten(S,G,g)
        [2, x, x*y*z + x*z^2]

    """
    mon_sign = R(g[1]).lm()
    newg=g
    for gg in G:
        if gg[0]==g[0] and R(gg[1]).lm().divides(mon_sign):
            glm = R(gg[1]).lm()
            mult = R.monomial_quotient(mon_sign,glm)
            newg = [gg[0],mult*gg[1],mult*gg[2]]
    return newg


def SymbolicPreprocessing(R,G, d, paires, monom):
    r"""
    
    Constructs the Macaulay matrix from the pairs in paires,
    while using the F5 criterion

    
    INPUT:

    - ``R`` -- a polynomial ring
    - ``G`` -- a list of polynomials with signature
    - ``d`` -- an integer, representing the degree of the processed polynomials
    - ``paires`` -- a list of pairs of polynomials with signature
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    
    OUTPUT:

    A Macaulay matrix, as a couple of a matrix built from polynomials and
    a list of signatures for the rows of this matrix 

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: G = [[0,y,x**2+2*z**2], [1,x,x*y+y**2], [2,1,y*z+z**2], [3,1,z**2]]
        sage: paires = []
        sage: paires.append([3,1,x,[1,x,x*y+y**2],y, [0,1,x**2+2*z**2]])
        sage: paires.append([3,3,y,[3,1,x*z+3*y*z+z**2],z, [3,1,x*y+2*y*z+z**2]])
        sage: d = 3
        sage: Mac = SymbolicPreprocessing(S, G, d, paires,[x**3, x**2*y, x*y**2, y**3, x**2*z, x*y*z, y**2*z, x*z**2, y*z**2, z**3])
        sage: Mac[0][0]
        [0 0 0 0 0 0 0 0 1 1]
        [0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 1 0]
        sage: Mac[0][1]
        [[2, z], [3, z], [3, y]]

    """

    paires_to_process = []
    paires_to_remove = []
    list_already_erased = []
    list_already_added = []
    # Might be faster with a list like list_already_passing_F5


    for paire in paires:
        i_acceptable = True
        j_acceptable = True
        if paire[0] == d :
            paires_to_remove.append(paire)
            signj = [paire[3][0],paire[3][1]*paire[2]]
            signi = [paire[5][0],paire[5][1]*paire[4]]
            if not(signj in list_already_erased) and not(signi in list_already_erased):
                if F5_criterion(signj,G):
                    j_acceptable = False
                    list_already_erased.append(signj)

                if F5_criterion(signi,G):
                    i_acceptable = False
                    list_already_erased.append(signi)

                if i_acceptable and j_acceptable:
                    paires_to_process.append(paire)



    #Erasing the pairs that we have considered
    for paire in paires_to_remove:
        paires.remove(paire)



    #Producing the list of polynomials to be processed
    # We produce a set of the leading monomials



    poly_to_do = []
    list_poly = []
    set_done_monom = set()
    for paire in paires_to_process:
        signj_polyj = [paire[3][0],paire[3][1]*paire[2],paire[2]*paire[3][2]]
        signj = [paire[3][0],paire[3][1]*paire[2]]
        if poly_to_do.count(signj_polyj)<1 and list_already_added.count(signj)<1:
            poly_to_do.append(signj_polyj)
            list_poly.append(signj_polyj)
            list_already_added.append(signj)
        
        #set_done_monom.add(signj_polyj[2].lm())
        
        signi_polyi = [paire[5][0],paire[5][1]*paire[4], paire[4]*paire[5][2]]
        signi = [paire[5][0],paire[5][1]*paire[4]]
        if poly_to_do.count(signi_polyi)<1 and list_already_added.count(signi)<1:
            poly_to_do.append(signi_polyi)
            list_poly.append(signi_polyi)
            list_already_added.append(signi)
        #set_done_monom.add(signi_polyi[2].lm())
    

    #Here, we add the initial polynomials that are of degree d to the matrix.
    #The goal is to prevent redundancy that might happen when an initial polynomial is far from being reduced
    for g in G:
        if g[2].degree() == d and g[1] == 1 and poly_to_do.count(g)<1 and list_already_added.count([g[0],g[1]])<1:
            poly_to_do.append(g)
            list_poly.append(g)
            list_already_added.append([g[0],g[1]])
    
    
    poly_to_do.sort()
    list_poly.sort()
    #Both are equal

    # Rewriting of the polynomials
    len_poly_to_do = len(poly_to_do)
    for ss in range(len_poly_to_do):
        poly_to_do[ss] = Rewritten(R,G, poly_to_do[ss])
        list_poly[ss] = poly_to_do[ss]



    while len(poly_to_do)>0:
        sign_et_poly = poly_to_do.pop(0)
        set_monom = set(sign_et_poly[2].monomials())-set_done_monom
        set_monom = list(set_monom)
        set_monom.sort()
        while len(set_monom)>0:
            mon = set_monom.pop()
            #print "test", mon
            # We do take for mon the biggest monomial available here
            # We look for a g in G and x**alpha such that  mon ==x**alpha* g.lm() with smallest signature
            # In other word, we look for a "reducer." Yet for signature reason, this reducer could be on more reduced than reducing
            smallest_g = 0
            smallest_sign = []
            for g in G:
                #if Reducible_resp_sign([sign_et_poly[0],sign_et_poly[1],mon],g):
                if g[2].lm().divides(mon):
                    glm = g[2].lm()
                    mult = R.monomial_quotient(mon,glm)
                    if len(smallest_sign)==0 and list_already_added.count([g[0],mult*g[1]])<1:
                        if not(F5_criterion([g[0],mult*g[1]],G)):
                            smallest_g = mult*g[2]
                            smallest_sign = [g[0],mult*g[1]]
                    else:
                        if [g[0],mult*g[1]] <= smallest_sign and list_already_added.count([g[0],mult*g[1]])<1:
                            if not(F5_criterion([g[0],mult*g[1]],G)):
                                smallest_g = mult*g[2]
                                smallest_sign = [g[0],mult*g[1]]
            if len(smallest_sign)>0:
                #new_poly = [g[0],mult*g[1],mult*g[2]]
                new_poly = smallest_sign+[smallest_g]
                poly_to_do.append(new_poly)
                poly_to_do.sort()
                list_poly.append(new_poly)
                list_poly.sort()
                list_already_added.append(smallest_sign)
            set_done_monom.add(mon)


    if len(list_poly)>0:
        #print("\nje suis ici...")
        #print(list_poly)
        Mac = Make_Macaulay_Matrix_sparse_dict(list_poly, monom)
        #print(list_monom(R,d))
        #print(Mac[0].str())
    else:
        #Mac = [Matrix(R.base_ring(),0,len(list_monom(R,d)),[]),[]]
        Mac = [Matrix(R.base_ring(),0,0,[]),[]]


    return Mac, paires


def list_indexes(Mat):
    r"""
    
    Computes the list of the indexes of the first non-zero coefficients
    of the rows of the matrix Mat

    
    INPUT:

    - ``Mat`` -- a matrix
    
    OUTPUT:

    a list of integers

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: mat = Matrix(QQ,4,4,[0,1,2,0,0,0,0,3,2,0,0,0,0,0,1,1])
        sage: list_indexes(mat)
        [1, 3, 0, 2]
    """
    n = Mat.nrows()
    m = Mat.ncols()
    listindexes = []
    for i in range(n):
        j=0
        while j<m:
            if Mat[i,j] <>0:
                break
            j=j+1
        listindexes.append(j)
    return listindexes

def RowEchelonMac_and_listpivots(Mac):
    r"""
    
    Computes the row-echelon form of the matrix Mac,
    with no choice of pivot: always the first non-zero entry in
    a column.

    
    INPUT:

    - ``Mac`` -- a matrix
    
    OUTPUT:

    a pair of Mactilde, a matrix under row-echelon form
    (up to permutation of its rows) and a list of the
    indexes of the pivots used for the computation of 
    Mactilde (number of rows plus one if there is no
    pivot on a row)

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: mat = Matrix(QQ,4,4,[0,1,2,0,-1,-2,0,3,2,0,0,0,2,1,3,1])
        sage: u = RowEchelonMac_and_listpivots(mat)
        sage: u[0]
        [  0   1   2   0]
        [  1   0  -4  -3]
        [  0   0   1 3/4]
        [  0   0   0   1]
        sage: u[1]
        [1, 0, 2, 3]
    """
    n = Mac.nrows()
    m = Mac.ncols()
    Mactilde = copy(Mac)
    listpivots = []
    for i in range(n):
        j=0
        while j<m:
            if Mactilde[i,j] <>0:
                break
            j=j+1
        listpivots.append(j)
        if j<m:
            Mactilde.rescale_row(i, Mactilde[i,j]**(-1), start_col=j)
            for l in range(i+1,n):
                if Mactilde[l,j]<>0:
                    Mactilde.add_multiple_of_row(l, i, -Mactilde[l,j] , start_col=j) 
    return Mactilde, listpivots

def Reconstruction(Mactilde,i,monom,sign):
    r"""
    
    Computes the polynomial with signature, corresponding to the row of Mactilde of index i

    
    INPUT:

    - ``Mactilde`` -- a Macaulay matrix (just the matrix)
    - ``i`` -- an integer
    - ``monom`` -- the list of the monomials corresponding to the columns of Mactilde
    - ``sign`` -- the signature corresponding to the row of Mactilde of index i
    
    OUTPUT:

    a polynomial with signature, corresponding to the row of Mactilde of index i

    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: mat = Matrix(QQ,2,3,[0,1,2,1,-1,-2])
        sage: monom = [x,y,z]
        sage: sign = [0,1]
        sage: Reconstruction(mat,1,monom,sign)
        [0, 1, x - y - 2*z]
    """
    r = Mactilde.row(i)
    m = Mactilde.ncols()
    pol = sum(r[j]*monom[j] for j in range(m))
    return sign+[pol]

def Eliminate_unreduced(G,d):
    r"""
    
    In the list of polynomials with signature G, remove the polynomials
    such that their signature is attained by another member of $G$ with smaller
    leading monomial. This is done only for polynomials of degree d.
    

    
    INPUT:

    - ``G`` -- a list of polynomials with signature, ordered by signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    
    OUTPUT:

    - ``newG`` -- an updated list of polynomials with signature, with no redundancy of signature

    EXAMPLES:: 

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: G = [[0,1,x],[1,1,x**2+2*x*y+y**2-z**2],[1,1,y**2-z**2]]
        sage: G2=Eliminate_unreduced(G,2)
        sage: G2
        [[0, 1, x], [1, 1, y^2 - z^2]]        
    """
    long = len(G)
    i = 0
    list_equal = []
    while i<long:
        list_equal_local = [G[i]]
        j = i+1
        while j<long and ([G[j][0],G[j][1]] == [G[i][0],G[i][1]]):
            list_equal_local.append(G[j])
            j=j+1
        list_equal.append(list_equal_local)
        i=j
    newG=[]
    for ll in list_equal:
        long2 = len(ll)
        if long2==1:
            newG.append(ll[0])
        else:
            kmin=0
            min_lm = ll[0][2].lm()
            for k in range(long2):
                if ll[k][2].lm()<min_lm:
                    min_lm=ll[k][2].lm()
                    kmin=k
            newG.append(ll[kmin])
    return newG


def UpdateGpaires(Mac, Mactilde, listpivots, G, paires,d, redundancy_test, monom):
    r"""
    
    In the F5 algorithm, updates the current list of polynomials with
    signatures G, according to the Macaulay matrix Mac,
    its row-echelon form (with no choice of pivot) Mactilde,
    and previous computation.
    paires, the list of the remaining S-pairs to proceed is also
    updated.
    G and paires are updated in place.

    
    INPUT:

    - ``Mac`` -- a Macaulay matrix: a list of a matrix and a list of the signatures corresponding to its rows
    - ``Mactilde`` -- a Macaulay matrix (just the matrix)
    - ``listpivots`` -- the list of the pivots used when computing the row-echelon form Mactilde of Mac
    - ``G`` -- a list of polynomials with signature
    - ``paires`` -- a list of S-pairs of polynomials with signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    - ``redundancy_test`` -- a boolean to check whether a redundancy test is required
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    
    OUTPUT:

    - ``G`` -- an updated list of polynomials with signature
    - ``paires`` -- an updated list of S-pairs of polynomials with signature

    EXAMPLES::  

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: mat = Matrix(QQ,2,6,[1,1,0,0,0,1,1,2,1,0,0,-1])
        sage: Mac = [mat,[[0,1],[1,1]]]
        sage: Mactilde = Matrix(QQ,2,6,[1,1,0,0,0,1,0,1,1,0,0,-2])
        sage: listpivots = [0,1]
        sage: G = [[0,1,x**2+x*y+z**2],[1,1,x**2+2*x*y+y**2-z**2]]
        sage: paires = []
        sage: monom = [x**2, x*y, y**2, x*z, y*z, z**2]
        sage: d = 2
        sage: G2, paires2=UpdateGpaires(Mac, Mactilde, listpivots, G, paires, d,true, monom)
        sage: G2
        [[0, 1, x^2 + x*y + z^2], [1, 1, x*y + y^2 - 2*z^2]]
        sage: paires2
        [[3, 1, x, [1, 1, x*y + y^2 - 2*z^2], y, [0, 1, x^2 + x*y + z^2]]]
    """
    n = Mac[0].nrows()
    m = Mac[0].ncols()

    if n == 0:
        return G,paires


    listindexes = list_indexes(Mac[0])
    R = G[0][2].parent()

    newG = []
    for i in range(n):
        if listindexes[i]<>listpivots[i] and listpivots[i]<m:
            #Beware, 20/01 modification of condition of reconstruction
            if not(Reducible_large_resp_sign_family([Mac[1][i][0],Mac[1][i][1],monom[listpivots[i]]],G)):
                newG.append(Reconstruction(Mactilde,i,monom,[Mac[1][i][0],Mac[1][i][1]]))

    #Beware, because we do not use Rewritting techniques, we recompute the necessary pairs



    G2 = G+newG
    G2.sort()
    if redundancy_test:
        paires = []
        G2=Eliminate_unreduced(G2,d)
        #can be simplified?
        G = [aa for aa in G if G2.count(aa)>0]
        newG = [aa for aa in newG if G2.count(aa)>0]

        long = len(G2)
        for j1 in range(long):
            for j2 in range(j1+1,long):
                g = G2[j1]
                g2 = G2[j2]
                u = g[2].lm().lcm(g2[2].lm())
                du = u.degree()
                if du > d :
                    paires.append([du,g2[0],R(u/g2[2].lm()),g2,R(u/g[2].lm()), g])
    
    else:

        for g in G:
            for g2 in newG:
                i1 = g[0]
                i2 = g2[0]
                gtemp1 = g
                gtemp2 = g2
                if i2< i1:
                    h = gtemp1
                    gtemp1 = gtemp2
                    gtemp2 = h
                u = gtemp1[2].lm().lcm(gtemp2[2].lm())
                du = u.degree()
                if du > d :
                    paires.append([du,gtemp2[0],R(u/gtemp2[2].lm()),gtemp2,R(u/gtemp1[2].lm()), gtemp1])
    
        long = len(newG)
        for j1 in range(long):
            for j2 in range(j1+1,long):
                g = newG[j1]
                g2 = newG[j2]
                u = g[2].lm().lcm(g2[2].lm())
                du = u.degree()
                if du > d :
                    paires.append([du,g2[0],R(u/g2[2].lm()),g2,R(u/g[2].lm()), g])



    if len(newG)>0 or redundancy_test:
        paires.sort()

    return G2, paires
    


def tentative_F5(list1):
    r"""
    Computes a Gröbner bases of the ideal generated by the list of polynomials list1.
    The monomial order is the ambiant monomial order in the polynomial ring containing
    the polynomials of list1

    
    INPUT:

    - ``list1`` -- a list of polynomials

    OUTPUT:

    A Gröbner basis of the ideal generated by list1. For now, it is not reduced


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z]
        sage: G = tentative_F5(list1)
        sage: G
        [x + y + z, y^2 + y*z + z^2, z^3]

    ALGORITHM:

    The following algorithm is adapted from [VY2017]_.

    REFERENCES:

    .. [VY2017] T. Vaccon, K.Yokoyama, A Tropical F5 Algorithm, proceedings of ISSAC 2017.

    """
    #Initialisation
    R = list1[0].parent()
    list_degrees = [aa.degree() for aa in list1]
    list_degrees.sort()
    list = copy(list1)
    s = len(list)
    listsign = []
    for i in range(s):
        listsign.append([i,1,list[i]])
    G = listsign
    l = len(G)
    paires = []
    for i in range(s):
        for j in range(i+1,l):
            u = list[i].lm().lcm(list[j].lm())
            d = u.degree()
            paires.append([d,j,R(u/list[j].lm()),listsign[j],R(u/list[i].lm()), listsign[i]])
            #watch out: j>i
    ######attention
    paires.sort()

    #Corps du programme   
    d = 1
    monom = list_monom_iteratif_no_repeat(R,d,[])
    while paires != [] :
        ###### No Rewriting
        # Beware for d, not necessarily increasing...
        # We might take d = paires[0][0], but this leads to (harmless ?) unnecessary computations
        #print "d=",d
        #print "remaining pairs=", len(paires)
        #print "paires", paires
        Mac, paires = SymbolicPreprocessing(R,G, d,paires, monom)
        #print "Macaulay matrix size=", [Mac[0].nrows(),Mac[0].ncols()]
        #print "je reduis"
        Mactilde, listpivots = RowEchelonMac_and_listpivots(Mac[0])
        # Beware, no pair of degree <=d will be considered anymore
        redundancy_check = (list_degrees.count(d)>0)
        #print "je fais l'update"
        G, paires = UpdateGpaires(Mac, Mactilde, listpivots, G, paires,d,redundancy_check, monom)
        #print "forme de G"
        #for g in G:
        #    print [g[0],g[1],g[2].lm()]
        d = d+1
        #print "calculs des monomes de degre d"
        monom = list_monom_iteratif_no_repeat(R,d,monom)
    G = [g[2] for g in G]
    return G


def tentative_F5_reduced(list1):
    r"""
    Computes the reduced Gröbner bases of the ideal generated by the list of polynomials list1.
    The monomial order is the ambiant monomial order in the polynomial ring containing
    the polynomials of list1

    
    INPUT:

    - ``list1`` -- a list of polynomials

    OUTPUT:

    The reduced Gröbner basis of the ideal generated by list1.


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z, z**4]
        sage: G = tentative_F5_reduced(list1)
        sage: G
        [z^3, y^2 + y*z + z^2, x + y + z]

    """
    G = tentative_F5(list1)
    G2 = Sequence(G)
    G3 = G2.reduced()
    return G3
