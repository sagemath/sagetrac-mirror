#This module implements variants of the F5 algorithm, with a view for implementation over finite-precision complete discrete valuation field.

#AUTHOR:

#- Tristan Vaccon (2016-2017)


#*****************************************************************************
#  Copyright (C) 2015 Tristan Vaccon <tristan.vaccon@univ-rennes1.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************




def list_monom(P,d):
    r"""

    Return an ordered list of the monomials of P of degree d
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import list_monom
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

    Return an ordered list of the monomials of P of degree d
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer
    - ``ld1`` -- the list of the monomials of degree d-1

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import list_monom
        sage: list_monom_dminus(S,2n,[x,y,z])
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

    Return the index of the last variable dividing mon

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import list_monom
        sage: index_last_variable(S,x*y)
        2
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

    Return an ordered list of the monomials of P of degree d
 
    INPUT:

    - ``P`` -- a polynomial ring
    - ``d`` -- an integer
    - ``ld1`` -- the list of the monomials of degree d-1

    OUTPUT:

    The list of monomials in P of degree d, ordered from the biggest to the smallest according to the monomial ordering on P

    EXAMPLES::

        sage: S.<x,y,z>=Qp(5)[]
        sage: from sage.rings.polynomial.padics.toy_MatrixF5 import list_monom
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

def Make_Macaulay_Matrix(list_poly):
    #Mac will be a couple of a matrix and a list of signatures
    # We assume list_poly is ordered by signature and compatible
    r"""

    Make the Macaulay matrix for the polynomials of degree d in list_poly
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials with signature

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix
        sage: list_poly = [[0,1,x^2+y^2], [1,1,x*y], [2,1,y*z], [3,1,z^2]]
        sage: Mac = Make_Macaulay_Matrix(list_poly)
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
    monom = list_monom(R,d)
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
    #Mac will be a couple of a matrix and a list of signatures
    # We assume list_poly is ordered by signature and compatible
    r"""

    Make the sparse Macaulay matrix for the polynomials of degree d in list_poly
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix
        sage: list_poly = [[0,1,x^2+y^2], [1,1,x*y], [2,1,y*z], [3,1,z^2]]
        sage: Mac = Make_Macaulay_Matrix_sparse_dict(list_poly,[x^2, x*y, y^2, x*z, y*z, z^2])
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

def Make_Macaulay_Matrix_no_sign(list_poly):
    r"""

    Make the Macaulay matrix for the polynomials of degree d in list_poly
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix
        sage: list_poly = [x^2+y^2,x*y,y*z,z^2]
        sage: Mac = Make_Macaulay_Matrix_no_sign(list_poly)
        sage: Mac
        [1 0 1 0 0 0]
        [0 1 0 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]        
    """

    R = list_poly[0].parent()
    d = list_poly[0].degree()
    monom = list_monom(R,d)
    l = len(monom)

    listsign = []
    Mac = Matrix(R.base_ring(), 0, l, [])
    for f in list_poly:
        list = f.monomials()
        flist = [ R(0) ] * l
        for mon in list:
            if mon in monom:
                i = monom.index(mon)
                flist[i] = f.monomial_coefficient(mon)
        fmat = Matrix(R.base_ring(), 1, l, flist)
        Mac = Mac.stack(fmat)
    return Mac

def Make_Macaulay_Matrix_sparse_dict_no_sign(list_poly, monom):
    #Mac will be a couple of a matrix and a list of signatures
    # We assume list_poly is ordered by signature and compatible
    r"""

    Make the sparse Macaulay matrix for the polynomials of degree d in list_poly
 
    INPUT:

    - ``list_poly`` -- a list of degree d polynomials
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly

    OUTPUT:

    The Macaulay matrix defined by list_poly

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5 import Make_Macaulay_Matrix
        sage: list_poly = [[0,1,x^2+y^2], [1,1,x*y], [2,1,y*z], [3,1,z^2]]
        sage: Mac = Make_Macaulay_Matrix_sparse_dict(list_poly,[x^2, x*y, y^2, x*z, y*z, z^2])
        sage: Mac[0]
        [1 0 1 0 0 0]
        [0 1 0 0 0 0]
        [0 0 0 0 1 0]
        [0 0 0 0 0 1]
        sage: Mac[1]
        [[0, 1], [1, 1], [2, 1], [3, 1]]
        
    """

    R = list_poly[0].parent()
    d = list_poly[0].degree()
    l = len(monom)

    dict_monom = {monom[aa]:aa for aa in range(l)}


    nrows = len(list_poly)


    Mac = MatrixSpace(R.base_ring(),nrows,l,sparse=True)(0)
    for u in range(nrows):
        f = list_poly[u]
        list = f.monomials()
        for mon in list:
            if dict_monom.has_key(mon):
                j = dict_monom[mon]
                Mac[u,j] = f.monomial_coefficient(mon)
    return Mac

def LT_trop(f,R,w):
    r"""
    
    Computes the leading term of f according to the tropical term order
    given by the monomial ordering on R and the weight w
 
    INPUT:

    - ``f`` -- a polynomial
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    the leading term of f  

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5trop import LM_trop
        sage: f = 2*x^2+3*y^2+z^2
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: LT_trop(f,R1,w1)
        3*y^2
        sage: LT_trop(f,R2,w2)
        z^2
    """
    lmon = f.monomials()
    lcoeff = f.coefficients()
    #print(lmon)
    nvar = len(f.parent().gens())
    K = R.base_ring()
    if len(lmon) == 0:
        return f
    else:
        imin = 0
        mon = lmon[0]
        trop_weight0 = K(lcoeff[0]).valuation()+sum([lmon[0].degrees()[vv]*w[vv] for vv in range(nvar)])
        for uu in range(len(lmon)):
            trop_weight = K(lcoeff[uu]).valuation()+sum([lmon[uu].degrees()[ww]*w[ww] for ww in range(nvar)])
            if trop_weight < trop_weight0 or (trop_weight == trop_weight0 and R(lmon[uu])> R(mon)):
                imin = uu
                mon = lmon[uu]
                trop_weight0 = trop_weight
    return lcoeff[imin]*lmon[imin]

def LM_trop(f,R,w):
    r"""
    
    Computes the monomial of the leading term of f according to the tropical term order
    given by the monomial ordering on R and the weight w
 
    INPUT:

    - ``f`` -- a polynomial
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a monomial of f, the monomial of the leading term of f  

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5trop import LM_trop
        sage: f = 2*x^2+3*y^2+z^2
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: LM_trop(f,R1,w1)
        y^2
        sage: LM_trop(f,R2,w2)
        z^2
    """
    lmon = f.monomials()
    lcoeff = f.coefficients()
    nvar = len(f.parent().gens())
    K = R.base_ring()
    if len(lmon) == 0:
        return f
    else:
        imin = 0
        mon = lmon[0]
        trop_weight0 = K(lcoeff[0]).valuation()+sum([lmon[0].degrees()[vv]*w[vv] for vv in range(nvar)])
        for uu in range(len(lmon)):
            trop_weight = K(lcoeff[uu]).valuation()+sum([lmon[uu].degrees()[ww]*w[ww] for ww in range(nvar)])
            if trop_weight < trop_weight0 or (trop_weight == trop_weight0 and R(lmon[uu])> R(mon)):
                imin = uu
                mon = lmon[uu]
                trop_weight0 = trop_weight
    return lmon[imin]


def LC_trop(f,R,w):
    r"""
    
    Computes the coefficient of the leading term of f according to the tropical term order
    given by the monomial ordering on R and the weight w
 
    INPUT:

    - ``f`` -- a polynomial
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a coefficient of f, the coefficient of the leading term of f  

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5trop import LM_trop
        sage: f = 2*x^2+3*y^2+z^2
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: LC_trop(f,R1,w1)
        3
        sage: LC_trop(f,R2,w2)
        1
    """
    lmon = f.monomials()
    lcoeff = f.coefficients()
    nvar = len(f.parent().gens())
    K = R.base_ring()
    if len(lmon) == 0:
        return f
    else:
        imin = 0
        mon = lmon[0]
        trop_weight0 = K(lcoeff[0]).valuation()+sum([lmon[0].degrees()[vv]*w[vv] for vv in range(nvar)])
        for uu in range(len(lmon)):
            trop_weight = K(lcoeff[uu]).valuation()+sum([lmon[uu].degrees()[ww]*w[ww] for ww in range(nvar)])
            if trop_weight < trop_weight0 or (trop_weight == trop_weight0 and R(lmon[uu])> R(mon)):
                imin = uu
                mon = lmon[uu]
                trop_weight0 = trop_weight
    return lcoeff[imin]


def F5_criterion_trop(sign,G,R,w):
    r"""
    
    Inside the F5 algorithm
    Check the F5 criterion for the signature sign = (i, x^alpha) and G
    an S-GB up to degree deg(f_i)+deg(x^alpha) and signature strictly less than (i,1)
 
    INPUT:

    - ``sign`` -- a (potential) signature of a polynomial
    - ``G`` -- a list of polynomials with signature
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a boolean, testing whether the F5 criterion for sign and G is passed

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.padics.toy_F5trop import F5_criterion
        sage: sign = [1,x^3]
        sage: G =  [[0,1,x^2], [1,1,x*y], [2,1,y*z], [3,1,2*x^2+3*y^2+z^2]]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: F5_criterion_trop(sign,G, R1,w1)
        True
        sage: sign2 = [0,x^2]
        sage: sign3 = [4,y^2]
        sage: F5_criterion_trop(sign2,G, R2, w2)
        False
        sage: F5_criterion_trop(sign3,G, R1, w1)
        True
        sage: F5_criterion_trop(sign3,G, R2, w2)
        False
        
    """
    for g in G:
        ei = g[0]
        glm = LM_trop(g[2],R,w)
        if ei < sign[0]:
            if parent(glm).monomial_divides(glm,sign[1]):
                return True
    return False 

def Reducible_resp_sign_trop(f,g,R,w):
    r"""
    
    Inside the F5 algorithm, f and g are of the form [i,x^alpha, poly].
    Tests whether g can reduce f while respecting their signatures with strict inequality
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``g`` -- a polynomial with signature
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a boolean, testing whether g can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: from sage.rings.polynomial.padics.toy_F5trop import Reducible_resp_sign
        sage: g = [0,x,x^2]
        sage: f = [1,x^2,x^3*y]
        sage: f2 = [0,x^2,x^3]
        sage: f3 = [0,1,x^3]
        sage: Reducible_resp_sign_trop(f,g,R1,w1)
        True
        sage: Reducible_resp_sign_trop(f2,g,R1,W1)
        False
        sage: Reducible_resp_sign_trop(f3,g,R2,w2)
        False
    """
    fpol = f[2]
    R0 = fpol.parent()
    flm = LM_trop(f[2],R,w)
    glm = LM_trop(g[2],R,w)
    nvar = len(fpol.parent().gens())
    trop_weight_sign = sum([f[1].degrees()[vv]*w[vv] for vv in range(nvar)])
    if parent(glm).monomial_divides(glm,flm):
        test_sign = R0.monomial_quotient(flm,glm)*g[1]
        trop_weight_test_sign = sum([test_sign.degrees()[vv]*w[vv] for vv in range(nvar)])
        bool_comp_sign =( trop_weight_test_sign > trop_weight_sign or (trop_weight_test_sign == trop_weight_sign and R(test_sign) < R(f[1])))
        return (f[0]> g[0] or (f[0] == g[0] and bool_comp_sign))
    return False

def Reducible_large_resp_sign_trop(f,g,R,w):
    r"""
    
    Inside the F5 algorithm, f and g are of the form [i,x^alpha, poly].
    Tests whether g can reduce f while respecting their signatures with non-strict inequality
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``g`` -- a polynomial with signature
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a boolean, testing whether g can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: from sage.rings.polynomial.padics.toy_F5 import Reducible_large_resp_sign
        sage: g = [0,x,x^2]
        sage: f = [1,x^2,x^3*y]
        sage: f2 = [0,x^2,x^3]
        sage: f3 = [0,1,x^3]
        sage: Reducible_large_resp_sign_trop(f,g,R1,w1)
        True
        sage: Reducible_large_resp_sign_trop(f2,g,R1,w1)
        True
        sage: Reducible_large_resp_sign_trop(f3,g,R2,w2)
        False
    """
    R0 = f[2].parent()
    flm = LM_trop(f[2],R,w)
    glm = LM_trop(g[2],R,w)
    nvar = len(f[2].parent().gens())
    trop_weight_sign = sum([f[1].degrees()[vv]*w[vv] for vv in range(nvar)])
    if parent(glm).monomial_divides(glm,flm):
        test_sign = R0.monomial_quotient(flm,glm)*g[1]
        trop_weight_test_sign = sum([test_sign.degrees()[vv]*w[vv] for vv in range(nvar)])
        bool_comp_sign = trop_weight_test_sign > trop_weight_sign or (trop_weight_test_sign == trop_weight_sign and R(test_sign) <= R(f[1]))
        return (f[0]> g[0] or (f[0] == g[0] and bool_comp_sign))
    return False


def Reducible_large_resp_sign_family_trop(f,G,R,w):
    r"""
    
    Inside the F5 algorithm, f and the elements of the list G are of the form [i,x^alpha, poly].
    Tests whether G can reduce f while respecting their signatures
    
    INPUT:

    - ``f`` -- a polynomial with signature
    - ``G`` -- a list of polynomials with signature
    - ``R`` -- a polynomial ring of dimension
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a boolean, testing whether G can reduce f while respecting their signatures

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: R2 = PolynomialRing(Qp(3),['x2','y2','z2'] , order = 'lex') 
        sage: w1 = [0,0,0]
        sage: w2 = [0,0,-1]
        sage: from sage.rings.polynomial.padics.toy_F5trop import all
        sage: G = [[0,1,x^2], [1,1,x*y], [2,1,y*z], [3,1,z^2]]
        sage: f = [1,x^2,x^3*y]
        sage: f2 = [0,x^2,x^3]
        sage: f3 = [0,1,x^3]
        sage: Reducible_large_resp_sign_family_trop(f,G,R1,w1)
        True
        sage: Reducible_large_resp_sign_family_trop(f2,G,R1,w1)
        True
        sage: Reducible_large_resp_sign_family_trop(f3,G,R2,w2)
        False
    """
    for g in G:
        if Reducible_large_resp_sign_trop(f,g,R,w):
            return True
    return False

def Rewritten(R,G,g, R1,w):
    r"""
    
    Find the latest element of G such it has a multiple with same signature as g.

    
    INPUT:

    - ``R`` -- a polynomial ring
    - ``G`` -- a list of polynomials with signature, ordered by signature
    - ``g`` -- a polynomial with signature
    - ``R1`` -- a polynomial ring
    - ``G`` -- a list of polynomials with signature

    
    OUTPUT:

    A multiple of an element of G with same signature and leading monomial
    as g. It is obtained with the latest element of g possible

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: w1 = [0,0,0]
        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: G = [[0,y,x^2+2*z^2], [1,x,x*y+y^2], [2,1,y*z+z^2], [3,1,z^2]]
        sage: g = [2,x,2*x*y*z+x*z^2+y^2*z]
        sage: Rewritten(G,g)
        [2,x,x*y*z+x*z^2]

    """
    mon_sign = LM_trop(R(g[1]),R1,w)
    newg=g
    for gg in G:
        if gg[0]==g[0] and R.monomial_divides(LM_trop(R(gg[1]),R1,w),mon_sign):
            glm = LM_trop(R(gg[1]),R1,w)
            mult = R.monomial_quotient(mon_sign,glm)
            newg = [gg[0],mult*gg[1],mult*gg[2]]
    return newg

def SymbolicPreprocessingTrop(R,w,G, d, paires, monom):
    r"""
    
    Constructs the Macaulay matrix from the pairs in paires,
    while using the F5 criterion

    
    INPUT:

    - ``R`` -- a polynomial ring
    - ``w`` -- a list of n rational numbers
    - ``G`` -- a list of polynomials with signature
    - ``d`` -- an integer, representing the degree of the processed polynomials
    - ``paires`` -- a list of pairs of polynomials with signature
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    
    OUTPUT:

    A Macaulay matrix, as a couple of a matrix built from polynomials and
    a list of signatures for the rows of this matrix 

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: R1.<x1,y1,z1> = Qp(2)[]
        sage: w1 = [0,0,0]
        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: G = [[0,y,x^2+2*z^2], [1,x,x*y+y^2], [2,S(1),y*z+z^2], [3,S(1),z^2]]
        sage: paires = []
        sage: paires.append([3,1,x,[1,x,x*y+y^2],y, [0,S(1),x^2+2*z^2]])
        sage: paires.append([3,3,y,[3,S(1),x*z+3*y*z+z^2],z, [3,S(1),x*y+2*y*z+z^2]])
        sage: d = 3
        sage: Mac = SymbolicPreprocessingTrop(R1,w1, G, d, paires)
        sage: Mac[0][0]
        [0 0 0 0 0 0 0 0 1 1]
        [0 0 0 0 0 0 1 0 1 0]
        [0 0 0 0 0 1 0 0 2 1]
        [0 0 0 0 0 1 3 0 1 0]
        sage: Mac[0][1]
        [[2, z], [2, y], [3, z], [3, y]]

    """
    R0 = G[0][2].parent()

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
                if F5_criterion_trop(signj,G,R,w):
                    j_acceptable = False
                    list_already_erased.append(signj)

                if F5_criterion_trop(signi,G,R,w):
                    i_acceptable = False
                    list_already_erased.append(signi)

                if i_acceptable and j_acceptable:
                    paires_to_process.append(paire)



    #Erasing the pairs that we have just considered
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
        
        set_done_monom.add(LM_trop(signj_polyj[2],R,w))
        
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
    
    
    # Refaire le sort??? Peut§etre pas necessaire
    poly_to_do.sort()
    list_poly.sort()

    #Rewritting

    len_poly_to_do = len(poly_to_do)
    for ss in range(len_poly_to_do):
        poly_to_do[ss] = Rewritten(R0,G, poly_to_do[ss],R,w)
        list_poly[ss] = poly_to_do[ss]
    

    while len(poly_to_do)>0:
        sign_et_poly = poly_to_do.pop(0)
        set_monom = set(sign_et_poly[2].monomials())-set_done_monom
        set_monom = list(set_monom)
        set_monom.sort()

        while len(set_monom)>0:
            mon = set_monom.pop()

            # We do take for mon the biggest monomial available here
            # We look for a g in G and x^alpha such that  mon ==x^alpha* g.lm() with smallest signature
            # In other word, we look for a "reducer." Yet for signature reason, this reducer could be on more reduced than reducing
            smallest_g = 0
            smallest_sign = []
            for g in G:
                #if Reducible_resp_sign([sign_et_poly[0],sign_et_poly[1],mon],g):
                if parent(mon).monomial_divides(LM_trop(g[2],R,w),mon):
                    glm =  LM_trop(g[2],R,w)
                    mult = parent(mon).monomial_quotient(mon,glm)
                    if len(smallest_sign)==0 and list_already_added.count([g[0],mult*g[1]])<1:
                        smallest_g = mult*g[2]
                        smallest_sign = [g[0],mult*g[1]]
                    else:
                        if [g[0],mult*g[1]] <= smallest_sign and list_already_added.count([g[0],mult*g[1]])<1:
                            smallest_g = mult*g[2]
                            smallest_sign = [g[0],mult*g[1]]

            if len(smallest_sign)>0:
                #new_poly = [g[0],mult*g[1],mult*g[2]]
                list_already_added.append(smallest_sign)
                new_poly = smallest_sign+[smallest_g]
                poly_to_do.append(new_poly)
                poly_to_do.sort()
                list_poly.append(new_poly)
                list_poly.sort()
            set_done_monom.add(mon)



    if len(list_poly)>0:
        #print("\nje suis ici...")
        #print(list_poly)
        Mac = Make_Macaulay_Matrix_sparse_dict(list_poly, monom)
        #print(list_monom(R,d))
        #print(Mac[0].str())
    else:
        Mac = [Matrix(R.base_ring(),0,0,[]),[]]

    return Mac, paires

def list_indexes_trop(Mat, monom, R,w):
    r"""
    
    Computes the list of the indexes of the leading coefficients
    of the rows of the Macaulay matrix Mat, given with
    the list of monomials monom indexing its columns
    and R and w defining a tropical term order on the
    polynomial ring in which the monomials of monom live.

    
    INPUT:

    - ``Mat`` -- a matrix
    - ``monom`` -- a list of monomials
    - ``R`` -- a polynomial ring of dimension n
    - ``w`` -- a list of n rational numbers
    
    OUTPUT:

    a list of integers

    EXAMPLES::

        sage: S.<x,y,z,t>=QQ[]
        sage: R1.<x1,y1,z1,t1> = Qp(2)[]
        sage: w1 = [0,0,0,0]
        sage: mat = Matrix(QQ,4,4,[0,1,2,0,0,0,0,3,2,0,0,0,0,0,1,1])
        sage: monom = [x,y,z,t]
        sage: list_indexes_trop(mat, monom, R1, w1)
        [1,3,0,2]
        sage: w2 = [2,9,1,2]
        sage: list_indexes_trop(mat, monom, R1, w2)
        [2,3,0,2]
    """
    n = Mat.nrows()
    m = Mat.ncols()
    listindexes = []
    T = R.base_ring()
    nvar = len(w)
    for i in range(n):
        #Looking for the LT for the tropical order
        max_term = None
        max_col = None
        entry_non_zero = None
        for j in range(m):
            if max_term is None and Mat[i,j] != 0:
                max_term = Mat[i,j]
                max_col = j
                max_monomial = monom[j]
                max_trop_weight = T(Mat[i,j]).valuation()+sum([monom[j].degrees()[uu]*w[uu] for uu in range(nvar)])
            elif Mat[i,j] != 0:
                term_trop_weight = T(Mat[i,j]).valuation()+sum([monom[j].degrees()[uu]*w[uu] for uu in range(nvar)])
                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and R(monom[j])> R(max_monomial)):
                    max_term = Mat[i,j]
                    max_col = j
                    max_trop_weight = term_trop_weight
                    max_monomial = monom[j]
        if max_term is None:
            listindexes.append(m)
        else:
            listindexes.append(max_col)
    return listindexes

def column_transposition(mat,list_monom,j1,j2):
    """
    performs the transposition of the columns j1 and j2 of the Macaulay matrix given by mat and list_monom.

    EXAMPLES::

        sage: T.<x,y,z> = QQ[]
        sage: Mac = Matrix(QQ,1,3,[1,2,3])
        sage: list_monom = [x,y,z]
        sage: column_transposition(Mac,list_monom,1,2)
        sage: Mac
        [1,3,2]
        sage: list_monom
        [x,z,y]
        
    """
    if j1 != j2:
        nrows = mat.nrows()
        ncols = mat.ncols
        temp = list_monom[j1]
        list_monom[j1] = list_monom[j2]
        list_monom[j2] = temp
        mat.swap_columns(j1,j2)


def Trop_RowEchelonMac_and_listpivots(Mac,list_monom, R, tropical_weight):
    r"""
    
    Computes the row-echelon form of the Macaulay matrix Mac,
    with no choice of pivot: always the LM of each row is used as pivot

    
    INPUT:

    - ``Mac`` -- a matrix
    - ``list_monom`` -- a list of the monomials of degree d in R
    - ``R`` -- a polynomial ring
    - ``tropical_weight`` -- a tropical weight, i.e. an n-tuple of real numbers
    
    OUTPUT:

    a list consisting of Mactilde, a matrix under row-echelon form
    (up to permutation of its rows), a list of the
    indexes of the pivots used for the computation of 
    Mactilde (number of rows plus one if there is no
    pivot on a row) and a list 

    EXAMPLES::
        sage: S.<x,y,z,t> = QQ[]
        sage: R.<xx,yy,zz,tt> = Qp(2,50)[]
        sage: w = [-1,0,2,-3]
        sage: list_monom = [x,y,z,t]
        sage: mat = Matrix(QQ,4,4,[0,1,2,0,-1,-2,0,3,2,0,0,0,2,1,3,1])
        sage: u = Trop_RowEchelonMac_and_listpivots(mat, list_monom,R,w)
        sage: u[0]
        [   1    0    0    2]
        [   0    1 -1/3  4/3]
        [   0    0    1    0]
        [   0    0    0    1]
        sage: u[1]
        [y, t, x, z]

    """
    n = Mac.nrows()
    m = Mac.ncols()
    Mactilde = copy(Mac)
    listpivots = []
    T = R.base_ring()
    nvar = len(tropical_weight)
    index_column = 0
    for i in range(n):
        #Looking for the LT for the tropical order
        max_term = None
        max_col = None
        entry_non_zero = None
        for j0 in range(index_column,m):
            if max_term is None and Mactilde[i,j0] != 0:
                max_term = Mactilde[i,j0]
                max_col = j0
                max_monomial = list_monom[j0]
                max_trop_weight = T(Mactilde[i,j0]).valuation()+sum([list_monom[j0].degrees()[uu]*tropical_weight[uu] for uu in range(nvar)])
            elif Mactilde[i,j0] != 0:
                term_trop_weight = T(Mactilde[i,j0]).valuation()+sum([list_monom[j0].degrees()[uu]*tropical_weight[uu] for uu in range(nvar)])
                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and R(list_monom[j0])> R(max_monomial)):
                    max_term = Mactilde[i,j0]
                    max_col = j0
                    max_trop_weight = term_trop_weight
                    max_monomial = list_monom[j0]
        if max_term is None:
            listpivots.append(m)
        else:
            listpivots.append(index_column)
            column_transposition(Mactilde, list_monom, index_column, max_col) 
            Mactilde.rescale_row(i, Mactilde[i,index_column]^(-1), start_col=index_column)
            for l in range(i+1,n):
                if Mactilde[l,index_column] != 0:
                    Mactilde.add_multiple_of_row(l, i, -Mactilde[l,index_column] , start_col=index_column)
            index_column +=1 
    return Mactilde, list_monom, listpivots

def Reconstruction(Mactilde,i,monom,sign):
    r"""
    
    Computes a polynomial with signature, 
    corresponding to the row of Mactilde of index i


    
    INPUT:

    - ``Mactilde`` -- a Macaulay matrix (just the matrix)
    - ``i`` -- an integer
    - ``monom`` -- the list of the monomials corresponding to the columns of Mactilde
    - ``sign`` -- the signature corresponding to the row of Mactilde of index i
    
    OUTPUT:

    a polynomial with signature, corresponding to the row of Mactilde of index i

    EXAMPLES::

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

def Eliminate_unreduced(G,F,d):
    r"""
    
    In the list of polynomials with signature G, remove the polynomials
    of degree d that are both in G and F, leading to a redundancy of signature    

    
    INPUT:

    - ``G`` -- a list of polynomials with signature, ordered by signature
    - ``F`` -- a list of polynomials with signature, ordered by signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    
    OUTPUT:

    - ``newG`` -- an updated list of polynomials with signature

    EXAMPLES:: 

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: F = [[0,1,x],[1,1,x^2+2*x*y+y^2-z^2]]
        sage: G = [[0,1,x],[1,1,x^2+2*x*y+y^2-z^2],[1,1,y^2-z^2]]
        sage: Eliminate_unreduced(G,F,2)
        sage: G
        [[0,1,x],[1,1,y^2-z^2]]        
    """

    newGd = []
    newG = []
    newGsignd = []
    for g in G:
         newG.append(g)
         if g[2].degree() == d:
             newGd.append(g)
             newGsignd.append([g[0],g[1]])
    tempF = []
    tempFsign = []
    for f in F:
         if f[2].degree() == d:
             tempF.append(f)
             tempFsign.append([f[0],f[1]])
    dl = len(tempF)
    for i in range(dl):
        if newGsignd.count(tempFsign[i])>1:
             if newGd.count(tempF[i])>0:
                 newG.remove(tempF[i])
    return newG       

def UpdateGpaires(Mac, monom, Mactilde, lmonom_mactilde, listpivots, G, paires,d, R1, w, redundancy_test, F):
    # Compute the new G and the new set of admissible pairs
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
    - ``monom`` -- the list of the monomials corresponding to the columns of Mac
    - ``Mactilde`` -- a Macaulay matrix (just the matrix), tropical row-echelon form of Mac
    - ``lmonom_mactilde`` -- the list of the monomials corresponding to the columns of Mac
    - ``listpivots`` -- the list of the pivots used when computing the tropical row-echelon form Mactilde of Mac
    - ``G`` -- a list of polynomials with signature
    - ``paires`` -- a list of S-pairs of polynomials with signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers
    - ``redundancy_test`` -- a boolean to check whether a redundancy test is required
    - ``F`` -- a list of polynomials with signature, ordered by signature 
   
    OUTPUT:

    - ``G`` -- an updated list of polynomials with signature
    - ``paires`` -- an updated list of S-pairs of polynomials with signature

    EXAMPLES::  

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: mat = Matrix(QQ,2,6,[2,1,0,0,0,4,0,1,1,0,0,-1])
        sage: monom = [x^2,x*y,y^2,x*z,y*z,z^2]
        sage: Mac = [mat,[[0,S(1)],[1,S(1)]]]
        sage: mat2 = Matrix(QQ,2,6,[1,0,2,0,0,4,0,1,-2,0,0,-9])
        sage: lmonom_mat2 = [x*y,y^2, x^2,x*z,y*z,z^2]
        sage: listpivots = [0,1]
        sage: G = [[0,1,2*x^2+x*y+4*z^2],[1,1,x*y+y^2-z^2]]
        sage: paires = []
        sage: d = 2
        sage: UpdateGpaires(Mac, monom, mat2, lmonom_mat2, listpivots, G, paires, d, R,w)
        sage: G
        [[0, 1, 2*x^2 + x*y + 4*z^2],
         [1, 1, x*y + y^2 - z^2],
         [1, 1, -2*x^2 + y^2 - 9*z^2]]
        sage: paires
        [[3, 1, x, [1, 1, -2*x^2 + y^2 - 9*z^2], y, [0, 1, 2*x^2 + x*y + 4*z^2]],
         [3, 1, x, [1, 1, -2*x^2 + y^2 - 9*z^2], y, [1, 1, x*y + y^2 - z^2]]]
        
    """
    n = Mac[0].nrows()
    m = Mac[0].ncols()

    if n == 0:
        return G, paires 

    listindexes = list_indexes_trop(Mac[0], monom,R1,w)
    R = G[0][2].parent()

    newG = []
    for i in range(n):
        if listpivots[i]<m and monom[listindexes[i]]<>lmonom_mactilde[listpivots[i]]:
            #Beware on whether this condition is necessary or not
            if not(Reducible_large_resp_sign_family_trop([Mac[1][i][0],Mac[1][i][1],lmonom_mactilde[listpivots[i]]],G,R1,w)):
                newG.append(Reconstruction(Mactilde,i,lmonom_mactilde,[Mac[1][i][0],Mac[1][i][1]]))

    G2 = G+newG
    G2.sort()
    if redundancy_test:
        paires = []
        G2=Eliminate_unreduced(G2,F,d)
        #can be simplified?
        G = [aa for aa in G if G2.count(aa)>0]
        newG = [aa for aa in newG if G2.count(aa)>0]

        long = len(G2)
        for j1 in range(long):
            for j2 in range(j1+1,long):
                g = G2[j1]
                g2 = G2[j2]
                u = R.monomial_lcm(LM_trop(g[2],R1,w),LM_trop(g2[2],R1,w))
                du = u.degree()
                if du > d :
                    paires.append([du,g2[0],R.monomial_quotient(u,LM_trop(g2[2],R1,w)),g2,R.monomial_quotient(u,LM_trop(g[2],R1,w)), g])

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
                u = R.monomial_lcm(LM_trop(gtemp1[2],R1,w),LM_trop(gtemp2[2],R1,w))
                du = u.degree()
                if du > d :
                    paires.append([du,gtemp2[0],R.monomial_quotient(u,LM_trop(gtemp2[2],R1,w)),gtemp2,R.monomial_quotient(u,LM_trop(gtemp1[2],R1,w)), gtemp1])
    
        long = len(newG)
        for j1 in range(long):
            for j2 in range(j1+1,long):
                g = newG[j1]
                g2 = newG[j2]
                u = R.monomial_lcm(LM_trop(g[2],R1,w),LM_trop(g2[2],R1,w))
                du = u.degree()
                if du > d :
                     paires.append([du,g2[0],R.monomial_quotient(u,LM_trop(g2[2],R1,w)),g2,R.monomial_quotient(u,LM_trop(g[2],R1,w)), g])
    
    if len(newG)>0 or redundancy_test:
        paires.sort()

    return G2, paires
    


def tentative_F5Trop(list1,R1,w):
    r"""
    Computes a tropical Gröbner bases of the ideal generated by the list of polynomials list1.
    The tropical term order is defined by the valuation on the base ring of R1, its ambiant
    monomial order, and the tropical weight w.

    
    INPUT:

    - ``list1`` -- a list of polynomials
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A Gröbner basis of the ideal generated by list1. For now, it is not reduced


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z]
        sage: G = tentative_F5Trop(list1,R,w)
        sage: G

    """
    #Initialisation
    R = list1[0].parent()
    list_degrees = [aa.degree() for aa in list1]
    list_degrees.sort()
    list = copy(list1)
    s = len(list)
    listsign = []
    for i in range(s):
        listsign.append([i,R(1),list[i]])
    G = listsign
    l = len(G)
    paires = []
    for i in range(s):
        for j in range(i+1,l):
            u = R.monomial_lcm(LM_trop(list[i],R1,w),LM_trop(list[j],R1,w))
            d = u.degree()
            paires.append([d,j,R.monomial_quotient(u,LM_trop(list[j],R1,w)),listsign[j],R.monomial_quotient(u,LM_trop(list[i],R1,w)), listsign[i]])
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
        print "d=",d
        #print "remaining pairs=", len(paires)
        Mac, paires = SymbolicPreprocessingTrop(R1,w,G, d, paires, monom)
        #print "Macaulay matrix size=", [Mac[0].nrows(),Mac[0].ncols()]
        #print "je reduis"
        Mactilde, list_monom_mactilde, listpivots = Trop_RowEchelonMac_and_listpivots(Mac[0],monom, R1, w)
        # Beware, no pair of degree <=d will be considered anymore
        redundancy_check = (list_degrees.count(d)>0)
        #print "je fais l'update"
        G, paires = UpdateGpaires(Mac, monom, Mactilde,list_monom_mactilde, listpivots, G, paires,d, R1,w,redundancy_check, listsign )
        #print "longueur de G", len(G)
        #print "calculs des monomes de degre d"
        d = d+1
        monom = list_monom_iteratif_no_repeat(R,d,monom)
    G = [g[2] for g in G]
    return G


###############################################
# Classical reduction for computed reduced GB #
###############################################

def minimalization(G,R1,w):
    r"""
    Computes a minimal tropical Gröbner basis of the ideal generated by the list of polynomials in G.
    The tropical term order is defined by the valuation on the base ring of R1, its ambiant
    monomial order, and the tropical weight w.
    G is a tropical Gröbner basis (for the same term order).

    
    INPUT:

    - ``G`` -- a list of polynomials, which is a tropical Gröbner basis
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A minimal tropical Gröbner basis of the ideal generated by G. For now, it is not reduced


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z]
        sage: G = tentative_F5Trop(list1,R,w)
        sage: minimalization(G,R1,w)
        [x + y + z, y^2 + y*z + z^2, z^3]

    """
    P = G[0].parent()
    list_minimal = []
    listlm = []
    for g in G:
        listlm.append(LM_trop(g,R1,w))
    i = 0
    while i < len(listlm):
        j = 0
        flag = True
        #since we proceed degree by degree, it is enough to only consider j<i, as dividers are necessarily of degree stricly less
        a = listlm[i]
        while j < i:
            b = listlm[j]
            if (i != j) and P.monomial_divides(b,a):
                flag = False
                break
            j = j + 1
        if flag:
            list_minimal.append(G[i])
        i = i + 1
    return list_minimal

###############################################
#          Reduction of a Trop GB             #
###############################################

def symb_preprocess_for_reduction(G,R1,w):
    r"""
    Computes a list of 3-tuples (integer, integer, Macaulay matrix) , starting from G a minimal tropical GB, in order to do the inter-reduction of G

    
    INPUT:

    - ``G`` -- a list of polynomials, which is a minimal tropical Gröbner basis
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A list of couples 3-tuples (integer, integer, Macaulay matrix)
    The first one is the degree of the Macaulay matrix, the second is the number of rows coming from
    polynomials of G in the matrix (placed at the top)
    a list indexed by d of the list of monomials of degree d, for all d up to the maximal required


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: list1 = [z^2, x^3 + 2*x^2*y + y^3, y^3 + y^2*z + y*z^2, y^2*z + z^3]
        sage: G = tentative_F5Trop(list1,R,w)
        sage: G = minimalization(G,R1,w)
        sage: symb_preprocess_for reduction(G,R1,w)
        [[2, 1, [0 0 0 0 0 1]],
         [
              [1 2 0 1 0 0 0 0 0 0]
              [0 0 0 1 0 0 1 0 1 0]
              [0 0 0 0 0 0 1 0 0 1]
              [0 0 0 0 0 0 0 0 1 0]
         3, 3, [0 0 0 0 0 0 0 0 0 1]
        ]]
    """ 
    listmat = []
    lmonom = []
    d0=1
    P = G[0].parent()
    lmonom.append(list_monom_iteratif_no_repeat(P,d0,[]))


    deg = set([])
    for g in G:
        deg.add(g.degree())
    for d in deg:
        for u in range(d0+1,d+1):
            lmonom.append(list_monom_iteratif_no_repeat(P,u,lmonom[u-2]))
        d0=d

        listd = []
        list_to_do = []
        set_done_monom = set([])
        nb_pol_degd = 0
        for g in G:
            if g.degree()==d:
                listd.append(g)
                list_to_do.append(g)
                set_done_monom.add(LM_trop(g,R1,w))
                nb_pol_degd +=1
        while len(list_to_do)>0:
            poly = list_to_do.pop(0)
            set_monom = set(poly.monomials())-set_done_monom
            set_monom = list(set_monom)
             
            while len(set_monom)>0:
                mon = set_monom.pop()
                for g in G:
                    glm =  LM_trop(g,R1,w)
                    if parent(mon).monomial_divides(glm,mon):
                        mult = parent(mon).monomial_quotient(mon,glm)
                        listd.append(mult*g)
                        list_to_do.append(mult*g)
                        break
                set_done_monom.add(mon)
        mat = Make_Macaulay_Matrix_sparse_dict_no_sign(listd,lmonom[d-1])
        listmat.append([d,nb_pol_degd,mat])
    
    return listmat, lmonom

def tropical_complete_row_reduction(Mat, lmonom, R,w):
    r"""
    Computes a tropical complete row reduction of the Macaulay matrix Mat, with columns indexed by lmonom
    and with term order defined by R1 and w

    
    INPUT:

    - ``Mat`` -- a matrix
    - ``lmonom`` -- a list of the monomials of a given degree
    - ``R`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A reduced tropical Gröbner basis of the ideal generated by G. 


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: Mat = Matrix(QQ,5,10,[1,2,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1])
        sage: lmonom = list_monom(S,3)
        sage: u = tropical_complete_row_reduction(Mat, lmonom, R,w)
        sage: u[0]
        [1 0 0 0 0 0 0 0 2 0]
        [0 1 0 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0]
        sage: u[1]
        [x^3, y^3, y^2*z, y*z^2, z^3, x*y*z, x*y^2, x*z^2, x^2*y, x^2*z]
 
    """
    n = Mat.nrows()
    m = Mat.ncols()
    Mactilde = copy(Mat)
    list_monom = copy(lmonom)
    T = R.base_ring()
    nvar = len(w)
    index_column = 0
    for i in range(n):
        #Looking for the LT for the tropical order
        max_term = None
        max_col = None
        entry_non_zero = None
        for j0 in range(index_column,m):
            if max_term is None and Mactilde[i,j0] != 0:
                max_term = Mactilde[i,j0]
                max_col = j0
                max_monomial = list_monom[j0]
                max_trop_weight = T(Mactilde[i,j0]).valuation()+sum([list_monom[j0].degrees()[uu]*w[uu] for uu in range(nvar)])
            elif Mactilde[i,j0] != 0:
                term_trop_weight = T(Mactilde[i,j0]).valuation()+sum([list_monom[j0].degrees()[uu]*w[uu] for uu in range(nvar)])
                if term_trop_weight < max_trop_weight or (term_trop_weight == max_trop_weight and R(list_monom[j0])> R(max_monomial)):
                    max_term = Mactilde[i,j0]
                    max_col = j0
                    max_trop_weight = term_trop_weight
                    max_monomial = list_monom[j0]
        if max_term is not None:
            column_transposition(Mactilde, list_monom, index_column, max_col) 
            Mactilde.rescale_row(i, Mactilde[i,index_column]^(-1), start_col=index_column)
            for l in range(n):
                if l <> i:
                    Mactilde.add_multiple_of_row(l, i, -Mactilde[l,index_column] , start_col=index_column)
            index_column +=1
    return Mactilde, list_monom



def reduction_trop_GB(G,R1,w):
    r"""
    Computes a reduced tropical Gröbner basis of the ideal generated by the list of polynomials in G.
    The tropical term order is defined by the valuation on the base ring of R1, its ambiant
    monomial order, and the tropical weight w.
    G is a minimal tropical Gröbner basis (for the same term order).

    
    INPUT:

    - ``G`` -- a list of polynomials, which is a tropical Gröbner basis
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A reduced tropical Gröbner basis of the ideal generated by G. 


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: list1 = [z^2, x^3 + 2*x^2*y + y^3, y^3 + y^2*z + y*z^2, y^2*z + z^3]
        sage: G = tentative_F5Trop(list1,R,w)
        sage: G = minimalization(G,R,w)
        sage: G = reduction_trop_GB(G,R,w)
        [z^2, x^3 + 2*x^2*y, y^3, y^2*z]

    """
    listmat, lmonom = symb_preprocess_for_reduction(G,R1,w)
    P = G[0].parent()
    G2 = []



    for u in listmat:
        mat = u[2]
        d = u[0]
        nb = u[1]
        Mactilde, lmonom2 = tropical_complete_row_reduction(mat, lmonom[d-1], R1,w)
        for i in range(nb):
            a = Reconstruction(Mactilde,i, lmonom2,  [])
            G2.append(a[0])
    return G2

###############################################
#              Reduced Trop GB                #
###############################################

def reduced_trop_GB(F,R1,w):
    r"""
    Computes a reduced tropical Gröbner basis of the ideal generated by the list of polynomials in F.
    The tropical term order is defined by the valuation on the base ring of R1, its ambiant
    monomial order, and the tropical weight w.

    
    INPUT:

    - ``F`` -- a list of polynomials
    - ``R1`` -- a polynomial ring
    - ``w`` -- a tropical weight, i.e. an n-tuple of real numbers

    OUTPUT:

    A reduced tropical Gröbner basis of the ideal generated by F. 


    EXAMPLES::

        sage: from sage.rings.polynomial.padics.toy_F5 import all
        sage: S.<x,y,z>=QQ[]
        sage: R.<x1,y1,z1>=Qp(2)[]
        sage: w = [0,0,0]
        sage: list1 = [z^2, x^3 + 2*x^2*y + y^3, y^3 + y^2*z + y*z^2, y^2*z + z^3]
        sage: G = reduced_trop_GB(G,R,w)
        [z^2, x^3 + 2*x^2*y, y^3, y^2*z]

    """    
    G = tentative_F5Trop(F,R1,w)
    G = minimalization(G,R1,w)
    G = reduction_trop_GB(G,R1,w)
    return G    

def sym_pol(n):
    #Computes the list of elementary symmetric polynomials in n variables (sorted by increasing degree)
    Sym = SymmetricFunctions(QQ)
    e = Sym.elementary()
    list1 = []
    for i in range(1,n+1):
        list1.append(e[i].expand(n))
    return list1


