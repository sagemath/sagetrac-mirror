# -*- coding: utf-8 -*-
"""

Educational Versions of the F5 Algorithm to compute Groebner bases.

We have followed [Fau02]_, [AP11]_, [EF17]_ and [VY17]_.

No attempt was made to optimize the algorithm as the emphasis of
these implementations is a clean and easy presentation. To compute a
Groebner basis in Sage efficiently use the
:meth:`sage.rings.polynomial.multi_polynomial_ideal.MPolynomialIdeal.groebner_basis()`
method on multivariate polynomial objects.

The main algorithm is ``educational_F5``; secondary algorithms of importance are ``SymbolicPreprocessing``,  
``RowEchelonMac_and_listpivots``  and ``UpdateGpairs``. 

The idea of the algorithm is to use an adapted Buchberger termination criterion
on S-pairs while reducing all computations to linear algebra
and eliminating useless polynomials with the F5 elimination criterion.
The core part of the main algorithm is a while loop
while there are still S-polynomials, where it applies:

* ``SymbolicPreprocessing`` to construct a matrix with the polynomials of the pairs of a given degree d and reducers to reduce them. It uses the F5 elimination criterion.
* ``RowEchelonMac_and_listpivots`` to echelonize the previous matrix.
* ``UpdateGpairs`` to update the list in construction of a Groebner basis (more precisely, an ``\mathfrak{S}``-Groebner basis), and obtain the new S-polynomials.

A necessary tool in order to use the F5 elimination criterion is the notion of signature.
Let ``f_0,\dots,f_s`` be homogeneous polynomials (ordered with increasing degree) in ``k[X_1,\dots,X_n]`` for ``k`` a field. 
Let ``I= \left\lbrace f_0,\dots,f_s \right\rbrace`` and ``f \in I`` be
a homogeneous polynmomial.
We can define the signature of ``f`` as the following.
Among all the possible ways of writing ``f=\sum_{j=0}^i a_j f_j``(``a_j`` homogeneous, all the non-zero ``a_j f_j`` of same degree), we select one such that:

*Firstly, ``i`` is minimal such that ``a_i \neq 0.``
*Secundly, among all such writing for this minimal ``i,`` we select one with ``x^\alpha:=a_i.lm()`` minimal.

Then we write that the signature of ``f`` is ``x^\alpha e_i,`` with ``(e_0,\dots,e_s)`` being the canonical basis of ``k[X_1,\dots,X_n]^s.``

It is possible to define a vector space filtration of ``I`` using increasing signature as filtration.
The F5 elimination criterion enables the detection of some of the vector spaces in this filtration that do not provide any new polynomial
to the ideal.

In order to achieve this detection, all the polynomial handled by the program are encapsulated into a list of three elements:
* an integer ``i.``
* a monomial ``x^\alpha.`` 
* the polynomial itself, whose signature is ``x^\alpha e_i``.

Only operations on polynomials that preserve the knowledge of the signature are performed (``\textit{i.e.}`` addition of two polynomials of distinct signature).
A Gröbner basis that is compatible with the filtration by signature is called an ``\mathfrak{S}``-Groebner basis. The F5 algorithm computes
``\mathfrak{S}``-Groebner bases.


Let us look at the following example coming from Christian Eder's talk at ISSAC 2017 in Kaiserslautern (reference: http://www.mathematik.uni-kl.de/~ederc/download/issac-2017.pdf).

In ``\mathbb{Q}[x,y,z]`` with degree-reverse lexicographical ordering as monomial ordering, we begin with ``f_1 = xy-z^2`` and ``f_2 = y^2-z^2.``

The algorithm will first produce the list ``G=[[0,1,xy-z^2],[1,1,y^2-z^2]]`` as Gröbner basis in construction.
It encodes two polynomials with signature: ``g_1=f_1`` of signature ``e_0`` and ``g_2=f_2`` of signature ``e_1.``

A corresponding S-pair is then produced: ``[3,1,x,[1,1,y^2-z^2],y,[0,1,xy-z^2]].``
It summarizes the fact that the S-polynomial ``S_{g_1,g_2}`` of ``g_1`` and ``g_2`` is obtained as `` xg_2-yg_1.``

The SymbolicPreprocessing will then produce a Macaulay matrix with two rows (in that order):
*one for ``yg_1=xy^2-yz^2`` with corresponding signature ``ye_0.``
*one for ``xg_2=xy^2-x*z^2`` with corresponding signature ``xe_1.``


This matrix is then echelonize by RowEchelonMac_and_listpivots.
The first row is used as pivot to eliminate ``xy^2`` on the second row.
A new polynomial is produced on this second row: ``g_3=-xz^2+yz^2,`` of signature ``xe_1.``
As it produces a new leading monomial (``xz^2``) for the polynomials of signature less or equal to ``xe_1,`` it is added to ``G``
as ``[1,x,-xz^2+yz^2].``

UpdateGpairs then produce two new S-pairs corresponding to the S-polynomials ``S_{g_1,g_3}`` and ``S_{g_2,g_3},`` respectively:
*``[4,1,y,[1,x,-xz^2+yz^2],z^2,[0,1,xy-z^2]].``
*``[5,1,y^2,[1,x,-xz^2+yz^2],xz^2,[1,1,y^2-z^2]].``

We now go again through the main while loop of the algorithm.
However, this time, we will see the F5 elimination criterion into action.
Indeed, for the first S-pair, the signature that would be produced on the matrix in degree 4 for the left polynomial is
``xye_1.`` Since ``g_1.lm()=f_0.lm()=xy,`` ``xy`` can be achieved as leading monomial of a polynomial of signature ``x^\beta e_j`` with ``j<1.``
The F5 elimination criterion then states exactly that both left and right polynomials of the S-pair can be discarded.

As a consequence, no matrix is built in degree 4.

We can then go again through the main while loop of the algorithm.
The F5 elimination criterion applies again.
Indeed, the only remaining S-pair is ``[5,1,y^2,[1,x,-xz^2+yz^2],xz^2,[1,1,y^2-z^2]].``
The signature that would be produced on the matrix in degree 5 for the right polynomial is
``xy^2e_1.`` Since ``g_1.lm()=f_0.lm()=xy,`` ``xy^2`` can be achieved as leading monomial of a polynomial of signature ``x^\beta e_j`` with ``j<1.``
The F5 elimination criterion then states again that both left and right polynomials of the S-pair can be discarded.

As a consequence, no matrix is built in degree 5.
No S-pair is remaining and the F5 algorithm then outputs 
``[y^2-z^2,xy-z^2,-xz^2+yz^2]`` as Gröbner basis of ``I= \left\lbrace f_0,f_1 \right\rbrace`` (indeed, in this case, this basis is
already reduced).

We should remark that the S-pair in degree 5 would have been ruled out by Buchberger's product criterion for elimination.
It is however not the case of the S-pair in degree 4, which here does require the F5 elimination criterion. Buchberger's chain criterion
would not have applied... note::
    The setting is strictly that of homogeneous initial polynomials.

EXAMPLES:

Consider Katsura-6 w.r.t. a ``degrevlex`` ordering.::

    sage: from sage.rings.polynomial.toy_F5 import *
    sage: P.<a,b,c,e,f,g,h,i,j,k> = PolynomialRing(GF(32003),10)
    sage: I = sage.rings.ideal.Katsura(P,6).homogenize()
    sage: F = I.gens()

    sage: g1 = F5_reduced(F)
    sage: g2 = I.groebner_basis()

All algorithms actually compute a Groebner basis::

    sage: Ideal(g1).basis_is_groebner()
    True
    sage: Ideal(g2).basis_is_groebner()
    True

The results are correct::

    sage: Ideal(g1) == Ideal(g2)
    True

REFERENCES:

- [Fau2002]_
- [AP2011]_
- [EF2017]_
- [VY2017]_

AUTHOR:

- Tristan Vaccon (2016--2017)
"""




from copy import copy
from sage.structure.sequence import Sequence
from sage.matrix.matrix_space import MatrixSpace
from sage.matrix.constructor import matrix, Matrix

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
        sage: from sage.rings.polynomial.toy_F5 import index_last_variable
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


def list_monom_no_repeat(P,d,ld1):
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
        sage: from sage.rings.polynomial.toy_F5 import list_monom_no_repeat
        sage: list_monom_no_repeat(S,2,[x,y,z])
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
        sage: from sage.rings.polynomial.toy_F5 import Make_Macaulay_Matrix_sparse_dict
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
        sage: from sage.rings.polynomial.toy_F5 import F5_criterion
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
        sage: from sage.rings.polynomial.toy_F5 import Reducible_large_resp_sign
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
        sage: from sage.rings.polynomial.toy_F5 import*
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
    - ``pairs`` -- a list of pairs of polynomials with signature

    
    OUTPUT:

    A multiple of an element of G with same signature and leading monomial
    as g. It is obtained with the latest element of G possible

    EXAMPLES::

        sage: S.<x,y,z>=QQ[]
        sage: from sage.rings.polynomial.toy_F5 import*
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


def SymbolicPreprocessing(R,G, d, pairs, monom):
    r"""
    
    Constructs the Macaulay matrix from the pairs in pairs,
    while using the F5 criterion

    
    INPUT:

    - ``R`` -- a polynomial ring
    - ``G`` -- a list of polynomials with signature
    - ``d`` -- an integer, representing the degree of the processed polynomials
    - ``pairs`` -- a list of pairs of polynomials with signature
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    
    OUTPUT:

    A Macaulay matrix, as a couple of a matrix built from polynomials and
    a list of signatures for the rows of this matrix 

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: G = [[0,y,x**2+2*z**2], [1,x,x*y+y**2], [2,1,y*z+z**2], [3,1,z**2]]
        sage: pairs = []
        sage: pairs.append([3,1,x,[1,x,x*y+y**2],y, [0,1,x**2+2*z**2]])
        sage: pairs.append([3,3,y,[3,1,x*z+3*y*z+z**2],z, [3,1,x*y+2*y*z+z**2]])
        sage: d = 3
        sage: Mac = SymbolicPreprocessing(S, G, d, pairs,[x**3, x**2*y, x*y**2, y**3, x**2*z, x*y*z, y**2*z, x*z**2, y*z**2, z**3])
        sage: Mac[0][0]
        [0 0 0 0 0 0 0 0 1 1]
        [0 0 0 0 0 0 0 0 0 1]
        [0 0 0 0 0 0 0 0 1 0]
        sage: Mac[0][1]
        [[2, z], [3, z], [3, y]]

    """

    pairs_to_process = []
    pairs_to_remove = []
    list_already_erased = []
    list_already_added = []


    for pair in pairs:
        i_acceptable = True
        j_acceptable = True
        if pair[0] == d :
            pairs_to_remove.append(pair)
            signj = [pair[3][0],pair[3][1]*pair[2]]
            signi = [pair[5][0],pair[5][1]*pair[4]]
            if not(signj in list_already_erased) and not(signi in list_already_erased):
                if F5_criterion(signj,G):
                    j_acceptable = False
                    list_already_erased.append(signj)

                if F5_criterion(signi,G):
                    i_acceptable = False
                    list_already_erased.append(signi)

                if i_acceptable and j_acceptable:
                    pairs_to_process.append(pair)



    #Erasing the pairs that we have considered
    for pair in pairs_to_remove:
        pairs.remove(pair)



    #Producing the list of polynomials to be processed
    # We produce a set of the leading monomials



    poly_to_do = []
    list_poly = []
    set_done_monom = set()
    for pair in pairs_to_process:
        signj_polyj = [pair[3][0],pair[3][1]*pair[2],pair[2]*pair[3][2]]
        signj = [pair[3][0],pair[3][1]*pair[2]]
        if poly_to_do.count(signj_polyj)<1 and list_already_added.count(signj)<1:
            poly_to_do.append(signj_polyj)
            list_poly.append(signj_polyj)
            list_already_added.append(signj)
        
        
        signi_polyi = [pair[5][0],pair[5][1]*pair[4], pair[4]*pair[5][2]]
        signi = [pair[5][0],pair[5][1]*pair[4]]
        if poly_to_do.count(signi_polyi)<1 and list_already_added.count(signi)<1:
            poly_to_do.append(signi_polyi)
            list_poly.append(signi_polyi)
            list_already_added.append(signi)

    

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
            # We do take for mon the biggest monomial available here
            # We look for a g in G and x**alpha such that  mon ==x**alpha* g.lm() with smallest signature
            # In other word, we look for a "reducer." Yet for signature reason, this reducer could be on more reduced than reducing
            smallest_g = 0
            smallest_sign = []
            for g in G:
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
                new_poly = smallest_sign+[smallest_g]
                poly_to_do.append(new_poly)
                poly_to_do.sort()
                list_poly.append(new_poly)
                list_poly.sort()
                list_already_added.append(smallest_sign)
            set_done_monom.add(mon)


    if len(list_poly)>0:
        Mac = Make_Macaulay_Matrix_sparse_dict(list_poly, monom)
    else:
        Mac = [Matrix(R.base_ring(),0,0,[]),[]]


    return Mac, pairs


def list_indexes(Mat):
    r"""
    
    Computes the list of the indexes of the first non-zero coefficients
    of the rows of the matrix Mat

    
    INPUT:

    - ``Mat`` -- a matrix
    
    OUTPUT:

    a list of integers

    EXAMPLES::

        sage: from sage.rings.polynomial.toy_F5 import*
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
            if Mat[i,j] !=0:
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

        sage: from sage.rings.polynomial.toy_F5 import*
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
            if Mactilde[i,j] !=0:
                break
            j=j+1
        listpivots.append(j)
        if j<m:
            Mactilde.rescale_row(i, Mactilde[i,j]**(-1), start_col=j)
            for l in range(i+1,n):
                if Mactilde[l,j]!=0:
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

        sage: from sage.rings.polynomial.toy_F5 import*
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
    such that their signature is attained by another member of ``G`` with smaller
    leading monomial. This is done only for polynomials of degree d.
    

    
    INPUT:

    - ``G`` -- a list of polynomials with signature, ordered by signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    
    OUTPUT:

    - ``newG`` -- an updated list of polynomials with signature, with no redundancy of signature

    EXAMPLES:: 

        sage: from sage.rings.polynomial.toy_F5 import*
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


def UpdateGpairs(Mac, Mactilde, listpivots, G, pairs,d, redundancy_test, monom):
    r"""
    
    In the F5 algorithm, updates the current list of polynomials with
    signatures G, according to the Macaulay matrix Mac,
    its row-echelon form (with no choice of pivot) Mactilde,
    and previous computation.
    pairs, the list of the remaining S-pairs to proceed is also
    updated.
    G and pairs are updated in place.

    
    INPUT:

    - ``Mac`` -- a Macaulay matrix: a list of a matrix and a list of the signatures corresponding to its rows
    - ``Mactilde`` -- a Macaulay matrix (just the matrix)
    - ``listpivots`` -- the list of the pivots used when computing the row-echelon form Mactilde of Mac
    - ``G`` -- a list of polynomials with signature
    - ``pairs`` -- a list of S-pairs of polynomials with signature
    - ``d`` -- an integer, to tell the degree of the polynomials that are processed
    - ``redundancy_test`` -- a boolean to check whether a redundancy test is required
    - ``monom`` -- the list of the monomials of degree d, ordered decreasingly
    
    OUTPUT:

    - ``G`` -- an updated list of polynomials with signature
    - ``pairs`` -- an updated list of S-pairs of polynomials with signature

    EXAMPLES::  

        sage: from sage.rings.polynomial.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: mat = Matrix(QQ,2,6,[1,1,0,0,0,1,1,2,1,0,0,-1])
        sage: Mac = [mat,[[0,1],[1,1]]]
        sage: Mactilde = Matrix(QQ,2,6,[1,1,0,0,0,1,0,1,1,0,0,-2])
        sage: listpivots = [0,1]
        sage: G = [[0,1,x**2+x*y+z**2],[1,1,x**2+2*x*y+y**2-z**2]]
        sage: pairs = []
        sage: monom = [x**2, x*y, y**2, x*z, y*z, z**2]
        sage: d = 2
        sage: G2, pairs2=UpdateGpairs(Mac, Mactilde, listpivots, G, pairs, d,true, monom)
        sage: G2
        [[0, 1, x^2 + x*y + z^2], [1, 1, x*y + y^2 - 2*z^2]]
        sage: pairs2
        [[3, 1, x, [1, 1, x*y + y^2 - 2*z^2], y, [0, 1, x^2 + x*y + z^2]]]
    """
    n = Mac[0].nrows()
    m = Mac[0].ncols()

    if n == 0:
        return G,pairs


    listindexes = list_indexes(Mac[0])
    R = G[0][2].parent()

    newG = []
    for i in range(n):
        if listindexes[i]!=listpivots[i] and listpivots[i]<m:
            if not(Reducible_large_resp_sign_family([Mac[1][i][0],Mac[1][i][1],monom[listpivots[i]]],G)):
                newG.append(Reconstruction(Mactilde,i,monom,[Mac[1][i][0],Mac[1][i][1]]))




    G2 = G+newG
    G2.sort()
    if redundancy_test:
        pairs = []
        G2=Eliminate_unreduced(G2,d)
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
                    pairs.append([du,g2[0],R(u/g2[2].lm()),g2,R(u/g[2].lm()), g])
    
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
                    pairs.append([du,gtemp2[0],R(u/gtemp2[2].lm()),gtemp2,R(u/gtemp1[2].lm()), gtemp1])
    
        long = len(newG)
        for j1 in range(long):
            for j2 in range(j1+1,long):
                g = newG[j1]
                g2 = newG[j2]
                u = g[2].lm().lcm(g2[2].lm())
                du = u.degree()
                if du > d :
                    pairs.append([du,g2[0],R(u/g2[2].lm()),g2,R(u/g[2].lm()), g])



    if len(newG)>0 or redundancy_test:
        pairs.sort()

    return G2, pairs
    


def educational_F5(list1):
    r"""
    Computes a Gröbner bases of the ideal generated by the list of polynomials list1.
    The monomial order is the ambiant monomial order in the polynomial ring containing
    the polynomials of list1

    
    INPUT:

    - ``list1`` -- a list of polynomials

    OUTPUT:

    A Gröbner basis of the ideal generated by list1. For now, it is not reduced


    EXAMPLES::

        sage: from sage.rings.polynomial.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z]
        sage: G = educational_F5(list1)
        sage: G
        [x + y + z, y^2 + y*z + z^2, z^3]

    ALGORITHM:

    The following algorithm is adapted from [VY17]_.

    REFERENCES:

    .. [VY17] T. Vaccon, K.Yokoyama, A Tropical F5 Algorithm, proceedings of ISSAC 2017.

    """
    #Initialization
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
    pairs = []
    for i in range(s):
        for j in range(i+1,l):
            u = list[i].lm().lcm(list[j].lm())
            d = u.degree()
            pairs.append([d,j,R(u/list[j].lm()),listsign[j],R(u/list[i].lm()), listsign[i]])
    pairs.sort()

    d = 1
    monom = list_monom_no_repeat(R,d,[])
    while pairs != [] :
        Mac, pairs = SymbolicPreprocessing(R,G, d,pairs, monom)
        Mactilde, listpivots = RowEchelonMac_and_listpivots(Mac[0])
        redundancy_check = (list_degrees.count(d)>0)
        G, pairs = UpdateGpairs(Mac, Mactilde, listpivots, G, pairs,d,redundancy_check, monom)
        d = d+1
        monom = list_monom_no_repeat(R,d,monom)
    G = [g[2] for g in G]
    return G


def F5_reduced(list1):
    r"""
    Computes the reduced Gröbner bases of the ideal generated by the list of polynomials list1.
    The monomial order is the ambiant monomial order in the polynomial ring containing
    the polynomials of list1

    
    INPUT:

    - ``list1`` -- a list of polynomials

    OUTPUT:

    The reduced Gröbner basis of the ideal generated by list1.


    EXAMPLES::

        sage: from sage.rings.polynomial.toy_F5 import*
        sage: S.<x,y,z>=QQ[]
        sage: list1 = [x+y+z, x*y+x*z+y*z,x*y*z, z**4]
        sage: G = F5_reduced(list1)
        sage: G
        [z^3, y^2 + y*z + z^2, x + y + z]

    """
    G = educational_F5(list1)
    G2 = Sequence(G)
    G3 = G2.reduced()
    return G3
