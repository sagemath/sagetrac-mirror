r"""
Computes the connectivity between (sets of) vertices in a graph given the 
probability matrix for that graph

The main function is "connectivity_probability(...)". It requires the input 
of a (symmetric) probability matrix for a graph. A probability matrix is a matrix
whose rows and columns represent the vertices of a graph and whose entries are 
the connectivity probability between those vertices. User may optionally input a 
subset of vertices for which to compute the connectivity, and/or a subset of vertices 
to compute the probability that all of the vertices in the subset are connected. 
The main output is a connectivity matrix for the graph given represented by 
the given probability matrix.  

AUTHORS:

- James Rudzinski (2015-12-31): initial version

EXAMPLES::


    This example shows connectivity_probability with input of a real-valued probability matrix. 
    The output matrix has the same type as the input matrix. Either value type (real or rational)
    can be used with any selection of optional inputs.
    
    ::
  
        sage: S
        [ 1.00000000000000 0.000000000000000 0.500000000000000 0.333333333333333]
        [0.000000000000000  1.00000000000000  1.00000000000000 0.750000000000000]
        [0.500000000000000  1.00000000000000  1.00000000000000 0.000000000000000]
        [0.333333333333333 0.750000000000000 0.000000000000000  1.00000000000000]
        sage: connectivity_probability(S)
        [ 1.00000000000000 0.625000000000000 0.625000000000000 0.583333333333333]
        [0.625000000000000  1.00000000000000  1.00000000000000 0.791666666666667]
        [0.625000000000000  1.00000000000000  1.00000000000000 0.791666666666667]
        [0.583333333333333 0.791666666666667 0.791666666666667  1.00000000000000]

    This example shows connectivity_probability run with a rational-valued probability matrix.

    ::

        sage: T
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: connectivity_probability(T)
        [    1   5/8   5/8  7/12]
        [  5/8     1     1 19/24]
        [  5/8     1     1 19/24]
        [ 7/12 19/24 19/24     1]

    We now show connectivity_probability with the probability matrix as above and the 
    vertex list containing the vertices 0 and 2. Only the submatrix for the 
    vertices 0 and 2 is calculated and output. The purpose of this is to speed up
    computation of the desired submatrix by only computing those values.

    ::

        sage: T
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: TE
        [0, 2]
        sage: connectivity_probability(T,TE)
        [  1 5/8]
        [5/8   1]

    This example shows the function with both optional inputs. 
    The output matrix is the submatrix for the connectivity for the vertices 
    0 and 2 and the second output is the probability that the vertices in 
    grouplist (vertices 0, 1, 2, and 3) are simultaneously all connected.
    Note that this simutaneous connectivity probability is not obtainable
    from the first output probability matrix, which only contains pairwise 
    connectivity probability information.

    ::

        sage: T,TE,TF
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 2], [0, 1, 2, 3]
        )
        sage: connectivity_probability(T,TE,TF)
        (
        [  1 5/8]
        [5/8   1], 13/24
        )

    This example shows connectivity_probability with the matrix T from above 
    and the optional input for grouplist without the input for vertexlist. 
    The output is the entire connectivity matrix and the probability that 
    the vertices 0,1,2, and 3 are all connected. Again, note that the 
    simultaneous connectivity probability is not obtainable from the first
    matrix of pairwise connectivity probabilities.

    ::
 
        sage: T,TF
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 1, 2, 3]
        )
        sage: connectivity_probability(T,grouplist=TF)
        (
        [    1   5/8   5/8  7/12]
        [  5/8     1     1 19/24]
        [  5/8     1     1 19/24]
        [ 7/12 19/24 19/24     1], 13/24
        )

"""

#*****************************************************************************
#       Copyright (C) 2015 James Rudzinski <jerudzin@uncg.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#internal function to update the current connectivity matrix, M, 
#with the added edge e
def _addedge(M,e):	
    if M[e[0]][e[1]] != 1:
        M[e[0]] = M[e[0]] + M[e[1]]
        for i in M.nonzero_positions_in_row(e[0]):
            M[i] = M[e[0]]
    return M

#updates a connectivity matrix, M, with the deleted edge e, 
#the graph G is the graph of the connectivity matrix with e already deleted
def _deleteedge(M,G,e):	#function to delete an edge from a connectivity matrix
    N = M
    L0 = G.connected_component_containing_vertex(e[0]) #connected component of e[0]
    if e[1] not in L0:
        L1 = [i for i in N.nonzero_positions_in_row(e[0]) if i not in L0]
        for i in L1:
            N[e[0],i] = 0
        for i in L0:
            N[i] = N[e[0]]
            N[e[1],i] = 0
        for i in L1:
            N[i] = N[e[1]]
    return N

#The function _probconn is called if only the required probability matrix is passed to 
#connectivity_probability. The algorithm used for computing the matrix is described in 
#the documentation for the connectivity_probability function.
def _probconn(P):
    
    from sage.combinat.gray_codes import product
    	
    n = P.nrows()    # n is the number of vertices of G
    
    ones = ones_matrix(n)
    
    Q = ones - P
    
    a = P.nonzero_positions()
    b = Q.nonzero_positions()

    A = copy(a)
    B = copy(b)

    for x in a:
        if x[0] >= x[1]:
            A.remove(x)

    for x in b:
        if x[0] >= x[1]:
            B.remove(x)

    E = sorted(list(set(A).intersection(B)))

    missing_edges = sorted(list(set(B)-set(E)))
    permanent_edges = sorted(list(set(A)-set(E)))
    
    k = len(E)    # k is the number of edges to percolate

    	
    l = [0 for i in [1..k]]    # l represents the percolated edges
        # when l[i] is 0 then E[i] will be missing
        # when l[i] is 1 then E[i] will be present
    	
    # H represents each individual "instance" of the percolation
    # initially H will include all vertices of G,
    # but only edges of G not included in E (those are permanent edges)
    
    h = zero_matrix(n,n)

    for x in permanent_edges:
        h[x[0],x[1]] = 1
        h[x[1],x[0]] = 1

    H = Graph(h)    	
    	
    # C is a "connectivity matrix" for H
    # initializing C
    C = matrix.identity(n)
    for e in permanent_edges:
        C = _addedge(C,e)
    
    d = 1
    for i in [0..k-1]:
        d = d*Q[E[i][0],E[i][1]]
    	
    # S is a summation matrix, keeping a total count of instances with 
    # connectivity between all pairs of vertices
    # initializing S
    S = d*C
    	
    for p,i in product([2 for i in [1..k]]):
        l[p] += i
        if i == 1:
            Ep = E[p]
            H.add_edge(Ep)
            C = _addedge(C,Ep)
            d = d * (P[Ep[0],Ep[1]]/Q[Ep[0],Ep[1]])
            S = S + d*C

        else:
            Ep = E[p]
            H.delete_edge(Ep)
            C = _deleteedge(C,H,Ep)
            d = d * (Q[Ep[0],Ep[1]]/P[Ep[0],Ep[1]])
            S = S + d*C
    
    return S

#The function _probconn2 is called if a vertexlist is passed to 
#connectivity_probability but a grouplist is not.
def _probconn2(P,F):
    
    from sage.combinat.gray_codes import product
    	
    n = P.nrows()    # n is the number of vertices of G
    
    ones = ones_matrix(n)
    
    Q = ones - P
    
    a = P.nonzero_positions()
    b = Q.nonzero_positions()

    A = copy(a)
    B = copy(b)

    for x in a:
        if x[0] >= x[1]:
            A.remove(x)

    for x in b:
        if x[0] >= x[1]:
            B.remove(x)

    E = sorted(list(set(A).intersection(B)))

    missing_edges = sorted(list(set(B)-set(E)))
    permanent_edges = sorted(list(set(A)-set(E)))
    
    k = len(E)    # k is the number of edges to percolate    	
    	
    l = [0 for i in [1..k]]    # l represents the percolated edges
        # when l[i] is 0 then E[i] will be missing
        # when l[i] is 1 then E[i] will be present
    	
    # H represents each individual "instance" of the percolation
    # initially H will include all vertices of G,
    # but only edges of G not included in E (those are permanent edges)
    
    h = zero_matrix(n,n)

    for x in permanent_edges:
        h[x[0],x[1]] = 1
        h[x[1],x[0]] = 1

    H = Graph(h)    	
    	
    # C is a "connectivity matrix" for H
    # initializing C
    C = matrix.identity(n)
    for e in permanent_edges:
        C = _addedge(C,e)
    
    d = 1

    for i in [0..k-1]:
        d = d*Q[E[i][0],E[i][1]]
    	
    # S is a summation matrix, keeping a total count of instances with 
    # connectivity between all pairs of vertices
    # initializing S
    dC = d*C
    
    f = len(F)-1
    S = zero_matrix(dC[0,0].parent(),f+1)
    for i in [0..f]:
        for j in [0..f]:
            S[i,j] = dC[F[i],F[j]]
    	
    for p,i in product([2 for i in [1..k]]):
        l[p] += i
        if i == 1:
            Ep = E[p]
            H.add_edge(Ep)
            C = _addedge(C,Ep)
            d = d * (P[Ep[0],Ep[1]]/Q[Ep[0],Ep[1]])
            for m in [0..f]:
                for j in [0..f]:
                    S[m,j] = S[m,j] + d*C[F[m],F[j]]

        else:
            Ep = E[p]
            H.delete_edge(Ep)
            C = _deleteedge(C,H,Ep)
            d = d * (Q[Ep[0],Ep[1]]/P[Ep[0],Ep[1]])
            for m in [0..f]:
                for j in [0..f]:
                    S[m,j] = S[m,j] + d*C[F[m],F[j]]
    
    return S


#The function _probconn3 is called when a grouplist (with or without a vertexlist)
# is passed to connectivity_probability. If no vertexlist is sent, the vertexlist 
# defaults to the entire set of vertices.
def _probconn3(P,F1,F2):
    
    from sage.combinat.gray_codes import product
    	
    n = P.nrows()    # n is the number of vertices of G

    ones = ones_matrix(n)
    
    Q = ones - P
    
    a = P.nonzero_positions()
    b = Q.nonzero_positions()

    A = copy(a)
    B = copy(b)

    for x in a:
        if x[0] >= x[1]:
            A.remove(x)

    for x in b:
        if x[0] >= x[1]:
            B.remove(x)

    E = sorted(list(set(A).intersection(B)))

    missing_edges = sorted(list(set(B)-set(E)))
    permanent_edges = sorted(list(set(A)-set(E)))
    
    k = len(E)    # k is the number of edges to percolate
    	
    l = [0 for i in [1..k]]    # l represents the percolated edges
        # when l[i] is 0 then E[i] will be missing
        # when l[i] is 1 then E[i] will be present
    	
    # H represents each individual "instance" of the percolation
    # initially H will include all vertices of G,
    # but only edges of G not included in E (those are permanent edges)

    h = zero_matrix(n,n)

    for x in permanent_edges:
        h[x[0],x[1]] = 1
        h[x[1],x[0]] = 1

    H = Graph(h)
    	
    # C is a "connectivity matrix" for H
    # initializing C
    C = matrix.identity(n)
    for e in permanent_edges:
        C = _addedge(C,e)
    
    d = 1

    for i in [0..k-1]:
        d = d*Q[E[i][0],E[i][1]]
    	
    # S is a summation matrix, keeping a total count of instances with 
    # connectivity between all pairs of vertices
    # initializing S
    dC = d*C
    
    f = len(F1)-1
    S = zero_matrix(dC[0,0].parent(),f+1)
    for i in [0..f]:
        for j in [0..f]:
            S[i,j] = dC[F1[i],F1[j]]
    T = 0

    z = 1
    for y in [0..len(F2)-1]:
        z = z*C[F2[0],F2[y]]
    if z == 1:
        T = d+T
  
    	
    for p,i in product([2 for i in [1..k]]):
        l[p] += i
        if i == 1:
            Ep = E[p]
            H.add_edge(Ep)
            C = _addedge(C,Ep)
            d = d * (P[Ep[0],Ep[1]]/Q[Ep[0],Ep[1]])
            for m in [0..f]:
                for j in [0..f]:
                    S[m,j] = S[m,j] + d*C[F1[m],F1[j]]
            z = 1
            for y in [0..len(F2)-1]:
                z = z*C[F2[0],F2[y]]
            if z == 1:
                T = d+T
        else:
            Ep = E[p]
            H.delete_edge(Ep)
            C = _deleteedge(C,H,Ep)
            d = d * (Q[Ep[0],Ep[1]]/P[Ep[0],Ep[1]])
            for m in [0..f]:
                for j in [0..f]:
                    S[m,j] = S[m,j] + d*C[F1[m],F1[j]]
            z = 1
            for y in [0..len(F2)-1]:
                z = z*C[F2[0],F2[y]]
            if z == 1:
                T = d+T

    return S, T


def connectivity_probability(probmatrix, vertexlist=[], grouplist=[]):
    r"""
    Returns the connectivity matrix for the graph. If vertexlist is present 
    returns the connectivity matrix for only the vertices in vertexlist.
    If grouplist is present, also returns the probability that all of the vertices
    in grouplist are simultaneously connected. 

    INPUT:

    - ``probmatrix`` -- matrix with entries from a subfield of the reals. 
      The probability matrix is an n x n matrix where n is the number of 
      vertices in the graph. Each entry in the matrix is the probability       
      that the corresponding edge is present in the graph (1 for always 
      present and 0 for never present). Directed and multi-graphs are not 
      supported at this time. Consequently, the probability matrix must be 
      symmetric and have ones on the diagonal. 
      
    - ``vertexlist`` -- list of integers (vertices) (default = []) 
      Optional input. If this parameter is entered, the function 
      sums and returns a subset of the connectivity matrix representing 
      only the vertices present in vertexlist

    - ``grouplist`` -- list of integers (vertices) (default = []) 
      Optional input. If this parameter is entered, the function 
      additionally returns the probability that all the vertices in
      grouplist are simultaneously connected in the graph. This information
      is not obtainable from the primary connectivity matrix output which 
      only returns pairwise connectivity information.

    OUTPUT:

    The connectivity matrix for the graph corresponding to the 
    probability matrix sent to the function. The connectivity matrix 
    represents the probability that the two index vertices of each entry
    are connected.
    If vertexlist is nonempty, only the submatrix of the connectivity 
    matrix represented by the vertices in vertexlist is fully computed 
    and returned.
    If grouplist is nonempty, the function also returns as the second output, 
    the probability that all the vertices in grouplist are simultaneously connected.

    .. SEEALSO::

        :func:`_addedge`
        :func:`_deleteedge`
        :func:`_probconn1`
        :func:`_probconn2`
        :func:`_probconn3`

    #########################################################################################


    EXAMPLES:


    This example shows connectivity_probability with input of a real-valued probability matrix. 
    The output matrix has the same type as the input matrix. Either value type (real or rational)
    can be used with any selection of optional inputs.
    
    ::
  
        sage: S
        [ 1.00000000000000 0.000000000000000 0.500000000000000 0.333333333333333]
        [0.000000000000000  1.00000000000000  1.00000000000000 0.750000000000000]
        [0.500000000000000  1.00000000000000  1.00000000000000 0.000000000000000]
        [0.333333333333333 0.750000000000000 0.000000000000000  1.00000000000000]
        sage: connectivity_probability(S)
        [ 1.00000000000000 0.625000000000000 0.625000000000000 0.583333333333333]
        [0.625000000000000  1.00000000000000  1.00000000000000 0.791666666666667]
        [0.625000000000000  1.00000000000000  1.00000000000000 0.791666666666667]
        [0.583333333333333 0.791666666666667 0.791666666666667  1.00000000000000]


    This example shows connectivity_probability run with a rational-valued probability matrix.

    ::

        sage: T
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: connectivity_probability(T)
        [    1   5/8   5/8  7/12]
        [  5/8     1     1 19/24]
        [  5/8     1     1 19/24]
        [ 7/12 19/24 19/24     1]


    We now show connectivity_probability with the probability matrix as above and the 
    vertex list containing the vertices 0 and 2. Only the submatrix for the 
    vertices 0 and 2 is calculated and output. The purpose of this is to speed up
    computation of the desired submatrix by only computing those values.


    ::

        sage: T
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: TE
        [0, 2]
        sage: connectivity_probability(T,TE)
        [  1 5/8]
        [5/8   1]


    This example shows the function with both optional inputs. 
    The output matrix is the submatrix for the connectivity for the vertices 0 and 2 
    and the second output is the probability that the vertices in 
    grouplist (vertices 0, 1, 2, and 3) are simultaneously connected.
    Note that this simutaneous connectivity probability is not obtainable
    from the first output probability matrix, which only contains pairwise 
    connectivity probability information.


    ::

        sage: T,TE,TF
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 2], [0, 1, 2, 3]
        )
        sage: connectivity_probability(T,TE,TF)
        (
        [  1 5/8]
        [5/8   1], 13/24
        )


    This example shows connectivity_probability with the matrix T from above 
    and the optional input for grouplist without the input for vertexlist. 
    The output is the entire connectivity matrix and the probability that 
    the vertices 0,1,2, and 3 are all connected. Again, note that the 
    simultaneous connectivity probability is not obtainable from the first
    matrix of pairwise connectivity probabilities.


    ::
 
        sage: T,TF
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 1, 2, 3]
        )
        sage: connectivity_probability(T,grouplist=TF)
        (
        [    1   5/8   5/8  7/12]
        [  5/8     1     1 19/24]
        [  5/8     1     1 19/24]
        [ 7/12 19/24 19/24     1], 13/24
        )


    
    #####################################################################################

    ERRORS:


    It is an error to input a matrix which is not a probability matrix. This means the matrix 
    must be a real valued (can be rational or integer) symmetric matrix with all entries 
    between 0 and 1 and one's along the diagonal. Each (i,j) entry of the matrix should
    represent the probability of the edge between vertex i and j existing.
    


    Our first error example shows the output when the input matrix is missing a one on the diagonal.

    ::
   
    
        TA = Matrix([[0,0,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])

        sage: TA
        [  0   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: connectivity_probability(TA)
        ---------------------------------------------------------------------------
        ValueError                                Traceback (most recent call last)
        <ipython-input-4-08109acec4d7> in <module>()
        ----> 1 connectivity_probability(TA)

        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/17930/connectivityprobability2.sageGBrxsq.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            703                 "must be in the range from 0 to 1. ")
            704         if not (probmatrix[i,i] == _sage_const_1 ):
        --> 705             raise ValueError("The diagonal entries of the probability matrix (probmatrix) "
            706             "must have a value of 1. ")
            707

        ValueError: The diagonal entries of the probability matrix (probmatrix) must have a value of 1.



    Our next error example shows the output when the input matrix is not symmetric

    ::


        TB = Matrix([[1,0,1/2,1/3],[0,1,1,3/4],[0,0,1,0],[1/5,1/4,0,1]])

        sage: TB
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [  0   0   1   0]
        [1/5 1/4   0   1]
        sage: connectivity_probability(TB)
        ---------------------------------------------------------------------------
        ValueError                                Traceback (most recent call last)
        <ipython-input-6-a27d36fc56e0> in <module>()
        ----> 1 connectivity_probability(TB)

        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/17930/connectivityprobability2.sageGBrxsq.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            707
            708     if not probmatrix.is_symmetric():
        --> 709         raise ValueError("The probability matrix (probmatrix) must be symmetric.")
            710
            711     #### end error checking

        ValueError: The probability matrix (probmatrix) must be symmetric.



    The next error example shows the output when the matrix entries are not from a subfield of the reals

    ::


        TC = Matrix([[i,i,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])  

        sage: TC
        [  I   I 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1]
        sage: connectivity_probability(TC)
        ---------------------------------------------------------------------------
        TypeError                                 Traceback (most recent call last)
        <ipython-input-8-7577613d86d8> in <module>()
        ----> 1 connectivity_probability(TC)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/17930/connectivityprobability2.sageGBrxsq.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            697         for j in range(n):
            698             if not (probmatrix[i,j] in RR):
        --> 699                 raise TypeError("The entries of the probability matrix (probmatrix) "
            700                 "must be real numbers. ")
            701             if not (_sage_const_0  <= probmatrix[i,j] <= _sage_const_1 ):

        TypeError: The entries of the probability matrix (probmatrix) must be real numbers.
    
    

    This error example shows the output when the matrix entries are not probabilities (between 0 and 1)

    ::

        TD = Matrix([[1,0,2,3],[0,1,1,4],[2,1,1,0],[3,4,0,1]])

        sage: TD
        [1 0 2 3]
        [0 1 1 4]
        [2 1 1 0]
        [3 4 0 1]
        sage: connectivity_probability(TD)
        ---------------------------------------------------------------------------
        ValueError                                Traceback (most recent call last)
        <ipython-input-10-807e430090f9> in <module>()
        ----> 1 connectivity_probability(TD)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/17930/connectivityprobability2.sageGBrxsq.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            700                 "must be real numbers. ")
            701             if not (_sage_const_0  <= probmatrix[i,j] <= _sage_const_1 ):
        --> 702                 raise ValueError("The entries of the probability matrix (probmatrix) "
            703                 "must be in the range from 0 to 1. ")
            704         if not (probmatrix[i,i] == _sage_const_1 ):
    
        ValueError: The entries of the probability matrix (probmatrix) must be in the range from 0 to 1.
    
    
    
    It will also cause an error if either one of the optional vertex lists are not lists of indices.
    This means that the optional vertex lists must be lists of integers between 0 and the number 
    of rows (or columns) of the input matrix.



    This error example shows the output when the optional vertex list is outside the range of the vertices

    ::

        TE = Matrix([[1,0,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])
        TE1 = [0,5]
     
        sage: TE,TE1
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 5]
        )
        sage: connectivity_probability(TE,TE1)
        ---------------------------------------------------------------------------
        IndexError                                Traceback (most recent call last)
        <ipython-input-33-2e2b1988c128> in <module>()
        ----> 1 connectivity_probability(TE,TE1)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            398             return probconn(probmatrix)
            399         else:
        --> 400             return probconn2(probmatrix, vertexlist)
            401     else:
            402         if not vertexlist:
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        probconn2(P, F)
            258     for i in (ellipsis_range(_sage_const_0 ,Ellipsis,f)):
            259         for j in (ellipsis_range(_sage_const_0 ,Ellipsis,f)):
        --> 260             S[i,j] = dC[F[i],F[j]]
            261
            262     for p,i in product([_sage_const_2  for i in (ellipsis_range(_sage_const_1 ,Ellipsis,k))]):
    
        /afs/uncg.edu/@sys/local/sage-6.6/src/sage/matrix/matrix0.pyx in 
        sage.matrix.matrix0.Matrix.__getitem__ (build/cythonized/sage/matrix/matrix0.c:5996)()
            918                     col += ncols
            919                 if col < 0 or col >= ncols:
        --> 920                     raise IndexError("matrix index out of range")
            921                 single_col = 1
            922
    
        IndexError: matrix index out of range
    

   
    This error example shows the output when the optional vertex list is not a list of integers.

    ::

        TF = Matrix([[1,0,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])
        TF1 = [0,2.1]
    
        sage: TF,TF1
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 2.10000000000000]
        )
        sage: connectivity_probability(TF,TF1)
        ---------------------------------------------------------------------------
        TypeError                                 Traceback (most recent call last)
        <ipython-input-37-c93b362cdbc2> in <module>()
        ----> 1 connectivity_probability(TF,TF1)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            398             return probconn(probmatrix)
            399         else:
        --> 400             return probconn2(probmatrix, vertexlist)
            401     else:
            402         if not vertexlist:
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        probconn2(P, F)
            258     for i in (ellipsis_range(_sage_const_0 ,Ellipsis,f)):
            259         for j in (ellipsis_range(_sage_const_0 ,Ellipsis,f)):
        --> 260             S[i,j] = dC[F[i],F[j]]
            261
            262     for p,i in product([_sage_const_2  for i in (ellipsis_range(_sage_const_1 ,Ellipsis,k))]):

        /afs/uncg.edu/@sys/local/sage-6.6/src/sage/matrix/matrix0.pyx in 
        sage.matrix.matrix0.Matrix.__getitem__ (build/cythonized/sage/matrix/matrix0.c:5932)()
            913             else:
            914                 if not PyIndex_Check(col_index):
        --> 915                     raise TypeError("index must be an integer")
            916                 col = col_index
            917                 if col < 0:

        TypeError: index must be an integer
    


    This error example shows the output when the optional group list is outside the range of vertices.

    ::

        TG = Matrix([[1,0,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])
        TG1 = [0,1,2,5]

        sage: TG,TG1
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 1, 2, 5]
        )
        sage: connectivity_probability(TG,grouplist=TG1)
        ---------------------------------------------------------------------------
        IndexError                                Traceback (most recent call last)
        <ipython-input-39-87f7562b9ee7> in <module>()
        ----> 1 connectivity_probability(TG,grouplist=TG1)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            402         if not vertexlist:
            403             vertexlist = (ellipsis_range(_sage_const_0 ,Ellipsis,probmatrix.nrows() - 
                            _sage_const_1 ))
        --> 404             return probconn3(probmatrix, vertexlist, grouplist)
            405         else:
            406             return probconn3(probmatrix, vertexlist, grouplist)
    
        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        probconn3(P, F1, F2)
            356     z = _sage_const_1
            357     for y in (ellipsis_range(_sage_const_0 ,Ellipsis,len(F2)-_sage_const_1 )):
        --> 358         z = z*C[F2[_sage_const_0 ],F2[y]]
            359     if z == _sage_const_1 :
            360         T = d+T
    
        /afs/uncg.edu/@sys/local/sage-6.6/src/sage/matrix/matrix0.pyx in 
        sage.matrix.matrix0.Matrix.__getitem__ (build/cythonized/sage/matrix/matrix0.c:5996)()
            918                     col += ncols
            919                 if col < 0 or col >= ncols:
        --> 920                     raise IndexError("matrix index out of range")
            921                 single_col = 1
            922
    
        IndexError: matrix index out of range



    This error example shows the ouput when the optional group list is not a list of integers.

    ::

        TH = Matrix([[1,0,1/2,1/3],[0,1,1,3/4],[1/2,1,1,0],[1/3,3/4,0,1]])
        TH1 = [0,1.0,2,3]

        sage: TH,TH1
        (
        [  1   0 1/2 1/3]
        [  0   1   1 3/4]
        [1/2   1   1   0]
        [1/3 3/4   0   1], [0, 1.00000000000000, 2, 3]
        )
        sage: connectivity_probability(TH,grouplist=TH1)
        ---------------------------------------------------------------------------
        TypeError                                 Traceback (most recent call last)
        <ipython-input-41-eae3b5a2d442> in <module>()
        ----> 1 connectivity_probability(TH,grouplist=TH1)

        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        connectivity_probability(probmatrix, vertexlist, grouplist)
            402         if not vertexlist:
            403             vertexlist = (ellipsis_range(_sage_const_0 ,Ellipsis,probmatrix.nrows() - 
                            _sage_const_1 ))
        --> 404             return probconn3(probmatrix, vertexlist, grouplist)
            405         else:
            406             return probconn3(probmatrix, vertexlist, grouplist)

        /afs/uncg.edu/user/j/jerudzin/.sage/temp/prdlinux1/27248/connectivityprobability.sageOCdU9W.py in 
        probconn3(P, F1, F2)
            356     z = _sage_const_1
            357     for y in (ellipsis_range(_sage_const_0 ,Ellipsis,len(F2)-_sage_const_1 )):
        --> 358         z = z*C[F2[_sage_const_0 ],F2[y]]
            359     if z == _sage_const_1 :
            360         T = d+T

        /afs/uncg.edu/@sys/local/sage-6.6/src/sage/matrix/matrix0.pyx in 
        sage.matrix.matrix0.Matrix.__getitem__ (build/cythonized/sage/matrix/matrix0.c:5932)()
            913             else:
            914                 if not PyIndex_Check(col_index):
        --> 915                     raise TypeError("index must be an integer")
            916                 col = col_index
            917                 if col < 0:

        TypeError: index must be an integer

    ############################################################################################


    .. NOTE::

        This function calls functions (_probconn1, _probconn2, _probconn3) 
        that step through the instances of the graph.

        The called functions use gray code (imported from sage.combinat.gray_codes) to 
        "step" through the iterations of a graph by changing only one edge each time.
        The probability is calculated by computing a matrix for each possible graph of 
        the probability matrix. The computed matrix is similar to an adjancency matrix 
        for a graph; it has a 1 in each entry where two vertices are contained in 
        the same connected component and a 0 if the two vertices are not in the same 
        component. This matrix is scaled by the probability that the current graph 
        would be obtained if edges were added (or deleted) from the original graph 
        according to the probabilities given by the probability matrix (``probmatrix``).
        The sum of all the scaled matrices from the possible graphs is the 
        output of the function (if ``vertexlist`` is empty). 

    ..  WARNING::
        
        This code will run for a very long time for graphs of size about 20 or larger. The
        algorithm time is exponential.  
    ...
    """

    


    #n is the number of rows (and columns) of the input matrix
    n = probmatrix.nrows()


    #### begin error checking

    for i in range(n):
        for j in range(n):
            if not (probmatrix[i,j] in RR):
                raise TypeError("The entries of the probability matrix (probmatrix) " 
                "must be real numbers. ")
            if not (0 <= probmatrix[i,j] <= 1):
                raise ValueError("The entries of the probability matrix (probmatrix) " 
                "must be in the range from 0 to 1. ")
        if not (probmatrix[i,i] == 1):
            raise ValueError("The diagonal entries of the probability matrix (probmatrix) " 
            "must have a value of 1. ")

    if not probmatrix.is_symmetric():
        raise ValueError("The probability matrix (probmatrix) must be symmetric.")

    #### end error checking
    ####
    #### At this point, matrices should be probability matrices, 
    #### that is, real-valued symmetric matrices with values 
    #### between 0 and 1 and all one's on the diagonal.
    


    if not grouplist:
        if not vertexlist:
            return _probconn(probmatrix)
        else:
            return _probconn2(probmatrix, vertexlist)
    else:
        if not vertexlist:
            vertexlist = [0..probmatrix.nrows() - 1]
            return _probconn3(probmatrix, vertexlist, grouplist)
        else:
            return _probconn3(probmatrix, vertexlist, grouplist)


    return

