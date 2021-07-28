"""
Stoll's algorithm for the reduction of point clusters in projective space.

A point cluster is a formal sum of points in projective space.

Given a stable point cluster, there exists a unique representative
that minimizes the covariant that Stoll defines in [Sto2009].
The ReduceCluster function computes an approximation to this
representative using steepest descent. The function below was provided
by Stoll.

AUTHORS:

- Alexander Galarraga (2021): initial implementation for Sage.
"""

from sage.matrix.constructor import matrix

from sage.modules.free_module_element import vector

from sage.rings.complex_mpfr import ComplexField
from sage.rings.real_mpfr import RealField


def ReduceCluster(S, eps=10**-6, c=1, precision=500):
    r"""
    Return an approximation to the minimal representative of a
    stable point cluster ``S``.

    A representative of a point cluster is minimal if it minimizes
    the covariant Stoll defines in [Sto2009].

    Note that this function will fail if the point cluster is not stable,
    and additionally it does NOT check if the point cluster ``S`` is not stable.

    This implementation is due to Stoll and based off of his implementation
    for Magma.

    INPUT:

    - ``S`` -- A list of points of projective space used to represent
      a point cluster.

    - ``eps`` -- (default: `10^{-6}`) The precision after which to treat
      a number as zero.

    - ``c`` -- (default: 1) The step size to use for steepest descent.

    - ``precision``-- (default: 500) The precision to use for the real and
      complex fields when approximating.

    OUTPUT:

    A tuple (``S``, ``mat1``, ``mat2``) where ``S`` is the approximated
    minimal representative of the point cluster and ``mat1`` is the
    matrix used to move the point cluster to its minimal representative.
    ``mat2`` is the inverse of ``mat1``.

    EXAMPLES::

        sage: from sage.schemes.projective.reduce_cluster import ReduceCluster
        sage: pnts = [P(1,0,0), P(0,1,0), P(0,0,1), P(1,1,1)]
        sage: ReduceCluster(pnts)
        (
                                                        [ 0  1  0]  [1 1 1]
                                                        [ 0  0  1]  [1 0 0]
        [(0, 1, 0), (0, 0, 1), (1, -1, -1), (1, 0, 0)], [ 1 -1 -1], [0 1 0]
        )
    """
    n = S[0].parent().codomain().dimension_relative() + 1
    cluster0 = [vector(list(pnt)) for pnt in S]
    C = ComplexField(prec = precision)
    R = RealField(prec = precision)
    for v in cluster0:
        v.change_ring(C)
    cluster = [v*1/R(v.norm()) for v in cluster0]
    Tr = matrix.identity(n)
    count = 0
    total_count = 0
    if S[1][1].parent().is_exact():
        prec = precision
    else:
        prec = S[1][1].prec()*2
    # ideally we would like to keep iterating indefinitely
    # however we stop once we lose enough precision
    while total_count < prec:
        total_count += 1
        cols = matrix(cluster).transpose()
        der1 = [R(col.hermitian_inner_product(col).real()) for col in cols]
        der2 = []
        for i in range(n):
            temp_list = []
            for j in range(i):
                inner_prod = R(cols[i].hermitian_inner_product(cols[j]).real())
                temp_list.append(2*(inner_prod))
            der2.append(temp_list)
        average_der1 = sum(der1)/n
        der1 = [r - average_der1 for r in der1]
        # absd is the absolute value of the derivative
        absd = 0
        der2_prime = [] # unpack the list der2
        for lst in der2:
            der2_prime += lst
        for d in der1 + der2_prime:
            absd += d**2
        if absd < eps**2:
            break
        for i in range(len(der2)):
            der2[i].append(der1[i])
        der3 = []
        for lst in der2:
            new_lst = lst + [C(0)]*(n - len(lst))
            der3.append(new_lst)
        # make symmetric matrix B1
        B1 = matrix(der3)
        for i in range(n-1):
            for j in range(i+1, n):
                B1[i, j] = B1[j, i]
        B = MatExp(-c*B1, eps=eps*(1**(-20))).change_ring(R)
        B *= B.determinant().nth_root(n).inverse_of_unit()
        Tr *= B
        cluster1 = [B*pnt for pnt in cluster]
        val = sum(pnt.hermitian_inner_product(pnt).log() for pnt in cluster1)
        if val > 0:
            c /= 2
        else:
            count += 1
            if count == 10:
                count = 0
                Tr *= R(Tr.det())**(-1/n)
                cluster1 = [pt*Tr for pt in cluster0]
            cluster = [v*1/R(v.norm()) for v in cluster1]
    Tr = Tr*Tr.transpose()
    Tr = Tr.change_ring(R)
    Tr = (Tr + Tr.transpose())/2
    U1 = Tr.LLL_gram()
    U1_inv = U1.inverse()
    return [v*U1_inv for v in cluster0], U1_inv, U1

def MatExp(mat, eps=1*10**-6):
    r"""
    An approximate matrix exponential, computed via power series.

    INPUT:

    - ``mat`` -- The matrix for which to compute `e^{\text{mat}}`

    - ``eps`` -- (default: `10^{-6}`) The precision after which
      to treat a number as zero

    OUPUT: A matrix which is an approximation of `e^{\text{mat}}`

    EXAMPLES::

        sage: from sage.schemes.projective.reduce_cluster import MatExp
        sage: m = matrix.identity(3)
        sage: MatExp(m)
        [9864101/3628800               0               0]
        [              0 9864101/3628800               0]
        [              0               0 9864101/3628800]
    """
    # computes the 'size' of a matrix
    def mat_prod(mat):
        lst = []
        for tup in mat:
            for ele in tup:
                lst.append(ele)
        vec = vector(lst)
        return vec.hermitian_inner_product(vec)

    n = len(mat[0])
    result = matrix.identity(n)
    k = 1
    m = result
    # approximate e^mat via power series
    while mat_prod(m) > eps**2:
        m *= 1/k * mat
        result += m
        k += 1
    return result