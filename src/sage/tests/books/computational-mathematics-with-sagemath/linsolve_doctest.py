## -*- encoding: utf-8 -*-
r"""
This file (./linsolve_doctest.sage) was *autogenerated* from ./linsolve.tex,
with sagetex.sty version 2011/05/27 v2.3.1.
It contains the contents of all the sageexample environments from this file.
You should be able to doctest this file with:
sage -t ./linsolve_doctest.sage
It is always safe to delete this file; it is not used in typesetting your
document.

Sage example in ./linsolve.tex, line 294::

  sage: def cond_hilbert(n):
  ....:     A = matrix(QQ, [[1/(i+j-1) for j in [1..n]] for i in [1..n]])
  ....:     return A.norm(Infinity) * (A^-1).norm(Infinity)

Sage example in ./linsolve.tex, line 346::

  sage: def diff_hilbert(n):
  ....:     x = vector(QQ,[1 for i in range(0,n)])
  ....:     A = matrix(QQ, [[1/(i+j-1) for j in [1..n]] for i in [1..n]])
  ....:     y = A*x
  ....:     A[n-1,n-1] = (1/(2*n-1))*(1+1/(10^5))   # modifies the matrix
  ....:     s = A\y
  ....:     return max(abs(float(s[i]-x[i])) for i in range(0,n))

Sage example in ./linsolve.tex, line 396::

  sage: def hilbert_diff(n):
  ....:     j = var("j")
  ....:     f = lambda i: sum(1/(i+j-1),j,1,n)
  ....:     y = vector(RDF, [f(i+1) for i in range(0,n)])
  ....:     A = matrix(RDF, [[1/(i+j-1) for i in [1..n]] for j in [1..n]])
  ....:     x = A.solve_right(y)
  ....:     return max(abs(x[i]-1.0) for i in range(0,n))

Sage example in ./linsolve.tex, line 543::

  sage: n = 20; cost = (n+1)*factorial(n); cost
  51090942171709440000

Sage example in ./linsolve.tex, line 559::

  sage: v = 3*10^9
  sage: print("%3.3f" % float(cost/v/3600/24/365))
  540.028

Sage example in ./linsolve.tex, line 653::

  sage: A = matrix(RDF, [[-1,2],[3,4]])
  sage: b = vector(RDF, [2,3])
  sage: x = A\b; x
  (-0.20000000000000018, 0.9000000000000001)

Sage example in ./linsolve.tex, line 666::

  sage: x = A.solve_right(b)

Sage example in ./linsolve.tex, line 678::

  sage: A = matrix(RDF, [[-1,2],[3,4]])
  sage: P, L, U = A.LU()

Sage example in ./linsolve.tex, line 723::

  sage: eps = 1e-16
  sage: y = (1-2*eps)/(1-eps)
  sage: x = (1-y)/eps
  sage: x, y
  (1.11022302462516, 1.00000000000000)

Sage example in ./linsolve.tex, line 736::

  sage: 1. + eps == 1.
  True

Sage example in ./linsolve.tex, line 758::

  sage: y = (1-2*eps)/(1-eps)
  sage: x = 2-y
  sage: x, y
  (1.00000000000000, 1.00000000000000)

Sage example in ./linsolve.tex, line 803::

  sage: A = random_matrix(RDF, 1000)
  sage: b = vector(RDF, range(1000))
  sage: c = vector(RDF, 2*list(range(500)))

Sage example in ./linsolve.tex, line 858::

  sage: m = random_matrix(RDF, 10)
  sage: A = transpose(m)*m
  sage: C = A.cholesky()

Sage example in ./linsolve.tex, line 946::

  sage: A = random_matrix(RDF,6,5)
  sage: Q, R = A.QR()

Sage example in ./linsolve.tex, line 1012::

  sage: A = matrix(RDF, [[1,3,2],[1,2,3],[0,5,2],[1,1,1]])
  sage: U, Sig, V = A.SVD()
  sage: A1 = A - U*Sig*transpose(V); A1 # abs tol 1e-9
  [ 2.220446049250313e-16                    0.0                    0.0]
  [3.3306690738754696e-16 -4.440892098500626e-16 -4.440892098500626e-16]
  [-9.298117831235686e-16 1.7763568394002505e-15 -4.440892098500626e-16]
  [ 4.440892098500626e-16 -8.881784197001252e-16 -4.440892098500626e-16]

Sage example in ./linsolve.tex, line 1123::

  sage: A = matrix(RDF, [[1,3,2],[1,4,2],[0,5,2],[1,3,2]])
  sage: b = vector(RDF, [1,2,3,4])
  sage: Z = transpose(A)*A
  sage: C = Z.cholesky()
  sage: R = transpose(A)*b
  sage: Z.solve_right(R) # abs tol 2e-13
  (-1.5000000000000135, -0.5000000000000085, 2.7500000000000213)

Sage example in ./linsolve.tex, line 1188::

  sage: A = matrix(RDF, [[1,3,2],[1,4,2],[0,5,2],[1,3,2]])
  sage: b = vector(RDF, [1,2,3,4])
  sage: Q, R = A.QR()
  sage: R1 = R[0:3,0:3]
  sage: b1 = transpose(Q)*b
  sage: c = b1[0:3]
  sage: R1.solve_right(c)  # abs tol 1e-13
  (-1.499999999999999, -0.49999999999999867, 2.749999999999997)

Sage example in ./linsolve.tex, line 1202::

  sage: Z = A.transpose()*A
  sage: Z.norm(Infinity)*(Z^-1).norm(Infinity)  # abs tol 2e-10
  1992.3750000000168

Sage example in ./linsolve.tex, line 1256::

  sage: A = matrix(RDF, [[1,3,2],[1,3,2],[0,5,2],[1,3,2]])
  sage: B = vector(RDF, [1,2,3,4])
  sage: U, Sig, V = A.SVD()
  sage: m = A.ncols()
  sage: x = vector(RDF, [0]*m)
  sage: lamb = vector(RDF, [0]*m)
  sage: for i in range(0,m):
  ....:     s = Sig[i,i]
  ....:     if s!=0.0:
  ....:         lamb[i]=U.column(i)*B/s
  sage: x = V*lamb; x  # random
  (0.2370370370370367, 0.4518518518518521, 0.3703703703703702)

Sage example in ./linsolve.tex, line 1289::

  sage: m = 3; [ Sig[i,i] for i in range(0,m) ]  # abs tol 1e-15
  [8.309316833256451, 1.3983038884881154, 0.0]

Sage example in ./linsolve.tex, line 1358::

  sage: A = matrix(RDF, [[1,2],[3,4],[5,6],[7,8]])

Sage example in ./linsolve.tex, line 1366::

  sage: set_random_seed(0) # to get reproducible values below

Sage example in ./linsolve.tex, line 1369::

  sage: th = 0.7
  sage: R = matrix(RDF, [[cos(th),sin(th)],[-sin(th),cos(th)]])
  sage: B = (A + 0.1*random_matrix(RDF,4,2)) * transpose(R)

Sage example in ./linsolve.tex, line 1374::

  sage: C = transpose(B)*A
  sage: U, Sigma, V = C.SVD()
  sage: Q = U*transpose(V)

Sage example in ./linsolve.tex, line 1381::

  sage: Q  # rel tol 1e-15
  [ 0.7612151656410957   0.648499399843978]
  [-0.6484993998439779  0.7612151656410955]
  sage: R
  [0.7648421872844885  0.644217687237691]
  [-0.644217687237691 0.7648421872844885]

Sage example in ./linsolve.tex, line 1562::

  sage: set_random_seed(0) # to get reproducible values below

Sage example in ./linsolve.tex, line 1675::

  sage: A = matrix(RDF, [[1,3,2],[1,2,3],[0,5,2]])

Sage example in ./linsolve.tex, line 1698::

  sage: A = matrix(RDF,[[1,3,2],[1,2,3],[0,5,2]])
  sage: mu = 0.56
  sage: AT = A - mu*identity_matrix(RDF,3)
  sage: X = vector(RDF,[1 for i in range(0,A.nrows())])
  sage: lam_old = 0
  sage: for i in range(1,1000): # abs tol 1e-9
  ....:     Z = AT.solve_right(X)
  ....:     X = Z/Z.norm()
  ....:     lam = X.dot_product(A*X)
  ....:     s = abs(lam - lam_old)
  ....:     print("{i} s={s} lambda={lam}".format(i=i, s=s, lam=lam))
  ....:     lam_old = lam
  ....:     if s<1.e-10:
  ....:         break
  1 s=0.56423627407 lambda=0.56423627407
  2 s=0.00371649959176 lambda=0.560519774478
  3 s=2.9833340176e-07 lambda=0.560519476145
  4 s=3.30288019157e-11 lambda=0.560519476112
  sage: X  # abs tol 1e-15
  (0.9276845629439007, 0.10329475725387141, -0.3587917847435305)

Sage example in ./linsolve.tex, line 1724::

  sage: A*X-lam*X # abs tol 1e-14
  (2.886579864025407e-15, 1.672273430841642e-15, 8.326672684688674e-15)

Sage example in ./linsolve.tex, line 1801::

  sage: m = matrix(RDF, [[1,2,3,4],[1,0,2,6],[1,8,4,-2],[1,5,-10,-20]])
  sage: Aref = transpose(m)*m
  sage: A = copy(Aref)
  sage: for i in range(0,20):
  ....:     Q, R = A.QR()
  ....:     A = R*Q

Sage example in ./linsolve.tex, line 1835::

  sage: Aref.eigenvalues()  # abs tol 1e-12
  [585.0305586200212, 92.91426499150643, 0.03226690899408103, 4.022909479477674]

Sage example in ./linsolve.tex, line 1895::

  sage: A = matrix(RDF, [[1,3,2],[1,2,3],[0,5,2]])
  sage: eigen_vals, eigen_vects = A.eigenmatrix_right()
  sage: eigen_vals  # abs tol 1e-12
  [   6.39294791648918                 0.0                 0.0]
  [                0.0   0.560519476111939                 0.0]
  [                0.0                 0.0 -1.9534673926011215]
  sage: eigen_vects # abs tol 2e-15
  [ 0.5424840601106511  0.9276845629439008 0.09834254667424457]
  [ 0.5544692861094349 0.10329475725386986  -0.617227053099068]
  [ 0.6310902116870117 -0.3587917847435306   0.780614827194734]

Sage example in ./linsolve.tex, line 1939::

  sage: def pol2companion(p):
  ....:     n = len(p)
  ....:     m = matrix(RDF,n)
  ....:     for i in range(1,n):
  ....:         m[i,i-1]=1
  ....:     m.set_column(n-1,-p)
  ....:     return m

Sage example in ./linsolve.tex, line 1965::

  sage: q = vector(RDF,[1,-1,2,3,5,-1,10,11])
  sage: comp = pol2companion(q); comp
  [  0.0   0.0   0.0   0.0   0.0   0.0   0.0  -1.0]
  [  1.0   0.0   0.0   0.0   0.0   0.0   0.0   1.0]
  [  0.0   1.0   0.0   0.0   0.0   0.0   0.0  -2.0]
  [  0.0   0.0   1.0   0.0   0.0   0.0   0.0  -3.0]
  [  0.0   0.0   0.0   1.0   0.0   0.0   0.0  -5.0]
  [  0.0   0.0   0.0   0.0   1.0   0.0   0.0   1.0]
  [  0.0   0.0   0.0   0.0   0.0   1.0   0.0 -10.0]
  [  0.0   0.0   0.0   0.0   0.0   0.0   1.0 -11.0]
  sage: roots = comp.eigenvalues(); roots  # abs tol 1e-12
  [0.3475215101190289 + 0.5665505533984981*I, 0.3475215101190289 - 0.5665505533984981*I,
   0.34502377696179265 + 0.43990870238588275*I, 0.34502377696179265 - 0.43990870238588275*I,
   -0.5172576143252197 + 0.5129582067889322*I, -0.5172576143252197 - 0.5129582067889322*I,
   -1.3669971645459291, -9.983578180965276]

Sage example in ./linsolve.tex, line 2121::

  sage: reset()

Sage example in ./linsolve.tex, line 2124::

  sage: def eval(P,x):
  ....:     if len(P) == 0:
  ....:         return 0
  ....:     else:
  ....:         return P[0]+x*eval(P[1:],x)

Sage example in ./linsolve.tex, line 2133::

  sage: def pscal(P,Q,lx):
  ....:     return float(sum(eval(P,s)*eval(Q,s) for s in lx))

Sage example in ./linsolve.tex, line 2139::

  sage: def padd(P,a,Q):
  ....:     for i in range(0,len(Q)):
  ....:         P[i] += a*Q[i]

Sage example in ./linsolve.tex, line 2149::

  sage: class BadParamsforOrthop(Exception):
  ....:     def __init__(self, degreeplusone, npoints):
  ....:         self.deg = degreeplusone - 1
  ....:         self.np = npoints
  ....:     def __str__(self):
  ....:         return "degree: " + str(self.deg) + \
  ....:                " nb. points: " + repr(self.np)

Sage example in ./linsolve.tex, line 2160::

  sage: def orthopoly(n,x):
  ....:     if n > len(x):
  ....:         raise BadParamsforOrthop(n, len(x))
  ....:     orth = [[1./sqrt(float(len(x)))]]
  ....:     for p in range(1,n):
  ....:         nextp = copy(orth[p-1])
  ....:         nextp.insert(0,0)
  ....:         s = []
  ....:         for i in range(p-1,max(p-3,-1),-1):
  ....:             s.append(pscal(nextp, orth[i], x))
  ....:         j = 0
  ....:         for i in range(p-1,max(p-3,-1),-1):
  ....:             padd(nextp, -s[j], orth[i])
  ....:             j += 1
  ....:         norm = sqrt(pscal(nextp, nextp, x))
  ....:         nextpn = [nextp[i]/norm for i in range(len(nextp))]
  ....:         orth.append(nextpn)
  ....:     return orth

Sage example in ./linsolve.tex, line 2209::

  sage: set_random_seed(3) # to get reproducible values below

Sage example in ./linsolve.tex, line 2212::

  sage: L = 40
  sage: X = [100*float(i)/L for i in range(40)]
  sage: Y = [float(1/(1+25*X[i]^2)+0.25*random()) for i in range(40)]
  sage: n = 15; orth = orthopoly(n, X)

Sage example in ./linsolve.tex, line 2222::

  sage: coeff = [sum(Y[j]*eval(orth[i],X[j]) for j in
  ....:         range(0,len(X))) for i in range(0,n)]

Sage example in ./linsolve.tex, line 2230::

  sage: polmin = [0 for i in range(0,n)]
  sage: for i in range(0,n):
  ....:     padd(polmin, coeff[i], orth[i])
  sage: p = lambda x: eval(polmin, x)
  sage: plot(p(x), x, 0, X[len(X)-1])
  Graphics object consisting of 1 graphics primitive

Sage example in ./linsolve.tex, line 2609::

  sage: from scipy.sparse.linalg.dsolve import *
  sage: from scipy.sparse import lil_matrix
  sage: from numpy import array
  sage: n = 200
  sage: n2 = n*n
  sage: A = lil_matrix((n2, n2))
  sage: h2 = 1./float((n+1)^2)
  sage: for i in range(0,n2):
  ....:    A[i,i]=4*h2+1.
  ....:    if i+1<n2: A[i,int(i+1)]=-h2
  ....:    if i>0:    A[i,int(i-1)]=-h2
  ....:    if i+n<n2: A[i,int(i+n)]=-h2
  ....:    if i-n>=0: A[i,int(i-n)]=-h2
  sage: Acsc = A.tocsc()
  sage: b = array([1 for i in range(0,n2)])
  sage: solve = factorized(Acsc) # LU factorization
  sage: S = solve(b)             # resolution

Sage example in ./linsolve.tex, line 2784::

  sage: from numpy.random import seed
  sage: seed(0) # to get reproducible values below

Sage example in ./linsolve.tex, line 2828::

  sage: from scipy import sparse
  sage: from numpy.linalg import *
  sage: from numpy import array
  sage: from numpy.random import rand
  sage: def power(A,x,N):             # power iteration
  ....:     for i in range(N):
  ....:         y = A*x
  ....:         z = y/norm(y)
  ....:         lam = sum(x*y)
  ....:         s = norm(x-z)
  ....:         print("{i} s={s} lambda={lam}".format(i=i, s=s, lam=lam))
  ....:         if s < 1e-3:
  ....:             break
  ....:         x = z
  ....:     return x
  sage: n = 1000
  sage: m = 5
  sage: # build a stochastic matrix of size n
  sage: # with m non-zero coefficients per row
  sage: A1 = sparse.lil_matrix((n, n))
  sage: for i in range(0,n):
  ....:     for j in range(0,m):
  ....:         l = int(n*rand())
  ....:         A1[l,i] = rand()
  sage: for i in range(0,n):
  ....:     s = sum(A1[i,0:n])
  ....:     A1[i,0:n] /= s
  sage: At = A1.transpose().tocsc()
  sage: x = array([rand() for i in range(0,n)])
  sage: # compute the dominant eigenvalue
  sage: # and the associated eigenvector
  sage: y = power(At, x, 5)  # rel tol 1e-10
  0 s=17.0241218112 lambda=235.567796432
  1 s=0.39337173784 lambda=0.908668201953
  2 s=0.230865716856 lambda=0.967356896036
  3 s=0.134156683993 lambda=0.986660315554
  4 s=0.0789423487458 lambda=0.995424635219
"""
