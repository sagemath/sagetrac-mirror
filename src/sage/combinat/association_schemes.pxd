# -*- coding: utf-8 -*-

class AssociationScheme:

    def _is_association_scheme(self):

        from itertools import product

        # check matrix size
        if self._matrix.ncols() != self._nX or \
           self._matrix.nrows() != self._matrix.ncols():
            raise ValueError("matrix has wrong size")

        # check R0
        for i in range(self._nX):
            if self._matrix[i][i] != 0:
                raise ValueError("identity is not R_0")

        # check symmetries
        r1 = self._r + 1
        for i in range(r1):
            if self._R[i] is None:
                self.R(i)

            RiT = set([(y,x) for (x,y) in self._R[i]])
            for j in range(self._r):
                if self._R[j] is None:
                    self.R(j)

                if RiT == self._R[j]:
                    break
            else:
                raise ValueError(f"There is no relation R_{i}^T")

        for i, j, k in product(range(r1), range(r1), range(r1)):
            if self.p(i, j, k) != self.p(j, i, k):
                raise ValueError("The intersection numbers are not symmetric")

        # chck intersection numbers
        # IMPORTANT: we rely on self._R[i] != None for all i
        #            this is safe due to the symmetry check
        for i, j, k in product(range(r1), range(r1), range(r1)):
            for (x,y) in self._R[k]:
                count = 0
                for z in range(self._nX):
                    if self._matrix[x, z] == i and self._matrix[z, y] == j:
                        count += 1

                if count != self.p(i,j,k):
                    raise ValueError("There are no valid intersection numbers")

    def __init__(self, points, relations=None, matrix=None, check=True):
        if relations is None and matrix is None:
            raise ValueError("Please specify one of relations or matrix")

        def relationToInt(R):
            return set([(self._XtoInt[x], self._XtoInt[y]) for (x,y) in R])

        self._X = list(points)
        self._nX = len(self._X)
        self._XtoInt = {x: i for i, x in enumerate(self._X)}

        if matrix is not None:
            self._matrix = Matrix(matrix)

            # compute number of classes
            self._r = 0
            for r in self._matrix:
                for c in r:
                    if c > self._r: self._r = c

            if relations is not None:
                self._R = list(map(relationToInt, relations))

                # now check consistency
                if len(self._R) != self._r:
                    raise ValueError(("matrix and relations are"
                                      "inconsistent with each others"))

                for k, s in enumerate(self._R):
                    for (x,y) in s:
                        if self._matrix[x, y] != k:
                            raise ValueError(("matrix and relations are"
                                              "inconsistent with each others"))
            else:
                self._R = [None]*self._r

        else:
            self._R = list(map(relationToInt, relations))
            self._r = len(self._R)

            # construct matrix
            self._matrix = Matrix([[None for i in range(self._nX)]
                                   for j in range(self._nX)])

            for k, s in enumerate(self._R):
                for (x,y) in s:
                    self._matrix[x, y] = k

        # empty intersection array
        r1 = self._r+1
        self._P = [[[None for k in range(r1)]
                    for j in range(r1)] for i in range(r1)]

        if check:
            self._is_association_scheme()

    def points(self):
        return self._X

    def num_points(self):
        return self._nX

    def matrix(self):
        return self._matrix

    def num_classes(self):
        return self._r

    def has_relation(self,x,y,i):
        return self.which_relation(x,y)  == i

    def which_relation(self,x,y):
        return self._matrix[self._XtoInt[x]][self._XtoInt[y]]

    def R(self,k):
        if self._R[k] is not None:
            return [(self._X[i], self._X[j]) for (i,j) in self._R[k]]

        Rk = []
        self._R[k] = set()
        for i in range(self._nX):
            for j in range(i,self._nX):
                if self._matrix[i][j] == k:
                    self._R[k].add((i,j))
                    Rk.append((self._X[i],self._X[j]))
                    Rk.append((self._X[j],self._X[i]))
        return Rk

    def p(self,i,j,k):
        if self._P[i][j][k] is not None: return self._P[i][j][k]

        nX = self._nX
        self._P[i][j][k] = 0

        # Find (x,y) in R_k
        found = False
        for x in range(nX):
            for y in range(nX):
                if self._matrix[x][y] == k:
                    found = True
                    break
            if found:
                break

        if self._matrix[x][y] != k:
            raise RuntimeError("something bad happend in code")

        # count z s.t. (x,z) in R_i and (z,y) in R_j
        for z in range(nX):
            if self._matrix[x][z] == i and self._matrix[z][y] == j:
                self._P[i][j][k] += 1

        return self._P[i][j][k]

    def is_symmetric(self):
        for i in range(self._r + 1):
            for (x,y) in self.R(i):
                if not self.has_relation(y,x,i):
                    return False

        return True

    def is_pseudocyclic(self):
        #we need p_{ii}^0 to be a constant f for 1<= i <= r
        #we need sum_{i} p_{ii}^k = f-1 for 1 <= k <= r
        r1 = self._r+1
        f = self.p(1,1,0)
        for i in range(2,r1):
            if self.p(i,i,0) != f:
                return False

        for k in range(1,r1):
            sumP = 0
            for i in range(r1):
                sumP += self.p(i,i,k)

            if sumP != f-1:
                return False
        return True
