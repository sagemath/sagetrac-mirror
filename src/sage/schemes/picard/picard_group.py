from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent

class BundledCurve(SageObject):
    def multiply_section_subspace(self, v, n, WD, m):
        pass # should have a generic implementation that uses a mul method.

    def multiply_subspaces(self, WD, n, WE, m, expected_dimension=None):
        E = None
        i = 0
        for v in WD:
            i += 1
            newrel = self.multiply_section_subspace(v, n, WE, m)
            E = newrel if E is None else E.stack(newrel)
            E.echelonize()
            r = E.rank()
            E = E.submatrix(0,0,r)
            if expected_dimension is not None and r >=expected_dimension:
                #print "mu_space_product: found expected %s-dim space after %s of %s multiplications of %s"%(expected_dimension,i,WD.nrows(),WE.nrows())
                break
        return E

    def divide_subspaces(self, V, m, restrictions, n, target, expected_codimension=None):
        KER = target.right_kernel_matrix()
        E = None
        i = 0
        for v in restrictions:
            i += 1
            if V is None:
                newrel = KER*self.mumat(v,n,m).transpose()
            else:
                newrel = KER*(V*self.mumat(v,n,m)).transpose()
            E = newrel if E is None else E.stack(newrel)
            E.echelonize()
            r = E.rank()
            E = E.submatrix(0,0,r)
            if expected_codimension is not None and r >= expected_codimension:
                #print "mu_preimage: found codim %s in %s tries of %s"%(expected_codimension,i,restrictions.nrows())
                break
        else:
            pass
            #print "mu_preimage: found codim %s instead of %s in %s tries"%(r,expected_codimension,i)
        return E.right_kernel_matrix()

class BundledCurve_plane(SageObject):
    def __init__(self, C):
        self.curve = C
        self.R = R = C.defining_polynomial().parent()
        self.k = R.base_ring()

    def _hom_deg_part(self, n):
        I = [tuple(v) for v in IntegerVectors(n, self.R.ngens())]
        B = [self.R({v:1}) for v in I]
        V = VectorSpace(self.k, len(B))
        def tovec(p):
            return V([p.coefficient(m) for m in B])
        def topol(v):
            return P(dict(zip(I,v)))
        return tovec, topol

    @cached_method
    def _hompart(self, n):
        tovec, topol = self._hom_deg_part(n)
        F = self.curve.defining_polynomial()
        dim = int((n+1)*(n+2)/2)
        if n < F.degree():
            reps = [topol(w) for w in matrix(self.k, dim, dim, 1)]
            l = dim
        else:
            In = matrix([tovec(m*F) for m in monomials_of_degree(self.R, n - F.degree())])
            toWn = In.right_kernel_matrix()
            l = toWn.nrows()
            idmat = matrix(self.k, dim, dim, 1)
            reps = [topol(idmat[i]) for i in toWn.pivots()]
            tovec = lambda g: toWn*tovec(g)
            topol = lambda v: sum(a*b for a,b in zip(reps, v))
        return l, tovec, topol, reps

    @cached_method
    def _W0(self,n):
        l,_,_,_=self._hompart(n)
        return matrix(self.k, l, l, 1)

    @cached_method
    def _mu(self,n,m):
        ln,_,_,repsn = self._hompart(n)
        lm,_,_,repsm = self._hompart(m)
        lt,tovec,_,_ = self._hompart(n+m)
        return (matrix([sum((tovec(a*b) for b in repsm),[]) for a in repsn]), lm, lt)

    def _mumat(self, v, n, m):
        M,lm,lsum=self._mu(n,m)
        return (matrix(v)*M)._matrices_from_rows(lm,lsum)[0]

    def multiply_section_subspace(self, v, n, WD, m):
        return WD * self._mumat(v, n, m)

    def point(self,pt,n):
        assert self.curve.defining_polynomial()(*pt) == 0
        l,tovecr,_,reps=self.hompart(n)
        WD=matrix([r(pt) for r in reps]).right_kernel_matrix()
        return divisor(self,n,WD,1,-n)

    def zeros(self,n):
        if n in self.zeros_cache:
            return self.zeros_cache[n]
        l,tovec,topol,reps=self.hompart(n)
        R=self.R
        b = self.curve.degree()
        g = self.curve.arithmetic_genus()
        k0=((n*b-2*g)/b).floor()
        z=R.gen(2)
        L=[]
        for k in range(k0+1):
            E=matrix([tovec(z^k*m) for m in monomials_of_degree(R,n-k)])
            E.echelonize()
            E=E.submatrix(0,0,E.rank())
            L.append(divisor(self,n,E,k*b,k-n))
        self.zeros_cache[n]=L
        return L

    def initpickvector(self,pt,n):
        l, tovec, topol, reps = self._hompart(n)
        G1=self.point(pt,n)
        L=[self.zeros(n)[0]]
        b = self.curve.degree()
        g = self.curve.arithmetic_genus()
        N=((n*b/2).ceil())
        for i in range(N):
            L.append(L[-1].add1(G1))
        spc=[l.WD for l in L]
        for i in range(N,n*b):
            spc.append(L[i//2].raw_add1(L[i -(i//2)]))
        spc.reverse()
        r=0
        B=None
        for s in spc:
            if s.rank() > r:
                assert s.rank() == r+1
                for v in s:
                    C = matrix(v) if B is None else B.stack(matrix(v))
                    if C.rank() > r:
                        B=C
                        r+=1
                        break
        C=C.rows()
        C.reverse()
        C=matrix(C)
        self.pickvector_cache[n]=(C,C^(-1))

    def pickvector(self,WD,n):
        if n in self.pickvector_cache:
            C,Cinv = self.pickvector_cache[n]
            B=WD*Cinv
            B.echelonize()
            return B[B.rank()-1]*C
        else:
            return WD[-1]

class BundledCurve_modular(SageObject):
    def __init__(self, Gamma):
        self.gamma = Gamma

class PicardGroup(Parent):
    def __init__(self, C):
        if isinstance(C, ProjectivePlaneCurve):
            self.curve = BundledCurve_plane(C)

