from sage.matrix.matrix_integer_2x2 import Matrix_integer_2x2 as mi2x2
import sage.matrix.all as matrix
#from libcpp.vector cimport vector

M2Z = matrix.MatrixSpace(ZZ,2)

@cached_function
def teich(a,p,M):
	R=pAdicField(p,M)
	return ZZ(R.teichmuller(a))

########################
### modsym utilities ###
########################

#### Hecke operator utilities ####

#cdef vector[mi2x2] unimod_matrices(int r, int s):
@cached_function
def unimod_matrices(r, s):
	"""connects the rational number r/s to 1/0"""
	if s != 0:
		list = convergents(r/s)
		v = [mi2x2(M2Z,[1, list[0].numerator(),0,list[0].denominator()], copy = False, coerce = False)]
		sig = -1
		for j in range(0,len(list)-1):
			a=list[j].numerator()
			c=list[j].denominator()
			b=list[j+1].numerator()
			d=list[j+1].denominator()
			v.append(mi2x2(M2Z, [sig * a, b, sig * c, d], copy = False, coerce = False))
			sig = -sig
		return v
	else:
		return []

@cached_function
def unimod_matrices2(r, s):
	"""connects the rational number 1/0 to r/s"""
	if s != 0:
		list = convergents(r/s)
		v = [mi2x2(M2Z,[-list[0].numerator(), 1, -list[0].denominator(), 0], copy = False, coerce = False)]
		sig = -1
		for j in range(0,len(list)-1):
			a=list[j].numerator()
			c=list[j].denominator()
			b=list[j+1].numerator()
			d=list[j+1].denominator()
			v.append(mi2x2(M2Z, [-b, sig * a, -d, sig * c], copy = False, coerce = False))
			sig = -sig
		return v
	else:
		return []

@cached_function
def prep_hecke_individual(ell,N,M,m):
	ans=[[] for a in range(len(M.mats))]
	for a in srange(ell):
		gama = mi2x2(M2Z, [1, a, 0, ell], copy = False, coerce = False)
		t = gama * M.mats[M.gens[m]]
		v = unimod_matrices2(t[0,0], t[1,0]) + unimod_matrices(t[0,1], t[1,1])
		for b in range(len(v)):
			A = v[b]
			i = M.P.index(A[1,0], A[1,1])
			j = M.P1_to_mats[i]
			B = M.mats[j]
			C = A.__invert__unit()
			gaminv = B * C
			ans[j] = ans[j] + [gaminv * gama]
	if N % ell != 0:
		gama=mi2x2(M2Z, [ell, 0, 0, 1], copy = False, coerce = False)
		t = gama * M.mats[M.gens[m]]
		v = unimod_matrices2(t[0,0], t[1,0]) + unimod_matrices(t[0,1], t[1,1])
		for b in range(len(v)):
			A = v[b]
			i = M.P.index(A[1,0], A[1,1])
			j = M.P1_to_mats[i]
			B = M.mats[j]
			C = A.__invert__unit()
			gaminv = B * C
			ans[j] = ans[j] + [gaminv * gama]

	return ans

@cached_function
def prep_hecke(ell,N,M):
	ans = []
	for m in range(len(M.gens)):
		ans = ans + [prep_hecke_individual(ell,N,M,m)]
	return ans

