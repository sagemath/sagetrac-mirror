def Skewify(A):
    """Input: A, where A is either a partition or a skew partition.
    Output: a skew partition associated to A if A was a partition and A itself if was a
    skew partition."""
    if isinstance(A[0],list)==False:
        return [A,[]]
    else:
        return A

def MaxContent(t = SkewTableau([[None,1,2],[3,4],[5],[6]])):
    """Input: t is a skew tableau of shape lambda/mu, where lambda is a partition (aka outer
    shape) and mu is a partition whose Young diagram contained in the Young diagram of lambda 
    (aka inner shape).
    Output: maximum value of j-i over all cells (i,j) in skew partition t.
    Strategy: Visits every row in the Young diagrams of lambda and mu and returns the maximum
    content on the first row where the parts differ in size."""
    A = t.outer_shape()
    B = t.inner_shape()
    if A == B:
        return None
    else:
        for i in range(len(B)):
            if A[i] > B[i]:
                return A[i] - (i + 1)
        return A[len(B)] - (len(B) + 1)

def MinContent(t = SkewTableau([[None,1,2],[3,4],[5],[6]])):
    """Input: t is a skew tableau with shape lambda/mu, where lambda is a partition 
    (aka outer shape) and mu is a partition whose Young diagram contained in the Young diagram
    of lambda (aka inner shape)
    Output: minimum value of j-i over all cells (i,j) in skew partition t. 
    Strategy: Computes the conjugates of the outer and inner shapes and mimics the MaxContent(t)
    function from there."""
    A = t.outer_shape().conjugate()
    B = t.inner_shape().conjugate()
    if A == B:
        return None
    else:
        for i in range(len(B)):
            if A[i] > B[i]:
                return - A[i] + (i + 1)
        return - A[len(B)] + (len(B) + 1)
# EXAMPLES:
# print MaxContent()
# print MinContent()
#
# t1=SkewTableau([[None,None,None,None,None,None],[None,None,None,None,None,None],[None,None,None,3,3,3],
# [None,None,None,4,4],[None,None,None],[None,None,6],[None,None],[None,8],[None,9],[10],[11],[12]])
# print MaxContent(t1)
# print MinContent(t1)
#
# t2=SkewTableau([[None,None,None,None,None,None],[None,None,None,None,None,None],[None,None,None,None,None,None],
# [None,None,None,None,None],[None,None,None],[None,None,None],[None,None],[None,None],[None,None],[None],[None],[None]])
# print MaxContent(t2)
# print MinContent(t2)
#
# t3=SkewTableau([[None,None,None,None,None,1],[None,None,None,None,None,2],[None,None,None,None,None,3],
# [None,None,None,None,None],[None,None,None],[None,None,6],[None,None],[None],[None]])
# print MaxContent(t3)
# print MinContent(t3)

def TableOfContents(t = SkewTableau([[None,1,2,3],[4,5,6,7],[8,9]])):
    """Input: t, a skew tableau with shape lambda/mu, where lambda is a partition 
    (aka outer shape) and mu is a partition whose Young diagram contained in the 
    Young diagram of lambda (aka inner shape).
    Output: a list of pairs (c,e) where c is an integer between MinContent(t) and 
    MaxContent(t) inclusive, and e is an entry of t whose content(cell containing e) = c."""
    return [[i,t.entries_by_content(i)] for i in range(MinContent(t),MaxContent(t)+1)]
# EXAMPLE:
# TableOfContents()

def DictOfToC(*T):
    """Input: T, an arbitrary tuple of skew tableaux t_1, t_2, ..., t_k of respective shapes
    lambda_k/mu_k, where lambda_k is a partition (aka outer shape) and mu_k is a partition whose
    Young diagram contained in the Young diagram of lambda_k (aka inner shape). Input could also
    be T, an arbitrary tuple of Young tableaux t_1, t_2, ..., t_k of respective shapes lambda_k 
    (Here, mu_k is taken to be empty for all k).
    Output: a dictionary of triples (c,i): e where i is the index of tableau t_i, c is an integer
    between MinContent(t_i) and MaxContent(t_i) inclusive, and e is a list of entries of t_i 
    whose content is c."""
    D={}
    for i in range(len(T)):
        t=SkewTableau(T[i])
        D.update( dict( ((c,i), t.entries_by_content(c) ) for c in range(MinContent(t),MaxContent(t)+1) ) )
    A=[b[0] for b in D.keys()]
    Max=max(A)
    Min=min(A)
    return D,Max,Min
# EXAMPLE:
# T=( SkewTableau([[None,None,1,2],[None,3,4,5],[6,7,8],[9]]), SkewTableau([[None,1,2,3],[4,5,6,7],[8,9]]), 
# SkewTableau([[None,None,1,2,3],[None,4,5],[None,6],[7,8],[9]]) )
# from pprint import pprint
# D,Max,Min=DictOfToC(*T)
# pprint(D)
# print Max,Min

def CountInv(A,B):
    if len(A)==0 or len(B)==0:
        return 0
    L=[(a,b) for a in A for b in B if a>b]
    return len(L)

# EXAMPLES:
# CountInv([],[])
# CountInv([],[1,2])
# CountInv([2,3],[])
# CountInv([1,3],[2,3,4])
# CountInv([2,3,4],[1,2])
 
def NumberOfInversions(*T):
    """Input: T, an arbitrary list of skew tableaux t_1, t_2, ..., t_k of respective shapes
    lambda_k/mu_k, where lambda_k is a partition (aka outer shape) and mu_k is a partition whose 
    Young diagram contained in the Young diagram of lambda_k (aka inner shape)
    Output: inv(*T), the number of inversions occuring in *T, that is the number of all pairs 
    ( (c(u),i,e(u)), (c(v),j,e(v)) ) where u and v are cells in t_i and t_j, have contents 
    c(u) and c(j), and have entries e(u) and e(v) respectively, where either:
    (i) c(u)=c(v) and i<j and e(u)>e(v)
    (ii) c(u)+1=c(v) and i>j and e(u)>e(v)"""
    D,Max,Min=DictOfToC(*T)
    n=len(T)
    inv=0
    for key in sorted(D.keys()):
        for j in range(n):
            k1=(key[0],j)
            k2=(key[0]+1,j)
            if key[1]<j and k1 in D.keys():
                inv+=CountInv(D[key],D[k1])
            elif key[1]>j and k2 in D.keys():
                inv+=CountInv(D[key],D[k2])
    return inv
# EXAMPLE:
# M1=( SkewTableau([[None,None,1,2],[None,3,4,5],[6,7,8],[9]]), SkewTableau([[None,1,2,3],[4,5,6,7],[8,9]]), 
# SkewTableau([[None,None,1,2,3],[None,4,5],[None,6],[7,8],[9]]) )
#
# from pprint import pprint
# D,Max,Min=DictOfToC(*T)
# pprint(D)
# print Max,Min
#
# NumberOfInversions(*M1)
# NumberOfInversions( SkewTableau([[None,None,1,2],[None,3,4,5],[6,7,8],[9]]), SkewTableau([[None,1,2,3],[4,5,6,7],[8,9]]), SkewTableau([[None,None,1,2,3],[None,4,5],[None,6],[7,8],[9]]) )
#
# M2=( Tableau([ [2,3],[4] ]), Tableau([ [3,3],[4,4] ]), Tableau([ [1,2,3],[4],[5] ]) ) 
# D2,Max2,Min2=DictOfToC(*T)
# pprint(D2)
# print Max2,Min2
# NumberOfInversions(*M2) NumberOfInversions(*c)

class LLTPolynomial:
    """Computes LLT Polynomial given the number of variables and a sequence of skew partitions
    (with inner and outer shapes specified)."""
    def __init__(self,max_alphabet,*nu):
        self.max_alphabet=max_alphabet
        self.nu=tuple(Skewify(i) for i in nu)
        for i in range(len(self.nu)):
            min_index = 0
            bool_min = False
            if len(self.nu[i][0])>min_index and len(self.nu[i][0])>self.max_alphabet:
                bool_min = True
                min_index = len(self.nu[i][0])
        if bool_min == True:
            print "max_alphabet should be at least",min_index,"!"

    def __repr__(self):
        return "<LLTPolynomial with max alphabet=%r and nu=%r>" % (self.max_alphabet, self.nu)
    def Nu(self):
        return self.nu
    def MaxAlphabet(self):
        return self.max_alphabet
    def LLT(self):
        """Input: max_alphabet, the number of variables n in x=(x_1,...,x_n) and S, an arbitrary tuple of shape of skew diagrams.
        Output: LLT polynomial p given by the sum of q^inv(T)x^T over all T in SSYT(S), S=(nu_1,nu_2,...,nu_k) and x=(x_1,...,x_n)."""
        Sym = SymmetricFunctions(QQ['q'])
        s = Sym.schur()
        P = PolynomialRing(QQ['q'], ['x%s'%i for i in range(1,self.max_alphabet+1)])
        x = P.gens()
        q = (P.base_ring()).gens()
        C = cartesian_product([SemistandardSkewTableaux(self.nu[i],max_entry=self.max_alphabet) for i in range(len(self.nu))])
        p = sum(q[0]**NumberOfInversions(*c)*prod(prod(x[i]**(c[j].weight()[i]) for i in range(len(c[j].weight()))) for j in range(len(c))) 
            for c in C)
        f = s.from_polynomial(p).restrict_partition_lengths(self.max_alphabet,False)
        return f
# EXAMPLES:
# S = ([3,3],[2,2,1])
# A = LLTPolynomial(3,*S)
# p = A.LLT()
# print"The LLT polynomial with max alphabet = ",A.MaxAlphabet(),"\n and nu = ",A.Nu()," is\n",p,"\n"
# The LLT polynomial with max alphabet =  3 
# and nu =  ([[2, 1], []], [[2, 1], []], [[1], []])  is 
# q^6*s[5, 4, 2] + q^5*s[5, 5, 1] 
#
# B = LLTPolynomial(3,[2,1],[2,1],[1])
# f = A1.LLT()
# print"The LLT polynomial with max alphabet = ",B.MaxAlphabet(),"\n and nu = ",B.Nu()," is \n",f,"\n"
# The LLT polynomial with max alphabet =  3 
# and nu =  ([[2, 1], []], [[2, 1], []], [[1], []])  is 
# (q^3+q^2)*s[4, 3] + q*s[5, 2] 

