from copy import deepcopy
from sage.all import sqrt, random, var, vector, RR, log
from numpy import array, mat, linalg, dot, concatenate

def projection(character, degree, genus, perm, v, n, lab):
    r"""
    Return the value of the projection p_i of v on the coordinate (lab, n)
    """
    sum = 0
    for indice in range(genus):
        sum += character[indice]*v.val(perm[indice].inverse()(n),lab)
    return (degree/genus*sum)


class VectPaths(object):
    r"""
    A vectors whose coordinates are referenced by a familly of paths on the corresponding translation surface
    which generate the relative homology group.
    """
    def __init__(self,labels,degrees):
        r"""
        INPUT:

        - ``labels`` -- labels of the intervals
        - ``degrees`` -- degrees of the the cover

        Initialise the vector to zero.
        """
        self._vect = {n:{a:0 for a in labels} for n in degrees}
        self._labels = labels
        self._degrees = degrees
        self._iter = [(n, a) for a in labels for n in degrees]

    def __repr__(self):
        rep = ""
        for n in self._degrees:
            for a in self._labels:
                rep += "%s "%(self._vect[n][a])
            rep +=" | "
        rep += "\n"
        return rep

    def copy(self, x):
        r"""
        copy x in each component of the vector
        """
        for n, a in self._iter:
            self._vect[n][a] = deepcopy(x)
        return self

    def random(self):
        r"""
        Give random coordinates to each components of the vector.

        EXAMPLES::
            sage: VectPaths(['a'], range(1,3)).random()
            0.111439293741  | 0.514347513419  |
        """
        for n, a in self._iter:
            self._vect[n][a] = random()
        return self

    def var(self,n_var,a_var):
        r"""
        Put variables in (n,a) coordinate of the vector and zeros elsewhere. 
        """
        for n, a in self._iter:
            self._vect[n][a] = 0
        self._vect[n_var][a_var] = var("%s_%s"%(a_var,n_var))

    def id(self,n_id,a_id):
        r"""
        Put 1 in (n,a) coordinate of the vector and zeros elsewhere. 
        """
        for n, a in self._iter:
            self._vect[n][a] = 0
        self._vect[n_id][a_id] = 1
        return self

    def dot_product(self, other):
        r"""
        Return the dot_product of the vector with the other.
        """
        assert self._labels == other._labels and self._degrees == self._degrees
        return(sum([self._vect[n][a]*other._vect[n][a] for n in self._degrees for a in self._labels]))

    def norm(self):
        return (sqrt(self.dot_product(self)))

    def map(self,f):
        image = VectPaths(self._labels, self._degrees)
        for n, a in self._iter:
            image._vect[n][a] = f(n,a,self._vect[n][a])
        return image

    def map_mod(self,f):
        for n, a in self._iter:
            self._vect[n][a] = f(n,a,self._vect[n][a])

    def substract(self, v):
        image = VectPaths(self._labels, self._degrees)
        for n, a in self._iter:
            image._vect[n][a] = self._vect[n][a] - v._vect[n][a]
        return image

    def apply(self,f):
        for n, a in self._iter:
            f(n,a,self._vect[n][a])

    def val(self,n,a):
        return(self._vect[n][a])

    def copy_vector(self, v):
        r"""
        EXAMPLES::
            sage: v = VectPaths(['a','c','b'], range(1, 3))
            sage: v.copy_vector([1, 2, 3, 4, 5, 6])
            1 2 3  | 4 5 6  |
            <BLANKLINE>
        """    
        d = len(self._labels)
        for i in range(len(self._labels)):
            for n in self._degrees:
                self._vect[n][self._labels[i]] = v[d*(n-1) + i]
        return self

    def to_list(self):
        r"""
        Convert PathVector to vector
        
        EXAMPLES::
            sage: v = VectPaths(['a','c','b'], range(1, 3)).random()
            sage: v
            0.111439293741 0.514347513419 0.0446896852482  | 0.332490606442 0.520600243542 0.182713505069  | 
            <BLANKLINE>

            sage: v.to_list()
            [0.111439293741037, 0.5143475134191677, 0.04468968524815642, 0.332490606442413, 0.5206002435424241, 0.18271350506851058]

        """
        return [self._vect[n][a] for n in self._degrees for a in self._labels]

    def to_vector(self):
        return vector(CC, self.to_list())

    def normalise(self):
        norm = self.norm()
        if norm == 0:
            print "Warning: Norm of vector egal to 0"
            return 1
        else :
            for n, a in self._iter:
                self._vect[n][a] = self._vect[n][a]/norm
        return norm

    def isotopic_projection(self, chi, degree, g, perm):
        r"""
        Return a new VectPaths corresponding to the projection of self to the isotropic space of the \chi_i character
        """
        res = VectPaths(self._labels, self._degrees)
        for n, a in res._iter:
            res._vect[n][a] = projection(chi, degree, g, perm, self, n, a)
        return(res)
        
        
def vect_paths_ortho(v) :
    N = len(v)
    aux = deepcopy(v)
    for i in range(N):
        for j in range(i):
            v_j = v[j] 
            u_ij = aux[i]
            v_squared_norm = v_j.dot_product(v_j)
            proj = u_ij.dot_product(v_j)/v_squared_norm
            aux[i] = u_ij.map(lambda n,a,u: u - proj*v_j.val(n,a))
        v[i] = aux[i]
    return [v[i].normalise() for i in range(N)]


def vect_paths_ortho_householder(v):
    r"""
    Orthonormalise the vector basis v, and return the normalisation coefficient of vectors

    EXAMPLES::
        sage: v_1, v_2 = VectPaths(['a'], range(1,4)).random(), VectPaths(['a'], range(1,4)).random()
        sage: v = [v_1, v_2]
        sage: vp = deepcopy(v)
        sage: v
        [0.111439293741  | 0.514347513419  | 0.0446896852482  | 
        , 0.332490606442  | 0.520600243542  | 0.182713505069  | 
        ]

        sage: vect_paths_ortho_householder(v)
        [0.52817539578965744, 0.2525967350241719]

        sage: vect_paths_ortho(vp)
        [0.5281753957896576, 0.25259673502417185]

        sage: v
        [-0.210989180165  | -0.973819525709  | -0.0846114484022  | 
        , 0.821317669973  | -0.22354971821  | 0.524845509154  | 
        ]

        sage: v[0].dot_product(v[1])
        -3.4694469519536142e-16
    """
    N = len(v)
    if N == 0:
        return []
    l = len(v[0].to_list())
    M = array([v[i].to_list() for i in range(N)]).transpose()
    q, r = linalg.qr(M)
    for i in range(N):
        v[i].copy_vector(q[:l,i])
    return [abs(r[i, i]) for i in range(N)]
