from interval_exchange import *
from path_vector import *
from sage.all import ceil, n

def lyap(ie, it) :
    t = 0
    N = ie.genus()
    if N == 0:
        return []
    v =  [VectPaths(ie.labels(), ie.degrees()).random() for i in xrange(N)]
    theta = [0 for i in xrange(N)]
    if ie._lengths == None:
        ie.rand_lg()
    ie.normalise()
    for i in xrange(it) :
        (A, B, c, perm_one, perm_two) = ie.rauzy()
        
        for j in xrange(N):                                                    #first apply the change of B
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for n, a in v[j]._iter:
                if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                    image._vect[n][a] = v[j].val(n,B) + c*v[j].val(perm_one(n),A)
                else:
                    image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        for j in xrange(N):                                                    #then apply the permutation
            image = VectPaths(v[j]._labels, v[j]._degrees)
            for n, a in v[j]._iter:
                if a == A:
                    image._vect[n][a] = v[j].val(perm_two(n),A)
                else:
                    image._vect[n][a] = v[j].val(n,a)
            v[j] = image
        
        if ie.diff_sum() > 2**(-precision+20):
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-6):
            nm = vect_paths_ortho_householder(v)
            t += -log(ie.normalise())
            for i in xrange(N):
                    theta[i] += log(nm[i])
    if t == 0:
        return [0 for i in xrange(N)]
    return [theta[i]/t for i in xrange(N)]


def lyap_isotopic_decomposed(ie, it) :
    t = 0
    isotopic_dimension = [ceil(ie.character_degree()[k]*ie.number_of_character_appeareance(k)/2) for k in xrange(ie.n_characters())]
    v_iter = [ (k, l) for k in xrange(ie.n_characters()) for l in xrange(isotopic_dimension[k])]
    v = [[ie.canonical_VectPaths().random() for l in xrange(isotopic_dimension[k])] for k in xrange(ie.n_characters())]
    theta = [[0 for l in xrange(isotopic_dimension[k])] for k in xrange(ie.n_characters())]

    if ie._lengths == None:
        ie.rand_lg()
    ie.normalise()

    def project():
        for k, l in v_iter:
            v[k][l] = ie.vector_isotopic_projection(k, v[k][l])
    for i in xrange(it) :
        project()

        (A, B, c, perm_one, perm_two) = ie.rauzy()
        
        for k, l in v_iter :
            w = v[k][l]                                                   #first apply the change of B
            image = VectPaths(w._labels, w._degrees)
            for n, a in w._iter:
                if a == B:                                                #return inversed and transposed matrix, thus  add (-1)
                    image._vect[n][a] = CC(w.val(n,B) + c*w.val(perm_one(n),A))
                else:
                    image._vect[n][a] = CC(w.val(n,a))
            v[k][l] = image
        for k,l in v_iter:                                                    #then apply the permutation
            w = v[k][l]
            image = VectPaths(w._labels, w._degrees)
            for n, a in w._iter:
                if a == A:
                    image._vect[n][a] = w.val(perm_two(n),A)
                else:
                    image._vect[n][a] = w.val(n,a)
            v[k][l] = image

        if ie.diff_sum() > diff_sum_seuil:
            t += -log(ie.normalise())
        if ie.min_length() < 2**(-10):
            t += -log(ie.normalise())
            for k in xrange(ie.n_characters()):
                nm = vect_paths_ortho_householder(v[k])
                for l in xrange(isotopic_dimension[k]):
                    theta[k][l] += log(nm[l])
    if t == 0:
        return [[0 for l in xrange(isotopic_dimension[k])] for k in xrange(ie.n_characters())]
    else:
        return [[theta[k][l]/t for l in xrange(isotopic_dimension[k])] for k in xrange(ie.n_characters())]
