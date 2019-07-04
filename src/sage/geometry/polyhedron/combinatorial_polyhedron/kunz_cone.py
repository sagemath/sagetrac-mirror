from sage.geometry.polyhedron.base import Polyhedron
from sage.rings.integer_ring import ZZ

def kunz_cone(m, backend='normaliz'):
    V = ZZ**(m-1)
    B = V.basis()
    #print B
    ieqs = [[0] +  list(B[i-1]+B[j-1]-B[(i+j-1) % m]) for j in range(1,m) for i in range(1,j+1) if i+j != m]
    #print len(ieqs)
    #return ieqs
    return Polyhedron(ieqs=ieqs, backend=backend)

def sort_matrix(M):
    a = M.rows()
    a.sort(key=sum, reverse=True)
    M1 = matrix(a)
    b = M1.columns()
    b.sort(key=sum, reverse=True)
    M2 = matrix(b)
    return M2.transpose()

