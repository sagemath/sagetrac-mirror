import time

def gauss(m, maxi=None):
    m = copy(m)
    n = m.nrows()
    if maxi is None:
        maxi=n
    for i in range(maxi):
        # Recherche du pivot
        for j in range(i,n):
            if m[j,i]:
                m.swap_rows(i,j)
                m.rescale_row(i, ~m[i,i])
                break
        if not m[i,i]:
            raise ValueError("non invertible matrix")

        for j in range(i+1, n):
            m.add_multiple_of_row(j, i, -m[j,i])
    return m

def matrice_inversible(n, corps=QQ):
    while True:
        m = random_matrix(corps, n, n)
        if m.det():
            return m

def hilbert(n):
    return matrix(QQ, n, lambda i,j: 1/QQ(i+j+1))

def temps(f, n, input=input):
    m = input(n)
    debut = time.time()
    f(m)
    t = time.time() - debut
    print n, t
    return t
