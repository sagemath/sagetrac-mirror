import time

def gauss(m):
    m = copy(m)
    n = m.nrows()
    for i in range(n):
        # Recherche du pivot
        for j in range(i,n):
            if m[j,i]:
                m.swap_rows(i,j)
                break
        if not m[i,i]:
            raise ValueError("non invertible matrix")

        for j in range(i+1, n):
            m.add_multiple_of_row(j, i, -m[j,i]/m[i,i])
    return m

def data(n, corps=QQ):
    # Renvoie une matrice inversible
    while True:
        m = random_matrix(corps, n, n)
        if m.det():
            return m

def temps(f, n, data=data, corps=QQ):
    print n
    m = data(n, corps)
    debut = time.time()
    f(m)
    return time.time() - debut
