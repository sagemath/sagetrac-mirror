from interval_exchange import *

def orientable_double_cover(intervals):
    r"""
    Return the interval exchange of the orientated cover
    EXEMPLES ::
        sage: intervals = iet.GeneralizedPermutation('a a b', 'b c c')
        sage: orientable_double_cover(intervals)
        Intervals    : [ a, -a,  b]
                       [ b,  c, -c]
        Lengths      : None
        Permutations : {'a': [2, 1], 'c': [2, 1], 'b': [1, 2]}
        <BLANKLINE>

    """
    def twin(intervals, i, j):
        for n in range(len(intervals[i])):
            if n <> j and intervals[i][n] == intervals[i][j]:
                return(i,n)
        for n in range(len(intervals[(i+1)%2])):
            if intervals[(i+1)%2][n] == intervals[i][j]:
                return((i+1)%2,n)

    def orientable_double_cover_permutation(i, j):
        from sage.combinat.permutation import Permutation, Permutations
        (i_twin, j_twin) = twin(intervals, i, j)
        if i_twin <> i :
            return Permutations(2).identity()
        else :
            return Permutation('(1,2)')
    return CoveringStratum(intervals, permutations = {intervals[i][j]: orientable_double_cover_permutation(i, j) for i in range(2) for j in range(len(intervals[i]))})


def square_tiled_surface(horizontal_permutation, vertical_permutation):
    return CoveringStratum([['a','b'],['b','a']], permutations = {'a' : horizontal_permutation, 'b' : vertical_permutation})

def square_tiled_surface_pillow_case(perm_a, perm_b, perm_c):
    return CoveringStratum([['a','a','b'],['b','c','c']], permutations = {'a' : perm_a, 'b' : perm_b, 'c' : perm_c})

def EKZ_example(n):
    from sage.combinat.permutation import Permutation
    horizontal_permutation = Permutation ([(k + 1) % n + 1 for k in range(n)])
    vertical_permutation = Permutation ([n - k for k in range(n)])
    return square_tiled_surface(horizontal_permutation, vertical_permutation)

def pillow_case_iet(N):
    from sage.combinat.permutation import Permutation
    perm_a = Permutation ([ (k + N - 2) % (2*N) + 1 for k in range(2*N) ])
    perm_b = Permutation ([ (k + 0) % (2*N) + 1 for k in range(2*N) ])
    perm_c = Permutation ([ (k + N) % (2*N) + 1 for k in range(2*N) ])
    return square_tiled_surface_pillow_case(perm_a, perm_b, perm_c)

def pillow_case_exact(N):
    a = range(4)
    if N % 2 == 1:
        a[0] = N + 2
        a[1] = N - 2
        a[2] = N
        a[3] = N
    if N % 2 == 0:
        a[0] = N-1
        a[1] = 1
        a[2] = N-1
        a[3] = 1
    return cyclic_cover(2*N, a)

def cyclic_cover_iet(N, a):
    from sage.combinat.permutation import Permutation
    if (a[0] + a[1] + a[2] + a[3]) % N <> 0 :
        raise NameError('Not square-tiled cyclic cover')
    n_a = a[1]
    n_b = -a[0] - a[1]
    n_c = a[2]
    perm_a = Permutation ([ (k + n_a) % N + 1 for k in range(N) ])
    perm_b = Permutation ([ (k + n_b) % N + 1 for k in range(N) ])
    perm_c = Permutation ([ (k + n_c) % N + 1 for k in range(N) ])
    return square_tiled_surface_pillow_case(perm_a, perm_b, perm_c)

def cyclic_cover_exact(N, a):
    from sage.rings.integer import Integer
    res = range(N-1)
    for k in range(1, N):
        t = range(4)
        for i in range(4):
            t[i] = RR(a[i]/(N+ 0.)*k).frac()
        res[k-1] = Integer(round(sum(t)))
    results = [(0,0)]
    for k in range(1, N):
        dim = res[k-1] + res[N - k - 1] - 2
        results += [(res[N - k - 1] - 1 , res[k - 1] - 1)]
    singularities_order = []
    from math import gcd
    for i in range(4):
        singularities_order += [N / gcd(a[i], N) - 2 for k in range(gcd(a[i], N))]
    return results

