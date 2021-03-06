# TODO:
# Turn things into iterators whenever possible
# tests tests tests tests

def positve_definite_by_orthogonal(genus,prime,e):
    r"""
    Returns all conjugacy classes of isometries of order ``prime^e``
    in ``genus``.

    The computation is using the enumeration of genera and
    computation of orthogonal groups.

    Input:

    - genus -- a positive definite genus
    - prime -- integer
    - e -- integer

    Output:

    - a
    """
    actions = []
    for R in genus.representatives():
        R = IntegralLattice(R)
        G = R.orthogonal_group()
        conj = [g for g in G.conjugacy_classes_representatives() if g.order()==prime^e]
        for f in conj:
            actions.append(LatticeWithIsometry(R,f.matrix()))
    return actions

def positve_definite_prime(genus, prime, e):
    r"""
    Returns all conjugacy classes of isometries of order ``prime^e``
    in ``genus``.

    The computation is using hermitian lattices and equivariant glueings.

    Input:

    - genus -- a positive definite genus
    - prime -- integer
    - e -- integer

    Output:



    """
    assert genus.signature_pair()[1] == 0, 'genus must be positive definite'
    rk = genus.rank()
    weights = [euler_phi(d) for d in (prime^e).divisors()]
    n = len(weights)
    floor = [0] + [0]*(n-2) + [1]
    ceiling = [rk // d for d in weights]
    P = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    P = [[p[k]*weights[k] for k in range(n)] for p in P]
    P = [p for p in P if sum(p)==rk]
    acts = []
    for ranks in P:
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        if prime == 2:
            signatures = [[0], [0]]
            signatures += [[0]*(weights[k]//2) for k in range(2,n)]
        else:
            signatures = [[0]]
            signatures += [[0]*(weights[k]//2) for k in range(1,n)]
        acts += prime_power_actions(genus,prime,ranks,signatures)
    results = []
    for act in acts:
        L = act[0]
        f = act[1]
        iso = LatticeWithIsometry(L,f)
        results.append(iso)
    return results


def check_prime_power(genus, prime, e):
    r"""
    Check that
    positive_definite_prime
    and
    positive_definite_by_orthogonal
    return the same number of results
    """
    expected = positve_definite_by_orthogonal(genus,prime,e)
    rep = positve_definite_prime(genus, prime, e)
    if len(rep) != len(expected):
        print(len(rep), expected)
        return False, rep, expected
    return True, rep, expected



def enriques_prime_power(genus, prime, e, verbose=False,rkT=None):
    r"""
    Returns all conjugacy classes of isometries of order `p^e` in `genus`
    such that the kernel of the `p^e`-part is of signature `(2,n)`.

    INPUT:

    - ``genus`` -- a genus of signature `(3,k)`.
    - ``p`` -- a prime number
    - ``e`` -- a natural number
    """
    assert genus.signature_pair()[0] == 2
    rk = genus.rank()
    weights = [euler_phi(d) for d in (prime^e).divisors()]
    n = len(weights)
    if prime == 2 and e==1:
        floor = [0, 2]     # minimal ranks  due to ample class and 2 form
    else:
        floor = [2] + [0]*(n-2) + [1]
    ceiling = [rk // d for d in weights]
    P = IntegerListsLex(max_sum=rk,length=n,floor=floor, ceiling=ceiling)
    P = [[p[k]*weights[k] for k in range(n)] for p in P]
    P = [p for p in P if sum(p)==rk]
    acts = []
    for ranks in P:
        if rkT is not None and ranks[-1]!=rkT:
            continue
        ranks_E = [ranks[k]//weights[k] for k in range(n)]
        if prime == 2 and e != 1:
            signatures = [[ranks_E[0]-2], [ranks_E[1]]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(2,n)]
        elif prime == 2 and e == 1:
            signatures = [[ranks_E[0]-2], [ranks_E[1]]]
        else:
            signatures = [[ranks_E[0]-2]]
            signatures += [[ranks_E[k]]*(weights[k]//2) for k in range(1,n)]
        for act in prime_power_actions(genus,prime,ranks,signatures,k3_unobstructed=True,verbose=verbose):
            yield act
