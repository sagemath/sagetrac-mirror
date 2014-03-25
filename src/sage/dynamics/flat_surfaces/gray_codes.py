r"""
Gray codes.

The file implements

1) A Gray code for k-subsets of {0,...,n-1}

For k-subsets of {0,...,n-1}, a Gray-code (called revolving door algorithm)  may
be obtained by

n = s + t
we want to generate all partitions of {0,...,n-1} into two sets of size
respectively s and t.

C(1,0) = 0 and C(1,1) = 1
C(n,k) = 0 C_(n-1)^k  u  1 rev(C_(n-1)^k-1)

It corresponds to algorithm R of Knuth (revolving-door combinations). It is a
genlex order for combinations (all sequences with the same prefix occurs
consecutively).

REFERENCES:

- Knuth, "The art of computer programming", fascicule 3A
Generating

- Jenkins McCarthy, Ars Combinatorica, 40, (1995)

- W. H. Payne, ACM Trans. Math Software 5 (1979), 163--172.
"""

###########################################
# The switches that generate combinations #
###########################################

def revolving_door_switch_iterator(n,t,verbose=False):
    r"""
    Iterator through the switches of the revolving-door algorithm.

    OUTPUT:

    - pairs (i,j) where item i is the one to unselect and j to select.

    EXAMPLES::

        sage: b = [1, 1, 1, 0, 0]
        sage: for i,j in revolving_door_switch_iterator(5,3):
        ...      b[i] = 0; b[j] = 1
        ...      print b
        [1, 0, 1, 1, 0]
        [0, 1, 1, 1, 0]
        [1, 1, 0, 1, 0]
        [1, 0, 0, 1, 1]
        [0, 1, 0, 1, 1]
        [0, 0, 1, 1, 1]
        [1, 0, 1, 0, 1]
        [0, 1, 1, 0, 1]
        [1, 1, 0, 0, 1]

        sage: b = [1,0,0,0]
        sage: for i,j in revolving_door_switch_iterator(4,1):
        ...       b[i] = 0; b[j] = 1
        ...       print b
        [0, 1, 0, 0]
        [0, 0, 1, 0]
        [0, 0, 0, 1]

        sage: list(revolving_door_switch_iterator(4,4))
        []
    """
    t = int(t)
    n = int(n)
    assert 0 <= t and t <= n, "Parameters not admissible"
    if t == 0 or t == n:
        return iter([])
    if t % 2:
        return _revolving_door_switch_iterator_odd(n,t,verbose)
    else:
        return _revolving_door_switch_iterator_even(n,t,verbose)

def _revolving_door_switch_iterator_odd(n,t,verbose=False):
    r"""
    Revolving door algorithm when t is odd.
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R2 visit
        if verbose: print ''.join(map(str,reversed(c[:-1])))

        # R3 : easy case
        if c[0] + 1 < c[1]:
            if verbose: print " R3"
            yield c[0], c[0]+1
            c[0] += 1
            continue

        j = 1
        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                if verbose: print " R4 with j = %d" %j
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                if verbose: print " R5 with j = %d" %j
                yield c[j-1], c[j]+1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break

def _revolving_door_switch_iterator_even(n,t,verbose=False):
    r"""
    Revolving door algorithm when t is even.
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R2 visit
        if verbose: print ''.join(map(str,reversed(c[:-1])))

        # R3 : easy case
        if c[0] > 0:
            if verbose: print " R3"
            yield c[0], c[0]-1
            c[0] -= 1
            continue

        j = 1
        # R5 : try to increase c[j]
        # at this point c[j-1] = j-1
        if c[j] + 1 < c[j+1]:
            if verbose: print " R5 with j = %d" %j
            yield c[j-1], c[j]+1
            c[j-1] = c[j]
            c[j] += 1
            continue
        j += 1

        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                if verbose: print " R4 with j = %d" %j
                yield c[j], j-1
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                if verbose: print " R5 with j = %d" %j
                yield c[j-1], c[j] + 1
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break

##############################
# Generation of combinations #
##############################

def revolving_door_iterator(n,t,verbose=False):
    r"""
    """
    t = int(t)
    n = int(n)
    assert 1 < t and t < n
    if t % 2:
        return _revolving_door_iterator_odd(n,t,verbose)
    else:
        return _revolving_door_iterator_even(n,t,verbose)

def _revolving_door_iterator_odd(n,t,verbose=False):
    r"""
    Revolving door algorithm when t is odd.
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R2 visit
        yield c[:t]
        if verbose: print ''.join(map(str,reversed(c[:-1])))

        # R3 : easy case
        if c[0] + 1 < c[1]:
            if verbose: print " R3"
            c[0] += 1
            continue

        j = 1
        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                if verbose: print " R4 with j = %d" %j
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                if verbose: print " R5 with j = %d" %j
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break

def _revolving_door_iterator_even(n,t,verbose=False):
    r"""
    Revolving door algorithm when t is even.
    """
    c = range(t) + [n]    # the combination (ordered list of numbers of length t+1)

    while True:
        # R2 visit
        yield c[:t]
        if verbose: print ''.join(map(str,reversed(c[:-1])))

        # R3 : easy case
        if c[0] > 0:
            if verbose: print " R3"
            c[0] -= 1
            continue

        j = 1
        # R5 : try to increase c[j]
        # at this point c[j-1] = j-1
        if c[j] + 1 < c[j+1]:
            if verbose: print " R5 with j = %d" %j
            c[j-1] = c[j]
            c[j] += 1
            continue
        j += 1

        while j < t:
            # R4 : try to decrease c[j]
            # at this point c[j] = c[j-1] + 1
            if c[j] > j:
                if verbose: print " R4 with j = %d" %j
                c[j] = c[j-1]
                c[j-1] = j-1
                break
            j += 1

            # R5 : try to increase c[j]
            # at this point c[j-1] = j-1
            if c[j] + 1 < c[j+1]:
                if verbose: print " R5 with j = %d" %j
                c[j-1] = c[j]
                c[j] += 1
                break
            j += 1

        else: # j == t
            break
