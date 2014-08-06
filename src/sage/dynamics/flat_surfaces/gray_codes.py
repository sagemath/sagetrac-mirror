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

- Knuth, "The art of computer programming", fascicules 2A (generating all
  n-tuples) and fascicule 3A (generating all combinations).

- Jenkins McCarthy, Ars Combinatorica, 40, (1995)

- W. H. Payne, ACM Trans. Math Software 5 (1979), 163--172.
"""
################################
# Switches for integer vectors #
################################

class TuplesGraySwitch:
    r"""
    Return the switch to perform for the enumeration of tuples of non negative
    integers up to a certain vector.

    This is an iterator over 2-tuples `(pos,incr)` which corresponds to a
    sequence of switches to perform in position ``pos`` by increment ``incr``.
    The increment ``incr`` is either ``+1`` or ``-1``.

    This is algorithm H in TAOCP 2a: loopless reflected mixed-radix Gray
    generation.

    EXAMPLES::

        sage: from sage.dynamics.flat_surfaces.separatrix_diagram import TuplesGraySwitch
        sage: l = [0,0,0]
        sage: for i,j in TuplesGraySwitch([3,3,3]):
        ....:     l[i] += j
        ....:     print l
        ....:     
        [1, 0, 0]
        [2, 0, 0]
        [2, 1, 0]
        [1, 1, 0]
        [0, 1, 0]
        [0, 2, 0]
        [1, 2, 0]
        [2, 2, 0]
        [2, 2, 1]
        [1, 2, 1]
        [0, 2, 1]
        [0, 1, 1]
        [1, 1, 1]
        [2, 1, 1]
        [2, 0, 1]
        [1, 0, 1]
        [0, 0, 1]
        [0, 0, 2]
        [1, 0, 2]
        [2, 0, 2]
        [2, 1, 2]
        [1, 1, 2]
        [0, 1, 2]
        [0, 2, 2]
        [1, 2, 2]
        [2, 2, 2]
        sage: l = [0,0]
        sage: for i,j in TuplesGraySwitch([2,1]):
        ....:     l[i] += j
        ....:     print l
        [1, 0]

    TESTS::

        sage: for t in [[2,2,2],[2,1,2],[3,2,1],[2,1,3]]:
        ....:     assert sum(1 for _ in TuplesGraySwitch(t)) == prod(t)-1
    """
    def __init__(self, m):
        self.m = []     # the limit (static)
        self.mm = []    # jump over entries with limit 1 that will not change
                        # during time (static)
        k = 0
        for i in m:
            if i <= 0:
                raise ValueError("accept only positive integers")
            if i == 1:
                k += 1
            else:
                self.m.append(i-1)
                self.mm.append(k)
        self.n = len(self.m)      # the actual number of elements on which the
                                  # algorithm perform the iteration (static)
        self.f = range(self.n+1)  # focus pointer (dynamic)
        self.o = [1]*self.n       # +1 or -1 for the switch (dynamic)
        self.a = [0]*self.n       # the current tuple (dynamic)
        
    def __iter__(self):
        r"""
        Return itself.

        TESTS::

            sage: from sage.dynamics.flat_surfaces.separatrix_diagramimport TuplesGraySwitch
            sage: G = TuplesGraySwitch([2,3,4])
            sage: iter(G) is G
            True
        """
        return self
         
    def next(self):
        r"""
        Return the next switch

        EXAMPLES::

            sage: from sage.dynamics.flat_surfaces.separatrix_diagram import TuplesGraySwitch
            sage: G = TuplesGraySwitch([2,2,2])
            sage: G.next()
            (0, 1)
            sage: G.next()
            (0, 1)
            sage: G.next()
            (1, 1)
        """
        j = self.f[0]
        if j == self.n: 
            raise StopIteration
        self.f[0] = 0            
        o = self.o[j]            
        self.a[j] += o
        if self.a[j] == 0 or self.a[j] == self.m[j]:
            self.f[j] = self.f[j+1]
            self.f[j+1] = j+1
            self.o[j] = -o
        return j+self.mm[j],o

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
