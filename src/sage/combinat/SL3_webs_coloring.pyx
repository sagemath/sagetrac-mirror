from cpython cimport array
import array


cdef int Sign(lie, array.array Hl):
    """
    This function calculates the sign of the term corresponding to
    the web colouring corresponding to Hl.

    At each vertex of the web, we are supposed to multiply by -1 if
    the labels of the edges are -1,0,1 in counterclockwise order.

    For each vertex, lie contains a pair of consecutive edges, so we
    check whether the corresponding labels differ by 1 or -1.

    INPUT:

    Hl is an array encoding a web colouring.
    lie is a list which, for each vertex, gives a pair of edges.

    OUTPUT:

    1 or -1, which the term should be multiplied by.

    Note that the output is only correct up to an overall sign,
    """
    cdef int prod
    prod = 1
    for x, y in lie:
        if Hl[x] % 3 == (Hl[y] - 1) % 3:
            prod = -prod
    return prod


cdef tuple CI(array.array Hl, He, int n):
    """
    This function calculates the exponent of the monomial
    corresponding to the colouring encoded by Hl.

    OUTPUT:

    The tuple encoding the exponent.
    """
    cdef array.array exp
    cdef array.array template = array.array('b', [])
    exp = array.clone(template, 3 * n, True)
    for k in range(len(He)):
        if He[k][0] <= n and He[k][1] > n:
            exp[3 * (He[k][0] - 1) + Hl[k] + 1] += 1
    return tuple(exp)


cpdef Color(array.array Hl, Hn, He, int E, dict Sum, lie, int n):
    """
    Compute all possible colorings of the web.

    INPUT:

    Hn is a list, for each edge, of the edges adjacent to it.

    He is the list of edges.  We need it because we need to pass it to CI.

    The edges 0 to E-1 are assumed to have been already
    coloured; these colours are recorded in Hl.  The entries in Hl at E
    and after are arbitrary.

    n is the number of vertices

    lie is needed because we need to pass it to Sign.

    Sum is a dictionary associating exponent tuples to coefficients.
    This is where we keep track of the terms we have calculated.
    """
    cdef int c      # a possible colour we could add to edge E
    cdef tuple exp  # the exponent of a colouring
    cdef int val
    cdef int coeff

    # we look at what colours already appear on edges we have already labelled.
    NeighborsColors = []
    for y in Hn[E]:
        if E > y:
            NeighborsColors.append(Hl[y])

    # the other colours are the colours we could colour edge E with.
    PossibleColors = [x for x in [1, 0, -1] if x not in NeighborsColors]

    # we loop through those possibilities
    for c in PossibleColors:
        Hl[E] = c
        # if we still have more edges to colour, call color recursively
        if E + 1 < len(Hl):
            Color(Hl, Hn, He, E + 1, Sum, lie, n)
        # otherwise, we have built an entire colouring, and we add the
        # corresponding term into Sum
        else:
            exp = CI(Hl, He, n)
            val = Sign(lie, Hl)
        # we look up exp in the dictionary; if it is not there, its coefficient is zero
            try:
                coeff = Sum[exp]
            except KeyError:
                coeff = 0
            Sum[exp] = coeff + val
    return Sum
