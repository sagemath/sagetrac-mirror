import enumerate_short_vectors_boost_python

def short_vectors(lattice, lower_bound, upper_bound, up_to_sign = False):
    r"""
    Enumerate vectors of minimal norm ``lower_bound`` and maximal norm
    ``upper_bound``, either up to sign or not.

    INPUT:
    
    - ``lattice`` - a list of list of integers that corresponds to the
                    Gram matrix of a binary quadratic form.

    - ``lower_bound`` - a positive integer.

    - ``upper_bound`` - a positive integer.

    - ``up_to_sign`` - a boolean (default: ``False``).

    OUPUT:

    - A dictionary mapping integers to a list of tuples.  Each tuple
      corresponds to a vector.

    EXAMPLES::

        sage: from sage.quadratic_forms.enumerate_short_vectors.enumerate_short_vectors_python import short_vectors
        sage: sorted(short_vectors([[2, 1], [1, 2]], 4, 10).items())
        [(4, []), (6, [(1, -2), (-1, 2), (-1, -1), (1, 1), (2, -1), (-2, 1)]), (8, [(0, -2), (0, 2), (2, -2), (-2, 2), (-2, 0), (2, 0)]), (10, [])]
    """
    return enumerate_short_vectors_boost_python \
        .enumerate_short_vectors_boost_python([[int(e) for e in r] for r in lattice],
                                              int(lower_bound), int(upper_bound),
                                              bool(up_to_sign))
