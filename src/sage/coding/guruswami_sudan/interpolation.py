from sage.functions.other import ceil, binomial
from sage.matrix.constructor import matrix

def _flatten_once(lstlst):
    r"""
    Flattens``lstlst`` only once and returns a generator.

    This is similar to Python's ``flatten`` method, except that here, if you
    provide a list of lists of lists (and so on), it returns a list of lists
    (and so on) and not a list.

    INPUT:

    - ``lstlst`` -- a list of lists.

    EXAMPLES::

    sage: from sage.coding.guruswami_sudan.interpolation import _flatten_once
    sage: ll = [[1,2], [3,4], [5,6]]
    sage: _flatten_once(ll) #random
    <generator object _flatten_once at 0x7fc1631cca50>
    """
    for lst in lstlst:
        for e in lst:
            yield e

def _monomial_list(maxdeg, l, wy):
    r"""
    Returns a list of the `(x,y)` powers of all monomials in `F[x,y]` whose
    (1,wy)-weighted degree is less than ``maxdeg`` and whose ``y-degree <= l``.

    INPUT:

    - ``maxdeg``, ``l``, ``wy`` -- integers.

    OUTPUT:

    - a list of pairs of integers.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _monomial_list
        sage: _monomial_list(8, 5, 4)
        [(0, 0),
         (1, 0),
         (2, 0),
         (3, 0),
         (4, 0),
         (5, 0),
         (6, 0),
         (7, 0),
         (0, 1),
         (1, 1),
         (2, 1),
         (3, 1)]
    """
    monomials = []
    for y in range(0, l+1):
        for x in range(0,  ceil(maxdeg - y*wy)):
            monomials.append((x, y))
    return monomials

def _interpol_matrix_by_mons(points, s, monomials):
    r"""
    Returns a generator of the interpolation matrix whose nullspace gives the coefficients
    for all interpolation polynomials, given the list of monomials allowed.

    Its ``i``-th column will be the coefficients on the ``i``-th monomial
    in ``monomials``.

    INPUT:

    - ``points`` -- a list of integers, the interpolation points.

    - ``s`` -- an integer, the multiplicity parameter from Guruswami-Sudan algorithm.

    - ``monomials`` -- a list of monomials.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _interpol_matrix_by_mons
        sage: points = [(0, 2), (1, 5), (2, 0), (3, 4), (4, 9), (5, 1), (6, 9), (7, 10)]
        sage: s = 1
        sage: monomials = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (0, 1), (1, 1)]
        sage: _interpol_matrix_by_mons(points, s, monomials) #random
        <generator object _flatten_once at 0x7fb5ff8cce10>
    """
    n = len(points)
    def eqs_affine(x0,y0):
        r"""
        Make equation for the affine point x0, y0. Return a list of
        equations, each equation being a list of coefficients corresponding to
        the monomials in mons.
        """
        eqs = []
        for i in range(0, s):
            for j in range(0, s-i):
                eq = dict()
                for monomial in monomials:
                    ihat = monomial[0]
                    jhat = monomial[1]
                    if ihat >= i and jhat >= j:
                        icoeff = binomial(ihat, i)*x0**(ihat-i) \
                                    if ihat > i else 1
                        jcoeff = binomial(jhat, j)*(y0**(jhat-j)) \
                                    if jhat > j else 1
                        eq[monomial] = jcoeff*icoeff
                eqs.append([eq.get(monomial, 0) for monomial in monomials])
        return eqs
    return _flatten_once([ eqs_affine(*point) for point in points ])

def _interpol_matrix_problem(points, tau, parameters, wy):
    r"""
    Returns the linear system of equations which ``Q`` should be a solution to.

    This linear system is returned as a matrix ``M`` and a list of monomials ``monomials``,
    where a vector in the right nullspace of ``M`` corresponds to an
    interpolation polynomial `Q`, by the `i`'th element
    being the coefficient of the `i`'th monomial in ``monomials`` of `Q`.

    INPUT:

    - ``points`` -- a list of interpolation points.

    - ``tau`` -- an integer, the number of errors one wants to decode.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _interpol_matrix_problem
        sage: points = [(0, 2), (1, 5), (2, 0), (3, 4), (4, 9), (5, 1), (6, 9), (7, 10)]
        sage: tau = 1
        sage: params = (1, 1)
        sage: wy = 1
        sage: _interpol_matrix_problem(points, tau, params, wy)
        (
        [     1      0      0      0      0      0      0      2      0      0      0      0      0]
        [     1      1      1      1      1      1      1      5      5      5      5      5      5]
        [     1      2      4      8     16     32     64      0      0      0      0      0      0]
        [     1      3      9     27     81    243    729      4     12     36    108    324    972]
        [     1      4     16     64    256   1024   4096      9     36    144    576   2304   9216]
        [     1      5     25    125    625   3125  15625      1      5     25    125    625   3125]
        [     1      6     36    216   1296   7776  46656      9     54    324   1944  11664  69984]
        [     1      7     49    343   2401  16807 117649     10     70    490   3430  24010 168070], [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (6, 0), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1), (5, 1)]
        )
    """
    s, l = parameters[0], parameters[1]
    monomials = _monomial_list((len(points)-tau)*s, l, wy)
    M = matrix(list(_interpol_matrix_by_mons(points, s, monomials)))
    return (M, monomials)

def _construct_Q_from_matrix(M, monomials):
    r"""
    Returns a satisfactory ``Q`` polynomial given the interpolation matrix problem
    and the corresponding list of monomials.

    IMPUT:

    - ``M`` -- a matrix.

    - ``monomials`` -- a list of monomials.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import _construct_Q_from_matrix
        sage: from sage.coding.guruswami_sudan.interpolation import _interpol_matrix_problem
        sage: points = [(0, 2), (1, 5), (2, 0), (3, 4), (4, 9), (5, 1), (6, 9), (7, 10)]
        sage: tau = 1
        sage: params = (1, 1)
        sage: wy = 1
        sage: res = _interpol_matrix_problem(points, tau, params, wy)
        sage: M, monomials = res[0], res[1]
        sage: _construct_Q_from_matrix(M, monomials)
        4202026*x^6 - 5614235*x^5*y - 29351399*x^5 + 64635986*x^4*y + 41894587*x^4 - 273534229*x^3*y + 3*x^3 + 508264978*x^2*y + x^2 - 297101040*x*y + 2*x - 840*y + 1680
    """
    if M.nrows() >= M.ncols():
        raise Exception("More rows than columns! This matrix is not satisfactory.")
    Sp = M.right_kernel()
    sol = Sp.an_element()
    #TODO: Option to pick out minimal element?
    while sol.is_zero():
        # Picking out e.g. element 1 directly seems to run into an infinite
        # loop for large matrices.
        sol = Sp.random_element()
    # Construct the Q polynomial
    PF = M.base_ring()['x', 'y'] #make that ring a ring in <x>
    x, y = PF.gens()
    Q = sum([x**monomials[i][0] * y**monomials[i][1] * sol[i] for i in range(0, len(monomials))])
    return Q

def construct_Q_linalg(points, tau, parameters, wy):
    r"""
    Returns an interpolation polynomial Q(x,y) for the given input.

    INPUT:

    - ``points`` -- a list of tuples ``(xi, yi)`` such that
    ``Q(xi,yi) = 0 ``with multiplicity ``s``.

    - ``tau`` -- an integer, the number of errors one wants to decode.

    - ``parameters`` -- (default: ``None``) a pair of integers, where:
        - the first integer is the multiplicity parameter of Guruswami-Sudan algorithm and
        - the second integer is the list size parameter.

    - ``wy`` -- an integer.

    EXAMPLES::

        sage: from sage.coding.guruswami_sudan.interpolation import construct_Q_linalg
        sage: points = [(0, 2), (1, 5), (2, 0), (3, 4), (4, 9), (5, 1), (6, 9), (7, 10)]
        sage: tau = 1
        sage: params = (1, 1)
        sage: wy = 1
        sage: construct_Q_linalg(points, tau, params, wy)
        4202026*x^6 - 5614235*x^5*y - 29351399*x^5 + 64635986*x^4*y + 41894587*x^4 - 273534229*x^3*y + 3*x^3 + 508264978*x^2*y + x^2 - 297101040*x*y + 2*x - 840*y + 1680
    """
    return _construct_Q_from_matrix(
                *_interpol_matrix_problem(points, tau, parameters, wy))
