"""
Fourier Coefficients

::

    sage: import sys
    sage: sys.path.append(".")

"""
import itertools

from sage.combinat.finite_state_machine import Transducer
from sage.functions.log import log
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc import srange, verbose
from sage.modules.free_module_element import vector
import sage.rings.arith
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_mpfi import min_RIF, max_RIF
from sage.structure.sage_object import SageObject


def infinity_vector_norm(v):
    """
    Compute the infinity norm of a vector of
    :class:`~sage.rings.real_interval` elements.

    INPUT:

    - ``v`` -- a vector of :class:`~sage.rings.real_interval` elements.

    OUTPUT:

    A :class:`~sage.rings.real_interval` element.

    EXAMPLES::

        sage: from fsm_fourier import infinity_vector_norm
        sage: a = RIF(1, 4)
        sage: b = RIF(-3, -2)
        sage: infinity_vector_norm(vector([a, b])).endpoints()
        (2.00000000000000, 4.00000000000000)
        sage: infinity_vector_norm(vector([b, a])).endpoints()
        (2.00000000000000, 4.00000000000000)

    Note that the
    :meth:`sage.rings.modules.free_module_element.FreeModuleElement.norm`
    method is inadequate, as it uses the Python max instead of
    :func:`max_RIF`::

        sage: vector([a, b]).norm(infinity).endpoints()
        (1.00000000000000, 4.00000000000000)
        sage: vector([b, a]).norm(infinity).endpoints()
        (2.00000000000000, 3.00000000000000)
    """
    return max_RIF(abs(r) for r in v)


def infinity_matrix_norm(A):
    """
    Compute the infinity norm of a matrix A, in the same ring as
    the original matrix A.

    INPUT:

    - ``A`` -- a matrix.

    OUTPUT:

    An entry.

    In contrast to :meth:`sage.matrix.matrix2.Matrix.norm`,
    sparse zero matrices are allowed and the result is not
    necessarily a
    :class:`~sage.rings.real_double.RealDoubleElement`.

    EXAMPLES::

        sage: from fsm_fourier import infinity_matrix_norm
        sage: M = matrix([[1, 2], [-5, -1]])
        sage: infinity_matrix_norm(M)
        6
        sage: M.norm(infinity)
        6.0
        sage: M = matrix([[0]], sparse=True)
        sage: infinity_matrix_norm(M)
        0
        sage: M.norm(infinity)
        Traceback (most recent call last):
        ...
        TypeError: base_ring (=Category of objects) must be a ring
    """
    if A.is_zero():
        return 0
    return max(sum(r) for r in A.apply_map(abs).rows())


def _hurwitz_zeta_(s, alpha,  m = 0):
    r"""
    Compute the truncated Hurwitz zeta function `\sum_{k\ge m} (k+\alpha)^{-s}`.

    INPUT:

    -   ``s`` -- a :class:`ComplexIntervalField` element

    -   ``alpha`` -- a :class:`RealIntervalField` element

    -   ``m`` -- a positive integer


    OUTPUT:

    A :class:`ComplexIntervalField` element in the same ring as ``s``.

    EXAMPLES:

    -   ::

            sage: from fsm_fourier import _hurwitz_zeta_
            sage: _hurwitz_zeta_(CIF(2), RIF(3/4), 10)
            0.097483848201852? + 0.?e-19*I

    -   Compare with well-known value `\zeta(2)=\zeta(2, 1)=\pi^2/6`::

            sage: _hurwitz_zeta_(CIF(2), 1)
            1.64493406684823? + 0.?e-17*I
            sage: (_hurwitz_zeta_(CIF(2), 1) - CIF(pi)^2/6).abs()<10^(-13)
            True
            sage: _hurwitz_zeta_(CIF(2*pi*I/log(2)), 1)
            1.598734526809? + 0.278338669639?*I

    -   There is a singularity at `s=1`. ::

            sage: _hurwitz_zeta_(CIF(RIF(0.9, 1.1), (-0.1, 0.1)), 1)
            Traceback (most recent call last):
            ...
            ZeroDivisionError: zeta is singular at 1.

    -   Checking the function for non-positive integers::

            sage: zeta(-1) in _hurwitz_zeta_(CIF(-1), 1)
            True
            sage: zeta(0) in _hurwitz_zeta_(CIF(0), 1)
            True

    -   Debugging output can be enabled using
        :func:`~sage.misc.misc.set_verbose`. To test it, We use a
        large imaginary part because convergence is worse in those
        cases::

            sage: set_verbose(2)
            sage: _hurwitz_zeta_(CIF(1+100/log(2)*I), 1)
            verbose 1 (...) _hurwitz_zeta_(1 + 144.2695040888963?*I, 1, 0): M = 172
            verbose 2 (...)     N = 2, error = 0.0352354068797?, acceptable_error = 2.2204460492503131?e-16, result = 2.125571548789? + 0.511221280470?*I
            verbose 2 (...)     N = 4, error = 0.000310532577681?, acceptable_error = 2.2204460492503131?e-16, result = 2.125575595864? + 0.51121897538?*I
            verbose 2 (...)     N = 6, error = 3.65215306101?e-6, acceptable_error = 2.2204460492503131?e-16, result = 2.125575660430? + 0.511218933060?*I
            verbose 2 (...)     N = 8, error = 4.83820940904?e-8, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661484? + 0.511218932225?*I
            verbose 2 (...)     N = 10, error = 6.84787778802?e-10, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.511218932208?*I
            verbose 2 (...)     N = 12, error = 1.011644165731?e-11, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.51121893221?*I
            verbose 2 (...)     N = 14, error = 1.54088223831?e-13, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.51121893221?*I
            verbose 2 (...)     N = 16, error = 2.40250320318?e-15, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.51121893221?*I
            verbose 2 (...)     N = 18, error = 3.81730778323?e-17, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.51121893221?*I
            verbose 1 (...)     N = 18, error = 3.81730778323?e-17, acceptable_error = 2.2204460492503131?e-16, result = 2.125575661501? + 0.51121893221?*I
            2.125575661501? + 0.51121893221?*I
            sage: set_verbose(0)

    -   The current implementation does not work well with negative real
        values, all precision is lost::

            sage: _hurwitz_zeta_(CIF(-15+I), 1)
            0.?e12 + 0.?e12*I
            sage: hurwitz_zeta(ComplexField(200)(-15 + I), 1)
            0.66621329305522618549073659441004805750538410627754288912677
            - 0.84614995218731390314834131849322502368334938802936485299779*I

        A work-around is to start with higher precision; however, we have
        to clear the cache first::

            sage: _hurwitz_zeta_(ComplexIntervalField(200)(-15 + I), 1)
            0.66621329305522618549073660?
            - 0.8461499521873139031483414?*I
    """
    from sage.rings.arith import bernoulli, falling_factorial

    CIF = s.parent()
    RIF = s.real().parent()

    if ZZ(1) in s:
        raise ZeroDivisionError("zeta is singular at 1.")

    # We rely on (2pi)^-N for convergence of the error term.
    # As a conservative estimate, 2pi is approximately 2^2,
    # so we will finally need N ~ s.prec()/2 to achieve an error which
    # is less than the resolution of s.
    # In order to have the falling factorial (-s)^\underline{N}
    # smaller than (M+a)^N, we choose M>|s|+N
    M = max(m, (s.prec()/ZZ(2)).ceil() + ZZ(s.abs().upper().ceil()))
    verbose("_hurwitz_zeta_(%s, %s, %s): M = %d" % (s, alpha, m, M),
            level=1)

    sigma = s.real()
    result = sum((r + alpha)**(-s) for r in reversed(srange(m, M)))
    result += (M + alpha)**(1-s) / (s-1)
    factor = (M + alpha)**(-s)
    result += factor/2

    N = 0
    error_factor = RIF(4)
    N_factorial = ZZ(1)

    while True:
        N += 2
        factor *= (-s - N + 2)/(M + alpha)
        #assert factor.overlaps(falling_factorial(-s, N - 1)/(M + alpha)**(s + N - 1))
        N_factorial *= (N - 1)*N
        #assert ZZ(N).factorial() == N_factorial
        result -= bernoulli(N)/N_factorial * factor

        factor *= (-s - N + 1)
        #assert factor.overlaps(falling_factorial(-s, N)/(M + alpha)**(s + N - 1))

        error_factor /= 4*RIF.pi()**2
        #assert error_factor.overlaps(4/(2*RIF.pi())**N)
        error_bound = error_factor / (sigma + N - 1) * factor.abs()

        if result.abs().upper().is_zero():
            error_acceptable = 0
        else:
            error_acceptable = RIF(2) ** (max(result.real().abs().upper().log2(),
                                             result.imag().abs().upper().log2()).floor()
                                         - result.prec())

        verbose("    N = %d, error = %s, acceptable_error = %s, result = %s" %
                (N, error_bound, error_acceptable, result), level=2)

        if error_bound.abs() < error_acceptable:
            error_real = RIF(-error_bound, error_bound)
            error = CIF(error_real, error_real)
            verbose("    N = %d, error = %s, acceptable_error = %s, result = %s" %
                    (N, error_bound, error_acceptable, result), level=1)
            return result + error

        factor /= (M + alpha)
        #assert factor.overlaps(falling_factorial(-s, N)/(M + alpha)**(s + N))

class FSM_Fourier_Component(SageObject):
    """Hold a final component and associated data."""
    def __init__(self, fsm, parent):
        self.fsm = fsm
        self.period = fsm.graph().period()
        self.n_states = len(self.fsm.states())
        self.parent = parent

    def mask(self, n):
        nrows = sum(c.n_states
                    for c in self.parent.components if c != self)
        mask = matrix(
            nrows, n,
            [self.parent.standard_basis[
                    self.parent.positions[state.label()]]
             for other in self.parent.components
             if other != self
             for state in other.fsm.iter_states()])
        return mask

    def eigenvectors(self, M):
        mask = self.mask(self.parent.M.nrows())
        def eigenvector(j):
            eigenvalue = self.parent.q * self.parent.alpha**(
                j * self.parent.common_period / self.period)
            S = matrix.block(
                [[M - eigenvalue*matrix.identity(M.nrows())],
                 [mask]],
                subdivide=False)
            kernel = S.right_kernel_matrix()
            assert kernel.nrows() == 1
            if j == 0:
                #normalize for positive eigenvector
                return kernel.row(0) / sum(kernel.row(0))
            else:
                return kernel.row(0)

        return [eigenvector(j) for j in range(self.period)]

    @cached_method()
    def right_eigenvectors(self):
        return self.eigenvectors(self.parent.M)

    @cached_method()
    def left_eigenvectors(self):
        left_eigenvectors = self.eigenvectors(self.parent.M.transpose())
        return [w/(v*w) for v, w
                in itertools.izip(self.right_eigenvectors(),
                                  left_eigenvectors)]

    @cached_method()
    def vectors_w(self):
        return [(self.parent.initial_vector*v)*w for v, w
                in itertools.izip(self.right_eigenvectors(),
                                  self.left_eigenvectors())]

    @cached_method()
    def coefficient_lambda(self):
        products = [w*self.parent.ones for w in self.vectors_w()]
        assert all(e.is_zero() for e in products[1:])
        return products[0]

    @cached_method()
    def mu_prime(self):
        Y = self.parent.Y
        p = self.fsm.adjacency_matrix(
            entry=lambda t:Y**sum(t.word_out)).charpoly('Z')
        Z = p.parent().gen()
        assert p(Y=1, Z=self.parent.q) == 0
        mu_prime_Z = (- p.derivative(Y)/p.derivative(Z))(
            Y=1, Z=self.parent.q)
        return self.parent.I*mu_prime_Z

    @cached_method()
    def a(self):
        return QQ(-self.parent.I * self.mu_prime()/self.parent.q)

    def w_ell(self, ell):
        if self.parent.common_period.divides(ell*self.period):
            k = self.period*ell/self.parent.common_period % self.period
            return vector(self.parent.field_to_CIF(c) for c in self.vectors_w()[k])
        else:
            return vector(0 for _ in self.parent.ones)

    def vector_v_prime(self, k):
        M = self.parent.M
        mask = self.mask(M.nrows())
        eigenvalue = self.parent.q * self.parent.alpha**(
                k * self.parent.common_period / self.period)
        S = matrix.block(
             [[M - eigenvalue*matrix.identity(M.nrows())],
              [matrix(self.parent.ones)],
              [mask]],
             subdivide=False)
        eigenvector_right = vector(self.right_eigenvectors()[k])
        M_prime = self.parent.I*self.parent.Delta
        right_side = - matrix.block(
                                    [[M_prime - self.mu_prime()*matrix.identity(M.nrows())],
                                     [0*matrix(self.parent.ones)],
                                     [0*mask]],
                                    subdivide=False) * eigenvector_right
        v_prime = S.solve_right(right_side)
        return v_prime

    def vector_w_prime(self, k):
        M = self.parent.M
        mask = self.mask(M.nrows())
        eigenvalue = self.parent.q * self.parent.alpha**(
                k * self.parent.common_period / self.period)
        eigenvector_right = vector(self.right_eigenvectors()[k])
        eigenvector_left = vector(self.left_eigenvectors()[k])
        M_prime = self.parent.I*self.parent.Delta
        S = matrix.block(
             [[M.transpose() - eigenvalue*matrix.identity(M.nrows())],
              [matrix(eigenvector_right)],
              [mask]],
             subdivide=False)
        right_side = - matrix.block(
                                    [[M_prime.transpose() - self.mu_prime()*matrix.identity(M.nrows())],
                                     [matrix(self.vector_v_prime(k))],
                                     [0*mask]],
                                    subdivide=False) * eigenvector_left
        left_prime = S.solve_right(right_side)
        w_prime = self.parent.initial_vector*self.vector_v_prime(k)*eigenvector_left \
            + self.parent.initial_vector*eigenvector_right*left_prime
        return w_prime


class FSMFourier(SageObject):
    """
    Fourier coefficients for the sum of output of transducers.
    """

    def __init__(self, transducer):
        r"""
        Return the common data needed for the computation of all
        Fourier coefficients of the periodic fluctuation of the sum of
        output.

        INPUT:

        Nothing.

        OUTPUT:

        A :class:`namedtuple` consisting of:

        - ``c`` -- number of final components.

        - ``periods`` -- list of periods of the final components.

        - ``period`` -- least common multiple of the periods.

        - ``T`` -- eigenvector matrix.

        - ``w`` -- list of lists of vectors `\mathbf{w}_{jk}`.

        - ``coefficient_lambda`` -- list of coefficients `\lambda_j`.

        - ``e_T`` -- constant `e_{\mathcal{T}}`, the coefficient of
          the main term of the expectation.

        - ``a`` -- list of constants `a_j`.

        - ``M`` -- adjacency matrix `M`.

        - ``M_epsilon`` -- list of partial adjacency matrices
          `M_\varepsilon`.

        - ``Delta`` -- output matrix `\Delta`.

        - ``Delta_epsilon`` -- list of partial output matrices
          `\Delta_\varepsilon`.

        - ``C_0`` -- `\max\{\|\mathbf{b}(r)\|_\infty: 0\le r<q\}`.

        - ``C_1`` -- `\max\{\|\Delta_\varepsilon\|_\infty: 0\le \varepsilon<q\}`.

        - ``components`` -- a list of :class:`FSM_Fourier_Component`, representing the final components.

        EXAMPLES:

        -   Binary sum of digits::

                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: function('f')
                f
                sage: var('n')
                n
                sage: from fsm_fourier import FSMFourier
                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [1]
                sage: F.common_period
                1
                sage: F.T
                [1]
                sage: F.w
                [[(1)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                1/2
                sage: F.a
                [1/2]
                sage: F.M
                [2]
                sage: F.M_epsilon
                [[1], [1]]
                sage: F.Delta
                [1]
                sage: F.Delta_epsilon
                [[0], [1]]
                sage: F.C_0
                1
                sage: F.C_1
                1

        -   NAF::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(4*n + 1) == f(n) + 1,
                ....:     f(4*n + 3) == f(n + 1) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [1]
                sage: F.common_period
                1
                sage: F.T
                [1/3   1   0]
                [1/3   0   1]
                [1/3  -1  -1]
                sage: F.w
                [[(1/3, 1/3, 1/3)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                1/3
                sage: F.a
                [1/3]
                sage: F.FourierCoefficient(0) # long time
                0.4478541793943? + 0.?e-16*I
                sage: F.FourierCoefficient(42) # long time
                -9.0701442?e-6 + 0.0001189561459?*I

        -   Abelian complexity of the paperfolding sequence::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(4*n) == f(2*n),
                ....:     f(4*n+2) == f(2*n+1)+1,
                ....:     f(16*n+1) == f(8*n+1),
                ....:     f(16*n+5) == f(4*n+1)+2,
                ....:     f(16*n+11) == f(4*n+3)+2,
                ....:     f(16*n+15) == f(2*n+2)+1,
                ....:     f(1) == 2, f(0) == 0]
                ....:     + [f(16*n+jj) == f(2*n+1)+2 for jj in [3,7,9,13]],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [1]
                sage: F.common_period
                1
                sage: F.T
                [1/10    1    0    0    0    0    0    0    0    0]
                [1/10    0    1    0    0    0    0    0    0    0]
                [1/10    0    0    1    0    0    0    0    0    0]
                [1/10    0    0    0    1    0    0    0    0    0]
                [1/10    0    0    0    0    1    0    0    0    0]
                [1/10    0    0    0    0    0    1    0    0    0]
                [1/10    0    0    0    0    0    0    1    0    0]
                [1/10    0    0    0    0    0    0    0    1    0]
                [1/10    0    0    0    0    0    0    0    0    1]
                [1/10    0    0   -3   -2   -2   -2   -1   -1   -1]
                sage: F.w
                [[(0, 0, 3/13, 2/13, 2/13, 2/13, 1/13, 1/13, 1/13, 1/13)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                8/13
                sage: F.a
                [8/13]
                sage: F.FourierCoefficient(0) # long time
                1.5308151287593? + 0.?e-15*I
                sage: F.FourierCoefficient(42) # long time
                9.380157?e-7 + 0.0001569568848?*I

        -   Artificial example, one-periodic, 2 states::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(4*n) == f(2*n)+0,
                ....:     f(4*n+2) == f(n)+1,
                ....:     f(2*n+1) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [1]
                sage: F.common_period
                1
                sage: F.T
                [1/2   1]
                [1/2  -1]
                sage: F.w
                [[(1/2, 1/2)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                1/4
                sage: F.a
                [1/4]
                sage: F.FourierCoefficient(0) # long time
                -0.20319811602320? + 0.?e-16*I
                sage: F.FourierCoefficient(42) # long time
                -0.0000280287200? + 0.0000741203739?*I

        -   Artificial example, period 3::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(8*n) == f(4*n+3)+3,
                ....:     f(8*n+4) == f(4*n+3)+1,
                ....:     f(8*n+2) == f(4*n+3)+2,
                ....:     f(8*n+6) == f(4*n+3)-1,
                ....:     f(8*n+1) == f(4*n)+5,
                ....:     f(8*n+5) == f(4*n+2)+1,
                ....:     f(8*n+3) == f(4*n+1)+2,
                ....:     f(8*n+7) == f(4*n+1),
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [3]
                sage: F.common_period
                3
                sage: F.T
                [         1/7               1               1            1            0            0            0]
                [         1/7  4*zeta12^2 - 4     -4*zeta12^2            0            1            0            0]
                [         1/7 -2*zeta12^2 + 2      2*zeta12^2            0            0            1            0]
                [         1/7     -4*zeta12^2  4*zeta12^2 - 4            0            0            0            1]
                [         1/7     -4*zeta12^2  4*zeta12^2 - 4            0            0            0           -1]
                [         1/7  4*zeta12^2 - 4     -4*zeta12^2            0            0            0            0]
                [         1/7               4               4            0            0            0            0]
                sage: F.w
                [[(0, 0, 0, 1/6, 1/6, 1/3, 1/3),
                (0, 0, 0, 1/24*zeta12^2 - 1/24, 1/24*zeta12^2 - 1/24, -1/12*zeta12^2, 1/12),
                (0, 0, 0, -1/24*zeta12^2, -1/24*zeta12^2, 1/12*zeta12^2 - 1/12, 1/12)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                7/4
                sage: F.a
                [7/4]
                sage: F.FourierCoefficient(0) # long time
                0.8794486773901? + 0.?e-14*I
                sage: F.FourierCoefficient(42) # long time
                0.000330400039? + 0.002040223215?*I
                sage: F.FourierCoefficient(43) # long time
                0.0001789628369? + 3.8774734?e-6*I

        -   Artificial example, period 2, vanishing w-vector::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(4*n) == f(2*n+1)+1,
                ....:     f(4*n+1) == f(2*n)+2,
                ....:     f(4*n+2) == f(2*n+1)+3,
                ....:     f(4*n+3) == f(2*n)-1,
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: F.c
                1
                sage: F.periods
                [2]
                sage: F.common_period
                2
                sage: F.T
                [1/3   0   1]
                [1/3   1   0]
                [1/3  -1   0]
                sage: F.w
                [[(0, 1/2, 1/2), (0, 0, 0)]]
                sage: F.coefficient_lambda
                [1]
                sage: F.e_T
                5/4
                sage: F.a
                [5/4]
                sage: F.FourierCoefficient(0) # long time
                -1.5912428334793? + 0.?e-15*I
                sage: F.FourierCoefficient(42) # long time
                -0.0002650957054? - 0.0001002936380?*I
                sage: F.FourierCoefficient(43) # long time
                0

        -   Artificial example with two final components of periods `2`
            and `3`, respectively::

                sage: F = FSMFourier(Transducer([(0, 1, 0, 1), (1, 2, 0, 1),
                ....:     (2, 1, 0, 2), (1, 2, 1, 0), (2, 1, 1, 2),
                ....:     (0, -1, 1, 1), (-1, -2, 1, 1), (-2, -3, 1, 1),
                ....:     (-3, -1, 1, 1), (-1, -2, 0, 2), (-2, -3, 0, 1),
                ....:     (-3, -1, 0, 3)],
                ....:     initial_states=[0],
                ....:     final_states=[0, 1, 2, -3, -2, -1]))
                sage: F.c
                2
                sage: F.periods
                [3, 2]
                sage: F.common_period
                6
                sage: F.T
                [        1/7              1              1         1/5           1           1]
                [          0              0              0         2/5          -2           0]
                [          0              0              0         2/5           2           0]
                [        2/7              2              2           0           0           0]
                [        2/7    -2*zeta12^2 2*zeta12^2 - 2           0           0           0]
                [        2/7 2*zeta12^2 - 2    -2*zeta12^2           0           0           0]
                sage: F.w
                [[(0, 0, 0, 1/6, 1/6, 1/6),
                  (0, 0, 0, 1/6, 1/6*zeta12^2 - 1/6, -1/6*zeta12^2),
                  (0, 0, 0, 1/6, -1/6*zeta12^2, 1/6*zeta12^2 - 1/6)],
                 [(0, 1/4, 1/4, 0, 0, 0), (0, -1/4, 1/4, 0, 0, 0)]]
                sage: F.coefficient_lambda
                [1/2, 1/2]
                sage: F.e_T
                11/8
                sage: F.a
                [5/4, 5/4]
                sage: F.FourierCoefficient(0) # long time
                -2.1863500631078? + 0.?e-14*I
                sage: F.FourierCoefficient(42) # long time
                -0.0024222079519? + 0.0010332514417?*I
                sage: F.FourierCoefficient(43) # long time
                0
                sage: F.FourierCoefficient(44) # long time
                0.0002142083798? - 0.0004539668075?*I
                sage: F.FourierCoefficient(45) # long time
                -0.0005573382643? - 0.0005531638345?*I

        -   Ternary sum of digits::

                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(3*n + 2) == f(n) + 2,
                ....:     f(3*n + 1) == f(n) + 1,
                ....:     f(3*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 3))
                sage: F.common_period
                1
                sage: F.e_T
                1
                sage: F.FourierCoefficient(0) # long time
                -0.2373314270632? + 0.?e-17*I
                sage: F.FourierCoefficient(42) # long time
                0.0001516409849? + 0.0000541593062?*I
        """


        import collections
        import operator

        from sage.calculus.var import var
        from sage.functions.log import exp
        from sage.modules.free_module import VectorSpace
        from sage.rings.arith import lcm
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

        self.transducer = transducer
        self.positions = dict((state.label(), j)
                         for j, state in enumerate(transducer.iter_states()))
        self.q = len(transducer.input_alphabet)
        self.Y = PolynomialRing(QQ, 'Y').gen()
        self.ones = vector(1 for _ in transducer.iter_states())



        self.components = [FSM_Fourier_Component(c, self)
                           for c in transducer.final_components()]
        self.common_period = lcm([c.period
                                  for c in self.components])
        self.field = CyclotomicField(lcm(4, self.common_period))
        self.alpha = self.field.zeta(self.common_period)
        self.I = self.field.zeta(4)
        self.field_to_CIF = self.field.hom(
            [ComplexIntervalField().zeta(lcm(4, self.common_period))], check=False)
        assert self.I**2 == -1
        assert self.alpha**self.common_period == 1
        assert all(self.alpha**j != 1 for j in range(1, self.common_period))
        assert self.field_to_CIF(self.I).overlaps(ComplexIntervalField()(0, 1))
        assert self.field_to_CIF(self.alpha).overlaps(ComplexIntervalField()(exp(2*ComplexIntervalField().pi()*ComplexIntervalField()(0, 1)/self.common_period)))
        self.M = transducer.adjacency_matrix(entry=lambda t: 1)
        self.standard_basis = VectorSpace(self.field, self.M.nrows()).basis()

        if len(transducer.initial_states()) != 1:
            raise NotImplementedError(
                "Transducer does not have a unique initial state.")
        self.initial_vector = self.standard_basis[self.positions[
                transducer.initial_states()[0].label()]]

        right_eigenvectors = list(itertools.chain(
                *(c.right_eigenvectors()
                  for c in self.components)))

        left_eigenvectors = list(itertools.chain(
                *(c.left_eigenvectors()
                  for c in self.components)))

        annihilated_by_left = matrix(left_eigenvectors).\
            right_kernel_matrix().transpose()

        self.T = matrix.block([[matrix.column(right_eigenvectors),
                                annihilated_by_left]],
                              subdivide=False)

        assert self.T.is_square()
        assert self.T.nrows() == self.M.nrows()
        assert self.T.is_invertible()

        check = self.T.inverse() * self.M * self.T
        eigenvalues = [self.q * self.alpha**(j * self.common_period/c.period)
                       for c in self.components
                       for j in range(c.period)]
        check_dont_care = check.submatrix(len(eigenvalues),
                                          len(eigenvalues))
        assert (matrix.block(
                [[matrix.diagonal(eigenvalues), ZZ(0)],
                 [ZZ(0), check_dont_care]],
                subdivide=False) - check).is_zero()

        assert (self.T.inverse().submatrix(nrows=len(left_eigenvectors))
                - matrix(left_eigenvectors)).is_zero()

        self.e_T = sum(c.a()*c.coefficient_lambda()
                       for c in self.components)

        var('n0')
        try:
            assert self.e_T == transducer.asymptotic_moments(n0)['expectation']\
                .coefficient(n0)
        except NotImplementedError:
            pass

        self.M_epsilon = [transducer.adjacency_matrix(input=epsilon,
                                                      entry=lambda t: 1)
                          for epsilon in range(self.q)]
        assert self.M == sum(self.M_epsilon)

        self.Delta_epsilon = [transducer.adjacency_matrix(
                input=epsilon,
                entry=lambda t: sum(t.word_out))
                              for epsilon in range(self.q)]
        self.Delta = transducer.adjacency_matrix(
            entry=lambda t: sum(t.word_out))
        assert self.Delta == sum(self.Delta_epsilon)

        self.C_0 = max(self._FC_b_direct_(r).norm(infinity)
                       for r in range(self.q))
        self.C_1 = max(infinity_matrix_norm(d)
                       for d in self.Delta_epsilon)

        self.c=len(self.components)
        self.periods=[c.period for c in self.components]
        self.w=[c.vectors_w() for c in self.components]
        self.coefficient_lambda=[c.coefficient_lambda()
                                for c in self.components]
        self.a=[c.a() for a in self.components]

    @cached_method
    def _FC_b_direct_(self, r):
        r"""

        Compute `\mathbf{b}(r)`, whose ``s`` th component is the
        output of ``self`` when reading the ``q``-ary expansion of
        ``r`` starting in state ``s``.

        INPUT:

        -   ``r`` -- non-negative integer, the input to read

        OUTPUT:

        A vector whose ``s``-th component is the sum of the output of
        ``self`` when reading the ``q``-ary expansion of ``r``
        starting in state ``s`` where ``q`` is the length of the input
        alphabet.

        In contrast to :meth:`_FC_b_recursive_`, the values are computed
        directly via :meth:`.__call__`.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: function('f')
            f
            sage: var('n')
            n
            sage: T = transducers.Recursion([
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: [FSMFourier(T)._FC_b_direct_(r) for r in range(3)]
            [(0, 1, 1), (1, 2, 1), (1, 2, 2)]

        .. SEEALSO::

            :meth:`_FC_b_recursive_`
        """
        from sage.modules.free_module_element import vector
        from sage.rings.integer_ring import ZZ

        expansion = ZZ(r).digits(self.q)
        return vector([sum(self.transducer(expansion, initial_state=s))
                       for s in self.transducer.iter_states()])

    @cached_method
    def _FC_b_recursive_(self, r):
        r"""
        Compute `\mathbf{b}(r)`, whose ``s`` th component is the
        output of ``self`` when reading the ``q``-ary expansion of
        ``r`` starting in state ``s``.

        INPUT:

        -   ``r`` -- non-negative integer, the input to read

        OUTPUT:

        A vector whose ``s``-th component is the sum of the output of
        ``self`` when reading the ``q``-ary expansion of ``r``
        starting in state ``s`` where ``q`` is the length of the input
        alphabet.

        In contrast to :meth:`_FC_b_direct_`, the values are computed
        recursively.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: function('f')
            f
            sage: var('n')
            n
            sage: F = FSMFourier(transducers.Recursion([
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2))
            sage: [F._FC_b_recursive_(r) for r in range(3)]
            [(0, 1, 1), (1, 2, 1), (1, 2, 2)]
            sage: all(F._FC_b_recursive_(r) ==
            ....:     F._FC_b_direct_(r) for r in range(8))
            True

        .. SEEALSO::

            :meth:`_FC_b_direct_`
        """
        from sage.modules.free_module_element import vector

        epsilon = r % self.q
        R = (r - epsilon)/ self.q
        if R == 0:
            return self._FC_b_direct_(r)

        return self.Delta_epsilon[epsilon]*self.ones +\
            self.M_epsilon[epsilon] * self._FC_b_recursive_(R)


    def _H_m_rhs_(self, s, m, remove_poles=False):
        r"""
        Compute `(I - q^{-s}M)\mathbf{H}_m(s)`.

        INPUT:

        -   ``s`` -- a
            :class:`sage.rings.complex_interval.ComplexInterval` element

        -   ``m`` -- non-negative integer

        -   ``remove_poles`` -- (default: ``False``) if ``True`` and ``s
            == 1``, the pole of the zeta functions is removed.

        The function has a pole for ``s == 1``. If ``remove_poles`` is ``True``,
        the zeroth coefficient of the Laurent series at ``s == 1`` is returned
        instead.

        OUTPUT:

        A :class:`sage.rings.complex_interval.ComplexInterval` element.

        EXAMPLES:

        -   One state, always ``1``: this corresponds to the Riemann zeta
            function (minus the first 99 summands, multiplied by 2). ::

                sage: function('f')
                f
                sage: var('n')
                n
                sage: F = FSMFourier(transducers.Recursion([
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     f, n, 2))
                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: F._H_m_rhs_(CIF(2), 100)
                (0.0050250833316668? + 0.?e-18*I)

            Evaluated at `s = 1 + 2k \pi i/\log 2` for `k\neq 0`, the
            result must be zero, as the Riemann zeta function does not
            have a pole there and the first factor therefore annihilates
            the result::

                sage: result = [F._H_m_rhs_(CIF(1+2*k*pi*I/log(2)), 30)[0] for
                ....:      k in range(1, 2)]
                sage: all(value.abs().absolute_diameter() < 1e-12
                ....:      for value in result)
                True
                sage: all(value.contains_zero() for value in result)
                True

        -   One state, one for every input digit, but subtract one at the
            end. This corresponds to the function `L`. Evaluated at `s =
            1 + 2k \pi i/\log 2` for `k\neq 0`, the result must equal
            `\log 2/(2 k \pi i)`, because the `1 \times 1`-matrix is
            singular so that we actually compute the residue at this
            point, multiplied by `\log 2`. The value for `m` does not
            matter, as the difference is analytic and thus does not
            contribute to the residue. ::

                sage: F = FSMFourier(transducers.Recursion([
                ....:         f(2*n + 1) == f(n) + 1,
                ....:         f(2*n) == f(n) + 1,
                ....:         f(0) == -1],
                ....:         f, n, 2))

            We first check that this is indeed the function `L`::

                sage: all(F._FC_b_recursive_(r)[0] ==
                ....:     floor(log(r, base=2)) for r in range(1, 9))
                True

            Next, we check that the result agrees with the known values::

                sage: all(F._H_m_rhs_(CIF(1 + 2*k*pi*I/log(2)),
                ....:         30)[0].overlaps(CIF(log(2)/(2*pi*I*k)))
                ....:     for k in range(1, 2))
                True

        TESTS::

            sage: F._H_m_rhs_(CIF(1), 20, remove_poles=False)
            Traceback (most recent call last):
            ...
            ValueError: remove_poles must be set if and only if s == 1.
            sage: F._H_m_rhs_(CIF(1), 20, remove_poles=False)
            Traceback (most recent call last):
            ...
            ValueError: remove_poles must be set if and only if s == 1.
        """
        verbose("_H_m_rhs_(%s, %s)" % (s, m), level=1)
        if remove_poles != (s == 1):
            raise ValueError("remove_poles must be set if "
                             "and only if s == 1.")

        sigma = s.real()
        CIF = s.parent()
        RIF = sigma.parent()

        def taylor_error(N, sigma, x):
            return min_RIF([
                max_RIF([RIF(1), 1/(1 + x)**(sigma + N)]),
                max_RIF([RIF(N), N/(1 + x)**(sigma + 1)])])


        q = ZZ(self.q)
        log_q = log(RIF(q))
        log_m_1 = log(RIF(m - 1))

        result = sum(self._FC_b_recursive_(r) * r**(-s)
                     for r in reversed(srange(m, q*m)))
        if remove_poles and s == 1:
            result += q**(-s) * sum(
                D * self.ones * (-RIF(ZZ(epsilon)/q + int(epsilon == 0)).psi()
                             - sum((k + ZZ(epsilon)/q + int(epsilon==0))**(-1)
                                   for k in range(0, m-int(epsilon==0))))
                for epsilon, D in enumerate(self.Delta_epsilon))
        else:
            result += q**(-s) * sum(
                D * self.ones * _hurwitz_zeta_(s, ZZ(epsilon)/q, m)
                for epsilon, D in enumerate(self.Delta_epsilon))

        N = 1
        factor = -s * q**(-s - 1)

        while True:
            #assert factor.overlaps(sage.rings.arith.binomial(-s, N) * q**(-s - N))

            if result.is_zero():
                error_acceptable = 0
            else:
                error_acceptable = RIF(2) ** (infinity_vector_norm(result).upper().log2()
                                             - result[0].prec())
            error_bound = factor.abs() * \
                (self.C_0 * (sigma + N -1) * log_q\
                     + self.C_1 * (sigma + N -1) * log_m_1\
                     + self.C_1) / \
                ((sigma + N - 1)**2 \
                     * (m - 1)**(sigma + N - 1)\
                     * log_q) * \
                sum(epsilon**N
                    * taylor_error(N, sigma, ZZ(epsilon)/ZZ(q*m))
                    for epsilon in range(q))

            verbose("    N = %d, error = %s, acceptable_error = %s, "
                    "result = %s" %
                    (N, error_bound, error_acceptable, result), level=2)

            if error_bound.abs() < error_acceptable:
                error_real = RIF(-error_bound, error_bound)
                error = CIF(error_real, error_real)
                error_vector = vector(error for _ in range(len(result)))
                verbose("    N = %d, error = %s, acceptable_error = %s, "
                        "result = %s" %
                        (N, error_bound, error_acceptable, result),
                        level=1)
                return result + error_vector

            result += factor * sum(
                epsilon**N * M
                for epsilon, M in enumerate(self.M_epsilon)) \
                * self._H_m_(s + N, m)

            factor *= (-s - N)/(N + 1) / q
            N += 1

    def _w_H_Res_(self, w, s):
        r"""
        Compute the Residue of `\mathbf{w}^\top\mathbf{H}` at `s`.

        INPUT:

        ``s`` -- a :class:`sage.rings.complex_interval.ComplexInterval` element

        OUTPUT:

        A :class:`sage.rings.complex_interval.ComplexInterval` element.

        EXAMPLES:
        """

        CIF = s.parent()
        log_q = log(CIF(self.q))
        m = s.abs().upper().ceil() + s.parent().precision()

        if s == 1:
            return w * self._H_m_rhs_(s, m, remove_poles=True)/log_q \
                - self.e_T/2
        else:
            return w * self._H_m_rhs_(s, m)/log_q

    # BEGIN_REMOVE_FOR_DOCUMENTATION
    @cached_method(key=lambda self, s, m: (s.real().lower(),
                                           s.real().upper(),
                                           s.imag().lower(),
                                           s.imag().upper(), m))
    # END_REMOVE_FOR_DOCUMENTATION
    def _H_m_(self, s, m):
        r"""
        Compute `\mathbf{H}_m(s)`.

        INPUT:

        -   ``s`` -- a :class:`sage.rings.complex_interval.ComplexInterval` element

        -   ``m`` -- non-negative integer

        OUTPUT:

        A :class:`sage.rings.complex_interval.ComplexInterval` element.

        EXAMPLES:

        -   One state, always ``1``: this corresponds to the Riemann
            zeta function (minus the first 99 summands). ::

                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     f, n, 2)
                sage: FSMFourier(T)._H_m_(CIF(2), 100)
                (0.0100501666633336? + 0.?e-18*I)
        """
        n = self.M.nrows()

        return (matrix.identity(n) - ZZ(self.q)**(-s)*self.M).solve_right(
            self._H_m_rhs_(s, m))

    def _H_(self, s):
        r"""
        Compute `\mathbf{H}_1(s)`.

        INPUT:

        -   ``s`` -- a :class:`sage.rings.complex_interval.ComplexInterval` element

        OUTPUT:

        A :class:`sage.rings.complex_interval.ComplexInterval` element.

        EXAMPLES:

        -   One state, always ``1``: this is the Riemann
            zeta function. ::

                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     f, n, 2)
                sage: result = FSMFourier(T)._H_(CIF(2))
                sage: result
                (1.644934066848226? + 0.?e-18*I)
                sage: result[0].overlaps(CIF(pi^2/6))
                True
        """
        m = s.abs().upper().ceil() + s.parent().precision()
        result = self._H_m_(s, m)
        result += sum(self._FC_b_recursive_(r) * r**(-s)
                      for r in reversed(srange(1, m)))

        return result

    @cached_method
    def FourierCoefficient(self, ell, CIF=ComplexIntervalField()):
        """
        Compute the `\ell`-th Fourier coefficient of the fluctuation
        of the summatory function of the sum of outputs.

        INPUT:

        -   `ell` -- an integer.

        OUTPUT:

        A :class:`ComplexIntervalField` element.

        EXAMPLES:

        -   Binary sum of digits::

                sage: function('f')
                f
                sage: var('n')
                n
                sage: from fsm_fourier import FSMFourier
                sage: T = FSMFourier(transducers.Recursion([
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2))
                sage: def FourierCoefficientDelange(k):
                ....:     if k == 0:
                ....:         return 1/(2*log(2))*(log(2*pi)-1)-3/4
                ....:     return CC(I/(2*k*pi)*(1 + 2*k*pi*I/log(2))^(-1) \
                ....:            *zeta(CC(2*k*pi*I/log(2))))
                sage: all(FourierCoefficientDelange(k)
                ....:     in T.FourierCoefficient(k)
                ....:     for k in range(0, 5)) # long time
                True
        """

        log_q = CIF(log(CIF(self.q)))
        chi_ell =  CIF(2*ell*CIF.pi()*self.I / (self.common_period*log_q))

        w = sum(c.w_ell(ell) for c in self.components)
        if any(w):
            result = 1/(1+chi_ell) * self._w_H_Res_(w, 1+chi_ell)
        else:
            result = CIF(0)

        if ell == 0:
            result += -self.e_T/log_q \
                - self.I*sum(c.vector_w_prime(0)*self.ones
                        for c in self.components)
        return result
