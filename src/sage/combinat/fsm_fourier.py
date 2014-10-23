"""
Fourier Coefficients

::

    sage: import sys
    sage: sys.path.append(".")

"""

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
from sage.rings.real_mpfi import min_RIF, max_RIF


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


class FSMFourier(Transducer):
    """
    Fourier coefficients for the sum of output of transducers.
    """

    @cached_method
    def _fourier_coefficient_data_(self):
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

        - ``components`` -- a list of :class:`FCComponent`, representing the final components.

        EXAMPLES:

        -   Binary sum of digits::

                sage: function('f')
                f
                sage: var('n')
                n
                sage: from fsm_fourier import FSMFourier
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2)
                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: FSMFourier(T)._fourier_coefficient_data_()
                FourierCoefficientData(c=1, periods=[1], period=1, T=[1],
                w=[[(1)]], coefficient_lambda=[1], e_T=1/2, a=[1/2],
                M=[2], M_epsilon=[[1], [1]], Delta=[1],
                Delta_epsilon=[[0], [1]], C_0=1, C_1=1,
                components=[<class 'fsm_fourier.FCComponent'>])

        -   NAF::

                sage: T = transducers.Recursion([
                ....:     f(4*n + 1) == f(n) + 1,
                ....:     f(4*n + 3) == f(n + 1) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2)
                sage: FSMFourier(T)._fourier_coefficient_data_()[:8]
                (
                           [1/3   1   0]
                           [1/3   0   1]
                1, [1], 1, [1/3  -1  -1], [[(1/3, 1/3, 1/3)]], [1],
                1/3, [1/3]
                )

        -   Abelian complexity of the paperfolding sequence::

                sage: T = transducers.Recursion([
                ....:     f(4*n) == f(2*n),
                ....:     f(4*n+2) == f(2*n+1)+1,
                ....:     f(16*n+1) == f(8*n+1),
                ....:     f(16*n+5) == f(4*n+1)+2,
                ....:     f(16*n+11) == f(4*n+3)+2,
                ....:     f(16*n+15) == f(2*n+2)+1,
                ....:     f(1) == 2, f(0) == 0]
                ....:     + [f(16*n+jj) == f(2*n+1)+2 for jj in [3,7,9,13]],
                ....:     f, n, 2)
                sage: FSMFourier(T)._fourier_coefficient_data_()[:8]
                (
                           [1/10    1    0    0    0    0    0    0    0    0]
                           [1/10    0    1    0    0    0    0    0    0    0]
                           [1/10    0    0    1    0    0    0    0    0    0]
                           [1/10    0    0    0    1    0    0    0    0    0]
                           [1/10    0    0    0    0    1    0    0    0    0]
                           [1/10    0    0    0    0    0    1    0    0    0]
                           [1/10    0    0    0    0    0    0    1    0    0]
                           [1/10    0    0    0    0    0    0    0    1    0]
                           [1/10    0    0    0    0    0    0    0    0    1]
                1, [1], 1, [1/10    0    0   -3   -2   -2   -2   -1   -1   -1],
                [[(0, 0, 3/13, 2/13, 2/13, 2/13, 1/13, 1/13, 1/13, 1/13)]],
                [1], 8/13, [8/13]
                )

        -   Artificial example, one-periodic, 2 states::

                sage: T = transducers.Recursion([
                ....:     f(4*n) == f(2*n)+0,
                ....:     f(4*n+2) == f(n)+1,
                ....:     f(2*n+1) == f(n),
                ....:     f(0) == 0],
                ....:     f, n, 2)
                sage: FSMFourier(T)._fourier_coefficient_data_()[:8]
                (
                           [1/2   1]
                1, [1], 1, [1/2  -1],
                [[(1/2, 1/2)]],
                [1], 1/4, [1/4]
                )

        -   Artificial example, period 3::

                sage: T = transducers.Recursion([
                ....:     f(8*n) == f(4*n+3)+3,
                ....:     f(8*n+4) == f(4*n+3)+1,
                ....:     f(8*n+2) == f(4*n+3)+2,
                ....:     f(8*n+6) == f(4*n+3)-1,
                ....:     f(8*n+1) == f(4*n)+5,
                ....:     f(8*n+5) == f(4*n+2)+1,
                ....:     f(8*n+3) == f(4*n+1)+2,
                ....:     f(8*n+7) == f(4*n+1),
                ....:     f(0) == 0],
                ....:     f, n, 2)
                sage: FSMFourier(T)._fourier_coefficient_data_()[:8]
                (
                              [         1/7            1            1            1            0            0            0]
                              [         1/7      4*zeta3 -4*zeta3 - 4            0            1            0            0]
                              [         1/7     -2*zeta3  2*zeta3 + 2            0            0            1            0]
                              [         1/7 -4*zeta3 - 4      4*zeta3            0            0            0            1]
                              [         1/7 -4*zeta3 - 4      4*zeta3            0            0            0           -1]
                              [         1/7      4*zeta3 -4*zeta3 - 4            0            0            0            0]
                   1, [3], 3, [         1/7            4            4            0            0            0            0],
                   [[(0, 0, 0, 1/6, 1/6, 1/3, 1/3),
                     (0, 0, 0, 1/24*zeta3, 1/24*zeta3, -1/12*zeta3 - 1/12, 1/12),
                     (0, 0, 0, -1/24*zeta3 - 1/24, -1/24*zeta3 - 1/24, 1/12*zeta3, 1/12)]],
                   [1], 7/4, [7/4]
                )

        -   Artificial example, period 2, vanishing w-vector::

                sage: T = transducers.Recursion([
                ....:     f(4*n) == f(2*n+1)+1,
                ....:     f(4*n+1) == f(2*n)+2,
                ....:     f(4*n+2) == f(2*n+1)+3,
                ....:     f(4*n+3) == f(2*n)-1,
                ....:     f(0) == 0],
                ....:     f, n, 2)
                sage: FSMFourier(T)._fourier_coefficient_data_()[:8]
                (
                           [1/3   0   1]
                           [1/3   1   0]
                1, [2], 2, [1/3  -1   0],
                [[(0, 1/2, 1/2), (0, 0, 0)]],
                [1], 5/4, [5/4]
                )

        -   Artificial example with two final components of periods `2`
            and `3`, respectively::

                sage: T = FSMFourier([(0, 1, 0, 1), (1, 2, 0, 1),
                ....:     (2, 1, 0, 2), (1, 2, 1, 0), (2, 1, 1, 2),
                ....:     (0, -1, 1, 1), (-1, -2, 1, 1), (-2, -3, 1, 1),
                ....:     (-3, -1, 1, 1), (-1, -2, 0, 2), (-2, -3, 0, 1),
                ....:     (-3, -1, 0, 3)],
                ....:     initial_states=[0],
                ....:     final_states=[0, 1, 2, -3, -2, -1])
                sage: T._fourier_coefficient_data_()[:8]
                (
                              [        1/7           1           1         1/5           1           1]
                              [          0           0           0         2/5          -2           0]
                              [          0           0           0         2/5           2           0]
                              [        2/7           2           2           0           0           0]
                              [        2/7    -2*zeta6 2*zeta6 - 2           0           0           0]
                2, [3, 2], 6, [        2/7 2*zeta6 - 2    -2*zeta6           0           0           0],
                [[(0, 0, 0, 1/6, 1/6, 1/6),
                  (0, 0, 0, 1/6, 1/6*zeta6 - 1/6, -1/6*zeta6),
                  (0, 0, 0, 1/6, -1/6*zeta6, 1/6*zeta6 - 1/6)],
                 [(0, 1/4, 1/4, 0, 0, 0), (0, -1/4, 1/4, 0, 0, 0)]],
                [1/2, 1/2], 11/8, [5/4, 5/4]
                )
        """


        import collections
        import itertools
        import operator

        from sage.calculus.var import var
        from sage.modules.free_module import VectorSpace
        from sage.rings.arith import lcm
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        from sage.structure.sage_object import SageObject
        from sage.symbolic.constants import I

        FourierCoefficientData = collections.namedtuple(
            "FourierCoefficientData",
            ["c", "periods", "period", "T", "w", "coefficient_lambda",
             "e_T", "a", "M", "M_epsilon", "Delta", "Delta_epsilon",
             "C_0", "C_1", "components"])

        positions = dict((state.label(), j)
                         for j, state in enumerate(self.iter_states()))
        q = len(self.input_alphabet)
        Y = PolynomialRing(QQ, 'Y').gen()

        class FCComponent(SageObject):
            """Hold a final component and associated data."""
            def __init__(self, fsm, parent):
                self.fsm = fsm
                self.period = fsm.graph().period()
                self.n_states = len(self.fsm.states())
                self.parent = parent

            def mask(self, n, components):
                nrows = sum(c.n_states for c in components if c != self)
                mask = matrix(
                    nrows, n,
                    [standard_basis[positions[state.label()]]
                     for other in components
                     if other != self
                     for state in other.fsm.iter_states()])
                return mask

            def eigenvectors(self, M, components):
                mask = self.mask(M.nrows(), components)
                def eigenvector(j):
                    eigenvalue = q * alpha**(
                        j * common_period / self.period)
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
                return self.eigenvectors(M, components)

            @cached_method()
            def left_eigenvectors(self):
                left_eigenvectors = self.eigenvectors(M.transpose(),
                                                      components)
                return [w/(v*w) for v, w
                        in itertools.izip(self.right_eigenvectors(),
                                          left_eigenvectors)]

            @cached_method()
            def vectors_w(self):
                return [(initial_vector*v)*w for v, w
                        in itertools.izip(self.right_eigenvectors(),
                                          self.left_eigenvectors())]

            @cached_method()
            def coefficient_lambda(self):
                ones = vector(1 for _ in range(M.nrows()))
                products = [w*ones for w in self.vectors_w()]
                assert all(e.is_zero() for e in products[1:])
                return products[0]

            @cached_method()
            def mu_prime(self):
                p = self.fsm.adjacency_matrix(
                    entry=lambda t:Y**sum(t.word_out)).charpoly('Z')
                Z = p.parent().gen()
                assert p(Y=1, Z=q) == 0
                mu_prime_Z = (- p.derivative(Y)/p.derivative(Z))(
                    Y=1, Z=q)
                I = CyclotomicField(4).gen()
                return I*mu_prime_Z

            @cached_method()
            def a(self):
                return QQ(-I * self.mu_prime()/q)

            def w_ell(self, ell):
                if common_period.divides(ell*self.period):
                    k = self.period*ell/common_period % self.period
                    return vector(field_to_CIF(c) for c in self.vectors_w()[k])
                else:
                    return vector(0 for _ in self.vectors_w()[0])

            def M_prime(self):
                I = CyclotomicField(4).gen()
                return self.parent.adjacency_matrix(entry=lambda t: sum(t.word_out))*I

            def vector_v_prime(self, k):
                mask = self.mask(M.nrows(), components)
                ones = matrix([1 for _ in self.vectors_w()[0]])
                eigenvalue = q * alpha**(
                        k * common_period / self.period)
                S = matrix.block(CyclotomicField(4*common_period),
                     [[M - eigenvalue*matrix.identity(M.nrows())],
                      [ones],
                      [mask]],
                     subdivide=False)
                eigenvector_right = vector(field, self.right_eigenvectors()[k])
                right_side = - matrix.block(CyclotomicField(4*common_period),
                                            [[self.M_prime() - self.mu_prime()*matrix.identity(M.nrows())],
                                             [0*ones],
                                             [0*mask]],
                                            subdivide=False) * eigenvector_right
                v_prime = S.solve_right(right_side)
                return v_prime

            def vector_w_prime(self, k):
                mask = self.mask(M.nrows(), components)
                eigenvalue = q * alpha**(
                        k * common_period / self.period)
                eigenvector_right = vector(field, self.right_eigenvectors()[k])
                eigenvector_left = vector(field, self.left_eigenvectors()[k])
                S = matrix.block(CyclotomicField(4*common_period),
                     [[M.transpose() - eigenvalue*matrix.identity(M.nrows())],
                      [matrix(eigenvector_right)],
                      [mask]],
                     subdivide=False)
                right_side = - matrix.block(CyclotomicField(4*common_period),
                                            [[self.M_prime().transpose() - self.mu_prime()*matrix.identity(M.nrows())],
                                             [matrix(self.vector_v_prime(k))],
                                             [0*mask]],
                                            subdivide=False) * eigenvector_left
                left_prime = S.solve_right(right_side)
                w_prime = initial_vector*self.vector_v_prime(k)*eigenvector_left \
                    + initial_vector*eigenvector_right*left_prime
                return w_prime

            def lambda_ell_prime(self, ell):
                ones = vector(1 for _ in self.vectors_w()[0])
                if common_period.divides(ell*self.period):
                    k = self.period*ell/common_period % self.period
                    lambda_prime = self.vector_w_prime(k)*ones
                    return lambda_prime
                else:
                    return 0
                    


        components = [FCComponent(c, self) for c in self.final_components()]
        common_period = lcm([c.period for c in components])
        field = CyclotomicField(common_period)
        alpha = field.gen()
        field_to_CIF = field.hom([ComplexIntervalField().zeta(common_period)], check=False)
        M = self.adjacency_matrix(entry=lambda t: 1)
        standard_basis = VectorSpace(field, M.nrows()).basis()

        if len(self.initial_states()) != 1:
            raise NotImplementedError(
                "Transducer does not have a unique initial state.")
        initial_vector = standard_basis[positions[
                self.initial_states()[0].label()]]

        right_eigenvectors = list(itertools.chain(
                *(c.right_eigenvectors()
                  for c in components)))

        left_eigenvectors = list(itertools.chain(
                *(c.left_eigenvectors()
                  for c in components)))

        annihilated_by_left = matrix(left_eigenvectors).\
            right_kernel_matrix().transpose()

        T = matrix.block([[matrix.column(right_eigenvectors),
                           annihilated_by_left]],
                         subdivide=False)

        assert T.is_square()
        assert T.nrows() == M.nrows()
        assert T.is_invertible()

        check = T.inverse() * M * T
        eigenvalues = [q * alpha**(j * common_period/c.period)
                       for c in components
                       for j in range(c.period)]
        check_dont_care = check.submatrix(len(eigenvalues),
                                          len(eigenvalues))
        assert (matrix.block(
                [[matrix.diagonal(eigenvalues), ZZ(0)],
                 [ZZ(0), check_dont_care]],
                     subdivide=False) - check).is_zero()

        assert (T.inverse().submatrix(nrows=len(left_eigenvectors))
                - matrix(left_eigenvectors)).is_zero()

        e_T = sum(c.a()*c.coefficient_lambda()
                  for c in components)

        var('n0')
        try:
            assert e_T == self.asymptotic_moments(n0)['expectation']\
                .coefficient(n0)
        except NotImplementedError:
            pass

        M_epsilon = [self.adjacency_matrix(input=epsilon,
                                           entry=lambda t: 1)
                     for epsilon in range(q)]
        assert M == sum(M_epsilon)

        Delta_epsilon = [self.adjacency_matrix(
                             input=epsilon,
                             entry=lambda t: sum(t.word_out))
                         for epsilon in range(q)]
        Delta = self.adjacency_matrix(
                             entry=lambda t: sum(t.word_out))
        assert Delta == sum(Delta_epsilon)

        C_0 = max(self._FC_b_direct_(r).norm(infinity)
                  for r in range(q))
        C_1 = max(infinity_matrix_norm(d)
                  for d in Delta_epsilon)

        return FourierCoefficientData(
            c=len(components),
            periods=[c.period for c in components],
            period=common_period,
            T=T,
            w=[c.vectors_w() for c in components],
            coefficient_lambda=[c.coefficient_lambda()
                                for c in components],
            e_T=e_T,
            a=[c.a() for a in components],
            M=M,
            M_epsilon=M_epsilon,
            Delta=Delta,
            Delta_epsilon=Delta_epsilon,
            C_0=C_0,
            C_1=C_1,
            components=components)

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

        q = len(self.input_alphabet)
        expansion = ZZ(r).digits(q)
        return vector([sum(self(expansion, initial_state=s))
                       for s in self.iter_states()])

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
            sage: T = transducers.Recursion([
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     f, n, 2)
            sage: [FSMFourier(T)._FC_b_recursive_(r) for r in range(3)]
            [(0, 1, 1), (1, 2, 1), (1, 2, 2)]
            sage: all(FSMFourier(T)._FC_b_recursive_(r) ==
            ....:     FSMFourier(T)._FC_b_direct_(r) for r in range(8))
            True

        .. SEEALSO::

            :meth:`_FC_b_direct_`
        """
        from sage.modules.free_module_element import vector

        q = len(self.input_alphabet)
        epsilon = r % q
        R = (r - epsilon)/ q
        if R == 0:
            return self._FC_b_direct_(r)

        d = self._fourier_coefficient_data_()
        ones = vector(1 for _ in range(d.M.nrows()))
        return d.Delta_epsilon[epsilon]*ones +\
            d.M_epsilon[epsilon] * self._FC_b_recursive_(R)


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
                sage: T = FSMFourier(transducers.Recursion([
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     f, n, 2))
                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: T._H_m_rhs_(CIF(2), 100)
                (0.0050250833316668? + 0.?e-18*I)

            Evaluated at `s = 1 + 2k \pi i/\log 2` for `k\neq 0`, the
            result must be zero, as the Riemann zeta function does not
            have a pole there and the first factor therefore annihilates
            the result::

                sage: result = [T._H_m_rhs_(CIF(1+2*k*pi*I/log(2)), 30)[0] for
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

                sage: T = FSMFourier(transducers.Recursion([
                ....:         f(2*n + 1) == f(n) + 1,
                ....:         f(2*n) == f(n) + 1,
                ....:         f(0) == -1],
                ....:         f, n, 2))

            We first check that this is indeed the function `L`::

                sage: all(T._FC_b_recursive_(r)[0] ==
                ....:     floor(log(r, base=2)) for r in range(1, 9))
                True

            Next, we check that the result agrees with the known values::

                sage: all(T._H_m_rhs_(CIF(1 + 2*k*pi*I/log(2)),
                ....:         30)[0].overlaps(CIF(log(2)/(2*pi*I*k)))
                ....:     for k in range(1, 2))
                True
        """
        verbose("_H_m_rhs_(%s, %s)" % (s, m), level=1)

        sigma = s.real()
        CIF = s.parent()
        RIF = sigma.parent()

        def taylor_error(N, sigma, x):
            return min_RIF([
                max_RIF([RIF(1), 1/(1 + x)**(sigma + N)]),
                max_RIF([RIF(N), N/(1 + x)**(sigma + 1)])])


        q = ZZ(len(self.input_alphabet))
        log_q = log(RIF(q))
        log_m_1 = log(RIF(m - 1))

        Delta_epsilon = self._fourier_coefficient_data_().Delta_epsilon
        M_epsilon = self._fourier_coefficient_data_().M_epsilon
        C_0 = self._fourier_coefficient_data_().C_0
        C_1 = self._fourier_coefficient_data_().C_1

        ones = vector(1 for _ in range(M_epsilon[0].ncols()))

        result = sum(self._FC_b_recursive_(r) * r**(-s)
                     for r in reversed(srange(m, q*m)))
        if remove_poles and s == 1:
            result += q**(-s) * sum(
                D * ones * (-RIF(ZZ(epsilon)/q + int(epsilon == 0)).psi()
                             - sum((k + ZZ(epsilon)/q + int(epsilon==0))**(-1)
                                   for k in range(0, m-int(epsilon==0))))
                for epsilon, D in enumerate(Delta_epsilon))
        else:
            result += q**(-s) * sum(
                D * ones * _hurwitz_zeta_(s, ZZ(epsilon)/q, m)
                for epsilon, D in enumerate(Delta_epsilon))

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
                (C_0 * (sigma + N -1) * log_q\
                     + C_1 * (sigma + N -1) * log_m_1\
                     + C_1) / \
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
                for epsilon, M in enumerate(M_epsilon)) \
                * self._H_m_(s + N, m)

            factor *= (-s - N)/(N + 1) / q
            N += 1

    def _H_Res_(self, s):
        r"""
        Compute the Residue of `\mathbf{H}` at `s`.

        INPUT:

        ``s`` -- a :class:`sage.rings.complex_interval.ComplexInterval` element

        OUTPUT:

        A :class:`sage.rings.complex_interval.ComplexInterval` element.

        EXAMPLES:
        """

        CIF = s.parent()
        q = len(self.input_alphabet)
        log_q = CIF(log(q))
        m = s.abs().upper().ceil() + s.parent().precision()

        if s == 1:
            return self._H_m_rhs_(s, m, remove_poles=True)/log_q
        else:
            return self._H_m_rhs_(s, m)/log_q

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
        q = ZZ(len(self.input_alphabet))
        M = self._fourier_coefficient_data_().M
        n = M.nrows()

        return (matrix.identity(n) - q**(-s)*M).solve_right(
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
                ....:     for k in range(0, 5))
                True
        """

        q = len(self.input_alphabet)
        log_q = CIF(log(q))
        data = self._fourier_coefficient_data_()
        I = CIF.gens()[0]
        chi_ell =  CIF(2*ell*CIF.pi()*I / (data.period*log_q))

        w = sum(c.w_ell(ell) for c in data.components)
        if any(w):
            result = 1/(1+chi_ell) * w * self._H_Res_(1+chi_ell)
        else:
            result = CIF(0)

        if ell == 0:
            result += -data.e_T/log_q - data.e_T/2 - I*sum(c.lambda_ell_prime(ell) for c in data.components)

        return result
