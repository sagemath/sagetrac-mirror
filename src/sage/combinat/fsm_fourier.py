"""
Fourier Coefficients

::

    sage: import sys
    sage: sys.path.append(".")

"""

from sage.combinat.finite_state_machine import Transducer
from sage.misc.cachefunc import cached_method
from sage.rings.integer_ring import ZZ

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
    from sage.misc.misc import srange, verbose
    from sage.rings.arith import bernoulli, falling_factorial

    CIF = s.parent()
    RIF = s.real().parent()

    if ZZ(1) in s:
        raise ZeroDivisionError("zeta is singular at 1.")

    # We rely on (2pi)^-N for convergence of the error term.
    # As a conservative estimate, 2pi is approximately 2^2,
    # so we will need N ~ s.prec()/2 to achieve an error which
    # is less than the resolution of s.
    # In order to have the falling factorial (-s)^\underline{N}
    # smaller than (M+a)^N, we choose M>|s|+N
    M = max(m, (s.prec()/ZZ(2)).ceil() + ZZ(s.abs().upper().ceil()))
    verbose("_hurwitz_zeta_(%s, %s, %s): M = %d" % (s, alpha, m, M),
            level=1)

    sigma = s.real()
    result = sum((r + alpha)**(-s) for r in reversed(srange(m, M)))
    result += (M + alpha)**(1-s)/(s-1)
    factor = (M + alpha)**(-s)
    result += factor/2

    N = 0
    error_factor = RIF(4)
    N_factorial = ZZ(1)

    while True:
        N += 2
        factor *= (-s - N + 2)/(M + alpha)
        assert factor.overlaps(falling_factorial(-s, N - 1)/(M + alpha)**(s + N - 1))
        N_factorial *= (N - 1)*N
        assert ZZ(N).factorial() == N_factorial
        result -= bernoulli(N)/N_factorial * factor

        factor *= (-s - N + 1)
        assert factor.overlaps(falling_factorial(-s, N)/(M + alpha)**(s + N - 1))

        error_factor /= 4*RIF.pi()**2
        assert error_factor.overlaps(4/(2*RIF.pi())**N)
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
        assert factor.overlaps(falling_factorial(-s, N)/(M + alpha)**(s + N))


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

        - ``Delta`` -- output matrix ``\Delta``.

        - ``Delta_epsilon`` -- list of partial output matrices
          `\Delta_\varepsilon`.

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
                sage: FSMFourier(T)._fourier_coefficient_data_()
                FourierCoefficientData(c=1, periods=[1], period=1, T=[1],
                w=[[(1)]], coefficient_lambda=[1], e_T=1/2, a=[1/2],
                M=[2], M_epsilon=[[1], [1]], Delta=[1],
                Delta_epsilon=[[0], [1]])

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
        from sage.matrix.constructor import matrix
        from sage.modules.free_module import VectorSpace
        from sage.modules.free_module_element import vector
        from sage.rings.arith import lcm
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
        from sage.rings.rational_field import QQ
        from sage.structure.sage_object import SageObject
        from sage.symbolic.constants import I

        FourierCoefficientData = collections.namedtuple(
            "FourierCoefficientData",
            ["c", "periods", "period", "T", "w", "coefficient_lambda",
             "e_T", "a", "M", "M_epsilon", "Delta", "Delta_epsilon"])

        positions = dict((state.label(), j)
                         for j, state in enumerate(self.iter_states()))
        q = len(self.input_alphabet)
        Y = PolynomialRing(QQ, 'Y').gen()

        class FCComponent(SageObject):
            """Hold a final component and associated data."""
            def __init__(self, fsm):
                self.fsm = fsm
                self.period = fsm.graph().period()
                self.n_states = len(self.fsm.states())

            def eigenvectors(self, M, components, common_period):
                nrows = sum(c.n_states for c in components if c != self)
                mask = matrix(
                    nrows, M.ncols(),
                    [standard_basis[positions[state.label()]]
                     for other in components
                     if other != self
                     for state in other.fsm.iter_states()])

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
                return self.eigenvectors(M, components, period)

            @cached_method()
            def left_eigenvectors(self):
                left_eigenvectors = self.eigenvectors(M.transpose(),
                                                      components,
                                                      period)
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
                return I*mu_prime_Z

            @cached_method()
            def a(self):
                return QQ(-I * self.mu_prime()/q)


        components = [FCComponent(c) for c in self.final_components()]
        period = lcm([c.period for c in components])
        field = CyclotomicField(period)
        alpha = field.gen()
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
        eigenvalues = [q * alpha**(j * period/c.period)
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

        return FourierCoefficientData(
            c=len(components),
            periods=[c.period for c in components],
            period=period,
            T=T,
            w=[c.vectors_w() for c in components],
            coefficient_lambda=[c.coefficient_lambda()
                                for c in components],
            e_T=e_T,
            a=[c.a() for a in components],
            M=M,
            M_epsilon=M_epsilon,
            Delta=Delta,
            Delta_epsilon=Delta_epsilon)

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
