r"""
Fourier Coefficients

Introduction
============

Let `T` be a complete deterministic transducer reading the `q`-ary
digit expansion of a random integer `n` with `0\le n<N`.

Then the expected sum of the output labels of the transducer is
`e_T \log_q N + \Phi(\log_q N) + o(1)` for some periodic fluctuation
`\Phi` and some constant `e_T`.

This module computes the Fourier Coefficients of `\Phi`, following
[HKP2015]_.

Contents
========

:class:`FSMFourier`
-------------------
The class :class:`FSMFourier` is the main class of this module, it
precomputes values needed for all Fourier coefficients and then offers
the method :meth:`~FSMFourier.FourierCoefficient`.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FSMFourier.FourierCoefficient` | Compute a Fourier Coefficient
    :meth:`~FSMFourier.plot_fluctuation` | Generate a plot of the fluctuation
    :meth:`~FSMFourier.plot_fluctuation_asymptote` | Generate an asymptote plot of the fluctuation
    :meth:`~FSMFourier.b0` | Compute the final output word

:class:`FSMFourierComponent`
------------------------------
Data corresponding to a final component of the transducer is stored in
:class:`FSMFourierComponent`. This mainly comprises eigenvectors.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FSMFourierComponent.mu_prime` | Derivative `\mu_{j0}'(0)`
    :meth:`~FSMFourierComponent.a` | Constant `a_j`
    :meth:`~FSMFourierComponent.coefficient_lambda` | Constant `\lambda_j`
    :meth:`~FSMFourierComponent.right_eigenvectors` | Right eigenvectors `\mathbf{v}_{jk}(0)`
    :meth:`~FSMFourierComponent.v_prime` | Derivative of the right eigenvector `\mathbf{v}'_{jk}(0)`
    :meth:`~FSMFourierComponent.left_eigenvectors` | Left eigenvectors
    :meth:`~FSMFourierComponent.w` | Scaled left eigenvectors `\mathbf{w}_{jk}(0)`
    :meth:`~FSMFourierComponent.w_prime` | Derivative of the left eigenvector `\mathbf{w}'_{jk}(0)`
    :meth:`~FSMFourierComponent.w_ell` | Left eigenvector to given eigenvalue.

:class:`FSMFourierCache`
------------------------
This is a mostly internal class speeding up some of the computation.

.. csv-table::
    :class: contentstable
    :widths: 30, 70
    :delim: |

    :meth:`~FSMFourierCache.b` | Compute `\mathbf{b}(r)`.
    :meth:`~FSMFourierCache.fluctuation_empirical` | Compute the fluctuation empirically.
    :meth:`~FSMFourierCache.fluctuation_fourier` | Compute the fluctuation via an the Fourier series.

Example
=======

::

    sage: function('f')
    f
    sage: var('n')
    n
    sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
    sage: T = transducers.Recursion([
    ....:     f(2*n + 1) == f(n) + 1,
    ....:     f(2*n) == f(n),
    ....:     f(0) == 0],
    ....:     2, f, n)
    sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
    sage: F = FSMFourier(T) # optional - arb
    sage: F.FourierCoefficient(0) # optional - arb
    -0.1455994557084? + 0.?e-16*I
    sage: F.FourierCoefficient(10) # optional - arb
    0.0016818573864? + 0.0003201978624?*I

AUTHORS:

- Clemens Heuberger (2014-08-14--2014-10-13): initial version
- Sara Kropf (2014-10-13--2014-10-21): corrections and 0th Fourier Coefficient
- Clemens Heuberger (2014-10-23--2014-10-26): Optimization and Preparation for Sage

ACKNOWLEDGEMENT:

- Clemens Heuberger and Sara Kropf are supported by the
  Austrian Science Fund (FWF): P 24644-N26.

.. TODO::

    The following private and special methods are included now for
    proofreading, they should be excluded. Furthermore, this module cannot
    be shown in the documentation as it depends on arb.

Classes and Methods
===================

.. autofunction:: _hurwitz_zeta_
.. autofunction:: _reduce_resolution_
.. automethod:: FSMFourierComponent.__init__
.. automethod:: FSMFourierComponent._eigenvectors_
.. automethod:: FSMFourierComponent._mask_
.. automethod:: FSMFourierCache.__init__
.. automethod:: FSMFourier.__init__
.. automethod:: FSMFourier._H_m_rhs_
.. automethod:: FSMFourier._w_H_Res_
.. automethod:: FSMFourier._H_m_
"""
#*****************************************************************************
# Copyright (C) 2014 Clemens Heuberger <clemens.heuberger@aau.at>
#               2014 Sara Kropf <sara.kropf@aau.at>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                http://www.gnu.org/licenses/
#*****************************************************************************
from libc.stdlib cimport malloc, realloc, free
import itertools

from sage.combinat.finite_state_machine import Transducer
from sage.functions.log import log
from sage.functions.other import floor
from sage.libs.arb.acb cimport *
from sage.libs.arb.acb_mat cimport *
from sage.matrix.constructor import matrix
from sage.matrix.matrix_complex_ball_dense cimport (
    matrix_to_acb_mat, acb_mat_to_matrix)
from sage.matrix.matrix_space import MatrixSpace
from sage.misc.cachefunc import cached_method
from sage.misc.misc import srange, verbose
from sage.modules.free_module_element import vector
from sage.modules.free_module_element cimport FreeModuleElement_generic_dense
import sage.rings.arith
from sage.rings.complex_interval cimport ComplexIntervalFieldElement
from sage.rings.complex_ball_acb cimport ComplexIntervalFieldElement_to_acb
from sage.rings.complex_interval_field import ComplexIntervalField
from sage.functions.other import ceil
from sage.rings.infinity import infinity
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.real_arb cimport arb_to_mpfi
from sage.rings.real_mpfi import RealIntervalField
from sage.rings.real_mpfi cimport RealIntervalFieldElement
from sage.rings.real_mpfi import RealIntervalFieldElement as RIF_E
from sage.structure.sage_object cimport SageObject



def infinity_vector_norm(v):
    """
    Compute the infinity norm of a vector of
    :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

    INPUT:

    - ``v`` -- a vector of
      :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`.

    OUTPUT:

    A :class:`~sage.rings.real_mpfi.RealIntervalFieldElement` element.

    EXAMPLES::

        sage: from sage.combinat.fsm_fourier import infinity_vector_norm # optional - arb
        sage: a = RIF(1, 4)
        sage: b = RIF(-3, -2)
        sage: infinity_vector_norm(vector([a, b])).endpoints() # optional - arb
        (2.00000000000000, 4.00000000000000)
        sage: infinity_vector_norm(vector([b, a])).endpoints() # optional - arb
        (2.00000000000000, 4.00000000000000)

    Note that the
    :meth:`sage.modules.free_module_element.FreeModuleElement.norm`
    method is inadequate, as it uses the Python max instead of
    :meth:`~sage.rings.real_mpfi.RealIntervalFieldElement.max`::

        sage: vector([a, b]).norm(infinity).endpoints()
        (1.00000000000000, 4.00000000000000)
        sage: vector([b, a]).norm(infinity).endpoints()
        (2.00000000000000, 3.00000000000000)
    """
    return RIF_E.max(*(abs(r) for r in v))


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

        sage: from sage.combinat.fsm_fourier import infinity_matrix_norm # optional - arb
        sage: M = matrix([[1, 2], [-5, -1]])
        sage: infinity_matrix_norm(M) # optional - arb
        6
        sage: M.norm(infinity)
        6.0
        sage: M = matrix([[0]], sparse=True)
        sage: infinity_matrix_norm(M) # optional - arb
        0
        sage: M.norm(infinity)
        Traceback (most recent call last):
        ...
        TypeError: base_ring (=Category of objects) must be a ring
    """
    if A.is_zero():
        return 0
    return max(sum(r) for r in A.apply_map(abs).rows())


def _hurwitz_zeta_(s, alpha,  m=0, max_approximation_error=0):
    r"""
    Compute the truncated Hurwitz zeta function `\sum_{k\ge m} (k+\alpha)^{-s}`.

    INPUT:

    -   ``s`` -- a :class:`ComplexIntervalFieldElement`

    -   ``alpha`` -- a :class:`~sage.rings.real_mpfi.RealIntervalFieldElement`

    -   ``m`` -- a positive integer

    -   ``max_approximation_error`` -- a non-negative number; an
        approximation error less than ``max_approximation_error`` is
        accepted even if it is not small with respect to the result.

    OUTPUT:

    A :class:`ComplexIntervalFieldElement` in the same field as ``s``.

    EXAMPLES:

    -   Simple example::

            sage: from sage.combinat.fsm_fourier import _hurwitz_zeta_ # optional - arb
            sage: _hurwitz_zeta_(CIF(2), RIF(3/4), 10) # optional - arb
            0.097483848201852? + 0.?e-17*I

    -   Compare with well-known value `\zeta(2)=\zeta(2, 1)=\pi^2/6`::

            sage: _hurwitz_zeta_(CIF(2), 1) # optional - arb
            1.64493406684823? + 0.?e-17*I
            sage: (_hurwitz_zeta_(CIF(2), 1) - CIF(pi)^2/6).abs()<10^(-13) # optional - arb
            True
            sage: _hurwitz_zeta_(CIF(2*pi*I/log(2)), 1) # optional - arb
            1.598734526809? + 0.278338669639?*I

    -   There is a singularity at `s=1`. ::

            sage: _hurwitz_zeta_(CIF(RIF(0.9, 1.1), (-0.1, 0.1)), 1) # optional - arb
            Traceback (most recent call last):
            ...
            ZeroDivisionError: zeta is singular at 1.

    -   Checking the function for non-positive integers::

            sage: zeta(-1) in _hurwitz_zeta_(CIF(-1), 1) # optional - arb
            True
            sage: zeta(0) in _hurwitz_zeta_(CIF(0), 1) # optional - arb
            True

    -   Debugging output can be enabled using
        :func:`~sage.misc.misc.set_verbose`. ::

            sage: set_verbose(2)
            sage: _hurwitz_zeta_(CIF(1+100/log(2)*I), 1) # optional - arb
            verbose 1 (...) _hurwitz_zeta_(1 + 144.2695040888963?*I, 1, 0): M = 172
            verbose 2 (...)     N = 2, error = 0.0352354068797?, acceptable_error = 3.41199291042927e-13, result = 2.125571548789? + 0.511221280470?*I
            verbose 2 (...)     N = 4, error = 0.000310532577681?, acceptable_error = 3.41254802194158e-13, result = 2.125575595864? + 0.51121897538?*I
            verbose 2 (...)     N = 6, error = 3.65215306101?e-6, acceptable_error = 3.41351946708813e-13, result = 2.125575660430? + 0.511218933060?*I
            verbose 2 (...)     N = 8, error = 4.83820940904?e-8, acceptable_error = 3.41407457860044e-13, result = 2.125575661484? + 0.511218932225?*I
            verbose 2 (...)     N = 10, error = 6.84787778802?e-10, acceptable_error = 3.41504602374699e-13, result = 2.125575661501? + 0.511218932208?*I
            verbose 2 (...)     N = 12, error = 1.011644165731?e-11, acceptable_error = 3.41560113525930e-13, result = 2.125575661501? + 0.51121893221?*I
            verbose 2 (...)     N = 14, error = 1.54088223831?e-13, acceptable_error = 3.41615624677161e-13, result = 2.125575661501? + 0.51121893221?*I
            verbose 1 (...)     N = 14, error = 1.54088223831?e-13, acceptable_error = 3.41615624677161e-13, result = 2.125575661501? + 0.51121893221?*I
            2.125575661501? + 0.51121893221?*I
            sage: set_verbose(0)

    -   The current implementation does not work well with values with negative real
        part, all precision is lost::

            sage: _hurwitz_zeta_(CIF(-15+I), 1) # optional - arb
            0.?e12 + 0.?e12*I
            sage: hurwitz_zeta(ComplexField(200)(-15 + I), 1) # optional - arb
            0.66621329305522618549073659441004805750538410627754288912677
            - 0.84614995218731390314834131849322502368334938802936485299779*I

        A work-around is to start with higher precision; however, we have
        to clear the cache first::

            sage: _hurwitz_zeta_(ComplexIntervalField(200)(-15 + I), 1) # optional - arb
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
            rounding_error = result.real().upper() + result.imag().upper() \
                - result.real().lower() - result.imag().lower()
            error_acceptable = max(max_approximation_error, rounding_error/8)

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

def _reduce_resolution_(data, x_min, x_max, resolution):
    r"""
    Reduce a list of pairs `(x, y)` to a list of triples `(x,
    y_\mathrm{min}, y_\mathrm{max})` corresponding to
    ``resolution`` equidistant `x` values.

    INPUT:

    - ``data`` -- a list (or iterable) of pairs of doubles.

    - ``x_min`` -- a double, start of the interval.

    - ``x_max`` -- a double, end of the interval.

    - ``resolution`` -- a positive integer, the number of points
      in `x` direction.

    OUTPUT:

    A list of triples of doubles.

    Each `(x, y)` is mapped to some
    `(x', y_\mathrm{min}, y_\mathrm{max})`
    such that `y_\mathrm{min}\le y\le y_\mathrm{max}` and
    `0\le x-x'< (x_\mathrm{max}-x_\mathrm{min})/\mathit{resolution}`.

    A list plot of the original list thus corresponds to plot of the
    vertical line segments defined by the output.

    EXAMPLES::

        sage: from sage.combinat.fsm_fourier import _reduce_resolution_ # optional - arb
        sage: _reduce_resolution_(((i/10, i) for i in range(10)), # optional - arb
        ....:                    0, 1, 2)
        [(0.0, 0, 4), (0.5, 5, 9)]
        sage: _reduce_resolution_(((1 - i/10, i) for i in range(1, 11)), # optional - arb
        ....:                    0, 1, 2)
        [(0.0, 6, 10), (0.5, 1, 5)]
        sage: _reduce_resolution_(((0, 2), (0.2, 1), (0.4, 0), # optional - arb
        ....:                      (0.6, 1.5), (0.8, 3)), 0, 1, 2)
        [(0.0, 0, 2), (0.5, 1.50000000000000, 3)]
        sage: _reduce_resolution_([(0, 10)],  # optional - arb
        ....:                    0, 1, 2)
        [(0.0, 10, 10)]

    TESTS::

        sage: _reduce_resolution_([(1, 10)],  # optional - arb
        ....:                    0, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: x values must be >= x_min and < x_max.
        sage: _reduce_resolution_([(-1, 10)],  # optional - arb
        ....:                    0, 1, 2)
        Traceback (most recent call last):
        ...
        ValueError: x values must be >= x_min and < x_max.
    """
    result = [None for _ in range(resolution)]
    f = (<double> resolution)/(x_max-x_min)

    for (x, y) in data:
        if x < x_min and x > x_min - 1/f:
            # Allow undershooting x_min by one pixel in order to be
            # more robust against numerical noise.
            x = x_min
        if x < x_min or x >= x_max:
            raise ValueError(
                "x values must be >= x_min and < x_max.")
        i = floor((x - x_min)*f)
        current = result[i]
        if current is None:
            result[i] = (y, y)
        else:
            result[i] = (min(y, current[0]), max(y, current[1]))
    return [(x_min + (<double> i)/f, y[0], y[1])
            for i, y in enumerate(result)
            if y is not None]


class FSMFourierComponent(SageObject):
    r"""
    Final component of a
    :class:`~sage.combinat.finite_state_machine.Transducer` and its
    associated data for computing the Fourier coefficients of the
    fluctuations of the sum of output.

    INPUT:

    - ``fsm`` -- the final component as a
      :class:`~sage.combinat.finite_state_machine.Transducer`.

    - ``parent`` -- an instance of
      :class:`~sage.combinat.fsm_fourier.FSMFourier` holding all
      relevant data for the full transducer.

    EXAMPLES::

        sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
        sage: function('f')
        f
        sage: var('n')
        n
        sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
        sage: F = FSMFourier(transducers.Recursion([ # indirect doctest; optional - arb
        ....:     f(2*n + 1) == f(n) + 1,
        ....:     f(2*n) == f(n),
        ....:     f(0) == 0],
        ....:     2, f, n))
        sage: F.components[0] # optional - arb
        <class 'sage.combinat.fsm_fourier.FSMFourierComponent'>
        sage: F.components[0].period # optional - arb
        1
    """

    def __init__(self, fsm, parent):
        r"""
        Initialize the :class:`FSMFourierComponent`.

        INPUT:

        - ``fsm`` -- the final component as a
          :class:`~sage.combinat.finite_state_machine.Transducer`.

        - ``parent`` -- an instance of :class:`FSMFourier` holding all
          relevant data for the full transducer.

        EXAMPLE::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: function('f')
            f
            sage: var('n')
            n
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: F = FSMFourier(transducers.Recursion([ # indirect doctest; optional - arb
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: F.components[0] # optional - arb
            <class 'sage.combinat.fsm_fourier.FSMFourierComponent'>
            sage: F.components[0].period # optional - arb
            1
        """
        self.fsm = fsm
        self.period = fsm.graph().period()
        self.n_states = len(self.fsm.states())
        self.parent = parent

    def _mask_(self):
        r"""
        Return a matrix whose rows consist of the incidence
        vectors of the states contained in the other final components
        of ``self.parent`` than self.

        INPUT:

        None.

        OUTPUT:

        A matrix.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: FC = FSMFourier(T).components[0] # optional - arb
            sage: FC._mask_() # optional - arb
            [0 0 1]
        """
        n = self.parent.M.nrows()
        nrows = sum(c.n_states
                    for c in self.parent.components if c != self)
        _mask_ = matrix(
            nrows, n,
            [self.parent.standard_basis[
                    self.parent.positions[state.label()]]
             for other in self.parent.components
             if other != self
             for state in other.fsm.iter_states()])
        return _mask_

    def _eigenvectors_(self, M):
        r"""
        Determine the dominant eigenvectors of the given full
        adjacency matrix where all coordinates corresponding to final
        components other than self vanish.

        INPUT:

        - ``M`` -- a matrix.

        OUTPUT:

        A list of vectors, the `k`-th entry is the eigenvector to
        the eigenvalue `q \exp(2\pi i k/p_j)` where `p_j` is
        the period of this component.

        The first eigenvector (to the positive eigenvalue) is
        normalized such that its entries sum up to one.

        This method is called by :meth:`.right_eigenvectors` with
        ``M = self.parent.M`` and by :meth:`.left_eigenvectors`
        with ``M = self.parent.M.transpose()``.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC._eigenvectors_(F.M) # optional - arb
            [(1/3, 2/3, 0)]

        .. SEEALSO::

            :meth:`.right_eigenvectors`, :meth:`.left_eigenvectors`
        """
        mask = self._mask_()
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
        r"""
        Determine the dominant right eigenvectors of the full
        adjacency matrix where all coordinates corresponding to final
        components other than self vanish.

        INPUT:

        None.

        OUTPUT:

        A list of vectors, the `k`-th entry is the eigenvector to
        the eigenvalue `q \exp(2\pi i k/p_j)` where `p_j` is
        the period of this component.

        The first eigenvector (to the positive eigenvalue) is
        normalized such that its entries sum up to one.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.right_eigenvectors() # optional - arb
            [(1/3, 2/3, 0)]

        .. SEEALSO::

            :meth:`.left_eigenvectors`
        """
        return self._eigenvectors_(self.parent.M)

    @cached_method()
    def left_eigenvectors(self):
        r"""
        Determine the dominant left eigenvectors of the full
        adjacency matrix where all coordinates corresponding to final
        components other than self vanish.

        INPUT:

        None.

        OUTPUT:

        A list of vectors, the `k`-th entry is the eigenvector to
        the eigenvalue `q \exp(2\pi i k/p_j)` where `p_j` is
        the period of this component.

        The eigenvectors are normalized such that their product
        with the corresponding right eigenvectors computed by
        :meth:`.right_eigenvectors` is `1`.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.left_eigenvectors() # optional - arb
            [(0, 3/2, 0)]

        .. SEEALSO::

            :meth:`.right_eigenvectors`
        """
        left_eigenvectors = self._eigenvectors_(self.parent.M.transpose())
        return [w/(v*w) for v, w
                in itertools.izip(self.right_eigenvectors(),
                                  left_eigenvectors)]

    @cached_method()
    def w(self):
        r"""
        Compute `w_{jk}` for `0\le k< p_j`.

        INPUT:

        None.

        OUTPUT:

        A list of vectors.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.w() # optional - arb
            [(0, 1/2, 0)]

        .. SEEALSO::

            :meth:`.left_eigenvectors`, :meth:`.right_eigenvectors`
        """
        return [(self.parent.initial_vector*v)*w for v, w
                in itertools.izip(self.right_eigenvectors(),
                                  self.left_eigenvectors())]

    @cached_method()
    def coefficient_lambda(self):
        r"""
        Compute `\lambda_j`.

        INPUT:

        None.

        OUTPUT:

        An element of a cyclotomic field representing the value
        `\lambda_j`.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.coefficient_lambda() # optional - arb
            1/2
        """
        products = [w*self.parent.ones for w in self.w()]
        assert all(e.is_zero() for e in products[1:])
        return products[0]

    @cached_method()
    def mu_prime(self):
        r"""
        Compute `\mu_{j0}'(0)`, the derivative of the dominant positive
        eigenvalue of `M(t)` corresponding to ``self`` at `t=0`.

        INPUT:

        None.

        OUTPUT:

        An element of a cyclotomic field representing the value
        `\mu'_{j0}(0)`.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.mu_prime() # optional - arb
            0
        """
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
        r"""
        Compute the constant `a_{j}`.

        INPUT:

        None.

        OUTPUT:

        An element of a cyclotomic field representing the value
        `a_j`.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 2, 0, 0), (2, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.a() # optional - arb
            0
        """
        return QQ(-self.parent.I * self.mu_prime()/self.parent.q)

    def w_ell(self, ell):
        r"""
        Compute the left eigenvector to the eigenvalue
        `\exp(2\pi i \ell/p)` to `M(0)` if the former is an eigenvalue
        of `M(0)`. Otherwise, return the zero vector.

        INPUT:

        - ``ell`` -- an integer.

        OUTPUT:

        A vector over a complex interval field or the zero vector.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 3, 0, 0), (2, 3, 1, 0),
            ....:                 (3, 2, 0, 0), (3, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2, 3]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.w_ell(0) # optional - arb
            (0, 0.50000000000000000?, 0, 0)
            sage: FC.w_ell(1) # optional - arb
            (0, 0, 0, 0)
        """
        if self.parent.common_period.divides(ell*self.period):
            k = self.period*ell/self.parent.common_period % self.period
            return vector(self.parent.field_to_CIF(c) for c in self.w()[k])
        else:
            return vector(0 for _ in self.parent.ones)

    def v_prime(self, k):
        r"""
        Compute the derivative `v_{jk}'(0)`, the right
        eigenvector to the `k`-th dominant eigenvalue.

        INPUT:

        - ``k`` -- a non-negative integer.

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 3, 0, 0), (2, 3, 1, 0),
            ....:                 (3, 2, 0, 0), (3, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2, 3]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.v_prime(0) # optional - arb
            (0, 0, 0, 0)
        """
        M = self.parent.M
        mask = self._mask_()
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

    def w_prime(self, k):
        r"""
        Compute the derivative `w_{jk}'(0)`, the left
        eigenvector to the `k`-th dominant eigenvalue.

        INPUT:

        - ``k`` -- a non-negative integer.

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 1, 0, 0), (0, 2, 1, 0),
            ....:                 (1, 1, 0, 0), (1, 1, 1, 0),
            ....:                 (2, 3, 0, 0), (2, 3, 1, 0),
            ....:                 (3, 2, 0, 0), (3, 2, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0, 1, 2, 3]) # optional - arb
            sage: F = FSMFourier(T) # optional - arb
            sage: FC = F.components[0] # optional - arb
            sage: FC.w_prime(0) # optional - arb
            (0, 0, 0, 0)
        """
        M = self.parent.M
        mask = self._mask_()
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
                                     [matrix(self.v_prime(k))],
                                     [0*mask]],
                                    subdivide=False) * eigenvector_left
        left_prime = S.solve_right(right_side)
        w_prime = self.parent.initial_vector*self.v_prime(k)*eigenvector_left \
            + self.parent.initial_vector*eigenvector_right*left_prime
        return w_prime

cdef class FSMFourierCache(SageObject):
    r"""
    Compute and cache `\mathbf{b}(r)` for increased performance
    and compute partial sums of the Dirichlet series `\mathbf{H}(s)`.

    INPUT:

    - ``parent`` -- a :class:`FSMFourier` instance.

    OUTPUT:

    None.

    EXAMPLES::

        sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
        sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
        sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, 1)],
        ....:                initial_states=[0],
        ....:                final_states=[0])
        sage: FSMFourier(T).cache # optional - arb
        <type 'sage.combinat.fsm_fourier.FSMFourierCache'>
    """
    cdef acb_mat_t *bb
    cdef unsigned long bb_computed
    cdef unsigned long bb_allocated
    cdef acb_mat_t *Delta_epsilon_ones
    cdef acb_mat_t *M_epsilon
    cdef unsigned long q
    cdef unsigned long precision
    cdef unsigned long n
    cdef object parent

    def __init__(self, parent):
        """
        Initialize the class.

        INPUT:

        - ``parent`` -- a :class:`FSMFourier` instance.

        OUTPUT:

        None.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: FSMFourier(T).cache # optional - arb
            <type 'sage.combinat.fsm_fourier.FSMFourierCache'>
        """
        cdef long epsilon

        self.parent = parent
        self.n = self.parent.M.nrows()
        self.precision = self.parent.CIF.precision()
        self.q = self.parent.q

        MS_n = MatrixSpace(self.parent.CIF, self.n, self.n)
        MS_1 = MatrixSpace(self.parent.CIF, self.n, 1)

        self.Delta_epsilon_ones = <acb_mat_t*> malloc(
            self.q * sizeof(acb_mat_t))
        self.M_epsilon = <acb_mat_t*> malloc(
            self.q * sizeof(acb_mat_t))

        for epsilon in range(self.q):
            acb_mat_init(self.M_epsilon[epsilon], self.n, self.n)
            acb_mat_init(self.Delta_epsilon_ones[epsilon], self.n, 1)

        for epsilon in range(self.q):
            matrix_to_acb_mat(self.M_epsilon[epsilon],
                              MS_n(parent.M_epsilon[epsilon]))

            matrix_to_acb_mat(self.Delta_epsilon_ones[epsilon],
                              MS_1(parent.Delta_epsilon[epsilon] *
                                   parent.ones))
        self.bb_allocated = self.q*1024
        self.bb = <acb_mat_t*> malloc(
            self.bb_allocated * sizeof(acb_mat_t))
        self.bb_computed = 0

    def __dealloc__(self):
        """
        Deallocate memory.
        """
        cdef unsigned long epsilon

        if self.M_epsilon != NULL:
            for epsilon in range(self.q):
                if self.M_epsilon[epsilon] != NULL:
                    acb_mat_clear(self.M_epsilon[epsilon])
            free(self.M_epsilon)

        if self.Delta_epsilon_ones != NULL:
            for epsilon in range(self.q):
                if self.Delta_epsilon_ones[epsilon] != NULL:
                    acb_mat_clear(self.Delta_epsilon_ones[epsilon])
            free(self.Delta_epsilon_ones)

        if self.bb != NULL:
            for epsilon in range(self.bb_computed):
                if self.bb[epsilon] != NULL:
                    acb_mat_clear(self.bb[epsilon])
            free(self.bb)

    cdef void compute_b(self, unsigned long M):
        r"""
        Compute and store `\mathbf{b}(r)` for `0\le r< M`.

        INPUT:

        - ``M`` -- a non-negative integer.

        OUTPUT:

        None.
        """
        cdef unsigned long r
        cdef unsigned long epsilon
        cdef unsigned long R
        cdef unsigned long new_M
        cdef unsigned long n
        cdef bint resize

        if M < self.bb_computed:
            return

        new_M = (M // self.q) * self.q
        if new_M < M:
            new_M += self.q

        resize = False
        while new_M > self.bb_allocated:
            self.bb_allocated *= self.q
            resize = True

        if resize:
            self.bb = <acb_mat_t *> realloc(
                self.bb,
                self.bb_allocated * sizeof(acb_mat_t))

        for r in range(self.bb_computed//self.q, new_M//self.q):
            for epsilon in range(self.q):
                R = r*self.q + epsilon
                acb_mat_init(self.bb[R], self.n, 1)
                if R == 0:
                    matrix_to_acb_mat(
                        self.bb[0],
                        matrix(self.parent.CIF,
                               self.parent.b0()).transpose())
                else:
                    acb_mat_mul(self.bb[R],
                                self.M_epsilon[epsilon],
                                self.bb[r],
                                self.precision)
                    acb_mat_add(self.bb[R],
                                self.bb[R],
                                self.Delta_epsilon_ones[epsilon],
                                self.precision)
        self.bb_computed = new_M

    cpdef b(self, unsigned long r):
        r"""
        Compute or retrieve `\mathbf{b}(r)`.

        INPUT:

        - ``r`` -- a non-negative integer.

        OUTPUT:

        A vector of :class:`ComplexIntervalFieldElement`.

        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: function('f')
            f
            sage: var('n')
            n
            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: cache = F.cache # optional - arb
            sage: [cache.b(r) for r in range(3)] # optional - arb
            [(0, 1, 1), (1, 2, 1), (1, 2, 2)]
        """
        self.compute_b(r+1)
        return acb_mat_to_matrix(self.bb[r], self.parent.CIF).column(0)

    cpdef fluctuation_empirical(self, unsigned long start,
                                unsigned long end):
        r"""
        Compute empirical fluctuation
        `\frac1n\sum_{0\le k<n} T(k) - e_T\log_q(n)` of the sum of
        outputs for `n` in the given interval.

        INPUT:

        - ``start`` -- a positive integer

        - ``end`` -- a non-negative integer

        OUTPUT:

        A list of pairs of doubles, the first entry is `\log_q(n)`,
        the second the empirical fluctuation at this point.

        EXAMPLES::

            sage: function('f')
            f
            sage: var('n')
            n
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = transducers.Recursion([
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n)
            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: F = FSMFourier(T) # optional - arb
            sage: F.cache.fluctuation_empirical(10, 12) # optional - arb; tolerance 1e-12
            [(3.3219280948873626, -0.1609640474436813),
             (3.4594316186372978, -0.18426126386410324)]

        The following transducer has sum of output `1` for all input. ::

            sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: T.state(0).final_word_out = 1
            sage: [T(i.bits()) for i in srange(3)]
            [[1], [0, 1], [0, 0, 1]]
            sage: F = FSMFourier(T) # optional - arb
            sage: F.cache.fluctuation_empirical(1, 4) # optional - arb; tolerance 1e-12
            [(0.0, 1.0), (1.0, 1.0), (1.5849625007211563, 1.0)]

        The following transducer has sum of output `\lfloor\log n\rfloor` for input `n`. ::

            sage: T = Transducer([(0, 0, 0, 1), (0, 0, 1, 1)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: T.state(0).final_word_out = -1
            sage: [sum(T(i.bits())) for i in srange(5)]
            [-1, 0, 1, 1, 2]
            sage: F = FSMFourier(T) # optional - arb
            sage: F.e_T # optional - arb
            1
            sage: F.cache.fluctuation_empirical(1, 4) # optional - arb; tolerance 1e-12
            [(0.0, -1.0),
             (1.0, -1.5),
             (1.5849625007211563, -1.5849625007211563)]

        The parameter ``start`` must be positive::

            sage: F.cache.fluctuation_empirical(0, 2) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: start must be positive.
        """
        cdef long initial
        cdef double *values
        cdef double sum
        cdef long i
        cdef RealIntervalFieldElement current
        cdef double e_T
        cdef list result
        cdef long precision
        cdef long q
        cdef double log_q

        initial = self.parent.positions[
            self.parent.transducer.initial_states()[0].label()]
        sum = 0
        e_T = self.parent.e_T
        q = self.parent.q
        self.compute_b(end)
        values = <double*> malloc((end-start) * sizeof(double))
        current = RealIntervalField()(0)
        precision = current.parent().precision()

        if start <= 0:
            raise ValueError("start must be positive.")

        for i in range(end):
            if i >= start:
                values[i - start] = sum/i - e_T * log(i, base=q)
            arb_to_mpfi(
                current.value,
                &acb_mat_entry(self.bb[i],
                              initial,
                              0).real,
                precision)
            sum += current.center()

        log_q = log(q)
        result = [(log(<double> i)/log_q, values[i-start])
                  for i in range(start, end)]

        free(values)
        return result

    cpdef fluctuation_fourier(self,
                              double start,
                              double end,
                              unsigned long resolution):
        r"""
        Computes the periodic fluctuation by using Fourier coefficients,
        using ``resolution`` points per period.

        INPUT:

        - ``start`` -- a double, start of the interval

        - ``end`` -- a double, end of the interval

        - ``resolution`` -- a positive integer, number of points per
          period interval.

        OUTPUT:

        A list of pairs of doubles, consisting of pairs `(x, \Phi_1(x))`.

        The periodic fluctuation is approximated by the ``resolution``-th
        Fourier polynomial.

        .. WARNING::

            The output of the transducer is assumed to be real, i.e.,
            the `-j`-th Fourier coefficient is the complex conjugate
            of the `j`-th Fourier coefficient.

        EXAMPLES:

        The following transducer has sum of output `1` for all input. ::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: T = Transducer([(0, 0, 0, 0), (0, 0, 1, 0)],
            ....:                initial_states=[0],
            ....:                final_states=[0])
            sage: T.state(0).final_word_out = 1
            sage: [T(i.bits()) for i in srange(3)]
            [[1], [0, 1], [0, 0, 1]]
            sage: F = FSMFourier(T) # optional - arb
            sage: F.cache.fluctuation_fourier(2, 3, 2) # optional - arb; rel tolerance 1e-10
            [(2.0, 0.9999999999999978), (2.5, 0.9999999999999991)]
        """
        cdef double c0
        cdef long i
        cdef unsigned long period

        import sage.parallel.ncpus
        import sage.gsl.fft

        period = self.parent.common_period

        c0 = self.parent.FourierCoefficient(0).real().center()
        self.parent.FourierCoefficient.precompute(
            range(resolution),
            sage.parallel.ncpus.ncpus())

        fft = sage.gsl.fft.FastFourierTransform(resolution)
        for j in range(resolution):
            fft[j] = self.parent.FourierCoefficient(j).center()

        fft.backward_transform()

        return [((<double>i*period)/resolution, 2*fft[i % resolution][0] - c0)
                for i in range(ceil(start/period * resolution),
                               ceil(end/period * resolution))]

    cdef FreeModuleElement_generic_dense partial_dirichlet(
        self,
        unsigned long start,
        unsigned long end,
        ComplexIntervalFieldElement s):
        r"""
        Compute `\sum_{\mathit{start}\le r<\mathit{end}} \mathbf{b}(r) r^{-s}`.

        INPUT:

        `start` -- a positive integer.

        `end` -- a positive integer.

        `s` -- a :class:`ComplexIntervalFieldElement`.

        OUTPUT:

        A vector of :class:`ComplexIntervalFieldElement`.
        """
        cdef unsigned long r;
        cdef acb_mat_t result;
        cdef acb_t scalar;
        cdef acb_t minuss;
        cdef FreeModuleElement_generic_dense result_vector

        self.compute_b(end)

        acb_mat_init(result, self.n, 1)
        acb_init(scalar)
        acb_init(minuss)

        ComplexIntervalFieldElement_to_acb(minuss, -s)
        for r in range(end-1, start-1, -1):
            acb_set_ui(scalar, r)
            acb_pow(scalar, scalar, minuss, self.precision)
            acb_mat_scalar_addmul_acb(result, self.bb[r],
                                      scalar, self.precision)
        result_vector = acb_mat_to_matrix(result, self.parent.CIF).column(0)

        acb_clear(minuss)
        acb_clear(scalar)
        acb_mat_clear(result)

        return result_vector

class FSMFourier(SageObject):
    """
    Compute Fourier coefficients of the periodic fluctuation of the
    summatory function of the sum of output of transducers.

    INPUT:

    - ``transducer`` -- a
      :class:`~sage.combinat.finite_state_machine.Transducer`.

    - ``CIF`` -- a
      :class:`~sage.rings.complex_interval_field.ComplexIntervalField_class`,
      indicating the precision for the floating point operations.

    The object stores various data (in particular, eigenvectors and
    their derivatives, final components) and provides the method
    :meth:`FourierCoefficient` to actually compute single Fourier
    coefficients.

    EXAMPLES:

    -   Binary sum of digits::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: function('f')
            f
            sage: var('n')
            n
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: F.common_period # optional - arb
            1
            sage: [FC] = F.components # optional - arb
            sage: FC.w() # optional - arb
            [(1)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            1/2
            sage: FC.a() # optional - arb
            1/2
            sage: F.M # optional - arb
            [2]
            sage: F.M_epsilon # optional - arb
            [[1], [1]]
            sage: F.Delta # optional - arb
            [1]
            sage: F.Delta_epsilon # optional - arb
            [[0], [1]]
            sage: F.C_0 # optional - arb
            1
            sage: F.C_1 # optional - arb
            1

    -   NAF::

            sage: F = FSMFourier(transducers.Recursion([  # optional - arb
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: F.common_period # optional - arb
            1
            sage: [FC] = F.components # optional - arb
            sage: FC.w() # optional - arb
            [(1/3, 1/3, 1/3)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            1/3
            sage: FC.a() # optional - arb
            1/3
            sage: F.FourierCoefficient(0) # optional - arb
            0.4478541793943? + 0.?e-14*I
            sage: F.FourierCoefficient(42) # optional - arb
            -9.0701442?e-6 + 0.0001189561459?*I

    -   Abelian complexity of the paperfolding sequence::

            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n+2) == f(2*n+1)+1,
            ....:     f(16*n+1) == f(8*n+1),
            ....:     f(16*n+5) == f(4*n+1)+2,
            ....:     f(16*n+11) == f(4*n+3)+2,
            ....:     f(16*n+15) == f(2*n+2)+1,
            ....:     f(1) == 2, f(0) == 0]
            ....:     + [f(16*n+jj) == f(2*n+1)+2 for jj in [3,7,9,13]],
            ....:     2, f, n))
            sage: [FC] = F.components # optional - arb
            sage: F.common_period # optional - arb
            1
            sage: FC.w() # optional - arb
            [(0, 0, 3/13, 2/13, 2/13, 2/13, 1/13, 1/13, 1/13, 1/13)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            8/13
            sage: FC.a() # optional - arb
            8/13
            sage: F.FourierCoefficient(0) # optional - arb
            1.5308151287593? + 0.?e-13*I
            sage: F.FourierCoefficient(42) # optional - arb
            9.38016?e-7 + 0.000156956885?*I

    -   Artificial example, one-periodic, 2 states::

            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(4*n) == f(2*n)+0,
            ....:     f(4*n+2) == f(n)+1,
            ....:     f(2*n+1) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: [FC] = F.components # optional - arb
            sage: F.common_period # optional - arb
            1
            sage: FC.w() # optional - arb
            [(1/2, 1/2)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            1/4
            sage: FC.a() # optional - arb
            1/4
            sage: F.FourierCoefficient(0) # optional - arb
            -0.2031981160232? + 0.?e-14*I
            sage: F.FourierCoefficient(42) # optional - arb
            -0.0000280287200? + 0.0000741203738?*I

    -   Artificial example, period 3::

            sage: T = Transducer([
            ....:     (0, 1, 0, []), (0, 2, 1, []), (1, 3, 0, []),
            ....:     (1, 4, 1, []), (2, 5, 0, []), (2, 6, 1, []),
            ....:     (3, 6, 0, 3), (3, 6, 1, 1), (4, 6, 0, 2),
            ....:     (4, 6, 1, -1), (5, 3, 0, 5), (5, 4, 1, 1),
            ....:     (6, 5, 0, 2), (6, 5, 1, 0)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2, 3, 4, 5, 6])
            sage: T.state(0).final_word_out = 0
            sage: T.state(1).final_word_out = 0
            sage: T.state(2).final_word_out = [5, 0]
            sage: T.state(3).final_word_out = 0
            sage: T.state(4).final_word_out = [2, 2, 5, 0]
            sage: T.state(5).final_word_out = [5, 0]
            sage: T.state(6).final_word_out = [2, 5, 0]
            sage: F = FSMFourier(T) # optional - arb
            sage: [FC] = F.components # optional - arb
            sage: F.common_period # optional - arb
            3
            sage: FC.w() # optional - arb
            [(0, 0, 0, 1/6, 1/6, 1/3, 1/3),
            (0, 0, 0, 1/24*zeta12^2 - 1/24, 1/24*zeta12^2 - 1/24, -1/12*zeta12^2, 1/12),
            (0, 0, 0, -1/24*zeta12^2, -1/24*zeta12^2, 1/12*zeta12^2 - 1/12, 1/12)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            7/4
            sage: FC.a() # optional - arb
            7/4
            sage: F.FourierCoefficient(0) # optional - arb
            0.879448677390? + 0.?e-13*I
            sage: F.FourierCoefficient(42) # optional - arb
            0.00033040004? + 0.002040223214?*I
            sage: F.FourierCoefficient(43) # optional - arb
            0.000178962837? + 3.877474?e-6*I

    -   Artificial example, period 2, vanishing `\mathbf{w}`-vector::

            sage: T = Transducer(
            ....:     [(0, 1, 0, []), (0, 2, 1, []), (1, 2, 0, 1),
            ....:      (1, 2, 1, 3), (2, 1, 0, 2), (2, 1, 1, -1)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2])
            sage: T.state(0).final_word_out = 0
            sage: T.state(1).final_word_out = 0
            sage: T.state(2).final_word_out = [2, 0]
            sage: F = FSMFourier(T) # optional - arb
            sage: [FC] = F.components # optional - arb
            sage: F.common_period # optional - arb
            2
            sage: FC.w() # optional - arb
            [(0, 1/2, 1/2), (0, 0, 0)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            5/4
            sage: FC.a() # optional - arb
            5/4
            sage: F.FourierCoefficient(0) # optional - arb
            -1.5912428334793? + 0.?e-13*I
            sage: F.FourierCoefficient(42) # optional - arb
            -0.000265095706? - 0.000100293638?*I
            sage: F.FourierCoefficient(43) # optional - arb
            0

    -   Artificial example with two final components of periods `2`
        and `3`, respectively::

            sage: F = FSMFourier(Transducer([(0, 1, 0, 1), (1, 2, 0, 1), # optional - arb
            ....:     (2, 1, 0, 2), (1, 2, 1, 0), (2, 1, 1, 2),
            ....:     (0, -1, 1, 1), (-1, -2, 1, 1), (-2, -3, 1, 1),
            ....:     (-3, -1, 1, 1), (-1, -2, 0, 2), (-2, -3, 0, 1),
            ....:     (-3, -1, 0, 3)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2, -3, -2, -1]))
            sage: [FC0, FC1] = F.components # optional - arb
            sage: FC0.period # optional - arb
            3
            sage: FC1.period # optional - arb
            2
            sage: F.common_period # optional - arb
            6
            sage: FC0.w() # optional - arb
            [(0, 0, 0, 1/6, 1/6, 1/6),
              (0, 0, 0, 1/6, 1/6*zeta12^2 - 1/6, -1/6*zeta12^2),
              (0, 0, 0, 1/6, -1/6*zeta12^2, 1/6*zeta12^2 - 1/6)]
            sage: FC1.w() # optional - arb
            [(0, 1/4, 1/4, 0, 0, 0), (0, -1/4, 1/4, 0, 0, 0)]
            sage: FC0.coefficient_lambda() # optional - arb
            1/2
            sage: FC1.coefficient_lambda() # optional - arb
            1/2
            sage: F.e_T # optional - arb
            11/8
            sage: [FC0.a(), FC1.a()] # optional - arb
            [3/2, 5/4]
            sage: F.FourierCoefficient(0) # optional - arb
            -2.1863500631078? + 0.?e-13*I
            sage: F.FourierCoefficient(42) # optional - arb
            -0.00242220796? + 0.00103325145?*I
            sage: F.FourierCoefficient(43) # optional - arb
            0
            sage: F.FourierCoefficient(44) # optional - arb
            0.000214208380? - 0.00045396681?*I
            sage: F.FourierCoefficient(45) # optional - arb
            -0.0005573382643? - 0.000553163835?*I

    -   Ternary sum of digits, with high precision::

            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(3*n + 2) == f(n) + 2,
            ....:     f(3*n + 1) == f(n) + 1,
            ....:     f(3*n) == f(n),
            ....:     f(0) == 0],
            ....:     3, f, n), ComplexIntervalField(200))
            sage: F.common_period # optional - arb
            1
            sage: F.e_T # optional - arb
            1
            sage: F.FourierCoefficient(0) # optional - arb
            -0.23733142706319409139996985552948551198247848839875752890?
            + 0.?e-57*I
            sage: F.FourierCoefficient(42) # optional - arb
            0.000151640984887988283825448153753614629500867823400070268?
            + 0.000054159306153387140043462960871097923443501309058799584?*I
    """

    common_period = None
    """`p`, the least common multiple of the periods
    of the final components"""

    e_T = None
    """The constant `e_{\mathcal{T}}`, the coefficient of
    the main term of the expectation."""

    M = None
    """`M(0)`, the adjacency matrix of the transducer."""

    M_epsilon = None
    """`(M_\varepsilon)_{0\le \varepsilon < q}`, the list of
    transition matrices."""

    Delta = None
    """`\Delta(0)`, the output matrix."""

    Delta_epsilon = None
    """`(Delta_\varepsilon)_{0\le\varepsilon <q}`, the list of
    partial of partial output matrices."""

    C_0 = None
    """`C_0 = \max\{\|\mathbf{b}(r)\|_\infty: 0\le r<q\}`, used
    in the bound of the Dirichlet series."""

    C_1 = None
    """`C_1 = \max\{\|\Delta_\varepsilon\|_\infty: 0\le \varepsilon<q\}`,
    used in the bound of the Dirichlet series."""

    components = None
    """A list of :class:`FSMFourierComponent`, representing the
    final components."""

    def __init__(self, transducer, CIF=ComplexIntervalField()):
        r"""
        Initialize the the common data needed for the computation of all
        Fourier coefficients of the periodic fluctuation of the sum of
        output.

        INPUT:

        - ``transducer`` -- a
          :class:`~sage.combinat.finite_state_machine.Transducer`.

        - ``CIF`` -- a
          :class:`~sage.rings.complex_interval_field.ComplexIntervalField_class`,
          indicating the precision for the floating point operations.

        OUTPUT:

        Nothing.

        EXAMPLES:

        Binary sum of digits::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: function('f')
            f
            sage: var('n')
            n
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: F = FSMFourier(transducers.Recursion([  # optional - arb
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: F.common_period # optional - arb
            1
            sage: [FC] = F.components # optional - arb
            sage: FC.w() # optional - arb
            [(1)]
            sage: FC.coefficient_lambda() # optional - arb
            1
            sage: F.e_T # optional - arb
            1/2
            sage: FC.a() # optional - arb
            1/2
            sage: F.M # optional - arb
            [2]
            sage: F.M_epsilon # optional - arb
            [[1], [1]]
            sage: F.Delta # optional - arb
            [1]
            sage: F.Delta_epsilon # optional - arb
            [[0], [1]]
            sage: F.C_0 # optional - arb
            1
            sage: F.C_1 # optional - arb
            1
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
        self.CIF = CIF
        self.positions = dict((state.label(), j)
                         for j, state in enumerate(transducer.iter_states()))
        self.q = len(transducer.input_alphabet)
        self.Y = PolynomialRing(QQ, 'Y').gen()
        self.ones = vector(1 for _ in transducer.iter_states())



        self.components = [FSMFourierComponent(c, self)
                           for c in transducer.final_components()]
        self.common_period = lcm([c.period
                                  for c in self.components])
        self.field = CyclotomicField(lcm(4, self.common_period))
        self.alpha = self.field.zeta(self.common_period)
        self.I = self.field.zeta(4)
        self.field_to_CIF = self.field.hom(
            [CIF.zeta(lcm(4, self.common_period))], check=False)
        assert self.I**2 == -1
        assert self.alpha**self.common_period == 1
        assert all(self.alpha**j != 1 for j in range(1, self.common_period))
        assert self.field_to_CIF(self.I).overlaps(CIF(0, 1))
        assert self.field_to_CIF(self.alpha).overlaps(CIF(exp(2*CIF.pi()*CIF(0, 1)/self.common_period)))
        self.M = transducer.adjacency_matrix(entry=lambda t: 1)
        self.standard_basis = VectorSpace(self.field, self.M.nrows()).basis()

        if len(transducer.initial_states()) != 1:
            raise NotImplementedError(
                "Transducer does not have a unique initial state.")
        self.initial_vector = self.standard_basis[self.positions[
                transducer.initial_states()[0].label()]]

        self.e_T = sum(c.a()*c.coefficient_lambda()
                       for c in self.components)

        n0 = var('n0')
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

        self.cache = FSMFourierCache(self)

        self.C_0 = max(self.cache.b(r).norm(infinity)
                       for r in range(self.q))
        self.C_1 = max(infinity_matrix_norm(d)
                       for d in self.Delta_epsilon)

        self.error_acceptable = ZZ(2)**(-2*CIF.precision())


    @cached_method
    def b0(self):
        r"""
        Compute `\mathbf{b}(0)`, the vector of final output sums.

        INPUT:

        None.

        OUTPUT:

        A vector whose ``s``-th component is the sum of the output of
        ``self`` when reading the empty word starting in the state ``s``.


        EXAMPLES::

            sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: function('f')
            f
            sage: var('n')
            n
            sage: T = transducers.Recursion([
            ....:     f(4*n + 1) == f(n) + 1,
            ....:     f(4*n + 3) == f(n + 1) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n)
            sage: FSMFourier(T).b0() # optional - arb
            (0, 1, 1)
        """
        from sage.modules.free_module_element import vector

        return vector([sum(self.transducer([], initial_state=s))
                       for s in self.transducer.iter_states()])

    def _H_m_rhs_(self, s, m, remove_poles=False):
        r"""
        Compute `(I - q^{-s}M)\mathbf{H}_m(s)`.

        INPUT:

        -   ``s`` -- a
            :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` element

        -   ``m`` -- non-negative integer

        -   ``remove_poles`` -- (default: ``False``) if ``True`` and ``s
            == 1``, the pole of the zeta functions is removed.

        The function has a pole for ``s == 1``. If ``remove_poles`` is ``True``,
        the zeroth coefficient of the Laurent series at ``s == 1`` is returned
        instead.

        OUTPUT:

        A
        :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement`
        element.

        EXAMPLES:

        -   One state, always ``1``: this corresponds to the Riemann zeta
            function (minus the first 99 summands, multiplied by 2). ::

                sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
                sage: function('f')
                f
                sage: var('n')
                n
                sage: F = FSMFourier(transducers.Recursion([ # optional - arb
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     2, f, n))
                sage: sage.combinat.finite_state_machine.FSMOldProcessOutput = False
                sage: F._H_m_rhs_(CIF(2), 100) # optional - arb
                (0.0050250833316668? + 0.?e-18*I)

            Evaluated at `s = 1 + 2k \pi i/\log 2` for `k\neq 0`, the
            result must be zero, as the Riemann zeta function does not
            have a pole there and the first factor therefore annihilates
            the result::

                sage: result = [F._H_m_rhs_(CIF(1+2*k*pi*I/log(2)), 30)[0] for  # optional - arb
                ....:      k in range(1, 2)]
                sage: all(value.abs().absolute_diameter() < 1e-12  # optional - arb
                ....:      for value in result)
                True
                sage: all(value.contains_zero() for value in result) # optional - arb
                True

        -   One state, one for every input digit, but subtract one at the
            end. This corresponds to the function `L`. Evaluated at `s =
            1 + 2k \pi i/\log 2` for `k\neq 0`, the result must equal
            `\log 2/(2 k \pi i)`, because the `1 \times 1`-matrix is
            singular so that we actually compute the residue at this
            point, multiplied by `\log 2`. The value for `m` does not
            matter, as the difference is analytic and thus does not
            contribute to the residue. ::

                sage: T = Transducer([(0, 0, 0, 1), (0, 0, 1, 1)],
                ....:     initial_states=[0], final_states=[0])
                sage: T.state(0).final_word_out = -1
                sage: F = FSMFourier(T) # optional - arb

            We check that the result agrees with the known values::

                sage: all(F._H_m_rhs_(CIF(1 + 2*k*pi*I/log(2)), # optional - arb
                ....:         30)[0].overlaps(CIF(log(2)/(2*pi*I*k)))
                ....:     for k in range(1, 2))
                True

        TESTS::

            sage: F._H_m_rhs_(CIF(1), 20, remove_poles=False) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: remove_poles must be set if and only if s == 1.
            sage: F._H_m_rhs_(CIF(2), 20, remove_poles=True) # optional - arb
            Traceback (most recent call last):
            ...
            ValueError: remove_poles must be set if and only if s == 1.
        """
        cdef FSMFourierCache cache

        verbose("_H_m_rhs_(%s, %s)" % (s, m), level=1)
        if remove_poles != (s == 1):
            raise ValueError("remove_poles must be set if "
                             "and only if s == 1.")

        sigma = s.real()
        CIF = s.parent()
        RIF = sigma.parent()

        def taylor_error(N, sigma, x):
            return RIF_E.min(
                RIF_E.max(RIF(1), 1/(1 + x)**(sigma + N)),
                RIF_E.max(RIF(N), N/(1 + x)**(sigma + 1)))


        q = ZZ(self.q)
        log_q = log(RIF(q))
        log_m_1 = log(RIF(m - 1))

        cache = self.cache
        result = cache.partial_dirichlet(
            m, q*m, s)

        if remove_poles and s == 1:
            result += q**(-s) * sum(
                D * self.ones * (-RIF(ZZ(epsilon)/q + int(epsilon == 0)).psi()
                             - sum((k + ZZ(epsilon)/q + int(epsilon==0))**(-1)
                                   for k in range(0, m-int(epsilon==0))))
                for epsilon, D in enumerate(self.Delta_epsilon))
        else:
            result += q**(-s) * sum(
                D * self.ones * _hurwitz_zeta_(s, ZZ(epsilon)/q, m,
                                               self.error_acceptable)
                for epsilon, D in enumerate(self.Delta_epsilon))

        N = 1
        factor = -s * q**(-s - 1)

        while True:
            #assert factor.overlaps(sage.rings.arith.binomial(-s, N) * q**(-s - N))

            if result.is_zero():
                error_acceptable = 0
            else:
                rounding_error = sum(
                    c.real().upper() + c.imag().upper()
                    -c.real().lower()-c.imag().lower()
                    for c in result)
                error_acceptable = max(self.error_acceptable,
                                       rounding_error/8)

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

        ``s`` -- a :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` element

        OUTPUT:

        A :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` element.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: function('f')
            f
            sage: var('n')
            n
            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(2*n + 1) == f(n),
            ....:     f(2*n) == f(n),
            ....:     f(0) == 1],
            ....:     2, f, n))
            sage: w = F.components[0].w(); w # optional - arb
            [(1)]
            sage: F._w_H_Res_(w[0], CIF(1)) # optional - arb
            1.00000000000000? + 0.?e-15*I
        """

        CIF = s.parent()
        log_q = log(CIF(self.q))
        m = s.abs().upper().ceil() + s.parent().precision()

        if s == 1:
            return w * self._H_m_rhs_(s, m, remove_poles=True)/log_q \
                - self.e_T/2
        else:
            return w * self._H_m_rhs_(s, m)/log_q

    @cached_method(key=lambda self, s, m: (s.real().lower(),
                                           s.real().upper(),
                                           s.imag().lower(),
                                           s.imag().upper(), m))
    def _H_m_(self, s, m):
        r"""
        Compute `\mathbf{H}_m(s)`.

        INPUT:

        -   ``s`` -- a :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` element

        -   ``m`` -- non-negative integer

        OUTPUT:

        A :class:`~sage.rings.complex_interval.ComplexIntervalFieldElement` element.

        EXAMPLES:

        -   One state, always ``1``: this corresponds to the Riemann
            zeta function (minus the first 99 summands). ::

                sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
                sage: function('f')
                f
                sage: var('n')
                n
                sage: T = transducers.Recursion([
                ....:     f(2*n + 1) == f(n),
                ....:     f(2*n) == f(n),
                ....:     f(0) == 1],
                ....:     2, f, n) # optional - arb
                sage: FSMFourier(T)._H_m_(CIF(2), 100) # optional - arb
                (0.010050166663334? + 0.?e-18*I)
        """
        n = self.M.nrows()

        return (matrix.identity(n) - ZZ(self.q)**(-s)*self.M).solve_right(
            self._H_m_rhs_(s, m))

    @cached_method
    def FourierCoefficient(self, ell):
        """
        Compute the `\ell`-th Fourier coefficient of the fluctuation
        of the summatory function of the sum of outputs.

        INPUT:

        -   ``ell`` -- an integer.

        OUTPUT:

        A :class:`ComplexIntervalFieldElement` element.

        EXAMPLES:

        -   Binary sum of digits::

                sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
                sage: function('f')
                f
                sage: var('n')
                n
                sage: F = FSMFourier(transducers.Recursion([ # optional - arb
                ....:     f(2*n + 1) == f(n) + 1,
                ....:     f(2*n) == f(n),
                ....:     f(0) == 0],
                ....:     2, f, n))
                sage: def FourierCoefficientDelange(k):
                ....:     if k == 0:
                ....:         return 1/(2*log(2))*(log(2*pi)-1)-3/4
                ....:     return CC(I/(2*k*pi)*(1 + 2*k*pi*I/log(2))^(-1)
                ....:            *zeta(CC(2*k*pi*I/log(2))))
                sage: all(FourierCoefficientDelange(k) # optional - arb
                ....:     in F.FourierCoefficient(k)
                ....:     for k in [0, 1, 42])
                True
        """

        CIF = self.CIF
        log_q = CIF(log(CIF(self.q)))
        chi_ell =  CIF(2*ell*CIF.pi()*self.I / (self.common_period*log_q))

        w = sum(c.w_ell(ell) for c in self.components)
        if any(w):
            result = 1/(1+chi_ell) * self._w_H_Res_(w, 1+chi_ell)
        else:
            result = CIF(0)

        if ell == 0:
            result += -self.e_T/log_q \
                - self.I*sum(c.w_prime(0)*self.ones
                        for c in self.components)
        return result

    def plot_fluctuation(self, start, end, resolution):
        r"""
        Plot periodic fluctuation.

        INPUT:

        - ``start`` -- a double, start of the plotting interval.

        - ``end`` -- a double, end of the plotting interval.

        - ``resolution`` -- a positive integer, number of points per
          period interval.

        OUTPUT:

        A plot.

        .. WARNING::

            The output of the transducer is assumed to be real, i.e.,
            the `-j`-th Fourier coefficient is the complex conjugate
            of the `j`-th Fourier coefficient.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: function('f')
            f
            sage: var('n')
            n
            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: F.plot_fluctuation(8, 9, 2) # optional - arb
            Graphics object consisting of 2 graphics primitives

        .. SEEALSO::

            :meth:`plot_fluctuation_asymptote`
        """
        from sage.plot.plot import list_plot

        empirical = self.cache.fluctuation_empirical(
            self.q ** start,
            self.q ** end)

        fourier = self.cache.fluctuation_fourier(start, end, resolution)

        return list_plot(empirical, color='blue') + \
            list_plot(fourier, plotjoined=True, color='red')


    def plot_fluctuation_asymptote(self, prefix, start, end,
                                   width=4, resolution=1200):
        r"""
        Provide `asymptote <http://asymptote.sourceforge.net>`_  files
        for the periodic fluctuation.

        INPUT:

        - ``prefix`` -- a string, filename for asymptote files. Three
          files are generated: a ``.asy``, an ``-empirical.dat`` and a
          ``-fourier.dat`` file; the latter two are read by asymptote
          when compiling the ``.asy`` file.

        - ``start`` -- a double, start of the plotting interval.

        - ``end`` -- a double, end of the plotting interval.

        - ``width`` -- (default: 4) a double giving
          the width in inches for the plotting area (excluding
          labels).

        - ``resolution`` -- (default: 1200) a positive integer, number
          of points per inch.

        OUTPUT:

        None.

        .. WARNING::

            The output of the transducer is assumed to be real, i.e.,
            the `-j`-th Fourier coefficient is the complex conjugate
            of the `j`-th Fourier coefficient.

        EXAMPLES::

            sage: from sage.combinat.fsm_fourier import FSMFourier # optional - arb
            sage: function('f')
            f
            sage: var('n')
            n
            sage: F = FSMFourier(transducers.Recursion([ # optional - arb
            ....:     f(2*n + 1) == f(n) + 1,
            ....:     f(2*n) == f(n),
            ....:     f(0) == 0],
            ....:     2, f, n))
            sage: dir = tmp_dir()
            sage: F.plot_fluctuation_asymptote( # optional - arb
            ....:     dir + "test", 8, 9, resolution=2)

        .. SEEALSO::

            :meth:`plot_fluctuation`
        """
        template = r'''
import graph;
locale("C");
size(0, 2inch, false);
unitsize(%ginch, 0);

file in = input("%s-empirical.dat").line();
real[][] a=in.dimension(0,0);
a=transpose(a);
real[] x=a[0];
real[] y_min=a[1];
real[] y_max=a[2];
int i;
for(i=0; i<x.length; ++i) {
  if (y_min[i] < y_max[i]) {
    draw((x[i],y_min[i])--(x[i],y_max[i]), blue);
  } else {
    draw((x[i],y_min[i]), blue);
  }
}

file in = input("%s-fourier.dat").line();
real[][] a=in.dimension(0,0);
a=transpose(a);
real[] x=a[0];
real[] y=a[1];
draw(graph(x,y), red+0.25bp);

xaxis(Bottom, RightTicks("$%%f$", Step=1, step=0.25));
yaxis(Left,RightTicks(trailingzero));
'''

        from sage.plot.plot import list_plot
        periods = (end-start)/self.common_period
        width_per_period = width/periods
        points_per_period = width_per_period * resolution

        empirical = _reduce_resolution_(
            self.cache.fluctuation_empirical(
                self.q ** start,
                self.q ** end),
            start,
            end,
            int(resolution * width))

        fourier = self.cache.fluctuation_fourier(start,
                                                 end,
                                                 points_per_period)

        with open("%s-empirical.dat" % prefix, "w") as f:
            for x, y, z in empirical:
                f.write("%g %g %g\n" % (x, y, z))
        with open("%s-fourier.dat" % prefix, "w") as f:
            for x, y in fourier:
                f.write("%g %g\n" % (x, y))
        with open("%s.asy" % prefix, "w") as f:
                f.write(template % ((<double> width)/(end-start),
                                    prefix, prefix))
