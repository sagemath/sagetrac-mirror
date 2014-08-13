"""
Fourier Coefficients

::

    sage: import sys
    sage: sys.path.append(".")

"""

from sage.combinat.finite_state_machine import Transducer
from sage.misc.cachefunc import cached_method
class FSMFourier(Transducer):
    """
    Fourier coefficients for the sum of output of transducers.
    """

    @cached_method
    def _fourier_coefficient_data_(self):
        """
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
            FourierCoefficientData(c=1, periods=[1], period=1, T=[1])

        -   Artificial example with two final components of periods `2`
            and `3`, respectively::

            sage: T = FSMFourier([(0, 1, 0, 1), (1, 2, 0, 1),
            ....:     (2, 1, 0, 2), (1, 2, 1, 0), (2, 1, 1, 2),
            ....:     (0, -1, 1, 1), (-1, -2, 1, 1), (-2, -3, 1, 1),
            ....:     (-3, -1, 1, 1), (-1, -2, 0, 2), (-2, -3, 0, 1),
            ....:     (-3, -1, 0, 3)],
            ....:     initial_states=[0],
            ....:     final_states=[0, 1, 2, -3, -2, -1])
            sage: T._fourier_coefficient_data_()
            (
                          [        1/7           1           1         1/5           1           1]
                          [          0           0           0         2/5          -2           0]
                          [          0           0           0         2/5           2           0]
                          [        2/7           2           2           0           0           0]
                          [        2/7    -2*zeta6 2*zeta6 - 2           0           0           0]
            2, [3, 2], 6, [        2/7 2*zeta6 - 2    -2*zeta6           0           0           0]
            )
        """
        import collections
        import operator

        from sage.matrix.constructor import matrix
        from sage.modules.free_module import VectorSpace
        from sage.rings.arith import lcm
        from sage.rings.integer_ring import ZZ
        from sage.rings.number_field.number_field import CyclotomicField
        from sage.structure.sage_object import SageObject

        FourierCoefficientData = collections.namedtuple(
            "FourierCoefficientData",
            ["c", "periods", "period", "T"])

        positions = dict((state.label(), j)
                         for j, state in enumerate(self.iter_states()))
        q = len(self.input_alphabet)

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
                    eigenvalue = q * alpha**(j * common_period / self.period)
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


        components = [FCComponent(c) for c in self.final_components()]
        period = lcm([c.period for c in components])
        field = CyclotomicField(period)
        alpha = field.gen()
        M = self.adjacency_matrix(entry=lambda t: 1)
        standard_basis = VectorSpace(field, M.nrows()).basis()

        right_eigenvectors = reduce(
            operator.add,
            [c.eigenvectors(M, components, period)
             for c in components],
            [])

        left_eigenvectors = reduce(
            operator.add,
            [c.eigenvectors(M.transpose(), components, period)
             for c in components],
            [])

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

        return FourierCoefficientData(
            c=len(components),
            periods=[c.period for c in components],
            period=period,
            T=T)
