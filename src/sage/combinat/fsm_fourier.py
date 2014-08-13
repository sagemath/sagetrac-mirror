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
            FourierCoefficientData(c=1, periods=[1], period=1)

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
            FourierCoefficientData(c=2, periods=[3, 2], period=6)
        """
        from sage.rings.arith import lcm
        import collections

        FourierCoefficientData = collections.namedtuple(
            "FourierCoefficientData",
            ["c", "periods", "period"])
        components = self.final_components()
        periods = [c.graph().period()
                   for c in components]
        period = lcm(periods)
        return FourierCoefficientData(
            c=len(components),
            periods=periods,
            period=period)
