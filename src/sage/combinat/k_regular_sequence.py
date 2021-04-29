r"""
`k`-regular Sequences

An introduction and formal definition of `k`-regular sequences can be
found, for example, on the :wikipedia:`k-regular_sequence` or in
[AS2003]_.


.. WARNING::

    As this code is experimental, warnings are thrown when a
    `k`-regular sequence space is created for the first time in a
    session (see :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        doctest:...: FutureWarning: This class/method/function is
        marked as experimental. It, its functionality or its interface
        might change without a formal deprecation.
        See http://trac.sagemath.org/21202 for details.


Examples
========

Binary sum of digits
--------------------

::

    sage: Seq2 = kRegularSequenceSpace(2, ZZ)
    sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
    ....:          left=vector([0, 1]), right=vector([1, 0]))
    sage: S
    2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
    sage: all(S[n] == sum(n.digits(2)) for n in srange(10))
    True

Number of odd entries in Pascal's triangle
------------------------------------------

::

    sage: @cached_function
    ....: def u(n):
    ....:     if n <= 1:
    ....:         return n
    ....:     return 2*u(floor(n/2)) + u(ceil(n/2))
    sage: tuple(u(n) for n in srange(10))
    (0, 1, 3, 5, 9, 11, 15, 19, 27, 29)

    sage: U = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
    ....:          left=vector([0, 1]), right=vector([1, 0]), transpose=True)
    sage: all(U[n] == u(n) for n in srange(30))
    True


Various
=======

.. SEEALSO::

    :mod:`recognizable series <sage.combinat.recognizable_series>`,
    :mod:`sage.rings.cfinite_sequence`,
    :mod:`sage.combinat.binary_recurrence_sequences`.

REFERENCES:

.. [AS2003] Jean-Paul Allouche, Jeffrey Shallit,
   *Automatic Sequences: Theory, Applications, Generalizations*,
   Cambridge University Press, 2003.
.. [HKL2021] Clemens Heuberger, Daniel Krenn, Gabriel Lipnik, *Asymptotic
   Analysis of `q`-Recursive Sequences*, 2021, unpublished manuscript.
.. [LRS2017] Julien Leroy, Michel Rigo, Manon Stipulanti, *Counting the
   number of non-zero coefficients in rows of generalized Pascal triangles*,
   Discrete Math. 340 (2017), no. 5, 862--881.

AUTHORS:

- Daniel Krenn (2016)
- Gabriel Lipnik (2021)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.
- Gabriel Lipnik is supported by the
  Austrian Science Fund (FWF): W 1230.


Classes and Methods
===================
"""
#*****************************************************************************
#       Copyright (C) 2016 Daniel Krenn <dev@danielkrenn.at>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import absolute_import

from .recognizable_series import RecognizableSeries
from .recognizable_series import RecognizableSeriesSpace
from sage.misc.cachefunc import cached_function, cached_method


class kRegularSequence(RecognizableSeries):

    def __init__(self, parent, mu, left=None, right=None):
        r"""
        A `k`-regular sequence.

        INPUT:

        - ``parent`` -- an instance of :class:`kRegularSequenceSpace`.

        - ``mu`` -- a family of square matrices, all of which have the
          same dimension. The indices of this family are `0,...,k-1`.
          ``mu`` may be a list or tuple of cardinality `k`
          as well. See
          :meth:`~sage.combinat.recognizable_series.RecognizableSeries.mu`
          for more details.

        - ``left`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        - ``right`` -- (default: ``None``) a vector.
          When evaluating the sequence, this vector is multiplied
          from the right to the matrix product. If ``None``, then this
          multiplication is skipped.

        When created via the parent :class:`kRegularSequenceSpace`, then
        the following option is available.

        - ``transpose`` -- (default: ``False``) a boolean. If set, then
          each of the matrices in
          :meth:`mu <sage.combinat.recognizable_series.RecognizableSeries.mu>`
          is transposed. Additionally the vectors
          :meth:`left <sage.combinat.recognizable_series.RecognizableSeries.left>`
          and
          :meth:`right <sage.combinat.recognizable_series.RecognizableSeries.right>`
          are switched.
          (This is done by calling :meth:`~sage.combinat.recognizable_series.RecognizableSeries.transposed`.)

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:      vector([0, 1]), vector([1, 0]),
            ....:      transpose=True)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequenceSpace`.
        """
        super(kRegularSequence, self).__init__(
            parent=parent, mu=mu, left=left, right=right)


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence.

        OUTPUT:

        A string.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: s = Seq2((Matrix([[3, 6], [0, 1]]), Matrix([[0, -6], [1, 5]])),
            ....:           vector([0, 1]), vector([1, 0]), transpose=True)
            sage: repr(s)  # indirect doctest
            '2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...'
        """
        from sage.misc.lazy_list import lazy_list_formatter
        return lazy_list_formatter(
            self,
            name='{}-regular sequence'.format(self.parent().k),
            opening_delimiter='', closing_delimiter='',
            preview=10)


    @cached_method
    def __getitem__(self, n, **kwds):
        r"""
        Return the `n`th entry of this sequence.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        An element of the universe of the sequence.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: S[7]
            3

        TESTS::

            sage: S[-1]
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned char

        ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: W = Seq2.indices()
            sage: M0 = Matrix([[1, 0], [0, 1]])
            sage: M1 = Matrix([[0, -1], [1, 2]])
            sage: S = Seq2((M0, M1), [0, 1], [1, 1])
            sage: S._mu_of_word_(W(0.digits(2))) == M0
            True
            sage: S._mu_of_word_(W(1.digits(2))) == M1
            True
            sage: S._mu_of_word_(W(3.digits(2))) == M1^2
            True
        """
        return self.coefficient_of_word(self.parent()._n_to_index_(n), **kwds)


    def __iter__(self):
        r"""
        Return an iterator.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([[1, 0], [0, 1]]), Matrix([[0, -1], [1, 2]])),
            ....:          left=vector([0, 1]), right=vector([1, 0]))
            sage: from itertools import islice
            sage: tuple(islice(S, 10))
             (0, 1, 1, 2, 1, 2, 2, 3, 1, 2)

        TESTS::

            sage: it = iter(S)
            sage: iter(it) is it
            True
            sage: iter(S) is not it
            True
        """
        from itertools import count
        return iter(self[n] for n in count())


class kRegularSequenceSpace(RecognizableSeriesSpace):
    r"""
    The space of `k`-regular Sequences over the given ``coefficients``.

    INPUT:

    - ``k`` -- an integer at least `2` specifying the base.

    - ``coefficients`` -- a (semi-)ring. If not specified (``None``),
      then the integer ring is used.

    - ``category`` -- (default: ``None``) the category of this
      space.

    EXAMPLES::

        sage: kRegularSequenceSpace(2, ZZ)
        Space of 2-regular sequences over Integer Ring
        sage: kRegularSequenceSpace(3, ZZ)
        Space of 3-regular sequences over Integer Ring

    .. SEEALSO::

        :doc:`k-regular sequence <k_regular_sequence>`,
        :class:`kRegularSequence`.
    """

    Element = kRegularSequence


    @classmethod
    def __normalize__(cls, k, coefficients=None, **kwds):
        r"""
        Normalizes the input in order to ensure a unique
        representation.

        For more information see :class:`kRegularSequenceSpace`.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2.category()
            Category of sets

        ::

            sage: Seq2 is kRegularSequenceSpace(2)
            True
        """
        from sage.arith.srange import srange
        from sage.rings.integer_ring import ZZ
        if coefficients is None:
            coefficients = ZZ
        nargs = super(kRegularSequenceSpace, cls).__normalize__(
            coefficients, alphabet=srange(k), **kwds)
        return (k,) + nargs


    def __init__(self, k, *args, **kwds):
        r"""
        See :class:`kRegularSequenceSpace` for details.

        INPUT:

        - ``k`` -- an integer at least `2` specifying the base.

        Other input arguments are passed on to
        :meth:`~sage.combinat.recognizable_series.RecognizableSeriesSpace.__init__`.

        TESTS::

            sage: kRegularSequenceSpace(2, ZZ)
            Space of 2-regular sequences over Integer Ring
            sage: kRegularSequenceSpace(3, ZZ)
            Space of 3-regular sequences over Integer Ring

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self.k = k
        super(kRegularSequenceSpace, self).__init__(*args, **kwds)


    def _repr_(self):
        r"""
        Return a representation string of this `k`-regular sequence space.

        OUTPUT:

        A string.

        TESTS::

            sage: repr(kRegularSequenceSpace(2, ZZ))  # indirect doctest
            'Space of 2-regular sequences over Integer Ring'
        """
        return 'Space of {}-regular sequences over {}'.format(self.k, self.base())


    def _n_to_index_(self, n):
        r"""
        Convert `n` to an index usable by the underlying
        recognizable series.

        INPUT:

        - ``n`` -- a nonnegative integer.

        OUTPUT:

        A word.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._n_to_index_(6)
            word: 011
            sage: Seq2._n_to_index_(-1)
            Traceback (most recent call last):
            ...
            OverflowError: can't convert negative value to unsigned char
        """
        from sage.rings.integer_ring import ZZ
        n = ZZ(n)
        W = self.indices()
        return W(n.digits(self.k))


    def _parse_recurrence_(self, equations, function, var):
        r"""
        Parse recurrence relations as admissible in :meth:`from_recurrence`.

        INPUT:

        - ``equations`` -- see :meth:`from_recurrence`

        - ``function`` -- see :meth:`from_recurrence`

        - ``var`` -- see :meth:`from_recurrence`

        OUTPUT:

        A vector consisting of

        - ``M``, ``m`` -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1

        - ``coeffs`` -- a dictionary mapping ``(r, j)`` to the coefficients
          `c_{r, j}` as given in [HKL2021]_, Equation (3.1)

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``th value of the sequence

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2._parse_recurrence_([
            ....:     f(4*n) == f(2*n) + 2*f(2*n + 1) + 3*f(2*n - 2),
            ....:     f(4*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n - 2),
            ....:     f(4*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n - 2),
            ....:     f(4*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n - 2),
            ....:     f(0) == 1, f(1) == 2, f(2) == 1], f, n)
            (2, 1, {(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12, (3, 0): 10,
            (3, 1): 11}, {0: 1, 1: 2, 2: 1})

        Stern--Brocot Sequence::

            sage: Seq2._parse_recurrence_([
            ....:    f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:    f(0) == 0, f(1) == 1], f, n)
            (1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 1})

        .. SEEALSO::

            :meth:`from_recurrence`

        TESTS:

        The following tests check that the equations are well-formed::

            sage: Seq2._parse_recurrence_([], f, n)
            Traceback (most recent call last):
            ...
            ValueError: List of recurrence equations is empty.

        ::

            sage: Seq2._parse_recurrence_([f(4*n + 1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(4*n + 1) is not an equation with ==.

        ::

            sage: Seq2._parse_recurrence_([42], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 42 is not a symbolic expression.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) + 1 == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(2*n) + 1 is not an evaluation of f.

        ::

            sage: Seq2._parse_recurrence_([f(2*n, 5) == 3], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(2*n, 5) does not have one argument.

        ::

            sage: Seq2._parse_recurrence_([f(1/n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ....:
            ValueError: 1/n + 1 is not a polynomial in n with integer coefficients.

        ::

            sage: Seq2._parse_recurrence_([f(2*n + 1/2) == f(n)], f, n)
            Traceback (most recent call last):
            ....:
            ValueError: 2*n + 1/2 is not a polynomial in n with integer coefficients.

        ::

            sage: Seq2._parse_recurrence_([f(4*n^2) == f(2*n^2)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 4*n^2 is not a polynomial of degree smaller 2.

        ::

            sage: Seq2._parse_recurrence_([f(42) == 1/2], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not in Integer Ring.

        ::

            sage: Seq2._parse_recurrence_([f(42) == 0, f(42) == 1], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value f(42) is given twice.

        ::

            sage: Seq2._parse_recurrence_([f(42) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(n) is not in Integer Ring.

        ::

            sage: Seq2._parse_recurrence_([f(4*n) == f(n), f(2*n) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2 does not equal 4.

        ::

            sage: Seq2._parse_recurrence_([f(3*n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a power of 2.

        ::

            sage: Seq2._parse_recurrence_([f(n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1 is less than 2.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(n), f(2*n) == 0], f, n)
            Traceback (most recent call last):
            ...
            ValueError: There are more than one recurrence relation for f(2*n).

        ::

            sage: Seq2._parse_recurrence_([f(2*n + 2) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2 is not smaller than 2.

        ::

            sage: Seq2._parse_recurrence_([f(2*n - 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: -1 is smaller than 0.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 2*n], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2*n does not contain f.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 1/2*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/2 is not a valid coefficient.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 1/f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/f(n) is not a valid right hand side.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 2*n*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2*n*f(n) is not a valid right hand side.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 2*f(n, 5)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f(n, 5) has more than one argument.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 2*f()], f, n)
            Traceback (most recent call last):
            ...
            ValueError: f() has no argument.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 1/f(n) + 2*f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/f(n) is not a valid summand.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == 2*f(1/n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/n is not a polynomial with integer coefficients.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(n + 1/2)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: n + 1/2 is not a polynomial with integer coefficients.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(1/2*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1/2*n is not a polynomial with integer coefficients.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(n^2 + 1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: n^2 + 1 does not have degree 1.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(1)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1 does not have degree 1.

        ::

            sage: Seq2._parse_recurrence_([f(4*n) == f(2*n) + f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1 does not equal 2.

        ::

            sage: Seq2._parse_recurrence_([f(4*n) == f(2*n), f(4*n + 1) == f(n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 1 does not equal 2.

        ::

            sage: Seq2._parse_recurrence_([f(4*n) == f(3*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 3 is not a power of 2.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(4*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 4 is not smaller than 2.

        ::

            sage: Seq2._parse_recurrence_([f(2*n) == f(2*n)], f, n)
            Traceback (most recent call last):
            ...
            ValueError: 2 is not smaller than 2.

        ::

            sage: Seq2._parse_recurrence_([f(42) == 0], f, n)
            Traceback (most recent call last):
            ...
            ValueError: No recurrence relations are given.

        Finally, also for the zero-sequence the output is as expected::

            sage: Seq2._parse_recurrence_([f(2*n) == 0, f(2*n + 1) == 0], f, n)
            (1, 0, {}, {})
        """
        from sage.arith.srange import srange
        from sage.functions.log import log
        from sage.rings.integer_ring import ZZ
        from sage.symbolic.operators import add_vararg, mul_vararg, operator

        k = self.k
        base_ring = self.base()
        M = None
        m = None
        coeffs = {}
        initial_values = {}
        remainders = []

        def _parse_multiplication_(op):
            operands = op.operands()
            assert op.operator() == mul_vararg and len(operands) == 2
            if operands[1].operator() == function:
                return [operands[0], operands[1]]
            elif operands[0].operator() == function:
                return [operands[1], operands[0]]
            else:
                raise ValueError('%s does not contain %s.' % (op, function))

        def _parse_one_summand_(summand):
            if summand.operator() == mul_vararg:
                coeff, op = _parse_multiplication_(summand)
            elif summand.operator() == function:
                coeff, op = 1, summand
            else:
                raise ValueError('%s is not a valid summand.' % (summand,))
            if len(op.operands()) > 1:
                raise ValueError('%s has more than one argument.' % (op,))
            elif len(op.operands()) == 0:
                raise ValueError('%s has no argument.' % (op,))
            try:
                poly = ZZ[var](op.operands()[0])
            except TypeError:
                raise ValueError('%s is not a polynomial with integer coefficients.'
                                 % (op.operands()[0],))
            if poly.degree() != 1:
                raise ValueError("%s does not have degree 1."
                                 % (poly,))
            d, base_power_m = list(poly)
            m = log(base_power_m, base=k)
            return [coeff, m, d]

        if not equations:
            raise ValueError("List of recurrence equations is empty.")

        for eq in equations:
            try:
                if eq.operator() != operator.eq:
                    raise ValueError("%s is not an equation with ==."  % eq)
            except AttributeError:
                raise ValueError("%s is not a symbolic expression."  % eq)
            left_side, right_side = eq.operands()
            if left_side.operator() != function:
                raise ValueError("%s is not an evaluation of %s."
                                 % (left_side, function))
            if  len(left_side.operands()) != 1:
                raise ValueError("%s does not have one argument." %
                                 (left_side,))
            try:
                polynomial_left = ZZ[var](left_side.operands()[0])
            except TypeError:
                raise ValueError("%s is not a polynomial in %s with "
                                 "integer coefficients."
                                 % (left_side.operands()[0], var))
            if polynomial_left.degree()  > 1:
                raise ValueError("%s is not a polynomial of degree smaller 2."
                                 % (polynomial_left,))
            if polynomial_left in ZZ:
                if right_side in base_ring:
                    if (polynomial_left in initial_values.keys() and
                        initial_values[polynomial_left] != right_side):
                        raise ValueError("Initial value %s is given twice."
                                         % (function(polynomial_left)))
                    initial_values.update({polynomial_left: right_side})
                else:
                    raise ValueError("%s is not in %s." % (right_side, base_ring))
            else:
                [r, base_power_M] = list(polynomial_left)
                M_new = log(base_power_M, base=k)
                if M and M != M_new:
                    raise ValueError("%s does not equal %s."
                                     % (base_power_M, k**M))
                elif not M:
                    M = M_new
                    if M not in ZZ:
                        raise ValueError("%s is not a power of %s."
                                         % (base_power_M, k))
                    if M < 1:
                        raise ValueError("%s is less than %s."
                                         % (base_power_M, k))
                if r in remainders:
                    raise ValueError("There are more than one recurrence relation for %s."
                                     % (left_side,))
                if r >= k**M:
                    raise ValueError("%s is not smaller than %s." % (r, k**M))
                elif r < 0:
                    raise ValueError("%s is smaller than 0." % (r,))
                else:
                    remainders.append(r)
                if right_side != 0:
                    if (len(right_side.operands()) == 1 and right_side.operator() == function
                        or right_side.operator() == mul_vararg and len(right_side.operands()) == 2):
                        summands = [right_side]
                    elif right_side.operator() == add_vararg:
                        summands = right_side.operands()
                    else:
                        raise ValueError("%s is not a valid right hand side."
                                         % (right_side,))
                    for summand in summands:
                        coeff, new_m, d = _parse_one_summand_(summand)
                        if coeff not in base_ring:
                            raise ValueError("%s is not a valid coefficient."
                                             % (coeff,))
                        if m and m != new_m:
                            raise ValueError("%s does not equal %s."
                                             % (k**new_m, k**m))
                        elif not m:
                            m = new_m
                            if m not in ZZ:
                                raise ValueError("%s is not a power of %s."
                                                 % (k**m, k))
                            if M <= m:
                                raise ValueError("%s is not smaller than %s."
                                                 % (k**m, k**M))
                        coeffs.update({(r, d): coeff})

        if not M:
            raise ValueError("No recurrence relations are given.")
        elif M and not m: # for the zero sequence
            m = M - 1

        return (M, m, coeffs, initial_values)


    def _get_parameters_from_recurrence_(self, M, m, coeffs, initial_values,
                                         offset):
        r"""
        Determine parameters from recurrence relations as admissible in
        :meth:`from_recurrence`.

        INPUT:

        - ``M``, ``m``, ``offset`` -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1 (see also :meth:`from_recurrence`)

        - ``coeffs`` -- a dictionary where ``coeffs[(r, j)]`` is the
          coefficient `c_{r,j}` as given in :meth:`from_recurrence`.
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``th value of the sequence

        OUTPUT:

        A namedtuple ``recurrence_rules`` consisting of

        - ``M``, ``m``, ``l``, ``u``  -- parameters of the recursive sequences,
          see [HKL2021]_, Definition 3.1

        - ``ll``, ``uu``, ``n1``, ``dim`` -- parameters and dimension of the
          resulting linear representation, see [HKL2021]_, Theorem A

        - ``coeffs`` -- a dictionary mapping ``(r, j)`` to the coefficients
          `c_{r, j}` as given in [HKL2021]_, Equation (3.1).
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``th value of the sequence

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._get_parameters_from_recurrence_(2, 1,
            ....: {(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            ....: (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12,
            ....: (3, 0): 10, (3, 1): 11}, {0: 1, 1: 2, 2: 1, 3: 4}, 0)
            recurrence_rules(M=2, m=1, l=-2, u=1, ll=-6, uu=3, dim=14,
            coeffs={(0, -2): 3, (0, 0): 1, (0, 1): 2, (1, -2): 6, (1, 0): 4,
            (1, 1): 5, (2, -2): 9, (2, 0): 7, (2, 1): 8, (3, -2): 12,
            (3, 0): 10, (3, 1): 11}, initial_values={0: 1, 1: 2, 2: 1, 3: 4,
            4: 12, 5: 30, 6: 48, 7: 66, 8: 75, 9: 204, 10: 333, 11: 462,
            12: 216, 13: 594, -6: 0, -5: 0, -4: 0, -3: 0, -2: 0, -1: 0},
            offset=1, n1=3)

        .. SEEALSO::

            :meth:`from_recurrence`

        TESTS::

            sage: Seq2._get_parameters_from_recurrence_(1, 0, {(0, 0): 1}, {}, 0)
            Traceback (most recent call last):
            ...
            ValueError: No initial values are given.

        ::

            sage: Seq2._get_parameters_from_recurrence_(1, 0, {(0, 0): 1},
            ....: {0: 1, 1: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={(0, 0): 1}, initial_values={0: 1, 1: 0}, offset=0, n1=0)

        Finally, also for the zero-sequence the output is as expected::

            sage: Seq2._get_parameters_from_recurrence_(1, 0, {}, {0: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={}, initial_values={0: 0}, offset=0, n1=0)

        ::

            sage: Seq2._get_parameters_from_recurrence_(1, 0,
            ....: {(0, 0): 0, (1, 1): 0}, {0: 0}, 0)
            recurrence_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
            coeffs={(0, 0): 0, (1, 1): 0}, initial_values={0: 0},
            offset=0, n1=0)
        """
        from collections import namedtuple

        from sage.arith.srange import srange
        from sage.functions.other import ceil, floor

        k = self.k
        keys = coeffs.keys()
        indices_right = [key[1] for key in keys if coeffs[key]]

        if not indices_right: # the sequence is the zero sequence
            l = 0
            u = 0
        else:
            l = min(indices_right)
            u = max(indices_right)

        if offset < max(0, -l/k**m):
            offset = max(0, -l/k**m)

        ll = (floor((l*k**(M-m) - k**M + 1)/(k**(M-m) - 1)) + 1)*(l < 0)
        uu = max([ceil((u*k**(M-m) + k**M - k**m)/(k**(M-m) - 1)) - 1, k**m - 1])
        n1 = offset - floor(ll/k**M)
        dim = (k**M - 1)/(k - 1) + (M - m)*(uu - ll - k**m + 1) + n1

        if not initial_values:
            raise ValueError("No initial values are given.")

        last_value_needed = max(
            k**(M-1) - k**m + uu + (n1 > 0)*k**(M-1)*(k*(n1 - 1) + k - 1), # for matrix W
            k**m*offset + u,
            max(initial_values.keys()))
        initial_values = self._get_values_from_recurrence_(
            M, m, l, u, ll, coeffs, initial_values, last_value_needed, offset)

        recurrence_rules = namedtuple('recurrence_rules',
                                      ['M', 'm', 'l', 'u', 'll', 'uu', 'dim',
                                       'coeffs', 'initial_values', 'offset', 'n1'])

        return recurrence_rules(M=M, m=m, l=l, u=u, ll=ll, uu=uu, dim=dim,
                                coeffs=coeffs, initial_values=initial_values,
                                offset=offset, n1=n1)


    def _get_values_from_recurrence_(self, M, m, l, u, ll, coeffs, initial_values,
                                     last_value_needed, offset):
        r"""
        Determine enough values of the corresponding recursive sequence by
        applying the recurrence relations given in :meth:`from_recurrence`
        to the values given in ``initial_values``.

        INPUT:

        - ``M``, ``m``, ``l``, ``u`` and ``offset`` -- parameters of the
          recursive sequences, see [HKL2021]_, Definition 3.1

        - ``ll`` -- parameter of the resulting linear representation,
          see [HKL2021]_, Theorem A

        - ``coeffs`` -- a dictionary where ``coeffs[(r, j)]`` is the
          coefficient `c_{r,j}` as given in :meth:`from_recurrence`.
          If ``coeffs[(r, j)]`` is not given for some ``r`` and ``j``,
          then it is assumed to be zero.

        - ``initial_values`` -- a dictionary mapping integers ``n`` to the
          ``n``th value of the sequence

        - ``last_value_needed`` -- last initial value which is needed to
          determine the linear representation

        OUTPUT:

        A dictionary mapping integers ``n`` to the ``n``th value of the
        sequence for all ``n`` up to ``last_value_needed``.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 1, 2: 1},
            ....: 20, 0)
            {0: 0, 1: 1, 2: 1, 3: 2, 4: 1, 5: 3, 6: 2, 7: 3, 8: 1, 9: 4, 10: 3,
            11: 5, 12: 2, 13: 5, 14: 3, 15: 4, 16: 1, 17: 5, 18: 4, 19: 7, 20: 3}

        .. SEEALSO::

            :meth:`from_recurrence`

        TESTS:

        For the equations `f(2n) = f(n)` and `f(2n + 1) = f(n) + f(n + 1)`::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 2}, 20, 0)
            {0: 0, 1: 2, 2: 2, 3: 4, 4: 2, 5: 6, 6: 4, 7: 6, 8: 2, 9: 8, 10: 6,
            11: 10, 12: 4, 13: 10, 14: 6, 15: 8, 16: 2, 17: 10, 18: 8, 19: 14,
            20: 6}

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 2*i}, 20, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial value for n = 1 is not in Integer Ring.

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {}, 20, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for n in [0, 1] are missing.

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0}, 20, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for n in [1] are missing.

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 2: 1}, 20, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for n in [1] are missing.

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 1, 0,
            ....: {(0, 0): 1, (1, 0): 1, (1, 1): 1}, {0: 0, 1: 2, 2:0}, 20, 0)
            Traceback (most recent call last):
            ...
            ValueError: Initial value for n = 2 does not match with the given
            recurrence relations.

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, -2, 2, -2,
            ....: {(0, -2): 1, (0, 2): 1, (1, -2): 1, (1, 2): 1},
            ....: {0: 0, 1: 2, 2: 4, 3: 3, 4: 2}, 20, 2)
            {-2: 0, -1: 0, 0: 0, 1: 2, 2: 4, 3: 3, 4: 2, 5: 2, 6: 4, 7: 4,
            8: 8, 9: 8, 10: 7, 11: 7, 12: 10, 13: 10, 14: 10, 15: 10, 16: 11,
            17: 11, 18: 11, 19: 11, 20: 18}

        Finally, also for the zero-sequence the output is as expected::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 0, 0, {}, {}, 10, 0)
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}

        ::

            sage: Seq2._get_values_from_recurrence_(1, 0, 0, 0, 0,
            ....: {(0, 0): 0, (1, 1): 0}, {}, 10, 0)
            {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0, 8: 0, 9: 0, 10: 0}
        """
        from sage.arith.srange import srange
        from sage.functions.other import ceil
        from sage.rings.integer_ring import ZZ

        k = self.k
        base_ring = self.base_ring()
        keys_initial = initial_values.keys()

        values = {n: None if n not in keys_initial else initial_values[n]
                  for n in srange(last_value_needed + 1)}
        missing_values = []

        @cached_function
        def _coeff_(r, k):
            try:
                return coeffs[(r, k)]
            except KeyError:
                return 0

        def _f_n_(n):
            f_n = values[n]
            if f_n is not None and f_n != "pending":
                if f_n not in base_ring:
                    raise ValueError("Initial value for n = %s is not in %s."
                                     % (n, base_ring))
                return f_n
            elif f_n == "pending":
                missing_values.append(n)
                return 0
            else:
                values.update({n: "pending"})
                q, r = ZZ(n).quo_rem(k**M)
                if q < offset:
                    missing_values.append(n)
                return sum([_coeff_(r, j)*_f_n_(k**m*q + j)
                            for j in srange(l, u + 1)
                            if _coeff_(r, j) != 0])

        for n in srange(last_value_needed + 1):
            values.update({n: _f_n_(n)})

        if missing_values:
            raise ValueError("Initial values for n in %s are missing."
                             % (missing_values,))

        for n in keys_initial:
            q, r = ZZ(n).quo_rem(k**M)
            if (q >= offset and
                values[n] != sum([_coeff_(r, j)*values[k**m*q + j]
                                  for j in srange(l, u + 1)])):
                raise ValueError("Initial value for n = %s does not match with "
                                 "the given recurrence relations." % (n,))

        values.update({n: 0 for n in srange(ll, 0)})

        return values

    @cached_method
    def _get_ind_from_recurrence_(self, M, m, ll, uu):
        r"""
        Determine the index operator corresponding to the recursive
        sequence given by ``recurrence_rules``, as defined in [HKL2021]_.

        INPUT:

        - ``recurrence_rules`` -- A namedtuple generated by
          :meth:`_parse_recurrence_`.

        OUTPUT:

        A dictionary which maps both row numbers to subsequence parameters and
        vice versa, i.e.,

        - ``ind[i]`` -- a pair ``(j, d)`` representing the sequence `x(k^j + d)`
          in the `i`th component (1-based) of the resulting linear representation

        - ``ind[(j, d)]`` -- the (1-based) row number of the sequence
          `x(k^j + d)` in the linear representation.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._get_ind_from_recurrence_(3, 1, -3, 3)
            {(0, 0): 1, (1, -1): 4, (1, -2): 3, (1, -3): 2,
            (1, 0): 5, (1, 1): 6, (1, 2): 7, (1, 3): 8, (2, -1): 11,
            (2, -2): 10, (2, -3): 9, (2, 0): 12, (2, 1): 13, (2, 2): 14,
            (2, 3): 15, (2, 4): 16, (2, 5): 17, 1: (0, 0), 10: (2, -2),
            11: (2, -1), 12: (2, 0), 13: (2, 1), 14: (2, 2), 15: (2, 3),
            16: (2, 4), 17: (2, 5), 2: (1, -3), 3: (1, -2), 4: (1, -1),
            5: (1, 0), 6: (1, 1), 7: (1, 2), 8: (1, 3), 9: (2, -3)}

        .. SEEALSO::

            :meth:`from_recurrence`
        """
        from sage.arith.srange import srange

        k = self.k
        ind = {}

        pos = 1
        for j in srange(m):
            for d in srange(k**j):
                ind.update({(j, d): pos, pos: (j, d)})
                pos += 1
        for j in srange(m, M):
            for d in srange(ll, k**j - k**m + uu + 1):
                ind.update({(j, d): pos, pos: (j, d)})
                pos += 1

        return ind

    def _v_eval_n_from_recurrence_(self, recurrence_rules, n):
        r"""
        Return the vector `v(n)` as given in [HKL2021]_, Theorem A.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`_parse_recurrence_`

        - ``n`` -- an integer

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: SB_rules = Seq2._get_parameters_from_recurrence_(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: Seq2._v_eval_n_from_recurrence_(SB_rules, 0)
            (0, 1, 1)

        .. SEEALSO::

            :meth:`from_recurrence`
        """
        from sage.arith.srange import srange
        from sage.modules.free_module_element import vector

        k = self.k
        M = recurrence_rules.M
        m = recurrence_rules.m
        ll = recurrence_rules.ll
        uu = recurrence_rules.uu
        initial_values = recurrence_rules.initial_values

        return vector(
            [initial_values[k**j*n + d] for j in srange(m)
             for d in srange(k**j)] + \
            [initial_values[k**j*n + d] for j in srange(m, M)
             for d in srange(ll, k**j - k**m + uu + 1)])

    def _get_matrix_from_recurrence_(self, recurrence_rules, rem,
                                     correct_offset=True):
        r"""
        Construct the matrix for remainder ``rem`` of the linear
        representation of the sequence represented by ``recurrence_rules``.

        INPUT:

        - ``recurrence_rules`` -- a namedtuple generated by
          :meth:`_parse_recurrence_`

        - ``rem`` -- an integer between ``0`` and ``k - 1``

        - ``function`` -- a function which represents the sequence

        - ``var`` -- a symbolic variable ``n``

        - ``correct_offset`` -- (default: ``True``) a boolean. If
          ``True``, then the resulting linear representation has no
          offset.  See [HKL2021]_ for more information.

        OUTPUT:

        A matrix.

        EXAMPLES:

        The following example illustrates how the coefficients in the
        right-hand sides of the recurrence relations correspond to the entries of
        the matrices. ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: M, m, coeffs, initial_values = Seq2._parse_recurrence_([
            ....:     f(8*n) == -1*f(2*n - 1) + 1*f(2*n + 1),
            ....:     f(8*n + 1) == -11*f(2*n - 1) + 10*f(2*n) + 11*f(2*n + 1),
            ....:     f(8*n + 2) == -21*f(2*n - 1) + 20*f(2*n) + 21*f(2*n + 1),
            ....:     f(8*n + 3) == -31*f(2*n - 1) + 30*f(2*n) + 31*f(2*n + 1),
            ....:     f(8*n + 4) == -41*f(2*n - 1) + 40*f(2*n) + 41*f(2*n + 1),
            ....:     f(8*n + 5) == -51*f(2*n - 1) + 50*f(2*n) + 51*f(2*n + 1),
            ....:     f(8*n + 6) == -61*f(2*n - 1) + 60*f(2*n) + 61*f(2*n + 1),
            ....:     f(8*n + 7) == -71*f(2*n - 1) + 70*f(2*n) + 71*f(2*n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7], f, n)
            sage: rules = Seq2._get_parameters_from_recurrence_(
            ....:     M, m, coeffs, initial_values, 0)
            sage: Seq2._get_matrix_from_recurrence_(rules, 0, False)
            [  0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            [  0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0 -61  60  61   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0 -71  70  71   0   0   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -21  20  21   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -31  30  31   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -41  40  41   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0]
            sage: Seq2._get_matrix_from_recurrence_(rules, 1, False)
            [  0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1   0]
            [  0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]
            [  0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -21  20  21   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -31  30  31   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -41  40  41   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -51  50  51   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -61  60  61   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0 -71  70  71   0   0   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0  -1   0   1   0   0   0   0   0   0   0   0   0]
            [  0   0   0   0   0 -11  10  11   0   0   0   0   0   0   0   0   0]

        Stern--Brocot Sequence::

            sage: SB_rules = Seq2._get_parameters_from_recurrence_(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: Seq2._get_matrix_from_recurrence_(SB_rules, 0)
            [1 0 0]
            [1 1 0]
            [0 1 0]
            sage: Seq2._get_matrix_from_recurrence_(SB_rules, 1)
            [1 1 0]
            [0 1 0]
            [0 1 1]

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: M, m, coeffs, initial_values = Seq2._parse_recurrence_([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n)
            sage: UB_rules = Seq2._get_parameters_from_recurrence_(
            ....:     M, m, coeffs, initial_values, 3)
            sage: Seq2._get_matrix_from_recurrence_(UB_rules, 0)
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  2  0  0  0  0  0  0  0  0  0 -1  0  0]
            [ 0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  1  0  1  0  0  0  0  0  0 -4  0  0]
            [ 0  0  0  0 -1  1  0  0  0  0  0  0  0  4  2  0]
            [ 0  0  0  0  0  2  0  0  0  0  0  0  0 -2  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0 -1  1  1  0  0  0  0  0  0  2  2  0]
            [ 0  0  0  0  2  0  1  0  0  0  0  0  0 -8 -4 -4]
            [ 0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0]
            sage: Seq2._get_matrix_from_recurrence_(UB_rules, 1)
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  2  0  0  0  0  0  0  0 -2  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0 -1  1  1  0  0  0  0  0  0  2  2  0]
            [ 0  0  0  0  2  0  1  0  0  0  0  0  0 -8 -4 -4]
            [ 0  0  0  0  0  0  0  2  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0 -1  1  0  0  0  2  0  0]
            [ 0  0  0  0  0  0  0  0  0  2  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0]

        .. SEEALSO::

            :meth:`from_recurrence`
        """
        from sage.arith.srange import srange
        from sage.matrix.constructor import Matrix
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.special import block_matrix, zero_matrix
        from sage.modules.free_module_element import vector

        base_ring = self.base_ring()
        k = self.k
        M = recurrence_rules.M
        m = recurrence_rules.m
        l = recurrence_rules.l
        ll = recurrence_rules.ll
        uu = recurrence_rules.uu
        dim = recurrence_rules.dim
        n1 = recurrence_rules.n1
        dim_without_corr = dim - n1
        coeffs = recurrence_rules.coeffs
        initial_values = recurrence_rules.initial_values
        ind = self._get_ind_from_recurrence_(M, m, ll, uu)

        mat = Matrix(base_ring, 0, dim_without_corr)

        @cached_function
        def _coeff_(r, k):
            try:
                return coeffs[(r, k)]
            except KeyError:
                return 0

        for i in srange(1, dim_without_corr + 1):
            j, d = ind[i]
            if j < M - 1:
                row = dim_without_corr*[0]
                row[ind[(j + 1, k**j*rem + d)] - 1] = 1
                mat = mat.stack(vector(row))
            else:
                rem_d = k**(M-1)*rem + (d%k**M)
                dd = d // k**M
                if rem_d < k**M:
                    lambd = l - ind[(m, (k**m)*dd + l)]
                    row = [_coeff_(rem_d, kk + lambd)
                           for kk in srange(1, dim_without_corr + 1)]
                    mat = mat.stack(vector(row))
                else:
                    lambd = l - ind[(m, k**m*dd + k**m + l)]
                    row = [_coeff_(rem_d - k**M, kk + lambd)
                           for kk in srange(1, dim_without_corr + 1)]
                    mat = mat.stack(vector(row))

        if n1 == 0 or not correct_offset:
            return mat
        else:
            W = Matrix(base_ring, dim_without_corr, 0)
            for i in srange(n1):
                W = W.augment(
                    self._v_eval_n_from_recurrence_(recurrence_rules, k*i + rem) -
                    mat*self._v_eval_n_from_recurrence_(recurrence_rules, i))

            J = Matrix(base_ring, 0, n1)
            for i in srange(n1):
                J = J.stack(vector([int(j*k == i - rem) for j in srange(n1)]))

            Mat = MatrixSpace(base_ring, dim, dim)
            return Mat(block_matrix([[mat, W],
                                     [zero_matrix(n1, dim_without_corr), J]]))


    def _get_left_from_recurrence_(self, dim):
        r"""
        Construct the vector ``left`` of the linear representation of
        recursive sequences.

        INPUT:

        - ``dim`` -- A positive integer.

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._get_left_from_recurrence_(5)
            (1, 0, 0, 0, 0)

        .. SEEALSO::

            :meth:`from_recurrence`
        """
        from sage.modules.free_module_element import vector

        return vector([1] + (dim - 1)*[0])


    def _get_right_from_recurrence_(self, recurrence_rules):
        r"""
        Construct the vector ``right`` of the linear
        representation of the sequence induced by ``recurrence_rules``.

        INPUT:

        - ``recurrence_rules`` -- A namedtuple generated by
          :meth:`_parse_recurrence_`.

        OUTPUT:

        A vector.

        TESTS:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: SB_rules = Seq2._get_parameters_from_recurrence_(
            ....:     1, 0, {(0, 0): 1, (1, 0): 1, (1, 1): 1},
            ....:     {0: 0, 1: 1, 2: 1}, 0)
            sage: Seq2._get_right_from_recurrence_(SB_rules)
            (0, 1, 1)

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: M, m, coeffs, initial_values = Seq2._parse_recurrence_([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n)
            sage: UB_rules = Seq2._get_parameters_from_recurrence_(
            ....:     M, m, coeffs, initial_values, 3)
            sage: Seq2._get_right_from_recurrence_(UB_rules)
            (1, 1, 2, 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, 1, 0, 0)

        .. SEEALSO::

            :meth:`from_recurrence`
        """
        from sage.arith.srange import srange
        from sage.modules.free_module_element import vector

        k = self.k
        M = recurrence_rules.M
        m = recurrence_rules.m
        ll = recurrence_rules.ll
        uu = recurrence_rules.uu
        initial_values = recurrence_rules.initial_values
        n1 = recurrence_rules.n1
        right = self._v_eval_n_from_recurrence_(recurrence_rules, 0)

        if n1 >= 1:
            right = vector(list(right) + [1] + (n1 - 1)*[0])

        return vector(right)


    def from_recurrence(self, equations, function, var, offset=0, minimize=False):
        r"""
        Construct a `k`-regular sequence that fulfills the recurrence relations
        given in ``equations``.

        INPUT:

        - ``equations`` -- A list of equations where the elements have
          either the form

          - `f(k^M n + r) = c_{r,l} f(k^m n + l) + c_{r,l + 1} f(k^m n
            + l + 1) + ... + c_{r,u} f(k^m n + u)` for some integers
            `0 \leq r < k^M`, `M > m \geq 0` and `l \leq u`, and some
            coefficients `c_{r,j}` from the (semi)ring ``coefficents``
            of the corresponding :class:`kRegularSequenceSpace`, valid
            for all integers `n \geq \text{offset}` for some integer
            `\text{offset} \geq \max(-l/k^m, 0)` (default: ``0``), and
            there is an equation of this form (with the same
            parameters `M` and `m`) for all `r`

          or the form

          - ``f(k) == t`` for some integer ``k`` and some ``t`` from the (semi)ring
            ``coefficients``.

          The recurrence relations above uniquely determine a `k`-regular sequence;
          see [HKL2021]_ for further information.

        - ``function`` -- symbolic function ``f`` occuring in the equations

        - ``var`` -- symbolic variable ``n`` occuring in the equations

        - ``minimize`` -- a boolean (default: ``False``). If ``True``, then
          :meth:`~sage.combinat.recognizable_series.RecognizableSeries.minimized`
          is called after the construction.

        OUTPUT: A :class:`kRegularSequence`.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2.from_recurrence([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            2-regular sequence 0, 1, 1, 2, 1, 3, 2, 3, 1, 4, ...

        Number of Odd Entries in Pascal's Triangle::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 3*f(n), f(2*n + 1) == 2*f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: Seq2.from_recurrence([
            ....:     f(8*n) == 2*f(4*n),
            ....:     f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2),
            ....:     f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14) == 4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8], f, n, 3)
            2-regular sequence 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, ...

        Number of Non-Zero Elements in the Generalized Pascal's Triangle (see [LRS2017]_)::

            sage: Seq2 = kRegularSequenceSpace(2, QQ)
            sage: Seq2.from_recurrence([
            ....:     f(4*n) == 5/3*f(2*n) - 1/3*f(2*n + 1),
            ....:     f(4*n + 1) == 4/3*f(2*n) + 1/3*f(2*n + 1),
            ....:     f(4*n + 2) == 1/3*f(2*n) + 4/3*f(2*n + 1),
            ....:     f(4*n + 3) == -1/3*f(2*n) + 5/3*f(2*n + 1),
            ....:     f(0) == 1, f(1) == 2], f, n)
            2-regular sequence 1, 2, 3, 3, 4, 5, 5, 4, 5, 7, ...

        TESTS::

            sage: Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n + 1024),
            ....:     f(0) == 1, f(1) == 1], f, n, 2)
            Traceback (most recent call last):
            ...
            ValueError: Initial values for n in [2, ..., 2044] are missing.

        ::

            sage: S = Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n + 16),
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7, f(16) == 4, f(18) == 4,
            ....:     f(20) == 4, f(22) == 4, f(24) == 6, f(26) == 6, f(28) == 6],
            ....:     f, n, 2)
            sage: all([S[4*i] == S[2*i] and
            ....:      S[4*i + 1] == S[2*i] and
            ....:      S[4*i + 2] == S[2*i] and
            ....:      S[4*i + 3] == S[2*i + 16] for i in srange(2, 100)])
            True

        ::

            sage: S = Seq2.from_recurrence([
            ....:     f(4*n) == f(2*n),
            ....:     f(4*n + 1) == f(2*n),
            ....:     f(4*n + 2) == f(2*n),
            ....:     f(4*n + 3) == f(2*n - 16),
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3, f(4) == 4,
            ....:     f(5) == 5, f(6) == 6, f(7) == 7, f(8) == 8, f(9) == 9,
            ....:     f(10) == 10, f(11) == 11, f(12) == 12, f(13) == 13,
            ....:     f(14) == 14, f(15) == 15, f(16) == 16, f(17) == 17,
            ....:     f(18) == 18, f(19) == 19, f(20) == 20, f(21) == 21,
            ....:     f(22) == 22, f(23) == 23, f(24) == 24, f(25) == 25,
            ....:     f(26) == 26, f(27) == 27, f(28) == 28, f(29) == 29,
            ....:     f(30) == 30, f(31) == 31], f, n, 8)
            sage: all([S[4*i] == S[2*i] and
            ....:      S[4*i + 1] == S[2*i] and
            ....:      S[4*i + 2] == S[2*i] and
            ....:      S[4*i + 3] == S[2*i - 16] for i in srange(8, 100)])
            True

        Zero-sequence with non-zero initial values::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 0, f(2*n + 1) == 0,
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3], f, n)
            Traceback (most recent call last):
            ...
            ValueError: Initial value for n = 0 does not match with the given recurrence relations.

        ::

            sage: Seq2.from_recurrence([
            ....:     f(2*n) == 0, f(2*n + 1) == 0,
            ....:     f(0) == 1, f(1) == 1, f(2) == 2, f(3) == 3], f, n, 2)
            2-regular sequence 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, ...
        """
        from sage.arith.srange import srange

        k = self.k
        M, m, coeffs, initial_values = self._parse_recurrence_(equations, function, var)
        recurrence_rules = self._get_parameters_from_recurrence_(
            M, m, coeffs, initial_values, offset)

        mu = [self._get_matrix_from_recurrence_(recurrence_rules, rem)
              for rem in srange(k)]

        seq = self(mu, self._get_left_from_recurrence_(recurrence_rules.dim),
                   self._get_right_from_recurrence_(recurrence_rules))

        if minimize:
            return seq.minimized()
        else:
            return seq
