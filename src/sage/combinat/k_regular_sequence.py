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

::

    sage: import logging
    sage: logging.basicConfig()

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

AUTHORS:

- Daniel Krenn (2016)

ACKNOWLEDGEMENT:

- Daniel Krenn is supported by the
  Austrian Science Fund (FWF): P 24644-N26.


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
from six import iteritems


def pad_right(T, length, zero=0):
    r"""
    Pad ``T`` to the right by ``zero``s to have
    at least the given ``length``.

    INPUT:

    - ``T`` -- A tuple, list or other iterable.

    - ``length`` -- a nonnegative integer.

    - ``zero`` -- (default: ``0``) the elements to pad with.

    OUTPUT:

    An object of the same type as ``T``.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import pad_right
        sage: pad_right((1,2,3), 10)
        (1, 2, 3, 0, 0, 0, 0, 0, 0, 0)
        sage: pad_right((1,2,3), 2)
        (1, 2, 3)

    TESTS::

        sage: pad_right([1,2,3], 10)
        [1, 2, 3, 0, 0, 0, 0, 0, 0, 0]
    """
    return T + type(T)(zero for _ in range(length - len(T)))


def value(D, k):
    r"""
    Return the value of the expansion with digits `D` in base `k`, i.e.

    .. MATH::

        \sum_{0\leq j < \operator{len}D} D[j] k^j.

    INPUT:

    - ``D`` -- a tuple or other iterable.

    - ``k`` -- the base.

    OUTPUT:

    An element in the common parent of the base `k` and of the entries
    of `D`.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import value
        sage: value(42.digits(7), 7)
        42
    """
    return sum(d * k**j for j, d in enumerate(D))


def split_interlace(n, k, p):
    r"""
    Split each digit in the `k`-ary expansion of `n` into `p` parts and
    return the value of the expansion obtained by each of these parts.

    INPUT:

    - ``n`` -- an integer.

    - ``k`` -- an integer specifying the base.

    - ``p`` -- a positive integer specifying in how many parts
      the input ``n`` is split. This has to be a divisor of ``k``.

    OUTPUT:

    A tuple of integers.

    EXAMPLES::

        sage: from sage.combinat.k_regular_sequence import split_interlace
        sage: [(n, split_interlace(n, 4, 2)) for n in srange(20)]
        [(0, (0, 0)), (1, (1, 0)), (2, (0, 1)), (3, (1, 1)),
         (4, (2, 0)), (5, (3, 0)), (6, (2, 1)), (7, (3, 1)),
         (8, (0, 2)), (9, (1, 2)), (10, (0, 3)), (11, (1, 3)),
         (12, (2, 2)), (13, (3, 2)), (14, (2, 3)), (15, (3, 3)),
         (16, (4, 0)), (17, (5, 0)), (18, (4, 1)), (19, (5, 1))]
        sage: [(n, split_interlace(n, 6, 3)) for n in srange(9)]
        [(0, (0, 0, 0)), (1, (1, 0, 0)), (2, (0, 1, 0)),
         (3, (1, 1, 0)), (4, (0, 0, 1)), (5, (1, 0, 1)),
         (6, (2, 0, 0)), (7, (3, 0, 0)), (8, (2, 1, 0))]

    TESTS::

        sage: split_interlace(42, 4, 3)
        Traceback (most recent call last):
        ...
        ValueError: p=3 is not a divisor of k=4.
    """
    if k % p != 0:
        raise ValueError('p={} is not a divisor of k={}.'.format(p, k))
    ki = k // p
    return tuple(value(D, ki)
                 for D in zip(*(d.digits(ki, padto=p)
                                for d in n.digits(k, padto=1))))


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
          from the left to the matrix product. If ``None``, then this
          multiplication is skipped.

        When created via the parent :class:`kRegularSequenceSpace`, then
        the following option is available.

        - ``transpose`` -- (default: ``False``) a boolean. If set, then
            each of the matrices in
            :meth:`mu <sage.combinat.recognizable_series.RecognizableSeries.mu>`
            is transposed. Additionally the vectors
            :meth`left <sage.combinat.recognizable_series.RecognizableSeries.left>`
            and
            :meth:`right <sage.combinat.recognizable_series.RecognizableSeries.right>`
            are switched.
            (This is done by calling :meth:`~sage.combinat.recognizable_series.RecognizableSeries.transposed`.)

        - ``heal`` -- (default: ``False``) a boolean. If set, then
          :meth:`healed` is called at the end of creating the element.

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
            sage: S = Seq2((M0, M1), vector([0, 1]), vector([1, 1]))
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


    @cached_method
    def is_healthy(self):
        r"""
        Return whether this `k`-regular sequence satisfies
        `\mu[0] \mathit{right} = \mathit{right}`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: S
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: S.is_healthy()
            False

        ::

            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C.is_healthy()
            True
        """
        from sage.rings.integer_ring import ZZ
        return (self.mu[ZZ(0)] * self.right) == self.right


    @cached_method
    def healed(self, minimize=True):
        r"""
        Return a `k`-regular sequence that satisfies
        `\mu[0] \mathit{right} = \mathit{right}`.

        INPUT:

        - ``minimize`` -- (default: ``True``) a boolean. If set, then
          :meth:`minimized` is called after the operation.

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)

        The following linear representation of `S` is chosen bad (is
        unhealty, see :meth:`is_healthy`), as `\mu(0)` applied on
        `\mathit{right}` does not equal `\mathit{right}`::

            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: S
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: S.is_healthy()
            False

        However, we can heal the sequence `S`::

            sage: H = S.healed()
            sage: H
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: H.mu[0], H.mu[1], H.left, H.right
            (
            [ 0  1]  [3 0]
            [-2  3], [6 0], (1, 0), (1, 1)
            )
            sage: H.is_healthy()
            True

        TESTS::

            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: H = S.healed(minimize=False)
            sage: H.mu[0], H.mu[1], H.left, H.right
            (
            [1 0]  [0 0]
            [0 2], [3 3], (1, 1), (1, 0)
            )
            sage: H.is_healthy()
            True

        ::

            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C.is_healthy()
            True
            sage: C.healed() is C
            True
        """
        if self.is_healthy():
            return self

        from sage.matrix.special import zero_matrix, identity_matrix
        from sage.modules.free_module_element import vector

        P = self.parent()
        dim = self.dimension()
        Z = zero_matrix(dim)
        I = identity_matrix(dim)

        itA = iter(P.alphabet())
        z = next(itA)
        mu = {z: I.augment(Z).stack(Z.augment(self.mu[z]))}
        mu.update((r, Z.augment(Z).stack(self.mu[r].augment(self.mu[r])))
                  for r in itA)

        result = P.element_class(
            P, mu,
            vector(2*tuple(self.left)),
            vector(tuple(self.right) + dim*(0,)))

        if minimize:
            return result.minimized()
        else:
            return result


    def subsequence(self, a, b, minimize=True):
        r"""
        Return the subsequence with indices `an+b` of this
        `k`-regular sequence.

        INPUT:

        - ``a`` -- a nonnegative integer.

        - ``b`` -- an integer.

          Alternatively, this is allowed to be a dictionary
          `b_j \mapsto c_j`. If so, the result will be the sum
          of all `c_j(an+b_j)`.

        - ``minimize`` -- (default: ``True``) a boolean. If set, then
          :meth:`minimized` is called after the operation.

        OUTPUT:

        A :class:`kRegularSequence`.

        .. NOTE::

            If `b` is negative (i.e., right-shift), then the
            coefficients when accessing negative indices are `0`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...

            sage: C.subsequence(2, 0)
            2-regular sequence 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, ...

            sage: S = C.subsequence(3, 1)
            sage: S
            2-regular sequence 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [ 0  1]  [ 6 -2]
            [-2  3], [10 -3], (1, 0), (1, 1)
            )

            sage: C.subsequence(3, 2)
            2-regular sequence 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, ...

        ::

            sage: S = C.subsequence(1, -1)
            sage: S
            2-regular sequence 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [ 0  1  0]  [ -2   2   0]
            [-2  3  0]  [  0   0   1]
            [-4  4  1], [ 12 -12   5], (1, 0, 0), (0, 0, 1)
            )

        We can build :meth:`backward_differences` manually by passing
        a dictionary for the parameter ``b``::

            sage: C.subsequence(1, {0: 1, -1: -1})
            2-regular sequence 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...

        TESTS::

            sage: C.subsequence(0, 4)
            2-regular sequence 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, ...
            sage: C.subsequence(1, 0) is C
            True

        The following test that the range for `c` in the code
        is sufficient::

            sage: C.subsequence(1, -1, minimize=False)
            2-regular sequence 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, ...
            sage: C.subsequence(1, -2, minimize=False)
            2-regular sequence 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, ...
            sage: C.subsequence(2, -1, minimize=False)
            2-regular sequence 0, 1, 3, 5, 7, 9, 11, 13, 15, 17, ...
            sage: C.subsequence(2, -2, minimize=False)
            2-regular sequence 0, 0, 2, 4, 6, 8, 10, 12, 14, 16, ...

            sage: C.subsequence(2, 21, minimize=False)
            2-regular sequence 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, ...
            sage: C.subsequence(2, 20, minimize=False)
            2-regular sequence 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, ...
            sage: C.subsequence(2, 19, minimize=False)
            2-regular sequence 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, ...
            sage: C.subsequence(2, -9, minimize=False)
            2-regular sequence 0, 0, 0, 0, 0, 1, 3, 5, 7, 9, ...

            sage: C.subsequence(3, 21, minimize=False)
            2-regular sequence 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, ...
            sage: C.subsequence(3, 20, minimize=False)
            2-regular sequence 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, ...
            sage: C.subsequence(3, 19, minimize=False)
            2-regular sequence 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, ...
            sage: C.subsequence(3, 18, minimize=False)
            2-regular sequence 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, ...

            sage: C.subsequence(10, 2, minimize=False)
            2-regular sequence 2, 12, 22, 32, 42, 52, 62, 72, 82, 92, ...
            sage: C.subsequence(10, 1, minimize=False)
            2-regular sequence 1, 11, 21, 31, 41, 51, 61, 71, 81, 91, ...
            sage: C.subsequence(10, 0, minimize=False)
            2-regular sequence 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, ...
            sage: C.subsequence(10, -1, minimize=False)
            2-regular sequence 0, 9, 19, 29, 39, 49, 59, 69, 79, 89, ...
            sage: C.subsequence(10, -2, minimize=False)
            2-regular sequence 0, 8, 18, 28, 38, 48, 58, 68, 78, 88, ...

        ::

            sage: C.subsequence(-1, 0)
            Traceback (most recent call last):
            ...
            ValueError: a=-1 is not nonnegative.

        The following linear representation of `S` is chosen bad (is
        unhealty, see :meth:`is_healthy`), as `\mu(0)` applied on
        `\mathit{right}` does not equal `\mathit{right}`::

            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: S
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...

        This leads to the wrong result
        ::

            sage: S.subsequence(1, -4)
            2-regular sequence 0, 0, 0, 0, 8, 12, 12, 18, 24, 36, ...

        We get the correct result by
        ::

            sage: S.healed().subsequence(1, -4)
            2-regular sequence 0, 0, 0, 0, 1, 3, 6, 9, 12, 18, ...
        """
        from sage.rings.integer_ring import ZZ
        zero = ZZ(0)
        a = ZZ(a)
        if not isinstance(b, dict):
            b = {ZZ(b): ZZ(1)}

        if a == 0:
            return sum(c_j * self[b_j] * self.parent().one_hadamard()
                       for b_j, c_j in iteritems(b))
        elif a == 1 and len(b) == 1 and zero in b:
            return b[zero] * self
        elif a < 0:
            raise ValueError('a={} is not nonnegative.'.format(a))

        from sage.arith.srange import srange
        from sage.matrix.constructor import Matrix
        from sage.modules.free_module_element import vector
        P = self.parent()
        A = P.alphabet()
        k = P.k
        dim = self.dimension()

        # Below, we use a dynamic approach to find the shifts of the
        # sequences in the kernel. According to [AS2003]_, the static range
        #    [min(b, 0), max(a, a + b))
        # suffices. However, it seems that the smaller set
        #    [min(b, 0), max(a, a + (b-1)//k + 1)) \cup {b}
        # suffices as well.
        kernel = list(b)

        def pad(T, d):
            di = kernel.index(d)
            return (di*dim)*(0,) + T
        def mu_line(r, i, c):
            d, f = (a*r + c).quo_rem(k)
            if d not in kernel:
                kernel.append(d)
            return pad(tuple(self.mu[f].rows()[i]), d)

        lines = dict((r, []) for r in A)
        ci = 0
        while ci < len(kernel):
            c = kernel[ci]
            for r in A:
                for i in srange(dim):
                    lines[r].append(mu_line(r, i, c))
            ci += 1

        ndim = len(kernel) * dim
        result = P.element_class(
            P,
            dict((r, Matrix([pad_right(row, ndim, zero=zero)
                             for row in lines[r]]))
                 for r in A),
            sum(c_j * vector(
                    pad_right(pad(tuple(self.left), b_j), ndim, zero=zero))
                for b_j, c_j in iteritems(b)),
            vector(sum((tuple(self.__getitem__(c, multiply_left=False))
                        if c >= 0 else dim*(zero,)
                        for c in kernel), tuple())))

        if minimize:
            return result.minimized()
        else:
            return result


    def backward_differences(self, **kwds):
        r"""
        Return the sequence of backward differences of this
        `k`-regular sequence.

        INPUT:

        - ``minimize`` -- (default: ``True``) a boolean. If set, then
          :meth:`minimized` is called after the operation.

        OUTPUT:

        A :class:`kRegularSequence`.

        .. NOTE::

            The coefficient to the index `-1` is `0`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
            sage: C.backward_differences()
            2-regular sequence 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...

        ::

            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...
            sage: E.backward_differences()
            2-regular sequence 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, ...
        """
        return self.subsequence(1, {0: 1, -1: -1}, **kwds)


    def forward_differences(self, **kwds):
        r"""
        Return the sequence of forward differences of this
        `k`-regular sequence.

        INPUT:

        - ``minimize`` -- (default: ``True``) a boolean. If set, then
          :meth:`minimized` is called after the operation.

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
            sage: C.forward_differences()
            2-regular sequence 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...

        ::

            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...
            sage: E.forward_differences()
            2-regular sequence -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, ...
        """
        return self.subsequence(1, {1: 1, 0: -1}, **kwds)


    def partial_sums(self, include_n=False, minimize=True):
        r"""
        Return the sequence of partial sums of this
        `k`-regular sequence. That is, the `n`th entry of the result
        is the sum of the first `n` entries in the original sequence.

        INPUT:

        - ``include_n`` -- (default: ``False``) a boolean. If set, then
          the `n`th entry of the result is the sum of the entries up
          to index `n` (included).

        - ``minimize`` -- (default: ``True``) a boolean. If set, then
          :meth:`minimized` is called after the operation.

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)

            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...
            sage: E.partial_sums()
            2-regular sequence 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, ...
            sage: E.partial_sums(include_n=True)
            2-regular sequence 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, ...

        ::

            sage: C = Seq2((Matrix([[2, 0], [2, 1]]), Matrix([[0, 1], [-2, 3]])),
            ....:          vector([1, 0]), vector([0, 1]))
            sage: C
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
            sage: C.partial_sums()
            2-regular sequence 0, 0, 1, 3, 6, 10, 15, 21, 28, 36, ...
            sage: C.partial_sums(include_n=True)
            2-regular sequence 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, ...

        The following linear representation of `S` is chosen bad (is
        unhealty, see :meth:`is_healthy`), as `\mu(0)` applied on
        `\mathit{right}` does not equal `\mathit{right}`::

            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: S
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...

        Therefore, building partial sums produces a wrong result::

            sage: H = S.partial_sums(include_n=True, minimize=False)
            sage: H
            2-regular sequence 1, 5, 16, 25, 62, 80, 98, 125, 274, 310, ...
            sage: H = S.partial_sums(minimize=False)
            sage: H
            2-regular sequence 0, 2, 10, 16, 50, 62, 80, 98, 250, 274, ...

        We can :meth:`~kRegularSequenceSpace.guess` the correct representation::

            sage: from itertools import islice
            sage: L = []; ps = 0
            sage: for s in islice(S, 110):
            ....:     ps += s
            ....:     L.append(ps)
            sage: G = Seq2.guess(lambda n: L[n])
            sage: G
            2-regular sequence 1, 4, 10, 19, 31, 49, 67, 94, 118, 154, ...
            sage: G.mu[0], G.mu[1], G.left, G.right
            (
            [  0   1   0   0]  [  0   0   1   0]
            [  0   0   0   1]  [ -5   3   3   0]
            [ -5   5   1   0]  [ -5   0   6   0]
            [ 10 -17   0   8], [-30  21  10   0], (1, 0, 0, 0), (1, 1, 4, 1)
            )
            sage: G.minimized().dimension() == G.dimension()
            True

        Or we heal the sequence `S` first::

            sage: S.healed().partial_sums(include_n=True, minimize=False)
            2-regular sequence 1, 4, 10, 19, 31, 49, 67, 94, 118, 154, ...
            sage: S.healed().partial_sums(minimize=False)
            2-regular sequence 0, 1, 4, 10, 19, 31, 49, 67, 94, 118, ...

        TESTS::

            sage: E = Seq2((Matrix([[0, 1], [0, 1]]), Matrix([[0, 0], [0, 1]])),
            ....:          vector([1, 0]), vector([1, 1]))
            sage: E.mu[0], E.mu[1], E.left, E.right
            (
            [0 1]  [0 0]
            [0 1], [0 1], (1, 0), (1, 1)
            )
            sage: P = E.partial_sums(minimize=False)
            sage: P.mu[0], P.mu[1], P.left, P.right
            (
            [ 0  1  0  0]  [0 1 0 0]
            [ 0  2  0 -1]  [0 2 0 0]
            [ 0  0  0  1]  [0 0 0 0]
            [ 0  0  0  1], [0 0 0 1], (1, 0, -1, 0), (1, 1, 1, 1)
            )
        """
        from sage.matrix.constructor import Matrix
        from sage.matrix.special import zero_matrix
        from sage.modules.free_module_element import vector

        P = self.parent()
        A = P.alphabet()
        k = P.k
        dim = self.dimension()

        B = dict((r, sum(self.mu[a] for a in A[r:])) for r in A)
        Z = zero_matrix(dim)
        B[k] = Z
        W = B[0].stack(Z)

        result = P.element_class(
            P,
            dict((r, W.augment((-B[r+1]).stack(self.mu[r])))
                 for r in A),
            vector(tuple(self.left) +
                   (dim*(0,) if include_n else tuple(-self.left))),
            vector(2*tuple(self.right)))

        if minimize:
            return result.minimized()
        else:
            return result


def _pickle_kRegularSequenceSpace(k, coefficients, category):
    r"""
    Pickle helper.

    TESTS::

        sage: Seq2 = kRegularSequenceSpace(2, ZZ)
        sage: from sage.combinat.k_regular_sequence import _pickle_kRegularSequenceSpace
        sage: _pickle_kRegularSequenceSpace(
        ....:     Seq2.k, Seq2.coefficients(), Seq2.category())
        Space of 2-regular sequences over Integer Ring
    """
    return kRegularSequenceSpace(k, coefficients, category=category)


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
            Category of modules over Integer Ring

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

        ::

            sage: from itertools import islice
            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: TestSuite(Seq2).run(  # long time
            ....:    verbose=True,
            ....:    elements=tuple(islice(Seq2.some_elements(), 4)))
            running ._test_additive_associativity() . . . pass
            running ._test_an_element() . . . pass
            running ._test_cardinality() . . . pass
            running ._test_category() . . . pass
            running ._test_elements() . . .
              Running the test suite of self.an_element()
              running ._test_category() . . . pass
              running ._test_eq() . . . pass
              running ._test_new() . . . pass
              running ._test_nonzero_equal() . . . pass
              running ._test_not_implemented_methods() . . . pass
              running ._test_pickling() . . . pass
              pass
            running ._test_elements_eq_reflexive() . . . pass
            running ._test_elements_eq_symmetric() . . . pass
            running ._test_elements_eq_transitive() . . . pass
            running ._test_elements_neq() . . . pass
            running ._test_eq() . . . pass
            running ._test_new() . . . pass
            running ._test_not_implemented_methods() . . . pass
            running ._test_pickling() . . . pass
            running ._test_some_elements() . . . pass
            running ._test_zero() . . . pass

        .. SEEALSO::

            :doc:`k-regular sequence <k_regular_sequence>`,
            :class:`kRegularSequence`.
        """
        self.k = k
        super(kRegularSequenceSpace, self).__init__(*args, **kwds)


    def __reduce__(self):
        r"""
        Pickling support.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: loads(dumps(Seq2))  # indirect doctest
            Space of 2-regular sequences over Integer Ring
        """
        return _pickle_kRegularSequenceSpace, \
            (self.k, self.coefficients(), self.category())


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


    def some_elements(self):
        r"""
        Return some elements of this `k`-regular sequence.

        See :class:`TestSuite` for a typical use case.

        OUTPUT:

        An iterator.

        EXAMPLES::

            sage: tuple(kRegularSequenceSpace(2, ZZ).some_elements())
            (2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...,
             2-regular sequence 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, -1, 0, 0, 1, -2, -1, ...,
             2-regular sequence 2, -1, 0, 0, 0, -1, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, 5, 0, 0, 1, -33, 5, ...,
             2-regular sequence -5, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence -59, -20, 0, -20, 0, 0, 0, -20, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, -10, 0, 0, 1, -60, -10, ...,
             2-regular sequence 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 142, -30, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, 44, 0, 0, 1, -1375, 44, ...,
             2-regular sequence -25, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence -275, -50, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, -36, 0, 0, 1, -432, -36, ...,
             2-regular sequence 46, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 478, 66, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, 178, 0, 0, 1, -15288, 178, ...,
             2-regular sequence -70, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence -756, -91, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, -88, 0, 0, 1, -1760, -88, ...,
             2-regular sequence 108, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1132, 112, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, 500, 0, 0, 1, -91260, 500, ...,
             2-regular sequence -150, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence -1608, -144, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 1, 1, 0, 1, -175, 0, 0, 1, -5250, -175, ...,
             2-regular sequence 210, 0, 0, 0, 0, 0, 0, 0, 0, 0, ...,
             2-regular sequence 2210, 170, 0, 0, 0, 0, 0, 0, 0, 0, ...)
        """
        return super(kRegularSequenceSpace, self).some_elements(heal=True)


    def _element_constructor_(self, *args, **kwds):
        r"""
        Return a `k`-regular sequence.

        See :class:`kRegularSequenceSpace` for details.

        TESTS::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]),
            ....:      heal=True)
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
        """
        heal = kwds.pop('heal', False)
        element = super(kRegularSequenceSpace, self)._element_constructor_(*args, **kwds)
        if heal:
            element = element.healed()
        elif not element.is_healthy():
            import logging
            logger = logging.getLogger(__name__)
            logger.warning('Unhealthy sequence: mu[0]*right != right. '
                           'Results might be wrong. '
                           'Use heal=True or method .healed() '
                           'for correcting this.')
        return element


    def guess(self, f, n_max=100, max_dimension=10, sequence=None):
        r"""
        Guess a `k`-regular sequence of `(f(n))_{n\geq0}`.

        INPUT:

        - ``f`` -- a function (callable) which determines the sequence.
          It takes nonnegative integers as an input.

        - ``n_max`` -- (default: ``100``) a positive integer. The resulting
          `k`-regular sequence coincides with `f` on the first ``n_max``
          terms.

        - ``max_dimension`` -- (default: ``10``) a positive integer specifying
          the maxium dimension which is tried when guessing the sequence.

        - ``sequence`` -- (default: ``None``) a `k`-regular sequence used
          for bootstrapping this guessing.

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES:

        Binary sum of digits::

            sage: @cached_function
            ....: def s(n):
            ....:     if n == 0:
            ....:         return 0
            ....:     return s(n//2) + ZZ(is_odd(n))
            sage: all(s(n) == sum(n.digits(2)) for n in srange(10))
            True
            sage: [s(n) for n in srange(10)]
            [0, 1, 1, 2, 1, 2, 2, 3, 1, 2]

        Variant 1::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: import logging
            sage: logging.getLogger().setLevel(logging.INFO)
            sage: S1 = Seq2.guess(s)
            INFO:...:including f_{1*m+0}
            INFO:...:M_0: f_{2*m+0} = (1) * X_m
            INFO:...:including f_{2*m+1}
            INFO:...:M_1: f_{2*m+1} = (0, 1) * X_m
            INFO:...:M_0: f_{4*m+1} = (0, 1) * X_m
            INFO:...:M_1: f_{4*m+3} = (-1, 2) * X_m
            sage: S1
            2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
            sage: S1.mu[0], S1.mu[1], S1.left, S1.right
            (
            [1 0]  [ 0  1]
            [0 1], [-1  2], (1, 0), (0, 1)
            )
            sage: logging.getLogger().setLevel(logging.WARN)

        Variant 2::

            sage: C = Seq2((Matrix([[1]]), Matrix([[1]])), vector([1]), vector([1])); C
            2-regular sequence 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ...
            sage: S2 = Seq2.guess(s, sequence=C)
            sage: S2
            2-regular sequence 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, ...
            sage: S2.mu[0], S2.mu[1], S2.left, S2.right
            (
            [1 0]  [1 0]
            [0 1], [1 1], (0, 1), (1, 0)
            )

        The sequence of all natural numbers::

            sage: S = Seq2.guess(lambda n: n)
            sage: S
            2-regular sequence 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [2 0]  [ 0  1]
            [2 1], [-2  3], (1, 0), (0, 1)
            )

        The indicator function of the even integers::

            sage: S = Seq2.guess(lambda n: ZZ(is_even(n)))
            sage: S
            2-regular sequence 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [0 1]  [0 0]
            [0 1], [0 1], (1, 0), (1, 1)
            )

        The indicator function of the odd integers::

            sage: S = Seq2.guess(lambda n: ZZ(is_odd(n)))
            sage: S
            2-regular sequence 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            (
            [0 0]  [0 1]
            [0 1], [0 1], (1, 0), (0, 1)
            )

        The following linear representation of `S` is chosen bad (is
        unhealty, see :meth:`is_healthy`), as `\mu(0)` applied on
        `\mathit{right}` does not equal `\mathit{right}`::

            sage: S = Seq2((Matrix([2]), Matrix([3])), vector([1]), vector([1]))
            WARNING:...:Unhealthy sequence: mu[0]*right != right.
                        Results might be wrong. Use heal=True or
                        method .healed() for correcting this.
            sage: S
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: S.is_healthy()
            False

        However, we can :meth:`~kRegularSequenceSpace.guess` a `2`-regular sequence of dimension `2`::

            sage: G = Seq2.guess(lambda n: S[n])
            sage: G
            2-regular sequence 1, 3, 6, 9, 12, 18, 18, 27, 24, 36, ...
            sage: G.mu[0], G.mu[1], G.left, G.right
            (
            [ 0  1]  [3 0]
            [-2  3], [6 0], (1, 0), (1, 1)
            )

            sage: G == S.healed()
            True

        TESTS::

            sage: S = Seq2.guess(lambda n: 2, sequence=C)
            sage: S
            2-regular sequence 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, ...
            sage: S.mu[0], S.mu[1], S.left, S.right
            ([1], [1], (2), (1))

        We :meth:`~kRegularSequenceSpace.guess` some partial sums sequences::

            sage: S = Seq2((Matrix([1]), Matrix([2])), vector([1]), vector([1]))
            sage: S
            2-regular sequence 1, 2, 2, 4, 2, 4, 4, 8, 2, 4, ...
            sage: from itertools import islice
            sage: L = []; ps = 0
            sage: for s in islice(S, 110):
            ....:     ps += s
            ....:     L.append(ps)
            sage: G = Seq2.guess(lambda n: L[n])
            sage: G
            2-regular sequence 1, 3, 5, 9, 11, 15, 19, 27, 29, 33, ...
            sage: G.mu[0], G.mu[1], G.left, G.right
            (
            [ 0  1]  [3 0]
            [-3  4], [3 2], (1, 0), (1, 1)
            )
            sage: G == S.partial_sums(include_n=True)
            True

        ::

            sage: Seq3 = kRegularSequenceSpace(3, QQ)
            sage: S = Seq3((Matrix([1]), Matrix([3]), Matrix([2])), vector([1]), vector([1]))
            sage: S
            3-regular sequence 1, 3, 2, 3, 9, 6, 2, 6, 4, 3, ...
            sage: from itertools import islice
            sage: L = []; ps = 0
            sage: for s in islice(S, 110):
            ....:     ps += s
            ....:     L.append(ps)
            sage: G = Seq3.guess(lambda n: L[n])
            sage: G
            3-regular sequence 1, 4, 6, 9, 18, 24, 26, 32, 36, 39, ...
            sage: G.mu[0], G.mu[1], G.mu[2], G.left, G.right
            (
            [ 0  1]  [18/5  2/5]  [ 6  0]
            [-6  7], [18/5 27/5], [24  2], (1, 0), (1, 1)
            )
            sage: G == S.partial_sums(include_n=True)
            True
        """
        import logging
        logger = logging.getLogger(__name__)

        from sage.arith.srange import srange, xsrange
        from sage.matrix.constructor import Matrix
        from sage.misc.mrange import cantor_product
        from sage.modules.free_module_element import vector

        k = self.k
        domain = self.coefficients()
        if sequence is None:
            mu = [[] for _ in srange(k)]
            seq = lambda m: tuple()
        else:
            mu = [M.rows() for M in sequence.mu]
            seq = lambda m: sequence.left * sequence._mu_of_word_(
                self._n_to_index_(m))

        zero = domain(0)
        one = domain(1)

        def values(m, lines):
            return tuple(seq(m)) + tuple(f(k**t_R * m + r_R) for t_R, r_R, s_R in lines)

        @cached_function(key=lambda lines: len(lines))  # we assume that existing lines are not changed (we allow appending of new lines)
        def some_inverse_U_matrix(lines):
            d = len(seq(0)) + len(lines)

            for m_indices in cantor_product(xsrange(n_max), repeat=d, min_slope=1):
                U = Matrix(domain, d, d, [values(m, lines) for m in m_indices]).transpose()
                try:
                    return U.inverse(), m_indices
                except ZeroDivisionError:
                    pass
            else:
                raise RuntimeError

        def guess_linear_dependence(t_L, r_L, lines):
            iU, m_indices = some_inverse_U_matrix(lines)
            X_L = vector(f(k**t_L * m + r_L) for m in m_indices)
            return X_L * iU

        def verify_linear_dependence(t_L, r_L, linear_dependence, lines):
            return all(f(k**t_L * m + r_L) ==
                       linear_dependence * vector(values(m, lines))
                       for m in xsrange(0, (n_max - r_L) // k**t_L + 1))

        def find_linear_dependence(t_L, r_L, lines):
            linear_dependence = guess_linear_dependence(t_L, r_L, lines)
            if not verify_linear_dependence(t_L, r_L, linear_dependence, lines):
                raise ValueError
            return linear_dependence

        left = None
        if seq(0):
            try:
                solution = find_linear_dependence(0, 0, [])
            except ValueError:
                pass
            else:
                left = vector(solution)

        to_branch = []
        lines = []
        def include(line):
            to_branch.append(line)
            lines.append(line)
            t, r, s = line
            logger.info('including f_{%s*m+%s}', k**t, r)

        if left is None:
            line_L = (0, 0, 0)  # entries (t, r, s) --> k**t * m + r, belong to M_s
            include(line_L)
            left = vector((len(seq(0)) + len(lines)-1)*(zero,) + (one,))

        while to_branch:
            line_R = to_branch.pop(0)
            t_R, r_R, s_R = line_R
            if t_R >= max_dimension:
                raise RuntimeError

            t_L = t_R + 1
            for s_L in srange(k):
                r_L = k**t_R * s_L + r_R
                line_L = t_L, r_L, s_L

                try:
                    solution = find_linear_dependence(t_L, r_L, lines)
                except ValueError:
                    include(line_L)
                    solution = (len(lines)-1)*(zero,) + (one,)
                logger.info('M_%s: f_{%s*m+%s} = %s * X_m',
                            s_L, k**t_L, r_L, solution)
                mu[s_L].append(solution)

        d = len(seq(0)) + len(lines)
        mu = tuple(Matrix(domain, [pad_right(tuple(row), d, zero=zero) for row in M])
                         for M in mu)
        right = vector(values(0, lines))
        left = vector(pad_right(tuple(left), d, zero=zero))
        return self(mu, left, right)


    def _parse_recursions_(self, equations, function, var, n0=0):
        r"""Parse recursion equations as admissible in :meth:`~.recursions`.

        INPUT:

        - ``equations`` -- see :meth:`~recursions`.

        - ``function`` -- see :meth:`~recursions`.

        - ``var`` -- see :meth:`~recursions`.

        OUTPUT:

        A namedtuple ``recursion_rules`` containing all the
        significant information obtained from ``equations``.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2._parse_recursions_([
            ....:     f(4*n) == f(2*n) + 2*f(2*n + 1) + 3*f(2*n - 2),
            ....:     f(4*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n - 2),
            ....:     f(4*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n - 2),
            ....:     f(4*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n - 2),
            ....:     f(0) == 1, f(1) == 2, f(2) == 1], f, n, 42)
            recursion_rules(M=2, m=1, l=-2, u=1, ll=-6, uu=3, dim=11,
            coeffs={(0, 1): 2, (0, 0): 1, (3, 1): 11, (3, 0): 10, (2, -2): 9,
            (2, 1): 8, (2, 0): 7, (3, -2): 12, (0, -2): 3, (1, 0): 4, (1, -2): 6,
            (1, 1): 5}, initial_values={0: 1, 1: 2, 2: 1}, n0=42)

        Stern--Brocot Sequence::

            sage: Seq2._parse_recursions_([f(2*n) == f(n),
            ....:    f(2*n + 1) == f(n) + f(n + 1), f(0) == 0,
            ....:    f(1) == 1, f(2) == 1], f, n)
            recursion_rules(M=1, m=0, l=0, u=1, ll=0, uu=2, dim=3,
            coeffs={(1, 0): 1, (0, 0): 1, (1, 1): 1},
            initial_values={0: 0, 1: 1, 2: 1}, n0=0)

        TESTS:

            The following tests check that the equations are well-formed::

                sage: Seq2._parse_recursions_([], f, n)
                Traceback (most recent call last):
                ...
                ValueError: List of recursion equations is empty.

            ::

                sage: Seq2._parse_recursions_([f(4*n + 1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(4*n + 1) is not an equation with ==.

            ::

                sage: Seq2._parse_recursions_([f(2*n) + 1 == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(2*n) + 1 is not an evaluation of f.

            ::

                sage: Seq2._parse_recursions_([f(2*n, 5) == 3], f, n)
                Traceback (most recent call last):
                ...
                ValueError: f(2*n, 5) does not have one argument.

            ::

                sage: Seq2._parse_recursions_([f(1/n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ....:
                ValueError: 1/n + 1 is not a polynomial in n.

            ::

                sage: Seq2._parse_recursions_([f(4*n^2) == f(2*n^2)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 4*n^2 is not a polynomial of degree smaller 2.

            ::

                sage: Seq2._parse_recursions_([f(3*n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: Seq2._parse_recursions_([f(n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 is less than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n), f(2*n) == 0], f, n)
                Traceback (most recent call last):
                ...
                ValueError: There are more than one recursions for f(2*n).

            ::

                sage: Seq2._parse_recursions_([f(2*n + 2) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n - 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: -1 is smaller than 0.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 2*n], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2*n does not contain f.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/2*f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/2 is not a valid coefficient.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/f(n) is not a valid right hand side.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == 1/f(n) + 2*f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/f(n) is not a valid summand.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n + 1/2)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: n + 1/2 does not have integer coefficients.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(1/2*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1/2*n does not have integer coefficients.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n^2 + 1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: n^2 + 1 does not have degree 1.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(1)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not have degree 1.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(2*n) + f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not equal 2.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(2*n), f(4*n + 1) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 1 does not equal 2.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(3*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 3 is not a power of 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(4*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 4 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(2*n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: 2 is not smaller than 2.

            ::

                sage: Seq2._parse_recursions_([f(2*n) == f(n)], f, n)
                Traceback (most recent call last):
                ...
                ValueError: Recursions for [f(2*n + 1)] are missing.

            ::

                sage: Seq2._parse_recursions_([f(4*n) == f(n), f(4*n + 3) == 0], f, n)
                Traceback (most recent call last):
                ...
                ValueError: Recursions for [f(4*n + 1), f(4*n + 2)] are missing.

            Finally, also for the zero-sequence the output is as expected::

                sage: Seq2._parse_recursions_([f(2*n) == 0, f(2*n + 1) == 0], f, n)
                recursion_rules(M=1, m=0, l=0, u=0, ll=0, uu=0, dim=1,
                coeffs={}, initial_values={}, n0=0)
        """
        from collections import namedtuple

        from sage.arith.srange import srange
        from sage.functions.log import log
        from sage.functions.other import ceil, floor
        from sage.rings.integer_ring import ZZ
        from sage.symbolic.operators import add_vararg, mul_vararg, operator

        k = self.k
        base_ring = self.base()
        indices_right = [0]
        coeffs = {}
        initial_values = {}
        remainders = []

        def _parse_multiplication_(op):
            if op.operator() != mul_vararg or len(op.operands()) != 2:
                raise ValueError("")
            operands = op.operands()
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
            try:
                poly = ZZ[var](op.operands()[0])
            except TypeError:
                raise ValueError('%s does not have integer coefficients.'
                                 % (op.operands()[0],))
            if poly.degree() != 1:
                raise ValueError("%s does not have degree 1."
                                 % (poly,))
            d, base_power_m = list(poly)
            m = log(base_power_m, base=k)
            return [coeff, m, d]

        if not equations:
            raise ValueError("List of recursion equations is empty.")

        for eq in equations:
            if eq.operator() != operator.eq:
                raise ValueError("%s is not an equation with ==."  % eq)
            left_side, right_side = eq.operands()
            if left_side.operator() != function:
                raise ValueError("%s is not an evaluation of %s."
                                 % (left_side, function))
            if  len(left_side.operands()) != 1:
                raise ValueError("%s does not have one argument." %
                                 (left_side,))
            try:
                polynomial_left = base_ring[var](left_side.operands()[0])
            except Exception:
                raise ValueError("%s is not a polynomial "
                                 "in %s." % (left_side.operands()[0], var))
            if polynomial_left.degree()  > 1:
                raise ValueError("%s is not a polynomial of degree smaller 2."
                                 % (polynomial_left,))
            if polynomial_left in base_ring and right_side in base_ring:
                initial_values.update({polynomial_left: right_side})
            else:
                poly_left = ZZ[var](left_side.operands()[0])
                [r, base_power_M] = list(poly_left)
                M_new = log(base_power_M, base=k)
                try:
                    if M != log(base_power_M, base=k):
                        raise ValueError("%s does not equal %s."
                                         % (base_power_M, k^M))
                except NameError:
                    M = M_new
                    if M not in ZZ:
                        raise ValueError("%s is not a power of %s."
                                         % (base_power_M, k))
                    if M < 1:
                        raise ValueError("%s is less than %s."
                                         % (base_power_M, k))
                if r in remainders:
                    raise ValueError("There are more than one recursions for %s."
                                     % (left_side,))
                if r not in ZZ:
                    raise ValueError("%s is not an integer." % (r,))
                if r >= k**M:
                    raise ValueError("%s is not smaller than %s." % (r, k**M))
                if r < 0:
                    raise ValueError("%s is smaller than 0." % (r,))
                remainders.append(r)

                if right_side != 0:
                    if (len(right_side.operands()) == 1 and right_side.operator() == function
                        or right_side.operator() == mul_vararg and len(right_side.operands()) == 2):
                        summands = [right_side]
                    elif right_side.operator() == add_vararg:
                        summands = right_side.operands()
                    else: # check this again
                        raise ValueError("%s is not a valid right hand side."
                                         % (right_side,))
                    for summand in summands:
                        coeff, new_m, d = _parse_one_summand_(summand)
                        if coeff not in base_ring:
                            raise ValueError("%s is not a valid coefficient."
                                             % (coeff,))
                        try:
                            if m != new_m:
                                raise ValueError("%s does not equal %s."
                                                 % (k**new_m, k**m))
                        except NameError:
                            m = new_m
                            if m not in ZZ:
                                raise ValueError("%s is not a power of %s."
                                                 % (k**m, k))
                            if M <= m:
                                raise ValueError("%s is not smaller than %s."
                                                 % (k**m, k**M))

                        indices_right.append(d)
                        coeffs.update({(r, d): coeff})

        remainders.sort()
        if remainders != srange(k**M):
            missing_equations = [function(k**M*var + r)
                                 for r in srange(k**M)
                                 if r not in remainders]
            raise ValueError("Recursions for %s are missing."
                             % missing_equations)

        if not coeffs:
            m = M - 1

        l = min(indices_right)
        u = max(indices_right)
        ll = (floor((l*k**(M-m) - k**M + 1)/(k**(M-m) - 1)) + 1)*(l < 0)
        uu = max([ceil((u*k**(M-m) + k**M - k**m)/(k**(M-m) - 1)) - 1, k**m - 1])
        n1 = n0 - floor(ll/q^M)
        dim = (k**M - 1)/(k - 1) + (M - m)*(uu - ll - k**m + 1)

        recursion_rules = namedtuple('recursion_rules',
                                     ['M', 'm', 'l', 'u',
                                      'll', 'uu', 'dim',
                                      'coeffs', 'initial_values', 'n0'])

        return recursion_rules(M=M, m=m, l=l, u=u, ll=ll, uu=uu, dim=dim,
                               coeffs=coeffs, initial_values=initial_values, n0=n0, n1=n1)


    def _get_matrix_from_recursions_(self, recursion_rules, rem, function, var):
        r"""
        Construct the matrix for remainder ``rem`` of the linear
        representation of the sequence induced by ``recursion_rules``.

        INPUT:

        - ``recursion_rules`` -- A namedtuple generated by
          :meth:`~_parse_recursions_`.

        - ``rem`` -- An integer between ``0`` and ``k - 1``.

        OUTPUT:

        A matrix.

        EXAMPLES:

        The following example illustrates how the coefficients in the
        right hand sides of the recursions correspond to the entries of
        the matrices. ::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: rules = Seq2._parse_recursions_([
            ....:     f(8*n) == 1*f(2*n) + 2*f(2*n + 1) + 3*f(2*n + 2),
            ....:     f(8*n + 1) == 4*f(2*n) + 5*f(2*n + 1) + 6*f(2*n + 2),
            ....:     f(8*n + 2) == 7*f(2*n) + 8*f(2*n + 1) + 9*f(2*n + 2),
            ....:     f(8*n + 3) == 10*f(2*n) + 11*f(2*n + 1) + 12*f(2*n + 2),
            ....:     f(8*n + 4) == 13*f(2*n) + 14*f(2*n + 1) + 15*f(2*n + 2),
            ....:     f(8*n + 5) == 16*f(2*n) + 17*f(2*n + 1) + 18*f(2*n + 2),
            ....:     f(8*n + 6) == 19*f(2*n) + 20*f(2*n + 1) + 21*f(2*n + 2),
            ....:     f(8*n + 7) == 22*f(2*n) + 23*f(2*n + 1) + 24*f(2*n + 2),],
            ....:     f, n)
            sage: Seq2._get_matrix_from_recursions_(rules, 0, f, n)
            [ 0  1  0  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  1  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  1  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  1  2  3  0  0  0  0  0  0  0  0  0]
            [ 0  4  5  6  0  0  0  0  0  0  0  0  0]
            [ 0  7  8  9  0  0  0  0  0  0  0  0  0]
            [ 0 10 11 12  0  0  0  0  0  0  0  0  0]
            [ 0 13 14 15  0  0  0  0  0  0  0  0  0]
            [ 0 16 17 18  0  0  0  0  0  0  0  0  0]
            [ 0 19 20 21  0  0  0  0  0  0  0  0  0]
            sage: Seq2._get_matrix_from_recursions_(rules, 1, f, n)
            [ 0  0  1  0  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  1  0  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  1  0  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  1  0  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  1  0]
            [ 0  0  0  0  0  0  0  0  0  0  0  0  1]
            [ 0 13 14 15  0  0  0  0  0  0  0  0  0]
            [ 0 16 17 18  0  0  0  0  0  0  0  0  0]
            [ 0 19 20 21  0  0  0  0  0  0  0  0  0]
            [ 0 22 23 24  0  0  0  0  0  0  0  0  0]
            [ 0  0  0  1  2  3  0  0  0  0  0  0  0]
            [ 0  0  0  4  5  6  0  0  0  0  0  0  0]
            [ 0  0  0  7  8  9  0  0  0  0  0  0  0]

        Stern--Brocot Sequence::

            sage: SB_rules = Seq2._parse_recursions_([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1], f, n)
            sage: Seq2._get_matrix_from_recursions_(SB_rules, 0, f, n)
            [1 0 0]
            [1 1 0]
            [0 1 0]
            sage: Seq2._get_matrix_from_recursions_(SB_rules, 1, f, n)
            [1 1 0]
            [0 1 0]
            [0 1 1]

        """
        from sage.arith.srange import srange
        from sage.matrix.constructor import Matrix
        from sage.matrix.matrix_space import MatrixSpace
        from sage.matrix.special import block_matrix, zero_matrix
        from sage.modules.free_module_element import vector

        k = self.k
        M = recursion_rules.M
        m = recursion_rules.m
        l = recursion_rules.l
        u = recursion_rules.u
        ll = recursion_rules.ll
        uu = recursion_rules.uu
        n0 = recursion_rules.n0
        dim = recursion_rules.dim
        coeffs = recursion_rules.coeffs
        initial_values = recursion_rules.initial_values
        n0 = recursion_rules.n0

        mat = []
        current_shift = 0

        for base_power in srange(m - 1):
            current_shift += k**base_power
            for d in srange(k**base_power):
                row = dim*[0]
                dd = k**base_power*rem + d
                index = current_shift + dd
                row[index] = 1
                mat.append(row)

        if m > 0:
            current_shift += k**(m-1)
            final_shift = current_shift - ll
            for d in srange(k**(m-1)):
                row = dim*[0]
                dd = k**(m-1)*rem + d
                index = current_shift + dd - ll
                row[index] = 1
                mat.append(row)
        else:
            final_shift = -ll

        for base_power in srange(m, M - 1):
            current_shift += k**base_power - k**m + uu - 2*ll + 1
            for d in srange(ll, k**base_power - k**m + uu + 1):
                row = dim*[0]
                dd = k**base_power*rem + d
                index = current_shift + dd
                row[index] = 1
                mat.append(row)

        for d in srange(ll, k**(M-1) - k**m + uu + 1):
            dd, r = d.quo_rem(k**M)
            rr = k**(M-1)*rem + r
            row = dim*[0]
            if rr < k**M:
                for i in range(l, u + 1):
                    try:
                        row[k**m*dd + i + final_shift] = coeffs[(rr, i)]
                    except KeyError:
                        pass
            else:
                for i in srange(l, u + 1):
                    try:
                        row[(dd + 1)*k**m + i + final_shift] = coeffs[(rr - k**M, i)]
                    except KeyError:
                        pass
            mat.append(row)

        mat = Matrix(mat)

        if n0 == 0:
            return mat
        else:
            arguments = [k**j*var + d for j in srange(m) for d in srange(k**j)] + \
                        [k**j*var + d for j in srange(m, M) for d in srange(ll, k**j - k**m + uu + 1)]
            W = []
            for i in srange(n0):
                v_eval_i = []
                v_eval_ki_plus_r = []
                for a in arguments:
                    try:
                        temp = a.substitute(var==i)
                        v_eval_i.append(initial_values[temp])
                        temp = a.substitute(var==k*i+rem)
                        v_eval_ki_plus_r.append(initial_values[temp])
                    except KeyError:
                        raise ValueError('Initial value %s is missing.'
                                         % (function(temp),))
                W.append(list(vector(v_eval_ki_plus_r) - mat*vector(v_eval_i)))

            J = []
            for i in srange(n0):
                J.append([int(i >= rem and i % k == rem and j*k == i - rem) for j in srange(n0)])

            Mat = MatrixSpace(self.base_ring(), dim + n0, dim + n0) 
            return Mat(block_matrix([[mat, Matrix(W).transpose()],
                                     [zero_matrix(n0, dim), Matrix(J)]]))

    
    def _get_left_from_recursions_(self, dim):
        r"""
        Construct the vector ``left`` of the linear representation of
        recursive sequences.

        INPUT:

        - ``dim`` -- A positive integer.

        OUTPUT:

        A vector.

        EXAMPLES::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: Seq2._get_left_from_recursions_(5)
            (1, 0, 0, 0, 0)
        """
        from sage.modules.free_module_element import vector

        return vector([1] + (dim - 1)*[0])


    def _get_right_from_recursions_(self, recursion_rules, function):
        r"""
        Construct the vector ``right`` of the linear
        representation of the sequence induced by ``recursion_rules``.

        INPUT:

        - ``recursion_rules`` -- A namedtuple generated by
          :meth:`~_parse_recursions_`.

        - ``function`` -- A function.

        OUTPUT:

        A vector.

        """
        from sage.arith.srange import srange
        from sage.modules.free_module_element import vector

        base = self.k
        M = recursion_rules.M
        m = recursion_rules.m
        ll = recursion_rules.ll
        uu = recursion_rules.uu
        initial_values = recursion_rules.initial_values
        n0 = recursion_rules.n0
        right = []

        for j in srange(m):
            for d in srange(base**j):
                try:
                    right.append(initial_values[d])
                except KeyError:
                    raise ValueError('Initial value %s is missing.'
                                     % (function(d),))

        for j in srange(m, M):
            for d in srange(ll, base**j - base**m + uu + 1):
                try:
                    right.append(initial_values[d])
                except KeyError:
                    raise ValueError('Initial value %s is missing.'
                                     % (function(d),))

        if n0 >= 1:
            right = right + [1] + (n0 - 1)*[0]

        return vector(right)


    def recursions(self, equations, function, var, n0=0, minimize=False):
        r"""
        Construct a `k`-regular sequence that fulfills the recursions
        given in ``equations``.

        INPUT:

        - ``equations`` -- A list of equations where the elements have
          either the form

          - ``f(k^M * n + r) == sum(f(k^m * n + k) for k in srange(l, u + 1)``
            for some integers ``0 <= r < k^M`` and ``M > m >= 0``
            and some ``l <= u`` -- valid for all non-negative
            integers ``n`` -- and there is an equation of this form
            for all ``r``

          or the form

          - ``f(k) == t`` for some integer ``k`` and some ``t``.

        - ``function`` -- symbolic function ``f`` occuring in the equations.

        - ``var`` -- symbolic variable ``n`` occuring in the equations.

        - ``minimize`` -- (default: ``False``) a boolean. If ``True``,
          the linear representation of the resulting sequence is
          minimized (see :meth:`~minimized`).

        OUTPUT:

        A :class:`kRegularSequence`.

        EXAMPLES:

        Stern--Brocot Sequence::

            sage: Seq2 = kRegularSequenceSpace(2, ZZ)
            sage: var('n')
            n
            sage: function('f')
            f
            sage: Seq2.recursions([
            ....:     f(2*n) == f(n), f(2*n + 1) == f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 1], f, n)
            2-regular sequence 0, 1, 1, 2, 1, 3, 2, 3, 1, 4, ...

        Number of Odd Entries in Pascal's Triangle::

            sage: Seq2.recursions([
            ....:     f(2*n) == 3*f(n), f(2*n + 1) == 2*f(n) + f(n + 1),
            ....:     f(0) == 0, f(1) == 1, f(2) == 3, f(3) == 5], f, n)
            2-regular sequence 0, 1, 3, 5, 9, 11, 15, 19, 27, 29, ...

        Number of Unbordered Factors in the Thue--Morse Sequence::

            sage: Seq2.recursions([
            ....:     f(8*n) == 2*f(4*n), f(8*n + 1) == f(4*n + 1),
            ....:     f(8*n + 2) == f(4*n + 1) + f(4*n + 3),
            ....:     f(8*n + 3) == -f(4*n + 1) + f(4*n + 2),
            ....:     f(8*n + 4) == 2*f(4*n + 2), f(8*n + 5) == f(4*n + 3),
            ....:     f(8*n + 6) == -f(4*n + 1) + f(4*n + 2) + f(4*n + 3),
            ....:     f(8*n + 7) == 2*f(4*n + 1) + f(4*n + 3),
            ....:     f(0) == 1, f(1) == 2, f(2) == 2, f(3) == 4, f(4) == 2,
            ....:     f(5) == 4, f(6) == 6, f(7) == 0, f(8) == 4, f(9) == 4,
            ....:     f(10) == 4, f(11) == 4, f(12) == 12, f(13) == 0, f(14)==4,
            ....:     f(15) == 4, f(16) == 8, f(17) == 4, f(18) == 8, f(19) == 0,
            ....:     f(20) == 8, f(21) == 4, f(22) == 4, f(23) == 8, f(24) == 24,
            ....:     f(25) == 0, f(26) == 4, f(27) == 4, f(28) == 8, f(29) == 4,
            ....:     f(30) == 8, f(31) == 4, f(32) == 16, f(33) == 4], f, n, 3)
            2-regular sequence 1, 2, 2, 4, 2, 4, 6, 0, 4, 4, ...

        Number of Non-Zero Elements in the a Generalized Pascal's Triangle [TODO: reference]::

            sage: Seq2 = kRegularSequenceSpace(2, QQ)
            sage: Seq2.recursions([
            ....:     f(4*n) == 5/3*f(2*n) - 1/3*f(2*n + 1),
            ....:     f(4*n + 1) == 4/3*f(2*n) + 1/3*f(2*n + 1),
            ....:     f(4*n + 2) == 1/3*f(2*n) + 4/3*f(2*n + 1),
            ....:     f(4*n + 3) == -1/3*f(2*n) + 5/3*f(2*n + 1),
            ....:     f(0) == 1, f(1) == 2, f(2) == 3, f(3) == 3], f, n)
            2-regular sequence 1, 2, 3, 3, 4, 5, 5, 4, 5, 7, ...
        """
        from sage.arith.srange import srange

        k = self.k
        mu = []

        recursion_rules = self._parse_recursions_(equations, function, var, n0)

        for rem in srange(k):
            mu.append(self._get_matrix_from_recursions_(recursion_rules, rem, function, var))

        seq = self(mu, self._get_left_from_recursions_(recursion_rules.dim + n0),
                   self._get_right_from_recursions_(recursion_rules, function))

        if minimize:
            return seq.minimized()
        else:
            return seq
