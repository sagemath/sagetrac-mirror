# -*- coding: utf-8 -*-
r"""
Piecewise-defined Functions

This module implement piecewise polynomial in a single variable. See
:mod:`sage.sets.real_set` for more information about how to construct
subsets of the real line for the domains.

EXAMPLES::

    sage: R.<t> = QQ[]
    sage: f1 = 1 - t
    sage: f2 = t^4 - t^2
    sage: D1 = RealSet([0, 1], [4, 5])
    sage: D2 = RealSet([2, 3], [6, 7])
    sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
    PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
    sage: p*p
    PiecewisePolynomial(t |--> t^2 - 2*t + 1 on [0, 1] ∪ [4, 5], t |--> t^8 - 2*t^6 + t^4 on [2, 3] ∪ [6, 7]; t)


AUTHORS:

- Yueqi Li, Yuan Zhou (2022-09): initial version
"""

# ****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 David Joyner <wdjoyner@gmail.com>
#                     2013 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.rings.polynomial.polynomial_element import Polynomial, is_Polynomial
from sage.sets.real_set import InternalRealInterval, RealSet
from heapq import merge
import bisect
from sage.structure.element import Element, ModuleElement


class Piecewise_Polynomial(ModuleElement):
    def __init__(self, function_pieces):
        """
        Piecewise polynomial

        INPUT:

        - ``function_pieces`` -- a list of pairs consisting of a
          domain and a polynomial function.

        OUTPUT:

        A piecewise-defined polynomial. A ``ValueError`` will be raised
        if the domains of the pieces are not pairwise disjoint or the input
        function is not polynomial

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)

        TESTS::
            sage: D3 = RealSet([1, 2])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2], [D3, f2]]); p
            Traceback (most recent call last):
            ...
            ValueError: domains must be pairwise disjoint
            sage: f3=1+x
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2], [D3, f3]]); p
            Traceback (most recent call last):
            ...
            ValueError: Invalid function. Function type must be Polynomial, got <class 'sage.symbolic.expression.Expression'> for x + 1.


        """
        self.domain_list = []
        self.func_list = []
        self.var = None
        input_func_type = None
        p = {}
        for i, (domain, func) in enumerate(function_pieces):
            if not isinstance(domain, RealSet):
                try:
                    domain = RealSet(domain)
                except:
                    raise ValueError(
                        "Invalid domain. Domain must be RealSet, but got {domain_type} for {domain_value}.".format(
                            domain_type=type(domain), domain_value=domain))
            if domain.is_empty():
                continue
            if not is_Polynomial(func):
                raise ValueError("Invalid function. Function type must be Polynomial, got {func_type} for {func_value}.".format(
                    func_type=type(func), func_value=func))
            if input_func_type is None:
                input_func_type = type(func)
            elif type(func) != input_func_type:
                raise ValueError("Inconsistant function types. Should be {first_type}, got {curr_type}".format(
                    first_type=input_func_type, curr_type=type(func)))
            if self.var is None:
                self.var = func.args()[0]
            elif func.args()[0] != self.var:
                raise ValueError("Inconsistant variable. Should be {0}, got {1}".format(self.var, func.args()[0]))

            non_point = []
            for interval in domain:
                if interval.is_point():
                    value = func(interval._lower)
                    # this value could not be Integer or Rational
                    # this could be LazyWrapper
                    if hasattr(value, "_value"):
                        value = value._value
                    func_ = value + 0 * func
                    if func_ in p:
                        i = p[func_]
                        self.domain_list[i] = self.domain_list[i].union(interval)
                    else:
                        self.domain_list.append(RealSet(interval))
                        self.func_list.append(func_)
                        p[func_] = len(self.domain_list) - 1
                else:
                    non_point.append(interval)

            if non_point:
                if func in p:
                    i = p[func]
                    self.domain_list[i] = RealSet().union(*([self.domain_list[i]] + non_point))
                else:
                    self.domain_list.append(RealSet().union(*non_point))
                    self.func_list.append(func)
                    p[func] = len(self.domain_list) - 1

        del p
        if not RealSet.are_pairwise_disjoint(*self.domain_list):
            raise ValueError("domains must be pairwise disjoint")

        self.support = RealSet().union(*self.domain_list)
        self._end_points = None
        self._end_points_list = None

    def _init_end_points(self):
        """
         Compute function value for all end point,  ` end_point^+`, function value at `end_point^-` . It's part of
         initialization.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
            sage: p._init_end_points()
            sage: p._end_points
            {0: [1, [None, 1]],
             1: [0, [0, None]],
             2: [12, [None, 12]],
             3: [72, [72, None]],
             4: [-3, [None, -3]],
             5: [-4, [-4, None]],
             6: [1260, [None, 1260]],
             7: [2352, [2352, None]]}
        """
        self._end_points = {}
        # [(x, epsilon), delta, i]
        self._end_points_list = [self._iterator(real_set, i) for i, real_set in enumerate(self.domain_list)]
        self._end_points_list = list(merge(*self._end_points_list))
        for i, (dom, func) in enumerate(zip(self.domain_list, self.func_list)):
            for (x, epsilon), delta in dom._scan():
                if x not in self._end_points: self._end_points[x] = [None, [None, None]]
                func_val = func(x)
                # this delta is from RealSet._scan()
                if delta < 0:
                    self._end_points[x][1][1] = func_val
                    if epsilon == 0:
                        self._end_points[x][0] = func_val
                else:
                    self._end_points[x][1][0] = func_val
                    if epsilon == 1:
                        self._end_points[x][0] = func_val

    def __repr__(self):
        """
        Return a string representation of piecewise polynomial

        OUTPUT:

        String.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: p = PiecewisePolynomial([[D1, f1]])
            sage: str(p)
            'PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5]; t)'
        """
        s = "PiecewisePolynomial("
        args = []
        for dom, func in zip(self.domain_list, self.func_list):
            args.append("{0} |--> {1} on {2}".format(self.var, func, dom))
        s += ", ".join(args) + '; {0})'.format(self.var)
        return s

    def __call__(self, *args, **kwds):
        r"""
        Piecewise polynomial

        OUTPUT:

        A piecewise-defined polynomial. A ``ValueError`` will be raised
        if more than 1 input

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: my_abs = PiecewisePolynomial([((-1, 0), -t), ([0, 1], t)]);  my_abs
            PiecewisePolynomial(t |--> -t on (-1, 0), t |--> t on [0, 1]; t)
            sage: [ my_abs(i/5) for i in range(-4, 5)]
            [4/5, 3/5, 2/5, 1/5, 0, 1/5, 2/5, 3/5, 4/5]
        """
        val = None
        if len(args) == 1:
            val = args[0]
            if len(kwds) > 0:
                raise ValueError("Invalid input: Got more than 1 inputs")
        elif len(args) == 0:
            if str(self.var) in kwds:
                val = kwds[str(self.var)]
                if len(kwds) > 1:
                    print("Warning: more than one input detected, only {0}={1} used.".format(str(self.var), val))
            else:
                raise ValueError("Invalid input: No positional input for {0} detected.".format(str(self.var)))
        else:
            raise ValueError("Invalid input: More than 1 arguments detected.")

        for dom, func in zip(self.domain_list, self.func_list):
            if val in dom:
                return func(val)
        print("Warning: {0}={1} does not fall in any domain defined. Returned None type.".format(str(self.var), val))
        return None

    def __len__(self):
        """
        Return the number of "pieces"

        OUTPUT:

        Integer

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]);
            sage: len(p)
            2
        """
        return len(self.func_list)

    def __iter__(self):
        """
        iter over piecewise polynomials

        OUTPUT:

        domains and functions

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t ^ 4 - t ^ 2
            sage: D1 = RealSet((-oo, 1), (4, 5))
            sage: D2 = RealSet([2, 3], x >= 6)
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]])
            sage: for dom, func in p:
            ....:     print(dom, func)
            ....:
            (-oo, 1) ∪ (4, 5) -t + 1
            [2, 3] ∪ [6, +oo) t^4 - t^2
        """
        for dom, func in zip(self.domain_list, self.func_list):
            yield dom, func

    def __add__(self, other):
        """
        Add two identical domain Piecewise polynomials

        INPUT:

        - ``other`` -- Another piecewise polynomial

        OUTPUT:

        A piecewise-defined polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: f3 = t^2
            sage: f4 = 1-t^7
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]])
            sage: q = PiecewisePolynomial([[D3, f3], [D4, f4]])
            sage: p+q
            PiecewisePolynomial(t |--> t^2 - t + 1 on [0, 1] ∪ [4, 5], t |--> -t^7 + t^4 - t^2 + 1 on [2, 3] ∪ [6, 7]; t)
        """
        if self.var != other.var:
            print("Invalid variables: Cannot add variable {0} with {1}".format(self.var, other.var))
            return None
        if type(other) == type(self):
            return self.PiecewisePolynomial_add([self, other])
        else:
            return Piecewise_Polynomial((dom, func + other) for dom, func in self.__iter__())

    __radd__ = __add__

    def __eq__(self, other):
        """
        Return if two piecewise polynomials are equal

        INPUT:

        - ``other`` -- Another piecewise polynomials

        OUTPUT:

        Boolean

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = t
            sage: f2 = 2 - t
            sage: f3 = t
            sage: f4 = 2 - t
            sage: D1 = RealSet([0, 1])
            sage: D2 = RealSet((1, 2))
            sage: D3 = RealSet(RealSet.closed_open(0, 1))
            sage: D4 = RealSet(RealSet.closed_open(1, 2))
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> t on [0, 1], t |--> -t + 2 on (1, 2); t)
            sage: q = PiecewisePolynomial([[D3, f3], [D4, f4]]); q
            PiecewisePolynomial(t |--> t on [0, 1), t |--> -t + 2 on [1, 2); t)
            sage: p==q
            True
        """

        if self.support != other.support or self.var != other.var or any(
                func != 0 * self.func_list[0] for _, func in (self - other)):
            return False
        return True

    def __neg__(self):
        """
        Return the negative of polynomials

        OUTPUT:

        A piecewise-defined polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = t
            sage: f2 = 2 - t
            sage: D1 = RealSet([0, 1])
            sage: D2 = RealSet((1, 2))
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> t on [0, 1], t |--> -t + 2 on (1, 2); t)
            sage: -p
            PiecewisePolynomial(t |--> -t on [0, 1], t |--> t - 2 on (1, 2); t)
        """

        return Piecewise_Polynomial((dom, -func) for dom, func in self.__iter__())

    def __sub__(self, other):
        """
        Subtraction another identical domain's piecewise polynomials

        INPUT:

        - ``other`` -- Another piecewise polynomial

        OUTPUT:

        A piecewise-defined polynomial

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: f3 = t^2
            sage: f4 = 1-t^7
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
            sage: q = PiecewisePolynomial([[D3, f3], [D4, f4]]); q
            PiecewisePolynomial(t |--> t^2 on [0, 1] ∪ [4, 5], t |--> -t^7 + 1 on [2, 3] ∪ [6, 7]; t)
            sage: p-q
            PiecewisePolynomial(t |--> -t^2 - t + 1 on [0, 1] ∪ [4, 5], t |--> t^7 + t^4 - t^2 - 1 on [2, 3] ∪ [6, 7]; t)
        """
        return self.__add__(other.__neg__())

    def __rsub__(self, other):
        return other.__add__(self.__neg__())

    def __mul__(self, other):
        r"""
        multiply two identical domain piecewise polynomial

        INPUT:

        - ``other`` -- Another piecewise polynomial

        OUTPUT:

        A piecewise-defined polynomial.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: f3 = t^2
            sage: f4 = 1-t^7
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
            sage: q = PiecewisePolynomial([[D3, f3], [D4, f4]]); q
            PiecewisePolynomial(t |--> t^2 on [0, 1] ∪ [4, 5], t |--> -t^7 + 1 on [2, 3] ∪ [6, 7]; t)
            sage: p*q
            PiecewisePolynomial(t |--> -t^3 + t^2 on [0, 1] ∪ [4, 5], t |--> -t^11 + t^9 + t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
        """
        if self.var != other.var:
            print("Invalid variables: Cannot add variable {0} with {1}".format(self.var, other.var))
            return None
        if type(other) == type(self):
            return self.PiecewisePolynomial_mul([self, other])
        else:
            return Piecewise_Polynomial((dom, func * other) for dom, func in self.__iter__())

    __rmul__ = __mul__

    # Other should be a number
    def __truediv__(self, other):
        """
        Divide by a scalar.

        INPUT:

        - ``other``: scalar

        OUTPUT:

        A piecewise-defined polynomial

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: f3 = t^2
            sage: f4 = 1-t^7
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
            sage: p/2
            PiecewisePolynomial(t |--> -1/2*t + 1/2 on [0, 1] ∪ [4, 5], t |--> 1/2*t^4 - 1/2*t^2 on [2, 3] ∪ [6, 7]; t)
        """
        return Piecewise_Polynomial((dom, func / other) for dom, func in self.__iter__())

    def __pow__(self, power):

        """
        Power of piecewise polynomial

        INPUT:

        - ``power``: scalar

        OUTPUT:

        A piecewise-defined polynomial

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t^4 - t^2
            sage: f3 = t^2
            sage: f4 = 1-t^7
            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on [0, 1] ∪ [4, 5], t |--> t^4 - t^2 on [2, 3] ∪ [6, 7]; t)
            sage: p^2
            PiecewisePolynomial(t |--> t^2 - 2*t + 1 on [0, 1] ∪ [4, 5], t |--> t^8 - 2*t^6 + t^4 on [2, 3] ∪ [6, 7]; t)
        """
        return Piecewise_Polynomial((dom, func ** power) for dom, func in self.__iter__())

    @staticmethod
    def _iterator(realset, i):
        """
        scan function for scan-line

        INPUT:

        - ``realset`` -- A collection of realsets
        - ``i`` -- integer

        OUTPUT:

        -  Generate events of the form ``((x, epsilon), delta), i``

        When ``x`` is the beginning of an interval ('on'):

        - ``epsilon`` is 0 if the interval is lower closed and 1 otherwise,
        - ``delta`` is 1

        When ``x`` is the end of an interval ('off'):

        - ``epsilon`` is 1 if the interval is upper closed and 0 otherwise,
        - ``delta`` is -1

        EXAMPLES::

            sage: D1 = RealSet([0, 2], [3, 5])
            sage: list(PiecewisePolynomial._iterator(D1, None))
            [(((0, 0), 1), None),
             (((2, 1), -1), None),
             (((3, 0), 1), None),
             (((5, 1), -1), None)]

        """
        for interval in realset:
            yield ((interval._lower, 1 - int(interval._lower_closed)), 1), i
            yield ((interval._upper, int(interval._upper_closed)), -1), i

    @staticmethod
    def finest_partitions(*real_set_collection):
        """
        Given a set of realsets, this function computes the set of intervals that could represent each realset by
        taking disjoint union of a subset of such intervals.

        INPUT:

        - ``*real_set_collection`` -- a list/tuple/iterable of :class:`RealSet`
          or data that defines one.

        OUTPUT:

        iterator of realsets and indicator that which input realsest contain the output realset

        EXAMPLES::

            sage: D1 = RealSet([0, 1], [4, 5])
            sage: D2 = RealSet([2, 3], [6, 7])
            sage: D3 = RealSet([0, 1], [4, 5])
            sage: D4 = RealSet([2, 3], [6, 7])
            sage: list( PiecewisePolynomial.finest_partitions(D1,D2,D3,D4))
            [([0, 1], (1, 0, 1, 0)),
             ([2, 3], (0, 1, 0, 1)),
             ([4, 5], (1, 0, 1, 0)),
             ([6, 7], (0, 1, 0, 1))]

        TESTS::

            sage: D1 = RealSet((0, 1))
            sage: D2 = RealSet([-1, 1])
            sage: D3 = RealSet((1, 2))
            sage: D4 = RealSet([1, 3])
            sage: list( PiecewisePolynomial.finest_partitions(D1,D2,D3,D4))
            [([-1, 0], (0, 1, 0, 0)),
             ((0, 1), (1, 1, 0, 0)),
             ({1}, (0, 1, 0, 1)),
             ((1, 2), (0, 0, 1, 1)),
             ([2, 3], (0, 0, 0, 1))]
        """
        scan = [Piecewise_Polynomial._iterator(real_set, i) for i, real_set in enumerate(real_set_collection)]
        indicator = [0] * len(scan)
        indicator_sum = 0
        scan = merge(*scan)
        prev_event = None

        for event, set_id in scan:
            (x, epsilon), delta = event
            if prev_event and event > prev_event:
                (x_, epsilon_), _ = prev_event
                yield RealSet(InternalRealInterval(x_, epsilon_ == 0, x, epsilon == 1)), tuple(indicator)
            indicator[set_id] += delta
            indicator_sum += delta

            prev_event = (x, epsilon), 1
            if indicator_sum == 0:
                prev_event = None

    @staticmethod
    def PiecewisePolynomial_add(PiecewisePolynomial_collection, union=None):
        """
        Compute the sum of PiecewisePolynomial collections

        INPUT:

        - ``PiecewisePolynomial_collection`` -- a list of PiecewisePolynomial polynomials
        - ``union``: If not given, then add same domains PiecewisePolynomial polynomial;
                    if True, take union of those PiecewisePolynomial collection, then add;
                    if False, take intersection of PiecewisePolynomial collection, then add
        OUTPUT:

        Sum of the input PiecewisePolynomial polynomials

        EXAMPLES::

           sage: R.<t> = QQ[]
           sage: p1 = PiecewisePolynomial([[[0, 2], t], [[3, 5], t^2], [[6, 7], 1 - t]])
           sage: p2 = PiecewisePolynomial([[[0, 2], t+1], [[5, 6], t], [[7, 7], -t]])
           sage: p3 = PiecewisePolynomial([[RealSet.closed_open(1, 2), (t+1)^2]])
           sage: PiecewisePolynomial.PiecewisePolynomial_add([p1,p2, p3],True)
           PiecewisePolynomial(t |--> 2*t + 1 on [0, 1), t |--> t^2 + 4*t + 2 on [1, 2), t |--> 5 on {2}, t |--> t^2 on [3, 5), t |--> 30 on {5}, t |--> t on (5, 6), t |--> 1 on {6}, t |--> -t + 1 on (6, 7), t |--> -13 on {7}; t)
           sage: PiecewisePolynomial.PiecewisePolynomial_add([p1,p2,p3],False)
           PiecewisePolynomial(t |--> t^2 + 4*t + 2 on [1, 2); t)
           sage: PiecewisePolynomial.PiecewisePolynomial_add([p1,p2,p3])
           Traceback (most recent call last):
           ...
           ValueError: Inconsistent domains. Please set union=True for union, or False for intersection.
        """

        realset_info = []
        for i, PiecewisePolynomial in enumerate(PiecewisePolynomial_collection):
            for j in range(len(PiecewisePolynomial)):
                realset_info.append((i, j))

        result_pairs = []
        p = Piecewise_Polynomial.finest_partitions(*[dom for PiecewisePolynomial in PiecewisePolynomial_collection for dom, _ in PiecewisePolynomial])
        for real_set, inds in p:
            new_pair = [real_set, 0]
            if sum(inds) < len(PiecewisePolynomial_collection):
                if union is None:
                    raise ValueError(
                        "Inconsistent domains. Please set union=True for union, or False for intersection.")
                elif union is False:
                    continue

            # i - index for realset, ind: True or False
            for i, ind in enumerate(inds):
                if ind > 0:
                    realset_id, func_id = realset_info[i]
                    new_pair[1] += PiecewisePolynomial_collection[realset_id].func_list[func_id]
            result_pairs.append(new_pair)

        return Piecewise_Polynomial(result_pairs)

    @staticmethod
    def PiecewisePolynomial_mul(PiecewisePolynomial_collection, union=None):
        """
        Compute the multiplication of PiecewisePolynomial collections

        INPUT:

        - ``PiecewisePolynomial_collection`` -- a list of PiecewisePolynomial polynomials
        - ``union``: If not given, then multiply same domains PiecewisePolynomial polynomial;
                    if True, take union of those PiecewisePolynomial collection, then multiply;
                    if False, take intersection of PiecewisePolynomial collection, then multiply;
        OUTPUT:

        Multiplication of the input PiecewisePolynomial polynomials

        EXAMPLES::

           sage: R.<t> = QQ[]
           sage: p1 = PiecewisePolynomial([[[0, 2], t], [[3, 5], t^2], [[6, 7], 1 - t]])
           sage: p2 = PiecewisePolynomial([[[0, 2], t+1], [[5, 6], t], [[7, 7], -t]])
           sage: p3 = PiecewisePolynomial([[RealSet.closed_open(1, 2), (t+1)^2]])
           sage: PiecewisePolynomial.PiecewisePolynomial_mul([p1,p2, p3],True)
           PiecewisePolynomial(t |--> t^2 + t on [0, 1), t |--> t^4 + 3*t^3 + 3*t^2 + t on [1, 2), t |--> 6 on {2}, t |--> t^2 on [3, 5), t |--> 125 on {5}, t |--> t on (5, 6), t |--> -30 on {6}, t |--> -t + 1 on (6, 7), t |--> 42 on {7}; t)
           sage: PiecewisePolynomial.PiecewisePolynomial_mul([p1,p2,p3],False)
           PiecewisePolynomial(t |--> t^4 + 3*t^3 + 3*t^2 + t on [1, 2); t)
           sage: PiecewisePolynomial.PiecewisePolynomial_mul([p1,p2,p3])
           Traceback (most recent call last):
           ...
           ValueError: Inconsistent domains. Please set union=True for union, or False for intersection.
        """

        realset_info = []
        for i, piecewise in enumerate(PiecewisePolynomial_collection):
            for j in range(len(piecewise)):
                realset_info.append((i, j))

        result_pairs = []
        p = Piecewise_Polynomial.finest_partitions(*[dom for piecewise in PiecewisePolynomial_collection for dom, _ in piecewise])
        for real_set, inds in p:
            new_pair = [real_set, 1]
            if sum(inds) < len(PiecewisePolynomial_collection):
                if union is None:
                    raise ValueError(
                        "Inconsistent domains. Please set union=True for union, or False for intersection.")
                elif union is False:
                    continue

            # i - index for realset, ind: True or False
            for i, ind in enumerate(inds):
                if ind > 0:
                    realset_id, func_id = realset_info[i]
                    new_pair[1] *= PiecewisePolynomial_collection[realset_id].func_list[func_id]
            result_pairs.append(new_pair)

        return Piecewise_Polynomial(result_pairs)

    def is_continuous(self):
        """
        Return if the function is continuous

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = t
            sage: f2 = (t-1) ^ 2 + 1
            sage: f3 = 1 - t
            sage: D1 = RealSet([0, 1], [2, 3])
            sage: D2 = RealSet((1, 2))
            sage: D3 = RealSet((3, +oo))
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2], [D3, f3]])
            sage: q = PiecewisePolynomial([[D1, f1], [D2, f2]])
            sage: p.is_continuous()
            False
            sage: q.is_continuous()
            True
        """
        if not self._end_points:
            self._init_end_points()

        for point in self._end_points:
            val, (left, right) = self._end_points[point]
            if val:
                if (left is not None and left != val) or (right is not None and right != val):
                    return False
        return True

    def is_continuous_defined(self, xmin=0, xmax=1):
        r"""
        Return if the function is continuous at certain interval or point

        INPUT:

        - ``xmin`` -- the lower bound of interval
        - ``xmax `` -- the upper bound of interval

        OUTPUT:

        Boolean.

        EXAMPLES::
            sage: R.<t> = QQ[]
            sage: f1 = t
            sage: f2 = (t - 1) ^ 2 + 1
            sage: f3 = 1 - t
            sage: D1 = RealSet([0, 1], RealSet.closed_open(2, 3))
            sage: D2 = RealSet((1, 2))
            sage: D3 = RealSet(x >= 3)
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2], [D3, f3]])
            sage: p.is_continuous_defined(0, 2)
            True
            sage: p.is_continuous_defined(1, 3)
            False
            sage: p.is_continuous_defined(2, 4)
            False
            sage: p.is_continuous_defined(3, 4)
            True
        """

        if not self._end_points_list:
            self._init_end_points()
        # Should discuss whether it is reasonable to force
        # [xmin, xmax] totally included in domain
        if not RealSet([xmin, xmax]).is_subset(self.support):
            return False

        # O(log n) implementation
        ind = bisect.bisect_left(self._end_points_list, (((xmin, 0), -1), 0))
        seen = set()
        for i in range(ind, len(self._end_points_list)):
            ((point, epsilon), delta), idx = self._end_points_list[i]
            if point in seen:
                continue
            seen.add(point)
            val, (left, right) = self._end_points[point]
            if point == xmin and val != right:
                return False
            elif point == xmax and val != left:
                return False
            elif xmin < point < xmax and (val != left or val != right):
                return False
            elif point > xmax:
                return True
        return True


    def which_pair(self, x0):
        """
        Find Input x0 in which function

        INPUT:

        - ``x0``: a number in domain of piecewise polynomial

        Returns:

        A piece in the piecewise polynomial

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t ^ 4 - t ^ 2
            sage: D1 = RealSet((-oo, 1), (4, 5))
            sage: D2 = RealSet([2, 3], x >= 6)
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]])
            sage: p.which_pair(0)
            ((-oo, 1) ∪ (4, 5), -t + 1)
            sage: p.which_pair(1)
            Invalid input: x0 not in domain
            sage: p.which_pair(2)
            ([2, 3] ∪ [6, +oo), t^4 - t^2)
            sage: p.which_pair(3)
            ([2, 3] ∪ [6, +oo), t^4 - t^2)
            sage: p.which_pair(4)
            Invalid input: x0 not in domain
        """

        if not self._end_points:
            self._init_end_points()
        idx = bisect.bisect_left(self._end_points_list, (((x0, 1), -1), 0))
        if idx <= 0 or idx >= len(self._end_points_list):
            print("Invalid input: x0 not in domain")
            return None
        ((x, epsilon), delta), func_idx = self._end_points_list[idx]
        if delta < 0:
            return self.domain_list[func_idx], self.func_list[func_idx]
        else:
            print("Invalid input: x0 not in domain")
            return None

    def limits(self, x0):
        """
        Returns [function value at `x_0`, function value at `x_0^+`, function value at `x_0^-`].

        INPUT:

        -``x0``: A number in domain of piecewise polynomial

        Returns:

        function value at `x_0`, function value at `x_0^+`, function value at `x_0^-`.

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t ^ 4 + 1
            sage: D1 = RealSet((-oo, 1), (5, 6))
            sage: D2 = RealSet([1, 2], x >= 7)
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]])
            sage: p.limits(1)
            (2, 0, 2)
            sage: p.limits(2)
            (17, 17, None)
            sage: p.limits(3)
            Invalid input: x0 not in domain
            (None, None, None)
            sage: p.limits(5)
            (None, None, -4)
            sage: p.limits(6)
            (None, -5, None)
            sage: p.limits(7)
            (2402, None, 2402)
            sage: p.limits(8)
            (4097, 4097, 4097)
        """

        if not self._end_points:
            self._init_end_points()

        if x0 in self._end_points:
            return self._end_points[x0][0], self._end_points[x0][1][0], self._end_points[x0][1][1]
        else:
            result = self.which_pair(x0)
            if result:
                _, func = result
                func_val = func(x0)
                return func_val, func_val, func_val
            else:
                return None, None, None

    def derivative(self, var=None):
        """
        Derivative of piecewise polynomial

        INPUT:

        - ``var``: the variable

        Returns:

        Piecewise function

        EXAMPLES::

            sage: R.<t> = QQ[]
            sage: f1 = 1 - t
            sage: f2 = t ^ 4 + 1
            sage: D1 = RealSet((-oo, 1), (5, 6))
            sage: D2 = RealSet([1, 2], x >= 7)
            sage: p = PiecewisePolynomial([[D1, f1], [D2, f2]]); p
            PiecewisePolynomial(t |--> -t + 1 on (-oo, 1) ∪ (5, 6), t |--> t^4 + 1 on [1, 2] ∪ [7, +oo); t)
            sage: p.derivative()
            PiecewisePolynomial(t |--> -1 on (-oo, 1) ∪ (5, 6), t |--> 4*t^3 on [1, 2] ∪ [7, +oo); t)
        """

        return Piecewise_Polynomial((dom, func._derivative(var)) for dom, func in self.__iter__())


PiecewisePolynomial = Piecewise_Polynomial
