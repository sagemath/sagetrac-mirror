# -*- coding: utf-8 -*-
r"""
Basic Linear Algebra Subroutines on dictionaries

This module provides functions for "level 1 basic linear algebra
subroutines" on sparse vectors represented as dictionaries. The API is
inspired from the standard BLAS API
:wikipedia:`Basic_Linear_Algebra_Subprograms`, but does not try to
follow it exactly.

For ``a, b, ...`` hashable objects, the dictionary ``{a: 2, b:3, ...}``
represents the formal linear combination `2 \cdot a + 3 \cdot b + \cdots`.
The values of the dictionary should all lie in the same parent `K`. For
simplicity, we call this formal linear combination a vector, as the
typical use case is `K` being a field. However the routines here can
be used with a ring `K` to represent free module elements, or a
semiring like `\NN` to represent elements of a commutative monoid, or
even just an additive semigroup. Of course not all operations are
meaningful in those cases. We are also assuming that ``-1 * x = -x``
and ``bool(x) == bool(-x)`` for all ``x`` in `K`.

Unless stated overwise, all values `v` in the dictionaries should be
non zero (as tested with `bool(v)`).

This is mostly used by :class:`CombinatorialFreeModule`.
"""
#*****************************************************************************
#       Copyright (C) 2010 Christian Stump <christian.stump@univie.ac.at>
#                     2016 Travis Scrimshaw <tscrimsh@umn.edu>
#                     2016 Nicolas M. Thiéry <nthiery at users.sf.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

cpdef int iaxpy(a, dict X, dict Y, bint remove_zeros=True, bint factor_on_left=True) except -1:
    r"""
    Mutate `Y` to represent `a X + Y`.

    INPUT:

    - ``a`` -- element of a parent `K` or `±1`
    - ``X,Y`` -- dictionaries representing a vector `X` over `K`
    - ``remove_zeros`` -- boolean (default: ``True``); whether to
      remove the keys whose values are zero after the addition has
      been performed
    - ``factor_on_left`` -- boolean (default: ``False``);
      whether to compute `a X + Y` or `X a + Y`

    OUTPUT:

    ``0`` when successful and ``-1`` on error. This is done because
    returning an ``int`` is faster than a (Python) ``None``.

    ``Y`` has been mutated to represent the vector `a \cdot X + Y`.
    If ``remove_zeros=True`` (and the input have no zero values
    themselves), then the output is guaranteed to have no zero
    values. Otherwise some zero values may be left in the output.
    Which ones is voluntarily left undefined. Use this in cases
    where you want to postpone clearing until the end of a long
    series of operations.

    The parent `K` should support addition unless `a = -1`, negation
    if `a = -1`, and multiplication if `a \neq ±1`. We are also
    assuming that ``-1 * x = -x`` and ``bool(x) == bool(-x)`` in `K`.

    See :mod:`sage.data_structures.blas_dict` for an overview.

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: X = {0: -1, 1: 1}
        sage: Y = {0: 1, 1: 1}

    Computing `X+Y`::

        sage: blas.iaxpy(1, X, Y)
        0
        sage: Y
        {1: 2}

    Reseting `Y` and computing `-X+Y`::

        sage: Y = {0: 1, 1: 1}
        sage: blas.iaxpy(-1, X, Y)
        0
        sage: Y
        {0: 2}

    Reseting `Y` and computing `2X+Y`::

        sage: Y = {0: 1, 1: 1}
        sage: blas.iaxpy(2, X, Y)
        0
        sage: Y
        {0: -1, 1: 3}

    Reseting `Y` and computing `-X+Y` without removing zeros::

        sage: Y = {0: 1, 1: 1}
        sage: blas.iaxpy(-1, X, Y, remove_zeros=False)
        0
        sage: Y
        {0: 2, 1: 0}
    """
    cdef int flag = 0
    if a == 1:
        flag = 1
    elif a == -1:
        flag = -1
    elif not a:
        return 0
    for (key, value) in X.iteritems():
        if flag == -1:
            if key in Y:
                Y[key]  -= value
            else:
                Y[key]  = -value
                continue # no need to check for zero
        else:
            if flag != 1:
                # If we had the guarantee that `a` is in the base
                # ring, we could use a._mul_(value), as suggested in
                # the documentation of sage.structure.element. However
                # we want to support elements that can act on values
                # in the base ring as well
                if factor_on_left:
                    value = a*value
                else:
                    value = value*a
                if not value:
                    continue # a is a zero divizor
            # assert value
            if key in Y:
                Y[key] += value
            else:
                Y[key]  = value
                continue # no need to check for zero
        if remove_zeros and not Y[key]:
            del Y[key]
    return 0

cpdef dict axpy(a, dict X, dict Y, bint factor_on_left=True):
    """
    Return `a X + Y`.

    All input dictionaries are supposed to have values in the same
    base ring `K` and have no nonzero values.

    INPUT:

    - ``a`` -- an element of `K` or `±1`
    - ``X``, ``Y`` -- dictionaries representing respectively vectors `X` and `Y`
    - ``factor_on_left`` -- boolean (default: ``False``);
      whether to compute `a X + Y` or `X a + Y`

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: Y = {0: 1, 1: 1}
        sage: X = {0: -1, 1: 1}

    Computing `X + Y`::

        sage: blas.axpy(1, X, Y)
        {1: 2}

    `X` and `Y` are not mutated:

        sage: X, Y
        ({0: -1, 1: 1}, {0: 1, 1: 1})

    Computing `-X + Y`::

        sage: blas.axpy(-1, X, Y)
        {0: 2}

    Computing `2X + Y`::

        sage: blas.axpy(2, X, Y)
        {0: -1, 1: 3}

    Testing the special cases of empty dictionaries::

        sage: Y = {}
        sage: Z = blas.axpy(2, X, Y)
        sage: Z is Y
        False
        sage: X, Y
        ({0: -1, 1: 1}, {})

        sage: Z = blas.axpy(2, Y, X)
        sage: Z is X
        True
        sage: X, Y
        ({0: -1, 1: 1}, {})
    """
    if X:
        Y = Y.copy()
        iaxpy(a, X, Y, True, factor_on_left)
    return Y

cpdef dict negate(dict D):
    r"""
    Return a dictionary representing the vector `-X`.

    INPUT:

    - ``X`` -- a dictionary representing a vector `X`

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: D1 = {0: 1, 1: 1}
        sage: blas.negate(D1)
        {0: -1, 1: -1}
    """
    return { key: -value for key, value in D.iteritems() }

cpdef dict scal(a, dict D, bint factor_on_left=True):
    r"""
    Return a dictionary representing the vector `a*X`.

    INPUT:

    - ``a`` -- an element of the base ring `K`
    - ``X`` -- a dictionary representing a vector `X`

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: R = IntegerModRing(12)              # a ring with zero divisors
        sage: D = {0: R(1), 1: R(2), 2: R(4)}
        sage: blas.scal(R(3), D)
        {0: 3, 1: 6}
    """
    # We could use a dict comprehension like for negate, but care
    # needs to be taken to remove products that cancel out.
    # So for now we just delegate to axpy.
    return axpy(a, D, {}, factor_on_left=factor_on_left)

cpdef dict add(dict D, dict D2):
    r"""
    Return the pointwise addition of dictionaries ``D`` and ``D2``.

    INPUT:

    - ``D``, ``D2`` -- dictionaries whose values are in a common ring
      and all values are non-zero

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: D = {0: 1, 1: 1}
        sage: D2 = {0: -1, 1: 1}
        sage: blas.add(D, D2)
        {1: 2}

    `D` and `D2` are left unchanged::

        sage: D, D2
        ({0: 1, 1: 1}, {0: -1, 1: 1})
    """
    # Optimization: swap the two dicts to ensure D is the largest
    if len(D) < len(D2):
        D, D2 = D2, D
    return axpy(1, D2, D)

cpdef dict sum(dict_iter):
    r"""
    Return the pointwise addition of dictionaries with coefficients.

    INPUT:

    - ``dict_iter`` -- iterator of dictionaries whose values are in
      a common ring and all values are non-zero

    OUTPUT:

    - a dictionary containing all keys of the dictionaries in ``dict_list``
      with values being the sum of the values in the different dictionaries
      (keys with zero value are omitted)

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: D = {0: 1, 1: 1}; D
        {0: 1, 1: 1}
        sage: blas.sum(D for x in range(5))
        {0: 5, 1: 5}

        sage: D1 = {0: 1, 1: 1}; D2 = {0: -1, 1: 1}
        sage: blas.sum([D1, D2])
        {1: 2}

        sage: blas.sum([{}, {}, D1, D2, {}])
        {1: 2}
    """
    cdef dict result = {}
    cdef D
    cdef list for_removal

    for D in dict_iter:
        if result:
            iaxpy(1, D, result, remove_zeros=False)
        elif D:
            result = D.copy()

    for_removal = [key for key in result if not result[key]]
    for key in for_removal:
        del result[key]
    return result

cpdef dict linear_combination(dict_factor_iter, bint factor_on_left=True):
    r"""
    Return the pointwise addition of dictionaries with coefficients.

    INPUT:

    - ``dict_factor_iter`` -- iterator of pairs ``D``, ``coeff``, where

      * the ``D``'s are dictionaries with values in a common ring
      * the ``coeff``'s are coefficients in this ring

    - ``factor_on_left`` -- boolean (default: ``True``); if ``True``,
      the coefficients are multiplied on the left, otherwise they are
      multiplied on the right

    OUTPUT:

    - a dictionary containing all keys of dictionaries in ``dict_list``
      with values being the sum of the values in the different
      dictionaries and each one first multiplied by the given factor
      (keys with zero value are omitted)

    EXAMPLES::

        sage: import sage.data_structures.blas_dict as blas
        sage: D = { 0:1, 1:1 }; D
        {0: 1, 1: 1}
        sage: blas.linear_combination( (D,i) for i in range(5) )
        {0: 10, 1: 10}
        sage: blas.linear_combination( [(D,1), (D,-1)] )
        {}
    """
    cdef dict result = {}
    cdef dict D
    cdef list for_removal

    for D, a in dict_factor_iter:
        if not a: # We multiply by 0, so nothing to do
            continue
        if not result and a == 1:
            result = D.copy()
        else:
            iaxpy(a, D, result, remove_zeros=False)

    for_removal = [key for key in result if not result[key]]
    for key in for_removal:
        del result[key]

    return result

