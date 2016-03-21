r"""
Using restriction to scalar indices, we compute Jacobi forms of arbitrary index.

AUTHOR:

- Martin Raum

REFERENCE:

.. [Ra] Martin Raum, Computation of Jacobi forms degree 1 and higher rank index.

EXAMPLES:

To compute a basis of weight `k` and index `m` Jacobi forms, we call
``higherrank_jacobi_forms``.  The index `m` is a quadratic form over
`\ZZ`.  This compute weight `10` Jacobi forms up to precision `20`.

::

    sage: from sage.modular.jacobi.all import *
    sage: m = QuadraticForm(ZZ, 2, [1,1,1])
    sage: jforms = higherrank_jacobi_forms(10, m, 5)
    sage: jforms
    [{(0, (0, 0)): 1,
      (1, (0, 0)): 0,
      (1, (1, 1)): -45,
      (2, (0, 0)): -59130,
      (2, (1, 1)): -12672,
      (3, (0, 0)): -1416960,
      (3, (1, 1)): -558045,
      (4, (0, 0)): -14346720,
      (4, (1, 1)): -7165440},
     {(0, (0, 0)): 0,
      (1, (0, 0)): 1,
      (1, (1, 1)): -1/6,
      (2, (0, 0)): -15,
      (2, (1, 1)): 5/3,
      (3, (0, 0)): 90,
      (3, (1, 1)): -4/3,
      (4, (0, 0)): -248,
      (4, (1, 1)): -155/3}]


Fourier expansions are represented as dictionaries.  They are similar to the ones used for classical Jacobi forms.  Fourier
coefficients of Jacobi forms are index by pairs `(n,r)`, where `n` is an integer and `r` is a tuple.  The notion of reduced `r` depends on a choice and is dictated by the following output::

    sage: (r_classes, _) = higherrank_jacobi_r_classes(m)
    sage: r_classes
    [[(0, 0)],
     [(1, 1), (-1, -1), (0, 1), (0, -1), (1, 0), (-1, 0)]]

This is a list of lists, the first elements of which are, by convention, reduced.  They are guarenteed to have minimal length in their class in `\ZZ^l / m \ZZ^l`, where `l` is the rank of `m` and `m` is identified with its Gram matrix::

    sage: ZZ^m.dim() / m.matrix().row_module()
    Finitely generated module V/W over Integer Ring with invariants (3)

Fourier expansions are given in terms of Fourier coefficients of
reduced index.  They can be accessed directly; If a reduced pair does
not occur as a key and `n` does not exceed the prescribed precision,
then the corresponding Fourier coefficient vanishes.

To access Fourier coefficients of non reduced index, we compute the attached reduction by ``classical_jacobi_reduce_fe_index``::

    sage: m_adj = QuadraticForm(2 * m.matrix().adjoint())
    sage: m_span = m.matrix().row_module()
    sage: higherrank_jacobi_reduce_fe_index((4, (2,1)), m, r_classes, m_adj, m_span)
    ((3, (0, 0)), 1)

As a result, we obtain `(n', r')` and a sign `s`.  The Fourier
coefficient at `(n,r)` equals `s^k` times the coefficient at `(n',
r')`.  In the specific examples (`k` is odd), the `(4, (2,1))`-th Fourier
coefficient of the first basis element is equal to::

   sage: jforms[0][(3,(0,0))]
   -1416960
"""

#===============================================================================
#
# Copyright (C) 2012-2014 Martin Raum
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

from sage.matrix.all import matrix, zero_matrix, identity_matrix
from sage.misc.cython import cython_lambda
from sage.misc.flatten import flatten
from sage.modular.jacobi.classical import (classical_jacobi_forms,
                    classical_jacobi_fe_indices, classical_jacobi_reduce_fe_index)
from sage.modular.jacobi.higherrank_dimension import jacobi_dimension
from sage.modules.all import FreeModule, vector, zero_vector
from sage.rings.all import ZZ, QQ
from sage.quadratic_forms.all import QuadraticForm
from sage.sets.all import Set

import operator
from random import Random


def higherrank_jacobi_reduce_fe_index((n, r), m, r_classes, m_adj, m_span):
    r"""
    Reduce a Fourier index `(n, r)`.

    INPUT:

    - `(n, r)` -- A pair of an integer and a tuple of integers.

    - `m` -- A quadratic form over `\ZZ`.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_adj`` -- A quadratic form over `\ZZ`.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A pair `((n', r'), s)` where `(n', r')` is the reduced index and
      `s = \pm 1` tells whether r or -r is equivalent to r modulo `m
      \ZZ^l`.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_reduce_fe_index
        sage: m = QuadraticForm(matrix([[2]]))
        sage: m_adj = QuadraticForm(2 * m.matrix().adjoint())
        sage: m_span = m.matrix().row_module()
        sage: r_classes = higherrank_jacobi_r_classes(m)[0]
        sage: higherrank_jacobi_reduce_fe_index((1,(2,)), m, r_classes, m_adj, m_span)
        ((0, (0,)), 1)

    TESTS:

    See ``test_higherrank.py``.
    """
    try:
        (rred, sgn) = _higherrank_jacobi_reduce_fe_index__r__cache[(m,r)]
    except KeyError:
        (rred, sgn) = _higherrank_jacobi_reduce_fe_index__r(r, r_classes, m_span)
        _higherrank_jacobi_reduce_fe_index__r__cache[(m,r)] = (rred, sgn)

    nred = n - (m_adj(r) - m_adj(rred)) // (2*m.det())

    return ((nred, rred), sgn)

_higherrank_jacobi_reduce_fe_index__r__cache = {}
def _higherrank_jacobi_reduce_fe_index__r(r, r_classes, m_span):
    r"""
    Find a representative in `r_classes` that is equivalent modulo `m
    \ZZ^l` and `\pm` to `r`.

    INPUT:

    - `r` -- A tuple of integers.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A pair `(r', s)`, where `r'` is reduced and `s = \pm 1` tells
      whether r or -r is equivalent to r modulo `m \ZZ^l`.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: from sage.modular.jacobi.higherrank import _higherrank_jacobi_reduce_fe_index__r
        sage: m = QuadraticForm(matrix([[2]]))
        sage: m_span = m.matrix().row_module()
        sage: r_classes = higherrank_jacobi_r_classes(m)[0]
        sage: _higherrank_jacobi_reduce_fe_index__r((2,), r_classes, m_span)
        ((0,), 1)

    TESTS:

    See ``test_higherrank.py``.
    """
    for r_class in r_classes:
        rred = r_class[0]

        if vector(r) - vector(rred) in m_span:
            return (rred, 1)
        if vector(r) + vector(rred) in m_span:
            return (rred, -1)
    else:
        raise RuntimeError( "Could not find reduced r" )

def higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=False):
    r"""
    Indices `(n, r)` of Fourier expansions of Jacobi forms of index
    `m`, where `n` is an integer and `r` is a tuple.

    INPUT:

    - `m` -- A quadratic form over `\ZZ`.

    - ``prec`` -- A nonnegative integer.

    - ``r_classes`` -- A list of list of integers.

    - ``reduced`` -- Boolean.  Default: ``False``.

    OUTPUT:

    - A generator of pairs `(n, r)`, where `n` is an integer and `r` is a tuple.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_fe_indices
        sage: m = QuadraticForm(matrix([[2]]))
        sage: r_classes = higherrank_jacobi_r_classes(m)[0]
        sage: list(higherrank_jacobi_fe_indices(m, 2, r_classes, reduced=True))
        [(0, (0,)), (1, (0,)), (1, (1,))]
        sage: list(higherrank_jacobi_fe_indices(m, 2, r_classes, reduced=False))
        [(0, (0,)), (1, (0,)), (1, (1,)), (1, (-1,)), (1, (2,)), (1, (-2,))]

    TESTS::

        sage: higherrank_jacobi_fe_indices(m, 2, r_classes, reduced=True)
        <generator object ...

    See also ``test_higherrank.py``.
    """
    m_adj = QuadraticForm(2 * m.matrix().adjoint())

    if reduced:
        for n in range(0,prec):
            for r_class in r_classes:
                r = r_class[0]
                if m_adj(r) <= 2*m.det()*n:
                    yield (n, r)
    else:
        short_vectors = m_adj.short_vector_list_up_to_length(2*m.det()*(prec - 1) + 1)
        for n in range(0, prec):
            for length in range(2 * m.det() * n + 1):
                for r in short_vectors[length]:
                    yield (n, tuple(r))

    raise StopIteration

def higherrank_jacobi_r_classes(m):
    r"""
    Let `l` be the dimension of `m`.  For each element of `(\ZZ^l / m
    \ZZ^l) \pm m` of minimal norm, we compute all representatives that
    minimize the norm with respect to the adjoint of `m`.

    INPUT:

    - `m` -- A quadratic form over `\ZZ`.

    OUTPUT:

    - A list of lists of vectors.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: m = QuadraticForm(matrix([[2]]))
        sage: higherrank_jacobi_r_classes(m)
        ([[(0,)], [(1,), (-1,)]], [[1], [1, 1]])

    TESTS:

    See ``test_higherrank.py``.
    """
    m_mat = m.matrix()
    m_span = m_mat.row_module()
    m_adj = QuadraticForm(2 * m_mat.adjoint())

    canonical_reps = [r.lift() for r in m_span.ambient_module() / m_span]
    max_norm = max(m_adj(r) for r in canonical_reps)

    current_max_length = 5
    short_vectors = m_adj.short_vector_list_up_to_length(current_max_length)

    r_classes = []
    r_classes_reduction_signs = []
    for r_can in canonical_reps:
        r_class_found = False
        for r_class in r_classes:
            if (vector(r_class[0]) - vector(r_can) in m_span
                or vector(r_class[0]) + vector(r_can) in m_span):
                r_class_found = True
                break
        if r_class_found:
            continue

        r_classes.append([])
        r_class = r_classes[-1]
        r_classes_reduction_signs.append([])
        r_class_reduction_signs = r_classes_reduction_signs[-1]

        # travers short vectors with respect to m_adj until we find a
        # vector that is equivalent to \pm r_can.modulo m \Z^l
        for length in range(0, max_norm + 1):
            if length >= current_max_length:
                current_max_length += 5
                short_vectors = m_adj.short_vector_list_up_to_length(current_max_length)

            for r in short_vectors[length]:
                if r_can - r in m_span:
                    r_class.append(tuple(r))
                    r_class_reduction_signs.append(1)
                elif r_can + r in m_span:
                    r_class.append(tuple(r))
                    r_class_reduction_signs.append(-1)

            if len(r_class) != 0:
                break

    for cl_ix in range(len(r_classes_reduction_signs)):
        if r_classes_reduction_signs[cl_ix][0] == -1:
            r_classes_reduction_signs[cl_ix] = map(operator.neg,
                                                   r_classes_reduction_signs[cl_ix])


    return (r_classes, r_classes_reduction_signs)

def higherrank_jacobi_forms(k, m, prec, algorithm="restriction"):
    r"""
    Compute the Fourier expansions of a basis of Jacobi forms (over
    `\QQ`) of weight `k` and index `m` (an quadratic form over `\ZZ`)
    up to given precision.

    ALGORITHM:

    See [Ra]_. The algorithm in [Ra]_ is applied for precision
    ``relation_prec``.  After this, the remaining Fourier coefficients
    are determined using as few restrictions as possible.

    INPUT:

    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``algorithm`` -- Only "restriction" is implemented.

    OUTPUT:

    A list of dictionaries, which describes the Fourier expansion of
    Jacobi forms.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_forms
        sage: k = 10
        sage: m = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: jforms = higherrank_jacobi_forms(k, m, 3)
        sage: Sequence(jforms, cr=True)
        [
        {(1, (1, 1)): -45, (2, (1, 1)): -12672, (1, (0, 0)): 0, (2, (0, 0)): -59130, (0, (0, 0)): 1},
        {(1, (1, 1)): -1/6, (2, (1, 1)): 5/3, (1, (0, 0)): 1, (2, (0, 0)): -15, (0, (0, 0)): 0}
        ]

    We access these the Fourier coefficients by means of the indices
    `n` and `r`, typical for Jacaobi forms.

    ::

        sage: n = 2; r = (1, 1)
        sage: jforms[0][(n,r)]
        -12672

    This works, since `r = (1,1)` is a reduced vector.  For general
    `r`, we have to invoke index reduction, to find a corresponding
    index of the dictionary.

    ::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_reduce_fe_index
        sage: n = 2; r = (2, 2)
        sage: m_adj = QuadraticForm(2 * m.matrix().adjoint())
        sage: r_classes = higherrank_jacobi_r_classes(m)[0]
        sage: m_span = m.matrix().row_module()
        sage: ((nred, rred), _) = higherrank_jacobi_reduce_fe_index((n,r), m, r_classes, m_adj, m_span)
        sage: jforms[0][(nred, rred)]
        -45

    TESTS:

    See ``test_higherrank.py``.
    """
    if algorithm != "restriction":
        raise NotImplementedError("Algorithm {} is not implemented.".format(algorithm))
    rand = Random()

    dim = jacobi_dimension(k, m)
    if dim == 0:
        return []

    (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
    m_span = m.matrix().row_module()

    rst_vectors_with_image = _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span)
    rst_vectors = [vector(s)
                   for s in Set(tuple(s) for (s, _) in rst_vectors_with_image)]

    max_rst_index = max([m(s) for s in rst_vectors])
    minimal_prec = 2 + (k + max_rst_index) // 12 + 5
    relation_prec = minimal_prec

    max_relation_rst_index = max_rst_index
    relation_rst_vectors = []


    while True:
        try:
            prec__max = max(prec, relation_prec)
            jforms = _higherrank_jacobi_forms__restriction(
                k, prec__max, relation_prec, dim,
                *_restriction_relation_matrices(k, m,
                                                prec__max, relation_prec,
                                                rst_vectors, relation_rst_vectors,
                                                r_classes, m_span)
            )

            if prec__max < prec:
                return jforms
            else:
                return [dict(((n,r), c) for ((n,r), c) in phi.items()
                             if n < prec)
                        for phi in jforms]

        except ValueError, err:
            if len(err.args) < 2 or err.args[1] != "INSUFFICIENT RELATIONS":
                raise

            relation_prec += minimal_prec
            max_relation_rst_index += 2

            if len(relation_rst_vectors):
                relation_rst_vectors = rst_vectors

            short_vectors = (
                flatten( m.short_vector_list_up_to_length(max_relation_rst_index+1, True)[max_relation_rst_index-2:],
                         max_level=1)
            )
            try:
                relation_rst_vectors += rand.sample(short_vectors, m.det())
            except ValueError:
                relation_rst_vectors += short_vectors


def _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span):
    r"""
    Given classes ``r_classes`` of elements in `\ZZ^l` find a complete
    set of restriction vectors.

    INPUT:

    - `m` -- A quadratic form.

    - ``r_classes`` -- A list of lists of tuples of integers.

    - ``r_classes_reduction_signs`` -- A list of lists of `\pm 1`.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A set of pairs, the first of which is a vector corresponding to
      an element in `\ZZ^l`, and the second of which is an integer.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import _complete_set_of_restriction_vectors
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: m = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: m_span = m.matrix().row_module()
        sage: (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
        sage: _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span)
        [((1, 0), 0), ((-2, 1), 1)]

    TESTS:

    See ``test_higherrank.py``
    """
    r_classes = [map(vector, r_class) for r_class in r_classes]

    length_inc = 5
    max_length = 5
    cur_length = 1
    short_vectors = m.short_vector_list_up_to_length(max_length+1, True)

    pm_fixed_point_indices = [ix for (ix, r_class) in enumerate(r_classes)
                               if 2*vector(r_class[0]) in m_span]

    rst_vectors = []
    nmb_rst_vectors_even = 0
    nmb_rst_vectors_odd = 0
    rst_kernel_even = FreeModule(QQ, len(r_classes))
    rst_kernel_odd_matrix = identity_matrix(QQ, len(r_classes))
    for ix in pm_fixed_point_indices:
        rst_kernel_odd_matrix[ix,ix] = 0
    rst_kernel_odd = rst_kernel_odd_matrix.row_module()


    while (nmb_rst_vectors_even < len(r_classes)):
        while len(short_vectors[cur_length]) == 0:
            cur_length += 1
            if max_length < cur_length:
                max_length += length_inc
                short_vectors = m.short_vector_list_up_to_length(max_length+1, True)

        s = vector( short_vectors[cur_length].pop() )


        rst_imgs_even = {}
        rst_imgs_odd = {}
        for (cl_ix, (r_class, r_signs)) \
            in enumerate(zip(r_classes, r_classes_reduction_signs)):

            for (r, r_sign) in zip(r_class, r_signs):
                r_rst = s.dot_product(r)

                if r_rst not in rst_imgs_even:
                    rst_imgs_even[r_rst] = zero_vector(len(r_classes))
                rst_imgs_even[r_rst][cl_ix] += 1

                if cl_ix not in pm_fixed_point_indices:
                    if r_rst not in rst_imgs_odd:
                        rst_imgs_odd[r_rst] = zero_vector(len(r_classes))
                    rst_imgs_odd[r_rst][cl_ix] += r_sign

        for r_rst in rst_imgs_even.keys():
            if (rst_kernel_even.basis_matrix() * rst_imgs_even[r_rst]).is_zero():
                continue
            if r_rst in rst_imgs_odd:
                contributes_to_odd = not (rst_kernel_odd.basis_matrix()
                                          * rst_imgs_odd[r_rst]).is_zero()

            if (nmb_rst_vectors_odd + len(pm_fixed_point_indices)
                <= nmb_rst_vectors_even
                and
                (r_rst not in rst_imgs_odd or not contributes_to_odd)):

                continue

            rst_vectors.append((s, r_rst))

            rst_kernel_even = rst_kernel_even.intersection(
                matrix(rst_imgs_even[r_rst]).right_kernel())
            nmb_rst_vectors_even += 1
            if r_rst in rst_imgs_odd and contributes_to_odd:
                rst_kernel_odd = rst_kernel_odd.intersection(
                    matrix(rst_imgs_odd[r_rst]).right_kernel())
                nmb_rst_vectors_odd += 1


    return rst_vectors

def _restriction_relation_matrices(k, m, prec, relation_prec,
                                   rst_vectors, relation_rst_vectors,
                                   r_classes, m_span):
    r"""
    INPUT:

    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``relation_prec`` -- A nonnegative integer.  Precision less than
                           or equal to ``prec`` up to which relations
                           of coefficients are computed.

    - ``rst_vectors`` -- A list of vectors.

    - ``relation_rst_vectors`` -- Compute relations for a give set of
                                  restriciton vectors.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A quintuple.  See `meth:_restriction_matrix` and
      `meth:_relation_matrix` for a more detailed description.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import _restriction_relation_matrices
        sage: from sage.modular.jacobi.higherrank import _complete_set_of_restriction_vectors
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: m = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: m_span = m.matrix().row_module()
        sage: prec = 2; relation_prec = 1
        sage: (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
        sage: rst_vectors_with_image = _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span)
        sage: rst_vectors = [vector(s) for s in Set(tuple(s) for (s, _) in rst_vectors_with_image)]
        sage: relation_rst_vectors = rst_vectors + [vector(ZZ, [2,0])]
        sage: Sequence(_restriction_relation_matrices(1, m, prec, relation_prec, rst_vectors, relation_rst_vectors, r_classes, m_span), cr=True)
        [
        [ 0  1  0]
        [ 2  0  0]
        [ 1  0  0]
        [ 2  1  0]
        [ 0  0 -2]
        [ 0  0  1]
        [ 2  0  0]
        [ 1  0  0], [((1, 0), 1, 0, 3), ((-2, 1), 3, 3, 5)], {1: {(1, 0): 0, (0, 0): 2, (1, 1): 1}, 3: {(1, 2): 2, (1, 0): 0, (1, 3): 3, (1, 1): 1, (0, 0): 4}},
        <BLANKLINE>
        [(0, (0, 0)), (1, (0, 0)), (1, (1, 1))], [1], [(0, (0, 0))]
        ]

    TESTS:

    Tested implicitely by ``meth:higherrank_jacobi_forms``.  See also ``test_higherrank.py``.
    """
    (restriction_matrix__big, row_groups, row_labels, column_labels) = \
        _restriction_matrix(k, m, prec, rst_vectors, False, r_classes, m_span)
    (relation_matrix, column_labels_relations) = \
        _relation_matrix(k, m, relation_prec, relation_rst_vectors, r_classes, m_span)

    restriction_matrix__big = restriction_matrix__big.change_ring(QQ)
    relation_matrix = relation_matrix.change_ring(QQ)

    return ( restriction_matrix__big, row_groups, row_labels, column_labels,
             relation_matrix, column_labels_relations )


def _restriction_matrix(k, m, prec, rst_vectors, find_relations, r_classes, m_span):
    r"""
    A matrix that maps the Fourier expansion of a Jacobi form of given precision
    to their restrictions with respect to the elements of S.

    INPUT:

    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``rst_vectors`` -- A list of vectors.

    - ``find_relation`` -- A boolean. If ``True``, then the restrictions to
                           nonreduced indices will also be computed.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    A quadruple ``(restriction_matrix, row_groups, row_labels, column_labels)``.

    - ``restriction_matrix`` -- A matrix that describes the
      restriction of Fourier expansion of Jacobi forms of index `m` so
      classical Jacobi forms.  It acts on column vectors.
      ``row_groups``, ``row_labels``, ``column_labels`` describe which
      entry corresponds to which Fourier index.

    - ``row_groups`` -- A list of of quadruples ``(s, m_rst, start,
      length)``.  This means that rows ``start`` to ``start + length``
      contain the image of Fourier expansions to `s z`.  The index of
      this image corresponds to classical Jacobi forms of index `m_rst`.

    - ``row_labels`` -- A dictionary, which assignes to each index `m_rst`
      of classical Jacobi forms that occur in ``row_groups`` a
      labelling.  The values of this dictionary are dictionaries
      themselves, which map pairs `(n,r)` of integers to integer
      indices ``ix``.  This means that for every restriction whose
      image has index `m` the row ``start + ix`` (``start`` was given
      above) corresponds to the Fourier coefficients of index `(n,r)`.

    - ``column_labels`` -- A dictionary that maps Fourier indices of
      Jacobi forms of index `m` to integer indices `ix`.  This means
      that that the `ix`-th column of the restrictio matrix
      corresponds to Fourier coefficients of index `(n,r)`.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import _restriction_matrix
        sage: from sage.modular.jacobi.higherrank import _complete_set_of_restriction_vectors
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: m = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: m_span = m.matrix().row_module()
        sage: prec = 2
        sage: (r_classes, r_classes_reduction_signs) = higherrank_jacobi_r_classes(m)
        sage: rst_vectors_with_image = _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs, m_span)
        sage: rst_vectors = [vector(s) for s in Set(tuple(s) for (s, _) in rst_vectors_with_image)]
        sage: Sequence(_restriction_matrix(0, m, prec, rst_vectors, False, r_classes, m_span), cr=True)
        [
        [0 1 2]
        [2 0 2]
        [1 0 0]
        [2 1 0]
        [0 0 2]
        [0 0 1]
        [2 0 0]
        [1 0 0],
        [((1, 0), 1, 0, 3), ((-2, 1), 3, 3, 5)], {1: {(1, 0): 0, (0, 0): 2, (1, 1): 1}, 3: {(1, 2): 2, (1, 0): 0, (1, 3): 3, (1, 1): 1, (0, 0): 4}},
        [(0, (0, 0)), (1, (0, 0)), (1, (1, 1))]
        ]

    TESTS:

    Tested implicitely by ``meth:higherrank_jacobi_forms``.  See also ``test_higherrank.py``.
    """
    k = k % 2
    m_adj = QuadraticForm(2 * m.matrix().adjoint())

    column_labels = list(higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=True))

    if len(rst_vectors) == 0:
        return (zero_matrix(ZZ, 0, len(column_labels)), [], {}, column_labels)

    rst_jacobi_indices = [m(s) for s in rst_vectors]
    rst_indices = dict((m_rst,
                        list(classical_jacobi_fe_indices(
                            m_rst, prec, reduced=not find_relations)))
                       for m_rst in Set(rst_jacobi_indices))

    row_groups = [len(rst_indices[mm_rst]) for mm_rst in rst_jacobi_indices]
    row_groups = [(s, mm_rst, sum(row_groups[:i]), row_groups[i])
                   for ((i, s), mm_rst) in zip(enumerate(rst_vectors), rst_jacobi_indices) ]
    row_labels = dict((m_rst, dict( (nr, i) for (i, nr) in enumerate(rst_indices[m_rst])))
                       for m_rst in Set(rst_jacobi_indices))

    reductions = dict((nr, []) for nr in column_labels)
    for nr in higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=False):
        (nrred, sgn) = higherrank_jacobi_reduce_fe_index(nr, m, r_classes, m_adj, m_span)
        reductions[nrred].append((nr, sgn))

    if sum(map(len, reductions.items())) > 10000:
        cython_dot_products = True
        dot_products = [ cython_lambda(
            ' , '.join([ 'int x{}'.format(ix) for ix in range(len(s)) ]),
            ' + '.join([ '{1}*x{0}'.format(*six) for six in enumerate(s) ]) )
                         for (s, _, _, _) in row_groups ]
    else:
        cython_dot_products = False

    restriction_matrix = zero_matrix(ZZ, row_groups[-1][2] + row_groups[-1][3],
                                     len(column_labels))
    for (col, nrred) in enumerate(column_labels):
        for ((n, r), sgn) in reductions[nrred]:
            for (ix, (s, m_rst, start, length)) in enumerate(row_groups):
                row_labels_dict = row_labels[m_rst]

                if cython_dot_products:
                    rst_r = dot_products[ix](*r)
                else:
                    rst_r = s.dot_product(vector(r))

                try:
                    restriction_matrix[start + row_labels_dict[(n, rst_r)], col] \
                        += 1 if k == 0 else sgn
                except KeyError:
                    assert not find_relations

    return (restriction_matrix, row_groups, row_labels, column_labels)


def _relation_matrix(k, m, prec, rst_vectors, r_classes, m_span):
    r"""
    Deduce relations of the coefficients of a Jacobi form from their
    specialization to Jacobi form of scalar index.

    INPUT:

    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``rst_vectors`` -- A list of vectors.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    A pair ``(relation_matrix, column_labels)``.

    - ``relation_matrix`` -- A matrix whose row space corresponds to
      relations of Fourier coefficients that are deduced from
      restrictions.  That is, linear combinations of Fourier
      coefficients described as the rows of this matrix must vanish.

    - ``column_labels`` -- See `meth:_restriction_matrix`.

    EXAMPLES::

        sage: from sage.modular.jacobi.higherrank import _relation_matrix
        sage: from sage.modular.jacobi.higherrank import _complete_set_of_restriction_vectors
        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_r_classes
        sage: m = QuadraticForm(matrix(2, [2,0,0,4]))
        sage: m_span = m.matrix().row_module()
        sage: prec = 2
        sage: (r_classes, _) = higherrank_jacobi_r_classes(m)
        sage: Sequence(_relation_matrix(1, m, prec, [], r_classes, m_span), cr=True)
        [
        [1 0 0 0 0 0 0]
        [0 1 0 0 0 0 0]
        [0 0 0 1 0 0 0]
        [0 0 0 0 1 0 0]
        [0 0 0 0 0 0 1],
        [(0, (0, 0)), (1, (0, 0)), (1, (0, 1)), (1, (0, 2)), (1, (1, 0)), (1, (1, 1)), (1, (1, 2))]
        ]

    TESTS:

    Tested implicitely by ``meth:higherrank_jacobi_forms``.  See also ``test_higherrank.py``.
    """
    k = k % 2
    # m_adj = QuadraticForm(2 * m.matrix().adjoint())

    ## relations computed from restrictions
    (mat, row_groups, row_labels, column_labels) = \
        _restriction_matrix(k, m, prec, rst_vectors, True, r_classes, m_span)

    relations = []
    for (s, m_rst, start, length) in row_groups:
        row_labels_dict = row_labels[m_rst]
        for (nr, ix) in row_labels_dict.items():
            (nrred, sgn) = classical_jacobi_reduce_fe_index(nr, m_rst)
            if nrred == nr:
                continue

            rel = (mat.row(start + row_labels_dict[nrred])
                   - (1 if k == 0 else sgn) * mat.row(start + ix))
            if rel != 0:
                relations.append(rel)

    ## Forces zeros in case of odd k
    if k == 1:
        for (ix, (n, r)) in enumerate(column_labels):
            if 2 * vector(r) in m_span:
                rel = zero_vector(len(column_labels))
                rel[ix] = 1
                relations.append(rel)

    return (matrix(len(relations), len(column_labels), relations), column_labels)


def _higherrank_jacobi_forms__restriction(
        k, prec, relation_prec, dim,
        restriction_matrix__big, row_groups, row_labels, column_labels,
        relation_matrix, column_labels_relations):
    r"""
    Compute the Fourier expansions of Jacobi forms (over `\Q`) of weight `k` and
    index `m` up to given precision.

    INPUT:

    - `k` -- An integer.  Only `k` modulo `2` is used.

    - ``prec`` -- A nonnegative integer.

    - ``relation_prec`` -- A nonnegative integer.

    - ``dim`` -- The dimension of the space of Jacobi forms.

    - ``restriction_matrix__big`` -- See output of
                                     `meth:_restriction_relation_matrices`.

    - ``row_groups`` -- See output of `meth:_restriction_relation_matrices`.

    - ``row_labels`` -- See output of `meth:_restriction_relation_matrices`.

    - ``column_labels`` -- See output of `meth:_restriction_relation_matrices`.

    - ``relation_matrix`` -- See output of `meth:_restriction_relation_matrices`.

    - ``column_labels_relations`` -- See output of
                                     `meth:_restriction_relation_matrices`.

    OUTPUT:

    A list of dictionaries.

    NOTE:

    A minimal example of how to call this function can be found in
    `meth:higherrank_jacobi_forms`.

    TESTS::

        sage: from sage.modular.jacobi.higherrank import higherrank_jacobi_forms
        sage: k = 8
        sage: m = QuadraticForm(matrix(2, [2,1,1,2]))
        sage: higherrank_jacobi_forms(k, m, 1)
        [{(0, (0, 0)): 1}]

    See also ``test_higherrank.py``.
    """
    assert relation_prec <= prec

    ## Construct a restriction matrix for relation_prec
    row_groups__small = [ len(filter( lambda (n,r): n < relation_prec, row_labels[m_rst].keys()))
                          for (_, m_rst, _, _) in row_groups ]
    row_groups__small = [ (sum(row_groups__small[:i]), row_groups__small[i])
                          for i in range(len(row_groups)) ]

    row_labels__small = dict()
    for (m_rst, row_labels_dict) in row_labels.items():
        row_labels_dict__small = {}

        label_nmb = 0
        for (nr, i) in row_labels_dict.items():
            if nr[0] < relation_prec:
                row_labels_dict__small[nr] = (label_nmb, i)
                label_nmb += 1

        row_labels__small[m_rst] = row_labels_dict__small

    row_indices__small = list()
    for ((s, m_rst, start, _), (start_small, length_small)) in zip(row_groups, row_groups__small):
        row_labels_dict = row_labels__small[m_rst]
        row_indices__sub = length_small * [None]
        for (_, (i, i_pre)) in row_labels_dict.items():
            row_indices__sub[i] = start + i_pre
        row_indices__small += row_indices__sub

    restriction_matrix = restriction_matrix__big \
        .matrix_from_rows_and_columns(
            row_indices__small,
            [column_labels.index(nr) for nr in column_labels_relations] )


    rst_jacobi_indices = [m_rst for (_, m_rst, _, _) in row_groups]
    # rst_indices = dict( (m_rst, list(classical_jacobi_fe_indices(m_rst, prec, reduced=True)))
    #                    for m_rst in Set(rst_jacobi_indices) )
    rst_jacobi_forms = dict( (m_rst, classical_jacobi_forms(k, m_rst, prec))
                             for m_rst in Set(rst_jacobi_indices) )

    rst_jacobi_vectors = []
    nmb_rst_coords = row_groups[-1][2] + row_groups[-1][3]
    for (s, m_rst, start, length) in row_groups:
        row_labels_dict = row_labels[m_rst]
        for phi in rst_jacobi_forms[m_rst]:
            v = vector(ZZ, len(row_labels_dict))
            for (nr, i) in row_labels_dict.items():
                (nrred, sgn) = classical_jacobi_reduce_fe_index(nr, m_rst)
                v[i] = ((1 if k % 2 == 0 else sgn) * phi[nr]) if nr in phi else 0

            rst_jacobi_vectors.append(vector(
                start * [0] + v.list() + (nmb_rst_coords - start - length) * [0] ))

    rst_jacobi_matrix__big = \
        matrix(len(rst_jacobi_vectors), nmb_rst_coords, rst_jacobi_vectors).transpose()
    rst_jacobi_matrix = rst_jacobi_matrix__big.matrix_from_rows(row_indices__small)
    rst_jacobi_space = rst_jacobi_matrix.column_module()

    rst_expansions_matrix = restriction_matrix.column_module() \
                                .intersection(rst_jacobi_space) \
                                .basis_matrix().transpose()
    rst_preimage = restriction_matrix.solve_right(rst_expansions_matrix).column_space() \
                   + restriction_matrix.right_kernel().change_ring(QQ)
    jacobi_expansions_space = rst_preimage.intersection(relation_matrix.right_kernel())


    if jacobi_expansions_space.dimension() < dim:
        raise RuntimeError( "There is a bug in the implementation of the restriction method. Dimensions: {} < {}!".format(jacobi_expansions_space.dimension(), dim) )
    if jacobi_expansions_space.dimension() > dim:
        raise ValueError( "Could not construct enough restrictions to determine Fourier expansions uniquely", "INSUFFICIENT RELATIONS" )


    ## reconstruct the whole Fourier expansion from partial ones
    restriction_coordinates = rst_jacobi_matrix.solve_right(
        restriction_matrix * jacobi_expansions_space.basis_matrix().transpose() )
    jacobi_expansions__big = restriction_matrix__big.solve_right( rst_jacobi_matrix__big * restriction_coordinates ).transpose()

    return [dict((nr, c) for (nr, c) in zip(column_labels, phi))
            for phi in jacobi_expansions__big.rows()]
