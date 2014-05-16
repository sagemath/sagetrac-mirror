r"""
Using restriction to scalar indices, we compute Jacobi forms of arbitrary index.

AUTHOR:

- Martin Raum

REFERENCE:

- [Ra] Martin Raum, Computation of Jacobi forms degree 1 and higher rank index.
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

from sage.matrix.constructor import matrix, zero_matrix
from sage.misc.cython import cython_lambda
from sage.misc.flatten import flatten
from sage.modular.jacobi.classical import (classical_jacobi_forms,
                    classical_jacobi_fe_indices, reduce_classical_jacobi_fe_index)
from sage.modular.jacobi.higherrank_dimension import dimension_jacobi
from sage.modules.all import FreeModule, vector, span
from sage.rings.all import ZZ, QQ
from sage.quadratic_forms.all import QuadraticForm
from sage.sets.all import Set

import operator
from random import Random


def reduce_higherrank_jacobi_fe_index((n, r), m, r_classes, m_span):
    r"""
    Reduce a Fourier index `(n, r)`.

    INPUT:

    - `(n, r)` -- A pair of an integer and a tuple of integers.

    - `m` -- A quadratic form over `\Z`.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A pair `((n', r'), s)` where `(n', r')` is the reduced index and
      `s = \pm 1` tells whether r or -r is equivalent to r modulo `m
      \Z^l`.
    """
    (rred, sgn) = _reduce_higherrank_jacobi_fe_index__r(r, r_classes, m_span)
    nred = n - (m(r) - m(rred)) // (2*m.det())

    return ((nred, rred), sgn)

def _reduce_higherrank_jacobi_fe_index__r(r, r_classes, m_span):
    r"""
    Find a representative in `r_classes` that is equivalent modulo `m
    \Z^l` and `\pm` to `r`.

    INPUT:

    - `r` -- A tuple of integers.

    - ``r_classes`` -- A list of lists of vectors.

    - ``m_span`` -- The row (or column) span `m`.

    OUTPUT:

    - A pair `(r', s)`, where `r'` is reduced and `s = \pm 1` tells
      whether r or -r is equivalent to r modulo `m \Z^l`.
    """
    for r_class in r_classes:
        rred = r_class[0]
        r_rred = vector(r) - vector(rred)
        if r_rred in m_sp:
            return (rred, 1)

        r_rred = vector(r) + vector(rred)
        if r_rred in m_span:
            return (rred, -1)
    else :
        raise RuntimeError( "Could not find reduced r" )

def higherrank_jacobi_fe_indices(m, prec, r_classes, reduced=False):
    r"""
    Indices `(n, r)` of Fourier expansions of Jacobi forms of index
    `m`, where `n` is an integer and `r` is a tuple.

    INPUT:

    - `m` -- A quadratic form over `\Z`.

    - ``prec`` -- A nonnegative integer.

    - ``r_classes`` -- A list of list of integers.

    - ``reduced`` -- Boolean.  Default: ``False``.

    OUTPUT:

    - A generator of pairs `(n, r)`, where `n` is an integer and `r` is a tupel.
    """
    m_adj = QuadraticForm(2 * m.matrix().adjoint())

    if reduced:
        for n in range(0,prec):
            for r_class in r_classes:
                r = r_class[0]
                if m_adj(r) <= 2*m.det()*n:
                    yield (n, r)
    else :
        short_vectors = m_adj.short_vector_list_up_to_length(2*m.det()*(prec - 1) + 1)
        for n in range(0, prec):
            for length in range(0, 2*m.det()*n + 1):
                for r in short_vectors[length] :
                    yield (n, r)

    raise StopIteration

def _higherrank_jacobi_r_classes(m):
    r"""
    Let `l` be the dimension of `m`.  For each element of `(\Z^l / m
    \Z^l) \pm m` of minimal norm, we compute all representatives that
    minimize the norm with respect to the adjoint of `m`.

    INPUT:

    - `m` -- A quadratic form over `\Z`.

    OUTPUT:

    - A list of lists of vectors.
    """
    m_mat = m.matrix()
    m_span = m_man.row_module()
    m_adj = QuadraticForm(2 * m_mat.adj())
    
    
    canonical_reps =  [r.lift() for r in m_span.ambient_module() / m_span]
    max_norm = max(m_adj(r) for r in canonical_reps)

    def recompute_short_vectors(current_max_length):
        return (current_max_length, short_vectors)
    (current_max_length, short_vectors) = recompute_short_vectors(0)

    r_classes = []
    r_classes_reduction_signs = []
    for r_can in canonical_reps:
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

            if len(r_class) != 0: break

    for cl_ix in range(len(r_classes_reduction_signs)):
        if r_classes_reduction_signs[cl_ix][0] == -1:
            r_classes_reduction_signs[cl_ix] = map(operator.neg, r_classes_reduction_signs[cl_ix])


    return (r_classes, r_classes_reduction_signs)

def higherrank_jacobi_forms(k, m, prec, algorithm="restriction"):
    r"""
    Compute the Fourier expansions of Jacobi forms (over `\QQ`) of weight `k` and 
    index `m` (an quadratic form over `\Z`) up to given precision.
    
    ALGORITHM:
    
    See [Ra]. The algorithm in [Ra] is applied for precision
    ``relation_prec``.  After this, the remaining Fourier coefficients
    are determined using as few restrictions as possible.
    
    INPUT:
    
    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``algorithm`` -- Only "restriction" is implemented.

    OUTPUT:
    
    - A list of dictionaries.
    
    TESTS::
    
        sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
        sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _coefficient_by_restriction
        sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
        sage: precision = indices.filter(20)
        sage: relation_precision = indices.filter(10)
        sage: _coefficient_by_restriction(precision, 10) == _coefficient_by_restriction(precision, 10, relation_precision) 
        True
    """
    if algorithm != "restriction":
        raise NotImplementedError("Algorithm {} is not implemented.".format(algorithm))
    rand = Random()


    dim = jacobi_dimension(k, m)
    if dim == 0: return []


    (r_classes, r_classes_reduction_signs) = _higherrank_jacobi_r_classes(k, m)
    m_span = m.matrix().row_module()


    rst_vectors_with_image = _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs)
    rst_vectors = Set(s for (s, _) in rst_vectors_with_image).list()


    max_rst_index = max([m(s) for s in rst_vectors])
    minimal_prec = 1 + (k + max_rst_index) // 12
    prec = max(minimal_prec, prec)
    relation_prec = minimal_prec


    relation_rst_index = max_rst_index - 1
    all_relation_rst_vectors = \
        flatten( m.short_vector_list_up_to_length(max_rst_index+1), max_level = 1 )
    if relation_rst_vectors is None:
        # We choose some s in order to compute relations, hoping that this
        # will be enough.  Using all just takes too long.
        relation_rst_vectors = []
        for _ in range(4*m.det()):
            s = rand.choice(all_relation_rst_vectors)
            if s not in relation_rst_vectors: relation_rst_vectors.append(s)

    while True:
        try:
            return _higherrank_jacobi_forms__restriction(
                k, prec, relation_prec, dim,
                *_restriction_relation_matrices(k, m, prec, relation_prec,
                                                rst_vectors, relation_rst_vectors,
                                                r_classes, m_span)
            )
        except ValueError, err:
            if len(err.args) == 1 and err.args[1] != "INSUFFICIENT RELATIONS":
                raise

            relation_prec += minimal_prec
            relation_rst_index += 1
            relation_rst_vectors = \
                flatten( m.short_vector_list_up_to_length(max_relation_rst_index+1), max_level = 1 )

def _complete_set_of_restriction_vectors(m, r_classes, r_classes_reduction_signs):
    r"""
    Given classes ``r_classes`` of elements in `\Z^l` find a complete
    set of restriction vectors.
    
    INPUT:
    
    - `m` -- A quadratic form.
    
    - ``r_classes`` -- A list of lists of tuples of integers.

    - ``r_classes_reduction_signs`` -- A list of lists of `\pm 1`.

    OUTPUT:
    
    - A set of pairs, the first of which is a vector corresponding to
      an element in `\ZZ^l`, and the second of which is an integer.
    
    # TESTS::
    
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _find_complete_set_of_restriction_vectors
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
    #     sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
    #     sage: _find_complete_set_of_restriction_vectors(indices.jacobi_index(), indices._r_representatives)
    #     [((-1, 0), 0), ((-1, 0), 1), ((2, -1), 1)]
    #     sage: _find_complete_set_of_restriction_vectors(indices.jacobi_index(), indices._r_representatives, reduction_function = indices.reduce_r)
    #     [((-1, 0), 0), ((-1, 0), 1), ((2, -1), 1)]
        
    # ::
     
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _local_restriction_matrix
    #     sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(4, [2,0,0,1, 0,2,0,1, 0,0,2,1, 1,1,1,2])))
    #     sage: S = _find_complete_set_of_restriction_vectors(indices.jacobi_index(), indices._r_representatives)
    #     sage: _local_restriction_matrix(indices._r_representatives, S).rank()
    #     4
    """
    r_classes = [map(vector, r_class) for r_class in r_classes]
    
    length_inc = 5
    max_length = 5
    cur_length = 1
    short_vectors = m.short_vector_list_up_to_length(max_length)
    
    rst_vectors = []
    rst_space = FreeModule(QQ, len(r_classes)).span([])
    
    while (len(rst_vectors) < len(r_classes)) :
        while len(short_vectors[2 * cur_length]) == 0:
            cur_length += 1
            if max_length < cur_length:
                max_length += length_inc
                short_vectors = m.short_vector_list_up_to_length(max_length)
        
        s = vector( short_vectors[2 * cur_length].pop() )
        
        restricted_r_candidates = Set([ s.dot_product(r)
                                        for r_class in r_classes for r in r_class ])
        
        for rst_r in restricted_r_candidates:
            v = vector([ sum(sgn for (r,sgn) in zip(r_classes, r_class_signs)
                             if s.dot_product(r) == rst_r)
                         for (r_class, r_class_signs) in zip(r_classes, r_classes_reduction_signs) ])
            if v not in rst_space :
                rst_vectors.append((s, rst_r))
                rst_space = rst_space + FreeModule(QQ, len(R)).span([v])
                
                if len(rst_vectors) == len(r_classes):
                    break
    
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

    - ``r_classes`` -- INSERT

    - ``m_span`` -- INSERT
    """
    (restriction_matrix__big, row_groups, row_labels, column_labels) = \
        _restriction_matrix(k, m, prec, rst_vectors, False, r_classes, m_span)
    (relation_matrix, column_labels_relations) = \
        _relation_matrix(k, m, relation_prec, relation_rst_vectors, r_classes, m_span)
    restriction_matrix__big.change_ring(QQ)
    relation_matrix.change_ring(QQ)

    return ( restriction_matrix__big, row_groups, row_labels, column_labels,
             relation_matrix, column_labels_relations )

def _restriction_matrix(k, m, prec, rst_vectors, find_relations, r_classes, m_span) :
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

    - ``r_classes`` -- INSERT

    - ``m_span`` -- INSERT
                           
    # TESTS::
    
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _global_restriction_matrix
    #     sage: precision = JacobiFormD1Filter(5, QuadraticForm(matrix(2, [2,1,1,2])))
    #     sage: (global_restriction_matrix, row_groups, row_labels, column_labels) = _global_restriction_matrix(precision, [vector((1,0))], 12)
    #     sage: global_restriction_matrix
    #     [1 0 0 0 0 0 0 0 0]
    #     [0 1 2 0 0 0 0 0 0]
    #     [2 0 2 0 0 0 0 0 0]
    #     [0 0 2 1 2 0 0 0 0]
    #     [0 2 0 0 2 0 0 0 0]
    #     [2 0 0 0 2 1 2 0 0]
    #     [0 0 2 2 0 0 2 0 0]
    #     [0 2 0 0 0 0 2 1 2]
    #     [0 0 0 0 2 2 0 0 2]
    #     sage: (row_groups, row_labels, column_labels)
    #     ([((1, 0), 1, 0, 9)], {1: {(0, 0): 0, (3, 0): 5, (3, 1): 6, (2, 1): 4, (2, 0): 3, (1, 0): 1, (4, 1): 8, (1, 1): 2, (4, 0): 7}}, [(0, (0, 0)), (1, (0, 0)), (1, (1, 1)), (2, (0, 0)), (2, (1, 1)), (3, (0, 0)), (3, (1, 1)), (4, (0, 0)), (4, (1, 1))])
    """
    k = k % 2

    rst_jacobi_indices = [ m(s) for s in rst_vectors ]
    rst_indices = dict( (m_rst, classical_jacobi_fe_indices(m, prec,
                                                        reduced = not find_relations))
                          for m_rst in Set(rst_jacobi_indices) )
    

    column_labels = higherrank_jacobi_fe_indices(m, prec, reduced=True)
    reductions = dict( (nr,[]) for nr in column_labels )
    for nr in higherrank_jacobi_fe_indices(m, prec, reduced=False):
        (nrred, sgn) = reduce_higherrank_jacobi_fe_index(nr, m, r_classes, m_span)
        reductions[nrred].append((nr, sgn))     

    row_groups = [ len(rst_indices[m_rst]) for m_rst in rst_jacobi_indices ]
    row_groups = [ (s, m_rst, sum(row_groups[:i]), row_groups[i])
                   for ((i, s), m_rst) in zip(enumerate(rst_vectors), rst_jacobi_indices) ]
    row_labels = dict( (m_rst, dict( (nr, i) for (i, nr) in enumerate(rst_indices[m_rst]) ))
                       for m_rst in Set(rst_jacobi_indices) )
    dot_products = [ cython_lambda(
        ' , '.join([ 'int x{}'.format(i) for i in range(len(s)) ]),
        ' + '.join([ '{}*x{}'.format(*si) for si in enumerate(s) ]) )
                     for (s, _, _, _) in row_groups ]
    
    restriction_matrix = zero_matrix(ZZ, row_groups[-1][2] + row_groups[-1][3],
                                     len(column_labels))
    
    for (col, nrred) in enumerate(column_labels):
        for ((n, r), sgn) in reductions[nrred]:
            for ((s, m_rst, start, length), dot_product) in zip(row_groups, dot_products):
                row_labels_dict = row_labels[m_rst]
                try :
                    rst_r = dot_product(*r)
                    restriction_matrix[start + row_labels_dict[(n, rst_r)], col] \
                      += 1 if k == 0 else sgn
                except KeyError :
                    pass

    return (restriction_matrix, row_groups, row_labels, column_labels)

def _relation_matrix(k, m, prec, rst_vectors, r_classes, m_span) :
    r"""
    Deduce relations of the coefficients of a Jacobi form from their
    specialization to Jacobi form of scalar index.
    
    INPUT:
    
    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.
    
    - ``rst_vectors`` -- A list of vectors.

    - ``r_classes`` -- INSERT

    - ``m_span`` -- INSERT
    """
    k = k % 2

    (mat, row_groups, row_labels, column_labels) = \
        _restriction_matrix(k, m, prec, rst_vectors, find_relations, r_classes, m_span)

    relations = list()
    for (s, m_rst, start, length) in row_groups :
        row_labels_dict = row_labels[m_rst]
        for (nr, i) in row_labels_dict.iteritems() :
            (nrred, sgn) = reduce_higherrank_jacobi_fe_index(nr, m, r_classes, m_span)
            if nrred == nr: continue
            
            relations.append(mat.row(start + row_labels_dict[nrred])
                             - (1 if k == 0 else sgn) * mat.row(start + i))

    return (matrix(len(relations), len(column_labels), relations), column_labels)

def _higherrank_jacobi_forms__restriction(
        k, prec, relation_prec, dim,
        restriction_matrix__big, row_groups, row_labels, column_labels,
        relation_matrix, column_labels_relations):
    r"""
    Compute the Fourier expansions of Jacobi forms (over `\QQ`) of weight `k` and 
    index `m` up to given precision.

    INPUT:
    
    - `k` -- An integer.  Only `k` modulo `2` is used.

    - `m` -- A quadratic form.

    - ``prec`` -- A nonnegative integer.

    - ``relation_prec`` -- A nonnegative integer.

    - ``dim`` -- The dimension of the space of Jacobi forms.
    
    OUTPUT:
    
    - A list of dictionaries.
    
    # TESTS::
    
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
    #     sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _coefficient_by_restriction
    #     sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
    #     sage: precision = indices.filter(20)
    #     sage: relation_precision = indices.filter(10)
    #     sage: _coefficient_by_restriction(precision, 10) == _coefficient_by_restriction(precision, 10, relation_precision) 
    #     True
    """
    assert relation_prec <= prec

    ## Construct a restriction matrix for relation_prec
    row_groups__small = [ len(filter( lambda (n,r): n < relation_prec, row_labels[m_rst].keys()))
                          for (_, m_rst, _, _) in row_groups ]
    row_groups__small = [ (sum(row_groups__small[:i]), row_groups__small[i])
                          for i in range(len(row_groups)) ]
        
    row_labels__small = dict()
    for (m_rst, row_labels_dict) in row_labels.iteritems():
        row_labels_dict__small = dict()
            
        label_nmb = 0
        for (nr, i) in row_labels_dict.items():
            if nr[0] < relation_prec :
                row_labels_dict__small[nr] = (label_nmb, i)
                label_nmb += 1
            
        row_labels__small[m_rst] = row_labels_dict__small
        
    row_indices__small = list()
    for ((s, m_rst, start, _), (start_small, length_small)) in zip(row_groups, row_groups__small):
        row_labels_dict = row_labels__small[m_rst]
        row_indices__sub = length_small * [None]
        for (_,(i, i_pre)) in row_labels_dict.items():
            row_indices__sub[i] = start + i_pre
        row_indices__small += row_indices__sub
        
    restriction_matrix = restriction_matrix__big \
        .matrix_from_rows_and_columns(
            row_indices__small,
            [column_labels.index(nr) for nr in column_labels_relations] )
    
    
    rst_jacobi_indices = [ m_rst for (_, m_rst, _, _) in row_groups ]
    rst_indices = dict( (m_rst, classical_jacobi_fe_indices(m_rst, prec, reduced=True))
                        for m_rst in Set(rst_jacobi_indices) )
    rst_jacobi_forms = dict( (m_rst, classical_jacobi_forms(k, m_rst, prec)) )
                         for m_rst in Set(rst_jacobi_indices) )
    
    rst_jacobi_vectors = []
    nmb_rst_coords = row_groups[-1][2] + row_groups[-1][3]
    for (s, m_rst, start, length) in row_groups:
        row_labels_dict = row_labels[m_rst]
        for phi in rst_jacobi_forms[m_rst]:
            v = vector(ZZ, len(row_labels_dict))
            for (nr, i) in row_labels_dict.iter():
                (sgn, nrred) = reduce_classical_jacobi_fe_index(nr, m_rst)
                v[i] = ((1 if k%2 == 0 else sgn) * f[nr]) if nr in phi else 0
    
            rst_jacobi_vectors.append(vector(
                start*[0] + v.list() + (nmb_rst_coords - start - length)*[0] ))


    rst_jacobi_matrix__big = \
        matrix(len(forms), nmb_rst_coords, rst_jacobi_vectors).transpose()
    rst_jacobi_matrix = rst_jacobi_matrix__big.matrix_from_rows(row_indices__small)
    rst_jacobi_space = rst_jacobi_matrix.column_module() 
    
    rst_expansions_matrix = restriction_matrix.column_module() \
                                .intersection(rst_jacobi_space) \
                                .basis_matrix().transpose()
    rst_preimage = restriction_matrix.solve_right(rst_expansions_matrix).column_space() \
                   + restriction_matrix.right_kernel().change_ring(QQ)
    jacobi_expansions_space = rst_preimage.intersection(relation_matrix.right_kernel())
    

    if jacobi_expansions_space.dimension() < dim:
        raise RuntimeError( "There is a bug in the implementation of the restriction method. Dimensions: {}, {}".format(jacobi_expansions_space.dimension(), dim) )
    if jacobi_expansions_space.dimension() > dim :
        raise ValueError( "Could not construct enough restrictions to determine Fourier expansions uniquely", "INSUFFICIENT RELATIONS" )

    
    ## reconstruct the whole Fourier expansion from partial ones
    restriction_coordinates = rst_jacobi_matrix.solve_right(
        restriction_matrix * jacobi_expansions_space.basis_matrix().transpose() )
    jacobi_expansions__big = restriction_matrix__big.solve_right( rst_jacobi_matrix__big * restriction_coordinates ).change_ring(QQ)


    return [dict((nr,c) for (nr,c) in zip(column_labels, phi))
            for phi in jacobi_expansions__big.rows()]
