r"""
Tests for higherrank.py.

AUTHOR:

- Martin Raum
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

from sage.modules.all import FreeModule
from sage.quadratic_forms.all import QuadraticForm
from sage.rings.all import ZZ


from sage.modular.jacobi.higherrank import (
    reduce_higher_rank_jacobi_fe_index
)


def _test_set__jacobi_m():
    return [QuadraticForm(matrix(1, [[2]])),
            QuadraticForm(matrix(2, [2, 0, 0, 2])),
            QuadraticForm(matrix(2, [2, 1, 1, 2])),
            QuadraticForm(matrix(3, [2,1,1, 1,2,1, 1,1,2]))]


def test__reduce_higher_rank_jacobi_fe_index():
    for m in _test_set__jacobi_m():
        r_classes = _higherrank_jacobi_r_classes(m)[0]
        m_span = m.matrix().row_module()
        fm = FreeModule(ZZ, m.dim())

        for _ in range(10):
            n = random_int()
            for _ in range(10):
                r = fm.random_element()
                yield (_test__reduce_higher_rank_jacobi_fe_index, ((n, r), m, r_classes, m_span))

def _test__reduce_higher_rank_jacobi_fe_index((n, r), m, r_classes, m_span):
    return NotImplemented()

def test__higherrank_jacobi_r_classes():
    for m in _test_set__jacobi_m():
        yield (_test__higherrank_jacobi_r_classes, (m,))

def _test__higherrank_jacobi_r_classes(m):
    (r_classes, r_classes_reduction_signs) = _higherrank_jacobi_r_classes(m)
    m_span = m.matrix().span()

    ## cardinality test
    assert (m.det()
            ==
            sum((2 if -1 in signs else 1) for signs in r_classes_reduction_signs)
    )

    ## equivalence in classes
    for (r_class, r_signs) in zip(r_classes, r_classes_reduction_signs):
        r0 = vector(r_class[0])
        for (r,s) in zip(map(vector, r_class), r_signs):
            assert (r - s*r0
                    in
                    m_span
            )

    ## distinction
    for ix in range(len(r_classes)):
        for r_class in r_classes[ix+1:]:
            assert (vector(r_classes[ix][0]) - vector(r_class[0]) not in m_span
                    and
                    vector(r_classes[ix][0]) + vector(r_class[0]) not in m_span
            )

# def _test__coefficient_by_restriction(precision, k, relation_precision = None, additional_lengths = 1 ) :
#     r"""
#     TESTS::
    
#         sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
#         sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _test__coefficient_by_restriction

#     ::

#         sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
#         sage: precision = indices.filter(10)
#         sage: _test__coefficient_by_restriction(precision, 10, additional_lengths = 10)
        
#     ::

#         sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [4,1,1,2])))
#         sage: precision = JacobiFormD1Filter(5, indices.jacobi_index())
#         sage: _test__coefficient_by_restriction(precision, 40, additional_lengths = 4) # long test
        
#     We use different precisions for relations and restrictions::

#         sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
#         sage: precision = indices.filter(20)
#         sage: relation_precision = indices.filter(2)
#         sage: _test__coefficient_by_restriction(precision, 10, relation_precision, additional_lengths = 4)
#     """
#     from sage.misc.misc import verbose
    
#     L = precision.jacobi_index()
    
#     if relation_precision is not None and not relation_precision <= precision :
#         raise ValueError( "Relation precision must be less than or equal to precision." )

#     expansions = _coefficient_by_restriction(precision, k, relation_precision)
#     verbose( "Start testing restrictions of {2} Jacobi forms of weight {0} and index {1}".format(k, L, len(expansions)) )
    
#     ch1 = JacobiFormD1WeightCharacter(k)
#     chL = JacobiFormD1WeightCharacter(k, L.matrix().nrows())
    
#     R = precision.monoid()._r_representatives
#     S_extended = _find_complete_set_of_restriction_vectors(L, R)

#     S = list()
#     for (s, _) in S_extended :
#         if s not in S :
#             S.append(s) 
#     max_S_length = max([L(s) for s in S])

#     Snew = flatten( enumerate_short_vectors__python( map(list, L.matrix().rows()), 2 * max_S_length + 2, max_S_length + additional_length, True).values(), max_level = 1 )
#     verbose( "Will use the following restriction vectors: {0}".format(Snew) )
    
#     jacobi_forms_dict = dict()
#     non_zero_expansions = list()
#     for s in Snew :
#         m = L(s)
#         verbose( "Restriction to index {0} via {1}".format(m, s) )
        
#         try :
#             jacobi_forms = jacobi_forms_dict[m]
#         except KeyError : 
#             jacobi_forms = JacobiFormsD1NN(QQ, JacobiFormD1NNGamma(k, m), JacobiFormD1NNFilter(precision.index(), m))
#             jacobi_forms_dict[m] = jacobi_forms
#         jacobi_forms_module = span([ vector( b[(ch1, k)] for k in jacobi_forms.fourier_expansion_precision() )
#                                      for b in map(lambda b: b.fourier_expansion(), jacobi_forms.graded_submodule(None).basis()) ])
        
#         fourier_expansion_module = jacobi_forms.fourier_expansion_ambient()
        
#         for (i, expansion) in enumerate(expansions) :
#             verbose( "Testing restriction of {0}-th form".format(i) )
#             restricted_expansion_dict = dict()
#             for (n,r) in precision.monoid_filter() :
#                 rres = s.dot_product(vector(r))
#                 try :
#                     restricted_expansion_dict[(n,rres)] += expansion[(chL,(n,r))]
#                 except KeyError :
#                     restricted_expansion_dict[(n,rres)] = expansion[(chL,(n,r))]
            
#             restricted_expansion = vector( restricted_expansion_dict.get(k, 0) for k in jacobi_forms.fourier_expansion_precision() )
#             if restricted_expansion not in jacobi_forms_module :
#                 raise RuntimeError( "{0}-th restricted via {1} is not a Jacobi form".format(i, s) )
            
#             if restricted_expansion != 0 :
#                 non_zero_expansions.append(i)
    
#     assert Set(non_zero_expansions) == Set(range(len(expansions)))
 
# def _local_restriction_matrix(r_classes, r_classes_reduction_signs, rst_vectors) :
#     r"""
#     Return a matrix whose rows correspond to the evaluations of the restriction
#     vectors (s, rst_r) in rst_vectors.
    
#     INPUT:
    
#     - ``r_classes`` -- A list of lists of tuples or vectors in L \otimes QQ (with
#                        given coordinates).  

#     - ``r_classes_reduction_signs`` -- A list of lists of `\pm 1`.
    
#     - ``rst_vectors`` -- A list of pairs `(s, r)`, where `s` is a vector,
#                          and `r` is an integer.
             
#     OUTPUT:
    
#     - A matrix with integer entries.
    
#     TESTS::

#         sage: from psage.modform.jacobiforms.jacobiformd1_fourierexpansion import *
#         sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _find_complete_set_of_restriction_vectors
#         sage: from psage.modform.jacobiforms.jacobiformd1_fegenerators import _local_restriction_matrix                
#         sage: indices = JacobiFormD1Indices(QuadraticForm(matrix(2, [2,1,1,2])))
#         sage: R = indices._r_representatives
#         sage: S = _find_complete_set_of_restriction_vectors(indices.jacobi_index(), R, 4)        
#         sage: _local_restriction_matrix(R, S)
#         [1 1 1]
#         [0 1 1]
#         [0 1 1]
#         [1 1 1]
#         [0 1 1]
#         [0 1 1]
#         [0 0 2]
#     """
#     R = [map(vector, rs) for rs in R]
    
#     return matrix([ _eval_restriction_vector(R, vector(s), r, reduction_function) for (s, r) in S ])
