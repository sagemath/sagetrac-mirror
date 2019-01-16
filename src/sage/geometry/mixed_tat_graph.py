r"""
Mixed t\^ete-\`a-t\^ete graphs.





- Pablo Portilla (2019)
"""

#*****************************************************************************
#       Copyright (C) 2016 Pablo Portilla  <p.portilla89@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.unique_representation import CachedRepresentation
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.rings.integer_ring import ZZ
from sage.misc.cachefunc import cached_method
from sage.misc.flatten import flatten
from copy import copy
from sage.geometry.ribbon_graph import *
from sage.functions.other import floor
from sage.rings.rational import Rational
from sage.matrix.matrix_space import MatrixSpace
from sage.geometry.tat_graph import *

def _find(l, k):
    r"""
    Return the two coordinates of the element ``k`` in the list of
    lists ``l``.

    INPUT:

    - ``l`` -- a list of lists
    - ``k`` -- a candidate to be in a list in ``l``

    OUTPUT:

    A list with two integers describing the position of the first
    instance of `k`` in ``l``.

    TESTS::

        sage: from sage.geometry.ribbon_graph import _find
        sage: A = [[2,3,4],[4,5,2],[8,7]]
        sage: _find(A,2)
        [0, 0]
        sage: _find(A,7)
        [2, 1]
        sage: _find(A,5)
        [1, 1]
        sage: _find(A,-1)
        Traceback (most recent call last):
        ...
        ValueError: element -1 not found
    """
    for i,lst in enumerate(l):
        if k in lst:
            return [i, lst.index(k)]
    raise ValueError("element {} not found".format(k))


class MixedTatGraph(SageObject):
    r"""


    INPUT:

    - ``list_ribbon`` -- a list of nested ribbon graphs.
    - ``metric`` -- a list of as many rational numbers as darts has the first 
      element of``ribbon``.
    - ``relative_boundary=[]`` -- a subset of ``ribbon.boundary()`` that 
      constitutes the relative boundary components of ``ribbon``. It is, by
      default, initialized to an empty list (for defining pure t\^ete-\`a-t\^ete graphs)

    Alternatively, one can pass in 2 integers and this will construct
    a bipartite t\^ete-\`a-t\^ete graph which realizes the corresponding 
    Brieskorn-Pham singularity.
    
    
    EXAMPLES:
    
        sage: T33 = bipartite_tat_graph(3,3); T33
        Tete-a-tete graph of order 3 on a ribbon graph of genus 1 and 3 boundary components.
        sage: T33.sigma()
        (1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)
        sage: T33.rho()
        (1,12)(2,15)(3,18)(4,11)(5,14)(6,17)(7,10)(8,13)(9,16)
        sage: print(T33.metric())
        [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]
        sage: T33._relative_boundary
        []
        sage: B33 = blow_up(T33, 1,1/8); B33 ; B33._relative_boundary
        Relative t\^ete-\`a-t\^ete graph of order 3 on a ribbon graph of genus 1 and 6 boundary components; where 3 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
        [[19, 24, 23, 22, 21, 20], [25, 30, 29, 28, 27, 26], [31, 36, 35, 34, 33, 32]]
        sage: T517=bipartite_tat_graph(5,17);T517
        Tete-a-tete graph of order 85 on a ribbon graph of genus 32 and 1 boundary components.
    """
    def __init__(self, list_ribbon, metric, relative_boundary=[]):
        r"""
        Initialize ``self``.
        """
        for i in range(len(relative_boundary)):
            assert relative_boundary[i] in ribbon.boundary()

        assert check_tat_property(ribbon, metric, relative_boundary) == True
        self._ribbon = ribbon
        self._metric = metric
        self._sigma = ribbon._sigma
        self._rho = ribbon._rho
        self._boundary = self._ribbon.boundary()
        self._relative_boundary = relative_boundary
        self._basis = self._ribbon.homology_basis()
        self._mu = 2*self._ribbon.genus() + self._ribbon.number_boundaries()-1


    def _repr_(self):
        r"""
        Return basic information about the mixed t\^ete-\`a-t\^ete graph in string
        format.

        EXAMPLES:

        Example of a relative t\^ete-\`a-t\^ete graph::

            sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
            sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)')
            sage: R0 = RibbonGraph(s0,r0)
            sage: m0 = 6*[1/8,3/8,1/8]; 
            sage: perm_bound = [[1, 9, 7, 6, 4, 3], [10, 18, 16, 15, 13, 12]]
            sage: T0 = TatGraph(R0,m0,relative_boundary = perm_bound); T0
            Relative t\^ete-\`a-t\^ete graph of order 6 on a ribbon graph of genus 1 and 3 boundary components; where 2 boundary components are part of the relative boundary and might be permuted by the automorphism induced.

        Example of a pure t\^ete-\`a-t\^ete graph::

            sage: T23 = bipartite_tat_graph(2,3); T23
            Tete-a-tete graph of order 6 on a ribbon graph of genus 1 and 1 boundary components.
        """
        if not self._relative_boundary:
            return "Tete-a-tete graph of order {} on a ribbon graph of genus {} and {} boundary components.".format(self.order(), self._ribbon.genus(), self._ribbon.number_boundaries())
        else:
            return "Relative t\^ete-\`a-t\^ete graph of order {} on a ribbon graph of genus {} and {} boundary components; where {} boundary components are part of the relative boundary and might be permuted by the automorphism induced.".format(self.order(), self._ribbon.genus(), self._ribbon.number_boundaries(), len(self._relative_boundary))

    def sigma(self):
        r"""
        Return the permutation `\sigma` of ``self._ribbon``.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R = RibbonGraph(s1, r1)
            sage: m = 6*[1/2]
            sage: T = TatGraph(R,m)
            sage: T.sigma()
            (1,3,5)(2,4,6)
        """
        return self._ribbon.sigma()

    def rho(self):
        r"""
        Return the permutation `\rho` of ``self._ribbon``.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R = RibbonGraph(s1, r1)
            sage: m = 6*[1/2]
            sage: T = TatGraph(R,m)
            sage: T.rho()
            (1,2)(3,4)(5,6)
        """
        return self._ribbon.rho()

    def metric(self):
        r"""
        Return a vector containing the metric of the graph where the `i`th
        (starting at `0`) value of the vector corresponds to the dart `i+1`.

        EXAMPLES::

            sage: T23 = bipartite_tat_graph(2,3); T23.metric()
            [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]
        """
        return self._metric

    def ribbon(self):
        r"""
        Return the underlying ribbon graph of the t\^ete-\`a-t\^ete graph ``self``.

        EXAMPLES::

            sage: T23 = bipartite_tat_graph(2,3); T23.ribbon(); T23.ribbon().sigma(); T23.ribbon().rho()
            Ribbon graph of genus 1 and 1 boundary components
            (1,2,3)(4,5,6)(7,8)(9,10)(11,12)
            (1,8)(2,10)(3,12)(4,7)(5,9)(6,11)
        """
        return self._ribbon


    def action_homology(self):
        r"""
        Return matrix representing the action of the t\^ete-\`a-t\^ete automorphism 
        on the first homology group.

        OUTPUT:

        - Return a `self._mu() \times self._mu()` matrix that represents
          the action of the t\^ete-\`a-t\^ete automorphism on the first homology group
          with respect to the basis self._ribbon.homology_basis(). This 
          matrix has only `0`, `1` and `-1` as entries.

        EXAMPLES::

            sage: T34 = bipartite_tat_graph(3,4); T34
            Tete-a-tete graph of order 12 on a ribbon graph of genus 3 and 1 boundary components.
            sage: T34._ribbon.homology_basis()
            [[[6, 17], [18, 2], [1, 15], [14, 5]],
            [[7, 20], [21, 3], [1, 15], [14, 5]],
            [[8, 23], [24, 4], [1, 15], [14, 5]],
            [[10, 16], [18, 2], [1, 15], [13, 9]],
            [[11, 19], [21, 3], [1, 15], [13, 9]],
            [[12, 22], [24, 4], [1, 15], [13, 9]]]
            sage: T34.action_homology(); T34.order(); T34.action_homology()**12
            [ 0  0  0  1 -1  0]
            [ 0  0  0  1  0 -1]
            [ 0  0  0  1  0  0]
            [-1  1  0  1 -1  0]
            [-1  0  1  1  0 -1]
            [-1  0  0  1  0  0]
            12
            [1 0 0 0 0 0]
            [0 1 0 0 0 0]
            [0 0 1 0 0 0]
            [0 0 0 1 0 0]
            [0 0 0 0 1 0]
            [0 0 0 0 0 1]

            sage: T33 = bipartite_tat_graph(3,3); T33
            Tete-a-tete graph of order 3 on a ribbon graph of genus 1 and 3 boundary components.
            sage: T33.action_homology(); T33.order(); T33.action_homology()**3
            [ 0  0  1 -1]
            [ 0  0  1  0]
            [-1  1  1 -1]
            [-1  0  1  0]
            3
            [1 0 0 0]
            [0 1 0 0]
            [0 0 1 0]
            [0 0 0 1]
            sage: BT33 = blow_up(T33, 0, 1/8)
            sage: BT33.action_homology(); BT33.action_homology()**3
            [ 0 -1  1  1  0  0  0]
            [ 0  0  1  0  0 -1  0]
            [ 0  0  0 -1  0  0  1]
            [ 0  0  1  1 -1  0  0]
            [ 0  0  1  0  0  0  0]
            [ 1  0  0 -1  0  0  0]
            [ 0  0  1  1  0  0  0]
            [1 0 0 0 0 0 0]
            [0 1 0 0 0 0 0]
            [0 0 1 0 0 0 0]
            [0 0 0 1 0 0 0]
            [0 0 0 0 1 0 0]
            [0 0 0 0 0 1 0]
            [0 0 0 0 0 0 1]

            
        """
        #Set the space of matrices that we will be working on and 
        #a copy  of the zero matrix that will be modified to get 
        #the action of the monodromy.
        M = MatrixSpace(ZZ, self._mu, self._mu) 
        T = copy(M.zero_matrix())

        #we run along all the darts of each of the elements of the 
        #basis of the homology. We check which is their image; if they
        #intersect with positive orientation (with respect to the fixed
        #orientation of homology_basis) then we set a +1 in the corresponding
        #column, if they intersect with negative orientation, a -1. 
        for i in range (self._mu):
            for j in range (self._mu):
                 for k in range (len(self._basis[i])):
                    if (safewalk(self._ribbon, 
                                 self._metric, 
                                 self._basis[i][k][0],
                                 self._relative_boundary) 
                        == self._basis[j][0][0]):
                        T[i,j] = 1
                    elif (safewalk(self._ribbon, 
                                   self._metric, 
                                   self._basis[i][k][0],
                                   self._relative_boundary) 
                          == self._basis[j][0][1]):

                        T[i,j] = -1
        return T


