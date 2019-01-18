# -*- coding: utf-8 -*-
"""
Tête-à-tête graphs



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
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.rings.integer_ring import ZZ
from sage.misc.flatten import flatten
from copy import copy
from sage.geometry.ribbon_graph import RibbonGraph
from sage.geometry.ribbon_graph import bipartite_ribbon_graph
from sage.rings.rational import Rational
from sage.matrix.matrix_space import MatrixSpace

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


def safewalk(ribbon_graph, metric, edge, relative_boundary = [], sf_length = 1):
    r"""
    Return the endpoint of the safe walk starting at 'edge'.

    INPUT:

    - ``ribbon_graph`` -- a ribbon graph.
    - ``metric`` --  a vector of as many positive rational numbers as
      darts has the graph.
    - ``edge`` -- is a integer from 1 to the number of darts.
    - ``relative_boundary`` -- a list of lists containing the permuted boundary
      components. It is empty by default if the ribbon graph is not
      relative.
    - ``sf_length`` -- the length of the safewalk. It is 1 by default.

    OUTPUT:

    - a number from 1 to the number of darts indicating the position of 
      the endpoint of the safewalk. It returns 0 if the image of the dart
      contains a vertex in its interior (in this case, the tat property
      can't be satisfied.

    EXAMPLES::

        sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
        sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)')
        sage: R0 = RibbonGraph(s0,r0); R0; R0.genus(); R0.number_boundaries(); print(R0.boundary())
        Ribbon graph of genus 1 and 3 boundary components
        1
        3
        [[1, 9, 7, 6, 4, 3], [2, 11, 12, 13, 14, 5, 6, 7, 8, 17, 18, 10, 11, 2, 3, 4, 5, 14, 15, 16, 17, 8, 9, 1], [10, 18, 16, 15, 13, 12]]
        sage: m0 = 6*[1/8,3/8,1/8]; print(m0)
        [1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8]
        sage: perm_bound = [[1,9,7,6,4,3], [10,18,16,15,13,12]]
        sage: safewalk(R0,m0,2,relative_boundary = perm_bound)
        14
        sage: safewalk(R0,m0,14,relative_boundary = perm_bound)
        8
        sage: safewalk(R0,m0,2,relative_boundary = perm_bound, sf_length= 2)
        8
        sage: safewalk(R0,m0,9,relative_boundary = perm_bound)
        12
    """
    #First we check if the given dart is in the relative boundary
    rel_b = edge in flatten(relative_boundary)
    end_point = edge

    #if the edge in not in the relative boundary we execute one algorithm
    #that takes the safewalk that starts by applying rho.
    if not rel_b:
            while (sf_length > 0):
                sf_length = sf_length - metric[ribbon_graph._rho(end_point)-1] \
                    - metric[ribbon_graph._sigma(ribbon_graph._rho(end_point))-1]
                end_point = ribbon_graph._sigma(ribbon_graph._rho(end_point))
            if (sf_length < 0):
                return 0
            elif sf_length == 0:
                return end_point

    #if the edge is in the relative boundary we execute another algorithm
    #that returns the end of the boundary safewalk.
    elif rel_b:
        aux_ind = _find(relative_boundary, edge)
        if (ribbon_graph._rho(edge) == \
           relative_boundary[aux_ind[0]][(aux_ind[1]-1) % 
                                         len(relative_boundary[aux_ind[0]])]):
            while (sf_length > 0):
                sf_length = sf_length - metric[ribbon_graph._rho(end_point)-1] \
                    - metric[ribbon_graph._sigma(ribbon_graph._rho(end_point))-1]
                end_point = ribbon_graph._sigma(ribbon_graph._rho(end_point))
            if (sf_length < 0):
                return 0
            elif sf_length == 0:
                return end_point
        elif (ribbon_graph._rho(edge) == \
           relative_boundary[aux_ind[0]][(aux_ind[1]+1) % 
                                         len(relative_boundary[aux_ind[0]])]):
            while (sf_length > 0):
                sf_length = sf_length - metric[ribbon_graph._sigma(end_point)-1] \
                    - metric[ribbon_graph._rho(ribbon_graph._sigma(end_point))-1]
                end_point = ribbon_graph._rho(ribbon_graph._sigma(end_point))
            if (sf_length < 0):
                return 0
            elif sf_length == 0:
                return end_point




def check_tat_property(ribbon_graph, metric, relative_boundary = [], sf_length=1):
    r"""
    Check wether the tat property is satisfied for the given
    ribbon graph and metric or not.

    INPUT:

    - ``ribbon_graph`` -- a ribbon graph.
    - ``metric`` --  a vector of as many positive rational numbers as
      darts has the graph.
    - ``relative_boundary`` -- a list of lists containing the permuted boundary
      components. It is empty by default if the ribbon graph is not
      relative.
    - ``sf_length`` -- a positive rational number. It is the length of the safe
      walks for which the tête-à-tête property has to be checked. By default, it
      is initialized to be `1`.

    OUTPUT:

    - True if the tat property is satisfied for the given
      ribbon graph and metric and False if not.

    EXAMPLES::

        sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
        sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)');
        sage: R0 = RibbonGraph(s0,r0); R0; R0.genus(); R0.number_boundaries(); print(R0.boundary());
        Ribbon graph of genus 1 and 3 boundary components
        1
        3
        [[1, 9, 7, 6, 4, 3], [2, 11, 12, 13, 14, 5, 6, 7, 8, 17, 18, 10, 11, 2, 3, 4, 5, 14, 15, 16, 17, 8, 9, 1], [10, 18, 16, 15, 13, 12]]
        sage: m0 = 6*[1/8,3/8,1/8]; print(m0);
        [1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8]
        sage: perm_bound = [[1,9,7,6,4,3], [10,18,16,15,13,12]]
        sage: check_tat_property(R0,m0,relative_boundary = perm_bound)
        True
        sage: m0[1]=2/3; print(m0);
        [1/8, 2/3, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8]
        sage: check_tat_property(R0,m0,relative_boundary = perm_bound)
        False
    """
    #auxiliary variable prop that holds if the tat property is true
    #or false.
    prop = True
    i = 0
    #store the length of the safe walks of the ribbon graph
    sl = sf_length
    #we run a while that checks the tat property on each edge (i.e. each
    #element of rho.
    while (i < len(ribbon_graph._rho.cycle_tuples()) and 
           prop):
        aux_edge = ribbon_graph._rho.cycle_tuples()[i]
        #see definition of safewalk to see that it returns 0
        #if at some point the remaining length to walk is negative.
        if safewalk(ribbon_graph,
                    metric,
                    aux_edge[0],
                    relative_boundary, sl)==0:
            prop = False
        #if the length is exactly 0 at some point, it means that the 
        #method safewalk ends at a dart, then this checks the tat
        #property.
        elif safewalk(ribbon_graph, 
                      metric, 
                      aux_edge[0],
                      relative_boundary,sl)!=0:
            prop = (ribbon_graph._rho(safewalk(ribbon_graph, 
                                          metric, 
                                          aux_edge[0],
                                          relative_boundary,sl)
                                    ) == \
                    safewalk(ribbon_graph, metric, aux_edge[1], 
                             relative_boundary,sl)
                   )
        i+=1

    return prop


def bipartite_tat_graph(p,q):
    r"""
    Return the bipartite_ribbon_graph(p,q) with a metric such that it models
    the Milnor fiber and the monodromy of the Brieskorn-Pham singularity
    `x^p+y^q`.

    INPUT:

    - ``p`` -- a positive integer.
    - ``q`` -- a positive integer.

    OUTPUT:

    - a tête-à-tête graph whose underlying ribbon graph is the complete
      bipartite graph of type `(p,q)` and with all darts of length 
      `1/4`. It models the Milnor fiber and monodromy of the Brieskorn-Pham
      singularity `x^p+y^q`.

    EXAMPLES::

        sage: T23 = bipartite_tat_graph(2,3); T32 = bipartite_tat_graph(3,2); T32; T32.sigma(); T32.rho(); T32.metric();
        Tete-a-tete graph of order 6 on a ribbon graph of genus 1 and 1 boundary components.
        (1,2)(3,4)(5,6)(7,8,9)(10,11,12)
        (1,9)(2,12)(3,8)(4,11)(5,7)(6,10)
        [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]
        sage: T23.ribbon().genus() ==  T32.ribbon().genus(); T23.ribbon().number_boundaries() ==  T32.ribbon().number_boundaries()
        True
        True
        sage: T32.action_homology()
        [ 0  1]
        [-1  1]
        sage: T23.action_homology()
        [ 1 -1]
        [ 1  0]

        sage: T46 = bipartite_tat_graph(4,6); T46; T46.sigma(); T46.rho(); print(T46.metric())
        Tete-a-tete graph of order 12 on a ribbon graph of genus 7 and 2 boundary components.
        (1,2,3,4,5,6)(7,8,9,10,11,12)(13,14,15,16,17,18)(19,20,21,22,23,24)(25,26,27,28)(29,30,31,32)(33,34,35,36)(37,38,39,40)(41,42,43,44)(45,46,47,48)
        (1,28)(2,32)(3,36)(4,40)(5,44)(6,48)(7,27)(8,31)(9,35)(10,39)(11,43)(12,47)(13,26)(14,30)(15,34)(16,38)(17,42)(18,46)(19,25)(20,29)(21,33)(22,37)(23,41)(24,45)
        [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]
    """
    metric = 2*p*q*[Rational('1/4')]
    return TatGraph(bipartite_ribbon_graph(p,q), metric)

class TatGraph(SageObject):
    r"""
    A tête-à-tête  graph codified as a ribbon graph and a list of rational 
    numbers.

    The original idea is due to A'Campo [AC2009]_. His original motivation was
    to model the monodromies of isolated plane curves in a combinatorial way.
    
    The theory  was further developed by C. Graf in [Gra2015]_ and later on by J. Fernandez de Bobadilla, M. Pe Pereira and P. Portilla in [BPP2017]_ .

    **Introduction**

    Let `\Sigma` be an orientable surface with non-empty boundary and let 
    `\Gamma` be the topological realization of a graph (a finite `1`-dimensional
    CW-cmplex) that is embedded in `\Sigma` in such a way that the graph is a 
    strong deformation retract of the surface. Suppose that `\Gamma` is endowed
    with a metric `m`, that is, an assignation of a rational number to each edge.

    A safe walk on `\Gamma` is an arc-length parametrized path starting at a 
    point `p \in \Gamma` with the properties:
    
    - When the path reaches  a vertex, it turns to the next edge indicated by 
      the permutation `\sigma` associated to that vertex (see Ribbon Graph
      documentation)
      
    - It has total length equal to `1`.
    
    Given a point `p \in \Gamma` in the interior of and edge, there exist
    exactly two safe walks starting at `p`. The tête-à-tête property says
    that "for all points in the interior of `\Gamma`, the two safe walks 
    starting at that point, end at the same point.
    
    As explained in the references that are mentioned above, a tête-à-tête
    graph models an oriented surface with boundary together with a mapping
    class of the surface which is freely periodic and has positive fractional
    Dehn twist coefficients at all boundary components. A quich summary of how 
    to build the automorphism is the following:
    
    (1) Cut the surface along the graph `\Gamma` to obtain a disjoint collection
        of cylinders.

    (2) Perform a boundary Dehn twist of length `1` on each of the cylinders
        around the boundary component that comes from cutting along the graph.
    
    (3) The tête-à-tête property exactly tells you that the boundary Dehn twists
        performed in step (2) are compatible with the gluing map that recovers
        `\Sigma`.

    This proccess defines an element in the mapping class group of the surface
    relative to the boundary.
    
    Another object which is also modeled in this package is a relative
    tête-à-tête graph. These are metric ribbon graphs with a marked set of 
    circles in the graph which correspond to certain boundary components. While
    in pure tête-à-tête, all boundary components are left invariant by
    the tête-à-tête automorphism, in their relative counterpart, the relative
    boundary components might be permuted by the relative tête-à-tête 
    mondromy.

    For a future implementation of mixed tête-à-tête graphs we include the 
    possibility of changing the length of the safe walk.
    INPUT:

    - ``ribbon`` -- a ribbon graph.
    - ``metric`` -- a list of as many rational numbers as darts has ``ribbon``.
    - ``relative_boundary=[]`` -- a subset of ``ribbon.boundary()`` that 
      constitutes the relative boundary components of ``ribbon``. It is, by
      default, initialized to an empty list (for defining pure tête-à-tête)
    - ``sf_length`` -- the length of the safe walks of the tête-à-tête graph. It
      is `1` by default.


    EXAMPLES::

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
        Relative tête-à-tête graph of order 3 on a ribbon graph of genus 1 and 6 boundary components; where 3 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
        [[19, 24, 23, 22, 21, 20], [25, 30, 29, 28, 27, 26], [31, 36, 35, 34, 33, 32]]
        sage: T517=bipartite_tat_graph(5,17);T517
        Tete-a-tete graph of order 85 on a ribbon graph of genus 32 and 1 boundary components.

        sage: T = bipartite_tat_graph(2,5)
        sage: T
        Tete-a-tete graph of order 10 on a ribbon graph of genus 2 and 1 boundary components.
        sage: B = blow_up(T, 2, 1/7)
        sage: B
        Relative tête-à-tête graph of order 10 on a ribbon graph of genus 2 and 6 boundary components; where 5 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
        sage: R = T.ribbon();R
        Ribbon graph of genus 2 and 1 boundary components
        sage: m = 20*[1]
        sage: S = TatGraph(R,m, sf_length = 4)
        sage: S
        Tete-a-tete graph of order 10 on a ribbon graph of genus 2 and 1 boundary components.
        sage: K=blow_up(S,2,1/2)
        sage: K
        Relative tête-à-tête graph of order 10 on a ribbon graph of genus 2 and 6 boundary components; where 5 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
        sage: K._sf_length
        4
        sage: K.action_homology()
        [ 0 -1  1  0  0  0  0  0  0]
        [ 0  0  1  0 -1  0  0  0  0]
        [ 0  0  1 -1  0  0  0  0  0]
        [ 0  0  1  0  0  0 -1  0  0]
        [ 0  0  1  0  0 -1  0  0  0]
        [ 0  0  1  0  0  0  0  0 -1]
        [ 0  0  1  0  0  0  0 -1  0]
        [-1  0  1  0  0  0  0  0  0]
        [ 0  0  1  0  0  0  0  0  0]
        sage: K.action_homology()^10
        [1 0 0 0 0 0 0 0 0]
        [0 1 0 0 0 0 0 0 0]
        [0 0 1 0 0 0 0 0 0]
        [0 0 0 1 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 1 0 0 0]
        [0 0 0 0 0 0 1 0 0]
        [0 0 0 0 0 0 0 1 0]
        [0 0 0 0 0 0 0 0 1]

    """
    def __init__(self, ribbon, metric, relative_boundary=[], sf_length=1):
        r"""
        Initialize ``self``.

        EXAMPLES::

            Example of a manual definition of a relative tete-a-tete graph:

            sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
            sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)');
            sage: R0 = RibbonGraph(s0,r0); R0; R0.genus(); R0.number_boundaries(); print(R0.boundary());
            Ribbon graph of genus 1 and 3 boundary components
            1
            3
            [[1, 9, 7, 6, 4, 3], [2, 11, 12, 13, 14, 5, 6, 7, 8, 17, 18, 10, 11, 2, 3, 4, 5, 14, 15, 16, 17, 8, 9, 1], [10, 18, 16, 15, 13, 12]]

            sage: m0 = 6*[1/8,3/8,1/8]; print(m0);
            [1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8, 1/8, 3/8, 1/8]
            sage: perm_bound = [[1,9,7,6,4,3], [10,18,16,15,13,12]]
            sage: T = TatGraph(R0,m0, relative_boundary = perm_bound);T;
            Relative tête-à-tête graph of order 6 on a ribbon graph of genus 1 and 3 boundary components; where 2 boundary components are part of the relative boundary and might be permuted by the automorphism induced.

            However, if we change the perm_bound to something that is not in the boundary, it raises an AssertionError:

            sage: perm_bound[0]= [1,9,31]; perm_bound
            [[1, 9, 31], [10, 18, 16, 15, 13, 12]]
            sage: T = TatGraph(R0,m0, relative_boundary = perm_bound)
            Traceback (most recent call last):
            ...
            AssertionError
            sage: T = bipartite_tat_graph(2,7)
            sage: T
            Tete-a-tete graph of order 14 on a ribbon graph of genus 3 and 1 boundary components.
            sage: R = T.ribbon()
            sage: R
            Ribbon graph of genus 3 and 1 boundary components
            sage: m = 28*[1]
            sage: check_tat_property(R,m)
            False
            sage: S = TatGraph(R,m,sf_length=4)
            sage: S
            Tete-a-tete graph of order 14 on a ribbon graph of genus 3 and 1 boundary components.

        """
        for i in range(len(relative_boundary)):
            assert relative_boundary[i] in ribbon.boundary()

        assert sf_length > 0
        assert check_tat_property(ribbon, metric, 
                                  relative_boundary, sf_length)
        self._ribbon = ribbon
        self._metric = metric
        self._sigma = ribbon._sigma
        self._rho = ribbon._rho
        self._boundary = self._ribbon.boundary()
        self._relative_boundary = relative_boundary
        self._basis = self._ribbon.homology_basis()
        self._mu = 2*self._ribbon.genus() + self._ribbon.number_boundaries()-1
        self._sf_length = sf_length


    def _repr_(self):
        r"""
        Return basic information about the tête-à-tête graph in string
        format.

        EXAMPLES:

        Example of a relative tête-à-tête graph::

            sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
            sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)')
            sage: R0 = RibbonGraph(s0,r0)
            sage: m0 = 6*[1/8,3/8,1/8]; 
            sage: perm_bound = [[1, 9, 7, 6, 4, 3], [10, 18, 16, 15, 13, 12]]
            sage: T0 = TatGraph(R0,m0,relative_boundary = perm_bound); T0
            Relative tête-à-tête graph of order 6 on a ribbon graph of genus 1 and 3 boundary components; where 2 boundary components are part of the relative boundary and might be permuted by the automorphism induced.

        Example of a pure tête-à-tête graph::

            sage: T23 = bipartite_tat_graph(2,3); T23
            Tete-a-tete graph of order 6 on a ribbon graph of genus 1 and 1 boundary components.
        """
        if not self._relative_boundary:
            return "Tete-a-tete graph of order {} on a ribbon graph of genus {} and {} boundary components.".format(self.order(), self._ribbon.genus(), self._ribbon.number_boundaries())
        else:
            return "Relative tête-à-tête graph of order {} on a ribbon graph of genus {} and {} boundary components; where {} boundary components are part of the relative boundary and might be permuted by the automorphism induced.".format(self.order(), self._ribbon.genus(), self._ribbon.number_boundaries(), len(self._relative_boundary))

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
        Return a vector containing the metric of the graph where the `i` -th
        (starting at `0`) value of the vector corresponds to the dart `i+1`.

        EXAMPLES::

            sage: T23 = bipartite_tat_graph(2,3); T23.metric()
            [1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4, 1/4]
        """
        return self._metric

    def ribbon(self):
        r"""
        Return the underlying ribbon graph of the tête-à-tête graph ``self``.

        EXAMPLES::

            sage: T23 = bipartite_tat_graph(2,3); T23.ribbon(); T23.ribbon().sigma(); T23.ribbon().rho()
            Ribbon graph of genus 1 and 1 boundary components
            (1,2,3)(4,5,6)(7,8)(9,10)(11,12)
            (1,8)(2,10)(3,12)(4,7)(5,9)(6,11)
        """
        return self._ribbon

    def rot_numbers(self):
        r"""
        Return a list containing the rotation numbers of the tat 
        automorphism at each boundary component which is not a relative
        boundary component. The list comes in the same order as the list
        resulting from removing relative_boundary from self.boundary().

        By definition, these are signed rational numbers formed by an 
        integer which tells how many times the safe walk winds around
        that boundary component plus a rational number in `[0,1)` that
        tells the rotation of the tête-à-tête automorphism around that boundary 
        component. In the case of tête-à-tête automorphisms, the rotation 
        numbers are always positive.

        EXAMPLES:

        Consider the complete bipartite graph of type (2,3) and length 1/4 for all darts::

            sage: T23 = bipartite_tat_graph(2,3); T23
            Tete-a-tete graph of order 6 on a ribbon graph of genus 1 and 1 boundary components.
            sage: T23.rot_numbers()
            [1/6]

        If we divide its metric by 20, we get the 20th power of its
        automorphism. This is an automorphism of order `6/gcd(6,20)=3`
        and with rotation number the previous one multiplied by 20::

            sage: mod_metric = [x/20 for x in T23._metric]; mod_metric
            [1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80, 1/80]
            sage: T23_mod = TatGraph(T23._ribbon, mod_metric); T23_mod
            Tete-a-tete graph of order 3 on a ribbon graph of genus 1 and 1 boundary components.
            sage: T23_mod.rot_numbers()
            [10/3]

        Other examples::

            sage: T33 = bipartite_tat_graph(3,3); T33
            Tete-a-tete graph of order 3 on a ribbon graph of genus 1 and 3 boundary components.
            sage: T33.rot_numbers()
            [1/3, 1/3, 1/3]

            sage: s0 = PermutationGroupElement('(1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)')
            sage: r0 = PermutationGroupElement('(1,9)(2,11)(3,4)(5,14)(6,7)(8,17)(10,18)(12,13)(15,16)')
            sage: R0 = RibbonGraph(s0,r0)
            sage: m0 = 6*[1/8,3/8,1/8]; 
            sage: perm_bound = [[1, 9, 7, 6, 4, 3], [10, 18, 16, 15, 13, 12]]
            sage: T0 = TatGraph(R0,m0,relative_boundary = perm_bound); T0
            Relative tête-à-tête graph of order 6 on a ribbon graph of genus 1 and 3 boundary components; where 2 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
            sage: print(T0._ribbon.boundary()); T0.rot_numbers()
            [[1, 9, 7, 6, 4, 3], [2, 11, 12, 13, 14, 5, 6, 7, 8, 17, 18, 10, 11, 2, 3, 4, 5, 14, 15, 16, 17, 8, 9, 1], [10, 18, 16, 15, 13, 12]]
            [1/6]
        """
        rot = []
        for i in range (self._ribbon.number_boundaries()):
            aux_rot = 0
            if self._ribbon.boundary()[i] not in self._relative_boundary:
                for j in range (len(self._ribbon.boundary()[i])):
                    aux_rot += self._metric[self._ribbon.boundary()[i][j]-1]

                rot.append(self._sf_length / aux_rot)
        return rot

    def action_homology(self):
        r"""
        Return matrix representing the action of the tête-à-tête automorphism 
        on the first homology group.

        OUTPUT:

        - Return a ``self._mu()`` by ``self._mu()`` matrix that represents
          the action of the tête-à-tête automorphism on the first homology group
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
                                 self._relative_boundary,
                                 self._sf_length) 
                        == self._basis[j][0][0]):
                        T[i,j] = 1
                    elif (safewalk(self._ribbon, 
                                   self._metric, 
                                   self._basis[i][k][0],
                                   self._relative_boundary,
                                   self._sf_length) 
                          == self._basis[j][0][1]):

                        T[i,j] = -1
        return T

    def order(self):
        r"""
        Return the order of the tête-à-tête automorphism.

        OUTPUT:

        - A positive integer which is the order of the tête-à-tête automorphism.

        EXAMPLES::

            sage: s1 = PermutationGroupElement('(1,3,5)(2,4,6)')
            sage: r1 = PermutationGroupElement('(1,2)(3,4)(5,6)')
            sage: R = RibbonGraph(s1, r1)
            sage: m = 6*[1/2]
            sage: T= TatGraph(R, m); T; T.order()
            Tete-a-tete graph of order 6 on a ribbon graph of genus 1 and 1 boundary components.
            6
        """
        return self.rot_numbers()[0].denom()


    def orbit_graph(self):
        r"""
        Return the orbits of the tête-à-tête automorphism together with a
        ribbon graph which is the orbit of the ribbon graph
        by the action of the tête-à-tête automorphism.

        OUTPUT:

        - A pair consisting of a list called ``orbits`` where each
          instance is a list containing all the edges of an orbit of
          the tête-à-tête automorphism; and a ribbon graph called ``orbit_graph``
          which is the orbit space of the action induced by that tat
          automorphism. They satisfy that the dart ``k`` in orbit_graph
          corresponds to the orbit ``orbits[k-1]``. (Observe that the
          smallest dart is labeled as `1`, hence the shifting).


        EXAMPLES:

        First example::

            sage: T34 = bipartite_tat_graph(3,4); T34
            Tete-a-tete graph of order 12 on a ribbon graph of genus 3 and 1 boundary components.
            sage: orb_vector, orb_graph  = T34.orbit_graph(); orb_vector; orb_graph; orb_graph.sigma(); orb_graph.rho();
            [[1, 10, 7, 4, 9, 6, 3, 12, 5, 2, 11, 8],
             [13, 17, 21, 22, 14, 18, 19, 23, 15, 16, 20, 24]]
            Ribbon graph of genus 0 and 1 boundary components
            ()
            (1,2)


        We see that in the previous example ``orb_graph.sigma()`` produces
        the empty permutation. This is because that permutation is formed by two
        cycles of length `1` and sage doesn't  show these cycles by default. We
        can force that these are show::

            sage: T34 = bipartite_tat_graph(3,4); T34
            Tete-a-tete graph of order 12 on a ribbon graph of genus 3 and 1 boundary components.
            sage: orb_vector, orb_graph  = T34.orbit_graph(); orb_vector; orb_graph; orb_graph.sigma().cycle_tuples(singletons = 1); orb_graph.rho();
            [[1, 10, 7, 4, 9, 6, 3, 12, 5, 2, 11, 8],
             [13, 17, 21, 22, 14, 18, 19, 23, 15, 16, 20, 24]]
            Ribbon graph of genus 0 and 1 boundary components
            [(1,), (2,)]
            (1,2)

        Other examples::

            sage: T33 = bipartite_tat_graph(3,3); T33
            Tete-a-tete graph of order 3 on a ribbon graph of genus 1 and 3 boundary components.
            sage: orb_vector, orb_graph = T33.orbit_graph(); orb_vector; orb_graph;
            [[1, 8, 6], [2, 9, 4], [3, 7, 5], [10, 14, 18], [11, 15, 16], [12, 13, 17]]
            Ribbon graph of genus 0 and 3 boundary components
            sage: BT33 = blow_up(T33, 0, 1/8)
            sage: orb_vector_blow, orb_graph_blow = BT33.orbit_graph(); orb_vector_blow; orb_graph_blow; orb_graph_blow.sigma(); orb_graph_blow.rho()
            [[1, 8, 6],
             [20, 28, 36],
             [19, 27, 35],
             [2, 9, 4],
             [22, 30, 32],
             [21, 29, 31],
             [3, 7, 5],
             [24, 26, 34],
             [23, 25, 33],
             [10, 14, 18],
             [11, 15, 16],
             [12, 13, 17]]
            Ribbon graph of genus 0 and 4 boundary components
            (1,2,3)(4,5,6)(7,8,9)(10,11,12)
            (1,12)(2,6)(3,8)(4,11)(5,9)(7,10)
            sage: edges = [[3*x+1,3*x+2,3*x+3] for x in range(28)]
            sage: edges += [[84+2*x+1, 84+2*x+2] for x in range(6)]; print(edges)
            [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15], [16, 17, 18], [19, 20, 21], [22, 23, 24], [25, 26, 27], [28, 29, 30], [31, 32, 33], [34, 35, 36], [37, 38, 39], [40, 41, 42], [43, 44, 45], [46, 47, 48], [49, 50, 51], [52, 53, 54], [55, 56, 57], [58, 59, 60], [61, 62, 63], [64, 65, 66], [67, 68, 69], [70, 71, 72], [73, 74, 75], [76, 77, 78], [79, 80, 81], [82, 83, 84], [85, 86], [87, 88], [89, 90], [91, 92], [93, 94], [95, 96]]
            sage: s_ex = PermutationGroupElement([tuple(x) for x in edges]); s_ex
            (1,2,3)(4,5,6)(7,8,9)(10,11,12)(13,14,15)(16,17,18)(19,20,21)(22,23,24)(25,26,27)(28,29,30)(31,32,33)(34,35,36)(37,38,39)(40,41,42)(43,44,45)(46,47,48)(49,50,51)(52,53,54)(55,56,57)(58,59,60)(61,62,63)(64,65,66)(67,68,69)(70,71,72)(73,74,75)(76,77,78)(79,80,81)(82,83,84)(85,86)(87,88)(89,90)(91,92)(93,94)(95,96)
            sage: r_ex = PermutationGroupElement('(1,15)(2,27)(3,39)(4,30)(5,42)(6,18)(7,33)(8,45)(9,21)(10,36)(11,48)(12,24)(50,17)(51,86)(49,31)(85,54)(53,35)(52,13)(56,20)(55,34)(57,88)(60,87)(58,16)(59,38)(62,23)(61,37)(63,90)(66,89)(65,41)(64,19)(69,92)(72,91)(67,40)(68,26)(70,22)(71,44)(75,94)(78,93)(74,29)(73,43)(77,47)(76,25)(81,96)(84,95)(79,46)(80,32)(82,28)(83,14)')
            sage: R_ex = RibbonGraph(s_ex,r_ex)
            sage: R_ex.genus(); R_ex.number_boundaries()
            7
            2
            sage: metric = 96*[1/88]
            sage: T_ex = TatGraph(R_ex, metric)
            sage: T_ex.orbit_graph()
            ([[1, 6, 9, 12, 2, 4, 7, 10, 3, 5, 8, 11],
              [13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46],
              [14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47],
              [15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48],
              [49, 55, 61, 67, 73, 79, 52, 58, 64, 70, 76, 82],
              [50, 56, 62, 68, 74, 80, 53, 59, 65, 71, 77, 83],
              [51, 57, 63, 69, 75, 81, 54, 60, 66, 72, 78, 84],
              [85, 87, 89, 91, 93, 95, 86, 88, 90, 92, 94, 96]],
             Ribbon graph of genus 0 and 2 boundary components)
        """
        aux_sigma = [list(x) for 
                     x in self._sigma.cycle_tuples(singletons = 1)]
        darts = [x for y in aux_sigma for x in y]

        n = self.order()
        orbits = []
        while (len(darts) > 0):
            aux_edgeorbit = [darts[0]]
            darts.remove(aux_edgeorbit[-1])

            for i in range (n - 1):
                aux_edgeorbit.append(
                                     safewalk(self._ribbon, 
                                              self._metric, 
                                              aux_edgeorbit[-1],
                                              relative_boundary = 
                                              self._relative_boundary,
                                              sf_length = self._sf_length
                                             )
                                    )
                darts.remove(aux_edgeorbit[-1])

            orbits.append(aux_edgeorbit)

        aux_orbits_sigma = copy(orbits)
        aux_orbits_rho = copy(orbits)
        orbit_sigma = []
        orbit_rho = []

        while aux_orbits_sigma:
            orbit_sigma += [[]]
            pos = _find(orbits,
                       aux_orbits_sigma[0][0]
                       )
            orbit_sigma[-1].append(pos[0]+1)
            del aux_orbits_sigma[0]
            finish_cycle = False
            while not finish_cycle:
                aux_pos = _find(orbits,
                               self._sigma(orbits
                                                 [
                                                    orbit_sigma[-1][-1]-1
                                                 ][0]
                                            )
                           )
                if (pos[0] == aux_pos[0]):
                    finish_cycle = True
                elif (pos[0] != aux_pos[0]):
                    del_pos = _find(aux_orbits_sigma,
                                self._sigma(orbits
                                                    [
                                                        orbit_sigma[-1][-1]-1
                                                    ][0]
                                                )
                            )
                    orbit_sigma[-1].append(aux_pos[0]+1)
                    del aux_orbits_sigma[del_pos[0]]

        while aux_orbits_rho:
            orbit_rho += [[]]
            pos = _find(orbits,
                       aux_orbits_rho[0][0]
                       )

            orbit_rho[-1].append(pos[0]+1)
            del aux_orbits_rho[0]
            finish_cycle = False

            while not finish_cycle:
                aux_pos = _find(orbits,
                               self._rho(orbits
                                                 [
                                                  orbit_rho[-1][-1]-1
                                                 ][0]
                                            )
                           )
                if (pos[0] == aux_pos[0]):
                    finish_cycle = True
                elif (pos[0] != aux_pos[0]):
                    del_pos = _find(aux_orbits_rho,
                                self._rho(orbits
                                                [
                                                    orbit_rho[-1][-1]-1
                                                ][0]
                                                )
                            )
                    orbit_rho[-1].append(aux_pos[0]+1)
                    del aux_orbits_rho[del_pos[0]]

        return orbits, RibbonGraph(
                                   PermutationGroupElement([tuple(x) for x in orbit_sigma]), 
                                   PermutationGroupElement([tuple(x) for x in orbit_rho])
                                   )

def blow_up(tat_graph, vertex, epsilon):
    r"""
    Return the result of performing a blow up on the orbit of the ``vertex``
    of length ``epsilon``.

    More concretely, the operation consists in taking the orbit of the given
    vertex performing a real blow up at each of these vertices while seeing
    the graph embedded in the thickenened surface. The result is a relative
    tête-à-tête graph whose relative part is the union of the previous relative
    tête-à-tête part (if any) with the new boundary components created by the
    blowing up (as many as vertices were in the orbit of ``vertex``. The epsilon
    input is taken to produce the final metric in the followin way: from each
    dart adjacent to the vertices in the orbit we substract a length of epsilon.
    Each new edge introduced has a total length of 2 epsilon in order to preserve
    the tête-à-tête property. For this reason, epsilon has to be strictly smaller
    than the lengths of all the darts adjacent to a vertex of the orbit.

    INPUT:

    - ``tat_graph`` -- a tête-à-tête graph or a relative tête-à-tête graph.
    - ``vertex`` -- a positive integer from 0 to len(tat_graph.sigma())
      indicating one of the vertices where the blow up will be performed
      (observe that the blow up is performed on all the vertices in the
      orbit of ``vertex``). Also we ask that ``vertex`` is not one vertex
      that is already on the relative boundary of ``tat_graph``.
    - ``epsilon`` -- a rational number that has to be smaller than all
      lengths of the darts adjacent to ``vertex``.

    OUTPUT:

    - A relative tête-à-tête graph resulting from the blow-up of ``tat_graph``
      at the indicated vertex of the indicated length.


    EXAMPLES:

    We define a bipartite tete-a-tete graph of type `(3,4)` and first we 
    blow up the orbit of 3 vertices; then we blow up the orbit formed by
    4 vertices. As we can see we obtain two differente relative tete-a-tete
    graphs::

        sage: T = bipartite_tat_graph(3,4); T
        Tete-a-tete graph of order 12 on a ribbon graph of genus 3 and 1 boundary components.
        sage: B1=blow_up(T,2,1/10); B1
        Relative tête-à-tête graph of order 12 on a ribbon graph of genus 3 and 4 boundary components; where 3 boundary components are part of the relative boundary and might be permuted by the automorphism induced.
        sage: B2 = blow_up(T,4,1/30); B2
        Relative tête-à-tête graph of order 12 on a ribbon graph of genus 3 and 5 boundary components; where 4 boundary components are part of the relative boundary and might be permuted by the automorphism induced.


    In the previous example we blew up with two different lengths, both
    smaller than `1/2` which is the length of each dart of T. Now we try to
    blow up with a length greater than that and it raises an AssertionError::

        sage: B3 = blow_up(T,2,2/3)
        Traceback (most recent call last):
        ...
        AssertionError
    """

    #we intialize the vector that is going to save the orbit of the 
    #vertex.
    orb_vertex = []

    aux_sigma = [list(x) for 
                     x in tat_graph._sigma.cycle_tuples()]

    aux_rho = [list(x) for 
                   x in tat_graph._rho.cycle_tuples()]
    orb_vector, orb_graph = tat_graph.orbit_graph()
    metric = copy(tat_graph._metric)
    relative_boundary = copy(tat_graph._relative_boundary)
    aux_dart = aux_sigma[vertex][0]
    aux_pos = _find(orb_vector, aux_dart)
    for i in range(len(orb_vector[aux_pos[0]])):
        new_v = _find(aux_sigma,orb_vector[aux_pos[0]][i])[0]
        if new_v not in orb_vertex:
            orb_vertex.append(new_v)

    #check that the lengths of the darts adjacent to the vertices in the orbit
    #are bigger than epsilon (it is enought to check one vertex of the orbit)
    for i in aux_sigma[orb_vertex[0]]:
         assert tat_graph._metric[i] > epsilon

    #we hold in one variable the maximun of the darts
    darts = [x for y in aux_sigma for x in y]
    m_dart = max(darts) 

    #finally we create the relative tat graph by modifying each vertex
    #in the orbit. Actually we substitute each vertex for as many vertices
    #as its valency. Each vertex produces a new relative component.
    for j in range(len(orb_vertex)):
        darts = [x for y in aux_sigma for x in y]
        m = max(darts)
        aux_ver = aux_sigma[orb_vertex[j]]
        for k in range(len(aux_ver)):
            aux_sigma.append([m+2*k+1,aux_ver[k],m+2*k+2])
            aux_rho.append([m+2*k+2, m + ((2*k+3) % (2*len(aux_ver)))])
            metric[aux_ver[k]-1] -= epsilon
            metric += [epsilon,epsilon]
    aux_sigma = [aux_sigma[s] for s in range(len(aux_sigma)) if s not in orb_vertex]
    aux_ribbon = RibbonGraph(
                                   PermutationGroupElement([tuple(x) for x in aux_sigma]), 
                                   PermutationGroupElement([tuple(x) for x in aux_rho])
                                   )
    relative_boundary = [x for x in aux_ribbon.boundary() if min(x) > m_dart]
    return TatGraph(aux_ribbon, metric, relative_boundary,tat_graph._sf_length)
