r"""
The network coding instance class

Network Coding is a paradigm for communication over networks assumed to
be directed acyclic hypergraphs with error-free links. An instance of
the network coding problem is described by tuple $(G,S,T,\beta)$
where $G=(V,E)$ is a directed acyclic hypergraph, $S,T\subset V$ are
subsets of source and sink nodes resp. and $\beta: T\rightarrow 2^S$ is
the sink demand function. A network code for a network coding instance
of size `n` is a collection of `n` discrete random variables (in the
most general sense). Entropy function associated with a collection
$\{X_1,\hdots,X_n\}$ of random variables is a set function
$h:2^{[n]}\rightarrow \mathbb R$ mapping each subset of indices of
random variables to respective Shannon entropy.

An arbitrary set function is a valid entropy function if it satisfies
certain inequalities known as information inequalities (not completely
characterized for
$n\geq 4$). For an arbitrary set function to have
arisen from a network code for a network coding instance, it must satisfy
additional linear constraints imposed by the network.  A set function which
is the rank
function of a representable matroid is always a valid entropy function. A
scalar linear network code can be then defined as a representable matroid
with ground set $E$ and a
map $f:E\rightarrow [n]$ s.t. the set function $rank(f(\cdot))$ satisfies
the constraints imposed by the network coding instance.

Basic Usage
===========
Given an instance $(G,S,T,\beta)$ of network coding, one can construct a
``NCinstance`` object using following parameters:

1. ``nsrc`` is $\vert S\vert$ and we label source messages with the set
$\{1,\hdots,\vert S\vert\}$
2. ``constraints`` is a list of network constraints (barring source
independence, which is implicit with the labels)
3. ``size`` is $\vert E\vert$ which is the total number of messages to be
transmitted

The simplest example is the 2-unicast butterfly network described as follows:

1. $G=(\{a,b,c,d,e\},\{(a,\{c,d\}), (b,\{c,e\}), (c,\{d,e\})\}) \}$ with 5
nodes and 3 hyperedges.
2. Sources $S=\{a,b\}$ with associated message labels $\{1,2\}$ and ``nsrc``
is $2$
3. Sinks $T=\{d,e\}$
3. Messages labeled $1,2,3$ associated with hyperedges
$(a,\{c,d\}), (b,\{c,e\}), (c,\{d,e\})$ respectively
4. Sink demand function $\beta$  s.t. $\beta(d)=b$ amd $\beta(e)=a$
5. Network constraints specified as a list of lists
``constraints=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[2,3],[1,2,3]]``
where each member ``[X,Y]`` enforces equality of the entropy of sets of
message random variables indexed by ``X`` and ``Y``. Independence of random
variables associated with
source messages is implied.
6. ``size`` is $3$

Construction
============
One can construct a ``NCinstance`` object by using
:func:`Matroid() <sage.matroids.network_coding.NCinstance>` with ``size``,
``nsrc`` and  ``constraints`` as input.

    sage: size=3
    sage: nsrc=2
    sage: constraints=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[2,3],[1,2,3]]
    sage: from sage.matroids.network_coding import *
    sage: NCinstance(size,nsrc,constraints)
    A network coding instance of size 3 with 2 sources

Class methods
=============

- :class:`NCinstance`

    - :func:`is_field_scalar_solvable() <sage.matroids.network_coding.NCinstance.is_field_scalar_solvable>`
    - :func:`is_achievable_rate_vector() <sage.matroids.linear_matroid.NCinstance.is_achievable_rate_vector>`

AUTHORS:

- Jayant Apte (2015-07-08): initial version

EXAMPLES::

    sage: size=3
    sage: nsrc=2
    sage: constraints=[[[1,2],[1,2,3]],[[1,3],[1,2,3]],[2,3],[1,2,3]]
    sage: from sage.matroids.network_coding import *
    sage: NCinstance(size,nsrc,constraints)
    A network coding instance of size 3 with 2 sources

REFERENCES
==========

..  [Yeung] Raymond W. Yeung. 2008. Information Theory and Network Coding (1 ed.). Springer Publishing Company, Incorporated.
..  [Oxley] James Oxley, "Matroid Theory, Second Edition". Oxford University Press, 2011.

"""

#*****************************************************************************
#       Copyright (C) 2015 Jayant Apte <jayant91089@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.structure.sage_object import SageObject
from sage.matroids.network_coding_helpers import *


class NCinstance(SageObject):

    r"""
    The network coding instance class

    """

    def __init__(self, size, nsrc, constraints):
        assert size > 1
        assert nsrc >= 1
        self.size = int(size)  # No. of random variables
        self.nsrc = int(nsrc)  # No. of Sources
        self.constraints = constraints  # constraints

    def _repr_(self):
        r"""
        Return a string representation of the network coding instance.

        EXAMPLES::

            sage: from sage.matroids.networks_catalog import *
            sage: V = VamosNet()
            sage: sage.matroids.network_coding.NCinstance._repr_(V)
            'A network coding instance of size 8 with 4 sources'

        """
        if len(self.constraints) == 0:
            return "An empty network coding instance of size %s" % (self.size)
        else:
            return "A network coding instance of size %s with %s sources" % (
                self.size, self.nsrc)

    def is_field_scalar_solvable(self, fieldsize):
        r"""
        Check if the network coding instance is scalar solvable over a field.

        Scalar solvability means that there exists a scalar linear network code
        over the finite field of size ``fieldsize`` such that a entropy of
        sources and edge messages is 1 bit each.

        INPUT:

        - ``fieldsize`` -- size of finite field over which to test solvability.

        OUTPUT:

        A 2-tuple containing:

        1. A boolean indicating whether ``self`` is scalar solvable over finite
        field of size ``fieldsize``.
        2. If 1. is ``True``, a pair ``[M,pcode]`` where ``M`` is a matroid of
        size ``self.size`` representable over finite field of size
        ``fieldsize`` and ``pcode`` is a dictionary mapping groundset elements
        to messages associated with the network

        EXAMPLES:

        sage: from sage.matroids.networks_catalog import *
        sage: V=VamosNet();
        sage: V.is_field_scalar_solvable(2) # Vamos matroid is not linearly representable
        (False, [])
        sage: F=FanoNet();
        sage: F.is_field_scalar_solvable(2) # Fano matroid is representable over GF(2)
        (True,
         [Binary matroid of rank 3 on 7 elements, type (3, 0),
          {0: 1, 1: 2, 2: 3, 3: 4, 4: 6, 5: 7, 6: 5}])
        sage: r,code=F.is_field_scalar_solvable(2)
        sage: r
        True
        sage: code
        [Binary matroid of rank 3 on 7 elements, type (3, 0),
         {0: 1, 1: 2, 2: 3, 3: 4, 4: 6, 5: 7, 6: 5}]
        sage: code[0].representation()
        [1 0 0 1 1 1 0]
        [0 1 0 1 0 1 1]
        [0 0 1 0 1 1 1]
        sage: code[0].is_isomorphic(matroids.named_matroids.Fano())
        True
        sage: F.is_field_scalar_solvable(3) # Fano matroid is not representable over GF(3)
        (False, [])

        """
        assert fieldsize.is_prime_power()
        return ratecertgen(self.constraints, [1] * self.size, self.nsrc, self.size, fieldsize)

    def is_achievable_rate_vector(self, rate_vector, fieldsize):
        r"""
        Check if a 0-1 rate vector is achievable over a given field.

        A 0-1 ``rate_vector`` is achievable means that there exists a scalar
        linear network code with message entropies specified by
        ``rate vector`` over finite field of size ``fieldsize``.

        INPUT:

        - ``rate_vector`` -- A 0-1 vector whose achievability we want to test.
        - ``fieldsize`` -- size of finite field over which to test solvability.

        OUTPUT:
        A 2-tuple containing:

        1. A boolean indicating whether ``rate_vector`` is achievable with
        scalar linear network coding over finite field of size ``fieldsize``.
        2. If 1. is ``True``, a pair ``[M,pcode]`` where ``M`` is a matroid of
        size ``self.size`` representable over finite field of size.
        ``fieldsize`` and ``pcode`` is a dictionary mapping groundset elements
        to messages associated with the network.

        EXAMPLES:

            sage: from sage.matroids.networks_catalog import *
            sage: V=VamosNet();
            sage: V.is_achievable_rate_vector([1,1,1,1,1,1,1,1],2)
            (False, [])
            sage: r,code=V.is_achievable_rate_vector([0,1,0,1,1,1,1,1],2) # delete sources 1 and 3
            sage: r
            True
            sage: code
            [Binary matroid of rank 2 on 8 elements, type (0, 2),
             {0: 2, 1: 4, 2: 1, 3: 3, 4: 5, 5: 6, 6: 7, 7: 8}]
            sage: code[0].representation()
            [1 0 0 0 1 1 1 1]
            [0 1 0 0 1 1 1 1]
        """
        assert fieldsize.is_prime_power()
        assert len(rate_vector) == self.size
        assert len([i for i in rate_vector if i == 1]) + len(
            [i for i in rate_vector if i == 0]) == self.size
        return ratecertgen(self.constraints, rate_vector, self.nsrc, self.size, fieldsize)
