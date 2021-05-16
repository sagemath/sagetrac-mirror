# -*- coding: utf-8 -*-
r"""
Generic structures for linear codes over the rank metric

Rank Metric
===========

In coding theory, the most common metric is the Hamming metric, where distance
between two codewords is  given by the number of positions in which they differ.
An alternative to this is the rank metric. Take two fields, `F_q` and `F_{q^m}`,
and define a code `C` to be a set of vectors of length `n` with entries from
`F_{q^m}`. Let `c` be a codeword. We can represent it as an `m \times n` matrix
`M` over `F_q`.

A detailed description on the relationship between the two representations can
be found in :meth:`sage.coding.linear_rank_metric.to_matrix_representation`
and :meth:`sage.coding.linear_rank_metric.from_matrix_representation`.

We can define a metric using the rank of the matrix representation of the
codewords. A distance between two codewords `a, b` is the rank of the matrix
representation of `a - b`. A weight of a codeword `c` is the rank of the matrix
representation of `c`.

This module allows representing rank metric codes which are linear over the
big field `F_{q^m}`, i.e. the usual linearity condition when the codewords are
considered in vector form. One can also consider rank metric codes which are
only linear over `F_q`, but these are not currently supported in SageMath.

Note that linear rank metric codes per the definition of this file are
mathematically just linear block codes, and so could be considered as a
:class:`sage.coding.linear_code.LinearCode`. However, since most of the
functionality of that class is specific to the Hamming metric, the two notions
are implemented as entirely different in SageMath. If you wish to investigate
Hamming-metric properties of a linear rank metric code ``C``, you can easily
convert it by calling ``C_hamm = LinearCode(C)``.

Linear Rank Metric Code and Gabidulin Codes
===========================================

The class :class:`sage.coding.linear_rank_metric.LinearRankMetricCode` is the
analog of :class:`sage.coding.linear_code.LinearCode`, i.e. it is a generator
matrix-based representation of a linear rank metric code without specific
knowledge on the structure of the code.

Gabidulin codes are the main family of structured linear rank metric codes.
These codes are the rank-metric analog of Reed-Solomon codes.

``AbstractLinearRankMetricCode``
--------------------------------

This is a base class designed to contain methods, features and parameters
shared by every linear rank metric code. For instance, generic algorithms for
computing the minimum distance, etc. Many of these algorithms are slow,
e.g. exponential in the code length. It also contains methods for swapping
between vector and matrix representation of elements.

``AbstractLinearCodeNoMetric`` is an abstract class for linear rank metric codes,
so any linear rank metric code  class should inherit from this class.
Also ``AbstractLinearCodeNoMetric`` should never itself be instantiated.

See :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
for details and examples.

``LinearRankMetricCode``
------------------------

This class is used to represent arbitrary and unstructured linear rank metric
codes. It mostly relies directly on generic methods provided by
``AbstractLinearRankMetricCode``, which means that basic operations on the code
(e.g. computation of the minimum distance) will use slow algorithms.

A ``LinearRankMetricCode`` is instantiated by providing a generator::

    sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
    sage: C = codes.LinearRankMetricCode(G, GF(4))
    sage: C
    [3, 2] linear rank metric code over GF(64)/GF(4)
    sage: C.generator_matrix()
    [1 1 0]
    [0 0 1]
    sage: c = vector(GF(64), (1, 1, 1))
    sage: c in C
    True

Further references
------------------

Read more about
`rank metric and Gabidulin codes <https://en.wikipedia.org/wiki/Rank_error-correcting_code>`_

AUTHORS:

- Marketa Slukova (2019-08-16): initial version

- Maaike van Leuken (2021-03-05): update, added basic functionalities and functions for GRS algorithm

TESTS::

    sage: MS = MatrixSpace(GF(2),4,7)
    sage: G = MS([[1,1,1,0,0,0,0], [1,0,0,1,1,0,0], [0,1,0,1,0,1,0], [1,1,0,1,0,0,1]])
    sage: C = LinearCode(G)
    sage: C == loads(dumps(C))
    True
"""

# ****************************************************************************
#       Copyright (C) 2019 MARKETA SLUKOVA <em.slukova@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from sage.categories.fields import Fields
from sage.matrix.constructor import Matrix
from sage.structure.element import is_Matrix, is_Vector
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity
from random import randint
from sage.modules.free_module import VectorSpace
from sage.matrix.constructor import matrix
from sage.modules.free_module import span

from sage.coding.linear_code_no_metric import AbstractLinearCodeNoMetric
from sage.coding.linear_code import LinearCodeGeneratorMatrixEncoder
from sage.coding.decoder import Decoder
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.combinat.q_analogues import q_binomial


def to_matrix_representation(v, sub_field=None, basis=None):
    r"""
    Return a matrix representation of ``v`` over ``sub_field`` in terms of
    ``basis``.

    Let `(b_1, b_2, \ldots, b_m)`, `b_i \in GF(q^m)`, be a basis of `GF(q^m)` as
    a vector space over `GF(q)`. Take an element `x \in GF(q^m)`. We can write
    `x` as `x = u_1 b_1 + u_2 b_2 + \ldots u_m b_m`, where `u_i \in GF(q)`. This
    way we can represent an element from `GF(q^m)` as a vector of length `m`
    over `GF(q)`.

    Given a vector ``v`` of length `n` over some field `F_{q^m}`, we can
    represent each entry as a vector of length `m`, yielding an `m \times n`
    matrix over ``sub_field``. In case ``sub_field`` is not given, we take the
    prime subfield `F_p` of `F_{q^m}`.

    INPUT:

    - ``v`` -- a vector over some field `F_{q^m}`

    - ``sub_field`` -- (default: ``None``) a sub field of `F_{q^m}`. If not
      specified, it is the prime subfield `F_p` of `F_{q^m}`.

    - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
      ``sub_field``. If not specified, given that `q = p^s`, let
      `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
      represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import to_matrix_representation
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: to_matrix_representation(a, GF(4))
        [1 1 1]
        [1 1 0]
        [0 0 0]

        sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
        sage: to_matrix_representation(m)
        Traceback (most recent call last):
        ...
        TypeError: Input must be a vector
    """
    if not is_Vector(v):
        raise TypeError("Input must be a vector")
    base_field = v.base_ring()
    if not sub_field:
        sub_field = base_field.prime_subfield()
    n = v.length()
    m = base_field.degree()//sub_field.degree()
    extension, to_big_field, from_big_field = base_field.vector_space(sub_field, basis, map=True)
    return Matrix(sub_field, m, n, lambda i, j: from_big_field(v[j])[i])

def from_matrix_representation(w, base_field=None, basis=None):
    r"""
    Return a vector representation of a matrix ``w`` over ``base_field`` in terms
    of ``basis``.

    Given an `m \times n` matrix over `F_q` and some ``basis`` of `F_{q^m}`
    over `F_q`, we can represent each of its columns as an element of `F_{q^m}`,
    yielding a vector of length `n` over `F_q`.

    In case ``base_field`` is not given, we take `F_{q^m}`, the field extension of
    `F_q` of degree `m`, the number of rows of ``w``.

    INPUT:

    - ``w`` -- a matrix over some field `F_q`

    - ``base_field`` -- (default: ``None``) an extension field of `F_q`. If not
      specified, it is the field `F_{q^m}`, where `m` is the number of rows of
      ``w``.

    - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
      ``F_q``. If not specified, given that `q = p^s`, let
      `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
      represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import from_matrix_representation
        sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
        sage: from_matrix_representation(m)
        (z6 + 1, z6 + 1, 1)

        sage: v = vector(GF(4), (1, 0, 0))
        sage: from_matrix_representation(v)
        Traceback (most recent call last):
        ...
        TypeError: Input must be a matrix
    """
    if not is_Matrix(w):
        raise TypeError("Input must be a matrix")
    sub_field = w.base_ring()
    if not base_field:
        base_field = sub_field.extension(w.nrows())
    v = []
    extension, to_big_field, from_big_field = base_field.vector_space(sub_field, basis, map=True)
    for i in range(w.ncols()):
        v.append(to_big_field(w.column(i)))
    return vector(v)

def rank_weight(c, sub_field=None, basis=None):
    r"""
    Return the rank of ``c`` as a matrix over ``sub_field``.

    If ``c`` is a vector over some field `F_{q^m}`, the function converts it
    into a matrix over `F_q`.

    INPUT:

    - ``c`` -- a vector over some field `F_{q^m}`; or a matrix over `F_q`

    - ``sub_field`` -- (default: ``None``) a sub field of `F_{q^m}`. If not
      specified, it is the prime subfield `F_p` of `F_{q^m}`.

    - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
      ``sub_field``. If not specified, given that `q = p^s`, let
      `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
      represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_weight
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: rank_weight(a, GF(4))
        2
    """
    if is_Vector(c):
        c = to_matrix_representation(c, sub_field, basis)
    return c.rank()

def rank_distance(a, b, sub_field=None, basis=None):
    r"""
    Return the rank of ``a`` - ``b`` as a matrix over ``sub_field``.

    Take two vectors ``a``, ``b`` over some field `F_{q^m}`. This function
    converts them to matrices over `F_q` and calculates the rank of their
    difference.

    If ``sub_field`` is not specified, we take the prime subfield `F_q` of
    `F_{q^m}`.

    INPUT:

    - ``a`` -- a vector over some field `F_{q^m}`

    - ``b`` -- a vector over some field `F_{q^m}`

    - ``sub_field`` -- (default: ``None``) a sub field of `F_{q^m}`. If not
      specified, it is the prime subfield `F_p` of `F_{q^m}`.

    - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
      ``sub_field``. If not specified, given that `q = p^s`, let
      `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
      represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

    EXAMPLES::

        sage: from sage.coding.linear_rank_metric import rank_distance
        sage: x = GF(64).gen()
        sage: a = vector(GF(64), (x + 1, x + 1, 1))
        sage: b = vector(GF(64), (1, 0, 0))
        sage: rank_distance(a, b, GF(4))
        2

        sage: c = vector(GF(4), (1, 0, 0))
        sage: rank_distance(a, c, GF(4))
        Traceback (most recent call last):
        ...
        ValueError: The base field of (z6 + 1, z6 + 1, 1) and (1, 0, 0) has to be the same

        sage: d = Matrix(GF(64), (1, 0, 0))
        sage: rank_distance(a, d, GF(64))
        Traceback (most recent call last):
        ...
        TypeError: Both inputs have to be vectors

        sage: e = vector(GF(64), (1, 0))
        sage: rank_distance(a, e, GF(64))
        Traceback (most recent call last):
        ...
        ValueError: The length of (z6 + 1, z6 + 1, 1) and (1, 0) has to be the same
    """
    if not (a.base_ring() == b.base_ring()):
        raise ValueError("The base field of {} and {} has to be the same".format(a, b))
    if not (is_Vector(a) and is_Vector(b)):
        raise TypeError("Both inputs have to be vectors")
    if not len(a) == len(b):
        raise ValueError("The length of {} and {} has to be the same".format(a, b))

    a = to_matrix_representation(a, sub_field, basis)
    b = to_matrix_representation(b, sub_field, basis)
    return (a - b).rank()


class AbstractLinearRankMetricCode(AbstractLinearCodeNoMetric):
    r"""
    Abstract class for linear rank metric codes.

    This class contains methods that can be used on families of linear rank
    metric codes. Every linear rank metric code class should inherit from this
    abstract class.

    This class is intended for codes which are linear over the ``base_field``.

    Codewords of rank metric codes have two representations. They can either be
    written as a vector of length `n` over `GF(q^m)`, or an `m \times n` matrix
    over `GF(q)`. This implementation principally uses the vector representation.
    However, one can always get the matrix representation using the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.to_matrix`
    method. To go back to a vector, use the
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.from_matrix`
    method.

    Instructions on how to make a new family of rank metric codes is analogous
    to making a new family of linear codes over the Hamming metric, instructions
    for which are in :class:`sage.coding.linear_code.AbstractLinearCode`. For an
    example on, see
    :meth:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__`

    .. WARNING::

        A lot of methods of the abstract class rely on the knowledge of a generator matrix.
        It is thus strongly recommended to set an encoder with a generator matrix implemented
        as a default encoder.
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, sub_field, length, default_encoder_name,
            default_decoder_name, basis=None):
        r"""
        Initialize mandatory parameters that every linear rank metric code has.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every linear rank metric code.
        The class :class:`sage.coding.linear_rank_metric.AbstractLinearRankMetricCode`
        should never be directly instantiated.

        INPUT:

        - ``base_field`` -- the base field of ``self``

        - ``sub_field`` -- the sub field of ``self``

        - ``length`` -- the length of ``self`` (a Python int or a Sage Integer),
          must be > 0 and at most the degree of the field extension

        - ``default_encoder_name`` -- the name of the default encoder of ``self``

        - ``default_decoder_name`` -- the name of the default decoder of ``self``

        - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
          ``sub_field``. If not specified, given that `q = p^s`, let
          `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
          represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

        EXAMPLES:

        The following example demonstrates how to use subclass
        `AbstractLinearRankMetricCode` for representing a new family of rank
        metric codes. The example is a rank repetition code::

             sage: from sage.coding.linear_rank_metric import AbstractLinearRankMetricCode
             sage: class RankRepetitionCode(AbstractLinearRankMetricCode):
             ....:   def __init__(self, base_field, sub_field, length):
             ....:       sage.coding.linear_rank_metric.AbstractLinearRankMetricCode.__init__(self, base_field, sub_field, length, "GeneratorMatrix", "NearestNeighbor")
             ....:       beta = base_field.gen()
             ....:       self._generator_matrix = matrix(base_field, [[ beta^i for i in range(length) ]])
             ....:   def generator_matrix(self):
             ....:       return self._generator_matrix
             ....:   def _repr_(self):
             ....:       return "[%d, %d] rank-metric repetition code over GF(%s)" % (self.length(), self.dimension(), self.base_field().cardinality())

        We now instantiate a member of our newly made code family::

            sage: C = RankRepetitionCode(GF(8), GF(2), 3)

        We can check its existence and parameters::

            sage: C
            [3, 1] rank-metric repetition code over GF(8)

        We can encode a vector::

            sage: word = vector(C.base_field(), [1])
            sage: E = codes.encoders.LinearCodeSystematicEncoder(C)
            sage: codeword = E(word)
            sage: codeword
            (1, z3, z3^2)

        We can get the matrix representation of the codeword::

            sage: C.matrix_form_of_vector(codeword)
            [1 0 0]
            [0 1 0]
            [0 0 1]

        We can decode the vector representation of the codeword::

            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decode_to_code(codeword)
            (1, z3, z3^2)
            sage: D.decode_to_message(codeword)
            (1)

        We can check that it is truly a part of the framework category::

            sage: C.parent()
            <class '__main__.RankRepetitionCode_with_category'>
            sage: C.category()
            Category of facade finite dimensional vector spaces with basis over Finite Field in z3 of size 2^3

        And any method that works on rank metric linear codes works for our new dummy code::

            sage: C.minimum_distance()
            3
            sage: C.metric()
            'rank'

        TESTS:

        If ``sub_field`` is not a field, an error is raised::

            sage: C = RankRepetitionCode(GF(8), ZZ, 3)
            Traceback (most recent call last):
            ...
            ValueError: 'sub_field' must be a field (and Integer Ring is not one)

        If ``sub_field`` is not a subfield of ``base_field``, an error is raised::

            sage: C = RankRepetitionCode(GF(8), GF(3), 2)
            Traceback (most recent call last):
            ...
            ValueError: 'sub_field' has to be a subfield of 'base_field'
        """
        self._registered_decoders["NearestNeighbor"] = LinearRankMetricCodeNearestNeighborDecoder
        if not sub_field.is_field():
            raise ValueError("'sub_field' must be a field (and {} is not one)".format(sub_field))
        if not (sub_field.degree().divides(base_field.degree()) and (sub_field.prime_subfield() == base_field.prime_subfield())):
            raise ValueError("'sub_field' has to be a subfield of 'base_field'")
        m = base_field.degree() // sub_field.degree()
        self._extension_degree = m
        self._sub_field = sub_field
        self._ranks = []
        self._indexed_ranks = []
        self._rank_distribution = []
        self._generic_constructor = LinearRankMetricCode
        
        super(AbstractLinearRankMetricCode, self).__init__(base_field, length, default_encoder_name, default_decoder_name, "rank")

    def sub_field(self):
        r"""
        Return the sub field of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.sub_field()
            Finite Field in z2 of size 2^2
        """
        return self._sub_field

    def extension_degree(self):
        r"""
        Return `m`, the degree of the field extension of ``self``.

        Let ``base_field`` be `GF(q^m)` and ``sub_field`` be `GF(q)`. Then this
        function returns `m`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.extension_degree()
            3
        """

        return self._extension_degree

    def field_extension(self):
        r"""
        Return the field extension of ``self``.

        Let ``base_field`` be some field `F_{q^m}` and ``sub_field`` `F_{q}`.
        This function returns the vector space of dimension `m` over `F_{q}`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.field_extension()
            Vector space of dimension 3 over Finite Field in z2 of size 2^2
        """
        return self.base_field().vector_space(self.sub_field(), map=False)

    def rank_distance_between_vectors(self, left, right):
        r"""
        Return the rank of the matrix of ``left`` - ``right``.

        INPUT:

        - ``left`` -- a vector over the ``base_field`` of ``self``

        - ``right`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: b = vector(GF(64), (1, 0, 0))
            sage: C.rank_distance_between_vectors(a, b)
            2
        """
        return rank_distance(left, right, self.sub_field())

    def minimum_distance(self):
        r"""
        Return the minimum distance of ``self``.

        This algorithm simply iterates over all the elements of the code and
        returns the minimum weight.

        EXAMPLES::

            sage: F.<a> = GF(8)
            sage: G = Matrix(F, [[1,a,a^2,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: C.minimum_distance()
            3
        """
        d = Infinity
        for c in self:
            if c == self.zero():
                continue
            d = min(self.rank_weight_of_vector(c), d)
        return d

    def rank_weight_of_vector(self, word):
        r"""
        Return the weight of the word, i.e. its rank.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.rank_weight_of_vector(a)
            2
        """
        return rank_weight(word, self.sub_field())

    def matrix_form_of_vector(self, word):
        r"""
        Return the matrix representation of a word.

        INPUT:

        - ``word`` -- a vector over the ``base_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: a = vector(GF(64), (x + 1, x + 1, 1))
            sage: C.matrix_form_of_vector(a)
            [1 1 1]
            [1 1 0]
            [0 0 0]
        """
        return to_matrix_representation(word, self.sub_field())

    def vector_form_of_matrix(self, word):
        r"""
        Return the vector representation of a word.

        INPUT:

        - ``word`` -- a matrix over the ``sub_field`` of ``self``

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: x = GF(64).gen()
            sage: m = Matrix(GF(4), [[1, 1, 1], [1, 1, 0], [0, 0, 0]])
            sage: C.vector_form_of_matrix(m)
            (z6 + 1, z6 + 1, 1)
        """
        return from_matrix_representation(word, self.base_field())
        
    def support(self, word, row_support = False):
        r"""
        Return the support of a codeword in ``self``.

        The support of a word in ``self.base_field()`` is the subspace over the not-extended 
        ``self.sub_field()`` generated by the word. This is the same as the span of that 
        word in ``self.sub_field()``. The dimension of this subspace generated by the word
        should be equal to the rank weight of that word. So for word ``word`` and LinearRankMetricCode
        ``C``, ``C.support(word).dimension() == C.rank_weight_of_vector(C.matrix_form_of_vector(word))``.

        INPUT:

        - ``word`` -- The word over which the support should be computed.

        - ``row_support`` -- Default the column support is computed, but one can also compute the row_support.

        EXAMPLES::
        
            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: word = C.random_element(); word # random                                                               
            (z6^5 + z6^4 + z6, z6^5 + z6^4 + z6, z6^4 + z6^2 + 1)
            sage: C.support(word) # random                                               
            Vector space of degree 3 and dimension 2 over Finite Field in z2 of size 2^2
            Basis matrix:
            [ 1  0 z2]
            [ 0  1  0]
        """
        M = self.matrix_form_of_vector(word)
        if row_support == False:
            basis = [M.column(i) for i in range(M.dimensions()[1])]
        else: 
            basis = [row for row in M]
        return span(basis, self.sub_field())

   
    def relative_distance(self):
        r"""
        Return the ratio of the minimum distance to the code length.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.relative_distance()
            0.3333333333333333
        """
        return self.minimum_distance() / self.length()

    
    def ranks(self, indexed=False):
        r"""      
        Return the rank of each element in ``self``.
        
        INPUT:

        - ``indexed`` -- boolean whether or not the resulting list should be indexed

        EXAMPLES::
        
            sage: C = codes.GabidulinCode(GF(8), 3, 1)                                               
            sage: C.ranks()                                                                 
            [0, 3, 3, 3, 3, 3, 3, 3]
            sage: C.ranks(indexed=True)                                                              
            [(0, 0), (3, 1), (3, 2), (3, 3), (3, 4), (3, 5), (3, 6), (3, 7)]

        .. WARNING::

            This computes the rank of each element in ``self``. This take a very
            long time, even for medium parameters.
        """
        if self._ranks == []:
            self._ranks = [self.rank_weight_of_vector(word) for word in self]
        if indexed == False: 
            return self._ranks
        elif self._indexed_ranks == []: 
            self._indexed_ranks = [(self.ranks()[i], i) for i in range(len(self.ranks()))]
            return self._indexed_ranks
        else:
            return self._indexed_ranks
        
    def rank_distribution(self):
        r"""
        Return the rank distribution of the ``self``, i.e. the number
        of words in ``self`` that have rank equal to the index of the list.

        Since `rank(a*c) = rank(c)` for all codewords `c` and for non-zero
        `a` from `F_{q^m}`, we only compute the rank over codewords that have
        `1` as first non-zero element. For each of those codewords `c'`, there are `q^m - 1` 
        elements that have equal rank as `c'`. 

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.rank_distribution()                                                     
            [1, 315, 3780]

            sage: C = codes.GabidulinCode(GF(2**6), 5, 2)
            sage: C.rank_distribution()
            [1, 0, 0, 0, 1953, 2142]
        """
        if self._rank_distribution == []:
            x = []                                                                    
            for c in range(len(self.sub_field())**(self.extension_degree() * self.dimension())):
                for i in range(self.length()): 
                    if self[c][i] == 1 and False not in [self[c][j] == 0 for j in range(i)]:
                        x += [self[c]]
            y = [self.rank_weight_of_vector(elem) for elem in x]
            temp = [1] + [(len(self.base_field()) - 1) * y.count(i) for i in range(1, self.extension_degree()+1)]
            last = [index for index, item in enumerate(temp) if item != 0][-1]
            self._rank_distribution = temp[:last + 1]
        return self._rank_distribution
     
    def random_element_of_min_rank(self):
        r"""
        Find a random codeword of minimal rank.
        
        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: x = C.random_element_of_min_rank(); x # random
            (z6^3, z6^3, z6^5 + z6 + 1)
            sage: C.rank_weight_of_vector(x) == C.minimal_rank_in_code()                    
            True             
        """
        return self.random_element_of_rank(self.minimal_rank_in_code())
            
    def random_element_of_rank(self, r):
        r""" 
        Find a random codeword of rank ``r``.
        
        INPUT:

        - ``r`` -- the rank the random codeword should have

        TESTS:

        The given rank should not exceed the maximal rank in the code:

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4));
            sage: C.random_element_of_rank(3)
            Traceback (most recent call last):
            ...
            ValueError: The given rank should be smaller than 2.

        The given rank should be nonnegative:

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4));
            sage: C.random_element_of_rank(-1)
            Traceback (most recent call last):
            ...
            ValueError: The given rank should be greater than 0.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.random_element_of_rank(0)                                              
            (0, 0, 0)
            sage: C.random_element_of_rank(1) # random                                             
            (z6^5 + z6^3 + z6^2 + z6 + 1, z6^5 + z6^3 + z6^2 + z6 + 1, z6^5 + z6^4 + z6^2)
            sage: C.random_element_of_rank(2) # random                                             
            (z6^3 + z6^2 + z6, z6^3 + z6^2 + z6, z6^5 + z6^4 + z6^3 + z6 + 1)
        """
        if r > self.maximal_rank_in_code():
            raise ValueError("The given rank should be smaller than " + str(self.maximal_rank_in_code()) + ".")
        if r < 0:
            raise ValueError("The given rank should be greater than 0.")
        if r == 0:
            return self[0] # Only the first codeword (0,0,0) has rank 0.
        if self.rank_distribution()[r] == 0:
            raise ValueError("There are 0 codewords with rank " + str(r) + ".")
        
        temp = [index for val, index in self.ranks(indexed=True) if val==r]
        return self[temp[randint(0, len(temp))]]
    
    def random_element_of_max_rank(self):
        r"""
        Find a random codeword of maximal rank.
        
        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: x = C.random_element_of_max_rank(); x # random
            (z6^4 + z6^3 + z6^2 + z6, z6^4 + z6^3 + z6^2 + z6, z6^4 + z6^3 + z6)
            sage: C.rank_weight_of_vector(x) == C.maximal_rank_in_code()                    
            True         
        """
        return self.random_element_of_rank(self.maximal_rank_in_code())
    
    def minimal_rank_in_code(self):
        r"""
        Return the minimal rank that occurs in ``self``.

        This is computed by finding the index of the first non-zero element
        in ``self.rank_distribution()``, where rank 0 is excluded.
        
        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.minimal_rank_in_code()                                                  
            1
        """
        return min([index for index, item in enumerate(self.rank_distribution()[1:]) if item != 0]) + 1
        
    def maximal_rank_in_code(self):
        r"""     
        Return the highest rank that occurs in the code.
       
        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.maximal_rank_in_code()                                                  
            2
        """
        return len(self.rank_distribution()) - 1
    
    def characteristic(self):
        r"""
        Return the characteristic of the base ring of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0],[0,0,1]])                                     
            sage: C = codes.LinearRankMetricCode(G, GF(4)); C
            [3, 2] linear rank metric code over GF(64)/GF(4)
            sage: C.characteristic()                                                        
            2
        """
        return (self.base_ring()).characteristic()
        
    def poly_form_of_vector(self, vector, base_field=None):
        r"""
        Return the polynomial form of ``vector``. This is the inverse
        of what the ``vector()`` function does for elements in `\mathbb{F}_{q^m}`.
        
        INPUT:

        - ``vector`` --the vector the polynomial form should be returned
        
        EXAMPLES::

            sage: q = 2; n = 4; k = 1; m = 5; w = 2; d = 2                                           
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); C                                           
            [4, 1, 2] LRPC code over GF(32)/GF(2)
            sage: x = list(C.base_field())[16]; x  
            z5^4 + z5^3 + z5 + 1
            sage: vector(x)                                                                          
            (1, 1, 0, 1, 1)
            sage: C.poly_form_of_vector(vector(x))                                       
            z5^4 + z5^3 + z5 + 1
        """
        if base_field == None:
            base_field = self.base_field()
        return sum([vector[i]*base_field.gen()**(i) for i in range(len(vector))])
 
    def rewrite_in_basis(self, x, basis):
        r"""
        Rewrite an element ``x`` `\in \mathbb{F}_q^m` in terms of ``basis``. 
        
        INPUT:

        - ``x`` -- the element which should be rewritten

        - ``basis`` -- the basis in which the element should be rewritten

        EXAMPLES::

            sage: C = codes.LRPCCode(GF(16), 4, 2, 2); C                                             
            [4, 2, 2] LRPC code over GF(16)/GF(2)
            sage: G = C.generator_matrix()                                                           
            sage: H = C.parity_check_matrix(); H # random                                                    
            [z4^3 + z4^2 + z4        z4^2 + z4 z4^3 + z4^2 + z4                0]
            [z4^3 + z4^2 + z4             z4^3 z4^3 + z4^2 + z4             z4^3]
            sage: F = C.support_of_parity(); F # random                                                 
            Vector space of degree 4 and dimension 2 over Finite Field of size 2
            Basis matrix:
            [0 1 1 0]
            [0 0 0 1]
            sage: hxij = [[C.rewrite_in_basis(elem, F.basis()) for elem in row] for row in H]; hxij # random
            [[(1, 1), (1, 0), (1, 1), (0, 0)], [(1, 1), (0, 1), (1, 1), (0, 1)]]
        """ 
        f = Matrix(basis).transpose()
        x = vector(x)
        return f.solve_right(x)

    def recover_from_basis(self, m, basis):
        r"""
        Recover the element that was rewritten in ``basis``. 
        
        INPUT:

        - ``m`` -- the matrix which needs to be recovered

        - ``basis`` -- the basis in which the element was rewritten

        EXAMPLES::

            sage: C = codes.LRPCCode(GF(16), 4, 2, 2); C                                             
            [4, 2, 2] LRPC code over GF(16)/GF(2)
            sage: G = C.generator_matrix()                                                           
            sage: H = C.parity_check_matrix(); H # random                                                    
            [z4^3 + z4^2 + z4        z4^2 + z4 z4^3 + z4^2 + z4                0]
            [z4^3 + z4^2 + z4             z4^3 z4^3 + z4^2 + z4             z4^3]
            sage: F = C.support_of_parity(); F # random                                                 
            Vector space of degree 4 and dimension 2 over Finite Field of size 2
            Basis matrix:
            [0 1 1 0]
            [0 0 0 1]
            sage: hxij = [[C.rewrite_in_basis(elem, F.basis()) for elem in row] for row in H]; hxij # random
            [[(1, 1), (1, 0), (1, 1), (0, 0)], [(1, 1), (0, 1), (1, 1), (0, 1)]]
            sage: C.recover_from_basis(hxij, F.basis()) == H
            True
        """
        basis = matrix(basis)
        return matrix([[self.poly_form_of_vector(elem*basis) for elem in row] for row in m])
        
####################### attacks ###############################
    
    def get_var_matrix(self, rows, columns):
        r"""
        Create a matrix with unknown variables which reside in the ring.
        This could be used to set up a system of equations.
        
        INPUT:

        - ``rows`` -- the amount of rows in the variable matrix

        - ``columns`` -- the amount of columns in the variable matrix
        
        EXAMPLES::

            sage: q = 2; n = 5; k = 1; m = 5; w = 2; d = 2                                  
            sage: C = codes.LRPCCode(GF(q**m), n, k, d);
            sage: X = C.get_var_matrix(3, m); X
            [ x0  x1  x2  x3  x4]
            [ x5  x6  x7  x8  x9]
            [x10 x11 x12 x13 x14]
        """
        K = self.sub_field()
        K_ext, _ = K.extension(self.extension_degree(), 'z', map=True)
        R = PolynomialRing(K_ext, 'x', rows*columns, order='lex')
        x = R.gens()
        return Matrix(R, rows, columns, x) 

    def find_e_of_rank(self, r):
        r"""
        Finding an error vector `e \in \mathbb{F}_{q^m}^n` of a certain rank ``r``.

        Find a random error support `E` of dimension ``r``. Then the error vector `e` can
        be constructed by picking a random element from `E` for every of the ``self.length()`` 
        elements. 

        INPUT:

        - ``r`` -- the rank the error vector should have

        EXAMPLES::

            sage: q = 2; n = 22; k = 11; m = 11; r = 5; d = 2;                              
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); C                                  
            [22, 11, 2] LRPC code over GF(2048)/GF(2)
            sage: e = C.find_e_of_rank(r);
            sage: C.rank_weight_of_vector(e) == r                                         
            True
        """
        E = self.find_subspace_of_dim(r)
        temp = []
        for i in range(self.length()):
            temp += [E.random_element()]
        e = self.vector_form_of_matrix(matrix(temp).transpose())
        return e
        
    def canonical_projection_base_element(self, elem, i):
        r"""
        Return the ``i``'th canonical projection of an element ``elem`` in ``self.base_field()``.
        This is defined as a projection from `\mathbb{F}_{q^m}` on `\mathbb{F}_{q}`
        such that ``elem`` can be rewritten in terms of basis `\beta_0, ..., \beta_{m-1}`
        as follows: for ``elem`` `= \sum_{i = 0}^{m - 1} x_i \beta_i, \varphi_i(elem) = x_i`.
        See [AGHT2018]_.

        The basis for `\mathbb{F}_{q^m}` over `\mathbb{F}_{q}` is often the identity matrix,
        but it can also be some other basis. The functions ``self.matrix_form_of_vector()`` and
        ``vector()`` take care of this.

        INPUT:

        - ``elem`` -- element in ``self.base_field()``

        - ``i`` -- the index of the projection

        EXAMPLE::

            sage: q = 2; n = 4; k = 1; m = 5; w = 2; d = 2                                           
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); C                                           
            [4, 1, 2] LRPC code over GF(32)/GF(2)
            sage: x = list(C.base_field())[16]; x                                                                                  
            z5^4 + z5^3 + z5 + 1
            sage: vector(x)                                                                          
            (1, 1, 0, 1, 1)
            sage: C.canonical_projection_base_element(x, 0)                                          
            1
            sage: C.canonical_projection_base_element(x, 1)                                          
            1
            sage: C.canonical_projection_base_element(x, 2)                                          
            0
            sage: C.canonical_projection_base_element(x, 3)                                          
            1
            sage: C.canonical_projection_base_element(x, 4)                                          
            1
        """
        from sage.modules.free_module_element import FreeModuleElement_generic_dense 
        from sage.rings.finite_rings.element_givaro import FiniteField_givaroElement

        if type(elem) == FreeModuleElement_generic_dense:
            return vector(self.matrix_form_of_vector(elem))[i]
        elif type(elem) == FiniteField_givaroElement:
            return vector(elem)[i]
        
    def find_subspace_of_dim(self, r, space = None, already_tried=None):
        r"""
        Find a subspace `F` of `F_{q^m}`, or ``space`` if given, of dimension `r`.
        A list can be kept such that this function returns a subspace `F` which is not
        present in the list ``already_tried``.

        INPUT:

        - ``w`` -- the dimension of the new subspace 

        - ``space`` -- (default: ``None``) the space of which we wish to find a subspace.
          If not given it is ``self.base_field()``.

        - ``already_tried`` -- (default: ``None``) a list of subspaces. The found subspace should
          be a new one, so it should not be present in this list. If not given it is the empty list.

        EXAMPLES::

            sage: q = 2; n = 22; k = 11; m = 11; r = 5; d = 2;
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix();
            sage: G = C.generator_matrix(); 
            sage: F = C.find_subspace_of_dim(r); F # random                                                                           
            Vector space of degree 11 and dimension 5 over Finite Field of size 2
            Basis matrix:
            [1 1 0 0 0 0 0 0 1 1 0]
            [0 0 1 0 0 0 0 1 1 1 0]
            [0 0 0 1 0 0 1 0 1 0 1]
            [0 0 0 0 1 0 0 1 1 1 0]
            [0 0 0 0 0 1 1 0 0 0 1]        
        """
        q = self.characteristic()
        if space == None:
            space = self.base_field()
            z = self.extension_degree()
        else:
            z = space.dimension()
        if already_tried == None:
            already_tried = []
        number_of_Fs = q_binomial(z,r,q)
        while True:
            if len(already_tried) >= number_of_Fs:
                return "No more different subspaces."
            F = span([vector(space.random_element())], self.sub_field())
            while F.dimension() < r:
                F += span([vector(space.random_element())], self.sub_field())
            if F not in already_tried:
                break;
        return F

    def error_support_attack_case_1(self, w , r, y):
        r"""
        Case 1 of the GRS algorithm, `m \leq n` in [AGHT2018]_. We search for a subspace F of `\mathbb{F}_{q^m}` of
        dimension ``r`` which contains the (column) support of the error.
        If this is the case, the solution found to the system of equations is of rank ``r`` and then 
        also a solution to the RSD problem.
        If this is not the case, we try another subspace `F` which we have not tried before.
        
        The complexity is `\mathcal{O}((n-k)^3 m^3 q^{w\left \lceil{\frac{km}{n}}\right \rceil})`.

        INPUT:

        - ``w`` -- the rank of the error vector

        - ``r`` -- the dimension of the subspace F

        - ``y`` -- the received vector

        EXAMPLES::

            sage: q = 2; n = 5; k = 1; m = 5; w = 1; d = 2;                                   
            sage: r = m - ceil(k*m/n); r                                                                                  
            4
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix();   
            sage: x = C.random_element(); e = C.find_e_of_rank(w); 
            sage: y = x + e;                       
            sage: e_found = C.error_support_attack_case_1(w, r, y)                                  
            Probability of choosing F correctly the first time is: 15/31
            The complexity of this attack is: O(2^14)
            sage: e_found == e                                                                       
            True
        """
        from sage.functions.other import ceil
        from sage.functions.log import log

        q = self.sub_field().characteristic()
        n = self.length()
        k = self.dimension() 
        m = self.extension_degree()
        H = self.parity_check_matrix()

        complexity = (n-k)**3*m**3*q**(w*ceil(k*m/n))
        print("Probability of choosing F correctly the first time is:", q_binomial(r, w, q)/q_binomial(m, w, q))
        print("The complexity of this attack is: O(2^" + str(ceil(log(complexity, 2).n())) + ")")

        s = H*y
        number_of_Fs = q_binomial(m,r,q)
        tried = []
        while len(tried) < number_of_Fs:
            
            # Step 1) Pick random subspace F of `\mathbb{F}_{q^m}` of dimension `r >= w` with basis `F_{0}, ..., F_{r-1}`.
            F = self.find_subspace_of_dim(r, already_tried=tried)
            tried += [F]

            # Step 2) Set up the system.
            right = vector([self.canonical_projection_base_element(elem, i) for i in range(m) for elem in s])
            left = self.system_left_case_1(F);
            
            # Step 3) Solve the system. 
            # The system often doesn't have a solution. Then we pick another F.
            try:
                sol = left.solve_right(vector(right))
            except ValueError:
                continue

            # Step 4) Recover the error e. The solution ``sol`` is a vector of coefficients in
            # `\mathbb{F}_{q}` of the basis of F. Find `e \in \mathbb{F}_{q^m}^n`.
            if type(sol) == list and len(sol) > 1:      # For multiple solutions to the system.
                e_rec = [self.recover_e_prime(elem, F) for elem in sol]
                e_found = [elem for elem in e_rec if self.rank_weight_of_vector(elem) == w]
                if len(e_found) > 0:
                    return e_found
            else: # For a single solution to the system
                e_found = self.recover_e_prime(sol, F) 
                if self.rank_weight_of_vector(e_found) == w:
                    return e_found

        return "Could not find a solution to the RSD problem."
        
    def system_left_case_1(self, F):
        r"""
        Construct the left side of the system of case 1 of the GRS algorithm.
        This results in an `(n-k)m` by `nr` matrix with elements in 
        `\mathbb{F}_q`. 

        INPUT: 

        - ``F`` -- the subspace of `\mathbb{F}_{q^m}` of dimension r

        EXAMPLE::

            sage: q = 2; n = 5; k = 1; m = 4; w = 1; d = 2;                                   
            sage: r = m - ceil(k*m/n); r                                                                                  
            3
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix();                                          
            sage: F = C.find_subspace_of_dim(r); F # random                                                   
            Vector space of degree 4 and dimension 3 over Finite Field of size 2
            Basis matrix:
            [1 0 0 0]
            [0 1 1 0]
            [0 0 0 1]
            sage: left = C.system_left_case_1(F); left # random                                         
            [1 0 1 0 0 1 0 0 1 0 0 0 1 0 1]
            [1 0 0 1 0 0 1 0 1 1 0 0 1 0 1]
            [1 0 1 1 0 0 0 0 0 1 0 0 0 0 0]
            [0 0 0 1 0 1 0 0 0 0 0 1 0 0 1]
            [1 1 1 1 0 1 1 0 1 0 0 0 1 1 1]
            [0 1 0 0 1 0 1 1 1 0 1 0 1 1 1]
            [1 1 1 0 1 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 1 1 1 0 0 0 1 0 1 1 0 1]
            [0 0 0 0 1 0 0 1 0 0 0 0 0 0 0]
            [0 1 0 0 1 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 1 0 0 0 0 0 1 0 0 0 0]
            [0 0 0 0 0 0 0 0 0 0 1 0 0 1 0]
            [0 1 1 0 1 0 0 1 0 0 0 0 0 1 1]
            [0 0 1 0 0 1 0 1 1 0 0 1 0 1 1]
            [0 1 1 0 0 1 0 0 0 0 0 1 0 0 0]
            [0 0 0 0 1 1 0 0 0 0 1 0 0 1 0]
        """
        H = self.parity_check_matrix()
        m = self.extension_degree()
        n = self.length()
        k = self.dimension()
        r = F.dimension()

        basis = F.basis()
        vals = []
        for i in range(m):
            for u in range(n-k):
                vals += [[self.canonical_projection_base_element(H[u][l]*self.poly_form_of_vector(basis[j]),i) \
                 for l in range(n) for j in range(r)]]
        return matrix(vals) 

    def recover_e_prime(self, sol, F):
        r"""
        Recover `e'` from the solution given. `sol` is written as the coefficients
        of `e` in F.

        INPUT:

        - ``sol`` -- the solution written in terms of the basis of F

        - ``F`` -- the subspace of `\mathbb{F}_{q^m}` of dimension r

        EXAMPLES::

            sage: q = 2; m = 5; n = 5; k = 1; d = 2; r = 4
            sage: C = codes.LRPCCode(GF(q**m), n, k, d);
            sage: sol = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0)
            sage: F = C.find_subspace_of_dim(r); F # random                                                                                                            
            Vector space of degree 5 and dimension 4 over Finite Field of size 2
            Basis matrix:
            [1 0 0 0 1]
            [0 1 0 0 0]
            [0 0 1 0 1]
            [0 0 0 1 0]
            sage: C.recover_e_prime(sol, F) # random                                                                                        
            (0, 0, 0, 0, z5^4 + z5^2 + z5)
        """
        r = F.dimension()
        n = self.length()
        
        e_rec = matrix(n, sol)
        F_poly = matrix(r, 1, [self.poly_form_of_vector(row) for row in matrix(F.basis())])
        return vector(e_rec * F_poly) 

    def error_support_attack_case_2(self, w, r, y, beta = None):
        r"""
        Case 2 of the GRS algorithm, `m > n` in [AGHT2018]_. We search for a subspace `F` of `\mathbb{F}_{q}^{n}` of 
        dimension ``r`` which contains the row support of the error. 
        If this is the case, the solution found to the system of equations is of rank ``r`` and then 
        also a solution to the RSD problem.
        If this is not the case, we try another subspace `F` which we have not tried before.

        The complexity is `\mathcal{O}((n-k)^3 m^3 q^{wk})`.

        INPUT:

        - ``w`` -- the rank of the error vector

        - ``r`` -- the dimension of the subspace F

        - ``y`` -- the received vector

        - ``beta`` -- (default: ``None``) the basis of `\mathbb{F}_{q^m}` over `\mathbb{F}_q`,
          if not specified, take the standard basis `I_m`.

        EXAMPLES::

            sage: q = 2; n = 4; k = 1; m = 5; w = 1; d = 2;                               
            sage: r = n-k; r                                                                         
            3
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix();  
            sage: e = C.find_e_of_rank(w);   
            sage: x = C.random_element(); 
            sage: y = x + e;              
            sage: e_found = C.error_support_attack_case_2(w, r, y)                          
            Probability of choosing F correctly the first time is:  7/15
            The complexity of this attack is: O(2^13)
            sage: e_found == e                                                              
            True
        """
        from sage.functions.other import ceil
        from sage.functions.log import log
        from sage.calculus.var import var 
        from sage.symbolic.relation import solve_mod 
        from sage.rings.finite_rings.integer_mod import IntegerMod
        from sage.symbolic.relation import solve
        import re   

        q = self.sub_field().characteristic()
        n = self.length()
        k = self.dimension() 
        m = self.extension_degree()
        H = self.parity_check_matrix()

        complexity = (n-k)**3*m**3*q**(w*k)
        print("Probability of choosing F correctly the first time is: ", q_binomial(r, w, q)/q_binomial(n, w, q))
        print("The complexity of this attack is: O(2^" + str(ceil(log(complexity, 2).n())) + ")")

        number_of_Fs = q_binomial(n,r,q) 
        tried = []  # This list is used to keep track of the subspaces F that were already tried.

        if beta == None:
            beta = matrix.identity(m)

        basis = [self.poly_form_of_vector(row) for row in beta]
        hxij = [[self.rewrite_in_basis(H[x][i], beta) for i in range(n)] for x in range(n-k)] 
        s = H*y
        while len(tried) < number_of_Fs:
            # Step 1) Find a random subspace F of `\mathbb{F}_{q}^{n}` of dimension ``r``. 
            VS = VectorSpace(self.sub_field(), n)
            F = self.find_subspace_of_dim(r, VS) 
            tried += [F]

            # Step 2) Set up the system.
            right = vector([vector(elem)[i] for elem in s for i in range(m)])
            left = self.system_left_case_2(F, hxij, basis)
            
            # Step 3) Solve the system. 
            # The system often doesn't have a solution. Then we pick another F.
            try:
                sol = left.solve_left(right)
            except ValueError:
                continue

            # Step 4) Recover the error e. The solution ``sol`` is e rewritten in the basis ``beta``
            # and in subspace F.
            if type(sol) == list and len(sol) > 1: 
                for attempt in sol:
                    e_in_beta = matrix(m, r, attempt)*matrix(F.basis())
                    e_found = vector([self.poly_form_of_vector(e_in_beta.column(i)) for i in range(n)])
                    if(self.rank_weight_of_vector(e_found) == w):
                        return e_found
            else:
                e_in_beta = matrix(m, r, sol)*matrix(F.basis())
                e_found = vector([self.poly_form_of_vector(e_in_beta.column(i)) for i in range(n)])
                if(self.rank_weight_of_vector(e_found) == w):
                    return e_found

        return "Could not find a solution to the RSD problem."

    def system_left_case_2(self, F, hxij, basis):
        r"""
        Construct the left side of the system of case 2 of the GRS algorithm.
        This results in a matrix of dimension `mr` unknowns by `m(n-k)` equations.

        INPUT: 

        - ``F`` -- the subspace of `\mathcal{F}_{q^m}` of dimension r

        - ``hxij`` -- the parity check matrix `H` rewritten in the basis `\beta`

        - ``basis`` -- the polynomial form of the basis `\beta`

        EXAMPLE::

            sage: q = 2; n = 4; k = 1; m = 5; w = 1; r = 3; d = 2;               
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix(); e = C.find_e_of_rank(w);  
            sage: G = C.generator_matrix(); x = C.random_element(); 
            sage: y = x + e;  
            sage: VS = VectorSpace(C.sub_field(), n) 
            sage: F = C.find_subspace_of_dim(r, VS)       
            sage: s = H*e       
            sage: beta = matrix.identity(m)  
            sage: basis = [C.poly_form_of_vector(row) for row in beta] 
            sage: hxij = [[C.rewrite_in_basis(H[x][i], beta) for i in range(n)] for x in range(n-k)]  
            sage: eil = C.matrix_form_of_vector(e)  
            sage: f = matrix(F.basis())                                                                                       
            sage: left = C.system_left_case_2(F, hxij, basis); left # random                               
            [1 1 0 1 1 1 1 0 1 1 1 1 0 1 1]
            [0 0 0 0 0 1 0 0 0 0 0 1 0 1 1]
            [1 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
            [1 1 0 0 1 1 1 0 0 1 1 1 0 0 1]
            [0 0 0 0 0 0 1 0 0 0 1 0 0 0 1]
            [0 1 0 0 0 0 0 0 0 0 0 1 0 0 0]
            [1 1 0 0 0 1 1 0 0 0 1 1 0 0 0]
            [0 0 0 0 0 0 0 1 0 0 1 1 1 0 0]
            [0 0 1 0 0 0 0 0 0 0 0 0 1 0 0]
            [0 1 1 0 0 0 1 1 0 0 0 1 1 0 0]
            [0 0 0 0 0 0 0 0 1 0 0 1 1 1 0]
            [0 0 0 1 0 0 0 0 0 0 0 0 0 1 0]
            [0 0 1 1 0 0 0 1 1 0 0 0 1 1 0]
            [0 0 0 0 0 0 0 0 0 1 0 0 1 1 1]
            [0 0 0 0 1 0 0 0 0 0 0 0 0 0 1]
            sage: #This describes the following system:                                     
            sage: X = matrix(1, m*r, vector(C.get_var_matrix(m, r))); X
            [ x0  x1  x2  x3  x4  x5  x6  x7  x8  x9 x10 x11 x12 x13 x14]
            sage: vector(X * left) # random                                                     
            (x0 + x2 + x3 + x6, x0 + x3 + x5 + x6 + x9, x8 + x9 + x12, x0 + x11 + x12, x0 + x3 + x14,
            x0 + x1 + x3 + x6, x0 + x3 + x4 + x6 + x9, x7 + x9 + x12, x0 + x10 + x12, x0 + x3 + x13, 
            x0 + x2 + x3 + x4 + x6 + x7, x0 + x1 + x3 + x5 + x6 + x7 + x9 + x10, 
            x7 + x8 + x9 + x10 + x12 + x13, x0 + x1 + x10 + x11 + x12 + x13, x0 + x1 + x3 + x4 + x13 + x14)
        """
        n = self.length()
        k = self.dimension() 
        m = self.extension_degree()
        f = matrix(F.basis())
        r = F.dimension()

        c = []
        for t in range(n-k):   
            for l in range(m):
                b = []     
                for i in range(n):     
                    for j in range(m):     
                        temp = vector(basis[j]*basis[l]) 
                        if(hxij[t][i][j]) == 0:
                            b += [vector([0 for g in range(m*r)])]   
                            continue
                        temp2 = [hxij[t][i][j]*[elem*f.column(i) for elem in temp]]
                        if temp2 != [[]]:
                            b += temp2
                c += [b]
        left = [[0 for j in range(m*(n-k))] for i in range(m*r)]
        for x in range(len(c)): 
            temp2 = [] 
            for i in range(m): 
                temp = 0 
                for j in range(len(c[x])): 
                    temp += c[x][j][i] 
                temp2 += [temp] 
            for u in range(r): 
                for j in range(m): 
                    left[(n-k)*(x%m)+u][(x // m)*m + j] = temp2[j][u] 

        return matrix(left)


    def error_support_attack(self, w , r, y):
        r"""
        Breaks the Rank Syndrome Decoding problem via the GRS algorithm, as described in [AGHT2018]_.

        Given an input vector `y = x + e`, one can recover this `e` of certain rank ``r`` as `e_{found}` and decode to `x = y + e_{found}`.
        The algorithm consists of two cases: 
        1) `m \leq n`
        2) `m > n`
        On of these cases is called based on the values of `m` and `n`. 

        INPUT:

        - ``w`` -- the rank of the error vector

        - ``r`` -- the dimension of the subspace F

        - ``y`` -- the received vector

        EXAMPLES::

            sage: q = 2; n = 5; k = 1; m = 4; w = 1; d = 2;                                   
            sage: r = m - ceil(k*m/n); r                                                                                  
            3
            sage: C = codes.LRPCCode(GF(q**m), n, k, d); H = C.parity_check_matrix(); e = C.find_e_of_rank(w);  
            sage: x = C.random_element(); 
            sage: y = x + e;                       
            sage: e_found = C.error_support_attack_case_1(w, r, y)                                  
            Probability of choosing F correctly the first time is:  7/15
            The complexity of this attack is: O(2^13)
            sage: e_found == e                                                                       
            True
        """
        n = self.length()
        m = self.extension_degree()

        if m <= n:
            return self.error_support_attack_case_1(w , r, y)
        else:
            return self.error_support_attack_case_2(w , r, y)


class LinearRankMetricCode(AbstractLinearRankMetricCode):
    r"""
    Linear rank metric codes over a finite field, represented using a
    generator matrix.

    This class should be used for arbitrary and unstructured linear rank metric
    codes. This means that basic operations on the code, such as the computation
    of the minimum distance, will use generic, slow algorithms.

    If you are looking for constructing a code from a more specific family, see
    if the family has been implemented by investigating ``codes.<tab>``. These
    more specific classes use properties particular to that family to allow
    faster algorithms, and could also have family-specific methods.

    EXAMPLES::

        sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
        sage: C = codes.LinearRankMetricCode(G, GF(4))
        sage: C
        [3, 2] linear rank metric code over GF(64)/GF(4)
        sage: C.base_field()
        Finite Field in z6 of size 2^6
        sage: C.sub_field()
        Finite Field in z2 of size 2^2
        sage: C.length()
        3
        sage: C.dimension()
        2
        sage: C[2]
        (z6, z6, 0)
        sage: E = codes.encoders.LinearCodeGeneratorMatrixEncoder(C)
        sage: word = vector(C.base_field(), [1, 0])
        sage: E(word)
        (1, 1, 0)
    """

    def __init__(self, generator, sub_field=None, basis=None):
        r"""
        See the docstring for :meth:`LinearRankMetricCode`.

        INPUT:

        - ``generator`` -- a generator matrix over the ``base_field`` with
          dimension `k \times n`, where `k` is the dimension of the code and
          `n` its length; or a code over a finite field

        - ``sub_field`` -- (default: ``None``) the sub field of ``self``, if not
          specified, it is the prime field of ``base_field``

        - ``basis`` -- (default: ``None``) a basis of `F_{q^m}` as a vector space over
          ``sub_field``. If not specified, given that `q = p^s`, let
          `1,\beta,\ldots,\beta^{sm}` be the power basis that SageMath uses to
          represent `F_{q^m}`. The default basis is then `1,\beta,\ldots,\beta^{m-1}`.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4)) # indirect doctest
            sage: C
            [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        base_field = generator.base_ring()
        if not base_field.is_field():
            raise ValueError("'generator' must be defined on a field (not a ring)")

        if not sub_field:
            sub_field = base_field.prime_subfield()

        try:
            gen_basis = None
            if hasattr(generator,"nrows"): # generator matrix case
                if generator.rank() < generator.nrows():
                    gen_basis = generator.row_space().basis()
            else:
                gen_basis = generator.basis() # vector space etc. case
            if not gen_basis is None:
                generator = matrix(base_field, gen_basis)
                if generator.nrows() == 0:
                    raise ValueError("this linear code contains no non-zero vector")
        except AttributeError:
            # Assume input is an AbstractLinearRankMetricCode, extract its generator matrix
            generator = generator.generator_matrix()

        self._generator_matrix = generator
        self._length = generator.ncols()
        super(LinearRankMetricCode, self).__init__(base_field, sub_field, self._length, "GeneratorMatrix", "NearestNeighbor", basis)

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C
            [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        R = self.base_field()
        S = self.sub_field()
        if R and S in Fields():
            return "[%s, %s] linear rank metric code over GF(%s)/GF(%s)"%(self.length(), self.dimension(), R.cardinality(), S.cardinality())
        else:
            return "[%s, %s] linear rank metric code over %s/%s"%(self.length(), self.dimension(), R, S)

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: latex(C)
            [3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}/\Bold{F}_{2^{2}}
        """
        return "[%s, %s]\\textnormal{ Linear rank metric code over }%s/%s"\
                % (self.length(), self.dimension(), self.base_field()._latex_(), self.sub_field()._latex_())

    def generator_matrix(self, encoder_name=None, **kwargs):
        r"""
        Return a generator matrix of ``self``.

        INPUT:

        - ``encoder_name`` -- (default: ``None``) name of the encoder which will be
          used to compute the generator matrix. ``self._generator_matrix``
          will be returned if default value is kept.

        - ``kwargs`` -- all additional arguments are forwarded to the construction of the
          encoder that is used.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: C.generator_matrix()
            [1 1 0]
            [0 0 1]
        """
        if encoder_name is None or encoder_name == 'GeneratorMatrix':
            g = self._generator_matrix
        else:
            g = super(LinearRankMetricCode, self).generator_matrix(encoder_name, **kwargs)
        g.set_immutable()
        return g


####################### decoders ###############################
class LinearRankMetricCodeNearestNeighborDecoder(Decoder):
    r"""
    Construct a decoder for Linear Rank Metric Codes.

    This decoder will decode to the nearest codeword found.
    """

    def __init__(self, code):
        r"""

        INPUT:

        - ``code`` -- A code associated to this decoder

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        super(LinearRankMetricCodeNearestNeighborDecoder, self).__init__(code, code.ambient_space(), \
                code._default_encoder_name)

    def __eq__(self, other):
        r"""
        Test equality between LinearRankMetricCodeNearestNeighborDecoder objects.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D1 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D2 = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D1 == D2
            True
        """
        return isinstance(other, LinearRankMetricCodeNearestNeighborDecoder)\
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Return a string representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D
            Nearest neighbor decoder for [3, 2] linear rank metric code over GF(64)/GF(4)
        """
        return "Nearest neighbor decoder for %s" % self.code()

    def _latex_(self):
        r"""
        Return a latex representation of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(64), [[1,1,0], [0,0,1]])
            sage: C = codes.LinearRankMetricCode(G, GF(4))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: latex(D)
            \textnormal{Nearest neighbor decoder for }[3, 2]\textnormal{ Linear rank metric code over }\Bold{F}_{2^{6}}/\Bold{F}_{2^{2}}
        """
        return "\\textnormal{Nearest neighbor decoder for }%s" % self.code()._latex_()

    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``word`` and returns a codeword.

        INPUT:

        - ``r`` -- a codeword of ``self``

        OUTPUT:

        - a vector of ``self``'s message space

        EXAMPLES::

            sage: F.<a> = GF(4)
            sage: G = Matrix(F, [[1,1,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decode_to_code(vector(F, [a, a, 1]))
            (a, a, 0)
        """
        C = self.code()
        c_min = C.zero()
        h_min = C.rank_weight_of_vector(r)
        for c in C:
            if C.rank_weight_of_vector(c-r) < h_min:
                h_min = C.rank_weight_of_vector(c-r)
                c_min = c
        c_min.set_immutable()
        return c_min

    def decoding_radius(self):
        r"""
        Return maximal number of errors ``self`` can decode.

        EXAMPLES::

            sage: F.<a> = GF(8)
            sage: G = Matrix(F, [[1,a,a^2,0]])
            sage: C = codes.LinearRankMetricCode(G, GF(2))
            sage: D = codes.decoders.LinearRankMetricCodeNearestNeighborDecoder(C)
            sage: D.decoding_radius()
            1
        """
        return (self.code().minimum_distance()-1) // 2

####################### registration ###############################

LinearRankMetricCode._registered_encoders["GeneratorMatrix"] = LinearCodeGeneratorMatrixEncoder
