r"""
Cyclic Code

Let `F` be a field. A `[n, k]` code `C` over `F` is called cyclic if every
cyclic shift of a codeword is also a codeword [R06]_:

    .. MATH::

        \forall c \in C,
        c = (c_{0}, c_{1}, \dots , c_{n-1}) \in C
        \Rightarrow (c_{n-1}, c_{0}, \dots , c_{n-2}) \in C

Let `c = (c_0, c_1, \dots, c_{n-1})` be a codeword of `C`.
This codeword can be seen as a polynomial over `F_q[x]` as follows:
`\Sigma_{i=0}^{n-1} c_i x^i`.
There is a unique monic polynomial `g(x)` such that for every
`c(x) \in F_q[x]` of degree less than `n-1`, we have
`c(x) \in C \Leftrightarrow g(x) | c(x)`.
This polynomial is called the generator polynomial of `C`.

For now, only single-root cyclic codes (i.e. whose length `n` and field order
`q` are coprimes) are implemented.


REFERENCES:

- [R06]_

TESTS:

This class uses the following experimental feature:
:class:`sage.coding.relative_finite_field_extension.RelativeFiniteFieldExtension`.
This test block is here only to trigger the experimental warning so it does not
interferes with doctests::

    sage: from sage.coding.relative_finite_field_extension import *
    sage: Fqm.<aa> = GF(16)
    sage: Fq.<a> = GF(4)
    sage: RelativeFiniteFieldExtension(Fqm, Fq)
    doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See http://trac.sagemath.org/20284 for details.
    Relative field extension between Finite Field in aa of size 2^4 and Finite Field in a of size 2^2
"""

#*****************************************************************************
#       Copyright (C) 2015 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from linear_code import (AbstractLinearCode,
                         LinearCodeSyndromeDecoder,
                         LinearCodeNearestNeighborDecoder)
from encoder import Encoder
from decoder import Decoder
from sage.rings.integer import Integer
from sage.arith.all import gcd
from sage.modules.free_module_element import vector
from sage.matrix.constructor import matrix
from sage.misc.cachefunc import cached_method
from sage.misc.misc_c import prod
from sage.rings.polynomial.polynomial_ring import PolynomialRing_general
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.all import Zmod
from sage.functions.log import log
from sage.categories.homset import Hom
from copy import copy
from sage.groups.generic import discrete_log
from sage.misc.functional import multiplicative_order
from relative_finite_field_extension import RelativeFiniteFieldExtension

def find_generator_polynomial(code, check=True):
    r"""
    Returns a possible generator polynomial for ``code``.

    If the code is cyclic, the generator polynomial is the gcd of all the
    polynomial forms of the codewords. Conversely, if this gcd exactly generates
    the code ``code``, then ``code`` is cyclic.

    If ``check`` is set to ``True``, then it also checks that the code is indeed
    cyclic. Otherwise it doesn't.

    INPUT:

    - ``code`` -- a linear code

    OUTPUT:

    - the generator polynomial (if the code is cyclic).

    EXAMPLES::

        sage: from sage.coding.cyclic_code import find_generator_polynomial
        sage: C = codes.GeneralizedReedSolomonCode(GF(2^3, 'a').list()[1:2^3], 2^2)
        sage: find_generator_polynomial(C)
        x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1
    """
    G = code.generator_matrix()
    n = code.length()
    k = code.dimension()
    F = code.base_ring()
    R = F['x']
    x = R.gen()
    g = R(G.row(0).list())

    i = 1
    while(g.degree() > n - k and i < k):
        p = R(G.row(i).list())
        g = gcd(g, p)
        i += 1

    if check:
        if g.degree() != n - k :
            raise ValueError("The code is not cyclic.")
        i = 0
        c = g.coefficients(sparse = False)
        c = _complete_list(c, n)
        while (vector(c) in code and i < k):
            c = [c[n-1]] + c[:n-1]
            i += 1
        if (i != k):
            raise ValueError("The code is not cyclic.")
    return g.monic()

def _complete_list(l, length):
    r"""
    Returns ``l`` with as many zeros appened to it as necessary to make
    it a list of size ``length``.

    INPUT:

    - ``l`` -- a list

    - ``length`` -- an integer

    OUTPUT:

    - a list formed by ``l`` and its completion to ``length``

    EXAMPLES::

        sage: l = [1, 2, 3, 4]
        sage: sage.coding.cyclic_code._complete_list(l, 6)
        [1, 2, 3, 4, 0, 0]
    """
    F = l[0].base_ring()
    l_len = len(l)
    if(l_len != length):
        l = l + ([F.zero()] * (length - l_len))
    return l

def bch_bound(n, D, arithmetic = False, bch_parameters = False):
    r"""
    Returns the BCH bound obtained for a cyclic code of length ``n`` and defining set ``D``.

    Considering a cyclic code `C`, with defining set `D`, length `n`, and minimum
    distance `d`. We have the following bound, called BCH bound, on `d`:
    `d \geq \delta + 1`, where `\delta` is the length of the longest arithmetic sequence
    (modulo `n`) of elements in `D`.

    That is, if `\exists c, \gcd(c,n) = 1` such that
    `\{l, l+c, \dots, l + (\delta - 1) \times c\} \subseteq D`,
    then `d \geq \delta + 1` [1]

    The BCH bound is often known in the particular case `c = 1`. The user can
    specify by setting ``arithmetic = False``.

    .. NOTE::

        As this is a specific use case of the BCH bound, it is *not* available if the global namespace.
        Call it by using ``sage.coding.cyclic_code.bch_bound``. You can also load it into the global
        namespace by typing ``from sage.coding.cyclic_code import bch_bound``.

    INPUT:

    - ``n`` -- an integer

    - ``D`` -- a list of integers

    - ``arithmetic`` -- (default: ``False``), if it is set to ``True``, then it
      computes the BCH bound using the longest arithmetic sequence definition

    - ``bch_parameters`` -- (default: ``False``), if it is set to ``True``, then
      the function also returns the parameters of the smallest BCH code that 
      contains a cyclic code having `D` as defining set. Only works when
      ``arithmetic = True``

    EXAMPLES::

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D)
        3

        sage: n = 15
        sage: D = [14,1,2,11,12]
        sage: sage.coding.cyclic_code.bch_bound(n, D, True)
        4
    """
    # Idea of the algorithm:
    # We first store a table representing D;
    # then, for each arithmetic step i:
    #   - we compute the size of the current sequence in D by
    #     incrementing a counter (cur_length)
    #   - we update the maximum size (max_length) if needed
    #   - we also keep the first and the last of these lengths in
    #     order to deal with a possible overlapping sequence (i.e.
    #     which contains 0)
    # Note: this algorithm has O(n^2) time and O(n) space complexity, assuming
    # that accessing a table is O(1)
    table_D = [ 1 if d in D else 0 for d in range(n) ]
    if arithmetic:
        max_loops = n
    else:
        max_loops = 2
    max_length = 0
    for i in range(1, max_loops):
        cur_length = 0
        if gcd(i, n) == 1:
            j = 0
            first_length = 0
            while table_D[i*j % n] == 1:
                first_length += 1
                j += 1
            cur_length = first_length
            if cur_length > max_length:
                max_length = cur_length
                bch_parameters = [0, j]
            cur_length = 0
            begin = j+1
            for j in range(begin, n):
                if table_D[i*j % n] == 1:
                    cur_length += 1
                else:
                    if cur_length > max_length:
                        max_length = cur_length
                        bch_parameters = [(j - cur_length * i) % n, i]
                    cur_length = 0
            overlap_length = first_length + cur_length
            if overlap_length > max_length:
                max_length = overlap_length
                bch_parameters = [ (- first_length * i) % n, i]
    if bch_parameters:
        return max_length + 1, bch_parameters
    return max_length + 1

class CyclicCode(AbstractLinearCode):
    r"""
    Representation of a cyclic code.

    We propose three different ways to create a new CyclicCode, either by providing:

    - the generator polynomial and the length (1)
    - an existing linear code. In that case, a generator polynomial will be computed
       from the provided linear code's parameters (2)
    - (a subset of) the defining set of the cyclic code (3)

    For now, only single-root cyclic codes are implemented. That is, only cyclic
    codes such that its length `n` and field order `q` are coprimes.

    Depending on which behaviour you want, you need to specify the names of the arguments to
    CyclicCode. See EXAMPLES section below for details.

    INPUT:

    - ``generator_pol`` -- (default: ``None``) the generator polynomial
      of ``self``. That is, the highest-degree monic polynomial which divides
      every polynomial representation of a codeword in ``self``.

    - ``length`` -- (default: ``None``) the length of ``self``. It has to be
      bigger than the degree of ``generator_pol``.

    - ``code`` -- (default: ``None``) a linear code.

    - ``check`` -- (default: ``False``) a boolean representing whether the
      cyclicity of ``self`` must be checked while finding the generator
      polynomial. See :meth:`find_generator_polynomial` for details.

    - ``D`` -- (default: ``None``) a list of integers between ``0`` and
      ``length-1``, corresponding to (a subset of) the defining set of the code.
      Will be modified if it is not cyclotomic-closed.

    - ``field`` -- (default: ``None``) the base field of ``self``.

    - ``primitive_root`` -- (default: ``None``) the primitive root of
      the splitting field which contains the roots of the generator polynomial.
      It has to be of multiplicative order ``length`` over this field.
      If the splitting field is not ``field``, it also have to be a polynomial
      in ``zx``, where ``x`` is the degree of the extension over the prime field.
      For instance, over ``GF(16)``, it has to be a polynomial in ``z4``.

    EXAMPLES:

    We can construct a CyclicCode object using three different methods.
    First (1), we provide a generator polynomial and a code length::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: C
        [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

    We can also provide a code (2). In that case, the program will try to
    extract a generator polynomial (see :meth:`find_generator_polynomial`
    for details)::

        sage: C = codes.GeneralizedReedSolomonCode(GF(2 ** 3, 'a').list()[1:2 ** 3], 2 ** 2)
        sage: Cc = codes.CyclicCode(code = C)
        sage: Cc
        [7, 4] Cyclic Code over Finite Field in a of size 2^3 with x^3 + (a^2 + 1)*x^2 + a*x + a^2 + 1 as generator polynomial

    Finally, we can give (a subset of) a defining set for the code (3).
    In this case, the generator polynomial will be computed::

        sage: F = GF(16, 'a')
        sage: n = 15
        sage: Cc = codes.CyclicCode(length = n, field = F, D = [1,2])
        sage: Cc
        [15, 13] Cyclic Code over Finite Field in a of size 2^4 with x^2 + (a^2 + a)*x + a^3 as generator polynomial
    """

    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, generator_pol=None, length=None, code=None, check=True, D = None, field = None, primitive_root = None):
        r"""
        TESTS:

        If one provides a generator polynomial and a length, we check that
        the length is bigger than the degree of the polynomial::

            sage: F.<x> = GF(2)[]
            sage: n = 2
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Only cyclic codes whose length and field order are coprimes are implemented.

        We also check that the polynomial is defined over a finite field::

            sage: F.<x> = RR[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Generator polynomial must be defined over a finite field

        And we check that the generator polynomial divides `x^{n} - 1`,
        where `n` is provided length::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 2 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            Traceback (most recent call last):
            ...
            ValueError: Provided polynomial must divide x^n - 1, where n is the provided length

        In the case of a code is passed as argument, if it's not possible
        to extract a generator polynomial, an exception is raised::

            sage: G = matrix(GF(2), [[1, 1, 1], [0, 1, 1]])
            sage: C = codes.LinearCode(G)
            sage: Cc = codes.CyclicCode(code = C)
            Traceback (most recent call last):
            ...
            ValueError: The code is not cyclic.
        """
        # Case (1) : generator polynomial and length are provided.
        if generator_pol is not None and length is not None:
            F = generator_pol.base_ring()
            if not F.is_finite() or not F.is_field():
                raise ValueError("Generator polynomial must be defined over a finite field")
            q = F.cardinality()
            if not gcd(length, q) == 1:
                raise ValueError("Only cyclic codes whose length and field order are coprimes are implemented.")
            R = generator_pol.parent()
            deg = generator_pol.degree()
            if not isinstance(length, Integer):
                length = Integer(length)
            if not generator_pol.divides(R.gen() ** length - 1):
                raise ValueError("Provided polynomial must divide x^n - 1, where n is the provided length")
            self._polynomial_ring = R
            self._dimension = length - deg
            if not generator_pol.is_monic():
                self._generator_polynomial = generator_pol.monic()
            else:
                self._generator_polynomial = generator_pol
            super(CyclicCode, self).__init__(F, length, "Vector", "Syndrome")

        # Case (2) : a code is provided.
        elif code is not None:
            if not isinstance(code, AbstractLinearCode):
                raise ValueError("code must be an AbstractLinearCode")
            F = code.base_ring()
            q = F.cardinality()
            n = code.length()
            if not gcd(n, q) == 1:
                raise ValueError("Only cyclic codes whose length and field order are coprimes are implemented.")
            try:
                g = find_generator_polynomial(code, check)
            except ValueError, e:
                raise ValueError(e)
            self._polynomial_ring = g.parent()
            self._generator_polynomial = g
            self._dimension = code.dimension()
            super(CyclicCode, self).__init__(code.base_ring(), n, "Vector", "Syndrome")

        # Case (3) : a defining set, a length and a field are provided
        elif D is not None and length is not None and field is not None:
            F = field
            if not F.is_finite() or not F.is_field():
                raise ValueError("A finite field must be given in complement to defining set.")
            n = length
            q = F.cardinality()
            if not gcd(n, q) == 1:
                raise ValueError("Only cyclic codes whose length and field order are coprimes are implemented.")
            
            R = F['x']
            x = R.gen()
            s = 1
            while not (q ** s - 1) % n == 0:
                s += 1

            if s == 1: # splitting field is F
                Fsplit = F
            else: # must compute a splitting field
                Fsplit, F_to_Fsplit = F.extension(Integer(s), map = True)
                FE = RelativeFiniteFieldExtension(Fsplit, F, embedding = F_to_Fsplit)
                Rsplit = Fsplit['xx']
                xx = Rsplit.gen()

            if primitive_root is not None:
                if (primitive_root not in Fsplit or
                    multiplicative_order(primitive_root) != n):
                    raise ValueError("primitive_root has to be an element of multiplicative order n in the extension field used to compute the generator polynomial")
                else:
                    alpha = primitive_root
            else:
                alpha = Fsplit.zeta(n)

            cosets = Zmod(n).cyclotomic_cosets(q, D)
            pows = [ item for l in cosets for item in l ]

            g = R.one()
            for J in cosets:
                if s == 1:
                    g *= x - alpha ** J[0]
                else:
                    pol = Rsplit.one()
                    for j in J:
                        pol *= xx - alpha**j
                    g *= R([FE.cast_into_relative_field(coeff) for coeff in pol])
            
            # we set class variables
            self._primitive_root = alpha
            self._defining_set = pows
            self._defining_set.sort()
            self._polynomial_ring = R
            self._generator_polynomial = g
            self._dimension = n - g.degree()
            super(CyclicCode, self).__init__(F, n, "Vector", "Syndrome")

        else:
            raise AttributeError("You must provide either a code, or a list of powers and the length and the field, or a generator polynomial and the code length")

    def __contains__(self, word):
        r"""
        Returns ``True`` if ``word`` belongs to ``self``, ``False`` otherwise.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: c in C
            True
        """
        g = self.generator_polynomial()
        R = self._polynomial_ring
        return (g.divides(R(word.list())) and word in self.ambient_space())

    def __eq__(self, other):
        r"""
        Tests equality between CyclicCode objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C1 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C2 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C1 == C2
            True
        """
        if not isinstance(other, CyclicCode):
            return False
        else:
            R = self._polynomial_ring
            return self.base_field() == other.base_field() \
                and self.length() == other.length() \
                and self.generator_polynomial() == R(other.generator_polynomial())

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C
            [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "[%s, %s] Cyclic Code over %s with %s as generator polynomial"\
                % (self.length(), self.dimension(),\
                self.base_field(), self.generator_polynomial())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: latex(C)
            [7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "[%s, %s] \\textnormal{ Cyclic Code over } %s \\textnormal{ with } %s \\textnormal{ as generator polynomial}"\
                % (self.length(), self.dimension(),\
                self.base_field()._latex_(), self.generator_polynomial()._latex_())

    def generator_polynomial(self):
        r"""
        Returns the generator polynomial of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.generator_polynomial()
            x^3 + x + 1
        """
        return self._generator_polynomial

    def defining_set(self, primitive_root = None):
        r"""
        Returns the set of powers of the root of ``self``'s generator polynomial
        over the extension field. It depends on the choice of the primitive
        root of the splitting field.

        EXAMPLES:

        We provide a defining set at construction time::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: C = codes.CyclicCode(length = n, field = F, D = [1,2])
            sage: C.defining_set()
            [1, 2]

        If the defining set was provided by the user, it might have been expanded
        at construction time. In this case, the expanded defining set will be returned::

            sage: C = codes.CyclicCode(length = 13, field = F, D = [1, 2])
            sage: C.defining_set()
            [1, 2, 3, 5, 6, 9]

        If a generator polynomial was passed at construction time,
        the defining set is computed using this polynomial::

            sage: F.<x> = GF(8, 'a')[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.defining_set()
            [1, 2, 4]

        Both operations give the same result::

            sage: C1 = codes.CyclicCode(length = n, field = GF(8, 'a'), D = [1, 2, 4])
            sage: C1.generator_polynomial() == g
            True

        Another one, in the revert order::

            sage: F = GF(16, 'a')
            sage: n = 13
            sage: C1 = codes.CyclicCode(length = n, field = F, D = [1, 2])
            sage: g = C1.generator_polynomial()
            sage: C2 = codes.CyclicCode(generator_pol = g, length = n)
            sage: C1.defining_set() == C2.defining_set()
            True
        """
        if (hasattr(self, "_defining_set") and
            (primitive_root is None or primitive_root == self._primitive_root)):
            return self._defining_set
        else:
            F = self.base_field()
            n = self.length()
            q = F.cardinality()
            g = self.generator_polynomial()

            s = 1
            while not (q ** s - 1) % n == 0:
                s += 1

            if primitive_root is None:
                Fsplit, F_to_Fsplit = F.extension(Integer(s), map = True)
                alpha = Fsplit.zeta(n)
            else:
                try:
                    alpha = primitive_root
                    Fsplit = alpha.parent()
                    FE = RelativeFiniteFieldExtension(Fsplit, F)
                    F_to_Fsplit = FE.embedding()
                except ValueError:
                    raise ValueError("primitive_root does not belong to the right splitting field")
                if alpha.multiplicative_order() != n:
                    raise ValueError("primitive_root must have multiplicative order n")
            self._primitive_root = alpha

            Rsplit = Fsplit['xx']
            gsplit = Rsplit([F_to_Fsplit(coeff) for coeff in g])
            roots = gsplit.roots(multiplicities = False)
            D = [discrete_log(root, alpha) for root in roots]
            D.sort()
            self._defining_set = D
            return D

    def primitive_root(self):
        r"""
        Returns the primitive root of the splitting field that is used
        to build the defining set of the code.

        If it has not been specified by the user, it is set by default with the
        output of the ``zeta`` method of the splitting field.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.primitive_root()
            z3

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: a = F.gen()
            sage: Cc = codes.CyclicCode(length = n, field = F, D = [1,2], primitive_root = a^2 + 1)
            sage: Cc.primitive_root()
            a^2 + 1
        """
        if hasattr(self, "_primitive_root"):
            return self._primitive_root
        else:
            _ = self.defining_set()
            return self._primitive_root

    def check_polynomial(self):
        r"""
        Returns the check polynomial of ``self``.

        Let `C` be a cyclic code of length `n` and `g` its
        generator polynomial.
        The following: `h = \frac{x^n - 1}{g(x)}` is called `C`'s
        check polynomial.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: h = C.check_polynomial()
            sage: h == (x**n - 1)/C.generator_polynomial()
            True
        """
        if hasattr(self, "_check_polynomial"):
            return self._check_polynomial
        R = self._polynomial_ring
        n = self.length()
        self._check_polynomial = (R.gen() ** n - 1) // self.generator_polynomial()
        return self._check_polynomial

    @cached_method
    def parity_check_matrix(self):
        r"""
        Returns the parity check matrix of ``self``.

        The parity check matrix of a linear code `C` corresponds to the
        generator matrix of the dual code of `C`.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: C.parity_check_matrix()
            [1 0 1 1 1 0 0]
            [0 1 0 1 1 1 0]
            [0 0 1 0 1 1 1]
        """
        k = self.dimension()
        n = self.length()
        R = self._polynomial_ring
        h = self.check_polynomial()
        l = h.coefficients(sparse = False)
        l.reverse()
        l = _complete_list(l, n)
        H = matrix(self.base_ring(), n - k, n)
        for i in range(n - k):
            H.set_row(i, l)
            l = l[n-1:] + l[:n-1]
        return H

    def bch_bound(self, arithmetic = False, bch_parameters = False):
        r"""
        Returns the BCH bound of ``self`` which is a bound on ``self``'s minimum distance.

        See :meth:`sage.coding.cyclic_code.bch_bound` for details.

        INPUT:

        - ``arithmetic`` -- (default: ``False``), if it is set to ``True``, then it
          computes the BCH bound using the longest arithmetic sequence definition

        - ``bch_parameters`` -- (default: ``False``), if it is set to ``True``, then
          the function also returns the parameters of the smallest BCH code that 
          contains a cyclic code having `D` as defining set. Only works when
          ``arithmetic = True``

        EXAMPLES::

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field = F, length = n, D = D)
            sage: C.bch_bound()
            3

            sage: F = GF(16, 'a')
            sage: n = 15
            sage: D = [14,1,2,11,12]
            sage: C = codes.CyclicCode(field = F, length = n, D = D)
            sage: C.bch_bound(True)
            4
        """
        return bch_bound(n = self.length(), D = self.defining_set(), arithmetic = arithmetic,\
                bch_parameters = bch_parameters)










class CyclicCodePolynomialEncoder(Encoder):
    r"""
    An encoder encoding polynomials into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    This encoder encodes any polynomial `p \in F[x]_{<k}` by computing
    `c = p \times g` and returning the vector of its coefficients.

    INPUT:

    - ``code`` -- The associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
        sage: E
        Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
    """

    def __init__(self, code):
        r"""
        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

        """
        if not isinstance(code, CyclicCode):
            raise ValueError("code has to be a CyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super(CyclicCodePolynomialEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between CyclicCodePolynomialEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E1 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E2 = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, CyclicCodePolynomialEncoder) \
            and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E
            Polynomial-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "Polynomial-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: latex(E)
            \textnormal{Polynomial-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "\\textnormal{Polynomial-style encoder for }%s" % self.code()._latex_()

    def encode(self, p):
        r"""
        Transforms ``p`` into an element of the associated code of ``self``.

        INPUT:

        - ``p`` -- A polynomial from ``self`` message space

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: m = x ** 2 + 1
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        C = self.code()
        k = C.dimension()
        n = C.length()
        if p.degree() >= k:
            raise ValueError("Degree of the message must be at most %s" % k-1)
        res = (p * C.generator_polynomial()).coefficients(sparse = False)
        res = _complete_list(res, n)
        return vector(C.base_field(), res)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            x^2 + 1
        """
        R = self.message_space()
        g = self.code().generator_polynomial()
        p = R(c.list())
        return p // g

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodePolynomialEncoder(C)
            sage: E.message_space()
            Univariate Polynomial Ring in x over Finite Field of size 2 (using NTL)
        """
        return self._polynomial_ring










class CyclicCodeVectorEncoder(Encoder):
    r"""
    An encoder which can encode vectors into codewords.

    Let `C` be a cyclic code over some finite field `F`,
    and let `g` be its generator polynomial.

    Let `m = (m_1, m_2, \dots, m_k)` be a vector in `F^{k}`.
    This codeword can be seen as a polynomial over `F[x]`, as follows:
    `P_m = \Sigma_{i=0}^{k-1} m_i \times x^i`.

    To encode `m`, this encoder does the following multiplication: `P_m \times g`.

    INPUT:

    - ``code`` -- The associated code of this encoder

    EXAMPLES::

        sage: F.<x> = GF(2)[]
        sage: n = 7
        sage: g = x ** 3 + x + 1
        sage: C = codes.CyclicCode(generator_pol = g, length = n)
        sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
        sage: E
        Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
    """

    def __init__(self, code):
        r"""

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial

        """
        if not isinstance(code, CyclicCode):
            raise ValueError("code has to be a CyclicCode")
        self._polynomial_ring = code._polynomial_ring
        super(CyclicCodeVectorEncoder, self).__init__(code)

    def __eq__(self, other):
        r"""
        Tests equality between CyclicCodeVectorEncoder objects.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E1 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E2 = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E1 == E2
            True
        """
        return isinstance(other, CyclicCodeVectorEncoder) \
                and self.code() == other.code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E
            Vector-style encoder for [7, 4] Cyclic Code over Finite Field of size 2 with x^3 + x + 1 as generator polynomial
        """
        return "Vector-style encoder for %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: latex(E)
            \textnormal{Vector-style encoder for }[7, 4] \textnormal{ Cyclic Code over } \Bold{F}_{2} \textnormal{ with } x^{3} + x + 1 \textnormal{ as generator polynomial}
        """
        return "\\textnormal{Vector-style encoder for }%s" % self.code()._latex_()

    def encode(self, m):
        r"""
        Transforms ``m`` into an element of the associated code of ``self``.

        INPUT:

        - ``m`` -- A element from ``self``'s message space

        OUTPUT:

        - A codeword in associated code of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: m = vector(GF(2), (1, 0, 1, 0))
            sage: E.encode(m)
            (1, 1, 1, 0, 0, 1, 0)
        """
        if hasattr(self, "_known_generator_matrix"):
            return super(CyclicCodeVectorEncoder, self).encode(m)

        k = self.code().dimension()
        n = self.code().length()
        F  = self.code().base_field()
        R = self._polynomial_ring
        p = R(m.list())
        if p.degree() >= k:
            raise ValueError("Degree of the message must be at most %s" % k-1)
        res = (p * self.code().generator_polynomial()).coefficients(sparse = False)
        res = _complete_list(res, n)
        return vector(F, res)

    def unencode_nocheck(self, c):
        r"""
        Returns the message corresponding to ``c``.
        Does not check if ``c`` belongs to the code.

        INPUT:

        - ``c`` -- A vector with the same length as the code

        OUTPUT:

        - An element of the message space

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: c = vector(GF(2), (1, 1, 1, 0, 0, 1, 0))
            sage: E.unencode_nocheck(c)
            (1, 0, 1, 0)
        """

        R = self._polynomial_ring
        g = self.code().generator_polynomial()
        p = R(c.list())
        l = (p//g).coefficients(sparse = False)
        l = _complete_list(l, self.message_space().dimension())
        return vector(self.code().base_field(), l)

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.generator_matrix()
            [1 1 0 1 0 0 0]
            [0 1 1 0 1 0 0]
            [0 0 1 1 0 1 0]
            [0 0 0 1 1 0 1]
        """

        self._known_generator_matrix = True
        C = self.code()
        k = C.dimension()
        n = C.length()
        l = C.generator_polynomial().coefficients(sparse = False)
        l = _complete_list(l, n)
        G = matrix(C.base_ring(), k, n)
        for i in range(k):
            G.set_row(i, l)
            l = l[n-1:] + l[:n-1]
        return G

    def message_space(self):
        r"""
        Returns the message space of ``self``

        EXAMPLES::

            sage: F.<x> = GF(2)[]
            sage: n = 7
            sage: g = x ** 3 + x + 1
            sage: C = codes.CyclicCode(generator_pol = g, length = n)
            sage: E = codes.encoders.CyclicCodeVectorEncoder(C)
            sage: E.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.code().base_ring() ** self.code().dimension()


####################### registration ###############################

CyclicCode._registered_encoders["Polynomial"] = CyclicCodePolynomialEncoder
CyclicCode._registered_encoders["Vector"] = CyclicCodeVectorEncoder
CyclicCode._registered_decoders["Syndrome"] = LinearCodeSyndromeDecoder
CyclicCode._registered_decoders["NearestNeighbor"] = LinearCodeNearestNeighborDecoder
