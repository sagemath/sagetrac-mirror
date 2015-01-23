# -*- coding: utf-8 -*-
r"""
Double Square tiles

If a polyomino P tiles the plane by translation, then there exists a
regular tiling of the plane by P [WVL1984]_, i.e., where the set of
translations forms a lattice. Such a polyomino was called *exact* by
Wijshoff and van Leeuven. There are two types of regular tiling of the
plane : square and hexagonal. These are characterized by the
Beauquier-Nivat condition [BN1991]_. Deciding whether a polyomino is exact
can be done efficiently from the boundary and in linear time for square
tiling [BFP2009]_. Brlek, Fédou, Provençal also remarked that there exist
polyominoes leading to more than one regular tilings but conjectured that
any polyomino produces at most two regular square tilings. This conjecture
was proved in [BBL2012]_. In [BBGL2011]_, two infinite families of *double
square tiles* were provided, that is polyominoes having exactly two
distinct regular square tilings of the plane, namely the Christoffel tiles
and the Fibonacci tiles. Finally, in [BGL2012]_, it was shown that any
double square tile can be constructed using two simple combinatorial rules:
EXTEND and SWAP.

This module is about double square tiles. Notations are chosen according to
[BGL2012]_. It allows to construct, study and show double square tiles.
Operations TRIM, SWAP and EXTEND are implemented. Double square tiles can
be shown using Sage 2D Graphics objects or using tikz.

REFERENCES:

.. [WVL1984] Wijshoff, H. A. G, et J. Van Leeuwen. Arbitrary versus periodic
   storage schemes and tessellations of the plane using one type of
   polyomino. INFORM. AND CONTROL 62 (1984): 1‑25.

.. [BN1991] Beauquier, D., and M. Nivat. On translating one polyomino to tile the
   plane. Discrete & Computational Geometry 6 (1991): 575‑592.
   :doi:`10.1007/BF02574705`

.. [BFP2009] S. Brlek, J.-M Fédou, X. Provençal, On the Tiling by Translation
   Problem, Discrete Applied Mathematics 157 Issue 3 (2009) 464-475.
   :doi:`10.1016/j.dam.2008.05.026`

.. [BBL2012] A. Blondin Massé, S. Brlek, S. Labbé, A parallelogram tile fills the
   plane by translation in at most two distinct ways, Discrete Applied
   Mathematics 160 (2012) 1011-1018. :doi:`10.1016/j.dam.2011.12.023`

.. [BBGL2011] A. Blondin Massé, S. Brlek, A. Garon, S. Labbé, Two
   infinite families of polyominoes that tile the plane by translation in two
   distinct ways, Theoret. Comput. Sci. 412 (2011) 4778-4786.
   :doi:`10.1016/j.tcs.2010.12.034`

.. [BGL2012] A. Blondin Massé, A. Garon, S. Labbé, Combinatorial properties
   of double square tiles, Theoretical Computer Science, Available online 2
   November 2012. :doi:`10.1016/j.tcs.2012.10.040`

AUTHORS:

    - Sébastien Labbé, 2008: initial version
    - Alexandre Blondin Massé, 2008: initial version
    - Sébastien Labbé, March 2013: rewrite for inclusion into Sage

EXAMPLES:

Double Square tile from the boundary word of a known double square::

    sage: from sage.combinat.double_square_tile import DoubleSquare
    sage: DoubleSquare(words.fibonacci_tile(2))
    Double Square Tile
      w0 = 32303010   w4 = 10121232
      w1 = 30323      w5 = 12101
      w2 = 21232303   w6 = 03010121
      w3 = 23212      w7 = 01030
    (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
    (d0, d1, d2, d3)         = (10, 16, 10, 16)
    (n0, n1, n2, n3)         = (0, 0, 0, 0)

::

    sage: DoubleSquare(words.christoffel_tile(4,7))
    Double Square Tile
      w0 = 03                          w4 = 21
      w1 = 0103010103010301010301030   w5 = 2321232321232123232123212
      w2 = 10103010                    w6 = 32321232
      w3 = 1                           w7 = 3
    (|w0|, |w1|, |w2|, |w3|) = (2, 25, 8, 1)
    (d0, d1, d2, d3)         = (26, 10, 26, 10)
    (n0, n1, n2, n3)         = (0, 2, 0, 0)

Double Square tile from the lengths of the `w_i`::

    sage: DoubleSquare((4,7,4,7))
    Double Square Tile
      w0 = 3232      w4 = 1010
      w1 = 1212323   w5 = 3030101
      w2 = 2121      w6 = 0303
      w3 = 0101212   w7 = 2323030
    (|w0|, |w1|, |w2|, |w3|) = (4, 7, 4, 7)
    (d0, d1, d2, d3)         = (14, 8, 14, 8)
    (n0, n1, n2, n3)         = (0, 0, 0, 0)

DoubleSquare tile from the words `(w_0, w_1, w_2, w_3)`::

    sage: DoubleSquare(([3,2], [3], [0,3], [0,1,0,3,0]))
    Double Square Tile
      w0 = 32              w4 = 10
      w1 = 3               w5 = 1
      w2 = 03              w6 = 21
      w3 = 01030           w7 = 23212
    (|w0|, |w1|, |w2|, |w3|) = (2, 1, 2, 5)
    (d0, d1, d2, d3)         = (6, 4, 6, 4)
    (n0, n1, n2, n3)         = (0, 0, 0, 1)

Reduction of a double square tile::

    sage: D = DoubleSquare(words.christoffel_tile(4,7))
    sage: D.reduction()
    ['TRIM_1', 'TRIM_1', 'TRIM_2', 'TRIM_1', 'TRIM_0', 'TRIM_2']
    sage: D.apply_reduction()
    Double Square Tile
      w0 =     w4 =
      w1 = 0   w5 = 2
      w2 =     w6 =
      w3 = 1   w7 = 3
    (|w0|, |w1|, |w2|, |w3|) = (0, 1, 0, 1)
    (d0, d1, d2, d3)         = (2, 0, 2, 0)
    (n0, n1, n2, n3)         = (0, NaN, 0, NaN)

The intermediate steps of the reduction of a double square tile::

    sage: E,op = D.reduce()
    sage: E
    Double Square Tile
      w0 = 03                w4 = 21
      w1 = 010301010301030   w5 = 232123232123212
      w2 = 10103010          w6 = 32321232
      w3 = 1                 w7 = 3
    (|w0|, |w1|, |w2|, |w3|) = (2, 15, 8, 1)
    (d0, d1, d2, d3)         = (16, 10, 16, 10)
    (n0, n1, n2, n3)         = (0, 1, 0, 0)
    sage: op
    'TRIM_1'

::

    sage: D.reduce_ntimes(3)
    Double Square Tile
      w0 = 03      w4 = 21
      w1 = 01030   w5 = 23212
      w2 = 10      w6 = 32
      w3 = 1       w7 = 3
    (|w0|, |w1|, |w2|, |w3|) = (2, 5, 2, 1)
    (d0, d1, d2, d3)         = (6, 4, 6, 4)
    (n0, n1, n2, n3)         = (0, 1, 0, 0)

Plot a double square tile and plot its reduction::

    sage: D = DoubleSquare((34,21,34,21))
    sage: D.plot()                 # long time (1s)
    sage: D.plot_reduction()       # long time (1s)

It is not said clear enough in the articles, but double square reduction
also works for double square tiles that are 8-connected polyominoes::

    sage: D = DoubleSquare((55,34,55,34))
    sage: D.plot()                 # long time (1s)
    sage: D.plot_reduction()       # long time (1s)

"""
#*****************************************************************************
#       Copyright (C) 2008-2013 Sebastien Labbe <slabqc@gmail.com>
#       Copyright (C) 2008-2012 Alexandre Blondin Massé
#
#  Distributed under the terms of the GNU General Public License version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from itertools import izip_longest
from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.real_mpfr import RR
from sage.rings.rational_field import QQ
from sage.modules.free_module_element import vector
from sage.symbolic.constants import NaN, pi
from sage.functions.trig import cos, sin
from sage.combinat.words.words import Words
from sage.combinat.words.word import Word
from sage.combinat.words.morphism import WordMorphism
from sage.plot.point import point
from sage.plot.plot import graphics_array
from sage.misc.misc_c import prod
from sage.misc.functional import numerical_approx
from sage.misc.table import table
from sage.misc.latex import LatexExpr
from sage.misc.lazy_attribute import lazy_attribute
from sage.misc.cachefunc import cached_method

###############################
# Class DoubleSquare
###############################
class DoubleSquare(SageObject):
    r"""
    A double square tile.

    We represent a double square tile by its boundary, that is a finite
    sequence on the alphabet `A=\{0,1,2,3\}` where `0` is a East step, `1` is
    a North step, `2` is a West step and `3` is a South step.

    INPUT:

    - ``data`` - can be one of the following:

      - word - word over over the alphabet A representing the boundary of a
        double square tile

      - tuple - tuple of 4 elements (w0,w1,w2,w3) or 8 elements
        (w0,w1,w2,w3,w4,w5,w6,w7) such that each wi is a sequence over the
        alphabet A. The condition `w_iw_{i+1} = hat(w_{i+4}w_{i+5})` must
        be verified for all `i` modulo 8.

      - tuple - tuple of 4 integers, the lengths of (w0,w1,w2,w3)

    - ``rot180`` - WordMorphism (default: None), involution on the alphabet
      A and representing a rotation of 180 degrees. If None, the morphism
      0->2, 1->3, 2->0, 3->1 is considered.

    - ``steps`` - dict (default: None), mapping letters of A to steps in
      the plane. If None, the corresondance 0->(1,0), 1->(0,1), 2->(-1,0),
      3->(0,-1) is considered.

    EXAMPLES:

    From a double square::

        sage: from sage.combinat.double_square_tile import DoubleSquare
        sage: DoubleSquare(words.fibonacci_tile(1))
        Double Square Tile
          w0 = 32   w4 = 10
          w1 = 3    w5 = 1
          w2 = 03   w6 = 21
          w3 = 0    w7 = 2
        (|w0|, |w1|, |w2|, |w3|) = (2, 1, 2, 1)
        (d0, d1, d2, d3)         = (2, 4, 2, 4)
        (n0, n1, n2, n3)         = (1, 0, 1, 0)
        sage: DoubleSquare(words.fibonacci_tile(2))
        Double Square Tile
          w0 = 32303010   w4 = 10121232
          w1 = 30323      w5 = 12101
          w2 = 21232303   w6 = 03010121
          w3 = 23212      w7 = 01030
        (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
        (d0, d1, d2, d3)         = (10, 16, 10, 16)
        (n0, n1, n2, n3)         = (0, 0, 0, 0)

    ::

        sage: DoubleSquare(words.christoffel_tile(9,7))
        Double Square Tile
          w0 = 03                        w4 = 21
          w1 = 0101030101030101030       w5 = 2323212323212323212
          w2 = 101010301010301010301010  w6 = 323232123232123232123232
          w3 = 1                         w7 = 3
        (|w0|, |w1|, |w2|, |w3|) = (2, 19, 24, 1)
        (d0, d1, d2, d3)         = (20, 26, 20, 26)
        (n0, n1, n2, n3)         = (0, 0, 1, 0)

    From the `w_i`::

        sage: D = DoubleSquare(([],[],[0,1,0,1],[0,1]))
        sage: D.rot180
        WordMorphism: 0->2, 1->3, 2->0, 3->1
        sage: D._steps
        {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}
        sage: D
        Double Square Tile
          w0 =        w4 =
          w1 =        w5 =
          w2 = 0101   w6 = 3232
          w3 = 01     w7 = 32
        (|w0|, |w1|, |w2|, |w3|) = (0, 0, 4, 2)
        (d0, d1, d2, d3)         = (2, 4, 2, 4)
        (n0, n1, n2, n3)         = (0, 0, 2, 0)

    One may also provide strings as long as other arguments are
    consistent::

        sage: steps = {'0':(1,0), '1':(0,1), '2':(-1,0), '3': (0,-1)}
        sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
        sage: DoubleSquare(('','','0101','01','','','3232','32'), rot180, steps)
        Double Square Tile
          w0 =        w4 =
          w1 =        w5 =
          w2 = 0101   w6 = 3232
          w3 = 01     w7 = 32
        (|w0|, |w1|, |w2|, |w3|) = (0, 0, 4, 2)
        (d0, d1, d2, d3)         = (2, 4, 2, 4)
        (n0, n1, n2, n3)         = (0, 0, 2, 0)

    The first four words wi are sufficient::

        sage: steps = {'0':(1,0), '1':(0,1), '2':(-1,0), '3': (0,-1)}
        sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
        sage: DoubleSquare(('','','0101','01'), rot180, steps)
        Double Square Tile
          w0 =        w4 =
          w1 =        w5 =
          w2 = 0101   w6 = 3232
          w3 = 01     w7 = 32
        (|w0|, |w1|, |w2|, |w3|) = (0, 0, 4, 2)
        (d0, d1, d2, d3)         = (2, 4, 2, 4)
        (n0, n1, n2, n3)         = (0, 0, 2, 0)

    """
    def __init__(self, data, rot180=None, steps=None):
        r"""
        Constructor.

        See :DoubleSquare: for documentation.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.christoffel_tile(5,6))
            sage: D
            Double Square Tile
              w0 = 03                           w4 = 21
              w1 = 01030                        w5 = 23212
              w2 = 10103010103010103010103010   w6 = 32321232321232321232321232
              w3 = 1                            w7 = 3
            (|w0|, |w1|, |w2|, |w3|) = (2, 5, 26, 1)
            (d0, d1, d2, d3)         = (6, 28, 6, 28)
            (n0, n1, n2, n3)         = (0, 0, 4, 0)
        """
        if isinstance(data, tuple):
            if not len(data) in (4, 8):
                raise ValueError("len(data) (=%s) must be 4 or 8" % len(data))
            if all(w in Words() for w in data):
                self._w = data
            elif all(isinstance(a, (int, Integer)) for a in data):
                self._w, rot180, steps = double_square_from_four_integers(*data)
            else:
                self._w = map(Word, data)
        elif data in Words():
            self._w, rot180, steps = double_square_from_boundary_word(data)
        else:
            raise TypeError, "Invalid arguments (=%s)" % data

        if rot180 is None:
            rot180 = WordMorphism({0:2,2:0,3:1,1:3})
        self.rot180 = rot180

        if steps is None:
            steps = {}
            steps[0] = vector((1,0))
            steps[1] = vector((0,1))
            steps[2] = vector((-1,0))
            steps[3] = vector((0,-1))
        self._steps = steps

        self.verify_definition()

    @lazy_attribute
    def hat(self):
        r"""
        Return the hat function returning the reversal of a word path.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.w(0)
            Path: 32303010
            sage: D.hat(D.w(0))
            Path: 23212101
        """
        return lambda x:self.rot180(x).reversal()

    def verify_definition(self):
        r"""
        Checks that the input verifies the definition.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.verify_definition()

        ::

            sage: DoubleSquare(([],[0],[0,1,0,1],[0,1]))
            Traceback (most recent call last):
            ...
            AssertionError: wiwi+1 = hat(wi+4,wi+5) is not verified for i=1
        """
        for i in range(4):
            msg = "wiwi+1 = hat(wi+4,wi+5) is not verified for i=%s"%i
            assert self.w(i) * self.w(i+1) == self.hat(self.w(i+4) * self.w((i+5))), msg

    @cached_method
    def w(self, i):
        r"""
        Return the factor w_i

        This corresponds to the new definition of configuration (solution).

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: [D.w(i) for i in range(8)]
            [Path: 32, Path: 3, Path: 03, Path: 0, Path: 10, Path: 1, Path: 21, Path: 2]

        ::

            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: [D.w(i) for i in range(8)]
            [Path: 32303010, Path: 30323, Path: 21232303, Path: 23212, Path: 10121232, Path: 12101, Path: 03010121, Path: 01030]
        """
        i = i % 8
        if 0 <= i < len(self._w):
            return self._w[i]
        elif i == 4:
            return self.hat(self._w[0] * self._w[1])[:len(self._w[0])]
        elif i == 5:
            return self.hat(self._w[0] * self._w[1])[len(self._w[0]):]
        elif i == 6:
            return self.hat(self._w[2] * self._w[3])[:len(self._w[2])]
        elif i == 7:
            return self.hat(self._w[2] * self._w[3])[len(self._w[2]):]
        else:
            raise ValueError, 'i (=%s) must be between 0 and 7.'%i

    @cached_method
    def d(self, i):
        r"""
        Return the integer d_i.

        The value of `d_i` is defined as `d_i=|w_{i-1}|+|w_{i+1}|`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: [D.d(i) for i in range(8)]
            [10, 16, 10, 16, 10, 16, 10, 16]
        """
        return len(self.w((i - 1))) + len(self.w((i + 1)))

    @cached_method
    def n(self, i):
        r"""
        Return the integer n_i.

        The value of `n_i` is defined as the quotient of `|w_i|` by `d_i`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: [D.n(i) for i in range(8)]
            [0, 0, 0, 0, 0, 0, 0, 0]

        ::

            sage: A = D.extend(1).extend(1).extend(1).extend(1)
            sage: [A.n(i) for i in range(8)]
            [0, 4, 0, 0, 0, 4, 0, 0]

        If `d_i=0` then `n_i` is not defined::

            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: B = D.reduce_ntimes(2)
            sage: [B.n(i) for i in range(8)]
            [0, NaN, 0, NaN, 0, NaN, 0, NaN]

        """
        d_i = self.d(i)
        if d_i != 0:
            return len(self.w(i)) // d_i
        else:
            return NaN

    @cached_method
    def u(self, i):
        r"""
        Return the word u_i.

        The word `u_i` is the unique word such that
        `w_i=(u_i*v_i)^{n_i}u_i` where `0\leq |u_i| < d_i`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.u(1)
            Path: 30323

        ::

            sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
            sage: steps = {'0':(1,0), '1':(0,1), '2':(-1,0), '3': (0,-1)}
            sage: D = DoubleSquare(('','03010','','1011'), rot180, steps)
            sage: D.u(0)
            word:
            sage: D.u(1)
            Traceback (most recent call last):
            ...
            ValueError: u_1 is not defined when d_1 == 0
        """
        p = self.hat(self.w((i-3))) * self.w((i-1))
        d_i = self.d(i)
        if d_i:
            len_ui = len(self.w(i)) % d_i
            return p[:len_ui]
        else:
            raise ValueError("u_%s is not defined when d_%s == 0" % (i,i))

    @cached_method
    def v(self, i):
        r"""
        Return the word v_i.

        The word `v_i` is the unique word such that
        `w_i=(u_i*v_i)^{n_i}u_i` where `0\leq |u_i| < d_i`,
        `\hat{w_{i-3}}w_{i-1}=u_iv_i` and `0 < |u_i| \leq d_i`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.v(1)
            Path: 21232303010

        ::

            sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
            sage: steps = {'0':(1,0), '1':(0,1), '2':(-1,0), '3': (0,-1)}
            sage: D = DoubleSquare(('','03010','','1011'), rot180, steps)
            sage: D.v(0)
            word: 030103323
            sage: D.v(1)
            Traceback (most recent call last):
            ...
            ValueError: v_1 is not defined when d_1 == 0

        """
        p = self.hat(self.w((i-3))) * self.w((i-1))
        d_i = self.d(i)
        if d_i:
            len_ui = len(self.w(i)) % d_i
            return p[len_ui:]
        else:
            raise ValueError("v_%s is not defined when d_%s == 0" % (i,i))

    def __eq__(self, other):
        r"""
        Returns True if w0, w1, w2 and w3 are the same.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D == D
            True
            sage: D == D.shift()
            False
            sage: D == D.extend(2).trim(2)
            True
        """
        return (isinstance(other, DoubleSquare) and
               self.w(0) == other.w(0) and
               self.w(1) == other.w(1) and
               self.w(2) == other.w(2) and
               self.w(3) == other.w(3))

    def __hash__(self):
        r"""
        Return the hash value.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: hash(D)                  # random
            398370403             # 32-bit
            3085492728729745145   # 64-bit
            sage: hash(D.swap(0).swap(0))  # random
            398370403             # 32-bit
            3085492728729745145   # 64-bit
        """
        return hash((self.w(0),self.w(1),self.w(2),self.w(3)))

    def _repr_(self):
        r"""
        Return the string representation.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D
            Double Square Tile
              w0 = 32303010   w4 = 10121232
              w1 = 30323      w5 = 12101
              w2 = 21232303   w6 = 03010121
              w3 = 23212      w7 = 01030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
            (d0, d1, d2, d3)         = (10, 16, 10, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
        """
        s = ""
        s += "Double Square Tile\n"
        left  = ["w%s = %s"%(i,self.w(i)) for i in range(4)]
        right = ["w%s = %s"%(i,self.w(i)) for i in range(4,8)]
        s += str(table(columns=[left,right]))
        s += "\n(|w0|, |w1|, |w2|, |w3|) = %s"%(tuple(map(len, (self.w(i) for i in range(4)))),)
        s += "\n(d0, d1, d2, d3)         = %s"%(tuple(self.d(i) for i in range(4)),)
        s += "\n(n0, n1, n2, n3)         = %s"%(tuple(self.n(i) for i in range(4)),)
        return s

    def _latex_(self):
        r"""
        Returns a 8-tuple representing the configuration in Latex

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: latex(D)
            (32303010,30323,21232303,23212,10121232,12101,03010121,01030)

        """
        w = [self.w(i).string_rep() for i in range(8)]
        return "(%s)" % ",".join(w)

    def alphabet(self):
        r"""
        Returns the python set of the letters that occurs in the boundary
        word.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.alphabet()
            set([0, 1, 2, 3])
        """
        return set(self.boundary_word())

    def is_singular(self):
        r"""
        Return whether self is singular.

        A double square is *singular* if there exists `i` such that
        `w_{i-1}` and `w_{i+1}` are empty, equivalently if `d_i=0`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.is_singular()
            False

        ::

            sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
            sage: D = DoubleSquare(('','03010','','1011'), rot180)
            sage: D.is_singular()
            True
        """
        return any(((self.d(0) == 0 and self.d(2) == 0),
                    (self.d(1) == 0 and self.d(3) == 0)))

    def is_flat(self):
        r"""
        Return whether self is flat.

        A double square is *flat* if one of the `w_iw_{i+1}` is empty.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.is_flat()
            False

        ::

            sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
            sage: D = DoubleSquare(('','','0101','01'), rot180)
            sage: D.is_flat()
            True

        """
        return any(len(self.w(i)*self.w(i+1)) == 0 for i in range(8))

    def is_degenerate(self):
        r"""
        Return whether self is degenerate.

        A double square is *degenerate* if one of the `w_i` is empty.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.is_degenerate()
            False

        ::

            sage: rot180 = WordMorphism('0->2,2->0,3->1,1->3')
            sage: D = DoubleSquare(('','0','10','1'), rot180)
            sage: D.is_degenerate()
            True
        """
        return any(len(self.w(i)) == 0 for i in range(8))

    def is_morphic_pentamino(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.is_morphic_pentamino()
            True
        """
        return ((len(self.u(0)) == 0 and len(self.u(2)) == 0) or
                (len(self.u(1)) == 0 and len(self.u(3)) == 0))

    def boundary_word(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.boundary_word()
            Path: 3230301030323212323032321210121232121010...
        """
        alphabet, steps = zip(*self._steps.items())
        from sage.combinat.words.paths import WordPaths
        WP = WordPaths(alphabet, steps)
        boundary = prod(self.w(i) for i in range(8))
        return WP(boundary)

    def turning_number(self):
        r"""
        Return the turning number of self.

        INPUT:

        - ``self`` - double square defined on the alphabet of integers
          `\{0,1,2,3\}`

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.turning_number()
            1
            sage: D.reverse().turning_number()
            -1

        Turning number of a degenerate double square::

            sage: D = DoubleSquare(([],[0],[1,0],[1]))
            sage: D.turning_number()
            1

        Turning number of a singular double square::

            sage: D = DoubleSquare(([],[0,3,0,1,0],[],[1,0,1,1]))
            sage: D.turning_number()
            1

        Turning number of a flat double square::

            sage: D = DoubleSquare(([],[],[0,1,0,1],[0,1]))
            sage: D.turning_number()
            0

        """
        boundary = self.boundary_word()
        boundary = boundary + boundary[:1]
        boundary = boundary.to_integer_word()
        turns = boundary.finite_differences(mod=4)
        ev = turns.evaluation_dict()
        return QQ(((ev[1] - ev[3]), 4))

    def shift(self):
        r"""
        Apply `SHIFT` on self.

        This replaces `w_i` by `w_{i+1}`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D
            Double Square Tile
              w0 = 32303010   w4 = 10121232
              w1 = 30323      w5 = 12101
              w2 = 21232303   w6 = 03010121
              w3 = 23212      w7 = 01030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
            (d0, d1, d2, d3)         = (10, 16, 10, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.shift()
            Double Square Tile
              w0 = 30323      w4 = 12101
              w1 = 21232303   w5 = 03010121
              w2 = 23212      w6 = 01030
              w3 = 10121232   w7 = 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)
            (d0, d1, d2, d3)         = (16, 10, 16, 10)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.shift().shift().shift().shift().shift().shift().shift().shift() == D
            True
            sage: D.shift().shift().shift().shift() == D
            False
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = tuple(self.w(i) for i in range(8))
        return DoubleSquare((w1,w2,w3,w4,w5,w6,w7,w0), self.rot180, self._steps)

    def reverse(self):
        r"""
        Apply `REVERSE` on self.

        This reverses the words `w_i`.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D
            Double Square Tile
              w0 = 32303010   w4 = 10121232
              w1 = 30323      w5 = 12101
              w2 = 21232303   w6 = 03010121
              w3 = 23212      w7 = 01030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
            (d0, d1, d2, d3)         = (10, 16, 10, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.reverse()
            Double Square Tile
              w0 = 21232      w4 = 03010
              w1 = 30323212   w5 = 12101030
              w2 = 32303      w6 = 10121
              w3 = 01030323   w7 = 23212101
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)
            (d0, d1, d2, d3)         = (16, 10, 16, 10)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.reverse().reverse() == D
            True
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = map(self.w, range(8))
        w = (self.hat(w7),self.hat(w6),self.hat(w5),self.hat(w4),
             self.hat(w3),self.hat(w2),self.hat(w1),self.hat(w0))
        return DoubleSquare(w, self.rot180, self._steps)

    def extend(self, i):
        r"""
        Apply `EXTEND_i` on self.

        This adds a period of length `d_i` to `w_i` and `w_{i+4}`.

        INPUT:

        - ``i`` - integer

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.extend(3)
            Double Square Tile
              w0 = 32303010                w4 = 10121232
              w1 = 30323                   w5 = 12101
              w2 = 21232303                w6 = 03010121
              w3 = 232121012123230323212   w7 = 010303230301012101030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 21)
            (d0, d1, d2, d3)         = (26, 16, 26, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 1)
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = [self.w(j) for j in range(i % 8, 8) + range(0, i % 8)]
        w = (w0*w1*self.hat(w3),w1,w2,w3,w4*w5*self.hat(w7),w5,w6,w7)
        w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
        return DoubleSquare(w, self.rot180, self._steps)

    def trim(self, i):
        r"""
        Apply `TRIM_i` on self.

        This removes a period of length `d_i` to `w_i` and `w_{i+4}`.

        INPUT:

        - ``i`` - integer

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare((3,6,3,2))
            sage: D.trim(1)
            Double Square Tile
              w0 = 212   w4 = 030
              w1 =       w5 =
              w2 = 303   w6 = 121
              w3 = 03    w7 = 21
            (|w0|, |w1|, |w2|, |w3|) = (3, 0, 3, 2)
            (d0, d1, d2, d3)         = (2, 6, 2, 6)
            (n0, n1, n2, n3)         = (1, 0, 1, 0)

        ::

            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.extend(3).trim(3)
            Double Square Tile
              w0 = 32303010   w4 = 10121232
              w1 = 30323      w5 = 12101
              w2 = 21232303   w6 = 03010121
              w3 = 23212      w7 = 01030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
            (d0, d1, d2, d3)         = (10, 16, 10, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
        """
        d = self.d(i)
        if len(self.w(i)) >= d:
            (w0,w1,w2,w3,w4,w5,w6,w7) = [self.w(j) for j in range(i % 8, 8) + range(0, i % 8)]
            w = (w0[d:],w1,w2,w3,w4[d:],w5,w6,w7)
            w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
            return DoubleSquare(w, self.rot180, self._steps)
        else:
            raise ValueError, 'trim cannot be applied on index %s of the following configuration\n%s'%(i,self)
    def swap(self, i):
        r"""
        Apply `SWAP_i` on self.

        This replaces `w_j` by `\hat{w_{j+4}}` for each `j=i,i+2,i+4,i+6`
        and `w_j=(u_j*v_j)^{n_j}u_j` by `(v_j*u_j)^{n_j}v_j` for each
        `j=i+1,i+3,i+5,i+7`. This is an involution if the `u_j` are non
        empty.

        INPUT:

        - ``i`` - integer

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D
            Double Square Tile
              w0 = 32303010   w4 = 10121232
              w1 = 30323      w5 = 12101
              w2 = 21232303   w6 = 03010121
              w3 = 23212      w7 = 01030
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)
            (d0, d1, d2, d3)         = (10, 16, 10, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.swap(1)
            Double Square Tile
              w0 = 30      w4 = 12
              w1 = 32303   w5 = 10121
              w2 = 23      w6 = 01
              w3 = 21232   w7 = 03010
            (|w0|, |w1|, |w2|, |w3|) = (2, 5, 2, 5)
            (d0, d1, d2, d3)         = (10, 4, 10, 4)
            (n0, n1, n2, n3)         = (0, 1, 0, 1)
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = [self.w(j) for j in range(i % 8, 8) + range(0, i % 8)]
        indexes = range(i % 8 + 1, 8, 2) + range((1 + i) % 2, i % 8, 2)
        (n1,n3,n5,n7) = [self.n(j) for j in indexes]
        (u1,u3,u5,u7) = [self.u(j) for j in indexes]
        (v1,v3,v5,v7) = [self.v(j) for j in indexes]
        # assertions for swap to be an involution:
        #assert len(u1) > 0, "one must have |u_%s| > 0 for swap_%s" % (i+1, i)
        #assert len(u3) > 0, "one must have |u_%s| > 0 for swap_%s" % (i+3, i)
        w = (self.hat(w4),(v1*u1)**n1*v1,self.hat(w6),(v3*u3)**n3*v3,
             self.hat(w0),(v5*u5)**n5*v5,self.hat(w2),(v7*u7)**n7*v7)
        w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
        return DoubleSquare(w, self.rot180, self._steps)

    def apply(self, L):
        r"""
        Return the double square obtained after the application of a list
        of operations.

        INPUT:

        - ``L`` - list, list of strings

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.apply(['SWAP_0', 'EXTEND_3', 'TRIM_3'])
            Double Square Tile
              w0 = 01030323      w4 = 23212101
              w1 = 21232303010   w5 = 03010121232
              w2 = 30323212      w6 = 12101030
              w3 = 10121232303   w7 = 32303010121
            (|w0|, |w1|, |w2|, |w3|) = (8, 11, 8, 11)
            (d0, d1, d2, d3)         = (22, 16, 22, 16)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)

        ::

            sage: D.apply(D.reduction())
            Double Square Tile
              w0 =     w4 =
              w1 = 3   w5 = 1
              w2 =     w6 =
              w3 = 2   w7 = 0
            (|w0|, |w1|, |w2|, |w3|) = (0, 1, 0, 1)
            (d0, d1, d2, d3)         = (2, 0, 2, 0)
            (n0, n1, n2, n3)         = (0, NaN, 0, NaN)
        """
        ds = self
        for op_i in L:
            op,i = op_i.split('_')
            i = int(i)
            if op == "TRIM":
                ds = ds.trim(i)
            elif op == "SWAP":
                ds = ds.swap(i)
            elif op == "EXTEND":
                ds = ds.extend(i)
            else:
                raise ValueError("Unknown operation (-%s)" % op)
        return ds

    def reduce(self):
        r"""
        Reduces self by the application of TRIM or otherwise SWAP.

        INPUT:

        - ``self`` - non singular double square tile on the alphabet
          `{0,1,2,3}` such that its turning number is +1 or -1.

        OUTPUT:

        - DoubleSquare - the reduced double square
        - string - the operation which was performed

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare((34,21,34,21))
            sage: E,op = D.reduce()
            sage: E
            Double Square Tile
              w0 = 32303010                w4 = 10121232
              w1 = 303232123230301030323   w5 = 121010301012123212101
              w2 = 21232303                w6 = 03010121
              w3 = 232121012123230323212   w7 = 010303230301012101030
            (|w0|, |w1|, |w2|, |w3|) = (8, 21, 8, 21)
            (d0, d1, d2, d3)         = (42, 16, 42, 16)
            (n0, n1, n2, n3)         = (0, 1, 0, 1)
            sage: op
            'SWAP_1'

        ::

            sage: D = DoubleSquare((1,2,2,1))
            sage: D
            Double Square Tile
              w0 = 1    w4 = 1
              w1 = 23   w5 = 03
              w2 = 12   w6 = 10
              w3 = 3    w7 = 3
            (|w0|, |w1|, |w2|, |w3|) = (1, 2, 2, 1)
            (d0, d1, d2, d3)         = (3, 3, 3, 3)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.reduce()
            Traceback (most recent call last):
            ...
            ValueError: not reducible, because self is nondegenerate and
            d_0 == d_1 == 3. Also, the turning number (=0) must be +1 or -1
            for the reduction to apply.

        TESTS::

            sage: D = DoubleSquare((5,4,3,4))
            sage: D
            Double Square Tile
              w0 = 90128   w4 = 40123
              w1 = 7659    w5 = 7654
              w2 = 012     w6 = 012
              w3 = 3765    w7 = 8765
            (|w0|, |w1|, |w2|, |w3|) = (5, 4, 3, 4)
            (d0, d1, d2, d3)         = (8, 8, 8, 8)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.reduce()
            Traceback (most recent call last):
            ...
            ValueError: not reducible, because self is nondegenerate and
            d_0 == d_1 == 8. Also, the turning number (=-1) must be +1 or
            -1 for the reduction to apply.
        """
        # We verify if we have a singular DS
        if self.is_singular():
            raise ValueError('not reducible, this is a singular DS-factorization')

        # Now we check if TRIM may be applied
        for i in range(8):
            if len(self.w(i)) >= self.d(i):
                return (self.trim(i), 'TRIM_%s' % i)

        # Now we check if SWAP may be applied
        if all(0 < len(self.w(i)) < self.d(i) for i in range(4)):
            # We try with SWAP
            for i in range(2):
                if len(self.v(i+1)) + len(self.v(i+3)) < len(self.u(i+1)) + len(self.u(i+3)):
                    return (self.swap(i), 'SWAP_%s' % i)

        # Now we check if SWAP may be applied
        for i in range(4):
            if all(( 0 == len(self.w(i)),
                    0 < len(self.w(i+1)) < self.d(i+1),
                    0 < len(self.w(i+2)) < self.d(i+2),
                    0 < len(self.w(i+3)) < self.d(i+3))):
                return (self.swap(i), 'SWAP_%s' % i)

        # Limit cases
        if not self.is_degenerate() and self.d(0) == self.d(1):
            if len(self.alphabet()) == 4:
                assert self.turning_number() not in [-1,1], "Lemma 16 of [BGL2012]"
            raise ValueError('not reducible, because self is nondegenerate and ' +
                    'd_0 == d_1 == %s. ' % self.d(0)+
                    "Also, the turning number (=%s) must be " % self.turning_number()+
                    "+1 or -1 for the reduction to apply.")

        # We should not get here
        raise ValueError("This double square seems to be not reducible,"
                         "it might be a counter example to a theorem !!!\n%s" % self)

    def reduce_ntimes(self, iteration=1):
        r"""
        Reduces the double square self until it is singular.

        INPUT:

        - ``iteration`` - integer (default: ``1``), number of iterations to
          perform

        OUTPUT:

            DoubleSquare

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare((34,21,34,21))
            sage: D.reduce_ntimes(10)
            Double Square Tile
              w0 =     w4 =
              w1 = 3   w5 = 1
              w2 =     w6 =
              w3 = 2   w7 = 0
            (|w0|, |w1|, |w2|, |w3|) = (0, 1, 0, 1)
            (d0, d1, d2, d3)         = (2, 0, 2, 0)
            (n0, n1, n2, n3)         = (0, NaN, 0, NaN)

        """
        ds = self
        for _ in range(iteration):
            if ds.is_singular():
                break
            ds, op = ds.reduce()
        return ds

    def reduction(self):
        r"""
        Return the list of operations to reduce self to a singular double
        square.

        OUTPUT:

        - list of strings

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare((34,21,34,21))
            sage: D.reduction()
            ['SWAP_1', 'TRIM_1', 'TRIM_3', 'SWAP_1', 'TRIM_1', 'TRIM_3', 'TRIM_0', 'TRIM_2']

        ::

            sage: D = DoubleSquare(words.christoffel_tile(9,7))
            sage: D.reduction()
            ['TRIM_2', 'TRIM_1', 'TRIM_1', 'TRIM_1', 'TRIM_0', 'TRIM_2', 'TRIM_2']
        """
        L = []
        ds = self
        while not ds.is_singular():
            ds, op = ds.reduce()
            L.append(op)
        return L

    def apply_reduction(self):
        r"""
        Apply the reduction algorithm on self.

        This is equivalent to ``self.apply(self.reduction())``.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.christoffel_tile(9,7))
            sage: D.apply_reduction()
            Double Square Tile
              w0 =     w4 =
              w1 = 0   w5 = 2
              w2 =     w6 =
              w3 = 1   w7 = 3
            (|w0|, |w1|, |w2|, |w3|) = (0, 1, 0, 1)
            (d0, d1, d2, d3)         = (2, 0, 2, 0)
            (n0, n1, n2, n3)         = (0, NaN, 0, NaN)

        ::

            sage: D = DoubleSquare((5,7,4,13))
            sage: D.apply_reduction()
            Double Square Tile
              w0 =     w4 =
              w1 =     w5 =
              w2 = 1   w6 = 0
              w3 =     w7 =
            (|w0|, |w1|, |w2|, |w3|) = (0, 0, 1, 0)
            (d0, d1, d2, d3)         = (0, 1, 0, 1)
            (n0, n1, n2, n3)         = (NaN, 0, NaN, 0)

        ::

            sage: D = DoubleSquare((5,2,4,13))
            sage: D.reduce_ntimes(3)
            Double Square Tile
              w0 = 0    w4 = 0
              w1 = 12   w5 = 32
              w2 = 01   w6 = 03
              w3 = 2    w7 = 2
            (|w0|, |w1|, |w2|, |w3|) = (1, 2, 2, 1)
            (d0, d1, d2, d3)         = (3, 3, 3, 3)
            (n0, n1, n2, n3)         = (0, 0, 0, 0)
            sage: D.apply_reduction()
            Traceback (most recent call last):
            ...
            ValueError: not reducible, because self is nondegenerate and
            d_0 == d_1 == 3. Also, the turning number (=0) must be +1 or -1
            for the reduction to apply.
        """
        ds = self
        while not ds.is_singular():
            ds, op = ds.reduce()
        return ds

    def factorization_points(self):
        r"""
        Returns the eight factorization points of this configuration

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.factorization_points()
            [0, 2, 3, 5, 6, 8, 9, 11]
        """
        lengths_of_w = map(len, self._w)
        return [sum(lengths_of_w[:i]) for i in range(8)]

    def apply_morphism(self, m):
        r"""
        INPUT:

        - ``m`` - a WordMorphism

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: m = WordMorphism({0:[0],1:[1,0,1],2:[2],3:[3,2,3]})
            sage: D.apply_morphism(m)
            Double Square Tile
              w0 = 3232   w4 = 1010
              w1 = 323    w5 = 101
              w2 = 0323   w6 = 2101
              w3 = 0      w7 = 2
            (|w0|, |w1|, |w2|, |w3|) = (4, 3, 4, 1)
            (d0, d1, d2, d3)         = (4, 8, 4, 8)
            (n0, n1, n2, n3)         = (1, 0, 1, 0)
        """
        w = tuple(m(self.w(j)) for j in range(8))
        return DoubleSquare(w, self.rot180, self._steps)

    def width(self):
        r"""
        Returns the width of this polyomino, i.e. the difference
        between its rightmost and leftmost coordinates

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.width()
            9

        ::

            sage: D = DoubleSquare((34,21,34,21))
            sage: D.width()
            23
        """
        points = list(self.boundary_word().points())
        return max(map(lambda p:p[0], points)) - min(map(lambda p:p[0], points))

    def height(self):
        r"""
        Returns the width of this polyomino, i.e. the difference
        between its uppermost and lowermost coordinates

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.height()
            9

        ::

            sage: D = DoubleSquare((34,21,34,21))
            sage: D.height()
            23
        """
        points = list(self.boundary_word().points())
        return max(map(lambda p:p[1], points)) - min(map(lambda p:p[1], points))

    def plot(self, pathoptions=dict(rgbcolor='black',thickness=3),
         fill=True, filloptions=dict(rgbcolor='black',alpha=0.2),
         startpoint=True, startoptions=dict(rgbcolor='black',pointsize=100),
         endarrow=True, arrowoptions=dict(rgbcolor='black',arrowsize=5,width=3),
         gridlines=False, gridoptions=dict(),
         axes=False):
        r"""
        Returns a 2d Graphics illustrating the double square tile associated to
        this configuration including the factorizations points.

        INPUT:

        - ``pathoptions`` - (dict,
          default:dict(rgbcolor='red',thickness=3)), options for the
          path drawing

        - ``fill`` - (boolean, default: True), if fill is True and if
          the path is closed, the inside is colored

        - ``filloptions`` - (dict,
          default:dict(rgbcolor='red',alpha=0.2)), ptions for the
          inside filling

        - ``startpoint`` - (boolean, default: True), draw the start point?

        - ``startoptions`` - (dict,
          default:dict(rgbcolor='red',pointsize=100)) options for the
          start point drawing

        - ``endarrow`` - (boolean, default: True), draw an arrow end at the end?

        - ``arrowoptions`` - (dict,
          default:dict(rgbcolor='red',arrowsize=20, width=3)) options
          for the end point arrow

        - ``gridlines``- (boolean, default: False), show gridlines?

        - ``gridoptions`` - (dict, default: {}), options for the gridlines

        - ``axes`` - (boolean, default: False), options for the axes


        EXAMPLES:

        The cross of area 5 together with its double square factorization points::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.plot()              # long time (1s)
        """
        path = self.boundary_word()
        points = list(path.points())
        points = [map(RR, x) for x in points]
        G = path.plot(pathoptions, fill, filloptions, startpoint, startoptions, endarrow, arrowoptions, gridlines, gridoptions)
        i = 0
        for p in self.factorization_points():
            if i % 2 == 0: G += point(points[p], pointsize=startoptions['pointsize'], rgbcolor="red")
            else: G += point(points[p], pointsize=startoptions['pointsize'], rgbcolor="blue")
            i += 1
        G.axes(axes)
        return G

    def plot_reduction(self, ncols=3, options={}):
        r"""
        Return a graphics array of the reduction.

        INPUT:

        - ``ncols`` - integer (default: ``3``), number of columns
        - ``options`` - dict (default: ``{}``), options given to the plot
          method of each double square

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.plot_reduction()          # long time (1s)

        Using the color options::

            sage: p = dict(rgbcolor='red', thickness=1)
            sage: q = dict(rgbcolor='blue', alpha=1)
            sage: options = dict(endarrow=False,startpoint=False,pathoptions=p,filloptions=q)
            sage: D.plot_reduction(options=options)      # long time (1s)
        """
        ds = self
        L = [ds]
        while not ds.is_singular():
            ds, op = ds.reduce()
            L.append(ds)
        nrows = (len(L)-1) / ncols + 1
        graphic_tiles = []
        for tile in L:
            graphic_tiles.append(tile.plot(**options))
        G = graphics_array(graphic_tiles, nrows, ncols)
        return G

    def latex_8_tuple(self):
        r"""
        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.latex_8_tuple()
            ('{\\bf 32303010}', '{\\bf 30323}', '{\\bf 21232303}', '{\\bf 23212}',
             '{\\bf 10121232}', '{\\bf 12101}', '{\\bf 03010121}', '{\\bf 01030}')
        """
        W = tuple("{\\bf %s}" % self.w(i) if self.w(i) else "\\varepsilon" for i in range(8))
        return W

    def latex_array(self):
        r"""
        Return a LaTeX array of self.

        This code was used to create Table 1 in [BGL2012]_.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: print D.latex_array()
            \begin{array}{lllllll}
            i & w_i & u_i & v_i & |w_i| & d_i & n_i
            \\
            \hline
            0 & {\bf 32} & {\bf } & {\bf 32} & 2 & 2 & 1\\
            1 & {\bf 3} & {\bf 3} & {\bf 032} & 1 & 4 & 0\\
            2 & {\bf 03} & {\bf } & {\bf 03} & 2 & 2 & 1\\
            3 & {\bf 0} & {\bf 0} & {\bf 103} & 1 & 4 & 0\\
            4 & {\bf 10} & {\bf } & {\bf 10} & 2 & 2 & 1\\
            5 & {\bf 1} & {\bf 1} & {\bf 210} & 1 & 4 & 0\\
            6 & {\bf 21} & {\bf } & {\bf 21} & 2 & 2 & 1\\
            7 & {\bf 2} & {\bf 2} & {\bf 321} & 1 & 4 & 0\\
            \hline
            \end{array}
        """
        s = "\\begin{array}{lllllll}\n"
        s += "i & w_i & u_i & v_i & |w_i| & d_i & n_i \n"
        s += "\\\\ \n\\hline\n"
        for i in range(8):
            s += "%s & " % i
            s += "{\\bf %s} & " % self.w(i)
            s += "{\\bf %s} & " % self.u(i)
            s += "{\\bf %s} & " % self.v(i)
            s += "%s & " % len(self.w(i))
            s += "%s & " % self.d(i)
            s += "%s" % self.n(i)
            s += "\\\\ \n"
        s += '\\hline\n\\end{array}\n'
        return s

    def latex_table(self):
        r"""
        Returns a Latex expression of a table containing the parameters of
        this double square.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(2))
            sage: D.latex_table()
            \begin{tabular}{|c|}
            \hline
            \\
            \begin{tikzpicture}
            [first/.style={circle,draw=black,fill=gray, inner sep=0pt, minimum size=3pt},
            second/.style={rectangle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]
            ...
            \end{tikzpicture} \\[1ex]
            \hline\\
            $(w_0,w_1,w_2,w_3) = (8,5,8,5)$ \\
            $u_0 = 32303010$\quad $u_1 = 30323$\\$u_2 = 21232303$\quad $u_3 = 23212$\\
            $v_0 = 30$\quad $v_1 = 21232303010$\\$v_2 = 23$\quad $v_3 = 10121232303$\\
            $(n_0,n_1,n_2,n_3) = (0,0,0,0)$ \\
            Turning number = -1\\
            Self-avoiding = True\\
            \hline
            \end{tabular}
        """
        remove_coma = lambda s:s.translate(None, ',')
        if_empty = lambda s:'\\varepsilon' if len(s) == 0 else s
        u = [if_empty(remove_coma(self.u(i).string_rep())) for i in range(4)]
        v = [if_empty(remove_coma(self.v(i).string_rep())) for i in range(4)]
        #uv = [self.u(i) for i in range(4)] + [self.v(i) for i in range(4)]

        s = '\\begin{tabular}{|c|}\n\\hline\n\\\\\n'
        s += '\\begin{tikzpicture}\n'
        s += '  [first/.style={circle,draw=black,fill=gray, inner sep=0pt, minimum size=3pt},\n'
        s += '   second/.style={rectangle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]\n'
        s += self.tikz_trajectory(step=5.0/max(self.width(), self.height()))
        s += '\n\\end{tikzpicture} \\\\[1ex] \n\\hline\\\\\n'
        #for i in range(4):
        #    s += '$w_{%s} = %s$\\\\\n'%(i, remove_coma(self.w(i).string_rep()))
        s += '$(w_0,w_1,w_2,w_3) = (%s,%s,%s,%s)$ \\\\\n'%tuple(len(self.w(i)) for i in range(4))
        s += '$u_0 = %s$\quad $u_1 = %s$\\\\$u_2 = %s$\quad $u_3 = %s$\\\\\n'%tuple(u[i] for i in range(4))
        s += '$v_0 = %s$\quad $v_1 = %s$\\\\$v_2 = %s$\quad $v_3 = %s$\\\\\n'%tuple(v[i] for i in range(4))
        s += '$(n_0,n_1,n_2,n_3) = (%s,%s,%s,%s)$ \\\\\n'%tuple(self.n(i) for i in range(4))
        s += 'Turning number = %s\\\\\n'%self.turning_number()
        s += 'Self-avoiding = %s\\\\\n'%self.boundary_word().is_simple()
        #s += 'Any $u_i$ or $v_i$ self-crossing ? %s\\\\\n'%any(map(lambda w:not w.is_simple(), uv))
        s += '\\hline\n\\end{tabular}\n'
        return LatexExpr(s)

    def tikz_trajectory(self, step=1, arrow='->'):
        r"""
        Returns a tikz string describing the double square induced by
        this configuration together with its factorization points

        The factorization points respectively get the tikz attribute 'first'
        and 'second' so that when including it in a tikzpicture environment,
        it is possible to modify the way those points appear.

        INPUT:

        - ``step`` - integer (default: ``1``)
        - ``arrow`` - string (default: ``->``), tikz arrow shape

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.tikz_trajectory()
            \filldraw[->, very thick, draw=black, fill=black!20] (0.000, 0.000)
            -- (0.000, -1.00) -- (-1.00, -1.00) -- (-1.00, -2.00) -- (0.000, -2.00) --
            (0.000, -3.00) -- (1.00, -3.00) -- (1.00, -2.00) -- (2.00, -2.00) -- (2.00,
            -1.00) -- (1.00, -1.00) -- (1.00, 0.000) -- (0.000, 0.000); \node[first] at
            (0.0000, 0.0000) {};
            \node[first] at (-1.000, -2.000) {};
            \node[first] at (1.000, -3.000) {};
            \node[first] at (2.000, -1.000) {};
            \node[second] at (-1.000, -1.000) {};
            \node[second] at (0.0000, -3.000) {};
            \node[second] at (2.000, -2.000) {};
            \node[second] at (1.000, 0.0000) {};
        """
        f = lambda x: numerical_approx(x,digits=3)
        step = numerical_approx(step, digits=4)
        points = map(lambda (x,y):(x*step,y*step), list(self.boundary_word().points()))
        l = [str(tuple(map(f, pt))) for pt in points]
        s = '\\filldraw[%s, very thick, draw=black, fill=black!20] ' %arrow
        s += ' -- '.join(l) + ';'
        [a1, b1, a2, b2, a3, b3, a4, b4] = self.factorization_points()
        s += '\n  \\node[first] at %s {};'% (points[a1],)
        s += '\n  \\node[first] at %s {};'% (points[a2],)
        s += '\n  \\node[first] at %s {};'% (points[a3],)
        s += '\n  \\node[first] at %s {};'% (points[a4],)
        s += '\n  \\node[second] at %s {};'% (points[b1],)
        s += '\n  \\node[second] at %s {};'% (points[b2],)
        s += '\n  \\node[second] at %s {};'% (points[b3],)
        s += '\n  \\node[second] at %s {};'% (points[b4],)
        return LatexExpr(s)

    def tikz_boxed(self, scale=1, boxsize=10):
        r"""
        Return a tikzpicture of self included in a box.

        INPUT:

        - ``scale`` - number (default: ``1``), tikz scale
        - ``boxsize`` - integer (default: ``10``), size of the box. If the
          width and height of the double square is less than the boxsize,
          then unit step are of size ``1`` and the `(w_i)` 8-tuple is added
          below the figure. Otherwise, if the width or height is larger
          than the boxsize, then the unit step are made smaller to fit the
          box and the `(w_i)` 8-tuple is not shown.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: D = DoubleSquare(words.fibonacci_tile(1))
            sage: D.tikz_boxed()
            \begin{tabular}{c}
            \begin{tikzpicture}
            [scale=1]
            \filldraw[-to, very thick, draw=black, fill=black!20] (0.000,
            0.000) -- (0.000, -1.00) -- (-1.00, -1.00) -- (-1.00, -2.00)
            -- (0.000, -2.00) -- (0.000, -3.00) -- (1.00, -3.00) -- (1.00,
            -2.00) -- (2.00, -2.00) -- (2.00, -1.00) -- (1.00, -1.00) --
            (1.00, 0.000) -- (0.000, 0.000);
            \node[first] at (0.0000, 0.0000) {};
            \node[first] at (-1.000, -2.000) {};
            \node[first] at (1.000, -3.000) {};
            \node[first] at (2.000, -1.000) {};
            \node[second] at (-1.000, -1.000) {};
            \node[second] at (0.0000, -3.000) {};
            \node[second] at (2.000, -2.000) {};
            \node[second] at (1.000, 0.0000) {};
            \end{tikzpicture}
            \\
            $({\bf 32},{\bf 3},{\bf 03},{\bf 0},$ \\
            $\phantom{((}{\bf 10},{\bf 1},{\bf 21},{\bf 2})$ \\
            \end{tabular}

        Smaller boxsize::

            sage: D.tikz_boxed(boxsize=1.5)
            \begin{tikzpicture}
            [scale=1]
            \filldraw[-to, very thick, draw=black, fill=black!20] (0.000, 0.000) --
            (0.000, -0.500) -- (-0.500, -0.500) -- (-0.500, -1.00) -- (0.000, -1.00) --
            (0.000, -1.50) -- (0.500, -1.50) -- (0.500, -1.00) -- (1.00, -1.00) --
            (1.00, -0.500) -- (0.500, -0.500) -- (0.500, 0.000) -- (0.000, 0.000);
            \node[first] at (0.0000, 0.0000) {};
            \node[first] at (-0.5000, -1.000) {};
            \node[first] at (0.5000, -1.500) {};
            \node[first] at (1.000, -0.5000) {};
            \node[second] at (-0.5000, -0.5000) {};
            \node[second] at (0.0000, -1.500) {};
            \node[second] at (1.000, -1.000) {};
            \node[second] at (0.5000, 0.0000) {};
            \end{tikzpicture}

        """
        (w, h) = (self.width(), self.height())
        M = max(w,h)
        step = boxsize / M if M > boxsize else 1
        s = ""
        if M <= boxsize:
            s += '\\begin{tabular}{c}\n'
        s += '\\begin{tikzpicture}\n'
        s += '[scale=%s]\n' % scale
        s += self.tikz_trajectory(step=step, arrow='-to')
        s += '\n\\end{tikzpicture}\n'
        if M <= boxsize:
            s += "\\\\\n"
            W = self.latex_8_tuple()
            s += '$(%s,%s,%s,%s,$ \\\\\n' % W[:4]
            s += '$\phantom{((}%s,%s,%s,%s)$ \\\\\n' % W[4:]
            s += '\\end{tabular}\n'
        return LatexExpr(s)

    def tikz_reduction(self, scale=1, ncols=3, gridstep=5, labels=True, newcommand=True):
        r"""
        INPUT:

        - ``scale`` - number
        - ``ncols`` - integer, number of columns displaying the
          reduction
        - ``gridstep`` - number (default: ``5``), the gridstep for the
          snake node positions
        - ``labels`` - arrow labels (default:``True``). It may take the
          following values:

            - ``True`` - prints TRIM, SWAP, etc.
            - ``'T'`` - prints T_i, etc.
            - ``False`` - print nothing

        - ``newcommand`` - bool (default: ``True``), whether newcommand
          which defines ``\SWAP``, ``\TRIM``, etc.

        EXAMPLES::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: fibo2 = words.fibonacci_tile(2)
            sage: cfibo2 = DoubleSquare(fibo2)
            sage: s = cfibo2.tikz_reduction(scale=0.5,ncols=4,labels=True)
            sage: s
            \newcommand{\TRIM}{\textsc{trim}}
            \newcommand{\EXTEND}{\textsc{extend}}
            \newcommand{\SWAP}{\textsc{swap}}
            \newcommand{\SHIFT}{\textsc{shift}}
            \newcommand{\REVERSE}{\textsc{reverse}}
            \begin{tikzpicture}
            [first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},
            second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]
            \node (q0) at (0, 0) {
            \begin{tikzpicture}
            [scale=0.500000000000000]
            ...
            \end{tikzpicture}
            \\
            $(\varepsilon,{\bf 3},\varepsilon,{\bf 2},$ \\
            $\phantom{((}\varepsilon,{\bf 1},\varepsilon,{\bf 0})$ \\
            \end{tabular}
            };
            \path[->] (q0) edge node[midway, rectangle, fill=white, rotate=90] {$\SWAP_1$} (q1);
            \path[->] (q1) edge node[midway, rectangle, fill=white, rotate=90] {$\TRIM_1$} (q2);
            \path[->] (q2) edge node[midway, rectangle, fill=white, rotate=90] {$\TRIM_3$} (q3);
            \path[->] (q3) edge node[midway, rectangle, fill=white] {$\TRIM_0$} (q4);
            \path[->] (q4) edge node[midway, rectangle, fill=white, rotate=90] {$\TRIM_2$} (q5);
            \end{tikzpicture}

        ::

            sage: S = WordMorphism({0:[0,0],1:[1,0,1],2:[2,2],3:[3,2,3]}, codomain=fibo2.parent())
            sage: cSfibo2 = cfibo2.apply_morphism(S)
            sage: s = cSfibo2.tikz_reduction(scale=0.15,ncols=4,labels='T')
            sage: s
            \newcommand{\TRIM}{\textsc{trim}}
            \newcommand{\EXTEND}{\textsc{extend}}
            \newcommand{\SWAP}{\textsc{swap}}
            \newcommand{\SHIFT}{\textsc{shift}}
            \newcommand{\REVERSE}{\textsc{reverse}}
            \begin{tikzpicture}
            [first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},
            second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]
            \node (q0) at (0, 0) {
            \begin{tikzpicture}
            [scale=0.150000000000000]
            ...
            \end{tikzpicture}
            \\
            $(\varepsilon,{\bf 323},\varepsilon,{\bf 22},$ \\
            $\phantom{((}\varepsilon,{\bf 101},\varepsilon,{\bf 00})$ \\
            \end{tabular}
            };
            \path[thick, ->] (q0) edge node[midway, above] {$T_1$} (q1);
            \path[thick, ->] (q1) edge node[midway, above] {$T_2$} (q2);
            \path[thick, ->] (q2) edge node[midway, above] {$T_3$} (q3);
            \path[thick, ->] (q3) edge node[midway, above] {$T_4$} (q4);
            \path[thick, ->] (q4) edge node[midway, above] {$T_5$} (q5);
            \end{tikzpicture}

        """
        s = ''
        if newcommand:
            s += "\\newcommand{\\TRIM}{\\textsc{trim}}\n"
            s += "\\newcommand{\\EXTEND}{\\textsc{extend}}\n"
            s += "\\newcommand{\\SWAP}{\\textsc{swap}}\n"
            s += "\\newcommand{\\SHIFT}{\\textsc{shift}}\n"
            s += "\\newcommand{\\REVERSE}{\\textsc{reverse}}\n"
        s += "\\begin{tikzpicture}\n"
        s += "[first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},\n"
        s += "second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]\n"
        D = self
        i = 0
        fns = []
        while True:
            s += '\\node (q%s) '%i
            x,y = snake(i, ncols)
            s += 'at %s '% ( (gridstep*x, gridstep*y), )
            s += '{\n%s};\n' % D.tikz_boxed(scale=scale, boxsize=0.8*gridstep)
            i += 1
            if D.is_singular(): break
            D,func = D.reduce()
            fns.append(func)
        if labels == 'T':
            for j in range(1, i):
                edge = "edge node[midway, above] "
                edge += "{$T_%s$}"%j
                s += "\\path[thick, ->] (q%s) %s (q%s);\n"%(j-1, edge, j)
        elif labels:
            for j,func in zip(range(1, i), fns):
                rotate90 = "" if j%ncols==0 else ", rotate=90"
                edge = "edge node[midway, rectangle, fill=white%s] "%rotate90
                edge += "{$\\%s$}" % func
                s += "\\path[->] (q%s) %s (q%s);\n"%(j-1, edge, j)
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

    def tikz_commutative_diagram(self, tile, N=1, scale=(1,1), labels=True,
            newcommand=True):
        r"""
        Return a tikz commutative diagram for the composition.

        INPUT:

        - ``tile`` - WordMorphism, a square tile
        - ``N`` - integer (default:``1``), length of the diagram
        - ``scale`` - tuple of number (default:``(1,1)``), one for each line
        - ``labels`` - arrow labels (default:``True``). It may take the
          following values:

            - ``True`` - prints TRIM, SWAP, etc.
            - ``'T'`` - prints T_i, etc.
            - ``False`` - print nothing

        - ``newcommand`` - bool (default: ``True``), whether newcommand
          which defines ``\SWAP``, ``\TRIM``, etc.

        EXAMPLES:

        The following command creates the tikz code for Figure 16 in
        [BGL2012]_::

            sage: from sage.combinat.double_square_tile import DoubleSquare
            sage: fibo2 = words.fibonacci_tile(2)
            sage: S = WordMorphism({0:[0,0],1:[1,0,1],2:[2,2],3:[3,2,3]}, codomain=fibo2.parent())
            sage: cfibo2 = DoubleSquare(fibo2)
            sage: options = dict(tile=S,N=3,scale=(0.25,0.15),labels=True,newcommand=True)
            sage: s = cfibo2.tikz_commutative_diagram(**options)     # long time (2s)
            sage: s                                                  # long time
            \newcommand{\TRIM}{\textsc{trim}}
            \newcommand{\EXTEND}{\textsc{extend}}
            \newcommand{\SWAP}{\textsc{swap}}
            \newcommand{\SHIFT}{\textsc{shift}}
            \newcommand{\REVERSE}{\textsc{reverse}}
            \begin{tikzpicture}
            [first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},
            second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]
            \node (q0) at (0, 0) {\begin{tikzpicture}
            [scale=0.250000000000000]
            ...
            \end{tikzpicture}};
            \path[thick, ->] (r0) edge node[midway, rectangle, fill=white, rotate=90] {$\SWAP_1$} (r1);
            \path[thick, ->] (r1) edge node[midway, rectangle, fill=white, rotate=90] {$\TRIM_1$} (r2);
            \path[thick, ->] (r2) edge node[midway, rectangle, fill=white, rotate=90] {$\TRIM_3$} (r3);
            \path[thick, ->] (q0) edge node[midway, left] {$\varphi$} (r0);
            \path[thick, ->] (q1) edge node[midway, left] {$\varphi$} (r1);
            \path[thick, ->] (q2) edge node[midway, left] {$\varphi$} (r2);
            \path[thick, ->] (q3) edge node[midway, left] {$\varphi$} (r3);
            \end{tikzpicture}

        """
        s = ''
        if newcommand:
            s += "\\newcommand{\\TRIM}{\\textsc{trim}}\n"
            s += "\\newcommand{\\EXTEND}{\\textsc{extend}}\n"
            s += "\\newcommand{\\SWAP}{\\textsc{swap}}\n"
            s += "\\newcommand{\\SHIFT}{\\textsc{shift}}\n"
            s += "\\newcommand{\\REVERSE}{\\textsc{reverse}}\n"
        s += "\\begin{tikzpicture}\n"
        s += "[first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},\n"
        s += "second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt}]\n"

        # top
        selfcomposed = DoubleSquare(tile(self.boundary_word()))
        for D,q,y,sc in zip([self, selfcomposed], 'qr', [0, -4], scale):
            i = 0
            fns = []
            while i <= N:
                s += '\\node (%s%s) '%(q,i)
                x,y = 4*i, y
                s += 'at %s '% ( (x, y), )
                s += '{\\begin{tikzpicture}\n'
                s += '[scale=%s]\n'%sc
                s += D.tikz_trajectory(arrow='')
                s += '\n\\end{tikzpicture}};\n'
                i += 1
                if D.is_singular(): break
                D,func = D.reduce()
                fns.append(func)
            if labels == 'T':
                for j in range(1, i):
                    edge = "edge node[midway, above] "
                    edge += "{$T_%s$}"%j
                    s += "\\path[thick, ->] (%s%s) %s (%s%s);\n"%(q,j-1,edge,q,j)
            elif labels:
                for j,func in zip(range(1, i), fns):
                    rotate90 = ", rotate=90"
                    edge = "edge node[midway, rectangle, fill=white%s] "%rotate90
                    edge += "{$\\%s$}"%func
                    s += "\\path[thick, ->] (%s%s) %s (%s%s);\n"%(q,j-1, edge,q,j)
        for j in range(i):
            edge = "edge node[midway, left] "
            edge += "{$\\varphi$}"
            s += "\\path[thick, ->] (q%s) %s (r%s);\n"%(j, edge, j)
        s += "\\end{tikzpicture}\n"
        return LatexExpr(s)

###############################
# Creation of Double Square from inputs
###############################
def find_square_factorisation(ds, factorisation=None, alternate=True):
    r"""
    Return a square factorisation of the double square ds, distinct from
    the given factorisation.

    INPUT:

    - ``ds`` - word, the boundary word of a square tile
    - ``factorisation`` - tuple (default: ``None``), a known factorisation
    - ``alternate`` - bool (default: ``True``), if True the search for the
      second factorisation is restricted to those who alternates with the
      first factorisation

    OUTPUT:

    tuple of four positions of a square factorisation

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import find_square_factorisation
        sage: find_square_factorisation(words.fibonacci_tile(0))
        (0, 1, 2, 3)
        sage: find_square_factorisation(words.fibonacci_tile(1))
        (0, 3, 6, 9)
        sage: find_square_factorisation(words.fibonacci_tile(2))
        (0, 13, 26, 39)
        sage: find_square_factorisation(words.fibonacci_tile(3))
        (0, 55, 110, 165)

    ::

        sage: f = find_square_factorisation(words.fibonacci_tile(3))
        sage: f
        (0, 55, 110, 165)
        sage: find_square_factorisation(words.fibonacci_tile(3),f)         # long time (6s)
        (34, 89, 144, 199)
        sage: find_square_factorisation(words.fibonacci_tile(3),f,False)   # long time (11s)
        (34, 89, 144, 199)

    ::

        sage: find_square_factorisation(words.christoffel_tile(4,5))
        (0, 7, 28, 35)
        sage: find_square_factorisation(words.christoffel_tile(4,5),_)
        (2, 27, 30, 55)

    TESTS::

        sage: find_square_factorisation(Words('abcd')('aaaaaa'))
        Traceback (most recent call last):
        ...
        ValueError: no square factorization found
        sage: find_square_factorisation(Words('abcd')('aaaaaa'),(1,2,3,4))
        Traceback (most recent call last):
        ...
        ValueError: no second square factorization found
    """
    e,n,w,s = ds.parent().alphabet()
    rot180 = WordMorphism({e:w,w:e,n:s,s:n},codomain=ds.parent())
    hat = lambda x:rot180(x).reversal()

    L = ds.length()
    half = L/2
    twice = ds * ds

    if factorisation and alternate:
        a,b,c,d = factorisation
        it = ((startA,endA) for startA in range(a+1,b) for endA in range(b+1,a+half))
    else:
        it = ((startA,endA) for startA in range(half) for endA in range(startA+1,startA+half+1) )

    for startA,endA in it:
        new = (startA,endA,(startA+half)%L,(endA+half)%L)
        if factorisation and set(factorisation) == set(new):
            continue
        A = twice[startA:endA]
        B = twice[endA:startA+half]
        A2 = twice[startA+half:endA+half]
        B2 = twice[endA+half:L+startA]
        if A == hat(A2) and B == hat(B2):
            return new

    if factorisation is None:
        raise ValueError, 'no square factorization found'
    else:
        raise ValueError, 'no second square factorization found'

def double_square_from_boundary_word(ds):
    r"""
    Creates a double square object from the boundary word of a double
    square tile.

    INPUT:

    - ``ds`` - word, the boundary of a double square. The parent alphabet
      is assumed to be in the order : East, North, West, South.

    OUTPUT:

    - tuple - tuple of 8 words over the alphabet A
    - WordMorphism, involution on the alphabet A and representing a
      rotation of 180 degrees.
    - dict - mapping letters of A to steps in the plane.

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import double_square_from_boundary_word
        sage: fibo = words.fibonacci_tile
        sage: W, rot180, steps = double_square_from_boundary_word(fibo(1))
        sage: map(len, W)
        [2, 1, 2, 1, 2, 1, 2, 1]
        sage: W, rot180, steps = double_square_from_boundary_word(fibo(2))
        sage: map(len, W)
        [8, 5, 8, 5, 8, 5, 8, 5]
        sage: W, rot180, steps = double_square_from_boundary_word(fibo(3))  # long time (6s)
        sage: map(len, W)                                                   # long time
        [34, 21, 34, 21, 34, 21, 34, 21]
        sage: rot180                                                        # long time
        WordMorphism: 0->2, 1->3, 2->0, 3->1

    """
    # Define rot180
    parent = ds.parent()
    alphabet = parent.alphabet()
    if not hasattr(alphabet,'cardinality') or alphabet.cardinality() != 4:
        raise ValueError, "The parent of ds must have a 4-letter alphabet."
    e,n,w,s = alphabet
    rot180 = WordMorphism({e:w,w:e,n:s,s:n},codomain=parent)

    # Define steps
    steps = {}
    steps[e] = vector((1,0))
    steps[n] = vector((0,1))
    steps[w] = vector((-1,0))
    steps[s] = vector((0,-1))

    # Compute the wi
    f = find_square_factorisation(ds)
    g = find_square_factorisation(ds, f, alternate=True)
    cuts = sorted(f+g)
    w0 = ds[cuts[0]:cuts[1]]
    w1 = ds[cuts[1]:cuts[2]]
    w2 = ds[cuts[2]:cuts[3]]
    w3 = ds[cuts[3]:cuts[4]]
    w4 = ds[cuts[4]:cuts[5]]
    w5 = ds[cuts[5]:cuts[6]]
    w6 = ds[cuts[6]:cuts[7]]
    w7 = ds[cuts[7]:] + ds[:cuts[0]]
    W = (w0,w1,w2,w3,w4,w5,w6,w7)

    return W, rot180, steps

def double_square_from_four_integers(l0, l1, l2, l3):
    r"""
    Creates a double square from the lengths of the `w_i`.

    INPUT:

    - ``l0`` - integer, length of `w_0`
    - ``l1`` - integer, length of `w_1`
    - ``l2`` - integer, length of `w_2`
    - ``l3`` - integer, length of `w_3`

    OUTPUT:

    - tuple - tuple of 8 words over alphabet A
    - WordMorphism, involution on the alphabet A and representing a
      rotation of 180 degrees.
    - dict - mapping letters of A to steps in the plane.

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import double_square_from_four_integers
        sage: w,rot180,steps = double_square_from_four_integers(2,1,2,1)
        sage: w
        (Path: 21, Path: 2, Path: 32, Path: 3, Path: 03, Path: 0, Path: 10, Path: 1)
        sage: rot180
        WordMorphism: 0->2, 1->3, 2->0, 3->1
        sage: sorted(steps.items())
        [(0, (1, 0)), (1, (0, 1)), (2, (-1, 0)), (3, (0, -1))]

    If the input integers do not define a double square uniquely, the
    alphabet might be larger than 8::

        sage: w,rot180,steps = double_square_from_four_integers(4,2,4,2)
        sage: w
        (Path: 7601,
         Path: 76,
         Path: 5476,
         Path: 54,
         Path: 2354,
         Path: 23,
         Path: 0123,
         Path: 01)
        sage: rot180
        WordMorphism: 0->4, 1->5, 2->6, 3->7, 4->0, 5->1, 6->2, 7->3
        sage: sorted(steps.items())
        [(0, (1, 0)),
         (1, (1/2*sqrt(2), 1/2*sqrt(2))),
         (2, (0, 1)),
         (3, (-1/2*sqrt(2), 1/2*sqrt(2))),
         (4, (-1, 0)),
         (5, (-1/2*sqrt(2), -1/2*sqrt(2))),
         (6, (0, -1)),
         (7, (1/2*sqrt(2), -1/2*sqrt(2)))]
    """
    # Initial words where every letter is assumed to different
    # The involution exchanges the sign
    somme = l0 + l1 + l2 + l3
    W = Words(range(1, somme+1) + range(-1, -(somme+1), -1))
    w0 = W(range(1, l0+1))
    w1 = W(range(l0+1,l0+l1+1))
    w2 = W(range(l0+l1+1,l0+l1+l2+1))
    w3 = W(range(l0+l1+l2+1,l0+l1+l2+l3+1))
    inv = lambda x : -x
    rot180 = WordMorphism(dict((x,inv(x)) for x in W.alphabet()), codomain=W)
    hat = lambda w: rot180(w).reversal()
    hatA = hat(w0 *w1)
    hatB = hat(w2 * w3)
    w4 = hatA[:l0]
    w5 = hatA[l0:]
    w6 = hatB[:l2]
    w7 = hatB[l2:]
    #print w0,w1,w2,w3,w4,w5,w6,w7

    # Creation of the overlap and of the resulting disjoint set
    p = hat(w1 * w2).overlap_partition(w5 * w6, 0, involution = inv)
    p = hat(w3 * w4).overlap_partition(w7 * w0, 0, p, involution = inv)

    # The alphabet representents and resulting new involution
    alphabet = p.root_to_elements_dict().keys()
    proj1 = p.element_to_root_dict()
    dbar = dict((a,proj1[inv(a)]) for a in alphabet)
    nbar = WordMorphism(dbar)
    A, B, C = nbar.partition_of_domain_alphabet()
    assert len(C) == 0
    len_A = len(A)
    steps = []
    for j in range(len_A):
        angle = pi * j / len_A
        steps.append( (cos(angle), sin(angle)) )
    A_ranged = range(len_A)
    AB_ranged = range(2 * len_A)

    # Projection of the initial words into words over the new alphabet
    from sage.combinat.words.paths import WordPaths
    P = WordPaths(AB_ranged, steps=steps)
    steps_dict = P.letters_to_steps()
    proj2 = dict(zip(A,A_ranged))
    proj2.update((dbar[a], proj2[a] + len_A) for a in A)
    Projection = WordMorphism(proj2, codomain=P) * WordMorphism(proj1)
    w0,w1,w2,w3,w4,w5,w6,w7 = map(Projection, (w0,w1,w2,w3,w4,w5,w6,w7))
    nbar = WordMorphism( dict((a,(a + len_A) % (2 * len_A)) for a in AB_ranged), codomain=P)
    #print w0,w1,w2,w3,w4,w5,w6,w7
    return ( (w0,w1,w2,w3,w4,w5,w6,w7), nbar, steps_dict)

###############################
# Creation of Figure 11 for [BGL2012]_
###############################
def figure_11_BGL2012(scale=0.5, boxsize=10, newcommand=True):
    r"""
    Return the tikz code of the Figure 11 for the article [BGL2012]_.

    INPUT:

    - ``scale`` - number (default: ``0.5``), tikz scale
    - ``boxsize`` - integer (default: ``10``), size of box the double
      squares must fit in
    - ``newcommand`` - bool (default: ``True``), whether to include
      latex newcommand for TRIM, EXTEND and SWAP

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import figure_11_BGL2012
        sage: s = figure_11_BGL2012()
        sage: s
        \newcommand{\TRIM}{\textsc{trim}}
        \newcommand{\EXTEND}{\textsc{extend}}
        \newcommand{\SWAP}{\textsc{swap}}
        \begin{tikzpicture}
        [first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},
        second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt},
        >=latex,
        node distance=3cm]
        \node (A) at (15,0)
        {
        \begin{tabular}{c}
        \begin{tikzpicture}
        [scale=0.5]
        ...
        \end{tikzpicture}
        };
        \path[<-] (A) edge node[midway, rectangle, fill=white] {$\TRIM_2$} (B);
        \path[<-] (B) edge node[midway, rectangle, fill=white] {$\TRIM_0$} (C);
        \path[<-] (C) edge node[midway, rectangle, fill=white] {$\TRIM_1$} (D);
        \path[<-] (D) edge node[midway, rectangle, fill=white] {$\TRIM_3$} (E);
        \path[<-] (E) edge node[midway, rectangle, fill=white] {$\SWAP_1$} (F);
        \path[<-] (F) edge node[midway, rectangle, fill=white] {$\TRIM_1$} (G);
        \path[<-] (G) edge node[midway, rectangle, fill=white,rotate=90] {$\TRIM_3$} (H);
        \path[<-] (H) edge node[midway, rectangle, fill=white,rotate=90] {$\SWAP_1$} (I);
        \path[<-] (D) edge node[midway, rectangle, fill=white] {$\TRIM_2$} (E2);
        \end{tikzpicture}

    """
    from sage.combinat.words.word_generators import words
    fibo = words.fibonacci_tile
    Dfibo2 = DoubleSquare(fibo(2))
    Dfibo2 = Dfibo2.reverse().shift().shift().shift()
    F = Dfibo2
    E = Dfibo2.swap(1)
    D = Dfibo2.swap(1).trim(3)
    C = Dfibo2.swap(1).trim(3).trim(1)
    B = Dfibo2.swap(1).trim(3).trim(1).trim(0)
    A = Dfibo2.swap(1).trim(3).trim(1).trim(0).trim(2)
    G = Dfibo2.extend(1)
    H = Dfibo2.extend(1).extend(3)
    I = Dfibo2.extend(1).extend(3).swap(1)
    E2 = D.extend(2)
    s = ''
    if newcommand:
        s += "\\newcommand{\\TRIM}{\\textsc{trim}}\n"
        s += "\\newcommand{\\EXTEND}{\\textsc{extend}}\n"
        s += "\\newcommand{\\SWAP}{\\textsc{swap}}\n"
    s += "\\begin{tikzpicture}\n"
    s += "[first/.style={circle,draw=black,fill=black, inner sep=0pt, minimum size=3pt},\n"
    s += "second/.style={circle,draw=black,fill=white, inner sep=0pt, minimum size=3pt},\n"
    s += ">=latex,\n"
    s += "node distance=3cm]\n"

    dx = 5
    dy = 5

    s += '\\node (A) at (%s,%s)\n' % (3*dx, 0)
    s += '{\n%s};\n' % A.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (B) at (%s,%s)\n' % (2*dx, 0)
    s += '{\n%s};\n' % B.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (C) at (%s,%s)\n' % (dx, 0)
    s += '{\n%s};\n' % C.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (D) at (0,0)\n'
    s += '{\n%s};\n' % D.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (E) at (%s,%s)\n' % (dx, dy)
    s += '{\n%s};\n' % E.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (E2) at (%s,%s)\n' % (0, dy)
    s += '{\n%s};\n' % E2.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (F) at (%s,%s)\n' % (2.6*dx, dy)
    s += '{\n%s};\n' % F.tikz_boxed(scale=scale, boxsize=boxsize)
    s += '\\node (G) at (%s,%s)\n' % (2.6*dx, 2.1*dy)
    s += '{\n%s};\n' % G.tikz_boxed(scale=scale, boxsize=0.7*boxsize)
    s += '\\node (H) at (%s,%s)\n' % (1.3*dx, 2.1*dy)
    s += '{\n%s};\n' % H.tikz_boxed(scale=scale, boxsize=0.7*boxsize)
    s += '\\node (I) at (%s,%s)\n' % (0.1*dx, 2.1*dy)
    s += '{\n%s};\n' % I.tikz_boxed(scale=scale, boxsize=0.7*boxsize)

    s += "\\path[<-] (A) edge node[midway, rectangle, fill=white] {$\\TRIM_2$} (B);\n"
    s += "\\path[<-] (B) edge node[midway, rectangle, fill=white] {$\\TRIM_0$} (C);\n"
    s += "\\path[<-] (C) edge node[midway, rectangle, fill=white] {$\\TRIM_1$} (D);\n"
    s += "\\path[<-] (D) edge node[midway, rectangle, fill=white] {$\\TRIM_3$} (E);\n"
    s += "\\path[<-] (E) edge node[midway, rectangle, fill=white] {$\\SWAP_1$} (F);\n"
    s += "\\path[<-] (F) edge node[midway, rectangle, fill=white] {$\\TRIM_1$} (G);\n"
    s += "\\path[<-] (G) edge node[midway, rectangle, fill=white,rotate=90] {$\\TRIM_3$} (H);\n"
    s += "\\path[<-] (H) edge node[midway, rectangle, fill=white,rotate=90] {$\\SWAP_1$} (I);\n"
    s += "\\path[<-] (D) edge node[midway, rectangle, fill=white] {$\\TRIM_2$} (E2);\n"
    s += "\\end{tikzpicture}\n"
    return LatexExpr(s)
###############################
# Usefull Functions
###############################
def snake(i, ncols=2):
    r"""
    Return the coordinate of the ith node of a snake.

    This is used for the tikz drawing of a double square reduction.

    INPUT:

    - ``i`` - integer, the ith node
    - ``ncols`` - integer (default ``2``), number of columns

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import snake
        sage: for i in range(8): snake(i, 3)
        (0, 0)
        (1, 0)
        (2, 0)
        (2, -1)
        (1, -1)
        (0, -1)
        (0, -2)
        (1, -2)
    """
    y = -(i // ncols)
    x = i % ncols
    if y % 2 == 1:
        x = ncols - 1  - x
    return (x, y)

###############################
# Triple square factorisation examples
###############################
def triple_square_example(i):
    r"""
    Return a triple square factorisation example.

    These words having three square factorisations were provided by Xavier
    Provençal.

    INPUT:

    - ``i`` - integer, accepted values are 1, 2 or 3.

    EXAMPLES::

        sage: from sage.combinat.double_square_tile import triple_square_example
        sage: triple_square_example(1)
        Path: abaBAbabaBAbabaBAbabABABabABABabABAB
        sage: triple_square_example(2)
        Path: abaBABaabaBABaabaBABaabABAAbabABAAbabABA...
        sage: triple_square_example(3)
        Path: aabAAbaabAAbaabAAbaaBAABaaBAABaaBAAB

    Triple square tile do not exist. Hence the example provided by Xavier
    Provençal can not be the boundary word of a tile. One can see it by
    ploting it or by the fact that the turning number is zero::

        sage: from sage.combinat.double_square_tile import DoubleSquare
        sage: D = DoubleSquare(triple_square_example(1))
        sage: D
        Double Square Tile
          w0 = a          w4 = a
          w1 = baBA       w5 = bABA
          w2 = babaB      w6 = BabAB
          w3 = AbabaBAb   w7 = ABabABAB
        (|w0|, |w1|, |w2|, |w3|) = (1, 4, 5, 8)
        (d0, d1, d2, d3)         = (12, 6, 12, 6)
        (n0, n1, n2, n3)         = (0, 0, 0, 1)
        sage: D.turning_number()
        0
    """
    from sage.combinat.words.paths import WordPaths
    WP = WordPaths('abAB')
    if i == 1:
        return WP('abaBAbabaBAbabaBAbabABABabABABabABAB')
    elif i == 2:
        return WP('abaBABaabaBABaabaBABaabABAAbabABAAbabABAAb')
    elif i == 3:
        return WP('aabAAbaabAAbaabAAbaaBAABaaBAABaaBAAB')
    else:
        raise ValueError("i (=%s) must be 1, 2 or 3" % i)

