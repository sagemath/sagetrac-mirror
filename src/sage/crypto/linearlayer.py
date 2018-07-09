r"""
Linear layer Representation

Linear layers play an important role in block cipher designs, especially in
substitution-permutation-networks (SPNs).

Available classes are the generic ``LinearLayer`` and ``AESLikeLinearLayer``.

Apart from that, there are also linear layers available which are used in the
literature, either in  the dictionary linearlayers, or as a seperate object.
This module provides the following linear layers:

    - AES ([DR2002]_)
    - Midori ([BBISHAR2015]_)
    - SKINNY ([BJKLMPSSS2016]_)
    - PRESENT (and SmallScalePRESENT) ([BKLPPRSV2007]_)

AUTHORS:

- Friedrich Wiemer (2018-07-02): initial version
"""

from six import integer_types

from sage.combinat.permutation import Permutation
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_mod2_dense import Matrix_mod2_dense
from sage.matrix.matrix_gf2e_dense import Matrix_gf2e_dense
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.structure.sage_object import SageObject

from sage.misc.superseded import experimental


def _branch_number(mtr):
    """
    Computes the branch number of the given matrix
    """
    from sage.coding.linear_code import LinearCode
    from sage.matrix.special import identity_matrix

    F = mtr.base_ring()
    n = mtr.nrows()
    I = identity_matrix(F, n)
    generator_matrix = Matrix(F, [list(a) + list(b) for a, b in zip(I, mtr)])
    return LinearCode(generator_matrix).minimum_distance()


def _ff_elem_to_binary(elem):
    """
    Converts a finite field element to the binary matrix carrying out the according multiplication
    """
    from sage.matrix.special import companion_matrix, identity_matrix, zero_matrix

    F = elem.parent()
    R = companion_matrix(F.modulus(), format='right')
    bits = ZZ(elem.integer_representation()).digits(base=2, padto=F.degree())
    result = zero_matrix(GF(2), F.degree())
    T = identity_matrix(GF(2), F.degree())
    for i in bits:
        if i==1:
            result += T
        T = T*R
    return result


def _ff_matrix_to_binary(mtr):
    """
    Converts a matrix over a finite field to the binary matrix carrying out the according multiplication
    """
    from sage.matrix.special import block_matrix

    result = []
    for row in mtr:
        for elem in row:
            result.append(_ff_elem_to_binary(elem))
    return block_matrix(mtr.nrows(), mtr.ncols(), result)


def _column_linear_layer(Ls):
    """
    Return the linear layer that applies a given linear layer ``L`` to ``ncols`` in parallel.

    The most common application for this is, to build the matrix that applies the
    AES MixColumns linear layer to each column of a state, while the state is represented
    as a vector.
    """
    from sage.matrix.special import block_matrix, diagonal_matrix

    F = Ls[0].base_ring()
    nblocks = len(Ls)
    n, m = Ls[0].dimensions()
    blockmtrs = [diagonal_matrix(F, nblocks, [Ls[i][k][l] for i in range(nblocks)]) for k in range(n) for l in range(m)]
    return LinearLayer.constructor(block_matrix(F, n, m, blockmtrs))


class LinearLayer:
    r"""
    Many modern block cipher constructions in symmetric cryptography use linear layers
    as one of their basic building blocks. A linear layer is typically used to spread
    locally diffused bits over the whole state of the block cipher. As linear layers
    are, well, linear, and operate on bits, nibbles, or bytes, we can represent it as
    a matrix over `\GF(2)^n`, `\GF(2^4)^n`, or `\GF(2^8)^n`. Application of the linear
    layer to x then just corresponds to left multiplication of this matrix: `A \cdot x`.

    EXAMPLES:

    We start with the simple identity linear layer::

        sage: from sage.crypto.linearlayer import LinearLayer
        doctest:warning
        ...
        FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See https://trac.sagemath.org/25735 for details.
        sage: id = LinearLayer.constructor(identity_matrix(GF(2), 2)); id
        LinearLayer of dimension 2 x 2 represented as
        [1 0]
        [0 1]

        sage: id.is_permutation()
        True

    A list (dictionary) of linear layers used in the literature is also available::

        sage: from sage.crypto.linearlayer import linearlayers
        sage: linearlayers['PRESENT']
        LinearLayer of dimension 64 x 64 represented as
        64 x 64 dense matrix over Finite Field of size 2 (use the '.str()' method to see the entries)
        sage: linearlayers['PRESENT'].is_permutation()
        True
    """
    @staticmethod
    @experimental(25735)
    def constructor(*args,  **kwargs):
        """
        Construct a linear layer for a given matrix `L` of dimension m x n.

        INPUT:

        - ``L`` - a matrix representing the linear layer part.
        """
        def LinearLayerFactory(K):
            if K.characteristic() == 2 and K.degree() == 1:
                return type("LinearLayerGF2", (LinearLayer, Matrix_mod2_dense), {})
            if K.characteristic() == 2 and K.degree() >= 1:
                return type("LinearLayerGF2E", (LinearLayer, Matrix_gf2e_dense), {})
            else:
                raise NotImplementedError

        if "L" in kwargs:
            L = kwargs["L"]
        elif len(args) == 1:
            L = args[0]
        else:
            raise TypeError("No matrix L as argument provided.")

        parent = L.parent()
        base_ring = L.base_ring()
        return LinearLayerFactory(base_ring)(parent, L)

    def _latex_(self):
        r"""
        Returns a `LaTeX` version of the operation table as a string,
        using a `LaTeX` ``array`` environment.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: id = LinearLayer.constructor(identity_matrix(GF(2), 2))
            sage: id._latex_()
            'x \\ {\\mapsto} \\left(\\begin{array}{rr}\n1 & 0 \\\\\n0 & 1\n\\end{array}\\right) {\\cdot} \\ x'
        """

        return "x \\ {\\mapsto} " + self.matrix()._latex_() + " {\\cdot} \\ x"

    def __str__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: ll = LinearLayer.constructor(matrix(GF(2), [[0,1,0],[1,0,1]]))
            sage: print(ll)
            LinearLayer of dimension 2 x 3 represented as
            [0 1 0]
            [1 0 1]

            sage: from sage.crypto.linearlayer import linearlayers
            sage: print(linearlayers['PRESENT'])
            LinearLayer of dimension 64 x 64 represented as
            [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
            ...
            [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
        """
        return "LinearLayer of dimension %d x %d represented as\n%s" \
                % (self.dimensions() + (self.matrix().str(),))

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: LinearLayer.constructor(matrix(GF(2), [[0,1,0],[1,0,1]]))  # indirect doctest
            LinearLayer of dimension 2 x 3 represented as
            [0 1 0]
            [1 0 1]

            sage: from sage.crypto.linearlayer import linearlayers
            sage: linearlayers['PRESENT']
            LinearLayer of dimension 64 x 64 represented as
            64 x 64 dense matrix over Finite Field of size 2 (use the '.str()' method to see the entries)
        """
        return "LinearLayer of dimension %d x %d represented as\n%s" \
                % (self.dimensions() + (self.matrix().__repr__(),))

    def matrix(self):
        """
        Returns the matrix representing this linear layer
        """
        return self.parent(self._matrix_())

    def binary_matrix(self):
        """
        Returns the matrix representing this linear layer in it's binary representation
        """
        if self.base_ring() is GF(2):
            return self._matrix_()
        return _ff_matrix_to_binary(self._matrix_())

    def __call__(self, x):
        """
        Apply linear layer to ``x``.

        INPUT:

        - ``x`` - either an integer, a tuple of `\GF{2}` elements of
          length ``len(self)`` or a finite field element in
          `\GF{2^n}`. As a last resort this function tries to convert
          ``x`` to an integer.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: L = LinearLayer.constructor(Matrix(GF(2), [[0,1,0],[1,0,1]]))
            sage: L(3)
            3

            sage: L((0, 1, 1))
            (1, 1)

            sage: L([0, 1, 1])
            [1, 1]

            sage: L(vector(GF(2), [0, 1, 1]))
            (1, 1)

            sage: L([0]*(L.ncols()+1))
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply LinearLayer to provided element, dimension mismatch.

            sage: L([0, 1/2, 1])
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply LinearLayer to provided element [0, 1/2, 1], conversion to vector failed.

            sage: L("failure")
            Traceback (most recent call last):
            ...
            TypeError: Unsupported type for input x: <type 'str'>.
        """
        from sage.modules.free_module_element import FreeModuleElement

        if isinstance(x, integer_types + (Integer,)):
            x = vector(GF(2), ZZ(x).digits(base=2, padto=self.ncols()))
            return ZZ(list(self.binary_matrix() * x), 2)

        elif isinstance(x, tuple):
            if len(x) != self.ncols():
                raise TypeError("Cannot apply LinearLayer to provided element, dimension mismatch.")

            try:
                x = vector(self.matrix().base_ring(), x)
            except (TypeError, ZeroDivisionError):
                raise TypeError("Cannot apply LinearLayer to provided element %r, conversion to vector failed." % (x,))

            return tuple(self * x)

        elif isinstance(x, list):
            if len(x) != self.ncols():
                raise TypeError("Cannot apply LinearLayer to provided element, dimension mismatch.")

            try:
                x = vector(self.base_ring(), x)
            except (TypeError, ZeroDivisionError):
                raise TypeError("Cannot apply LinearLayer to provided element %r, conversion to vector failed." % (x,))

            return list(self * x)

        elif isinstance(x, FreeModuleElement):
            return self * x

        else:
            raise TypeError("Unsupported type for input x: %s." % type(x))

    def is_permutation(self):
        r"""
        Check if the linear layer is a permutation, i.e., if the representing matrix is a permutation.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: L1 = LinearLayer.constructor(identity_matrix(GF(2), 4, 4))
            sage: L1.is_permutation()
            True

            sage: L2 = LinearLayer.constructor(Matrix(GF(2), [[0,1,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]]))
            sage: L2.is_permutation()
            False

            sage: L3 = LinearLayer.constructor(Matrix(GF(2), [[0,1,0], [1,0,1]]))
            sage: L3.is_permutation()
            False
        """
        # a permutation matrix has to be a square
        if not self.is_square():
            return False

        # in each row, all entries execpt one should be 0, the other should be 1
        for r in self.rows():
            if r.hamming_weight() != 1:
                return False
            if r[r.nonzero_positions()[0]] != 1:
                return False

        # in each column, all entries execpt one should be 0, the other should be 1
        for c in self.columns():
            if c.hamming_weight() != 1:
                return False
            if c[c.nonzero_positions()[0]] != 1:
                return False

        return True

    def naive_xor_count(self, algorithm="naive"):
        """
        Counts the number of xor operations needed for a naive implementation

        INPUT:

        - ``algorithm`` - string choosing which algorithm to use to compute
            the xor count. Available algorithms

            - 'naive' - non-optimized implementation (default)
        """
        avail_algs = ["naive"]
        if not algorithm in avail_algs:
            raise TypeError("algorithm must be one of %s" % avail_algs)

        if algorithm is "naive":
            mtr = self.binary_matrix()
            n, m = mtr.dimensions()
            return mtr.density() * n * m - n

    @cached_method
    def differential_branch_number(self):
        """
        Computes the differential branch number of the linear layer

        EXAMPLES::

            sage: from sage.crypto.linearlayer import PRESENT
            sage: PRESENT.differential_branch_number()
            2
        """
        if self.is_permutation():
            return 2
        return _branch_number(self.matrix())

    @cached_method
    def linear_branch_number(self):
        """
        Computes the linear branch number of the linear layer

        EXAMPLES::

            sage: from sage.crypto.linearlayer import PRESENT
            sage: PRESENT.linear_branch_number()
            2
        """
        if self.is_permutation():
            return 2
        return _branch_number(self.matrix().transpose())


class AESLikeLinearLayer(LinearLayer, Matrix_gf2e_dense):
    r"""
    An AES like linear layer consists of ShiftRows and MixColumns matrices,
    corresponding to the linear layer structure used in the AES.

    Here, ShiftRows is either a Permutation object, or a permutation matrix,
    and MixColumns any arbitrary matrix.

    EXAMPLES::

        sage: from sage.crypto.linearlayer import linearlayers
        sage: linearlayers['AES']
        AES like LinearLayer of dimension 128 x 128 represented by ShiftRows
        [1, 6, 11, 16, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12]
        and MixColumns
        [    x x + 1     1     1]
        [    1     x x + 1     1]
        [    1     1     x x + 1]
        [x + 1     1     1     x]
    """
    @staticmethod
    @experimental(25735)
    def constructor(sr, mc, apply_columns=True):
        """
        Construct an AES linear layer for a given matrix `L` of dimension m x n.

        INPUT:

        - ``ShiftRows`` - a permutation or matrix representing the shift rows part
        - ``MixColumns`` - a matrix representing the mix columns part
        - ``apply_columns`` - Bool; wether to apply MixColumns column-wise (default)
            or row-wise

        EXAMPLES::

            sage: from sage.crypto.linearlayer import linearlayers
            sage: linearlayers['AES']
            AES like LinearLayer of dimension 128 x 128 represented by ShiftRows
            [1, 6, 11, 16, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12]
            and MixColumns
            [    x x + 1     1     1]
            [    1     x x + 1     1]
            [    1     1     x x + 1]
            [x + 1     1     1     x]
        """
        K = mc.base_ring()
        if not (K.characteristic() == 2 and K.degree() >= 2):
            raise NotImplementedError

        if isinstance(sr, Permutation):
            sr = Matrix(K, sr.to_matrix())

        m = Matrix(K, _column_linear_layer([mc]*mc.ncols()))*sr
        if not apply_columns:
            raise NotImplementedError

        ll = AESLikeLinearLayer(m.parent(), m)
        ll._sr = sr
        ll._mc = mc

        return ll

    def _latex_(self):
        raise NotImplementedError

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.crypto.linearlayer import Midori
            sage: Midori  # indirect doctest
            AES like LinearLayer of dimension 64 x 64 represented by ShiftRows
            [1, 11, 6, 16, 15, 5, 12, 2, 10, 4, 13, 7, 8, 14, 3, 9]
            and MixColumns
            [0 1 1 1]
            [1 0 1 1]
            [1 1 0 1]
            [1 1 1 0]
        """
        # convert ShiftRows permutation matrix to a permutation object
        perm_sr = Permutation(map(lambda x: 1+list(x).index(1), self._sr.columns()))
        n, m = self.dimensions()
        degree = self.base_ring().degree()
        return "AES like LinearLayer of dimension %d x %d represented by ShiftRows\n%s\nand MixColumns\n%s" \
               % (n*degree, m*degree, perm_sr, self._mc)

    @cached_method
    def differential_branch_number(self):
        """
        Computes the differential branch number of the MixColumn Matrix

        EXAMPLES::

            sage: from sage.crypto.linearlayer import Midori
            sage: Midori.differential_branch_number()
            4
        """
        M = self._mc
        return _branch_number(M)

    @cached_method
    def linear_branch_number(self):
        """
        Computes the linear branch number of the MixColumn Matrix

        EXAMPLES::

            sage: from sage.crypto.linearlayer import Midori
            sage: Midori.linear_branch_number()
            4
        """
        #if self.is_permutation():
        #    return 2
        M = self._mc
        M = M.transpose()
        return _branch_number(M)


Left_ShiftRows = Permutation([1, 6, 11, 16, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12])
Right_ShiftRows = Permutation([1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3, 16, 13, 10, 7, 4])

_AES_irreducible_polynomial = PolynomialRing(GF(2), name="alpha")("alpha^8 + alpha^4 + alpha^3 + alpha + 1")
_AES_field = GF(2**8, name="x", modulus=_AES_irreducible_polynomial)
AES_ShiftRows = Left_ShiftRows
AES_MixColumns = Matrix(_AES_field, 4, 4, map(_AES_field.fetch_int, [2, 3, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 3, 1, 1, 2]))
AES = AESLikeLinearLayer.constructor(AES_ShiftRows, AES_MixColumns)

Midori_ShuffelCells = Permutation([1, 11, 6, 16, 15, 5, 12, 2, 10, 4, 13, 7, 8, 14, 3, 9])
Midori_MixColumns = Matrix(GF(2**4), [[0, 1, 1, 1] ,[1, 0, 1, 1] ,[1, 1, 0, 1] ,[1, 1, 1, 0]])
Midori = AESLikeLinearLayer.constructor(Midori_ShuffelCells, Midori_MixColumns)

SKINNY_ShiftRows = Right_ShiftRows
SKINNY_4_MixColumns = Matrix(GF(2**4), [[1, 0, 1, 1] ,[1, 0, 0, 0] ,[0, 1, 1, 0] ,[1, 0, 1, 0]])
SKINNY_8_MixColumns = Matrix(GF(2**8), [[1, 0, 1, 1] ,[1, 0, 0, 0] ,[0, 1, 1, 0] ,[1, 0, 1, 0]])
SKINNY_4 = AESLikeLinearLayer.constructor(SKINNY_ShiftRows, SKINNY_4_MixColumns)
SKINNY_8 = AESLikeLinearLayer.constructor(SKINNY_ShiftRows, SKINNY_8_MixColumns)


def smallscale_present_linearlayer(nsboxes=16):
    """
    The matrix representing SmallPRESENT with nsboxes many S-boxes.

    Original PRESENT corresponds to nsboxes=16.

    INPUT:

    - ``nsboxes`` - integer, number of sboxes the linear layer operates on
      (default: 16).

    TESTS:

        sage: from sage.crypto.linearlayer import smallscale_present_linearlayer
        sage: smallscale_present_linearlayer(4)
        LinearLayer of dimension 16 x 16 represented as
        [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]
        [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0]
        [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0]
        [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]

    """
    from sage.modules.free_module import VectorSpace
    from sage.modules.free_module_element import vector

    def present_llayer(n, x):
        dim = 4*n
        y = [0]*dim
        for i in range(dim-1):
            y[i] = x[(n * i) % (dim - 1)]
        y[dim-1] = x[dim-1]
        return vector(GF(2), y)

    m = Matrix(GF(2), [present_llayer(nsboxes, ei)
                       for ei in VectorSpace(GF(2), 4*nsboxes).basis()])
    return LinearLayer.constructor(m)

PRESENT = smallscale_present_linearlayer(nsboxes=16)


# Dictionary of all available linear layers
linearlayers = {}
import sys
for k in dir(sys.modules[__name__]):
    v = getattr(sys.modules[__name__], k)
    if isinstance(v, (LinearLayer, AESLikeLinearLayer)):
        linearlayers[k] = v
