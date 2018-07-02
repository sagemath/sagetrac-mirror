r"""
LinearLayer Representation

LinearLayers play an important role in block cipher designs, especially in
substitution-permutation-networks (SPNs).
"""

from six import integer_types

from sage.combinat.permutation import Permutation
from sage.matrix.constructor import matrix
from sage.matrix.special import block_matrix, block_diagonal_matrix, companion_matrix, \
    diagonal_matrix, identity_matrix, zero_matrix
from sage.misc.cachefunc import cached_method
from sage.modules.free_module_element import vector
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.structure.sage_object import SageObject

from sage.misc.superseded import experimental


def _branch_number(mtr):
    """
    Computes the branch number of the given matrix
    """
    from sage.coding.linear_code import LinearCode

    F = mtr.base_ring()
    n = mtr.nrows()
    I = identity_matrix(F, n)
    generator_matrix = matrix(F, [list(a) + list(b) for a, b in zip(I, mtr)])
    return LinearCode(generator_matrix).minimum_distance()


def _ff_elem_to_binary(elem):
    """
    Converts a finite field element to the binary matrix carrying out the according multiplication
    """
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
    result = []
    for row in mtr:
        for elem in row:
            result.append(_ff_elem_to_binary(elem))
    return block_matrix(mtr.nrows(), mtr.ncols(), result)


class LinearLayer(SageObject):
    r"""
    Many modern block cipher constructions in symmetric cryptography use linear layers
    as one of their basic building blocks. A linear layer is typically used to spread
    locally diffused bits over the whole state of the block cipher. As linear layers
    are, well, linear, and operate on bits, nibbles, or bytes, we can represent it as
    a matrix over `\GF(2)^n`, `\GF(2^4)^n`, or `\GF(2^8)^n`. Application of the linear
    layer to x then just corresponds to left multiplication of this matrix: `A*x`.

    EXAMPLES:

    We start with the simple identity linear layer::

        sage: from.sage.crypto.linearlayer import LinearLayer
        sage: id = LinearLayer(identity_matrix(GF(2), 2)); id
        LinearLayer of dimension 2 x 2 represented as
        [1 0]
        [0 1]

        sage: id.is_permutation()
        True

    A list of linear layers used in the literature is also available::

        sage: from.sage.crypto.linearlayers import linearlayers
        sage: linearlayers['PRESENT']
        LinearLayer of dimension 64 x 64 represented as
        [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ...
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]
        sage: PRESENT.is_permutation()
        True

    AUTHORS:

    - Friedrich Wiemer (2018-07-02): initial version
    """

    @experimental(25735)
    def __init__(self, *args,  **kwargs):
        """
        Construct a linear layer for a given matrix `L` of dimension m times n.

        INPUT:

        - ``L`` - a matrix representing the linear layer with finite field elements

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: LinearLayer(matrix(GF(2), [[0,1,0],[1,0,1]]))  # indirect doctest
            LinearLayer of dimension 2 x 3 represented as
            [0 1 0]
            [1 0 1]
        """
        if "L" in kwargs:
            L = kwargs["L"]
        elif len(args) == 1:
            L = args[0]
        else:
            raise TypeError("No matrix L as argument provided.")

        self._L = L
        self.m = L.nrows()
        self.n = L.ncols()

    def _latex_(self):
        raise NotImplementedError

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: LinearLayer(matrix(GF(2), [[0,1,0],[1,0,1]]))  # indirect doctest
            LinearLayer of dimension 2 x 3 represented as
            [0 1 0]
            [1 0 1]
        """
        return "LinearLayer of dimension %d x %d represented as\n%s" % (self.m, self.n, self._L)

    def matrix(self):
        """
        Returns the matrix representing this linear layer
        """
        return self._L

    def binary_matrix(self):
        """
        Returns the matrix representing this linear layer in it's binary representation
        """
        if self._L.base_ring() is GF(2):
            return self._L
        return _ff_matrix_to_binary(self._L)

    def __eq__(self, other):
        """
        LinearLayers are considered to be equal if they are represented by the
        same matrix.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: L = LinearLayer(matrix(GF(2), [[0,1,0],[1,0,1]]))
            sage: loads(dumps(L)) == L
            True
        """
        return (self._L) == (other._L)

    def __ne__(self, other):
        """
        LinearLayers are considered to be equal if they are represented by the
        same matrix.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: L = LinearLayer(matrix(GF(2), [[0,1,0],[1,0,1]]))
            sage: L != L
            False
        """
        return not self.__eq__(other)

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
            sage: L = LinearLayer(matrix(GF(2), [[0,1,0],[1,0,1]]))
            sage: L(3)
            3

            sage: L((0, 1, 1))
            (1, 1)

            sage: L([0, 1, 1])
            [1, 1]

            sage: L(vector(GF(2), [0, 1, 1]))
            (1, 1)

            sage: L([0]*(L.n+1))
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply LinearLayer to provided element, dimension mismatch.

            sage: L([0, 1/2, 1])
            Traceback (most recent call last):
            ...
            TypeError: Cannot apply LinearLayer to provided element, conversion to vector failed.

            sage: L("failure")
            Traceback (most recent call last):
            ...
            TypeError: Unsupported type for input x: <type 'str'>.
        """
        from sage.modules.free_module_element import FreeModuleElement

        if isinstance(x, integer_types + (Integer,)):
            x = vector(GF(2), ZZ(x).digits(base=2, padto=self.n))
            return ZZ(list(self.binary_matrix() * x), 2)

        elif isinstance(x, tuple):
            if len(x) != self.n:
                raise TypeError("Cannot apply LinearLayer to provided element, dimension mismatch.")

            try:
                x = vector(self.matrix().base_ring(), x)
            except TypeError, ZeroDivisionError:
                raise TypeError("Cannot apply LinearLayer to provided element, conversion to vector failed.")

            return tuple(self.matrix() * x)

        elif isinstance(x, list):
            if len(x) != self.n:
                raise TypeError("Cannot apply LinearLayer to provided element, dimension mismatch.")

            try:
                x = vector(self.matrix().base_ring(), x)
            except TypeError, ZeroDivisionError:
                raise TypeError("Cannot apply LinearLayer to provided element, conversion to vector failed.")

            return list(self.matrix() * x)

        elif isinstance(x, FreeModuleElement):
            return self.matrix() * x

        else:
            raise TypeError("Unsupported type for input x: %s." % type(x))

    def is_permutation(self):
        r"""
        Check if the linear layer is a permutation, i.e., if the representing matrix is a permutation.

        EXAMPLES::

            sage: from sage.crypto.linearlayer import LinearLayer
            sage: L1 = LinearLayer(identity_matrix(4, 4))
            sage: L1.is_permutation()
            True

            sage: L2 = LinearLayer(matrix([[0,1,1,0], [1,0,0,0], [0,1,0,0], [0,0,0,1]]))
            sage: L2.is_permutation()
            False

            sage: L3 = LinearLayer(matrix([[0,1,0], [1,0,1]]))
            sage: L3.is_permutation()
            False
        """
        # a permutation matrix has to be a square
        if not self._L.is_square():
            return False

        # in each row, all entries execpt one should be 0, the other should be 1
        for r in self._L.rows():
            if r.hamming_weight() != 1:
                return False
            if r[r.nonzero_positions()[0]] != 1:
                return False

        # in each column, all entries execpt one should be 0, the other should be 1
        for c in self._L.columns():
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

            sage: from sage.crypto.linearlayers import PRESENT
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

            sage: from sage.crypto.linearlayers import PRESENT
            sage: PRESENT.linear_branch_number()
            2
        """
        if self.is_permutation():
            return 2
        return _branch_number(self.matrix().transpose())


class AESLikeLinearLayer(LinearLayer):
    r"""
    An AES like linear layer consists of ShiftRows and MixColumns matrices,
    corresponding to the linear layer structure used in the AES.

    Here, ShiftRows is either a Permutation object, or a permutation matrix,
    and MixColumns any arbitrary matrix.

    EXAMPLES::

        sage: from.sage.crypto.linearlayers import linearlayers
        sage: linearlayers['AES']
        AES like LinearLayer of dimension 128 x 128 represented by ShiftRows
        [1, 6, 11, 16, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12]
        and MixColumns
        [    x x + 1     1     1]
        [    1     x x + 1     1]
        [    1     1     x x + 1]
        [x + 1     1     1     x]

    AUTHORS:

    - Friedrich Wiemer (2018-07-02): initial version
    """

    @experimental(25735)
    def __init__(self, *args,  **kwargs):
        """
        Construct an AES linear layer for a given matrix `L` of dimension m x n.

        INPUT:

        - ``ShiftRows`` - a permutation or matrix representing the shift rows part.
        - ``MixColumns`` - a matrix representing the mix columns part.

        """
        if "ShiftRows" in kwargs and "MixColumns" in kwargs:
            SR = kwargs["ShiftRows"]
            if isinstance(SR, Permutation):
                SR = SR.to_matrix()

            MC = kwargs["MixColumns"]
        elif len(args) == 2:
            SR = args[0]
            if isinstance(SR, Permutation):
                SR = SR.to_matrix()

            MC = args[1]
        else:
            raise TypeError("No ShiftRows and MixColumns matrices as arguments provided.")

        # TODO: check if ShiftRows is permutation object
        self._L = _ff_matrix_to_binary(column_linear_layer(MC).matrix() * SR)
        self.m = self._L.nrows()
        self.n = self._L.ncols()
        #additionaly set Mixcolumn and ShiftRows
        self._mc = MC[0]
        self._sr = SR

    def _latex_(self):
        raise NotImplementedError

    def _repr_(self):
        """
        EXAMPLES::

            sage: from sage.crypto.linearlayers import Midori
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
        return "AES like LinearLayer of dimension %d x %d represented by ShiftRows\n%s\nand MixColumns\n%s" \
               % (self.m, self.n, perm_sr, self._mc)

    @cached_method
    def differential_branch_number(self):
        """
        Computes the differential branch number of the MixColumn Matrix

        EXAMPLES::

            sage: from sage.crypto.linearlayers import Midori
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

            sage: from sage.crypto.linearlayers import Midori
            sage: Midori.linear_branch_number()
            4
        """
        #if self.is_permutation():
        #    return 2
        M = self._mc
        M = M.transpose()
        return _branch_number(M)
