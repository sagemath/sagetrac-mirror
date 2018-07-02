r"""
LinearLayers used in cryptographic schemes.

All linear layers are available in the dictionary linearlayers, or as a seperate object.
This module provides the following linear layers:

    - AES ([DR2002]_)
    - Midori ([BBISHAR2015]_)
    - SKINNY ([BJKLMPSSS2016]_)
    - PRESENT (and SmallScalePRESENT) ([BKLPPRSV2007]_)

EXAMPLES::

    sage: from sage.crypto.linearlayers import PRESENT
    doctest:warning
    ...
    FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
    See https://trac.sagemath.org/25735 for details.
    sage: PRESENT
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

    - Friedrich Wiemer (2018-04-06): initial version
"""

from sage.crypto.linearlayer import LinearLayer
from sage.crypto.linearlayer import AESLikeLinearLayer

from sage.combinat.permutation import Permutation
from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

Left_ShiftRows = Permutation([1, 6, 11, 16, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12])
Right_ShiftRows = Permutation([1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3, 16, 13, 10, 7, 4])

_AES_irreducible_polynomial = PolynomialRing(GF(2), name="alpha")("alpha^8 + alpha^4 + alpha^3 + alpha + 1")
_AES_field = GF(2**8, name="x", modulus=_AES_irreducible_polynomial)
AES_ShiftRows = Left_ShiftRows
AES_MixColumns = [matrix(_AES_field, 4, 4, map(_AES_field.fetch_int, [2, 3, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 3, 1, 1, 2]))]*4
AES = AESLikeLinearLayer(AES_ShiftRows, AES_MixColumns)

#Midori_ShuffelCells = matrix(GF(2**4), Permutation([1, 6, 16, 11, 14, 9, 3, 8, 12, 15, 5, 2, 7, 4, 10, 13]).to_matrix())
Midori_ShuffelCells = Permutation([1, 11, 6, 16, 15, 5, 12, 2, 10, 4, 13, 7, 8, 14, 3, 9])
Midori_MixColumns = [matrix(GF(2**4), [[0, 1, 1, 1] ,[1, 0, 1, 1] ,[1, 1, 0, 1] ,[1, 1, 1, 0]])]*4
Midori = AESLikeLinearLayer(Midori_ShuffelCells, Midori_MixColumns)

SKINNY_ShiftRows = Right_ShiftRows
SKINNY_4_MixColumns = [matrix(GF(2**4), [[1, 0, 1, 1] ,[1, 0, 0, 0] ,[0, 1, 1, 0] ,[1, 0, 1, 0]])]*4
SKINNY_8_MixColumns = [matrix(GF(2**8), [[1, 0, 1, 1] ,[1, 0, 0, 0] ,[0, 1, 1, 0] ,[1, 0, 1, 0]])]*4
SKINNY_4 = AESLikeLinearLayer(SKINNY_ShiftRows, SKINNY_4_MixColumns)
SKINNY_8 = AESLikeLinearLayer(SKINNY_ShiftRows, SKINNY_8_MixColumns)

def smallscale_present_linearlayer(nsboxes=16):
    """
    The matrix representing SmallPRESENT with nsboxes many S-boxes.

    Original PRESENT corresponds to nsboxes=16.

    INPUT:

    - ``nsboxes`` - integer, number of sboxes the linear layer operates on
      (default: 16).

    TESTS:

        sage: from sage.crypto.linearlayers import smallscale_present_linearlayer
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

    return LinearLayer(matrix(GF(2), [present_llayer(nsboxes, ei)
                                      for ei in VectorSpace(GF(2), 4*nsboxes).basis()
                                     ]))

PRESENT = smallscale_present_linearlayer(nsboxes=16)

# Dictionary of all available linear layers
linearlayers = {}
import sys
for k in dir(sys.modules[__name__]):
    v = getattr(sys.modules[__name__], k)
    if isinstance(v, (LinearLayer, AESLikeLinearLayer)):
        linearlayers[k] = v

