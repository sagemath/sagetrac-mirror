r"""
Punctured code

Let `C` be a linear code. Let `C_i` be the set of all words of `C` with the
`i`-th coordinate being removed. `C_i` is the punctured code of `C`
on the `i`-th position.
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

from linear_code import AbstractLinearCode
from encoder import Encoder
from decoder import Decoder, DecodingError
from sage.misc.cachefunc import cached_method
from sage.rings.integer import Integer
from sage.modules.free_module import VectorSpace
from sage.rings.finite_rings.constructor import GF
from copy import copy

def puncture(v, points, code):
    r"""
    Returns v punctured as the positions listed in ``points``.

    INPUT:

    - ``v`` -- a vector

    - ``points`` -- a list of integers

    - ``code`` -- the code in which ``v`` lives

    EXAMPLES::

        sage: C = codes.RandomLinearCode(11, 5, GF(7))
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: v = vector(GF(7), (2,3,0,2,1,5,1,5,6,5,3))
        sage: sage.coding.punctured_code.puncture(v, Cp.punctured_positions(), Cp)
        (2, 3, 0, 1, 5, 1, 5, 6, 5, 3)
    """
    S = code.ambient_space()
    vl = v.list()
    v_final = []
    start = 0
    for i in points:
        v_final += vl[start:i]
        start = i + 1
    v_final += vl[start:len(vl)]
    return S(v_final)


class PuncturedCode(AbstractLinearCode):
    r"""
    Representation of a punctured code.

    - ``C`` -- A linear code

    - ``positions`` -- the positions where ``C`` will be punctured. It can be either an integer
      if one need to puncture only one position, or a list of positions to puncture.
      If the same position is passed several times, it will be considered only once.

    EXAMPLES::

        sage: C = codes.RandomLinearCode(11, 5, GF(7))
        sage: Cp = codes.PuncturedCode(C, 3)
        sage: Cp
        Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]

        sage: Cp = codes.PuncturedCode(C, [3, 5])
        sage: Cp
        Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3, 5]
    """
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, C, positions):
        r"""
        TESTS:

        If one of the positions to puncture is bigger than the length of ``C``, an exception will be raised::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, [4,8,15])
            Traceback (most recent call last):
            ...
            ValueError: Positions to puncture must be positive integers smaller than the length of the provided code
        """
        if not isinstance(positions, (Integer, int, tuple, list)):
            raise TypeError("positions must be either a Sage Integer, a Python int, a tuple or a list")
        if isinstance(positions, (Integer, int)):
            positions = [positions]
        if not isinstance(C, AbstractLinearCode):
            raise ValueError("Provided code must be a linear code")
        if not all (i in range(0, C.length()) for i in positions):
            raise ValueError("Positions to puncture must be positive integers smaller than the length of the provided code")
        unique_positions = set()
        for i in positions:
            unique_positions.add(i)
        positions = []
        for i in unique_positions:
            positions.append(i)
        super(PuncturedCode, self).__init__(C.base_ring(), C.length() - len(positions), \
                "PuncturedMatrix", "OriginalCode")
        positions.sort()
        self._original_code = C
        self._positions = positions

    def __eq__(self, other):
        r"""
        Tests equality between two Punctured codes.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp1 = codes.PuncturedCode(C, 2)
            sage: Cp2 = codes.PuncturedCode(C, 2)
            sage: Cp1 == Cp2
            True
        """
        return isinstance(other, PuncturedCode) \
                and self.punctured_positions() == other.punctured_positions() \
                and self.original_code() == other.original_code()

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp
            Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        return "Punctured code coming from %s punctured on position(s) %s"\
                % (self.original_code(), self.punctured_positions())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: latex(Cp)
            \textnormal{Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) } [3]
        """
        return "\\textnormal{Punctured code coming from %s punctured on position(s) } %s"\
                % (self.original_code(), self.punctured_positions())

    def punctured_positions(self):
        r"""
        Returns the list of positions which were punctured on the original code.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.punctured_positions()
            [3]
        """
        return self._positions

    def original_code(self):
        r"""
        Returns the linear code which was punctured to get ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.original_code()
            Linear code of length 11, dimension 5 over Finite Field of size 7
        """
        return self._original_code

    def dimension(self):
        r"""
        Returns the dimension of ``self``.

        EXAMPLES::

            sage: set_random_seed(42)
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: Cp.dimension()
            5
        """
        if hasattr(self, '_dimension'):
            return self._dimension
        self._dimension = self.generator_matrix().rank()
        return self._dimension

    def random_element(self, *args, **kwds):
        r"""
        Returns a random codeword of ``self``.

        This methods does not trigger the computation of
        ``self``'s :meth:`sage.coding.linear_code.generator_matrix`.

        INPUT:

        - ``agrs``, ``kwds`` - extra positional arguments passed to
          :meth:`sage.modules.free_module.random_element`.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: set_random_seed(10)
            sage: Cp.random_element()
            (2, 0, 1, 3, 3, 3, 2, 6, 0, 5)
        """
        C_original = self.original_code()
        m = (C_original.base_ring() ** C_original.dimension()).random_element()
        c = C_original.encode(m)
        return puncture(c, self.punctured_positions(), self)

    def encode(self, m, original_encode=False, encoder_name=None, **kwargs):
        r"""
        Transforms an element of the message space into an element of the code.

        INPUT:

        - ``m`` -- a vector of the message space of the code.

        - ``original_encode`` -- (default: ``False``) if this is set to ``True``,
          ``m`` will be encoded using an Encoder of ``self``'s :meth:`original_code`.
          This allow to avoid the computation of a generator matrix for ``self``.

        - ``encoder_name`` -- (default: ``None``) Name of the encoder which will be used
          to encode ``word``. The default encoder of ``self`` will be used if
          default value is kept

        OUTPUT:

        - an element of ``self``

        EXAMPLES::

           sage: M = matrix(GF(7), [[1, 0, 0, 0, 3, 4, 6], [0, 1, 0, 6, 1, 6, 4], [0, 0, 1, 5, 2, 2, 4]])
           sage: C_original = LinearCode(M)
           sage: Cp = codes.PuncturedCode(C_original, 2)
           sage: m = vector(GF(7), [1, 3, 5])
           sage: Cp.encode(m)
            (1, 3, 5, 5, 0, 2)
        """
        if original_encode:
            c = self.original_code().encode(m, encoder_name, **kwargs)
            return puncture(c, self.punctured_positions, self)
        return self.encoder(encoder_name, **kwargs).encode(m)

    @cached_method
    def structured_representation(self):
        r"""
        Returns ``self`` as a structured code object.

        If ``self`` has a specific structured representation (e.g. a punctured GRS code is
        a GRS code too), it will return this representation, else it returns a
        :class:`sage.coding.linear_code.LinearCode`.

        EXAMPLES:

        We consider a GRS code::

            sage: C_grs = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)

        A punctured GRS code is still a punctured code::

            sage: Cp_grs = codes.PuncturedCode(C_grs, 3)
            sage: Cp_grs.structured_representation()
            [39, 12, 28] Generalized Reed-Solomon Code over Finite Field of size 59

        Which is not the case for generic linear codes::

            sage: set_random_seed(42)
            sage: C_lin  = codes.RandomLinearCode(10, 5, GF(2))
            sage: Cp_lin = codes.PuncturedCode(C_lin, 2)
            sage: Cp_lin.structured_representation()
            Linear code of length 9, dimension 5 over Finite Field of size 2
        """
        C = self.original_code()
        pts = copy(self.punctured_positions())
        cpt = 1
        while(isinstance(C, PuncturedCode)):
            cur_pts = C.punctured_positions()
            for i in cur_pts:
                pts.append(i + cpt)
                cpt += 1
            C = C.original_code()
        return C._punctured_form(pts)








class PuncturedCodePuncturedMatrixEncoder(Encoder):
    r"""
    Encoder using original code generator matrix to compute the punctured code's one.

    INPUT:

    - ``code`` -- The associated code of this encoder.

        EXAMPLES::
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
    """

    def __init__(self, code):
        r"""
        EXAMPLES::
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        super(PuncturedCodePuncturedMatrixEncoder, self).__init__(code)

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E
            Punctured matrix-based encoder for the Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3]
        """
        return "Punctured matrix-based encoder for the %s" % self.code()

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: latex(E)
            \textnormal{Punctured matrix-based encoder for the }\textnormal{Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) } [3]
        """
        return "\\textnormal{Punctured matrix-based encoder for the }%s" % self.code()._latex_()

    @cached_method
    def generator_matrix(self):
        r"""
        Returns a generator matrix of the associated code of ``self``.

        EXAMPLES::

            sage: set_random_seed(10)
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: E = codes.encoders.PuncturedCodePuncturedMatrixEncoder(Cp)
            sage: E.generator_matrix()
            [1 0 0 0 0 5 2 6 0 6]
            [0 1 0 0 0 5 2 2 1 1]
            [0 0 1 0 0 6 2 4 0 4]
            [0 0 0 1 0 0 6 3 3 3]
            [0 0 0 0 1 0 1 3 4 3]
        """
        C = self.code().original_code()
        pos = self.code().punctured_positions()
        M = C.generator_matrix()
        G = M.delete_columns(pos)
        G = G.echelon_form()
        delete = []
        cpt = 0
        for i in G.rows():
            if i.is_zero():
                delete.append(cpt)
            cpt += 1
        return G.delete_rows(delete)









class PuncturedCodeOriginalCodeDecoder(Decoder):
    r"""
    Decoder decoding through a decoder over the original code of its punctured code.

    INPUT:

    - ``code`` -- The associated code of this encoder

    - ``strategy`` -- (dafault: ``None``) the strategy used to decode.
      The available strategies are:

        * ``'error-erasure'`` -- uses an error-erasure decoder over the original code if available,
           fails otherwise.

        * ``'random-values'`` -- fills the punctured positions with random elements
           in ``code``'s base field and tries to decode using
           the default decoder of the original code

        * ``'try-all'`` -- fills the punctured positions with every possible combination of
           symbols until decoding succeeds, or until every combination have been tried

        * ``None`` -- uses ``error-erasure`` if an error-erasure decoder is available,
           switch to ``random-values`` behaviour otherwise

    - ``original_decoder`` -- (default: ``None``) the decoder that will be used over the original code.
      It has to be a decoder object over the original code.
      This argument takes precedence over ``strategy``: if both ``original_decoder`` and ``strategy``
      are filled, ``self`` will use the ``original_decoder`` to decode over the original code.
      If ``original_decoder`` is set to ``None``, it will use the decoder picked by ``strategy``.

    - ``**kwargs`` -- all extra arguments are forwarded to original code's decoder

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D
            Decoder of Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3] through Syndrome decoder for Linear code of length 11, dimension 5 over Finite Field of size 7
    """

    def __init__(self, code, strategy = None, original_decoder = None, **kwargs):
        r"""
        TESTS:

        If one tries to use ``'error-erasure'`` strategy when the original code has no such
        decoder, it returns an error::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp, strategy = 'error-erasure')
            Traceback (most recent call last):
            ...
            ValueError: Original code has no error-erasure decoder
        """
        original_code = code.original_code()
        if original_decoder is not None:
            if not isinstance(original_decoder, Decoder):
                raise TypeError("original_decoder must be a decoder object")
            if not original_decoder.code() == original_code:
                raise ValueError("original decoder must have the original code as associated code")
            if 'error-erasure' in original_decoder.decoder_type():
                strategy = 'error-erasure'
            self._original_decoder = original_decoder
        elif strategy == 'error-erasure':
            error_erasure = 0
            for D in original_code._registered_decoders.values():
                if 'error-erasure' in D._decoder_type:
                    error_erasure = 1
                    self._original_decoder = D(original_code, **kwargs)
                    break
            if not error_erasure:
                raise ValueError("Original code has no error-erasure decoder")
        elif strategy == 'random-values' or strategy == 'try-all':
            self._original_decoder = code.original_code().decoder(**kwargs)
        else:
            error_erasure = 0
            for D in original_code._registered_decoders.values():
                if 'error-erasure' in D._decoder_type:
                    error_erasure = 1
                    self._original_decoder = D(original_code, **kwargs)
                    break
            if not error_erasure:
                self._original_decoder = original_code.decoder(**kwargs)
        self._strategy = strategy
        self._decoder_type = copy(self._decoder_type)
        self._decoder_type.remove("dynamic")
        self._decoder_type = self._original_decoder.decoder_type()
        super(PuncturedCodeOriginalCodeDecoder, self).__init__(code, code.ambient_space(),\
                self._original_decoder.connected_encoder())

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D
            Decoder of Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3] through Syndrome decoder for Linear code of length 11, dimension 5 over Finite Field of size 7
        """
        return "Decoder of %s through %s" % (self.code(), self.original_decoder())


    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: latex(D)
            \textnormal{Decoder of } Punctured code coming from Linear code of length 11, dimension 5 over Finite Field of size 7 punctured on position(s) [3] \textnormal{ through } Syndrome decoder for Linear code of length 11, dimension 5 over Finite Field of size 7
        """
        return "\\textnormal{Decoder of } %s \\textnormal{ through } %s" % (self.code(), self.original_decoder())

    def original_decoder(self):
        r"""
        Returns the decoder over the original code that will be used to decode words of
        :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D.original_decoder()
            Syndrome decoder for Linear code of length 11, dimension 5 over Finite Field of size 7
        """
        return self._original_decoder

    def decode_to_code(self, y):
        r"""
        Decodes ``y`` to an element in :meth:`sage.coding.decoder.Decoder.code`.

        EXAMPLES::

            sage: set_random_seed(42)
            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: c = vector(GF(7), [6, 5, 5, 6, 2, 3, 4, 4, 1, 5])
            sage: c in Cp
            True
            sage: y = vector(GF(7), [6, 5, 5, 6, 0, 3, 4, 4, 1, 5])
            sage: y in Cp
            False
            sage: D.decode_to_code(y) == c
            True
        """
        D = self.original_decoder()
        C = self.code()
        A = C.original_code().ambient_space()
        Cor = C.original_code()
        pts = C.punctured_positions()
        F = self.code().base_field()
        zero, one = F.zero(), F.one()
        if "error-erasure" in D.decoder_type():
            if isinstance(y, (tuple, list)):
                y, e = y[0], y[1]
                e_list = e.list()
                shift = 0
                for i in range(Cor.length()):
                    if i in pts:
                        e_list.insert(i + shift, one)
            else:
                e_list = []
                for i in range(Cor.length()):
                    if i in pts:
                        e_list.append(one)
                    else:
                        e_list.append(zero)
            e = vector(GF(2), e_list)
            yl = y.list()
            shift = 0
            for i in pts:
                yl.insert(i + shift, F.zero())
                shift += 1
            y = A(yl)
            return puncture(D.decode_to_code((y, e)), pts, C)
        elif self._strategy == 'try-all':
            end = False
            yl = y.list()
            I = iter(VectorSpace(F, len(pts)))
            shift = 0
            for i in pts:
                yl.insert(i + shift, F.zero())
                shift += 1
            values = I.next()
            while not end:
                try:
                    shift = 0
                    for i in pts:
                        yl[i + shift] =  values[shift]
                        shift += 1
                    y = A(yl)
                    values = I.next()
                    try:
                        c_or = self.original_decoder().decode_to_code(y)
                        end = True
                        break
                    except:
                        pass
                except StopIteration:
                    raise DecodingError
            return puncture(c_or, pts, C)
        A = C.original_code().ambient_space()
        yl = y.list()
        shift = 0
        for i in pts:
            yl.insert(i + shift, F.random_element())
            shift += 1
        y = A(yl)
        return puncture(self.original_decoder().decode_to_code(y), pts, C)

    def decoding_radius(self, number_erasures = None):
        r"""
        Returns maximal number of errors that ``self`` can decode.

        EXAMPLES::

            sage: C = codes.RandomLinearCode(11, 5, GF(7))
            sage: Cp = codes.PuncturedCode(C, 3)
            sage: D = codes.decoders.PuncturedCodeOriginalCodeDecoder(Cp)
            sage: D.decoding_radius()
            1
        """
        punctured = len(self.code().punctured_positions())
        D = self.original_decoder()
        if self._strategy != 'try-all' and "error-erasure" not in D.decoder_type():
            return D.decoding_radius() - punctured
        elif "error-erasure" in D.decoder_type() and number_erasures is not None:
            diff = self.code().original_code().minimum_distance() - number_erasures
            if diff <= 0:
                raise ValueError("The number of erasures exceed decoding capability")
            return (diff - punctured) // 2
        elif "error-erasure" in D.decoder_type() and number_erasures is None:
            raise ValueError("You must provide the number of erasures")
        return D.decoding_radius()

####################### registration ###############################

PuncturedCode._registered_encoders["PuncturedMatrix"] = PuncturedCodePuncturedMatrixEncoder
PuncturedCode._registered_decoders["OriginalCode"] = PuncturedCodeOriginalCodeDecoder
PuncturedCodeOriginalCodeDecoder._decoder_type = {"dynamic"}
