r"""
Decoder

Representation of an error-correction algorithm for a code.

AUTHORS:

- David Joyner (2009-02-01): initial version
- David Lucas (2015-06-29): abstract class version

"""
#*****************************************************************************
#       Copyright (C) 2009 David Joyner <wdjoyner@gmail.com>
#                     2015 David Lucas <david.lucas@inria.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or later (at your preference).
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.abstract_method import abstract_method
from sage.structure.sage_object import SageObject

class Decoder(SageObject):
    r"""
    Abstract top-class for :class:`Decoder` objects.

    Every decoder class should inherit from this abstract class.

    To implement an decoder, you need to:

    - inherit from :class:`Decoder`

    - call ``Decoder.__init__`` in the subclass constructor.
      Example: ``super(SubclassName, self).__init__(code, input_space,
      connected_encoder_name)``.
      By doing that, your subclass will have all the parameters described above initialized.

    - Then, you need to override one of decoding methods, either :meth:`decode_to_code` or
      :meth:`decode_to_message`. You can also override the optional method :meth:`decoding_radius`.

    - By default, comparison of :class:`Decoder` (using methods ``__eq__`` and ``__ne__`` ) are
      by memory reference: if you build the same decoder twice, they will be different. If you
      need something more clever, override ``__eq__`` and ``__ne__`` in your subclass.

    - As :class:`Decoder` is not designed to be instantiated, it does not have any representation
      methods. You should implement ``_repr_`` and ``_latex_`` methods in the subclass.
    """

    @classmethod
    def decoder_type(cls):
        r"""
        Returns the set of types of ``self``. These types describe the nature of ``self``
        and its decoding algorithm.

        This method can be called on both an uninstantiated decoder class,
        or on an instance of a decoder class.

        EXAMPLES:

        We call it on a class::

            sage: codes.decoders.LinearCodeSyndromeDecoder.decoder_type()
            {'dynamic', 'hard-decision', 'unique'}

        We can also call it on a instance of a Decoder class::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.decoder_type()
            {'complete', 'hard-decision', 'might-error', 'unique'}
        """
        return cls._decoder_type

    def _instance_decoder_type(self):
        r"""
        Returns the set if types of ``self``.
        These types describe the nature of ``self`` and its decoding algorithm.

        This method is used as a copy of :meth:`decoder_type` when instantiating a
        class, so the returned set of types will be the one specific to the instance,
        and no longer the one for the uninstantiated class.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D._instance_decoder_type()
            {'complete', 'hard-decision', 'might-error', 'unique'}
        """
        return self._decoder_type

    def __init__(self, code, input_space, connected_encoder_name):
        r"""
        Initializes mandatory parameters for :class:`Decoder` objects.

        This method only exists for inheritance purposes as it initializes
        parameters that need to be known by every decoder. An abstract
        decoder object should never be created.

        INPUT:

        - ``code`` -- the associated code of ``self``

        - ``input_space`` -- the input space of ``self``, which is the ambient space
          of ``self``'s ``code``

        - ``connected_encoder_name`` -- the associated encoder, which will be
          used by ``self`` to recover elements from the message space

        EXAMPLES:

        We first create a new :class:`Decoder` subclass::

            sage: from sage.coding.decoder import Decoder
            sage: class DecoderExample(Decoder):
            ....:   def __init__(self, code):
            ....:       in_space = code.ambient_space()
            ....:       connected_enc = "GeneratorMatrix"
            ....:       super(DecoderExample, self).__init__(code, in_space, connected_enc)

        We now create a member of our brand new class::

            sage: G = Matrix(GF(2), [[1, 0, 0, 1], [0, 1, 1, 1]])
            sage: C = LinearCode(G)
            sage: D = DecoderExample(C)

        We can check its parameters::

            sage: D.input_space()
            Vector space of dimension 4 over Finite Field of size 2
            sage: D.connected_encoder()
            Generator matrix-based encoder for [4, 2] linear code over GF(2)
            sage: D.code()
            [4, 2] linear code over GF(2)
        """
        self.decoder_type = self._instance_decoder_type
        self._code = code
        self._input_space = input_space
        self._connected_encoder_name = connected_encoder_name

    def __hash__(self):
        r"""
        Returns the hash value of ``self``.

        This is a generic implementation which should be overwritten on decoders
        with extra arguments.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: hash(D) #random
            7575380076354998465
        """
        C = self.code()
        Str = str(C)
        return hash((C, Str)) ^ hash(Str) ^ hash(C)

    def __ne__(self, other):
        r"""
        Tests inequality of ``self`` and ``other``.

        This is a generic implementation, which returns the inverse of ``__eq__`` for self.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: D1 = LinearCode(G).decoder()
            sage: D2 = LinearCode(G).decoder()
            sage: D1 != D2
            False
            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,1,1]])
            sage: D2 = LinearCode(G).decoder()
            sage: D1 != D2
            True
        """
        return not self == other



    def decode_to_code(self, r):
        r"""
        Corrects the errors in ``r`` and returns a codeword.

        This is a default implementation which assumes that the method
        :meth:`decode_to_message` has been implemented, else it returns an exception.

        INPUT:

        - ``r`` -- a element of the input space of ``self``.

        OUTPUT:

        - a vector of :meth:`code`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: word in C
            True
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: w_err in C
            False
            sage: D = C.decoder()
            sage: D.decode_to_code(w_err)
            (1, 1, 0, 0, 1, 1, 0)
        """
        if hasattr(self, "defaulting_decode_to_message"):
            raise NotImplementedError
        else:
            word = self.decode_to_message(r)
            return self.connected_encoder().encode(word)

    def connected_encoder(self):
        r"""
        Returns the connected encoder of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.connected_encoder()
            Generator matrix-based encoder for [7, 4] linear code over GF(2)
        """
        return self.code().encoder(encoder_name=self._connected_encoder_name)

    def decode_to_message(self, r):
        r"""
        Decodes ``r`` to the message space of meth:`connected_encoder`.

        This is a default implementation, which assumes that the method
        :meth:`decode_to_code` has been implemented, else it returns an exception.

        INPUT:

        - ``r`` -- a element of the input space of ``self``.

        OUTPUT:

        - a vector of :meth:`message_space`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: word = vector(GF(2), (1, 1, 0, 0, 1, 1, 0))
            sage: w_err = word + vector(GF(2), (1, 0, 0, 0, 0, 0, 0))
            sage: D = C.decoder()
            sage: D.decode_to_message(w_err)
            (0, 1, 1, 0)
        """
        self.defaulting_decode_to_message = True
        return self.code().unencode(self.decode_to_code(r))

    def code(self):
        r"""
        Returns the code for this :class:`Decoder`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.code()
            [7, 4] linear code over GF(2)
        """
        return self._code

    def message_space(self):
        r"""
        Returns the message space of ``self``'s :meth:`connected_encoder`.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.message_space()
            Vector space of dimension 4 over Finite Field of size 2
        """
        return self.connected_encoder().message_space()

    def input_space(self):
        r"""
        Returns the input space of ``self``.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = C.decoder()
            sage: D.input_space()
            Vector space of dimension 7 over Finite Field of size 2
        """
        if hasattr(self, "_input_space"):
            return self._input_space
        else:
            raise NotImplementedError("Decoder does not have an _input_space parameter")

    @abstract_method(optional = True)
    def decoding_radius(self, **kwargs):
        r"""
        Returns the maximal number of errors that ``self`` is able to correct.

        This is an abstract method and it should be implemented in subclasses.

        EXAMPLES::

            sage: G = Matrix(GF(2), [[1,1,1,0,0,0,0],[1,0,0,1,1,0,0],[0,1,0,1,0,1,0],[1,1,0,1,0,0,1]])
            sage: C = LinearCode(G)
            sage: D = codes.decoders.LinearCodeSyndromeDecoder(C)
            sage: D.decoding_radius()
            1
        """
        raise NotImplementedError

class DecodingError(Exception):
    r"""
    Special exception class to indicate an error during decoding.
    """
    pass
