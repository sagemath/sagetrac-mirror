r"""
Channels

Given an input space and an output space, a channel takes element from
the input space (the message) and transforms it into an element of the output space
(the transmitted message).

This file contains the following elements:

    - *AbstractChannel*, the abstract class for Channels
    - *StaticErrorRateChannel*, which creates a specific number of errors in each
      transmitted message
    - *ErrorErasureChannel*, which creates a specific number of errors and a
      specific number of erasures in each transmitted message
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

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.rings.finite_rings.constructor import GF
from sage.misc.prandom import randint, random, sample
from sage.modules.free_module_element import vector
from sage.misc.abstract_method import abstract_method
from sage.combinat.cartesian_product import CartesianProduct
from sage.modules.free_module import VectorSpace
from copy import copy

def random_error_vector(n, F, error_positions):
    r"""
    Return a vector of length ``n`` over ``F`` filled with random non-zero coefficients
    at the positions given by ``error_positions``.

    .. NOTE::

        This is a helper function, which should only be used when implementing new channels.

    INPUT:

    - ``n`` -- the length of the vector

    - ``F`` -- the field over which the vector is defined

    - ``error_positions`` -- the non-zero positions of the vector

    OUTPUT:

    - a vector of ``F``

    AUTHORS:

    This function is taken from codinglib (https://bitbucket.org/jsrn/codinglib/)
    and was written by Johan Nielsen.

    EXAMPLES::

        sage: sage.coding.channel_constructions.random_error_vector(5, GF(2), [1,3])
        (0, 1, 0, 1, 0)
    """
    vect = [F.zero()]*n
    for i in error_positions:
        vect[i] = F._random_nonzero_element()
    return vector(F, vect)

class AbstractChannel(SageObject):
    r"""
    Abstract top-class for Channel objects.

    All channel objects must inherit from this class. To implement a channel subclass, one should
    do the following:

    - inherit from this class,

    - call the super constructor,

    - override :meth:`transmit_unsafe`.

    While not being mandatory, it might be useful to reimplement representation methods (``__repr__`` and
    ``_latex_``), along with a comparison method (``__eq__``).

    This abstract class provides the following parameters:

        - ``input_space`` -- the space of the words to transmit

        - ``output_space`` -- the space of the transmitted words
    """

    def __init__(self, input_space, output_space):
        r"""
        Initializes parameters for a Channel object.

        This is a private method, which should be called by the constructor
        of every encoder, as it automatically initializes the mandatory
        parameters of a Channel object.

        INPUT:

        - ``input_space`` -- the space of the words to transmit

        - ``output_space`` -- the space of the transmitted words

        EXAMPLES:

        We first create a new Channel subclass::

            sage: class ChannelExample(sage.coding.channel_constructions.AbstractChannel):
            ....:   def __init__(self, input_space, output_space):
            ....:       super(ChannelExample, self).__init__(input_space, output_space)

        We now create a member of our newly made class::

            sage: input = VectorSpace(GF(7), 6)
            sage: output = VectorSpace(GF(7), 5)
            sage: Chan = ChannelExample(input, output)

        We can check its parameters::

            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 7
            sage: Chan.output_space()
            Vector space of dimension 5 over Finite Field of size 7
        """
        self._input_space = input_space
        self._output_space = output_space

    def transmit(self, message):
        r"""
        Returns ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.
        Checks if ``message`` belongs to the input space, and returns an exception if not.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space of ``self``

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: Chan.transmit(msg) # random
            (4, 14, 15, 16, 17, 42)

        If we transmit a vector which is not in the input space of ``self``::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = (4, 8, 15, 16, 23, 42)
            sage: Chan.transmit(msg)
            Traceback (most recent call last):
            ...
            TypeError: Message must be an element of the input space for the given channel

        .. NOTE::

            One can also call directly ``Chan(message)``, which does the same as ``Chan.transmit(message)``
        """
        if message in self.input_space():
            return self.transmit_unsafe(message)
        else :
            raise TypeError("Message must be an element of the input space for the given channel")

    def __call__(self, message):
        r"""
        This is an alias for ``transmit`` method.
        With this method, one can transmit a message ``msg`` through a Channel ``Chan``
        by calling ``Chan(msg)`` instead of having to write ``Chan.transmit(msg)``.
        See :meth:`transmit` for details.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: set_random_seed(10)
            sage: m1 = Chan(msg)
            sage: set_random_seed(10)
            sage: m2 = Chan.transmit(msg)
            sage: m1 == m2
            True
        """
        return self.transmit(message)

    def input_space(self):
        r"""
        Returns the input space of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 59

        """
        return self._input_space

    def output_space(self):
        r"""
        Returns the output space of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan.output_space()
            Vector space of dimension 6 over Finite Field of size 59
        """
        return self._output_space

    @abstract_method
    def transmit_unsafe(self, message):
        r"""
        Returns ``message``, modified accordingly with the algorithm of the channel it was
        transmitted through.
        This method does not check if ``message`` belongs to the input space of``self``.

        This is an abstract method which should be reimplemented in all the subclasses of
        Channel.
        """
        raise NotImplementedError










class StaticErrorRateChannel(AbstractChannel):
    r"""
    Constructs a channel which adds a static number of errors to each message
    it transmits. The input space and the output space of this channel
    are the same.

    INPUT:

    - ``space`` -- the space of both input and output

    - ``number_errors`` -- the number of errors added to each transmitted message
      It can be either an integer of a tuple. If a tuple is passed as
      argument, the number of errors will be a random integer between the
      two bounds of the tuple.

    EXAMPLES:

    We construct a StaticErrorRateChannel which adds 2 errors
    to any transmitted message::

        sage: F = VectorSpace(GF(59), 40)
        sage: n_err = 2
        sage: Chan = channels.StaticErrorRateChannel(F, n_err)
        sage: Chan
        Static error rate channel creating 2 error(s)

    We can also pass a tuple for the number of errors::

        sage: F = VectorSpace(GF(59), 40)
        sage: n_err = (1, 10)
        sage: Chan = channels.StaticErrorRateChannel(F, n_err)
        sage: Chan
        Static error rate channel creating between 1 and 10 errors
    """

    def __init__(self, space, number_errors):
        r"""
        TESTS:

        If the number of errors exceeds the dimension of the input space,
        it will return an error::

            sage: F = VectorSpace(GF(59), 40)
            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            Traceback (most recent call last):
            ...
            ValueError: There might be more errors than the dimension of the input space
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (number_errors, number_errors)
        if not isinstance(number_errors, (tuple, list)):
            raise ValueError("number_errors must be a tuple, a list, an Integer or a Python int")
        super(StaticErrorRateChannel, self).__init__(space, space)
        if number_errors[1] > space.dimension():
            raise ValueError("There might be more errors than the dimension of the input space")
        self._number_errors = number_errors

    def __repr__(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan
            Static error rate channel creating 42 error(s)
        """
        no_err = self.number_errors()
        if no_err[0] == no_err[1]:
            return "Static error rate channel creating %s error(s)"\
                    % no_err[0]
        else:
            return "Static error rate channel creating between %s and %s errors"\
                    % (no_err[0], no_err[1])

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan._latex_()
            '\\textnormal{Static error rate channel, creating }42 \\textnormal{ error(s)}'
        """
        no_err = self.number_errors()
        if no_err[0] == no_err[1]:
            return "\\textnormal{Static error rate channel, creating }%s \\textnormal{ error(s)}"\
                    % no_err[0]
        else:
            return "\\textnormal{Static error rate channel, creating between %s and %s errors}"\
                    % (no_err[0], no_err[1])

    def __eq__(self, other):
        r"""
        Checks if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err = 42
            sage: Chan1 = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan2 = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan1 == Chan2
            True
        """
        return isinstance(other, StaticErrorRateChannel)\
                and self.input_space() == other.input_space()\
                and self.number_errors() == other.number_errors()

    def transmit_unsafe(self, message):
        r"""
        Returns ``message`` with as many errors as ``self._number_errors`` in it.
        If ``self._number_errors`` was passed as a tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.
        This method does not check if ``message`` belongs to the input space of``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: Chan.transmit_unsafe(msg) # random
            (4, 14, 15, 16, 17, 42)
        """
        w = copy(message)
        number_errors = randint(*self.number_errors())
        V = self.input_space()
        for i in sample(xrange(V.dimension()), number_errors):
            w[i] = V.base_ring().random_element()
        return w

    def number_errors(self):
        r"""
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err = 3
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors










class ErrorErasureChannel(AbstractChannel):
    r"""
    Constructs a channel which adds errors to any message it transmits. It also erases several positions
    in the transmitted message. The output space of this channel is a cartesian product
    between its input space and a VectorSpace of the same dimension over GF(2)

    INPUT:

    - ``space`` -- the input and output space

    - ``number_errors`` -- the number of errors created in each transmitted
      message. It can be either an integer of a tuple. If an tuple is passed as
      an argument, the number of errors will be a random integer between the
      two bounds of this tuple.

    - ``number_erasures`` -- the number of erasures created in each transmitted
      message. It can be either an integer of a tuple. If an tuple is passed as an
      argument, the number of erasures will be a random integer between the
      two bounds of this tuple.

    EXAMPLES:

    We construct a ErrorErasureChannel which adds 2 errors
    and 2 erasures to any transmitted message::

        sage: F = VectorSpace(GF(59), 40)
        sage: n_err, n_era = 2, 2
        sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating 2 error(s) and 2 erasure(s)

    We can also pass the number of errors and erasures as a couple of integers::

        sage: F = VectorSpace(GF(59), 40)
        sage: n_err, n_era = (1, 10), (1, 10)
        sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating between 1 and 10 errors and between 1 and 10 erasures
    """

    def __init__(self, space, number_errors, number_erasures):
        r"""


        TESTS:

        If the sum of number of errors and number of erasures
        exceeds (or may exceed, in the case of tuples) the dimension of the input space,
        it will return an error::

            sage: F = VectorSpace(GF(59), 40)
            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            Traceback (most recent call last):
            ...
            ValueError: The total number of errors and erasures can exceed the dimension of the input space
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (number_errors, number_errors)
        if not isinstance(number_errors, (tuple, list)):
            raise ValueError("number_errors must be a tuple, a list, an Integer or a Python int")

        if isinstance(number_erasures, (Integer, int)):
            number_erasures = (number_erasures, number_erasures)
        if not isinstance(number_erasures, (tuple, list)):
            raise ValueError("number_erasures must be a tuple, a list, an Integer or a Python int")

        output_space = CartesianProduct(space, VectorSpace(GF(2), space.dimension()))
        super(ErrorErasureChannel, self).__init__(space, output_space)
        if number_errors[1] + number_erasures[1] > space.dimension():
            raise ValueError("The total number of errors and erasures can exceed the dimension of the input space")
        self._number_errors = number_errors
        self._number_erasures = number_erasures

    def __repr__(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: Chan
            Error-and-erasure channel creating 21 error(s) and 21 erasure(s)
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        if no_err[0] == no_err[1] and no_era[0] == no_era[1]:
            return "Error-and-erasure channel creating %s error(s) and %s erasure(s)"\
                    %(no_err[0], no_era[0])
        elif no_err[0] != no_err[1] and no_era[0] == no_era[1]:
            return "Error-and-erasure channel creating between %s and %s errors and %s erasure(s)"\
                    % (no_err[0], no_err[1], no_era[0])
        elif no_err[0] == no_err[1] and no_era[0] != no_era[1]:
            return "Error-and-erasure channel creating %s error(s) and between %s and %s erasures"\
                    % (no_err[0], no_era[0], no_era[1])
        else:
            return "Error-and-erasure channel creating between %s and %s errors and between %s and %s erasures"\
                    % (no_err[0], no_err[1], no_era[0], no_era[1])

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: latex(Chan)
            \textnormal{Error-and-erasure channel creating 21 error(s) and 21 erasure(s)}
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        if no_err[0] == no_err[1] and no_era[0] == no_era[1]:
            return "\\textnormal{Error-and-erasure channel creating %s error(s) and %s erasure(s)}"\
                    %(no_err[0], no_era[0])
        elif no_err[0] != no_err[1] and no_era[0] == no_era[1]:
            return "\\textnormal{Error-and-erasure channel creating between %s and %s error(s) and %s erasure(s)}"\
                    % (no_err[0], no_err[1], no_era[0])
        elif no_err[0] == no_err[1] and no_era[0] != no_era[1]:
            return "\\textnormal{Error-and-erasure channel creating %s error(s) and between %s and %s erasure(s)}"\
                    % (no_err[0], no_era[0], no_era[1])
        else:
            return "\\textnormal{Error-and-erasure channel creating between %s and %s error(s) and between %s and %s erasure(s)}"\
                    % (no_err[0], no_err[1], no_era[0], no_era[1])

    def __eq__(self, other):
        r"""
        Checks if ``self`` is equal to ``other``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 50)
            sage: n_err = 17
            sage: n_era = 7
            sage: Chan1 = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: Chan2 = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: Chan1 == Chan2
            True
        """
        return isinstance(other, ErrorErasureChannel)\
                and self.input_space() == other.input_space()\
                and self.number_errors() == other.number_errors()\
                and self.number_erasures() == other.number_erasures()

    def transmit_unsafe(self, message):
        r"""
        Returns ``message`` with as many errors as ``self._number_errors`` in it, and as many erasures
        as ``self._number_erasures`` in it.
        If ``self._number_errors`` was passed as an tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.
        The same method applies with ``self._number_erasures``.

        All erased positions are set to 0 in the transmitted message.
        It is guaranteed that the erasures and the errors will never overlap:
        the received message will always contains exactly as many errors and erasures
        as expected.

        This method does not check if ``message`` belongs to the input space of``self``.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a couple of vectors, namely:

            - the transmitted message, which is ``message`` with erroneous and erased positions
            - the erasure vector, which contains ``1`` at the erased positions of the transmitted message
              , 0 elsewhere.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 11)
            sage: n_err, n_era = 2, 2
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: msg = F((3, 14, 15, 9, 26, 53, 58, 9, 7, 9, 3))
            sage: Chan.transmit_unsafe(msg) # random
            ((0, 14, 15, 0, 26, 53, 45, 9, 7, 14, 3), (1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0))
        """
        number_errors = randint(*self.number_errors())
        number_erasures = randint(*self.number_erasures())
        V = self.input_space()
        n = V.dimension()

        erroneous_positions = sample(xrange(n),\
                number_errors + number_erasures)
        error_split = sample(xrange(number_errors + number_erasures),\
                number_errors)
        error_positions = [erroneous_positions[i] for i in\
                range(number_errors + number_erasures) if i in error_split]
        erasure_positions = [erroneous_positions[i] for i in\
                range(number_errors + number_erasures) if i not in error_split]

        error_vector = random_error_vector(n, V.base_ring(), error_positions)
        erasure_vector = random_error_vector(n , GF(2), erasure_positions)

        message = message + error_vector

        for i in erasure_positions:
            message[i] = 0
        return message, erasure_vector

    def number_errors(self):
        r"""
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err, n_era = 3, 0
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors

    def number_erasures(self):
        r"""
        Returns the number of erasures created by ``self``.

        EXAMPLES::

            sage: F = VectorSpace(GF(59), 6)
            sage: n_err, n_era = 0, 3
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: Chan.number_erasures()
            (3, 3)
        """
        return self._number_erasures
