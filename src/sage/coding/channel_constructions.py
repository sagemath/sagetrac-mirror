r"""
Channels

Given an input space and an output space, a channel takes element from
the input space (the message) and transforms it into an element of the output space
(the transmitted message).

This file contains the following elements:

    - :class:`Channel`, the abstract class for Channels
    - :class:`StaticErrorRateChannel`, which creates a specific number of errors in each
      transmitted message
    - :class:`ErrorErasureChannel`, which creates a specific number of errors and a
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
from sage.categories.cartesian_product import cartesian_product
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

        sage: from sage.coding.channel_constructions import random_error_vector
        sage: random_error_vector(5, GF(2), [1,3])
        (0, 1, 0, 1, 0)
    """
    vect = [F.zero()]*n
    for i in error_positions:
        vect[i] = F._random_nonzero_element()
    return vector(F, vect)

def format_interval(t):
    r"""
    Returns a formatted string representation of ``t``.

    This method should be called by any representation function in Channel classes.

    .. NOTE::

        This is a helper function, which should only be used when implementing new channels.

    INPUT:

    - ``t`` -- a list or a tuple

    OUTPUT:

    - a string

    TESTS::

        sage: from sage.coding.channel_constructions import format_interval
        sage: t = (5, 5)
        sage: format_interval(t)
        '5'

        sage: t = (2, 10)
        sage: format_interval(t)
        'between 2 and 10'

    """
    return str(t[0]) if t[0] == t[1] else 'between %s and %s' % (t[0], t[1])

class Channel(SageObject):
    r"""
    Abstract top-class for Channel objects.

    All channel objects must inherit from this class. To implement a channel subclass, one should
    do the following:

    - inherit from this class,

    - call the super constructor,

    - override :meth:`transmit_unsafe`.

    While not being mandatory, it might be useful to reimplement representation methods (``_repr_`` and
    ``_latex_``).

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

            sage: from sage.coding.channel_constructions import Channel
            sage: class ChannelExample(Channel):
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
        Note that ``message`` itself is never modified by the channel.

        INPUT:

        - ``message`` -- a vector

        OUTPUT:

        - a vector of the output space of ``self``

        EXAMPLES::

            sage: F = GF(59)^6
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: set_random_seed(10)
            sage: Chan.transmit(msg)
            (4, 8, 4, 16, 23, 53)

        We can check that the input ``msg`` is not modified::

            sage: msg
            (4, 8, 15, 16, 23, 42)

        If we transmit a vector which is not in the input space of ``self``::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
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

    #Alias for transmit method
    __call__ = transmit

    def input_space(self):
        r"""
        Returns the input space of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: Chan.input_space()
            Vector space of dimension 6 over Finite Field of size 59

        """
        return self._input_space

    def output_space(self):
        r"""
        Returns the output space of ``self``.

        EXAMPLES::

            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
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










class StaticErrorRateChannel(Channel):
    r"""
    Channel which adds a static number of errors to each message it transmits.

    The input space and the output space of this channel are the same.

    The main purpose of communication channels is to transmit messages, which can be achieved with
    two methods:

    - with the method :meth:`Channel.transmit`. Considering a channel ``Chan``
      and a message ``msg``, transmitting
      ``msg`` with ``Chan`` can be done this way::

        Chan.transmit(msg)

      It can also be written in a more convenient way::

        Chan(msg)

    - with the method :meth:`transmit_unsafe`. It does the exact same thing as :meth:`transmit` except
      that it does not check if ``msg`` belongs to the input space of ``Chan``::

        Chan.transmit_unsafe(msg)

    INPUT:

    - ``space`` -- the space of both input and output

    - ``number_errors`` -- the number of errors added to each transmitted message
      It can be either an integer of a tuple. If a tuple is passed as
      argument, the number of errors will be a random integer between the
      two bounds of the tuple.

    EXAMPLES:

    We construct a StaticErrorRateChannel which adds 2 errors
    to any transmitted message::

        sage: n_err = 2
        sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
        sage: Chan
        Static error rate channel creating 2 errors, of input and output space
        Vector space of dimension 40 over Finite Field of size 59

    We can also pass a tuple for the number of errors::

        sage: n_err = (1, 10)
        sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
        sage: Chan
        Static error rate channel creating between 1 and 10 errors,
        of input and output space Vector space of dimension 40 over Finite Field of size 59
    """

    def __init__(self, space, number_errors):
        r"""
        TESTS:

        If the number of errors exceeds the dimension of the input space,
        it will return an error::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^40, n_err)
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

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^50, n_err)
            sage: Chan
            Static error rate channel creating 42 errors, of input and output space
            Vector space of dimension 50 over Finite Field of size 59
        """
        no_err = self.number_errors()
        return "Static error rate channel creating %s errors, of input and output space %s"\
                    % (format_interval(no_err), self.input_space())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: n_err = 42
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^50, n_err)
            sage: latex(Chan)
            \textnormal{Static error rate channel creating 42 errors, of
            input and output space Vector space of dimension 50 over Finite Field of size 59}
        """
        no_err = self.number_errors()
        return "\\textnormal{Static error rate channel creating %s errors, of input and output space %s}"\
                % (format_interval(no_err), self.input_space())

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

            sage: F = GF(59)^6
            sage: n_err = 2
            sage: Chan = channels.StaticErrorRateChannel(F, n_err)
            sage: msg = F((4, 8, 15, 16, 23, 42))
            sage: set_random_seed(10)
            sage: Chan.transmit_unsafe(msg)
            (4, 8, 4, 16, 23, 53)

        This checks that trac #19863 is fixed::

            sage: V = VectorSpace(GF(2), 1000)
            sage: Chan = channels.StaticErrorRateChannel(V, 367)
            sage: c = V.random_element()
            sage: (c - Chan(c)).hamming_weight()
            367
        """
        w = copy(message)
        number_errors = randint(*self.number_errors())
        V = self.input_space()
        R = V.base_ring()
        for i in sample(xrange(V.dimension()), number_errors):
            err = R.random_element()
            while (w[i] == err):
                err = R.random_element()
            w[i] = err
        return w

    def number_errors(self):
        r"""
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err = 3
            sage: Chan = channels.StaticErrorRateChannel(GF(59)^6, n_err)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors










class ErrorErasureChannel(Channel):
    r"""
    Channel which adds errors and erases several positions in any message it transmits.

    The output space of this channel is a Cartesian product
    between its input space and a VectorSpace of the same dimension over GF(2)

    The main purpose of communication channels is to transmit messages, which can be achieved with
    two methods:

    - with the method :meth:`Channel.transmit`. Considering a channel ``Chan``
      and a message ``msg``, transmitting
      ``msg`` with ``Chan`` can be done this way::

        Chan.transmit(msg)

      It can also be written in a more convenient way::

        Chan(msg)

    - with the method :meth:`transmit_unsafe`. It does the exact same thing as :meth:`transmit` except
      that it does not check if ``msg`` belongs to the input space of ``Chan``::

        Chan.transmit_unsafe(msg)

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

        sage: n_err, n_era = 2, 2
        sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating 2 errors and 2 erasures
        of input space Vector space of dimension 40 over Finite Field of size 59
        and output space The Cartesian product of (Vector space of dimension 40
        over Finite Field of size 59, Vector space of dimension 40 over Finite Field of size 2)

    We can also pass the number of errors and erasures as a couple of integers::

        sage: n_err, n_era = (1, 10), (1, 10)
        sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
        sage: Chan
        Error-and-erasure channel creating between 1 and 10 errors and
        between 1 and 10 erasures of input space Vector space of dimension 40
        over Finite Field of size 59 and output space The Cartesian product of
        (Vector space of dimension 40 over Finite Field of size 59,
        Vector space of dimension 40 over Finite Field of size 2)
    """

    def __init__(self, space, number_errors, number_erasures):
        r"""


        TESTS:

        If the sum of number of errors and number of erasures
        exceeds (or may exceed, in the case of tuples) the dimension of the input space,
        it will return an error::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^40, n_err, n_era)
            Traceback (most recent call last):
            ...
            ValueError: The total number of errors and erasures can not exceed the dimension of the input space
        """
        if isinstance(number_errors, (Integer, int)):
            number_errors = (number_errors, number_errors)
        if not isinstance(number_errors, (tuple, list)):
            raise ValueError("number_errors must be a tuple, a list, an Integer or a Python int")

        if isinstance(number_erasures, (Integer, int)):
            number_erasures = (number_erasures, number_erasures)
        if not isinstance(number_erasures, (tuple, list)):
            raise ValueError("number_erasures must be a tuple, a list, an Integer or a Python int")

        output_space = cartesian_product([space, VectorSpace(GF(2), space.dimension())])
        super(ErrorErasureChannel, self).__init__(space, output_space)
        if number_errors[1] + number_erasures[1] > space.dimension():
            raise ValueError("The total number of errors and erasures can not exceed the dimension of the input space")
        self._number_errors = number_errors
        self._number_erasures = number_erasures

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^50, n_err, n_era)
            sage: Chan
            Error-and-erasure channel creating 21 errors and 21 erasures
            of input space Vector space of dimension 50 over Finite Field of size 59
            and output space The Cartesian product of (Vector space of dimension 50
            over Finite Field of size 59, Vector space of dimension 50 over Finite Field of size 2)
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        return "Error-and-erasure channel creating %s errors and %s erasures of input space %s and output space %s"\
                % (format_interval(no_err), format_interval(no_era), self.input_space(), self.output_space())

    def _latex_(self):
        r"""
        Returns a latex representation of ``self``.

        EXAMPLES::

            sage: n_err, n_era = 21, 21
            sage: Chan = channels.ErrorErasureChannel(GF(59)^50, n_err, n_era)
            sage: latex(Chan)
            \textnormal{Error-and-erasure channel creating 21 errors and 21 erasures
            of input space Vector space of dimension 50 over Finite Field of size 59
            and output space The Cartesian product of (Vector space of dimension 50
            over Finite Field of size 59, Vector space of dimension 50 over Finite Field of size 2)}
        """
        no_err = self.number_errors()
        no_era = self.number_erasures()
        return "\\textnormal{Error-and-erasure channel creating %s errors and %s erasures of input space %s and output space %s}"\
                % (format_interval(no_err), format_interval(no_era), self.input_space(), self.output_space())

    def transmit_unsafe(self, message):
        r"""
        Returns ``message`` with as many errors as ``self._number_errors`` in it, and as many erasures
        as ``self._number_erasures`` in it.

        If ``self._number_errors`` was passed as an tuple for the number of errors, it will
        pick a random integer between the bounds of the tuple and use it as the number of errors.
        It does the same with ``self._number_erasures``.

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

            sage: F = GF(59)^11
            sage: n_err, n_era = 2, 2
            sage: Chan = channels.ErrorErasureChannel(F, n_err, n_era)
            sage: msg = F((3, 14, 15, 9, 26, 53, 58, 9, 7, 9, 3))
            sage: set_random_seed(10)
            sage: Chan.transmit_unsafe(msg)
            ((31, 0, 15, 9, 38, 53, 58, 9, 0, 9, 3), (0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0))
        """
        number_errors = randint(*self.number_errors())
        number_erasures = randint(*self.number_erasures())
        V = self.input_space()
        n = V.dimension()
        zero = V.base_ring().zero()

        errors = sample(xrange(n), number_errors + number_erasures)
        error_positions   = errors[:number_errors]
        erasure_positions = errors[number_errors:]

        error_vector = random_error_vector(n, V.base_ring(), error_positions)
        erasure_vector = random_error_vector(n , GF(2), erasure_positions)

        message = message + error_vector

        for i in erasure_positions:
            message[i] = zero
        return message, erasure_vector

    def number_errors(self):
        r"""
        Returns the number of errors created by ``self``.

        EXAMPLES::

            sage: n_err, n_era = 3, 0
            sage: Chan = channels.ErrorErasureChannel(GF(59)^6, n_err, n_era)
            sage: Chan.number_errors()
            (3, 3)
        """
        return self._number_errors

    def number_erasures(self):
        r"""
        Returns the number of erasures created by ``self``.

        EXAMPLES::

            sage: n_err, n_era = 0, 3
            sage: Chan = channels.ErrorErasureChannel(GF(59)^6, n_err, n_era)
            sage: Chan.number_erasures()
            (3, 3)
        """
        return self._number_erasures
