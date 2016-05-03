r"""
Benchmarking mechanism for Linear codes

This module provides a structure to easily run experiments
on linear codes, gather and manipulate results from these experiments.

Overview
========

A ``Benchmark`` takes an ``Encoder``, a ``Decoder`` and a ``Channel`` as
input::

    sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
    sage: E = C.encoder()
    sage: D = C.decoder()
    sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
    sage: B = codes.Benchmark(E, D, Chan)

Once the benchmark is created, it is possible to perform experiments using it::

    sage: B.run()

This command will trigger the following behaviour:

    - generate a random message ``m`` from the encoder's message space,
    - encode ``m``, this will return a codeword ``c``,
    - use the channel to scramble ``c``, giving a new word ``y``, and
    - decode ``y`` into ``y_dec`` using the decoder.

At every step, some data is stored to keep track of every event which
occured during a run. In the end, ``m``, ``c``, ``y``, ``y_dec`` are
stored. Furthermore, the encoding time, the decoding time and the result
of the decoding attempt (success or failure) are stored as well.
For more details on stored data, see
:meth:`sage.coding.benchmark.Benchmark.experimental_data`.
For more details on how a run works, see
:meth:`sage.coding.benchmark.Benchmark.run`

After a run, it is possible to use utility methods to get an overview
of the results. For instance, to access the median decoding time::

    sage: B.compute_timings("decoding", median) #random
    {"_0": 0.0008940000000000059}

For more details on these utility methods, see
:meth:`sage.coding.benchmark.Benchmark.compute_timings` and
:meth:`sage.coding.benchmark.Benchmark.decoding_success_rate`.

Combining benchmarks
====================

As it was presented above, this structure seems of litle interest.
Indeed, what if I want to increase the number of errors and see how
it impacts the decoding speed?

To cover this case, and all similar cases, it is possible to add
Benchmark objects.

Let's consider the following codes:

    sage: C1 = codes.GeneralizedReedSolomonCode(GF(101).list()[:100], 50)
    sage: C2 = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 125)
    sage: C3 = codes.GeneralizedReedSolomonCode(GF(503).list()[:502], 251)
    sage: C4 = codes.GeneralizedReedSolomonCode(GF(751).list()[:750], 375)
    sage: C5 = codes.GeneralizedReedSolomonCode(GF(1009).list()[:1008], 504)

We build benchmarks based on them:

    sage: chan = channels.StaticErrorRateChannel
    sage: B1 = codes.Benchmark(C1.encoder(), C1.decoder(), chan(C1.ambient_space(), C1.decoder().decoding_radius()))
    sage: B2 = codes.Benchmark(C2.encoder(), C2.decoder(), chan(C2.ambient_space(), C2.decoder().decoding_radius()))
    sage: B3 = codes.Benchmark(C3.encoder(), C3.decoder(), chan(C3.ambient_space(), C3.decoder().decoding_radius()))
    sage: B4 = codes.Benchmark(C4.encoder(), C4.decoder(), chan(C4.ambient_space(), C4.decoder().decoding_radius()))
    sage: B5 = codes.Benchmark(C5.encoder(), C5.decoder(), chan(C5.ambient_space(), C5.decoder().decoding_radius()))

We can add those benchmarks to create a big benchmark structure:

    sage: B_final = B1 + B2 + B3 + B4 + B5
    sage: B_final
    Benchmarking structure consisting of 5 benchmark primitives

Before going up to data exploitation (with plotting!), let's stop on some properties of
the addition for benchmarks:

    - it is possible to add any number of benchmarks
    - if those benchmarks identifiers are the default one (``"_0"``), Sage will
      rename them by incrementing their id by one for each benchmark but the
      first one. So, in the example above, ``B1`` will get the id ``"_0"``,
      ``B2`` will get ``"_1"``, ``B3`` will get ``"_2"`` etc.
    - if one changed their ids, it of course possible to add them as long as their ids
      are all different. If two benchmarks share the same custom id, an exception will be
      raised.
    - in the final benchmark object, any call to getter methods (:meth:`encoder`, :meth:`decoder` etc)
      will return a dictionary whose keys are the ids and values the associated object. It is possible to
      specify the id as input (e.g. ``B_final.encoder("_3") to get a specific value).
    - when the addition is performed, if some benchmark already had experimental data computed, this
      data will be propagated.

Plotting data
=============

Once one gets experimental data, it's possible to plot this data, using the provided
:meth:`plot` method.

To use it, one has to pass to it a method which generates the points to plot.
See :meth:`plot` for details.

    sage:  B_final.plot(lambda i:(B_final.decoder(i).decoding_radius(), B_final.compute_timings("decoding", median, i)), legend_label = 'number errors vs. median decoding time', plotjoined=True)

Saving benchmarks
=================

It's possible to save your experiments on your hard drive by using
the following interface provided by Sage:


    C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
    E = C.encoder()
    D = C.decoder()
    Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
    from sage.coding.benchmark import Benchmark
    B = Benchmark(E, D, Chan)
    B.start_experiments()
    f = os.path.join("path_to_file", 'bench.sobj')
    save(B, f)
    load(f)

.. WARNING::

    As this code is experimental, a warning is thrown when an
    benchmark is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
        sage: E = C.encoder()
        sage: D = C.decoder()
        sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
        sage: B = codes.Benchmark(E, D, Chan)
        doctest:...: FutureWarning: This class/method/function is marked as
        experimental. It, its functionality or its interface might change
        without a formal deprecation.
        See http://trac.sagemath.org/20526 for details.

"""

#*****************************************************************************
#       Copyright (C) 2016 David Lucas <david.lucas@inria.fr>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import time
from copy import copy
from encoder import Encoder
from decoder import Decoder, DecodingError
from channel_constructions import Channel
from sage.rings.integer import Integer
from sage.stats.basic_stats import mean, median
from sage.misc.superseded import experimental
from sage.misc.cachefunc import cached_method
from sage.misc.functional import numerical_approx
from sage.structure.sage_object import SageObject
from sage.plot.plot import list_plot

class Benchmark(SageObject):
    r"""
    Creates a benchmarking structure for linear codes.

    See documentation above to learn more on this structure.

    INPUT:

    - ``encoder`` -- (default: ``None``) an instance of :class:`sage.coding.encoder.Encoder`.
                     If set to ``None``, it creates an empty benchmark.

    - ``decoder`` -- (default: ``None``) an instance of :class:`sage.coding.decoder.Decoder`.
                     If set to ``None``, it creates an empty benchmark.

    - ``channel`` -- (default: ``None``) an instance of :class:`sage.coding.channel_constructions.Channel`.
                     If set to ``None``, it creates an empty benchmark.

    - ``identifier`` -- (default: ``None``) the name used to identify this specific benchmark.
                     If set to ``None``, this benchmark will be identified as ``"_0"``.

    .. NOTE::

        Leaving both ``encoder``, ``decoder`` and ``channel`` set to ``None`` will generate an empty
        benchmark. If only several of these fields are filled and the others left to ``None``, an
        error will be raised.

    .. WARNING::

        The default identifiers ``"_integer"`` are used when adding several benchmarks together.
        Never use these names as your own identifiers or it will break the internal behaviour
        of the addition for Benchmark objects! See above for more details on the addition of
        Benchmark objects.

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
        sage: E = C.encoder()
        sage: D = C.decoder()
        sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
        sage: codes.Benchmark(E, D, Chan)



    """

    _super_id = 0

    @experimental(trac_number=20526)
    def __init__(self, encoder=None, decoder=None, channel=None, identifier=None):
        r"""
        TESTS:

        If one passes an encoder and a decoder which are not associated with the
        same code, an exception is raised:

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: Ch = channels.StaticErrorRateChannel(C1.ambient_space(), 5)
            sage: codes.Benchmark(C1.encoder(), C2.decoder(), Ch)
            Traceback (most recent call last):
            ...
            ValueError: encoder and decoder have to be associated with the same code

        If one passes a channel which input space is not the ambient space of
        the code associated with ``encoder`` and ``decoder``, an exception
        is raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: Ch = channels.StaticErrorRateChannel(GF(59) ** 42, 5)
            sage: codes.Benchmark(C.encoder(), C.decoder(), Ch)
            Traceback (most recent call last):
            ...
            ValueError: channel's input space has to be the ambient space of the code
        """

        inp = [encoder, decoder, channel]
        self._encoders = dict()
        self._decoders = dict()
        self._channels = dict()
        self._codes = dict()
        self._ids = set()
        if identifier is None:
            self._ids.add("_0")
            local_id = "_0"
        else:
            self._ids.add(identifier)
            local_id = identifier

        if all(i is None for i in inp):
            self._encoders[local_id] = None
            self._decoders[local_id] = None
            self._codes[local_id] = None
            self._channels[local_id] = None
            self._experimental_data = dict()
            return
        if None in inp:
            raise ValueError("encoder, decoder and channel have to be all set to None or all specified")
        elif not isinstance(encoder, Encoder):
            raise ValueError("encoder has to be an Encoder object")
        elif not isinstance(decoder, Decoder):
            raise ValueError("decoder has to be a Decoder object")
        elif not isinstance(channel, Channel):
            raise ValueError("channel has to be a Channel object")
        elif not encoder.code() == decoder.code():
            raise ValueError("encoder and decoder have to be associated with the same code")
        self._encoders[local_id] = encoder
        self._decoders[local_id] = decoder
        self._codes[local_id] = encoder.code()
        if not channel.input_space() == encoder.code().ambient_space():
            raise ValueError("channel's input space has to be the ambient space of the code")
        self._channels[local_id] = channel
        self._experimental_data = dict()

    def __eq__(self, other):
        r"""
        Tests equality between Benchmark objects.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B1 = codes.Benchmark(E, D, Chan)
            sage: B2 = codes.Benchmark(E, D, Chan)
            sage: B1.__eq__(B2)
            True
        """
        return isinstance(other, Benchmark)\
                and self.encoder() == other.encoder()\
                and self.decoder() == other.decoder()\
                and self.channel() == other.channel()\
                and self.identifier() == other.identifier()

    def __ne__(self, other):
        r"""
        Tests difference between Benchmark objects.

        EXAMPLES::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: E1 = C1.encoder()
            sage: E2 = C2.encoder()
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B1 = codes.Benchmark(E1, D1, Chan1)
            sage: B2 = codes.Benchmark(E2, D2, Chan2)
            sage: B1.__ne__(B2)
            True
        """
        return not self == other

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: codes.Benchmark(E, D, Chan)

        """
        l = len(self.identifier())
        if l > 1:
            return "Benchmarking structure consisting of %s benchmark primitives" % l
        return "Benchmarking structure using %s and %s, over %s"\
                % (self.encoder(), self.decoder(), self.channel())

    def __add__(self, other):
        r"""
        Adds ``other`` to ``self`` and forms a new benchmark object.

        See documentation above to learn more on Benchmark addition.

        INPUT:

        - ``other`` -- a Benchmark object

        EXAMPLES::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: C3 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 16)
            sage: E1 = C1.encoder()
            sage: E2 = C2.encoder()
            sage: E3 = C3.encoder()
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: D3 = C3.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: Chan3 = channels.StaticErrorRateChannel(C3.ambient_space(), D3.decoding_radius())
            sage: B1 = codes.Benchmark(E1, D1, Chan1)
            sage: B2 = codes.Benchmark(E2, D2, Chan2)
            sage: B3 = codes.Benchmark(E3, D3, Chan3)
            sage: B4 = B1 + B2 + B3
            sage: B4
            Benchmarking structure consisting of 3 benchmark primitives

        TESTS:

        Adding two empty benchmarks is impossible::

            sage: B1 = codes.Benchmark()
            sage: B2 = codes.Benchmark()
            sage: B1 + B2
            Traceback (most recent call last):
            ...
            ValueError: Cannot add two enpty benchmarks

        Adding two benchmarks with the same custom id is impossible::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: E1 = C1.encoder()
            sage: E2 = C2.encoder()
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B1 = codes.Benchmark(E1, D1, Chan1, identifier = "my_id")
            sage: B2 = codes.Benchmark(E2, D2, Chan2, identifier = "my_id")
            sage: B1 + B2
            Traceback (most recent call last):
            ...
            ValueError("Benchmark with custom id 'my_id' is already in self" % i
        """
        if not isinstance(other, Benchmark):
            raise TypeError("%s must be a Benchmark object" % other)
        #if both self and other are empty, do nothing
        if (len(self.identifier()) == 1 and self.encoder() is None
                and len(other.identifier()) == 1 and other.encoder() is None):
            raise ValueError("Cannot add two empty benchmarks")
        if self == other:
            return self
        import re
        pattern_default_id = re.compile('_[0-9]+')
        b = Benchmark()
        #if self is None, we put other in b and immediately return this object
        if len(self.identifier()) == 1 and self.encoder() is None:
            b._encoders = copy(other._encoders)
            b._decoders = copy(other._decoders)
            b._channels = copy(other._channels)
            b._codes = copy(other._codes)
            b._ids = copy(other.identifier())
            b._experimental_data = copy(other.experimental_data())
            return b
        else:
            b._encoders = copy(self._encoders)
            b._decoders = copy(self._decoders)
            b._channels = copy(self._channels)
            b._codes = copy(self._codes)
            b._ids = copy(self.identifier())
            b._experimental_data = copy(self.experimental_data())
        for i in other.identifier():
            #if i is a custom id, which is already known by self, raise an Exception
            if i in b.identifier() and isinstance(i, str) and pattern_default_id.match(i) is None:
                raise ValueError("Benchmark with custom id %s is already in self" % i)
            #if other's encoder, decoder, channel and code are already know by self, do nothing
            elif (i in b.identifier() and self.encoder(i) == other.encoder(i)
                    and self.decoder(i) == other.decoder(i)
                    and self.channel(i) == other.channel(i)
                    and self.code(i) == self.code(i)):
                pass
            #if i is a default id, change it and add other to self
            elif i in b.identifier():
                if other.encoder(identifier = i) is None:
                    pass
                else:
                    Benchmark._super_id += 1
                    new_id = "_"+str(Benchmark._super_id)
                    b._ids.add(new_id)
                    b._encoders[new_id] = other._encoders[i]
                    b._decoders[new_id] = other._decoders[i]
                    b._channels[new_id] = other._channels[i]
                    b._codes[new_id] = other._codes[i]
                    other_exp_data = other.experimental_data()
                    if len(other_exp_data) is not 0:
                        run = 0
                        while True:
                            try:
                                b._experimental_data[new_id, run] = other_exp_data[i, run]
                                run += 1
                            except KeyError:
                                break
            #if i in a custom id, add other as is to self
            else:
                if other.encoder(identifier = i) is None:
                    pass
                else:
                    b._ids.add(i)
                    b._encoders[i] = other._encoders[i]
                    b._decoders[i] = other._decoders[i]
                    b._channels[i] = other._channels[i]
                    b._codes[i] = other._codes[i]
                    other_exp_data = other.experimental_data()
                    if len(other_exp_data) is not 0:
                        run = 0
                        while True:
                            try:
                                b._experimental_data[i, run] = other_exp_data[i, run]
                                run += 1
                            except KeyError:
                                break
        return b

    def encoder(self, identifier=None):
        r"""
        Returns the encoder used by ``self`` to perform experiments.

        If ``self`` is a combination of several benchmarks, it possesses a list
        of several encoders. In that case, one can specify the ``identifier`` of
        a sub-benchmark to access its encoder.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark which
          encoder will be returned. If it is set to ``None``, the complete dictionary
          of encoders will be returned.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B.encoder()
            Evaluation vector-style encoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        enc = self._encoders
        if identifier is None:
            return (enc.values()[0] if len(enc) == 1 else enc)
        else:
            if identifier not in self.identifier():
                raise ValueError("Given identifier is not known by self")
            return enc[identifier]

    def decoder(self, identifier=None):
        r"""
        Returns the decoder used by ``self`` to perform experiments.

        If ``self`` is a combination of several benchmarks, it possesses a list
        of several decoders. In that case, one can specify the ``identifier`` of
        a sub-benchmark to access its decoder.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark which
          decoder will be returned. If it is set to ``None``, the complete dictionary
          of decoders will be returned.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B.decoder()
            Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        dec = self._decoders
        if identifier is None:
            return (dec.values()[0] if len(dec) == 1 else dec)
        else:
            if identifier not in self.identifier():
                raise ValueError("Given identifier is not known by self")
            return dec[identifier]

    def channel(self, identifier=None):
        r"""
        Returns the channel used by ``self`` to perform experiments.

        If ``self`` is a combination of several benchmarks, it possesses a list
        of several channels. In that case, one can specify the ``identifier`` of
        a sub-benchmark to access its channel.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark which
          channel will be returned. If it is set to ``None``, the complete dictionary
          of channels will be returned.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B.channel()
            Static error rate channel creating 14 errors, of input and output space Vector space of dimension 40 over Finite Field of size 59
        """
        chan = self._channels
        if identifier is None:
            return (chan.values()[0] if len(chan) == 1 else chan)
        else:
            if identifier not in self.identifier():
                raise ValueError("Given identifier is not known by self")
            return chan[identifier]

    def code(self, identifier=None):
        r"""
        Returns the code used by ``self`` to perform experiments.

        If ``self`` is a combination of several benchmarks, it possesses a list
        of several codes. In that case, one can specify the ``identifier`` of
        a sub-benchmark to access its code.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark which
          code will be returned. If it is set to ``None``, the complete dictionary
          of codes will be returned.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B.code()
            [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59
        """
        c = self._codes
        if identifier is None:
            return (c.values()[0] if len(c) == 1 else c)
        else:
            if identifier not in self.identifier():
                raise ValueError("Given identifier is not known by self")
            return c[identifier]

    def experimental_data(self):
        r"""
        Returns the dictionary which contains all the experimental data generated.

        EXAMPLES::
            sage: C = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B.start_experiments(no_tests = 1)
            sage: B.experimental_data() #random
        """
        return self._experimental_data

    def identifier(self):
        r"""
        Returns the set of identifiers of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.identifier()
            {"_0"}
        """
        return self._ids

    def clear_experimental_data(self):
        r"""
        Erases all experimental data gathered during previous runs.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.start_experiments()
            sage: len(B.experimental_data())
            100
            sage: B.clear_experimental_data()
            sage: len(B.experimental_data())
            0
        """
        self.experimental_data().clear()

    def _perform_experiments_for_single_id(self, identifier, no_tests, verbosity):
        r"""
        Runs a ``no_tests`` rounds of experimental data computation using ``self``'s
        parameters.

        See :meth:`run` for details.

        INPUT:

        - ``no_tests`` -- the number of rounds to perform during this run.

        - ``verbosity_level`` -- a number which indicates how verbose this run will
          be:
            - ``0`` means completely silent, nothing will be written on the output
            - ``1`` means small verbosity, one message will be written every time 25% of the tests are done
            - ``2`` means very verbose, one message will be written every time 10% of the tests are done


        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: from sage.coding.benchmark import Benchmark
            sage: B = Benchmark(E, D, Chan)
            sage: B._perform_experiments_for_single_id('_0', 100, 0)
        """
        #setting local variables and checking validity of data to use
        data = self.experimental_data()
        try:
            data[identifier, 0]
            print "Some experimental data already exist for benchmark %s.\nIf you want to perform new experiments, please run clear_experimental_data first" % identifier
            return
        except KeyError:
            pass
        C = self.code(identifier)
        if C is None:
            raise ValueError("Impossible to run tests on an empty benchmark")
        E = self.encoder(identifier)
        D = self.decoder(identifier)
        Chan = self.channel(identifier)
        previous = 0 #used for verbosity
        if verbosity != 0:
            print "Starting run for benchmark %s" % identifier

        for i in range(no_tests):
            m = E.message_space().random_element()
            start = time.clock()
            c = E(m)
            encode_time = time.clock() - start
            y = Chan(c)
            try:
                start = time.clock()
                dec_y = D.decode_to_code(y)
                decode_time = time.clock() - start
                decoding_success = True
            except DecodingError:
                decode_time = time.clock() - start
                dec_y = None #as decoding failed
                decoding_success = False
            data[identifier, i] = (m, c, encode_time, y, dec_y, decode_time, decoding_success)

            #verbosity
            if verbosity ==1:
                cur = (i*4) // no_tests
                if (i*4) // no_tests != previous:
                    previous = cur
                    print "%s percent complete" % (cur * 25)
            if verbosity ==2:
                cur = (i*10) // no_tests
                if (i*10) // no_tests != previous:
                    previous = cur
                    print "%s percent complete" % (cur * 10)

        if verbosity != 0:
            print "Run complete for benchmark %s" % identifier

    def run(self, no_tests = 100, verbosity_level = 0):
        r"""
        Runs a ``no_tests`` rounds of experimental data computation using ``self``'s
        parameters. The user can be informed about the current run's status by setting
        ``verbosity_level`` to the appropriate value.

        A round of experimental data computation consists in:

        - generating a random message ``m`` from the ``self``'s encoder's message space,
        - encoding ``m``, this will return a codeword ``c``,
        - using the ``self``'s channel to scramble ``c``, giving a new word ``y``, and
        - decoding ``y`` into ``y_dec`` using ``self``'s decoder.

        INPUT:

        - ``no_tests`` -- (default: ``100``) the number of rounds to perform during this run.

        - ``verbosity_level`` -- (default: ``0``) a number which indicates how verbose this run will
          be:
            - ``0`` means completely silent, nothing will be written on the output,
            - ``1`` means small verbosity, one message will be written every time 25% of the tests are done,
            - ``2`` means very verbose, one message will be written every time 10% of the tests are done.

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark
          whose success rate will be computed. If it is let to ``None``, success
          rate of every subbenchmark will be computed.

        .. NOTE::

            This methods has a safeguard to prevent erasing experimental data by accident:
            if :meth:`experimental_data` is not empty, nothing will be done.
            See :meth:`clear_experimental_data` to delete data done in previous runs.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.run()
        Another one in verbose::

            sage: B.clear_experimental_data()
            sage: B.run(100, 1)
            Starting run for benchmark _0
            25.0 percent complete
            50.0 percent complete
            75.0 percent complete
            Run complete for benchmark _0
        """
        if not isinstance(no_tests, (int, Integer)):
            raise ValueError("no_tests has to be a positive Python int or Sage Integer")
        if not no_tests > 0:
            raise ValueError("no_tests has to be a positive Python int or Sage Integer")
        if not verbosity_level in {0,1,2}:
            raise ValueError("verbosity_level has to be either 0, 1 or 2")
        for i in self.identifier():
            self._perform_experiments_for_single_id(i, no_tests, verbosity_level)

    def _process_experimental_data(self, identifier, target, processing):
        r"""
        Returns processed ``target``'s data using ``processing`` for ``self``
        designated by ``identifier``.

        INPUT:

        - ``identifier`` -- the identifier of the benchmark
          whose success rate will be computed.

        - ``target`` -- the operation whose timings will be computed, can be
          either ``"encoding"`` or ``"decoding"``

        - ``processing`` -- the processing to apply on ``target``'s timings.
          It has to be a Sage/Python method. It can be either a native one
          (e.g. ``mean`` or ``median``) or one implemented by the user.
          If it's a custom one, it has to take a Python list as input.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.run()
            sage: B._process_experimental_data('_0', "decoding", median) #random
            0.0008940000000000059
        """
        if target is "encoding":
            pos = 2
        elif target is "decoding":
            pos = 5
        else:
            raise ValueError("target has to be set to either \"encoding\" or \"decoding\"")
        i = 0
        timings = []
        data = self.experimental_data()
        while True:
            try:
                timings.append(data[identifier, i][pos])
                i += 1
            except KeyError:
                break
        return processing(timings)

    @cached_method
    def compute_timings(self, target, processing, identifier = None):
        r"""
        Returns ``target``'s timings after applying ``processing`` on them.

        INPUT:

        - ``target`` -- the operation whose timings will be computed, can be
          either ``"encoding"`` or ``"decoding"``

        - ``processing`` -- the processing to apply on ``target``'s timings.
          It has to be a Sage/Python method. It can be either a native one
          (e.g. ``mean`` or ``median``) or one implemented by the user.
          If it's a custom one, it has to take a Python list as input.

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark
          whose success rate will be computed. If it is let to ``None``, success
          rate of every subbenchmark will be computed.

        OUTPUT:

        - a dictionary whose keys are the benchmarks' identifiers and values the associated
          timings if ``identifier`` is ``None``, or

        - the timing of the benchmark designated by ``identifier`` if it isn't ``None``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.run()
            sage: B.compute_timings("decoding", median) #random
            {'_0', 0.0008940000000000059}

        It is also possible to use your own processing method::

            sage: def my_dumb_processing(l):
            ....:    return -1
            sage: B.compute_timings("decoding", my_dumb_processing)
            {'_0', -1}
        """
        if not target in ("decoding", "encoding"):
            raise ValueError("target has to be set to either \"encoding\" or \"decoding\"")
        if not hasattr(processing, '__call__'):
            raise ValueError("processing has to be a method")
        data = self.experimental_data()
        ids = self.identifier()
        if len(data) == 0:
            raise ValueError("No experimental data has been computed!")
        if identifier is not None and identifier not in ids:
            raise ValueError("Unknown identifier")
        elif identifier is None:
            timings = dict()
            for i in ids:
                timings[i] = self._process_experimental_data(i, target, processing)
            return timings
        elif identifier is not None:
            return self._process_experimental_data(identifier,target, processing)

    @cached_method
    def decoding_success_rate(self, identifier = None):
        r"""
        Returns the success rate of the decoding attempts.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark
          whose success rate will be computed. If it is let to ``None``, success
          rate of every subbenchmark will be computed.

        OUTPUT:

        - a dictionary whose keys are the benchmarks' identifiers and values the associated
          success rates if ``identifier`` is ``None``, or

        - the success rate of the benchmark designated by ``identifier`` if it isn't ``None``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: E = C.encoder()
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(E, D, Chan)
            sage: B.run()
            sage: B.decoding_success_rate()
            {'_0': 1.00}
        """
        def compute_success_rate(identifier):
            i = 0
            success_sum = 0
            while True:
                try:
                    success_sum += data[identifier, i][6]
                    i += 1
                except KeyError:
                    break
            return numerical_approx(success_sum/i, digits = 3)

        data = self.experimental_data()
        ids = self.identifier()
        if len(data) == 0:
            raise ValueError("No experimental data has been computed!")
        if identifier is not None and identifier not in ids:
            raise ValueError("Unknown identifier")
        elif identifier is None:
            rates = dict()
            for i in ids:
                rates[i] = compute_success_rate(i)
            return rates
        elif identifier is not None:
            return compute_success_rate(i)

    def plot(self, points_generator, **kwargs):
        r"""
        Returns a plot of points generator by ``points_generator``.

        INPUT:

        - ``points_generator`` -- a method which takes an identifier in input
          and returns a tuple

        - ``**kwargs`` -- all extra arguments are forwarded to :meth:`sage.plot.plot.list_plot`

        EXAMPLES::


            sage: C1 = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 125)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 150)
            sage: C3 = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 175)
            sage: E1 = C1.encoder()
            sage: E2 = C2.encoder()
            sage: E3 = C3.encoder()
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: D3 = C3.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: Chan3 = channels.StaticErrorRateChannel(C3.ambient_space(), D3.decoding_radius())
            sage: B1 = codes.Benchmark(E1, D1, Chan1)
            sage: B2 = codes.Benchmark(E2, D2, Chan2)
            sage: B3 = codes.Benchmark(E3, D3, Chan3)
            sage: B4 = B1 + B2 + B3
            sage: B4.run()
            sage: B4.plot(lambda i: (B4.decoder(i).decoding_radius(), B4.compute_timings("decoding", median, i)))
        """
        if not hasattr(points_generator, '__call__'):
            raise ValueError("points_generator has to be a method")
        l = [points_generator(i) for i in self.identifier()]
        return list_plot(l, **kwargs)
