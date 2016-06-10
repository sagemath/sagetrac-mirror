r"""
Benchmarking mechanism for Linear codes

This module provides a structure to easily run experiments
on linear codes, gather and manipulate results from these experiments.

.. WARNING::

    As this code is experimental, a warning is thrown when a
    benchmark is created for the first time in a session (see
    :class:`sage.misc.superseded.experimental`).

    TESTS::

        sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
        sage: D = C.decoder()
        sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
        sage: B = codes.Benchmark(C, D, Chan)
        doctest:...: FutureWarning: This class/method/function is marked as experimental. It, its functionality or its interface might change without a formal deprecation.
        See http://trac.sagemath.org/20526 for details.

Overview
========

A ``Benchmark`` takes a ``Code``, a ``Decoder`` and a ``Channel`` as
input::

    sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
    sage: D = C.decoder()
    sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
    sage: B = codes.Benchmark(C, D, Chan)

Once the benchmark is created, it is possible to perform experiments using it::

    sage: B.run()

This command will trigger the following behaviour:

    - generate a random codeword ``c`` in ``self``'s code,
    - use the channel to scramble ``c``, giving a new word ``y``, and
    - decode ``y`` into ``y_dec`` using the decoder.

At every step, some data is stored to keep track of every event which
occured during a run. In the end, ``m``, ``c``, ``y``, ``y_dec`` are
stored. Furthermore, the time used to generate ``c``, the decoding time
and the result
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

Let's consider the following codes::

    sage: C1 = codes.GeneralizedReedSolomonCode(GF(101).list()[:100], 50)
    sage: C2 = codes.GeneralizedReedSolomonCode(GF(251).list()[:250], 125)
    sage: C3 = codes.GeneralizedReedSolomonCode(GF(503).list()[:502], 251)
    sage: C4 = codes.GeneralizedReedSolomonCode(GF(751).list()[:750], 375)
    sage: C5 = codes.GeneralizedReedSolomonCode(GF(1009).list()[:1008], 504)

We build benchmarks based on them::

    sage: chan = channels.StaticErrorRateChannel
    sage: B1 = codes.Benchmark(C1, C1.decoder(), chan(C1.ambient_space(), C1.decoder().decoding_radius()))
    sage: B2 = codes.Benchmark(C2, C2.decoder(), chan(C2.ambient_space(), C2.decoder().decoding_radius()))
    sage: B3 = codes.Benchmark(C3, C3.decoder(), chan(C3.ambient_space(), C3.decoder().decoding_radius()))
    sage: B4 = codes.Benchmark(C4, C4.decoder(), chan(C4.ambient_space(), C4.decoder().decoding_radius()))
    sage: B5 = codes.Benchmark(C5, C5.decoder(), chan(C5.ambient_space(), C5.decoder().decoding_radius()))

We can add those benchmarks to create a big benchmark structure::

    sage: B_final = B1 + B2 + B3 + B4 + B5
    sage: B_final
    Benchmarking structure consisting of 5 benchmark primitives

Before going up to data exploitation (with plotting!), let's stop on some properties of
the addition for benchmarks:

    - It is possible to add any number of benchmarks.

    - If those benchmarks identifiers are the default one (``"_0"``), Sage will
      rename them by incrementing their id by one for each benchmark but the
      first one. So, in the example above, ``B1`` will get the id ``"_0"``,
      ``B2`` will get ``"_1"``, ``B3`` will get ``"_2"`` etc.

    - If one changed their ids, it is of course possible to add them as long as their ids
      are all different. If two benchmarks share the same custom id, an exception will be
      raised.

    - In the final benchmark object, any call to getter methods
      (:meth:`sage.coding.benchmark.Benchmark.code`, :meth:`sage.coding.benchmark.Benchmark.decoder` etc)
      will return a dictionary whose keys are the ids and values the associated object. It is possible to
      specify the id as input (e.g. ``B_final.decoder("_3")`` to get a specific value).

    - When the addition is performed, if some benchmark already had experimental data computed, this
      data will be propagated.

Plotting data
=============

Once one gets experimental data, it's possible to plot this data, using the provided
:meth:`plot` method.

To use it, one has to pass to it a method which generates the points to plot.
See :meth:`plot` for details.
We give here an exemple of what one can do::

    B_final.plot(lambda i:(B_final.decoder(i).decoding_radius(), B_final.compute_timings("decoding", median, i)), legend_label = 'number errors vs. median decoding time', plotjoined=True)


Saving benchmarks
=================

It's possible to save your experiments on your hard drive by using
the following interface provided by Sage::


    C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
    D = C.decoder()
    Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
    B = codes.Benchmark(C, D, Chan)
    B.start_experiments()
    f = os.path.join("path_to_file", 'bench.sobj')
    save(B, f)
    load(f)
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
from linear_code import AbstractLinearCode
from decoder import Decoder, DecodingError
from channel_constructions import Channel
from sage.rings.integer import Integer
from sage.stats.basic_stats import mean, median
from sage.misc.superseded import experimental
from sage.misc.cachefunc import cached_method
from sage.misc.functional import numerical_approx
from sage.structure.sage_object import SageObject
from sage.plot.plot import list_plot
from sage.parallel.decorate import parallel
from sage.functions.other import ceil

class Benchmark(SageObject):
    r"""
    Creates a benchmarking structure for linear codes.

    See documentation above to learn more on this structure.

    INPUT:

    - ``code`` -- (default: ``None``) a linear code.
                     If set to ``None``, it creates an empty benchmark.

    - ``decoder`` -- (default: ``None``) an instance of :class:`sage.coding.decoder.Decoder`.
                     If set to ``None``, it creates an empty benchmark.

    - ``channel`` -- (default: ``None``) an instance of :class:`sage.coding.channel_constructions.Channel`.
                     If set to ``None``, it creates an empty benchmark.

    - ``identifier`` -- (default: ``None``) the name used to identify this specific benchmark.
                     If set to ``None``, this benchmark will be identified as ``"_0"``.

    - ``no_tests`` -- (default: ``100``) the number of tests to perform on ``self`` when doing
                      experiments.

    .. NOTE::

        Leaving both ``code``, ``decoder`` and ``channel`` set to ``None`` will generate an empty
        benchmark. If only several of these fields are filled and the others left to ``None``, an
        error will be raised.

    .. WARNING::

        The default identifiers ``"_integer"`` are used when adding several benchmarks together.
        Never use these names as your own identifiers or it will break the internal behaviour
        of the addition for Benchmark objects! See above for more details on the addition of
        Benchmark objects.

    EXAMPLES::

        sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
        sage: D = C.decoder()
        sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
        sage: codes.Benchmark(C, D, Chan)
        Benchmarking structure using [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59 and Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59, over Static error rate channel creating 14 errors, of input and output space Vector space of dimension 40 over Finite Field of size 59
    """

    _super_id = 0

    @experimental(trac_number=20526)
    def __init__(self, code=None, decoder=None, channel=None, identifier=None, no_tests=100):
        r"""
        TESTS:

        If ``decoder`` is not associated with ``code``, an exception is raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 18)
            sage: Ch = channels.StaticErrorRateChannel(GF(59) ** 42, 5)
            sage: codes.Benchmark(C, C1.decoder(), Ch)
            Traceback (most recent call last):
            ...
            ValueError: decoder has to be associated with code


        If one passes a channel whose input space is not the ambient space of
        ``code``, an exception is raised::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: Ch = channels.StaticErrorRateChannel(GF(59) ** 42, 5)
            sage: codes.Benchmark(C, C.decoder(), Ch)
            Traceback (most recent call last):
            ...
            ValueError: channel's input space has to be the ambient space of the code
        """

        inp = [code, decoder, channel]
        self._codes = dict()
        self._decoders = dict()
        self._channels = dict()
        self._no_tests = dict()
        self._ids = set()
        if identifier is None:
            self._ids.add("_0")
            local_id = "_0"
        else:
            self._ids.add(identifier)
            local_id = identifier

        if all(i is None for i in inp):
            self._codes[local_id] = None
            self._decoders[local_id] = None
            self._channels[local_id] = None
            self._experimental_data = dict()
            return
        if None in inp:
            raise ValueError("code, decoder and channel have to be all set to None or all specified")
        elif not isinstance(code, AbstractLinearCode):
            raise ValueError("code has to be an AbstractLinearCode object")
        elif not isinstance(decoder, Decoder):
            raise ValueError("decoder has to be a Decoder object")
        elif not isinstance(channel, Channel):
            raise ValueError("channel has to be a Channel object")
        elif not code == decoder.code():
            raise ValueError("decoder has to be associated with code")
        if not isinstance(no_tests, (int, Integer)):
            raise ValueError("no_tests has to be a positive Python int or Sage Integer")
        if not no_tests > 0:
            raise ValueError("no_tests has to be a positive Python int or Sage Integer")
        #Puts generator matrix and parity check matrix in the local cache
        G = C.generator_matrix()
        H = C.parity_check_matrix()
        self._codes[local_id] = code
        self._decoders[local_id] = decoder
        self._no_tests[local_id] = no_tests
        if not channel.input_space() == code.ambient_space():
            raise ValueError("channel's input space has to be the ambient space of the code")
        self._channels[local_id] = channel
        self._experimental_data = dict()

    def __eq__(self, other):
        r"""
        Tests equality between Benchmark objects.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B1 = codes.Benchmark(C, D, Chan)
            sage: B2 = codes.Benchmark(C, D, Chan)
            sage: B1.__eq__(B2)
            True
        """
        return isinstance(other, Benchmark)\
                and self.code() == other.code()\
                and self.decoder() == other.decoder()\
                and self.channel() == other.channel()\
                and self.identifier() == other.identifier()

    def __ne__(self, other):
        r"""
        Tests difference between Benchmark objects.

        EXAMPLES::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1)
            sage: B2 = codes.Benchmark(C2, D2, Chan2)
            sage: B1.__ne__(B2)
            True
        """
        return not self == other

    def _repr_(self):
        r"""
        Returns a string representation of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: codes.Benchmark(C, D, Chan)
            Benchmarking structure using [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59 and Gao decoder for [40, 12, 29] Generalized Reed-Solomon Code over Finite Field of size 59, over Static error rate channel creating 14 errors, of input and output space Vector space of dimension 40 over Finite Field of size 59

        If ``self`` is an empty benchmark, a specific string representation
        is given::

            sage: codes.Benchmark()
            Empty Benchmark object

        If ``self`` is the result of an addition between several benchmarks,
        a specific string representation is given::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1)
            sage: B2 = codes.Benchmark(C2, D2, Chan2)
            sage: B1 + B2
            Benchmarking structure consisting of 2 benchmark primitives
        """
        if self.code() is None:
            return "Empty Benchmark object"
        l = len(self.identifier())
        if l > 1:
            return "Benchmarking structure consisting of %s benchmark primitives" % l
        return "Benchmarking structure using %s and %s, over %s"\
                % (self.code(), self.decoder(), self.channel())

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
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: D3 = C3.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: Chan3 = channels.StaticErrorRateChannel(C3.ambient_space(), D3.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1)
            sage: B2 = codes.Benchmark(C2, D2, Chan2)
            sage: B3 = codes.Benchmark(C3, D3, Chan3)
            sage: B4 = B1 + B2 + B3
            sage: B4
            Benchmarking structure consisting of 3 benchmark primitives

        TESTS:

        Adding two benchmarks with the same custom id is impossible::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 14)
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1, identifier = "my_id")
            sage: B2 = codes.Benchmark(C2, D2, Chan2, identifier = "my_id")
            sage: B1 + B2
            Traceback (most recent call last):
            ...
            ValueError: Benchmark with custom id my_id is already in self
        """
        if not isinstance(other, Benchmark):
            raise TypeError("%s must be a Benchmark object" % other)
        #if both self and other are empty, return an empty benchmark
        #if they are the same, return self
        if (len(self.identifier()) == 1 and self.code() is None
                and len(other.identifier()) == 1 and other.code() is None
                or self == other):
            return self
        import re
        pattern_default_id = re.compile('_[0-9]+')
        b = Benchmark()
        #if self is None, we put other in b and immediately return this object
        if len(self.identifier()) == 1 and self.code() is None:
            b._codes = copy(other._codes)
            b._decoders = copy(other._decoders)
            b._channels = copy(other._channels)
            b._ids = copy(other.identifier())
            b._experimental_data = copy(other.experimental_data())
            b._no_tests = copy(other._no_tests)
            return b
        else:
            b._codes = copy(self._codes)
            b._decoders = copy(self._decoders)
            b._channels = copy(self._channels)
            b._ids = copy(self.identifier())
            b._experimental_data = copy(self.experimental_data())
            b._no_tests = copy(self._no_tests)
        for i in other.identifier():
            #if i is a custom id, which is already known by self, raise an Exception
            if i in b.identifier() and isinstance(i, str) and pattern_default_id.match(i) is None:
                raise ValueError("Benchmark with custom id %s is already in self" % i)
            #if other's codes, decoder, channel and code are already know by self, do nothing
            elif (i in b.identifier() and self.code(i) == other.code(i)
                    and self.decoder(i) == other.decoder(i)
                    and self.channel(i) == other.channel(i)):
                pass
            #if i is a default id, change it and add other to self
            elif i in b.identifier():
                if other.code(identifier = i) is None:
                    pass
                else:
                    Benchmark._super_id += 1
                    new_id = "_"+str(Benchmark._super_id)
                    b._ids.add(new_id)
                    b._codes[new_id] = other._codes[i]
                    b._decoders[new_id] = other._decoders[i]
                    b._channels[new_id] = other._channels[i]
                    b._no_tests[new_id] = other.number_of_tests(i)
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
                if other.code(identifier = i) is None:
                    pass
                else:
                    b._ids.add(i)
                    b._codes[i] = other._codes[i]
                    b._decoders[i] = other._decoders[i]
                    b._channels[i] = other._channels[i]
                    b._no_tests[new_id] = other.number_of_tests(i)
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
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
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
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
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
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
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

    def number_of_tests(self, identifier=None):
        r"""
        Returns the number of tests done by ``self`` while performing experiments.

        If ``self`` is a combination of several benchmarks, it possesses a list
        of several number of tests. In that case, one can specify the ``identifier`` of
        a sub-benchmark to access its specific number of tests argument.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark whose
          number of tests will be returned. If it is set to ``None``, the complete dictionary
          of number of tests will be returned.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.number_of_tests()
            100
        """
        nt = self._no_tests
        if identifier is None:
            return (nt.values()[0] if len(nt) == 1 else nt)
        else:
            if identifier not in self.identifier():
                raise ValueError("Given identifier is not known by self")
            return nt[identifier]

    def experimental_data(self, identifier=None, data=None):
        r"""
        Returns the dictionary which contains all the experimental data generated.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark whose
          experimental_data will be returned. If it is set to ``None``, the complete dictionary
          of experimental_data will be returned.

        - ``data`` -- (default: ``None``) the specific data to return.
          It can be either:

          - ``"codeword"``, which will return the whole list of codewords,
          - ``"generate_codeword_time"``, which will return the whole list of timings measured while
            generating the codeword,
          - ``"scrambled_codeword"``, which will return the whole list of codewords as they were
            after being scrambled by ``self``'s channel,
          - ``"decoded_word"``, which will return the whole list of scrambled codewords as they were
            after being decoded by ``self``'s decoder,
          - ``"decoding_time"``, which will return the whole list of timings measured while performing
            decoding attempts,
          - ``"decoding_error"``, which will return the whole list of decoding errors registered while
            performing decoding attempts,
          - ``"decoding_failure"``, which will return the whole list of decoding failures registered while
            performing decoding attempts, or
          - ``None``, which will return the whole dictionary of experimental data.

        .. NOTE::

            If ``identifier`` is set to ``None``, and ``data`` is set to a specific keyword

        EXAMPLES::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: D1 = C1.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1, no_tests = 2)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: D2 = C2.decoder()
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: B2 = codes.Benchmark(C2, D2, Chan2, no_tests = 2)
            sage: B = B1 + B2
            sage: B.run()
            sage: B.experimental_data() #random
            {('_0', 0): ((5, 3, 6, 0, 6, 3),
              0.0006960000000000299,
              (5, 5, 6, 0, 6, 3),
              (5, 3, 6, 0, 6, 3),
              0.0025980000000001002,
              False,
              False),
             ('_0', 1): ((5, 6, 5, 2, 4, 4),
              0.0002469999999998862,
              (5, 1, 5, 2, 4, 4),
              (5, 6, 5, 2, 4, 4),
              0.0006909999999999972,
              False,
              False),
             ('_1', 0): ((2, 4, 4, 2, 5, 6),
              0.021160000000000068,
              (2, 4, 4, 1, 5, 6),
              (2, 4, 4, 2, 5, 6),
              0.026135999999999937,
              False,
              False),
             ('_1', 1): ((5, 4, 4, 5, 0, 3),
              0.0007790000000000852,
              (5, 4, 4, 1, 0, 3),
              (5, 4, 4, 5, 0, 3),
              0.0008300000000001084,
              False,
              False)}

        To access directly some specific data, one can set the data to fetch::

            sage: B.experimental_data(data="decoding_time") #random
            {('_0', 0): 0.0025980000000001002,
             ('_0', 1): 0.0006909999999999972,
             ('_1', 0): 0.026135999999999937,
             ('_1', 1): 0.0008300000000001084}

        And even set the id of the sub-benchmark::

            sage: B.experimental_data(identifier="_0", data="decoding_time") #random
            {('_0', 0): 0.0025980000000001002,
             ('_0', 1): 0.0006909999999999972}
        """
        exp_data = self._experimental_data
        ids = self.identifier()
        allowed_input = ["codeword", "generate_codeword_time", "scrambled_codeword",
                "decoded_word", "decoding_time", "decoding_error", "decoding_failure", None]
        if data not in allowed_input:
            raise ValueError("data must be one of the allowed values. See documentation for details")

        if data is not None:
            data_id = allowed_input.index(data)

        if identifier is None and data is None:
            return exp_data
        elif identifier is None:
            selected_data = dict()
            for i in ids:
                if exp_data.has_key((i, 0)):
                    no_tests = self.number_of_tests(i)
                    for j in range(no_tests):
                        try:
                            selected_data[i, j] = exp_data[i, j][data_id]
                        except KeyError:
                            break
            return selected_data

        if identifier not in self.identifier():
            raise ValueError("Given identifier is not known by self")

        selected_data = dict()
        no_tests = self.number_of_tests(identifier)
        if exp_data.has_key((identifier, 0)):
            if data is None:
                for n in range(no_tests):
                    try:
                        selected_data[identifier, n] = exp_data[identifier, n]
                    except KeyError:
                        break
            else:
                for n in range(no_tests):
                    try:
                        selected_data[identifier, n] = exp_data[identifier, n][data_id]
                    except KeyError:
                        break
        return selected_data

    def identifier(self):
        r"""
        Returns the set of identifiers of ``self``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.identifier()
            {'_0'}
        """
        return self._ids

    def clear_experimental_data(self):
        r"""
        Erases all experimental data gathered during previous runs.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(7).list()[:6], 3)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.run()
            sage: len(B.experimental_data())
            100
            sage: B.clear_experimental_data()
            sage: len(B.experimental_data())
            0
        """
        self.experimental_data().clear()

    @parallel
    def _perform_parallel_experiments_for_single_id(self, identifier, no_tests):
        r"""
        Performs rounds of experimental data computation using ``self``'s
        parameters.

        This method is used when :meth:`run` is asked to be run in its multithreaded version.

        See :meth:`run` for details.

        INPUT:

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark on which experiments
          will be performed.

        - ``no_tests`` -- the number of experiments to run

        OUTPUT:

        - the experimental data computed during this run

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B._perform_parallel_experiments_for_single_id('_0', 1)#random
            [((30, 35, 30, 25, 36, 52, 31, 4, 15, 23, 8, 0, 26, 30, 32, 12, 11, 45, 47, 23, 35, 52, 54, 33, 37, 17, 46, 45, 29, 30, 49, 58, 36, 1, 34, 34, 20, 16, 18, 45),
            0.0020709999999999895,
            (30, 54, 30, 6, 36, 52, 23, 4, 49, 23, 8, 0, 57, 30, 32, 12, 32, 45, 47, 23, 48, 52, 53, 3, 37, 17, 46, 45, 55, 30, 49, 58, 40, 47, 45, 10, 20, 16, 18, 45),
            (30, 35, 30, 25, 36, 52, 31, 4, 15, 23, 8, 0, 26, 30, 32, 12, 11, 45, 47, 23, 35, 52, 54, 33, 37, 17, 46, 45, 29, 30, 49, 58, 36, 1, 34, 34, 20, 16, 18, 45),
            0.004525000000000001,
            False,
            False)]
        """
        cur_no_tests = 0
        C = self.code(identifier)
        D = self.decoder(identifier)
        Chan = self.channel(identifier)
        results = []
        for i in range(no_tests):
            start = time.clock()
            c = C.random_element()
            generate_c_time = time.clock() - start
            y = Chan(c)
            decoding_error = True
            try:
                start = time.clock()
                dec_y = D.decode_to_code(y)
                decode_time = time.clock() - start
                decoding_failure = False
                if "list-decoder" in D.decoder_type():
                    decoding_error = (c not in dec_y)
                else:
                    decoding_error = (c != dec_y)
            except DecodingError:
                decode_time = time.clock() - start
                dec_y = None #as decoding failed
                decoding_failure = True
            results.append((c, generate_c_time, y, dec_y, decode_time, decoding_error, decoding_failure))

        return results

    def _task_master(self, verbosity):
        r"""
        Setups and runs experiments for the multithreaded version of :meth:`run`.

        The multithreaded run as follows:

            - It begins with a setup phase. During the setup phase, 4 tests are
              runned per sub-benchmark of ``self``.

            - Then, it takes the median of the decoding timings measured on the
              previous step. It uses this timings to generate tasks which takes
              approximately 3 seconds to be performed.

            - And it runs all the tasks generated above.

        INPUT:

        - ``verbosity`` -- a boolean. If set to ``True``, it will send an information
          to the user everytime a task has been completed by a thread.


        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B._task_master(False)
        """
        def preparse_tasks(tasks_preparsing):
            r"""
            Prepares the list of tasks to execute
            according to ``tasks_preparsing``.
            """
            tasks = []
            while len(tasks_preparsing) != 0:
                for i in tasks_preparsing:
                    remaining_no_tests = i[1]
                    task_size = i[2]
                    if remaining_no_tests - task_size <= 0:
                        tasks_preparsing.remove(i)
                        task_size = remaining_no_tests
                    else:
                        i[1] = remaining_no_tests - task_size
                    task = (i[0], task_size)
                    tasks.append(task)
            return tasks

        def register_results(results, global_no_tests_dict):
            r"""
            Registers the experimental results computed
            by a thread into ``self``'s dictionary of
            experimental data.
            """
            data = self.experimental_data()
            if verbosity:
                no_tasks = len(tasks)
                cur_task = 1
            for (task, local_results) in results:
                if verbosity:
                    print "Finished task %s of %s" % (cur_task, no_tasks)
                    cur_task += 1
                bench = task[0][0]
                no_test = global_no_test_dict[bench]
                for res in local_results:
                    data[bench, no_test] = res
                    no_test += 1
                global_no_test_dict[bench] = no_test

        #Setting variables
        setup_phase_tests = 4
        tasks_preparsing= []
        tasks = []
        global_no_test_dict = dict()
        if verbosity:
            print "Starting setup run"

        #Setup run phase
        for b in self.identifier():
            if self.code(b) is None:
                pass
            else:
                already_performed_tests = len(self.experimental_data(b))
                total_tests_to_perform = self.number_of_tests(b)
                remaining_tests = total_tests_to_perform - already_performed_tests
                if remaining_tests <= setup_phase_tests and remaining_tests > 0:
                    tasks.append((b, remaining_tests))
                    global_no_test_dict[b] = remaining_tests
                elif remaining_tests <= 0:
                    pass
                else:
                    global_no_test_dict[b] = already_performed_tests
                    tasks_preparsing.append([b, remaining_tests-setup_phase_tests, setup_phase_tests])
                    tasks.append((b, setup_phase_tests))

        results = self._perform_parallel_experiments_for_single_id(tasks)
        register_results(results, global_no_test_dict)

        if verbosity:
            print "Starting main run"
        #Main run

        #Adjusting chunk size
        target_task_time = 3.0
        for t in tasks_preparsing:
            bench = t[0]
            median_time = self.compute_timings("decoding",
                    median, identifier = bench)
            worst_time = self.experimental_data()[(b, 0)][4]

            ind = tasks_preparsing.index(t)
            tasks_preparsing[ind][2] = ceil(target_task_time/median_time)
            #ceil(worst_time) + ceil(target_task_time/median_time)

        tasks = preparse_tasks(tasks_preparsing)
        results = self._perform_parallel_experiments_for_single_id(tasks)
        register_results(results, global_no_test_dict)

    def _perform_experiments_for_single_id(self, identifier, verbosity):
        r"""
        Performs rounds of experimental data computation using ``self``'s
        parameters.

        See :meth:`run` for details.

        INPUT:

        - ``verbosity_level`` -- a number which indicates how verbose this run will
          be:
            - ``0`` means completely silent, nothing will be written on the output
            - ``1`` means small verbosity, one message will be written every time 25% of the tests are done
            - ``2`` means very verbose, one message will be written every time 10% of the tests are done

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark on which experiments
          will be performed.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B._perform_experiments_for_single_id('_0', 0)
        """
        #setting local variables and checking validity of data to use
        data = self.experimental_data(identifier = identifier)
        no_tests = self.number_of_tests(identifier = identifier)
        cur_no_tests = len(data)
        if cur_no_tests == no_tests:
            print "Experimental data already exist for benchmark %s.\nIf you want to perform new experiments, please run clear_experimental_data first" % identifier
            return
        if cur_no_tests == 0:
            data = self._experimental_data
        C = self.code(identifier)
        if C is None:
            print "Benchmark %s is empty, skipping computation" % identifier
            return
        D = self.decoder(identifier)
        Chan = self.channel(identifier)
        previous = 0 #used for verbosity
        if verbosity != 0:
            print "Starting run for benchmark %s" % identifier
        while(cur_no_tests < no_tests):
            start = time.clock()
            c = C.random_element()
            generate_c_time = time.clock() - start
            y = Chan(c)
            decoding_error = True
            try:
                start = time.clock()
                dec_y = D.decode_to_code(y)
                decode_time = time.clock() - start
                decoding_failure = False
                if "list-decoder" in D.decoder_type():
                    decoding_error = (c != dec_y)
                else:
                    decoding_error = (c != dec_y)
            except DecodingError:
                decode_time = time.clock() - start
                dec_y = None #as decoding failed
                decoding_failure = True
            data[identifier, cur_no_tests] = (c, generate_c_time, y, dec_y, decode_time, decoding_error, decoding_failure)

            #verbosity
            if verbosity ==1:
                cur = (cur_no_tests*4) // no_tests
                if (cur_no_tests*4) // no_tests != previous:
                    previous = cur
                    print "Benchmark %s: %s percent complete" % (identifier, cur * 25)
            if verbosity ==2:
                cur = (cur_no_tests*10) // no_tests
                if (cur_no_tests*10) // no_tests != previous:
                    previous = cur
                    print "Benchmark %s: %s percent complete" % (identifier, cur * 10)
            cur_no_tests += 1

        if verbosity != 0:
            print "Run complete for benchmark %s" % identifier

    def run(self, verbosity_level = 0, parallel = False):
        r"""
        Runs rounds of experimental data computation using ``self``'s
        parameters. The user can be informed about the current run's status by setting
        ``verbosity_level`` to the appropriate value.

        A round of experimental data computation consists in:

        - generating a random codeword in ``self``'s code,
        - using the ``self``'s channel to scramble ``c``, giving a new word ``y``, and
        - decoding ``y`` into ``y_dec`` using ``self``'s decoder.

        INPUT:

        - ``verbosity_level`` -- (default: ``0``) a number which indicates how verbose this run will
          be:

            - ``0`` means completely silent, nothing will be written on the output,
            - ``1`` means small verbosity, one message will be written every time 25% of the tests are done,
            - ``2`` means very verbose, one message will be written every time 10% of the tests are done.


        .. NOTE::

            This methods has a safeguard to prevent erasing experimental data by accident:
            if :meth:`experimental_data` is not empty, nothing will be done.
            See :meth:`clear_experimental_data` to delete data done in previous runs.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.run()

        Another one in verbose mode::

            sage: B.clear_experimental_data()
            sage: B.run(1)
            Starting run for benchmark _0
            Benchmark  _0: 25 percent complete
            Benchmark  _0: 50 percent complete
            Benchmark  _0: 75 percent complete
            Run complete for benchmark _0
        """
        if not verbosity_level in {0,1,2}:
            raise ValueError("verbosity_level has to be either 0, 1 or 2")
        if parallel == True:
            if verbosity_level != 0:
                verbosity_level = True
            else:
                verbosity_level = False
            self._task_master(verbosity_level)
        else:
            for i in self.identifier():
                self._perform_experiments_for_single_id(i, verbosity_level)

    def _process_experimental_data(self, identifier, target, processing):
        r"""
        Returns processed ``target``'s data using ``processing`` for ``self``
        designated by ``identifier``.

        INPUT:

        - ``identifier`` -- the identifier of the benchmark
          whose success rate will be computed.

        - ``target`` -- the operation whose timings will be computed, can be
          either ``"codeword_generation"`` or ``"decoding"``

        - ``processing`` -- the processing to apply on ``target``'s timings.
          It has to be a Sage/Python method. It can be either a native one
          (e.g. ``mean`` or ``median``) or one implemented by the user.
          If it's a custom one, it has to take a Python list as input.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.run()
            sage: B._process_experimental_data('_0', "decoding", median) #random
            0.0008940000000000059
        """
        if target is "codeword_generation":
            pos = 1
        elif target is "decoding":
            pos = 4
        else:
            raise ValueError("target has to be set to either \"codeword_generation\" or \"decoding\"")
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
          either ``"codeword_generation"`` or ``"decoding"``

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
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.run()
            sage: B.compute_timings("decoding", median) #random
            {'_0', 0.0008940000000000059}

        It is also possible to use your own processing method::

            sage: def my_silly_processing(l):
            ....:    return -1
            sage: B.compute_timings("decoding", my_silly_processing)
            {'_0': -1}
        """
        if not target in ("decoding", "codeword_generation"):
            raise ValueError("target has to be set to either \"codeword_generation\" or \"decoding\"")
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
    def decoding_success_rate(self, target, identifier = None):
        r"""
        Returns the success rate of the decoding attempts.

        INPUT:

        - ``target`` -- the operation whose success rate will be computed, can be
          either:

            - ``"failure``, which will compute a ratio on the number of decoding attempts
              which returned a word, whichever word it was
            - ``"error"``, which will compute a ratio on the number of decoding attempts
              which returned the original codeword

        - ``identifier`` -- (default: ``None``) the identifier of the benchmark
          whose success rate will be computed. If it is let to ``None``, success
          rate of every subbenchmark will be computed.

        OUTPUT:

        - a dictionary whose keys are the benchmarks' identifiers and values the associated
          success rates if ``identifier`` is ``None``, or

        - the success rate of the benchmark designated by ``identifier`` if it isn't ``None``.

        EXAMPLES::

            sage: C = codes.GeneralizedReedSolomonCode(GF(59).list()[:40], 12)
            sage: D = C.decoder()
            sage: Chan = channels.StaticErrorRateChannel(C.ambient_space(), D.decoding_radius())
            sage: B = codes.Benchmark(C, D, Chan)
            sage: B.run()
            sage: B.decoding_success_rate("failure")
            {'_0': 1.00}
        """
        if not target in ("error", "failure"):
            raise ValueError("target has to be set to either \"error\" or \"failure\"")
        t = (5 if target is "error" else 6)
        def compute_success_rate(identifier, target):
            i = 0
            success_sum = 0
            while True:
                try:
                    success_sum += not(data[identifier, i][target])
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
                rates[i] = compute_success_rate(i, t)
            return rates
        elif identifier is not None:
            return compute_success_rate(i, t)

    def plot(self, points_generator, **kwargs):
        r"""
        Returns a plot of points generator by ``points_generator``.

        INPUT:

        - ``points_generator`` -- a method which takes an identifier in input
          and returns a tuple

        - ``**kwargs`` -- all extra arguments are forwarded to :meth:`sage.plot.plot.list_plot`

        EXAMPLES::

            sage: C1 = codes.GeneralizedReedSolomonCode(GF(251).list()[:200], 125)
            sage: C2 = codes.GeneralizedReedSolomonCode(GF(251).list()[:200], 150)
            sage: C3 = codes.GeneralizedReedSolomonCode(GF(251).list()[:200], 175)
            sage: D1 = C1.decoder()
            sage: D2 = C2.decoder()
            sage: D3 = C3.decoder()
            sage: Chan1 = channels.StaticErrorRateChannel(C1.ambient_space(), D1.decoding_radius())
            sage: Chan2 = channels.StaticErrorRateChannel(C2.ambient_space(), D2.decoding_radius())
            sage: Chan3 = channels.StaticErrorRateChannel(C3.ambient_space(), D3.decoding_radius())
            sage: B1 = codes.Benchmark(C1, D1, Chan1)
            sage: B2 = codes.Benchmark(C2, D2, Chan2)
            sage: B3 = codes.Benchmark(C3, D3, Chan3)
            sage: B4 = B1 + B2 + B3
            sage: B4.run()
            sage: B4.plot(lambda i: (B4.decoder(i).decoding_radius(), B4.compute_timings("decoding", median, i)))
            Graphics object consisting of 1 graphics primitive
        """
        if not hasattr(points_generator, '__call__'):
            raise ValueError("points_generator has to be a method")
        l = [points_generator(i) for i in self.identifier()]
        return list_plot(l, **kwargs)
