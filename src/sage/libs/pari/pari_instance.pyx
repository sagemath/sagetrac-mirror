"""
PARI C-library interface

AUTHORS:

- William Stein (2006-03-01): updated to work with PARI 2.2.12-beta

- William Stein (2006-03-06): added newtonpoly

- Justin Walker: contributed some of the function definitions

- Gonzalo Tornaria: improvements to conversions; much better error
  handling.

- Robert Bradshaw, Jeroen Demeyer, William Stein (2010-08-15):
  Upgrade to PARI 2.4.3 (:trac:`9343`)

- Jeroen Demeyer (2011-11-12): rewrite various conversion routines
  (:trac:`11611`, :trac:`11854`, :trac:`11952`)

- Peter Bruin (2013-11-17): split off this file from gen.pyx
  (:trac:`15185`)

- Jeroen Demeyer (2014-02-09): upgrade to PARI 2.7 (:trac:`15767`)

- Jeroen Demeyer (2014-09-19): upgrade to PARI 2.8 (:trac:`16997`)

- Jeroen Demeyer (2015-03-17): automatically generate methods from
  ``pari.desc`` (:trac:`17631` and :trac:`17860`)

EXAMPLES::

    sage: pari('5! + 10/x')
    (120*x + 10)/x
    sage: pari('intnum(x=0,13,sin(x)+sin(x^2) + x)')
    83.8179442684285  # 32-bit
    84.1818153922297  # 64-bit
    sage: f = pari('x^3-1')
    sage: v = f.factor(); v
    [x - 1, 1; x^2 + x + 1, 1]
    sage: v[0]   # indexing is 0-based unlike in GP.
    [x - 1, x^2 + x + 1]~
    sage: v[1]
    [1, 1]~

Arithmetic obeys the usual coercion rules::

    sage: type(pari(1) + 1)
    <type 'sage.libs.pari.gen.gen'>
    sage: type(1 + pari(1))
    <type 'sage.libs.pari.gen.gen'>

GUIDE TO REAL PRECISION AND THE PARI LIBRARY

The default real precision in communicating with the PARI library
is the same as the default Sage real precision, which is 53 bits.
Inexact Pari objects are therefore printed by default to 15 decimal
digits (even if they are actually more precise).

Default precision example (53 bits, 15 significant decimals)::

    sage: a = pari(1.23); a
    1.23000000000000
    sage: a.sin()
    0.942488801931698

Example with custom precision of 200 bits (60 significant
decimals)::

    sage: R = RealField(200)
    sage: a = pari(R(1.23)); a   # only 15 significant digits printed
    1.23000000000000
    sage: R(a)         # but the number is known to precision of 200 bits
    1.2300000000000000000000000000000000000000000000000000000000
    sage: a.sin()      # only 15 significant digits printed
    0.942488801931698
    sage: R(a.sin())   # but the number is known to precision of 200 bits
    0.94248880193169751002382356538924454146128740562765030213504

It is possible to change the number of printed decimals::

    sage: R = RealField(200)    # 200 bits of precision in computations
    sage: old_prec = pari.set_real_precision(60)  # 60 decimals printed
    sage: a = pari(R(1.23)); a
    1.23000000000000000000000000000000000000000000000000000000000
    sage: a.sin()
    0.942488801931697510023823565389244541461287405627650302135038
    sage: pari.set_real_precision(old_prec)  # restore the default printing behavior
    60

Unless otherwise indicated in the docstring, most Pari functions
that return inexact objects use the precision of their arguments to
decide the precision of the computation. However, if some of these
arguments happen to be exact numbers (integers, rationals, etc.),
an optional parameter indicates the precision (in bits) to which
these arguments should be converted before the computation. If this
precision parameter is missing, the default precision of 53 bits is
used. The following first converts 2 into a real with 53-bit
precision::

    sage: R = RealField()
    sage: R(pari(2).sin())
    0.909297426825682

We can ask for a better precision using the optional parameter::

    sage: R = RealField(150)
    sage: R(pari(2).sin(precision=150))
    0.90929742682568169539601986591174484270225497

Warning regarding conversions Sage - Pari - Sage: Some care must be
taken when juggling inexact types back and forth between Sage and
Pari. In theory, calling p=pari(s) creates a Pari object p with the
same precision as s; in practice, the Pari library's precision is
word-based, so it will go up to the next word. For example, a
default 53-bit Sage real s will be bumped up to 64 bits by adding
bogus 11 bits. The function p.python() returns a Sage object with
exactly the same precision as the Pari object p. So
pari(s).python() is definitely not equal to s, since it has 64 bits
of precision, including the bogus 11 bits. The correct way of
avoiding this is to coerce pari(s).python() back into a domain with
the right precision. This has to be done by the user (or by Sage
functions that use Pari library functions in gen.pyx). For
instance, if we want to use the Pari library to compute sqrt(pi)
with a precision of 100 bits::

    sage: R = RealField(100)
    sage: s = R(pi); s
    3.1415926535897932384626433833
    sage: p = pari(s).sqrt()
    sage: x = p.python(); x  # wow, more digits than I expected!
    1.7724538509055160272981674833410973484
    sage: x.prec()           # has precision 'improved' from 100 to 128?
    128
    sage: x == RealField(128)(pi).sqrt()  # sadly, no!
    False
    sage: R(x)               # x should be brought back to precision 100
    1.7724538509055160272981674833
    sage: R(x) == s.sqrt()
    True

Elliptic curves and precision: If you are working with elliptic
curves, you should set the precision for each method::

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: eta1 = e.elleta(precision=100)[0]
    sage: eta1.sage()
    3.6054636014326520859158205642077267748
    sage: eta1 = e.elleta(precision=180)[0]
    sage: eta1.sage()
    3.60546360143265208591582056420772677481026899659802474544

Number fields and precision: TODO

TESTS:

Check that output from PARI's print command is actually seen by
Sage (:trac:`9636`)::

    sage: pari('print("test")')
    test

"""

#*****************************************************************************
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from .cypari_instance cimport PariInstance

# from .paridecl cimport *
# from .paripriv cimport *
include "cysignals/signals.pxi"
# cdef extern from *:
#     int sig_on_count "cysigs.sig_on_count"

# import sys

# cimport libc.stdlib
# from libc.stdio cimport *
cimport cython

# include "cysignals/memory.pxi"
from sage.ext.memory import init_memory_functions
from sage.structure.parent cimport Parent
from sage.libs.gmp.all cimport *
from sage.libs.flint.fmpz cimport fmpz_get_mpz, COEFF_IS_MPZ, COEFF_TO_PTR
from sage.libs.flint.fmpz_mat cimport *

from .gen cimport gen, objtogen
# from .handle_error cimport _pari_init_error_handling
from sage.misc.superseded import deprecation, deprecated_function_alias

# # real precision in decimal digits: see documentation for
# # get_real_precision() and set_real_precision().  This variable is used
# # in gp to set the precision of input quantities (e.g. sqrt(2)), and for
# # determining the number of digits to be printed.  It is *not* used as
# # a "default precision" for internal computations, which always use
# # the actual precision of arguments together (where relevant) with a
# # "prec" parameter.  In ALL cases (for real computations) the prec
# # parameter is a WORD precision and NOT decimal precision.  Pari reals
# # with word precision w have bit precision (of the mantissa) equal to
# # 32*(w-2) or 64*(w-2).
# #
# # Hence the only relevance of this parameter in Sage is (1) for the
# # output format of components of objects of type
# # 'sage.libs.pari.gen.gen'; (2) for setting the precision of pari
# # variables created from strings (e.g. via sage: pari('1.2')).
# #
# # WARNING: Many pari library functions take a last parameter "prec"
# # which should be a words precision.  In many cases this is redundant
# # and is simply ignored.  In our wrapping of these functions we use
# # the variable prec here for convenience only.
# cdef long prec

# #################################################################
# # conversions between various real precision models
# #################################################################

# def prec_bits_to_dec(unsigned long prec_in_bits):
#     r"""
#     Convert from precision expressed in bits to precision expressed in
#     decimal.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_bits_to_dec
#         sage: prec_bits_to_dec(53)
#         15
#         sage: [(32*n, prec_bits_to_dec(32*n)) for n in range(1, 9)]
#         [(32, 9),
#         (64, 19),
#         (96, 28),
#         (128, 38),
#         (160, 48),
#         (192, 57),
#         (224, 67),
#         (256, 77)]
#     """
#     cdef double log_2 = 0.301029995663981
#     return int(prec_in_bits*log_2)

# def prec_dec_to_bits(unsigned long prec_in_dec):
#     r"""
#     Convert from precision expressed in decimal to precision expressed
#     in bits.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_dec_to_bits
#         sage: prec_dec_to_bits(15)
#         49
#         sage: [(n, prec_dec_to_bits(n)) for n in range(10, 100, 10)]
#         [(10, 33),
#         (20, 66),
#         (30, 99),
#         (40, 132),
#         (50, 166),
#         (60, 199),
#         (70, 232),
#         (80, 265),
#         (90, 298)]
#     """
#     cdef double log_10 = 3.32192809488736
#     return int(prec_in_dec*log_10)

# cpdef long prec_bits_to_words(unsigned long prec_in_bits):
#     r"""
#     Convert from precision expressed in bits to pari real precision
#     expressed in words. Note: this rounds up to the nearest word,
#     adjusts for the two codewords of a pari real, and is
#     architecture-dependent.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_bits_to_words
#         sage: prec_bits_to_words(70)
#         5   # 32-bit
#         4   # 64-bit

#     ::

#         sage: [(32*n, prec_bits_to_words(32*n)) for n in range(1, 9)]
#         [(32, 3), (64, 4), (96, 5), (128, 6), (160, 7), (192, 8), (224, 9), (256, 10)] # 32-bit
#         [(32, 3), (64, 3), (96, 4), (128, 4), (160, 5), (192, 5), (224, 6), (256, 6)] # 64-bit
#     """
#     if not prec_in_bits:
#         return prec
#     cdef unsigned long wordsize = BITS_IN_LONG

#     # This equals ceil(prec_in_bits/wordsize) + 2
#     return (prec_in_bits - 1)//wordsize + 3

# cpdef long prec_words_to_bits(long prec_in_words):
#     r"""
#     Convert from pari real precision expressed in words to precision
#     expressed in bits. Note: this adjusts for the two codewords of a
#     pari real, and is architecture-dependent.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_words_to_bits
#         sage: prec_words_to_bits(10)
#         256   # 32-bit
#         512   # 64-bit
#         sage: [(n, prec_words_to_bits(n)) for n in range(3, 10)]
#         [(3, 32), (4, 64), (5, 96), (6, 128), (7, 160), (8, 192), (9, 224)]  # 32-bit
#         [(3, 64), (4, 128), (5, 192), (6, 256), (7, 320), (8, 384), (9, 448)] # 64-bit
#     """
#     # see user's guide to the pari library, page 10
#     return (prec_in_words - 2) * BITS_IN_LONG

# cpdef long default_bitprec():
#     r"""
#     Return the default precision in bits.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import default_bitprec
#         sage: default_bitprec()
#         64
#     """
#     return (prec - 2) * BITS_IN_LONG

# def prec_dec_to_words(long prec_in_dec):
#     r"""
#     Convert from precision expressed in decimal to precision expressed
#     in words. Note: this rounds up to the nearest word, adjusts for the
#     two codewords of a pari real, and is architecture-dependent.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_dec_to_words
#         sage: prec_dec_to_words(38)
#         6   # 32-bit
#         4   # 64-bit
#         sage: [(n, prec_dec_to_words(n)) for n in range(10, 90, 10)]
#         [(10, 4), (20, 5), (30, 6), (40, 7), (50, 8), (60, 9), (70, 10), (80, 11)] # 32-bit
#         [(10, 3), (20, 4), (30, 4), (40, 5), (50, 5), (60, 6), (70, 6), (80, 7)] # 64-bit
#     """
#     return prec_bits_to_words(prec_dec_to_bits(prec_in_dec))

# def prec_words_to_dec(long prec_in_words):
#     r"""
#     Convert from precision expressed in words to precision expressed in
#     decimal. Note: this adjusts for the two codewords of a pari real,
#     and is architecture-dependent.

#     EXAMPLES::

#         sage: from sage.libs.pari.pari_instance import prec_words_to_dec
#         sage: prec_words_to_dec(5)
#         28   # 32-bit
#         57   # 64-bit
#         sage: [(n, prec_words_to_dec(n)) for n in range(3, 10)]
#         [(3, 9), (4, 19), (5, 28), (6, 38), (7, 48), (8, 57), (9, 67)] # 32-bit
#         [(3, 19), (4, 38), (5, 57), (6, 77), (7, 96), (8, 115), (9, 134)] # 64-bit
#     """
#     return prec_bits_to_dec(prec_words_to_bits(prec_in_words))


# The unique running Pari instance.
cdef SagePariInstance pari_instance, P
pari_instance = SagePariInstance()
P = pari_instance   # shorthand notation

# # PariInstance.__init__ must not create gen objects because their parent is not constructed yet
# sig_on()
# pari_instance.PARI_ZERO = pari_instance.new_gen_noclear(gen_0)
# pari_instance.PARI_ONE  = pari_instance.new_gen_noclear(gen_1)
# pari_instance.PARI_TWO  = pari_instance.new_gen_noclear(gen_2)
# sig_off()

# Also a copy of PARI accessible from external pure python code.
pari = pari_instance


# # Callbacks from PARI to print stuff using sys.stdout.write() instead
# # of C library functions like puts().
# cdef PariOUT sage_pariOut

# cdef void sage_putchar(char c):
#     cdef char s[2]
#     s[0] = c
#     s[1] = 0
#     sys.stdout.write(s)
#     # Let PARI think the last character was a newline,
#     # so it doesn't print one when an error occurs.
#     pari_set_last_newline(1)

# cdef void sage_puts(char* s):
#     sys.stdout.write(s)
#     pari_set_last_newline(1)

# cdef void sage_flush():
#     sys.stdout.flush()


# include 'auto_instance.pxi'

@cython.final
cdef class SagePariInstance(PariInstance):
    def __init__(self, long size=1000000, long sizemax=0, unsigned long maxprime=500000):
        """
        Initialize the PARI system.

        INPUT:


        -  ``size`` -- long, the number of bytes for the initial
           PARI stack (see note below)

        -  ``sizemax`` -- long, the number of bytes for the virtual
           PARI stack, or 0 for default (see note below)
        
        -  ``maxprime`` -- unsigned long, upper limit on a
           precomputed prime number table (default: 500000)


        .. note::

           In Sage, the PARI stack is different than in GP or the
           PARI C library. In Sage, instead of the PARI stack
           holding the results of all computations, it *only* holds
           the results of an individual computation. Each time a new
           Python/PARI object is computed, it it copied to its own
           space in the Python heap, and the memory it occupied on the
           PARI stack is freed. Thus it is not necessary to make the
           stack very large. Also, unlike in PARI, if the stack does
           overflow, in most cases the PARI stack is automatically
           increased and the relevant step of the computation rerun.

           This design obviously involves some performance penalties
           over the way PARI works, but it scales much better and is
           far more robust for large projects.

        .. note::

           If you do not want prime numbers, put ``maxprime=2``, but be
           careful because many PARI functions require this table. If
           you get the error message "not enough precomputed primes",
           increase this parameter.
        """
        # PARI has a "real" stack size (parisize) and a "virtual" stack
        # size (parisizemax). The idea is that the real stack will be
        # used if possible, but the stack might be increased up to
        # the complete virtual stack. Therefore, it is not a problem to
        # set the virtual stack size to a large value. There are two
        # constraints for the virtual stack size:
        # 1) on 32-bit systems, even virtual memory can be a scarce
        #    resource since it is limited by 4GB (of which the kernel
        #    needs a significant part)
        # 2) the system should actually be able to handle a stack size
        #    as large as the complete virtual stack.
        # As a simple heuristic, we set the virtual stack to 1/4 of the
        # virtual memory.

        from sage.misc.memory_info import MemoryInfo
        mem = MemoryInfo()
        
        super(PariInstance, self).__init__(size,  mem.virtual_memory_limit() // 4, maxprime)
        
        # pari_init_opts() overrides MPIR's memory allocation functions,
        # so we need to reset them.
        init_memory_functions()

    cpdef _coerce_map_from_(self, x):
        """
        Return ``True`` if ``x`` admits a coercion map into the
        PARI interface.

        This currently always returns ``True``.

        EXAMPLES::

            sage: pari._coerce_map_from_(ZZ)
            True
            sage: pari.coerce_map_from(ZZ)
            Call morphism:
              From: Integer Ring
              To:   Interface to the PARI C library
        """
        return True

    # def __richcmp__(left, right, int op):
    #     """
    #     EXAMPLES::

    #         sage: pari == pari
    #         True
    #         sage: pari == gp
    #         False
    #         sage: pari == 5
    #         False
    #     """
    #     return (<Parent>left)._richcmp(right, op)

    cdef gen new_gen_from_mpz_t(self, mpz_t value):
        """
        Create a new gen from a given MPIR-integer ``value``.

        EXAMPLES::

            sage: pari(42)       # indirect doctest
            42

        TESTS:

        Check that the hash of an integer does not depend on existing
        garbage on the stack (:trac:`11611`)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        sig_on()
        return self.new_gen(self._new_GEN_from_mpz_t(value))

    cdef inline GEN _new_GEN_from_mpz_t(self, mpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``mpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef unsigned long limbs = mpz_size(value)

        cdef GEN z = cgeti(limbs + 2)
        # Set sign and "effective length"
        z[1] = evalsigne(mpz_sgn(value)) + evallgefint(limbs + 2)
        mpz_export(int_LSW(z), NULL, -1, sizeof(long), 0, 0, value)

        return z

    cdef inline GEN _new_GEN_from_fmpz_t(self, fmpz_t value):
        r"""
        Create a new PARI ``t_INT`` from a ``fmpz_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        if COEFF_IS_MPZ(value[0]):
            return self._new_GEN_from_mpz_t(COEFF_TO_PTR(value[0]))
        else:
            return stoi(value[0])

    cdef gen new_gen_from_int(self, int value):
        sig_on()
        return self.new_gen(stoi(value))

    cdef gen new_gen_from_mpq_t(self, mpq_t value):
        """
        Create a new gen from a given MPIR-rational ``value``.

        EXAMPLES::

            sage: pari(-2/3)
            -2/3
            sage: pari(QQ(42))
            42
            sage: pari(QQ(42)).type()
            't_INT'
            sage: pari(QQ(1/42)).type()
            't_FRAC'

        TESTS:

        Check that the hash of a rational does not depend on existing
        garbage on the stack (:trac:`11854`)::

            sage: foo = pari(2^(32*1024));  # Create large integer to put PARI stack in known state
            sage: a5 = pari(5/7);
            sage: foo = pari(0xDEADBEEF * (2^(32*1024)-1)//(2^32 - 1));  # Dirty PARI stack
            sage: b5 = pari(5/7);
            sage: a5.__hash__() == b5.__hash__()
            True
        """
        sig_on()
        return self.new_gen(self._new_GEN_from_mpq_t(value))

    cdef inline GEN _new_GEN_from_mpq_t(self, mpq_t value):
        r"""
        Create a new PARI ``t_INT`` or ``t_FRAC`` from a ``mpq_t``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN num = self._new_GEN_from_mpz_t(mpq_numref(value))
        if mpz_cmpabs_ui(mpq_denref(value), 1) == 0:
            # Denominator is 1, return the numerator (an integer)
            return num
        cdef GEN denom = self._new_GEN_from_mpz_t(mpq_denref(value))
        return mkfrac(num, denom)

    cdef gen new_t_POL_from_int_star(self, int *vals, int length, long varnum):
        """
        Note that degree + 1 = length, so that recognizing 0 is easier.

        varnum = 0 is the general choice (creates a variable in x).
        """
        cdef GEN z
        cdef int i

        sig_on()
        z = cgetg(length + 2, t_POL)
        z[1] = evalvarn(varnum)
        if length != 0:
            setsigne(z,1)
            for i from 0 <= i < length:
                set_gel(z,i+2, stoi(vals[i]))
        else:
            ## polynomial is zero
            setsigne(z,0)

        return self.new_gen(z)

    cdef gen new_gen_from_padic(self, long ordp, long relprec,
                                mpz_t prime, mpz_t p_pow, mpz_t unit):
        cdef GEN z
        sig_on()
        z = cgetg(5, t_PADIC)
        z[1] = evalprecp(relprec) + evalvalp(ordp)
        set_gel(z, 2, self._new_GEN_from_mpz_t(prime))
        set_gel(z, 3, self._new_GEN_from_mpz_t(p_pow))
        set_gel(z, 4, self._new_GEN_from_mpz_t(unit))
        return self.new_gen(z)

    def double_to_gen(self, x):
        cdef double dx
        dx = float(x)
        return self.double_to_gen_c(dx)

    cdef gen double_to_gen_c(self, double x):
        """
        Create a new gen with the value of the double x, using Pari's
        dbltor.

        EXAMPLES::

            sage: pari.double_to_gen(1)
            1.00000000000000
            sage: pari.double_to_gen(1e30)
            1.00000000000000 E30
            sage: pari.double_to_gen(0)
            0.E-15
            sage: pari.double_to_gen(-sqrt(RDF(2)))
            -1.41421356237310
        """
        # Pari has an odd concept where it attempts to track the accuracy
        # of floating-point 0; a floating-point zero might be 0.0e-20
        # (meaning roughly that it might represent any number in the
        # range -1e-20 <= x <= 1e20).

        # Pari's dbltor converts a floating-point 0 into the Pari real
        # 0.0e-307; Pari treats this as an extremely precise 0.  This
        # can cause problems; for instance, the Pari incgam() function can
        # be very slow if the first argument is very precise.

        # So we translate 0 into a floating-point 0 with 53 bits
        # of precision (that's the number of mantissa bits in an IEEE
        # double).

        sig_on()
        if x == 0:
            return self.new_gen(real_0_bit(-53))
        else:
            return self.new_gen(dbltor(x))

    cdef GEN double_to_GEN(self, double x):
        if x == 0:
            return real_0_bit(-53)
        else:
            return dbltor(x)

    # def complex(self, re, im):
    #     """
    #     Create a new complex number, initialized from re and im.
    #     """
    #     cdef gen t0 = self(re)
    #     cdef gen t1 = self(im)
    #     sig_on()
    #     return self.new_gen(mkcomplex(t0.g, t1.g))


    cdef GEN _new_GEN_from_fmpz_mat_t(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**``.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,j))
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef GEN _new_GEN_from_fmpz_mat_t_rotate90(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc):
        r"""
        Create a new PARI ``t_MAT`` with ``nr`` rows and ``nc`` columns
        from a ``mpz_t**`` and rotate the matrix 90 degrees
        counterclockwise.  So the resulting matrix will have ``nc`` rows
        and ``nr`` columns.  This is useful for computing the Hermite
        Normal Form because Sage and PARI use different definitions.

        For internal use only; this directly uses the PARI stack.
        One should call ``sig_on()`` before and ``sig_off()`` after.
        """
        cdef GEN x
        cdef GEN A = zeromatcopy(nc, nr)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_fmpz_t(fmpz_mat_entry(B,i,nc-j-1))
                set_gcoeff(A, j+1, i+1, x)  # A[j+1, i+1] = x (using 1-based indexing)
        return A

    cdef gen integer_matrix(self, fmpz_mat_t B, Py_ssize_t nr, Py_ssize_t nc, bint permute_for_hnf):
        """
        EXAMPLES::

            sage: matrix(ZZ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        sig_on()
        cdef GEN g
        if permute_for_hnf:
            g = self._new_GEN_from_fmpz_mat_t_rotate90(B, nr, nc)
        else:
            g = self._new_GEN_from_fmpz_mat_t(B, nr, nc)
        return self.new_gen(g)

    cdef GEN _new_GEN_from_mpq_t_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        cdef GEN x
        # Allocate zero matrix
        cdef GEN A = zeromatcopy(nr, nc)
        cdef Py_ssize_t i, j
        for i in range(nr):
            for j in range(nc):
                x = self._new_GEN_from_mpq_t(B[i][j])
                set_gcoeff(A, i+1, j+1, x)  # A[i+1, j+1] = x (using 1-based indexing)
        return A

    cdef gen rational_matrix(self, mpq_t** B, Py_ssize_t nr, Py_ssize_t nc):
        """
        EXAMPLES::

            sage: matrix(QQ,2,[1..6])._pari_()   # indirect doctest
            [1, 2, 3; 4, 5, 6]
        """
        sig_on()
        cdef GEN g = self._new_GEN_from_mpq_t_matrix(B, nr, nc)
        return self.new_gen(g)

    cdef _coerce_c_impl(self, x):
        """
        Implicit canonical coercion into a PARI object.
        """
        try:
            return self(x)
        except (TypeError, AttributeError):
            raise TypeError("no canonical coercion of %s into PARI" % x)

    cdef _an_element_c_impl(self):  # override this in Cython
        return self.PARI_ZERO

     def genus2red(self, P, P0=None):
        """
        Let `P` be a polynomial with integer coefficients.
        Determines the reduction of the (proper, smooth) genus 2
        curve `C/\QQ`, defined by the hyperelliptic equation `y^2 = P`.
        The special syntax ``genus2red([P,Q])`` is also allowed, where
        the polynomials `P` and `Q` have integer coefficients, to
        represent the model `y^2 + Q(x)y = P(x)`.

        EXAMPLES::

            sage: x = polygen(QQ)
            sage: pari.genus2red([-5*x^5, x^3 - 2*x^2 - 2*x + 1])
            [1416875, [2, -1; 5, 4; 2267, 1], x^6 - 240*x^4 - 2550*x^3 - 11400*x^2 - 24100*x - 19855, [[2, [2, [Mod(1, 2)]], []], [5, [1, []], ["[V] page 156", [3]]], [2267, [2, [Mod(432, 2267)]], ["[I{1-0-0}] page 170", []]]]]

        This is the old deprecated syntax::

            sage: pari.genus2red(x^3 - 2*x^2 - 2*x + 1, -5*x^5)
            doctest:...: DeprecationWarning: The 2-argument version of genus2red() is deprecated, use genus2red(P) or genus2red([P,Q]) instead
            See http://trac.sagemath.org/16997 for details.
            [1416875, [2, -1; 5, 4; 2267, 1], x^6 - 240*x^4 - 2550*x^3 - 11400*x^2 - 24100*x - 19855, [[2, [2, [Mod(1, 2)]], []], [5, [1, []], ["[V] page 156", [3]]], [2267, [2, [Mod(432, 2267)]], ["[I{1-0-0}] page 170", []]]]]
        """
        if P0 is not None:
            from sage.misc.superseded import deprecation
            deprecation(16997, 'The 2-argument version of genus2red() is deprecated, use genus2red(P) or genus2red([P,Q]) instead')
            P = [P0, P]
        cdef gen t0 = objtogen(P)
        sig_on()
        return self.new_gen(genus2red(t0.g, NULL))


cdef inline void INT_to_mpz(mpz_ptr value, GEN g):
    """
    Store a PARI ``t_INT`` as an ``mpz_t``.
    """
    if typ(g) != t_INT:
        pari_err(e_TYPE, <char*>"conversion to mpz", g)

    cdef long size = lgefint(g) - 2
    mpz_import(value, size, -1, sizeof(long), 0, 0, int_LSW(g))

    if signe(g) < 0:
        mpz_neg(value, value)

cdef void INTFRAC_to_mpq(mpq_ptr value, GEN g):
    """
    Store a PARI ``t_INT`` or ``t_FRAC`` as an ``mpq_t``.
    """
    if typ(g) == t_FRAC:
        INT_to_mpz(mpq_numref(value), gel(g, 1))
        INT_to_mpz(mpq_denref(value), gel(g, 2))
    elif typ(g) == t_INT:
        INT_to_mpz(mpq_numref(value), g)
        mpz_set_ui(mpq_denref(value), 1)
    else:
        pari_err(e_TYPE, <char*>"conversion to mpq", g)
