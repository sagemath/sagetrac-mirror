"""
Interface between Sage and PARI

Guide to real precision in the PARI interface
=============================================

In the PARI interface, "real precision" refers to the precision of real
numbers, so it is the floating-point precision. This is a non-trivial
issue, since there are various interfaces for different things.

Internal representation and conversion between Sage and PARI
------------------------------------------------------------

Real numbers in PARI have a precision associated to them, which is
always a multiple of the CPU wordsize. So, it is a multiple of 32
of 64 bits. When converting from Sage to PARI, the precision is rounded
up to the nearest multiple of the wordsize::

    sage: x = 1.0
    sage: x.precision()
    53
    sage: pari(x)
    1.00000000000000
    sage: pari(x).bitprecision()
    64

With a higher precision::

    sage: x = RealField(100).pi()
    sage: x.precision()
    100
    sage: pari(x).bitprecision()
    128

When converting back to Sage, the precision from PARI is taken::

    sage: x = RealField(100).pi()
    sage: y = pari(x).sage()
    sage: y
    3.1415926535897932384626433832793333156
    sage: parent(y)
    Real Field with 128 bits of precision

So ``pari(x).sage()`` is definitely not equal to ``x`` since it has
28 bogus bits.

Therefore, some care must be taken when juggling reals back and forth
between Sage and PARI. The correct way of avoiding this is to convert
``pari(x).sage()`` back into a domain with the right precision. This has
to be done by the user (or by Sage functions that use PARI library
functions). For instance, if we want to use the PARI library to compute
``sqrt(pi)`` with a precision of 100 bits::

    sage: R = RealField(100)
    sage: s = R(pi); s
    3.1415926535897932384626433833
    sage: p = pari(s).sqrt()
    sage: x = p.sage(); x    # wow, more digits than I expected!
    1.7724538509055160272981674833410973484
    sage: x.prec()           # has precision 'improved' from 100 to 128?
    128
    sage: x == RealField(128)(pi).sqrt()  # sadly, no!
    False
    sage: R(x)               # x should be brought back to precision 100
    1.7724538509055160272981674833
    sage: R(x) == s.sqrt()
    True

Output precision for printing
-----------------------------

Even though PARI reals have a precision, not all significant bits are
printed by default. The maximum number of digits when printing a PARI
real can be set using the methods
:meth:`Pari.set_real_precision_bits` or
:meth:`Pari.set_real_precision`.

We create a very precise approximation of pi and see how it is printed
in PARI::

    sage: pi = pari(RealField(1000).pi())

The default precision is 15 digits::

    sage: pi
    3.14159265358979

With a different precision::

    sage: _ = pari.set_real_precision(50)
    sage: pi
    3.1415926535897932384626433832795028841971693993751

Back to the default::

    sage: _ = pari.set_real_precision(15)
    sage: pi
    3.14159265358979

Input precision for function calls
----------------------------------

When we talk about precision for PARI functions, we need to distinguish
three kinds of calls:

1. Using the string interface, for example ``pari("sin(1)")``.

2. Using the library interface with exact inputs, for example
   ``pari(1).sin()``.

3. Using the library interface with inexact inputs, for example
   ``pari(1.0).sin()``.

In the first case, the relevant precision is the one set by the methods
:meth:`Pari.set_real_precision_bits` or
:meth:`Pari.set_real_precision`::

    sage: pari.set_real_precision_bits(150)
    sage: pari("sin(1)")
    0.841470984807896506652502321630298999622563061
    sage: pari.set_real_precision_bits(53)
    sage: pari("sin(1)")
    0.841470984807897

In the second case, the precision can be given as the argument
``precision`` in the function call, with a default of 53 bits.
The real precision set by
:meth:`Pari.set_real_precision_bits` or
:meth:`Pari.set_real_precision` is irrelevant.

In these examples, we convert to Sage to ensure that PARI's real
precision is not used when printing the numbers. As explained before,
this artificially increases the precision to a multiple of the
wordsize. ::

    sage: s = pari(1).sin(precision=180).sage(); print(s); print(parent(s))
    0.841470984807896506652502321630298999622563060798371065673
    Real Field with 192 bits of precision
    sage: s = pari(1).sin(precision=40).sage(); print(s); print(parent(s))
    0.841470984807896507
    Real Field with 64 bits of precision
    sage: s = pari(1).sin().sage(); print(s); print(parent(s))
    0.841470984807896507
    Real Field with 64 bits of precision

In the third case, the precision is determined only by the inexact
inputs and the ``precision`` argument is ignored::

    sage: pari(1.0).sin(precision=180).sage()
    0.841470984807896507
    sage: pari(1.0).sin(precision=40).sage()
    0.841470984807896507
    sage: pari(RealField(100).one()).sin().sage()
    0.84147098480789650665250232163029899962

Elliptic curve functions
------------------------

An elliptic curve given with exact `a`-invariants is considered an
exact object. Therefore, you should set the precision for each method
call individually::

    sage: e = pari([0,0,0,-82,0]).ellinit()
    sage: eta1 = e.elleta(precision=100)[0]
    sage: eta1.sage()
    3.6054636014326520859158205642077267748
    sage: eta1 = e.elleta(precision=180)[0]
    sage: eta1.sage()
    3.60546360143265208591582056420772677481026899659802474544

"""

def _get_pari_instance():
    # There are two constraints for the virtual stack size:
    # 1) on 32-bit systems, even virtual memory can be a scarce
    #    resource since it is limited by 4GB (of which the kernel
    #    needs a significant part)
    # 2) the system should actually be able to handle a stack size
    #    as large as the complete virtual stack.
    # As a simple heuristic, we set the virtual stack to 1/4 of the
    # virtual memory.
    from sage.misc.getusage import virtual_memory_limit

    sizemax = virtual_memory_limit() // 4

    from sage.env import CYGWIN_VERSION
    if CYGWIN_VERSION and CYGWIN_VERSION < (2, 5, 2):
        # Cygwin's mmap is broken for large NORESERVE mmaps (>~ 4GB) See
        # http://trac.sagemath.org/ticket/20463 So we set the max stack
        # size to a little below 4GB (putting it right on the margin proves
        # too fragile)
        #
        # The underlying issue is fixed in Cygwin v2.5.2
        sizemax = min(sizemax, 0xf0000000)

    from cypari2 import Pari
    P = Pari(1000000, sizemax)

    # pari_init_opts() overrides MPIR's memory allocation functions,
    # so we need to reset them.
    from sage.ext.memory import init_memory_functions
    init_memory_functions()

    return P

pari = _get_pari_instance()
