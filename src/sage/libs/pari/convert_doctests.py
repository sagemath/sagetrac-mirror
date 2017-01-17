"""
This file is a direct copy of the content of the docstrings of the file
convert.pyx that used to be part of the PARI interface in Sage. After
:trac:`20238` the PARI interface will be distributed as an independent
module but the integrality of the doctests are kept for consistency
and testing.
"""
# doctest from convert
######################
# sage.libs.cypari2.convert
r"""Convert PARI objects to/from Python/C native types

This modules contains the following conversion routines:

- integers, long integers <-> PARI integers
- list of integegers -> PARI polynomials
- doubles -> PARI reals
- pairs of doubles -> PARI complex numbers

PARI integers are stored as an array of limbs of type ``pari_ulong``
(which are 32-bit or 64-bit integers). Depending on the kernel
(GMP or native), this array is stored little-endian or big-endian.
This is encapsulated in macros like ``int_W()``:
see section 4.5.1 of the
`PARI library manual <http://pari.math.u-bordeaux.fr/pub/pari/manuals/2.7.0/libpari.pdf>`_.

Python integers of type ``int`` are just C longs. Python integers of
type ``long`` are stored as a little-endian array of type ``digit``
with 15 or 30 bits used per digit. The internal format of a ``long`` is
not documented, but there is some information in
`longintrepr.h <https://github.com/python-git/python/blob/master/Include/longintrepr.h>`_.

Because of this difference in bit lengths, converting integers involves
some bit shuffling.
"""
# sage.libs.cypari2.convert.integer_to_gen
r"""Convert a Python ``int`` or ``long`` to a PARI ``gen`` of type
``t_INT``.

EXAMPLES::

    sage: from sage.libs.cypari2.convert import integer_to_gen
    sage: a = integer_to_gen(int(12345)); a; type(a)
    12345
    <type 'sage.libs.cypari2.gen.gen'>
    sage: a = integer_to_gen(long(12345)); a; type(a)
    12345
    <type 'sage.libs.cypari2.gen.gen'>
    sage: integer_to_gen(float(12345))
    Traceback (most recent call last):
    ...
    TypeError: integer_to_gen() needs an int or long argument, not float

TESTS::

    sage: for i in range(10000):
    ....:     x = 3**i
    ....:     if pari(long(x)) != pari(x):
    ....:         print(x)
"""
# sage.libs.cypari2.convert.gen_to_integer
r"""Convert a PARI ``gen`` to a Python ``int`` or ``long``.

INPUT:

- ``x`` -- a PARI ``t_INT``, ``t_FRAC``, ``t_REAL``, a purely
  real ``t_COMPLEX``, a ``t_INTMOD`` or ``t_PADIC`` (which are
  lifted).

EXAMPLES::

    sage: from sage.libs.cypari2.convert import gen_to_integer
    sage: a = gen_to_integer(pari("12345")); a; type(a)
    12345
    <type 'int'>
    sage: a = gen_to_integer(pari("10^30")); a; type(a)
    1000000000000000000000000000000L
    <type 'long'>
    sage: gen_to_integer(pari("19/5"))
    3
    sage: gen_to_integer(pari("1 + 0.0*I"))
    1
    sage: gen_to_integer(pari("3/2 + 0.0*I"))
    1
    sage: gen_to_integer(pari("Mod(-1, 11)"))
    10
    sage: gen_to_integer(pari("5 + O(5^10)"))
    5
    sage: gen_to_integer(pari("Pol(42)"))
    42
    sage: gen_to_integer(pari("x"))
    Traceback (most recent call last):
    ...
    TypeError: unable to convert PARI object x of type t_POL to an integer
    sage: gen_to_integer(pari("x + O(x^2)"))
    Traceback (most recent call last):
    ...
    TypeError: unable to convert PARI object x + O(x^2) of type t_SER to an integer
    sage: gen_to_integer(pari("1 + I"))
    Traceback (most recent call last):
    ...
    TypeError: unable to convert PARI object 1 + I of type t_COMPLEX to an integer

TESTS::

    sage: for i in range(10000):
    ....:     x = 3**i
    ....:     if long(pari(x)) != long(x):
    ....:         print(x)
    sage: gen_to_integer(pari("1.0 - 2^64"))
    -18446744073709551615L
    sage: gen_to_integer(pari("1 - 2^64"))
    -18446744073709551615L

Check some corner cases::

    sage: for s in [1, -1]:
    ....:     for a in [1, 2^31, 2^32, 2^63, 2^64]:
    ....:         for b in [-1, 0, 1]:
    ....:             Nstr = str(s * (a + b))
    ....:             N1 = gen_to_integer(pari(Nstr))  # Convert via PARI
    ....:             N2 = int(Nstr)                   # Convert via Python
    ....:             if N1 != N2:
    ....:                 print(Nstr, N1, N2)
    ....:             if type(N1) is not type(N2):
    ....:                 print(N1, type(N1), N2, type(N2))
"""
# sage.libs.cypari2.convert.NULL
r"""Convert a PARI object to a PARI integer.

This function is shallow and not stack-clean.
"""
# sage.libs.cypari2.convert.new_t_POL_from_int_star
r"""Note that degree + 1 = length, so that recognizing 0 is easier.

varnum = 0 is the general choice (creates a variable in x).
"""
