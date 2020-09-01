# -*- coding: utf-8 -*-
"""
The Sage Preparser

AUTHORS:

    - William Stein (2006-02-19)

      - Fixed bug when loading .py files.

    - William Stein (2006-03-09)

      - Fixed crash in parsing exponentials.
      - Precision of real literals now determined by digits of input
        (like Mathematica).

    - Joe Wetherell (2006-04-14)

      - Added MAGMA-style constructor preparsing.

    - Bobby Moretti (2007-01-25)

      - Added preliminary function assignment notation.

    - Robert Bradshaw (2007-09-19)

      - Added strip_string_literals, containing_block utility
        functions. Arrr!
      - Added [1,2,..,n] notation.

    - Robert Bradshaw (2008-01-04)

      - Implicit multiplication (off by default).

    - Robert Bradshaw (2008-09-23)

      - Factor out constants.

    - Robert Bradshaw (2000-01)

      - Simplify preparser by making it modular and using regular
        expressions.
      - Bug fixes, complex numbers, and binary input.

EXAMPLES:

Preparsing::

    sage: preparse('2/3')
    'Integer(2)/Integer(3)'
    sage: preparse('2.5')
    "RealNumber('2.5')"
    sage: preparse('2^3')
    'Integer(2)**Integer(3)'
    sage: preparse('a^b')            # exponent
    'a**b'
    sage: preparse('a**b')
    'a**b'
    sage: preparse('G.0')            # generator
    'G.gen(0)'
    sage: preparse('a = 939393R')    # raw
    'a = 939393'
    sage: implicit_multiplication(True)
    sage: preparse('a b c in L')     # implicit multiplication
    'a*b*c in L'
    sage: preparse('2e3x + 3exp(y)')
    "RealNumber('2e3')*x + Integer(3)*exp(y)"

A string with escaped quotes in it (the point here is that the
preparser doesn't get confused by the internal quotes)::

    sage: "\"Yes,\" he said."
    '"Yes," he said.'
    sage: s = "\\"; s
    '\\'

A hex literal::

    sage: preparse('0x2e3')
    'Integer(0x2e3)'
    sage: 0xA
    10
    sage: 0xe
    14

Raw and hex work correctly::

    sage: type(0xa1)
    <type 'sage.rings.integer.Integer'>
    sage: type(0xa1r)
    <type 'int'>
    sage: type(0Xa1R)
    <type 'int'>

The preparser can handle PEP 515 (see :trac:`28490`)::

    sage: 1_000_000 + 3_000 # py3
    1003000

In Sage, methods can also be called on integer and real literals (note
that in pure Python this would be a syntax error)::

    sage: 16.sqrt()
    4
    sage: 87.factor()
    3 * 29
    sage: 15.10.sqrt()
    3.88587184554509
    sage: preparse('87.sqrt()')
    'Integer(87).sqrt()'
    sage: preparse('15.10.sqrt()')
    "RealNumber('15.10').sqrt()"

Note that calling methods on int literals in pure Python is a syntax
error, but Sage allows this for Sage integers and reals, because users
frequently request it::

    sage: eval('4.__add__(3)')
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax

Symbolic functional notation::

    sage: a=10; f(theta, beta) = theta + beta; b = x^2 + theta
    sage: f
    (theta, beta) |--> beta + theta
    sage: a
    10
    sage: b
    x^2 + theta
    sage: f(theta,theta)
    2*theta

    sage: a = 5; f(x,y) = x*y*sqrt(a)
    sage: f
    (x, y) |--> sqrt(5)*x*y

This involves an =-, but should still be turned into a symbolic
expression::

    sage: preparse('a(x) =- 5')
    '__tmp__=var("x"); a = symbolic_expression(- Integer(5)).function(x)'
    sage: f(x)=-x
    sage: f(10)
    -10

This involves -=, which should not be turned into a symbolic
expression (of course a(x) isn't an identifier, so this will never be
valid)::

    sage: preparse('a(x) -= 5')
    'a(x) -= Integer(5)'

Raw literals:

Raw literals are not preparsed, which can be useful from an efficiency
point of view.  Just like Python ints are denoted by an L, in Sage raw
integer and floating literals are followed by an"r" (or "R") for raw,
meaning not preparsed.

We create a raw integer::

    sage: a = 393939r
    sage: a
    393939
    sage: type(a)
    <type 'int'>

We create a raw float::

    sage: z = 1.5949r
    sage: z
    1.5949
    sage: type(z)
    <type 'float'>

You can also use an upper case letter::

    sage: z = 3.1415R
    sage: z
    3.1415
    sage: type(z)
    <type 'float'>

This next example illustrates how raw literals can be very useful in
certain cases.  We make a list of even integers up to 10000::

    sage: v = [ 2*i for i in range(10000)]

This takes a noticeable fraction of a second (e.g., 0.25
seconds). After preparsing, what Python is really executing is the
following::

    sage: preparse('v = [ 2*i for i in range(10000)]')
    'v = [ Integer(2)*i for i in range(Integer(10000))]'

If instead we use a raw 2 we get execution that is *instant* (0.00
seconds)::

    sage: v = [ 2r * i for i in range(10000r)]

Behind the scenes what happens is the following::

    sage: preparse('v = [ 2r * i for i in range(10000r)]')
    'v = [ 2 * i for i in range(10000)]'

.. WARNING::

   The results of the above two expressions are different.  The
   first one computes a list of Sage integers, whereas the second
   creates a list of Python integers.  Python integers are typically
   much more efficient than Sage integers when they are very small;
   large Sage integers are much more efficient than Python integers,
   since they are implemented using the GMP C library.
"""

#*****************************************************************************
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function

import collections
import os
import re

from sage.repl.load import load_wrap

implicit_mul_level = False
numeric_literal_prefix = '_sage_const_'

def implicit_multiplication(level=None):
    r"""
    Turn implicit multiplication on or off, optionally setting a
    specific ``level``.

    INPUT:

    - ``level`` -- a boolean or integer (default: 5); how aggressive to be in
      placing \*'s

      -  0 - Do nothing
      -  1 - Numeric followed by alphanumeric
      -  2 - Closing parentheses followed by alphanumeric
      -  3 - Spaces between alphanumeric
      - 10 - Adjacent parentheses (may mangle call statements)

    OUTPUT:

    The current ``level`` if no argument is given.

    EXAMPLES::

      sage: implicit_multiplication(True)
      sage: implicit_multiplication()
      5
      sage: preparse('2x')
      'Integer(2)*x'
      sage: implicit_multiplication(False)
      sage: preparse('2x')
      '2x'

    Note that the `IPython automagic
    <https://ipython.readthedocs.io/en/stable/interactive/magics.html#magic-automagic>`_
    feature cannot be used if ``level >= 3``::

        sage: implicit_multiplication(3)
        sage: preparse('cd Documents')
        'cd*Documents'
        sage: implicit_multiplication(2)
        sage: preparse('cd Documents')
        'cd Documents'
        sage: implicit_multiplication(False)

    In this case, one can use the explicit syntax for IPython magics such as
    ``%cd Documents``.
    """
    global implicit_mul_level
    if level is None:
        return implicit_mul_level
    elif level is True:
        implicit_mul_level = 5
    else:
        implicit_mul_level = level

def isalphadigit_(s):
    """
    Return ``True`` if ``s`` is a non-empty string of alphabetic characters
    or a non-empty string of digits or just a single ``_``

    EXAMPLES::

        sage: from sage.repl.preparse import isalphadigit_
        sage: isalphadigit_('abc')
        True
        sage: isalphadigit_('123')
        True
        sage: isalphadigit_('_')
        True
        sage: isalphadigit_('a123')
        False
    """
    return s.isalpha() or s.isdigit() or s == "_"

in_single_quote = False
in_double_quote = False
in_triple_quote = False

def in_quote():
    return in_single_quote or in_double_quote or in_triple_quote

QuoteStackFrame = collections.namedtuple("QuoteStackFrame", [
    'delim',    # the quote characters used: ', ", ''', or """
    'raw',      # whether we're in a raw string
    'f_string', # whether we're in an F-string (PEP 498)
    'braces'    # in an F-string, how many unclosed {'s have we encountered?
])

def strip_string_literals(code, state=None):
    r"""
    Returns a string with all literal quotes replaced with labels and
    a dictionary of labels for re-substitution.  This makes parsing
    easier.

    INPUT:

    - ``code`` - a string; the input

    - ``state`` - a list of :class:`QuoteStackFrame` (default: None);
      state with which to continue processing, e.g., across multiple calls
      to this function

    OUTPUT:

    - a 4-tuple of the processed code, the dictionary of labels,
      any accumulated state, and the list of string delimiters
      (', ", ''', \"\"\") that are safe to insert into the processed code.

    EXAMPLES::

        sage: from sage.repl.preparse import strip_string_literals
        sage: s, literals, _, _ = strip_string_literals(r'''['a', "b", 'c', "d\""]''')
        sage: s
        '[%(L1)s, %(L2)s, %(L3)s, %(L4)s]'
        sage: literals
        {'L1': "'a'", 'L2': '"b"', 'L3': "'c'", 'L4': '"d\\""'}
        sage: print(s % literals)
        ['a', "b", 'c', "d\""]
        sage: print(strip_string_literals(r'-"\\\""-"\\"-')[0])
        -%(L1)s-%(L2)s-

    Triple-quotes are handled as well::

        sage: s, literals, _, _ = strip_string_literals("[a, '''b''', c, '']")
        sage: s
        '[a, %(L1)s, c, %(L2)s]'
        sage: print(s % literals)
        [a, '''b''', c, '']

    Comments are substitute too::

        sage: s, literals, _, _ = strip_string_literals("code '#' # ccc 't'"); s
        'code %(L1)s #%(L2)s'
        sage: s % literals
        "code '#' # ccc 't'"

    A state is returned so one can break strings across multiple calls to
    this function::

        sage: s, _, state, _ = strip_string_literals('s = "some'); s
        's = %(L1)s'
        sage: s, _, _, _ = strip_string_literals('thing" * 5', state); s
        '%(L1)s * 5'

    TESTS:

    Even for raw strings, a backslash can escape a following quote::

        sage: s, _, _, _ = strip_string_literals(r"r'somethin\' funny'"); s
        'r%(L1)s'
        sage: dep_regex = r'^ *(?:(?:cimport +([\w\. ,]+))|(?:from +(\w+) +cimport)|(?:include *[\'"]([^\'"]+)[\'"])|(?:cdef *extern *from *[\'"]([^\'"]+)[\'"]))' # Ticket 5821
    """
    new_code = []
    literals = {}
    counter = 0
    start = q = 0
    quote_stack = state or []
    safe_delims = collections.OrderedDict.fromkeys(["'", '"', "'''", '"""'])
    for quote in quote_stack:
        if quote.f_string:
            safe_delims.pop(quote.delim, None)
    while True:
        quote = quote_stack and quote_stack[-1]
        in_quote = quote and (not quote.f_string or not quote.braces)
        sig_q = code.find("'", q)
        dbl_q = code.find('"', q)
        hash_q = -1 if quote else code.find('#', q)
        if quote and quote.f_string:
            lbrace_q = code.find('{', q)
            rbrace_q = code.find('}', q)
            brace_q = min(lbrace_q, rbrace_q)
            if brace_q == -1:
                brace_q = max(lbrace_q, rbrace_q)
        else:
            brace_q = -1
        q = min(sig_q, dbl_q)
        if q == -1:
            q = max(sig_q, dbl_q)
        if hash_q != -1 and (q == -1 or hash_q < q):
            # it's a comment
            newline = code.find('\n', hash_q)
            if newline == -1:
                newline = len(code)
            counter += 1
            label = "L%s" % counter
            literals[label] = code[hash_q+1:newline]
            new_code.append(code[start:hash_q].replace('%','%%'))
            new_code.append("#%%(%s)s" % label)
            start = q = newline
        elif brace_q != -1 and (q == -1 or brace_q < q):
            # Inside the replacement section of an F-string?
            if quote.braces:
                # Treat the preceding substring as code.
                new_code.append(code[start:brace_q].replace('%','%%'))
            else:
                # Skip over {{ and }} escape sequences in the literal portion.
                if brace_q+1 < len(code) and code[brace_q+1] == code[brace_q]:
                    q = brace_q+2
                    continue
                # Treat the preceding substring as literal.
                counter += 1
                label = "L%s" % counter
                literals[label] = code[start:brace_q]
                new_code.append("%%(%s)s" % label)
            # Treat the brace itself as code.
            new_code.append(code[brace_q])
            # Increment/decrement brace count.
            if code[brace_q] == '{':
                quote_stack[-1] = quote._replace(braces=quote.braces+1)
            elif quote.braces > 0:
                quote_stack[-1] = quote._replace(braces=quote.braces-1)
            # Skip ahead just past the brace.
            start = q = brace_q+1
        elif q == -1:
            if in_quote:
                counter += 1
                label = "L%s" % counter
                literals[label] = code[start:]
                new_code.append("%%(%s)s" % label)
            else:
                new_code.append(code[start:].replace('%','%%'))
            break
        elif in_quote:
            if code[q-1] == '\\':
                k = 2
                while code[q-k] == '\\':
                    k += 1
                if k % 2 == 0:
                    q += 1
            if code[q:q+len(quote.delim)] == quote.delim:
                counter += 1
                label = "L%s" % counter
                literals[label] = code[start:q+len(quote.delim)]
                new_code.append("%%(%s)s" % label)
                q += len(quote.delim)
                start = q
                quote_stack.pop()
            else:
                q += 1
        else:
            if q>0 and code[q-1] in 'rR':
                raw = True
                f_string = q>1 and code[q-2] in 'fF'
            elif q>0 and code[q-1] in 'fF':
                f_string = True
                raw = q>1 and code[q-2] in 'rR'
            else:
                raw = f_string = False
            if len(code) >= q+3 and (code[q+1] == code[q] == code[q+2]):
                delim = code[q]*3
            else:
                delim = code[q]
            quote_stack.append(QuoteStackFrame(
                delim=delim, raw=raw, f_string=f_string, braces=0))
            if f_string:
                safe_delims.pop(delim, None)
            new_code.append(code[start:q].replace('%', '%%'))
            start = q
            q += len(delim)

    return "".join(new_code), literals, quote_stack, list(safe_delims)


def containing_block(code, idx, delimiters=['()','[]','{}'], require_delim=True):
    """
    Find the code block given by balanced delimiters that contains the position ``idx``.

    INPUT:

    - ``code`` - a string

    - ``idx`` - an integer; a starting position

    - ``delimiters`` - a list of strings (default: ['()', '[]',
      '{}']); the delimiters to balance. A delimiter must be a single
      character and no character can at the same time be opening and
      closing delimiter.

    - ``require_delim`` - a boolean (default: True); whether to raise
      a SyntaxError if delimiters are present. If the delimiters are
      unbalanced, an error will be raised in any case.

    OUTPUT:

    - a 2-tuple ``(a,b)`` of integers, such that ``code[a:b]`` is
      delimited by balanced delimiters, ``a<=idx<b``, and ``a``
      is maximal and ``b`` is minimal with that property. If that
      does not exist, a ``SyntaxError`` is raised.

    - If ``require_delim`` is false and ``a,b`` as above can not be
      found, then ``0, len(code)`` is returned.

    EXAMPLES::

        sage: from sage.repl.preparse import containing_block
        sage: s = "factor(next_prime(L[5]+1))"
        sage: s[22]
        '+'
        sage: start, end = containing_block(s, 22)
        sage: start, end
        (17, 25)
        sage: s[start:end]
        '(L[5]+1)'
        sage: s[20]
        '5'
        sage: start, end = containing_block(s, 20); s[start:end]
        '[5]'
        sage: start, end = containing_block(s, 20, delimiters=['()']); s[start:end]
        '(L[5]+1)'
        sage: start, end = containing_block(s, 10); s[start:end]
        '(next_prime(L[5]+1))'

    TESTS::

        sage: containing_block('((a{))',0)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('((a{))',1)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('((a{))',2)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('((a{))',3)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('((a{))',4)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('((a{))',5)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('(()()',1)
        (1, 3)
        sage: containing_block('(()()',3)
        (3, 5)
        sage: containing_block('(()()',4)
        (3, 5)
        sage: containing_block('(()()',0)
        Traceback (most recent call last):
        ...
        SyntaxError: Unbalanced delimiters
        sage: containing_block('(()()',0, require_delim=False)
        (0, 5)
        sage: containing_block('((})()',1, require_delim=False)
        (0, 6)
        sage: containing_block('abc',1, require_delim=False)
        (0, 3)

    """
    openings = "".join([d[0] for d in delimiters])
    closings = "".join([d[-1] for d in delimiters])
    levels = [0] * len(openings)
    p = 0
    start = idx
    while start >= 0:
        if code[start] in openings:
            p = openings.index(code[start])
            levels[p] -= 1
            if levels[p] == -1:
                break
        elif code[start] in closings and start < idx:
            p = closings.index(code[start])
            levels[p] += 1
        start -= 1
    if start == -1:
        if require_delim:
            raise SyntaxError("Unbalanced or missing delimiters")
        else:
            return 0, len(code)
    if levels.count(0) != len(levels)-1:
        if require_delim:
            raise SyntaxError("Unbalanced delimiters")
        else:
            return 0, len(code)
    p0 = p
    # We now have levels[p0]==-1. We go to the right hand side
    # till we find a closing delimiter of type p0 that makes
    # levels[p0]==0.
    end = idx
    while end < len(code):
        if code[end] in closings:
            p = closings.index(code[end])
            levels[p] += 1
            if p==p0 and levels[p] == 0:
                break
        elif code[end] in openings and end > idx:
            p = openings.index(code[end])
            levels[p] -= 1
        end += 1
    if levels.count(0) != len(levels):
        # This also occurs when end==len(code) without finding a closing delimiter
        if require_delim:
            raise SyntaxError("Unbalanced delimiters")
        else:
            return 0, len(code)
    return start, end+1


def parse_ellipsis(code, preparse_step=True):
    """
    Preparses [0,2,..,n] notation.

    INPUT:

    - ``code`` - a string

    - ``preparse_step`` - a boolean (default: True)

    OUTPUT:

    - a string

    EXAMPLES::

        sage: from sage.repl.preparse import parse_ellipsis
        sage: parse_ellipsis("[1,2,..,n]")
        '(ellipsis_range(1,2,Ellipsis,n))'
        sage: parse_ellipsis("for i in (f(x) .. L[10]):")
        'for i in (ellipsis_iter(f(x) ,Ellipsis, L[10])):'
        sage: [1.0..2.0]
        [1.00000000000000, 2.00000000000000]

    TESTS:

    Check that nested ellipsis is processed correctly (:trac:`17378`)::

        sage: preparse('[1,..,2,..,len([1..3])]')
        '(ellipsis_range(Integer(1),Ellipsis,Integer(2),Ellipsis,len((ellipsis_range(Integer(1),Ellipsis,Integer(3))))))'

    """
    ix = code.find('..')
    while ix != -1:
        if ix == 0:
            raise SyntaxError("Cannot start line with ellipsis.")
        elif code[ix-1]=='.':
            # '...' be valid Python in index slices
            code = code[:ix-1] + "Ellipsis" + code[ix+2:]
        elif len(code) >= ix+3 and code[ix+2]=='.':
            # '...' be valid Python in index slices
            code = code[:ix] + "Ellipsis" + code[ix+3:]
        else:
            start_list, end_list = containing_block(code, ix, ['()','[]'])

            #search the current containing block for other '..' occurrences that may
            #be contained in proper subblocks. Those need to be processed before
            #we can deal with the present level of ellipses.
            ix = code.find('..',ix+2,end_list)
            while ix != -1:
                if code[ix-1]!='.' and code[ix+2]!='.':
                    start_list,end_list = containing_block(code,ix,['()','[]'])
                ix = code.find('..',ix+2,end_list)

            arguments = code[start_list+1:end_list-1].replace('...', ',Ellipsis,').replace('..', ',Ellipsis,')
            arguments = re.sub(r',\s*,', ',', arguments)
            if preparse_step:
                arguments = arguments.replace(';', ', step=')
            range_or_iter = 'range' if code[start_list]=='[' else 'iter'
            code = "%s(ellipsis_%s(%s))%s" %  (code[:start_list],
                                               range_or_iter,
                                               arguments,
                                               code[end_list:])
        ix = code.find('..')
    return code

def extract_numeric_literals(code):
    """
    Pulls out numeric literals and assigns them to global variables.
    This eliminates the need to re-parse and create the literals,
    e.g., during every iteration of a loop.

    INPUT:

    - ``code`` - a string; a block of code

    OUTPUT:

    - a (string, string:string dictionary) 2-tuple; the block with
      literals replaced by variable names and a mapping from names to
      the new variables

    EXAMPLES::

        sage: from sage.repl.preparse import extract_numeric_literals
        sage: code, nums = extract_numeric_literals("1.2 + 5")
        sage: print(code)
        _sage_const_1p2  + _sage_const_5
        sage: print(nums)
        {'_sage_const_1p2': "RealNumber('1.2')", '_sage_const_5': 'Integer(5)'}

        sage: extract_numeric_literals("[1, 1.1, 1e1, -1e-1, 1.]")[0]
        '[_sage_const_1 , _sage_const_1p1 , _sage_const_1e1 , -_sage_const_1en1 , _sage_const_1p ]'

        sage: extract_numeric_literals("[1.sqrt(), 1.2.sqrt(), 1r, 1.2r, R.1, R0.1, (1..5)]")[0]
        '[_sage_const_1 .sqrt(), _sage_const_1p2 .sqrt(), 1 , 1.2 , R.1, R0.1, (_sage_const_1 .._sage_const_5 )]'
    """
    return preparse_numeric_literals(code, True)

all_num_regex = None

def preparse_numeric_literals(code, extract=False, quotes="'"):
    """
    This preparses numerical literals into their Sage counterparts,
    e.g. Integer, RealNumber, and ComplexNumber.

    INPUT:

    - ``code`` - a string; a code block to preparse

    - ``extract`` - a boolean (default: False); whether to create
      names for the literals and return a dictionary of
      name-construction pairs

    - ``quotes`` - a string (default: "'"); used to surround string
      arguments to Integer, RealNumber, and ComplexNumber

    OUTPUT:

    - a string or (string, string:string dictionary) 2-tuple; the
      preparsed block and, if ``extract`` is True, the
      name-construction mapping

    EXAMPLES::

        sage: from sage.repl.preparse import preparse_numeric_literals
        sage: preparse_numeric_literals("5")
        'Integer(5)'
        sage: preparse_numeric_literals("5j")
        "ComplexNumber(0, '5')"
        sage: preparse_numeric_literals("5jr")
        '5J'
        sage: preparse_numeric_literals("5l")
        '5l'
        sage: preparse_numeric_literals("5L")
        '5L'
        sage: preparse_numeric_literals("1.5")
        "RealNumber('1.5')"
        sage: preparse_numeric_literals("1.5j")
        "ComplexNumber(0, '1.5')"
        sage: preparse_numeric_literals(".5j")
        "ComplexNumber(0, '.5')"
        sage: preparse_numeric_literals("5e9j")
        "ComplexNumber(0, '5e9')"
        sage: preparse_numeric_literals("5.")
        "RealNumber('5.')"
        sage: preparse_numeric_literals("5.j")
        "ComplexNumber(0, '5.')"
        sage: preparse_numeric_literals("5.foo()")
        'Integer(5).foo()'
        sage: preparse_numeric_literals("5.5.foo()")
        "RealNumber('5.5').foo()"
        sage: preparse_numeric_literals("5.5j.foo()")
        "ComplexNumber(0, '5.5').foo()"
        sage: preparse_numeric_literals("5j.foo()")
        "ComplexNumber(0, '5').foo()"
        sage: preparse_numeric_literals("1.exp()")
        'Integer(1).exp()'
        sage: preparse_numeric_literals("1e+10")
        "RealNumber('1e+10')"
        sage: preparse_numeric_literals("0x0af")
        'Integer(0x0af)'
        sage: preparse_numeric_literals("0x10.sqrt()")
        'Integer(0x10).sqrt()'
        sage: preparse_numeric_literals('0o100')
        'Integer(0o100)'
        sage: preparse_numeric_literals('0b111001')
        'Integer(0b111001)'
        sage: preparse_numeric_literals('0xe')
        'Integer(0xe)'
        sage: preparse_numeric_literals('0xEAR')
        '0xEA'
        sage: preparse_numeric_literals('0x1012Fae')
        'Integer(0x1012Fae)'

    Test underscores as digit separators (PEP 515,
    https://www.python.org/dev/peps/pep-0515/)::

        sage: preparse_numeric_literals('123_456')
        'Integer(123_456)'
        sage: preparse_numeric_literals('123_456.78_9_0')
        "RealNumber('123_456.78_9_0')"
        sage: preparse_numeric_literals('0b11_011')
        'Integer(0b11_011)'
        sage: preparse_numeric_literals('0o76_321')
        'Integer(0o76_321)'
        sage: preparse_numeric_literals('0xaa_aaa')
        'Integer(0xaa_aaa)'
        sage: preparse_numeric_literals('1_3.2_5e-2_2')
        "RealNumber('1_3.2_5e-2_2')"

        sage: for f in ["1_1.", "11_2.", "1.1_1", "1_1.1_1", ".1_1", ".1_1e1_1", ".1e1_1",
        ....:           "1e12_3", "1_1e1_1", "1.1_3e1_2", "1_1e1_1", "1e1", "1.e1_1",
        ....:           "1.0", "1_1.0"]:
        ....:     preparse_numeric_literals(f)
        ....:     assert preparse(f) == preparse_numeric_literals(f), f
        "RealNumber('1_1.')"
        "RealNumber('11_2.')"
        "RealNumber('1.1_1')"
        "RealNumber('1_1.1_1')"
        "RealNumber('.1_1')"
        "RealNumber('.1_1e1_1')"
        "RealNumber('.1e1_1')"
        "RealNumber('1e12_3')"
        "RealNumber('1_1e1_1')"
        "RealNumber('1.1_3e1_2')"
        "RealNumber('1_1e1_1')"
        "RealNumber('1e1')"
        "RealNumber('1.e1_1')"
        "RealNumber('1.0')"
        "RealNumber('1_1.0')"

    Having consecutive underscores is not valid Python syntax, so
    it is not preparsed, and similarly with a trailing underscore::

        sage: preparse_numeric_literals('123__45')
        '123__45'
        sage: 123__45
        Traceback (most recent call last):
        ...
        SyntaxError: invalid ...

        sage: preparse_numeric_literals('3040_1_')
        '3040_1_'
        sage: 3040_1_
        Traceback (most recent call last):
        ...
        SyntaxError: invalid ...

    Using the ``quotes`` parameter::

        sage: preparse_numeric_literals('5j', quotes='"')
        'ComplexNumber(0, "5")'
        sage: preparse_numeric_literals('3.14', quotes="'''")
        "RealNumber('''3.14''')"
        sage: preparse_numeric_literals('01', quotes='~')
        'Integer(~01~)'

    """
    literals = {}
    last = 0
    new_code = []

    global all_num_regex
    if all_num_regex is None:
        hex_num = r"\b0x[0-9a-f]+(_[0-9a-f]+)*"
        oct_num = r"\b0o[0-7]+(_[0-7]+)*"
        bin_num = r"\b0b[01]+(_[01]+)*"
        # This is slightly annoying as floating point numbers may start
        # with a decimal point, but if they do the \b will not match.
        float_num = r"((\b\d+(_\d+)*([.](\d+(_\d+)*)?)?)|([.]\d+(_\d+)*))(e[-+]?\d+(_\d+)*)?"
        all_num = r"((%s)|(%s)|(%s)|(%s))(rj|rL|jr|Lr|j|L|r|)\b" % (hex_num, oct_num, bin_num, float_num)
        all_num_regex = re.compile(all_num, re.I)

    for m in all_num_regex.finditer(code):
        start, end = m.start(), m.end()
        num = m.group(1)
        postfix = m.groups()[-1].upper()

        if 'R' in postfix:
            postfix = postfix.replace('L', '')
            num_name = num_make = num + postfix.replace('R', '')
        elif 'L' in postfix:
            num_name = num_make = num + postfix.replace('L', '')
        else:

            # The Sage preparser does extra things with numbers, which we need to handle here.
            if '.' in num:
                if start > 0 and num[0] == '.':
                    if code[start-1] == '.':
                        # handle Ellipsis
                        start += 1
                        num = num[1:]
                    elif re.match(r'[\w\])]', code[start-1]):
                        # handle R.0
                        continue
                elif end < len(code) and num[-1] == '.':
                    if re.match(r'[^\W\d]', code[end]):
                        # handle 4.sqrt()
                        end -= 1
                        num = num[:-1]
            elif end < len(code) and code[end] == '.' and not postfix and re.match(r'\d+(_\d+)*$', num):
                # \b does not match after the . for floating point
                # two dots in a row would be an ellipsis
                if end+1 == len(code) or code[end+1] != '.':
                    end += 1
                    num += '.'

            num_name = numeric_literal_prefix + num.replace('.', 'p').replace('-', 'n').replace('+', '')

            if 'J' in postfix:
                num_make = "ComplexNumber(0, %s%s%s)" % (quotes, num, quotes)
                num_name += 'j'
            elif len(num) < 2 or num[1] in 'oObBxX':
                num_make = "Integer(%s)" % num
            elif '.' in num or 'e' in num or 'E' in num:
                num_make = "RealNumber(%s%s%s)" % (quotes, num, quotes)
            elif num[0] == "0":
                num_make = "Integer(%s%s%s)" % (quotes, num, quotes)
            else:
                num_make = "Integer(%s)" % num

            literals[num_name] = num_make

        new_code.append(code[last:start])
        if extract:
            new_code.append(num_name+' ')
        else:
            new_code.append(num_make)
        last = end

    new_code.append(code[last:])
    code = ''.join(new_code)
    if extract:
        return code, literals
    else:
        return code


def strip_prompts(line):
    r"""
    Removes leading sage: and >>> prompts so that pasting of examples
    from the documentation works.

    INPUT:

    - ``line`` - a string to process

    OUTPUT:

    - a string stripped of leading prompts

    EXAMPLES::

        sage: from sage.repl.preparse import strip_prompts
        sage: strip_prompts("sage: 2 + 2")
        '2 + 2'
        sage: strip_prompts(">>>   3 + 2")
        '3 + 2'
        sage: strip_prompts("  2 + 4")
        '  2 + 4'
    """
    for prompt in ['sage:', '>>>']:
        if line.startswith(prompt):
            line = line[len(prompt):].lstrip()
            break
    return line


def preparse_calculus(code):
    r"""
    Supports calculus-like function assignment, e.g., transforms::

       f(x,y,z) = sin(x^3 - 4*y) + y^x

    into::

       __tmp__=var("x,y,z")
       f = symbolic_expression(sin(x**3 - 4*y) + y**x).function(x,y,z)

    AUTHORS:

    - Bobby Moretti

      - Initial version - 02/2007

    - William Stein

      - Make variables become defined if they aren't already defined.

    - Robert Bradshaw

      - Rewrite using regular expressions (01/2009)

    EXAMPLES::

        sage: preparse("f(x) = x^3-x")
        '__tmp__=var("x"); f = symbolic_expression(x**Integer(3)-x).function(x)'
        sage: preparse("f(u,v) = u - v")
        '__tmp__=var("u,v"); f = symbolic_expression(u - v).function(u,v)'
        sage: preparse("f(x) =-5")
        '__tmp__=var("x"); f = symbolic_expression(-Integer(5)).function(x)'
        sage: preparse("f(x) -= 5")
        'f(x) -= Integer(5)'
        sage: preparse("f(x_1, x_2) = x_1^2 - x_2^2")
        '__tmp__=var("x_1,x_2"); f = symbolic_expression(x_1**Integer(2) - x_2**Integer(2)).function(x_1,x_2)'

    For simplicity, this function assumes all statements begin and end
    with a semicolon::

        sage: from sage.repl.preparse import preparse_calculus
        sage: preparse_calculus(";f(t,s)=t^2;")
        ';__tmp__=var("t,s"); f = symbolic_expression(t^2).function(t,s);'
        sage: preparse_calculus(";f( t , s ) = t^2;")
        ';__tmp__=var("t,s"); f = symbolic_expression(t^2).function(t,s);'

    TESTS:

    The arguments in the definition must be symbolic variables (:trac:`10747`)::

        sage: preparse_calculus(";f(_sage_const_)=x;")
        Traceback (most recent call last):
        ...
        ValueError: Argument names should be valid python identifiers.

    Although preparse_calculus returns something for f(1)=x, when
    preparsing a file an exception is raised because it is invalid python::

        sage: preparse_calculus(";f(1)=x;")
        ';__tmp__=var("1"); f = symbolic_expression(x).function(1);'

        sage: from sage.repl.preparse import preparse_file
        sage: preparse_file("f(1)=x")
        Traceback (most recent call last):
        ...
        ValueError: Argument names should be valid python identifiers.

        sage: from sage.repl.preparse import preparse_file
        sage: preparse_file("f(x,1)=2")
        Traceback (most recent call last):
        ...
        ValueError: Argument names should be valid python identifiers.

    Check support for unicode characters (:trac:`29278`)::

        sage: preparse("μ(x) = x^2")
        '__tmp__=var("x"); μ = symbolic_expression(x**Integer(2)).function(x)'
    """
    new_code = []
    last_end = 0
    #                                 f         (  vars  )   =      expr
    for m in re.finditer(r";(\s*)([^\W\d]\w*) *\(([^()]+)\) *= *([^;#=][^;#]*)", code):
        ident, func, vars, expr = m.groups()
        stripped_vars = [v.strip() for v in vars.split(',')]
        # if the variable name starts with numeric_literal_prefix
        # the argument name for the symbolic expression is a numeric literal
        # such as f(2)=5
        if any(n.startswith(numeric_literal_prefix) for n in stripped_vars):
            raise ValueError("Argument names should be valid python identifiers.")
        vars = ','.join(stripped_vars)

        new_code.append(code[last_end:m.start()])
        new_code.append(';%s__tmp__=var("%s"); %s = symbolic_expression(%s).function(%s)' %
                        (ident, vars, func, expr, vars))
        last_end = m.end()

    if last_end == 0:
        return code
    else:
        new_code.append(code[m.end():])
        return ''.join(new_code)


def preparse_generators(code):
    r"""
    Parses generator syntax, converting::

        obj.<gen0,gen1,...,genN> = objConstructor(...)

    into::

        obj = objConstructor(..., names=("gen0", "gen1", ..., "genN"))
        (gen0, gen1, ..., genN,) = obj.gens()

    and::

        obj.<gen0,gen1,...,genN> = R[interior]

    into::

        obj = R[interior]; (gen0, gen1, ..., genN,) = obj.gens()

    INPUT:

    - ``code`` - a string

    OUTPUT:

    - a string

    LIMITATIONS:

       - The entire constructor *must* be on one line.

    AUTHORS:

    - 2006-04-14: Joe Wetherell (jlwether@alum.mit.edu)

      - Initial version.

    - 2006-04-17: William Stein

      - Improvements to allow multiple statements.

    - 2006-05-01: William

      - Fix bug that Joe found.

    - 2006-10-31: William

      - Fix so obj doesn't have to be mutated.

    - 2009-01-27: Robert Bradshaw

      - Rewrite using regular expressions

    TESTS::

        sage: from sage.repl.preparse import preparse, preparse_generators

    Vanilla::

        sage: preparse("R.<x> = ZZ['x']")
        "R = ZZ['x']; (x,) = R._first_ngens(1)"
        sage: preparse("R.<x,y> = ZZ['x,y']")
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(2)"

    No square brackets::

        sage: preparse("R.<x> = PolynomialRing(ZZ, 'x')")
        "R = PolynomialRing(ZZ, 'x', names=('x',)); (x,) = R._first_ngens(1)"
        sage: preparse("R.<x,y> = PolynomialRing(ZZ, 'x,y')")
        "R = PolynomialRing(ZZ, 'x,y', names=('x', 'y',)); (x, y,) = R._first_ngens(2)"

    Names filled in::

        sage: preparse("R.<x> = ZZ[]")
        "R = ZZ['x']; (x,) = R._first_ngens(1)"
        sage: preparse("R.<x,y> = ZZ[]")
        "R = ZZ['x, y']; (x, y,) = R._first_ngens(2)"

    Names given not the same as generator names::

        sage: preparse("R.<x> = ZZ['y']")
        "R = ZZ['y']; (x,) = R._first_ngens(1)"
        sage: preparse("R.<x,y> = ZZ['u,v']")
        "R = ZZ['u,v']; (x, y,) = R._first_ngens(2)"

    Number fields::

        sage: preparse("K.<a> = QQ[2^(1/3)]")
        'K = QQ[Integer(2)**(Integer(1)/Integer(3))]; (a,) = K._first_ngens(1)'
        sage: preparse("K.<a, b> = QQ[2^(1/3), 2^(1/2)]")
        'K = QQ[Integer(2)**(Integer(1)/Integer(3)), Integer(2)**(Integer(1)/Integer(2))]; (a, b,) = K._first_ngens(2)'

    Just the .<> notation::

        sage: preparse("R.<x> = ZZx")
        'R = ZZx; (x,) = R._first_ngens(1)'
        sage: preparse("R.<x, y> = a+b")
        'R = a+b; (x, y,) = R._first_ngens(2)'
        sage: preparse("A.<x,y,z>=FreeAlgebra(ZZ,3)")
        "A = FreeAlgebra(ZZ,Integer(3), names=('x', 'y', 'z',)); (x, y, z,) = A._first_ngens(3)"

    Ensure we don't eat too much::

        sage: preparse("R.<x, y> = ZZ;2")
        'R = ZZ; (x, y,) = R._first_ngens(2);Integer(2)'
        sage: preparse("R.<x, y> = ZZ['x,y'];2")
        "R = ZZ['x,y']; (x, y,) = R._first_ngens(2);Integer(2)"
        sage: preparse("F.<b>, f, g = S.field_extension()")
        "F, f, g  = S.field_extension(names=('b',)); (b,) = F._first_ngens(1)"

    For simplicity, this function assumes all statements begin and end
    with a semicolon::

        sage: preparse_generators(";  R.<x>=ZZ[];")
        ";  R = ZZ['x']; (x,) = R._first_ngens(1);"

    See :trac:`16731`::

        sage: preparse_generators('R.<x> = ')
        'R.<x> = '

    Check support for unicode characters (:trac:`29278`)::

        sage: preparse('Ω.<λ,μ> = QQ[]')
        "Ω = QQ['λ, μ']; (λ, μ,) = Ω._first_ngens(2)"
    """
    new_code = []
    last_end = 0
    #                                obj       .< gens >      ,  other   =   constructor
    for m in re.finditer(r";(\s*)([^\W\d]\w*)\.<([^>]+)> *((?:,[\w, ]+)?)= *([^;#]+)", code):
        ident, obj, gens, other_objs, constructor = m.groups()
        gens = [v.strip() for v in gens.split(',')]
        constructor = constructor.rstrip()
        if len(constructor) == 0:
            pass   # SyntaxError will be raised by Python later
        elif constructor[-1] == ')':
            if '(' not in constructor:
                raise SyntaxError("Mismatched ')'")
            opening = constructor.rindex('(')
            # Only use comma if there are already arguments to the constructor
            comma = ', ' if constructor[opening+1:-1].strip() != '' else ''
            names = "('%s',)" % "', '".join(gens)
            constructor = constructor[:-1] + comma + "names=%s)" % names
        elif constructor[-1] == ']':
            # Could be nested.
            if '[' not in constructor:
                raise SyntaxError("Mismatched ']'")
            opening = constructor.rindex('[')
            closing = constructor.index(']', opening)
            if constructor[opening+1:closing].strip() == '':
                names = "'" + ', '.join(gens) + "'"
                constructor = constructor[:opening+1] + names + constructor[closing:]
        else:
            pass
        gens_tuple = "(%s,)" % ', '.join(gens)
        new_code.append(code[last_end:m.start()])
        new_code.append(";%s%s%s = %s; %s = %s._first_ngens(%s)" %
                        (ident, obj, other_objs, constructor, gens_tuple, obj, len(gens)))
        last_end = m.end()

    if last_end == 0:
        return code
    else:
        new_code.append(code[m.end():])
        return ''.join(new_code)


quote_state = None

def preparse(line, reset=True, do_time=False, ignore_prompts=False,
             numeric_literals=True):
    r"""
    Preparses a line of input.

    INPUT:

    - ``line`` - a string

    - ``reset`` - a boolean (default: True)

    - ``do_time`` - a boolean (default: False)

    - ``ignore_prompts`` - a boolean (default: False)

    - ``numeric_literals`` - a boolean (default: True)

    OUTPUT:

    - a string

    EXAMPLES::

        sage: preparse("ZZ.<x> = ZZ['x']")
        "ZZ = ZZ['x']; (x,) = ZZ._first_ngens(1)"
        sage: preparse("ZZ.<x> = ZZ['y']")
        "ZZ = ZZ['y']; (x,) = ZZ._first_ngens(1)"
        sage: preparse("ZZ.<x,y> = ZZ[]")
        "ZZ = ZZ['x, y']; (x, y,) = ZZ._first_ngens(2)"
        sage: preparse("ZZ.<x,y> = ZZ['u,v']")
        "ZZ = ZZ['u,v']; (x, y,) = ZZ._first_ngens(2)"
        sage: preparse("ZZ.<x> = QQ[2^(1/3)]")
        'ZZ = QQ[Integer(2)**(Integer(1)/Integer(3))]; (x,) = ZZ._first_ngens(1)'
        sage: QQ[2^(1/3)]
        Number Field in a with defining polynomial x^3 - 2 with a = 1.259921049894873?

        sage: preparse("a^b")
        'a**b'
        sage: preparse("a^^b")
        'a^b'
        sage: 8^1
        8
        sage: 8^^1
        9
        sage: 9^^1
        8

        sage: preparse("A \\ B")
        'A  * BackslashOperator() * B'
        sage: preparse("A^2 \\ B + C")
        'A**Integer(2)  * BackslashOperator() * B + C'
        sage: preparse("a \\ b \\") # There is really only one backslash here, it's just being escaped.
        'a  * BackslashOperator() * b \\'

        sage: preparse("time R.<x> = ZZ[]", do_time=True)
        '__time__=misc.cputime(); __wall__=misc.walltime(); R = ZZ[\'x\']; print("Time: CPU %.2f s, Wall: %.2f s"%(misc.cputime(__time__), misc.walltime(__wall__))); (x,) = R._first_ngens(1)'

    TESTS:

    Check support for unicode characters (:trac:`29278`)::

        sage: preparse("Ω.0")
        'Ω.gen(0)'
    """
    global quote_state
    if reset:
        quote_state = None

    L = line.lstrip()
    if len(L) > 0 and L[0] in ['#', '!']:
        return line

    if L.startswith('...'):
        i = line.find('...')
        return line[:i+3] + preparse(line[i+3:], reset=reset, do_time=do_time, ignore_prompts=ignore_prompts)

    if ignore_prompts:
        # Get rid of leading sage: and >>> so that pasting of examples from
        # the documentation works.
        line = strip_prompts(line)

    # This part handles lines with semi-colons all at once
    # Then can also handle multiple lines more efficiently, but
    # that optimization can be done later.
    L, literals, quote_state, safe_delims = strip_string_literals(line, quote_state)

    # Ellipsis Range
    # [1..n]
    try:
        L = parse_ellipsis(L, preparse_step=False)
    except SyntaxError:
        pass

    if implicit_mul_level:
        # Implicit Multiplication
        # 2x -> 2*x
        L = implicit_mul(L, level = implicit_mul_level)

    if numeric_literals and safe_delims:
        # Wrapping
        # 1 + 0.5 -> Integer(1) + RealNumber('0.5')
        L = preparse_numeric_literals(L, quotes=safe_delims[0])

    # Generators
    # R.0 -> R.gen(0)
    L = re.sub(r'(\b[^\W\d]\w*|[)\]])\.(\d+)', r'\1.gen(\2)', L)

    # Use ^ for exponentiation and ^^ for xor
    # (A side effect is that **** becomes xor as well.)
    L = L.replace('^', '**').replace('****', '^')

    # Make it easy to match statement ends
    L = ';%s;' % L.replace('\n', ';\n;')

    if do_time:
        # Separate time statement
        L = re.sub(r';(\s*)time +(\w)', r';time;\1\2', L)

    # Construction with generators
    # R.<...> = obj()
    # R.<...> = R[]
    L = preparse_generators(L)

    # Calculus functions
    # f(x,y) = x^3 - sin(y)
    L = preparse_calculus(L)

    # Backslash
    L = re.sub(r'''\\\s*([^\t ;#])''', r' * BackslashOperator() * \1', L)

    if do_time:
        # Time keyword
        L = re.sub(r';time;(\s*)(\S[^;]*)',
                   r';\1__time__=misc.cputime(); __wall__=misc.walltime(); \2; print(' +
                        '"Time: CPU %%.2f s, Wall: %%.2f s"%%(misc.cputime(__time__), misc.walltime(__wall__)))',
                   L)

    # Remove extra ;'s
    L = L.replace(';\n;', '\n')[1:-1]

    line = L % literals

    return line


######################################################
## Apply the preparser to an entire file
######################################################

def preparse_file(contents, globals=None, numeric_literals=True):
    """
    Preparses input, attending to numeric literals and load/attach
    file directives.

    .. note:: Temporarily, if @parallel is in the input, then
       numeric_literals is always set to False.

    INPUT:

    - ``contents`` - a string

    - ``globals`` - dict or None (default: None); if given, then
      arguments to load/attach are evaluated in the namespace of this
      dict.

    - ``numeric_literals`` - bool (default: True), whether to factor
      out wrapping of integers and floats, so they don't get created
      repeatedly inside loops

    OUTPUT:

    - a string

    TESTS::

        sage: from sage.repl.preparse import preparse_file
        sage: lots_of_numbers = "[%s]" % ", ".join(str(i) for i in range(3000))
        sage: _ = preparse_file(lots_of_numbers)
        sage: print(preparse_file("type(100r), type(100)"))
        _sage_const_100 = Integer(100)
        type(100 ), type(_sage_const_100 )
    """
    if not isinstance(contents, str):
        raise TypeError("contents must be a string")

    if globals is None:
        globals = {}

    # This is a hack, since when we use @parallel to parallelize code,
    # the numeric literals that are factored out do not get copied
    # to the subprocesses properly.  See trac #4545.
    if '@parallel' in contents:
        numeric_literals = False

    if numeric_literals:
        contents, literals, _, _ = strip_string_literals(contents)
        contents, nums = extract_numeric_literals(contents)
        contents = contents % literals
        if nums:
            # Stick the assignments at the top, trying not to shift
            # the lines down.
            ix = contents.find('\n')
            if ix == -1:
                ix = len(contents)
            if not re.match(r"^ *(#.*)?$", contents[:ix]):
                contents = "\n"+contents
            assignments = ["%s = %s" % x for x in nums.items()]
            # the preparser recurses on semicolons, so we only attempt
            # to preserve line numbers if there are a few
            if len(assignments) < 500:
                contents = "; ".join(assignments) + contents
            else:
                contents = "\n".join(assignments) + "\n\n" + contents

    # The list F contains the preparsed lines so far.
    F = []
    # A is the input, as a list of lines.
    A = contents.splitlines()
    # We are currently parsing the i-th input line.
    i = 0
    while i < len(A):
        L = A[i]
        do_preparse = True
        for cmd in ['load', 'attach']:
            if L.lstrip().startswith(cmd+' '):
                j = L.find(cmd+' ')
                s = L[j+len(cmd)+1:].strip()
                if not s.startswith('('):
                    F.append(' '*j + load_wrap(s, cmd=='attach'))
                    do_preparse = False
                    continue
        if do_preparse:
            F.append(preparse(L, reset=(i==0), do_time=True, ignore_prompts=False,
                              numeric_literals=not numeric_literals))
        i += 1

    return '\n'.join(F)


def implicit_mul(code, level=5):
    r"""
    Insert \*'s to make implicit multiplication explicit.

    INPUT:

    - ``code``  -- a string; the code with missing \*'s

    - ``level`` -- an integer (default: 5); see :func:`implicit_multiplication`
      for a list

    OUTPUT:

    - a string

    EXAMPLES::

        sage: from sage.repl.preparse import implicit_mul
        sage: implicit_mul('(2x^2-4x+3)a0')
        '(2*x^2-4*x+3)*a0'
        sage: implicit_mul('a b c in L')
        'a*b*c in L'
        sage: implicit_mul('1r + 1e3 + 5exp(2)')
        '1r + 1e3 + 5*exp(2)'
        sage: implicit_mul('f(a)(b)', level=10)
        'f(a)*(b)'

    TESTS:

    Check handling of Python 3 keywords (:trac:`29391`)::

        sage: implicit_mul('nonlocal a')
        'nonlocal a'

    Although these are not keywords in Python 3, we explicitly avoid implicit
    multiplication in these cases because the error message will be more
    helpful (:trac:`29391`)::

        sage: implicit_mul('print 2')
        'print 2'
        sage: implicit_mul('exec s')
        'exec s'

    Check support for unicode characters (:trac:`29278`)::

        sage: implicit_mul('3λ')
        '3*λ'
    """
    from keyword import iskeyword
    keywords_py2 = ['print', 'exec']

    def re_no_keyword(pattern, code):
        for _ in range(2): # do it twice in because matches don't overlap
            for m in reversed(list(re.finditer(pattern, code))):
                left, right = m.groups()
                if not iskeyword(left) and not iskeyword(right) \
                   and left not in keywords_py2:
                    code = "%s%s*%s%s" % (code[:m.start()],
                                          left,
                                          right,
                                          code[m.end():])
        return code

    code, literals, _, _ = strip_string_literals(code)
    if level >= 1:
        no_mul_token = " '''_no_mult_token_''' "
        code = re.sub(r'\b0x', r'0%sx' % no_mul_token, code)  # hex digits
        code = re.sub(r'( *)time ', r'\1time %s' % no_mul_token, code)  # first word may be magic 'time'
        code = re.sub(r'\b(\d+(?:\.\d+)?(?:e\d+)?)([rR]\b)', r'\1%s\2' % no_mul_token, code)  # exclude such things as 10r
        code = re.sub(r'\b(\d+(?:\.\d+)?)e([-\d])', r'\1%se%s\2' % (no_mul_token, no_mul_token), code)  # exclude such things as 1e5
        code = re_no_keyword(r'\b((?:\d+(?:\.\d+)?)|(?:%s[0-9eEpn]*\b)) *([^\W\d(]\w*)\b' % numeric_literal_prefix, code)
    if level >= 2:
        code = re.sub(r'(\%\(L\d+\))s', r'\1%ss%s' % (no_mul_token, no_mul_token), code) # literal strings
        code = re_no_keyword(r'(\)) *(\w+)', code)
    if level >= 3:
        code = re_no_keyword(r'(\w+) +(\w+)', code)
    if level >= 10:
        code = re.sub(r'\) *\(', ')*(', code)
    code = code.replace(no_mul_token, '')
    return code % literals



def _strip_quotes(s):
    """
    Strips one set of outer quotes.

    INPUT:

    - ``s`` - a string

    OUTPUT:

    - a string with any single and double quotes on either side of
      ``s`` removed

    EXAMPLES:

    Both types of quotes work::

        sage: import sage.repl.preparse
        sage: sage.repl.preparse._strip_quotes('"foo.sage"')
        'foo.sage'
        sage: sage.repl.preparse._strip_quotes("'foo.sage'")
        'foo.sage'

    The only thing that is stripped is at most one set of outer quotes::

        sage: sage.repl.preparse._strip_quotes('""foo".sage""')
        '"foo".sage"'
    """
    if len(s) == 0:
        return s
    if s[0] in ["'", '"']:
        s = s[1:]
    if s[-1] in ["'", '"']:
        s = s[:-1]
    return s


def handle_encoding_declaration(contents, out):
    r"""
    Find a PEP 263-style Python encoding declaration in the first or
    second line of ``contents``. If found, output it to ``out`` and return
    ``contents`` without the encoding line; otherwise output a default
    UTF-8 declaration and return ``contents``.

    EXAMPLES::

        sage: from sage.repl.preparse import handle_encoding_declaration
        sage: import sys
        sage: c1='# -*- coding: latin-1 -*-\nimport os, sys\n...'
        sage: c2='# -*- coding: iso-8859-15 -*-\nimport os, sys\n...'
        sage: c3='# -*- coding: ascii -*-\nimport os, sys\n...'
        sage: c4='import os, sys\n...'
        sage: handle_encoding_declaration(c1, sys.stdout)
        # -*- coding: latin-1 -*-
        'import os, sys\n...'
        sage: handle_encoding_declaration(c2, sys.stdout)
        # -*- coding: iso-8859-15 -*-
        'import os, sys\n...'
        sage: handle_encoding_declaration(c3, sys.stdout)
        # -*- coding: ascii -*-
        'import os, sys\n...'
        sage: handle_encoding_declaration(c4, sys.stdout)
        # -*- coding: utf-8 -*-
        'import os, sys\n...'

    TESTS:

    These are some of the tests listed in PEP 263::

        sage: contents = '#!/usr/bin/python\n# -*- coding: latin-1 -*-\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # -*- coding: latin-1 -*-
        '#!/usr/bin/python\nimport os, sys'

        sage: contents = '# This Python file uses the following encoding: utf-8\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # This Python file uses the following encoding: utf-8
        'import os, sys'

        sage: contents = '#!/usr/local/bin/python\n# coding: latin-1\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # coding: latin-1
        '#!/usr/local/bin/python\nimport os, sys'

    Two hash marks are okay; this shows up in SageTeX-generated scripts::

        sage: contents = '## -*- coding: utf-8 -*-\nimport os, sys\nprint(x)'
        sage: handle_encoding_declaration(contents, sys.stdout)
        ## -*- coding: utf-8 -*-
        'import os, sys\nprint(x)'

    When the encoding declaration doesn't match the specification, we
    spit out a default UTF-8 encoding.

    Incorrect coding line::

        sage: contents = '#!/usr/local/bin/python\n# latin-1\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # -*- coding: utf-8 -*-
        '#!/usr/local/bin/python\n# latin-1\nimport os, sys'

    Encoding declaration not on first or second line::

        sage: contents ='#!/usr/local/bin/python\n#\n# -*- coding: latin-1 -*-\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # -*- coding: utf-8 -*-
        '#!/usr/local/bin/python\n#\n# -*- coding: latin-1 -*-\nimport os, sys'

    We don't check for legal encoding names; that's Python's job::

        sage: contents ='#!/usr/local/bin/python\n# -*- coding: utf-42 -*-\nimport os, sys'
        sage: handle_encoding_declaration(contents, sys.stdout)
        # -*- coding: utf-42 -*-
        '#!/usr/local/bin/python\nimport os, sys'


    .. NOTE::

        - :pep:`263` says that Python will interpret a UTF-8
          byte order mark as a declaration of UTF-8 encoding, but I don't
          think we do that; this function only sees a Python string so it
          can't account for a BOM.

        - We default to UTF-8 encoding even though PEP 263 says that
          Python files should default to ASCII.

        - Also see https://docs.python.org/ref/encodings.html.

    AUTHORS:

    - Lars Fischer
    - Dan Drake (2010-12-08, rewrite for :trac:`10440`)
    """
    lines = contents.splitlines()
    for num, line in enumerate(lines[:2]):
        if re.search(r"coding[:=]\s*([-\w.]+)", line):
            out.write(line + '\n')
            return '\n'.join(lines[:num] + lines[(num+1):])

    # If we didn't find any encoding hints, use utf-8. This is not in
    # conformance with PEP 263, which says that Python files default to
    # ascii encoding.
    out.write("# -*- coding: utf-8 -*-\n")
    return contents

def preparse_file_named_to_stream(name, out):
    r"""
    Preparse file named \code{name} (presumably a .sage file), outputting to
    stream \code{out}.
    """
    name = os.path.abspath(name)
    with open(name) as f:
        contents = f.read()
    contents = handle_encoding_declaration(contents, out)
    parsed = preparse_file(contents)
    out.write('#'*70+'\n')
    out.write('# This file was *autogenerated* from the file %s.\n' % name)
    out.write('#'*70+'\n')
    out.write(parsed)

def preparse_file_named(name):
    r"""
    Preparse file named \code{name} (presumably a .sage file), outputting to a
    temporary file.  Returns name of temporary file.
    """
    from sage.misc.temporary_file import tmp_filename
    tmpfilename = tmp_filename(os.path.basename(name)) + '.py'
    with open(tmpfilename, 'w') as out:
        preparse_file_named_to_stream(name, out)
    return tmpfilename
