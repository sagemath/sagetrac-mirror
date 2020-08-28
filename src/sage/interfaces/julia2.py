r"""
Interface to Julia interpreter.

Julia is a high-level, high-performance dynamic programming language for
numerical computing, see the official website <JuliaLang https://julialang.org/>_
for installation and further information.

The commands in this section only work if you have "julia" installed and
available in your PATH. It's not necessary to install any special Sage packages.

EXAMPLES::

    sage: julia.eval('2+2')                 # optional - julia
    '4'

::

    sage: a = julia(10)                     # optional - julia
    sage: a**10                             # optional - julia
    10000000000

we use Nemo and Hecke automatically

::
    sage: V = julia('VectorSpace(GF(3),2)')           # optional - julia
    sage: V                                           # optional - julia
    Vector space of dimension 2 over Galois field with characteristic 3

AUTHORS:

- Marcelo Forets (2017-08-20) : initial implementation, strongly based on the
  matlab.py interface by W. Stein and D. Joyner.

Tutorial
--------

EXAMPLES::

    sage: julia('4+10')                     # optional - julia
    14
    sage: julia('5*10 + 6')                 # optional - julia
    56
    sage: julia('(6+6)/3')                  # optional - julia
    4.0
    sage: julia('9')^2                      # optional - julia
    81
    sage: a = julia(10); b = julia(20); c = julia(30)    # optional - julia
    sage: #avg = (a+b+c)/3; avg              # optional - julia
    sage: #parent(avg)                       # optional - julia

::

    sage: my_scalar = julia('3.1415')        # optional - julia
    sage: my_scalar                          # optional - julia
    3.1415
    sage: row_vector = julia('[1 5 7]')      # optional - julia
    sage: row_vector                         # optional - julia
    1×3 Array{Int64,2}:
     1  5  7

    sage: column_vector = julia('[1; 5; 7]')   # optional - julia
    sage: column_vector
    3-element Array{Int64,1}:
     1
     5
     7

    sage: row_vector * column_vector         # optional - julia
    1-element Array{Int64,1}:
    75

    sage: column_vector * row_vector         # optional - julia
    3×3 Array{Int64,2}:
     1   5   7
     5  25  35
     7  35  49

::

    sage: row_vector1 = julia('[1 2 3]')             # optional - julia
    sage: row_vector2 = julia('[3 2 1]')             # optional - julia
    sage: matrix_from_row_vec = julia('[%s; %s]'%(row_vector1.name(), row_vector2.name()))     # optional - julia
    sage: matrix_from_row_vec                            # optional - julia
    2×3 Array{Int64,2}:
     1  2  3
     3  2  1

::

    sage: column_vector1 = julia('[1;3]')               # optional - julia
    sage: column_vector2 = julia('[2;8]')               # optional - julia
    sage: matrix_from_col_vec = julia('[%s %s]'%(column_vector1.name(), column_vector2.name()))                                    # optional - julia
    sage: matrix_from_col_vec                            # optional - julia
    2×2 Array{Int64,2}:
     1  2
     3  8

::

    sage: my_matrix = julia('[8 12 19; 7 3 2; 12 4 23; 8 1 1]')    # optional - julia
    sage: my_matrix                                      # optional - julia
    4×3 Array{Int64,2}:
      8  12  19
      7   3   2
     12   4  23
      8   1   1

::

    sage: combined_matrix = julia('[%s %s]'%(my_matrix.name(), my_matrix.name()))                                        # optional - julia
    sage: combined_matrix                               # optional - julia
    4×6 Array{Int64,2}:
      8  12  19   8  12  19
      7   3   2   7   3   2
     12   4  23  12   4  23
      8   1   1   8   1   1

::

    sage: tm = julia('collect(0.5:2:10)')               # optional - julia
    sage: tm                                            # optional - julia
    5-element Array{Float64,1}:
     0.5
     2.5
     4.5
     6.5
     8.5

::

    sage: my_vector1 = julia('[1 5 7]')                # optional - julia
    sage: my_vector1[1]                                 # optional - julia
    1
    sage: my_vector1[2]                                 # optional - julia
    5
    sage: my_vector1[3]                                 # optional - julia
    7

Unicode characters work fine as long as they are valid python identifiers::

    sage: julia.eval('function å() return 1 end')
    'å (generic function with 1 method)'
    sage: julia.å()
    1

Unfortunately the exclaimation mark character (used to denote functions that modify their
arguments in Julia) is not a valid python identifier so we cannot do::
    
    sage: julia.pop!([1,2])
    Traceback (most recent call last):
    ...
    SyntaxError: invalid syntax

Instead we use BANG in place of ! to call such methods ::

    sage: julia.popBANG([1,2])
    2

Matrix indexing works as follows::

    sage: my_matrix = julia('[8 12 19; 7 3 2; 12 4 23; 8 1 1]')     # optional - julia
    sage: my_matrix[3,2]                                # optional - julia
    4

Setting using parenthesis cannot work (because of how the Python
language works). Use square brackets or the set function::

    sage: my_matrix = julia('[8 12 19; 7 3 2; 12 4 23; 8 1 1]')    # optional - julia
    sage: #my_matrix[2,3] = 1999    # optional - julia
    sage: #my_matrix         # optional - julia
    sage: my_matrix.setindexBANG(2000, int(2), int(3))                          # optional - julia
    4×3 Array{Int64,2}:
      8  12    19
      7   3  2000
     12   4    23
      8   1     1
"""

#*****************************************************************************
#       Copyright (C) 2017 Marcelo Forets <mforets@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

import re
import sys
import os

from sage.interfaces.expect import Expect, ExpectElement
from sage.docs.instancedoc import instancedoc

SAGE_REF = "_sage_ref"
SAGE_REF_RE = re.compile(r'%s\d+' % SAGE_REF)

class Julia(Expect):
    """
    Interface to the julia interpreter.

    EXAMPLES::

        sage: a = julia('[1 1 2; 3 5 8; 13 21 33]')             # optional - julia
        sage: b = julia('[1; 3; 13]')                           # optional - julia
        sage: c = a * b                                         # optional - julia
        sage: print(c)                                          # optional - julia
        3-element Array{Int64,1}:
          30
         122
         505
    """
    def __init__(self, maxread=None, script_subdirectory=None,
                 logfile=None, server=None,server_tmpdir=None):
        Expect.__init__(self,
                        name = 'julia',
                        prompt = 'julia> ',
                        command = "sage-native-execute julia -q",
                        server = server,
                        server_tmpdir = server_tmpdir,
                        script_subdirectory = script_subdirectory,
                        restart_on_ctrlc = False,
                        verbose_start = False,
                        logfile = logfile,
                        eval_using_file_cutoff=100)

    def __reduce__(self):
        return reduce_load_julia, tuple([])

    def _an_element_impl(self):
        """
        EXAMPLES::

            sage: julia._an_element_impl()
            0
        """
        return self(0)

    def _coerce_impl(self, x, use_special=True):
        if isinstance(x, int):
            return self(str(x))
        return super()._coerce_impl(self, x, use_special)


    def _read_in_file_command(self, filename):
        """
        Returns the command used to read in and execute a file in julia.

        EXAMPLES::

            sage: julia = Julia()
            sage: julia._read_in_file_command('/tmp/julia_file')
            'include("/tmp/julia_file");'

        Here is an indirect doctest to check that it does indeed work::

            sage: m = identity_matrix(ZZ, 10) # optional - julia
            sage: sm = julia.sage2julia_matrix_string(m) # optional - julia
            sage: m = julia(sm)                             # optional - julia
        """
        return 'include("{0}");'.format(filename)

    def _quit_string(self):
        return 'exit()'

    def _install_hints(self):
        return """
        You must obtain the program Julia in order to use Julia from Sage.
        You can read all about julia at <julialang.org http://julialang.org/>_.
        """

    def _start(self):
        Expect._start(self)
        self.eval("using Nemo, Hecke")

    def set(self, var, value):
        """
        Set the variable var to the given value.
        """
        cmd = '{0}={1};'.format(var, value)
        out = self.eval(cmd)
        if out.find("error") != -1:
            raise TypeError("Error executing code in julia\nCODE:\n\t{0}\njulia ERROR:\n\t{1}".format(cmd, out))

    def get(self, var):
        """
        Get the value of the variable var.

        EXAMPLES::

            sage: s = julia.eval('a = 2')      # optional - julia
            sage: julia.get('a')               # optional - julia
            '2'
        """
        s = self.eval('{0}'.format(var))
        return self.strip_answer(s)

    def strip_answer(self, s):
        r"""
        Returns the string s with Julia's answer prompt removed.

        EXAMPLES::

        """
        return s

    def eval(self, x, strip=True, **kwds):
        """
        Evaluate the given block x of code in Julia and return the output
        as a string.

        INPUT:

        -  ``x`` - string of code

        -  ``strip`` - ignored


        OUTPUT: string

        EXAMPLES:

        We evaluate a string that involves assigning to a
        variable and printing.

        ::

            sage: julia.eval("a = 10;print(2+a);")      # optional - julia
            '12'

        We evaluate a large input line (note that no weird output appears
        and that this works quickly).

        ::

            sage: julia.eval("a = %s;"%(10^10000))    # optional - julia
            ''

        A random example using Nemo::

            sage: nl=chr(10) # newline character
            sage: julia.eval(  # optional - julia
            ....: "R,x = PolynomialRing(QQ,\"x\")"+nl+
            ....: "g = R(0)"+nl+
            ....: "while g != 0 && roots(g) != []"+nl+
            ....: "  b=rand(-10:10)"+nl+
            ....: "  c=rand(-10:10)" + nl +
            ....: "  global g =3*b*x^4+18*c*x^3-6*b^2*x^2-6*b*c*x-b^3-9*c^2" + nl + 
            ....: "end "+nl+
            ....: "print(\"success\")")
            'success'

        Verify that :trac:`11401` is fixed::

            sage: nl=chr(10) # newline character
            sage: julia.eval("a=3;"+nl+"b=5;")  # optional - julia
            ''
            sage: julia("[a,b]")                  # optional - julia
            2-element Array{Int64,1}:
             3
             5

        """
        #x = self._preparse(x)
        x = str(x).rstrip()
        #if len(x) == 0 or x[len(x) - 1] != ';':
        #    x += ';'
        ans = Expect.eval(self, x, **kwds).rstrip()
        if 'ERROR:' in ans:
            raise RuntimeError("Error evaluating Julia code.\nIN:%s\nOUT:%s" % (x, ans))
        return ans

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: julia.function_call('sin', ['2'])
            0.9092974268256817
            sage: julia.sin(int(2))
            0.9092974268256817
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        function = function.replace("BANG","!")
        return self.new("%s(%s)" % (function,
                                    ",".join([s.name() for s in args])))

    def _true_symbol(self):
        """
        EXAMPLES::

            sage: julia._true_symbol()
            'true'
        """
        return 'true'

    def _false_symbol(self):
        """
        EXAMPLES::

            sage: julia._false_symbol()
            'false'
        """
        return 'false'

    def _equality_symbol(self):
        """
        """
        return "=="



    def console(self):
        julia_console()

    def version(self):
        return julia_version()[1:]

    def chdir(self, directory):
        """
        Change julia's current working directory.

        EXAMPLES::

            sage: julia.chdir('/')          # optional - julia
            sage: julia.pwd()               # optional - julia
            "/"

        """
        self.eval('cd("{0}")'.format(directory))

    def sage2julia_matrix_string(self, A):
        """
        Return a Julia matrix from a Sage matrix.

        INPUT: 

        - ``A`` -- a Sage matrix with entries in the rationals or reals

        OUTPUT:

        A string that evaluates to an julia matrix.

        EXAMPLES::

            sage: M33 = MatrixSpace(QQ,3,3)
            sage: A = M33([1,2,3,4,5,6,7,8,0])
            sage: julia.sage2julia_matrix_string(A)   # optional - julia
            '[1 2 3; 4 5 6; 7 8 0]'
        """
        return str(A.rows()).replace('), (', '; ').replace('(', '').replace(')','').replace(',', ' ')

    def _object_class(self):
        return JuliaElement


@instancedoc
class JuliaElement(ExpectElement):
    def __cmp__(self, other):
        """
        EXAMPLES::

            sage: one = julia(1); two = julia(2)
            sage: one == one
            True
            sage: one != two
            True
            sage: one < two
            True
            sage: two > one
            True
            sage: #one < 1
            sage: #two == 2
        """
        P = self._check_valid()
        if not hasattr(other, 'parent') or P is not other.parent():
            other = P(other)

        if P.eval('%s == %s' % (self.name(), other.name())) == P._true_symbol():
            return 0
        elif P.eval('%s < %s' % (self.name(), other.name())) == P._true_symbol():
            return -1
        else:
            return 1

    def bool(self):
        """
        EXAMPLES::

            sage: julia(int(1)).bool()
            True
            sage: julia(int(0)).bool()
            False
            sage: bool(julia(int(1)))
            True
        """
        P = self._check_valid()
        return P.eval("Bool(%s)" % self.name()) == P._true_symbol()

    def _add_(self, right):
        """
        EXAMPLES::

            sage: a = julia(1); b = julia(2)
            sage: a + b
            3
        """
        P = self._check_valid()
        return P.new('%s + %s' % (self.name(), right.name()))

    def _sub_(self, right):
        """
        EXAMPLES::

            sage: a = julia(1); b = julia(2)
            sage: a - b
            -1
        """
        P = self._check_valid()
        return P.new('%s - %s' % (self.name(), right.name()))

    def _mul_(self, right):
        """
        EXAMPLES::

            sage: a = julia(1); b = julia(2)
            sage: a * b
            2
        """
        P = self._check_valid()
        return P.new('%s * %s' % (self.name(), right.name()))

    def _div_(self, right):
        """
        EXAMPLES::

            sage: a = julia(1.0); b = julia(2.0)
            sage: a / b
            0.5
        """
        P = self._check_valid()
        return P.new('%s / %s' % (self.name(), right.name()))

    def __pow__(self, n):
        """
        EXAMPLES::

            sage: a = julia(3)
            sage: a^3
            27
        """
        P = self._check_valid()
        right = P(n)
        return P.new('%s ^ %s' % (self.name(), right.name()))

    def __len__(self):
        r"""
        Return length of this Julia element.

        This is the same as ``length(self)`` in Julia.

        EXAMPLES::

            sage: len(julia("1:9"))                                      # optional - julia
            9
        """
        P = self._check_valid()
        return int(P.eval('length(%s)' % self.name()))

    def __iter__(self):
        """
        EXAMPLES::

            sage: A = julia("[1,2,3]")
            sage: list(iter(A))
            [1, 2, 3]
            sage: A = julia("(1,2,3)")
            sage: list(iter(A))
            [1, 2, 3]
        """
        l = len(self)

        for i in range(1, int(l + 1)):
            yield self[i]


# An instance
julia = Julia(logfile="logo")

def reduce_load_Julia():
    return julia


def julia_console():
    """
    This requires that the optional julia program be installed and in
    your PATH, but no optional Sage packages need be installed.

    EXAMPLES::

        sage: julia_console()                # optional - julia; not tested
                                       < M A T L A B >
                           Copyright 1984-2006 The MathWorks, Inc.
        ...
        >> 2+3

    ans =

    5

    exit()

    Typing quit exits the Julia console and returns you to Sage.
    Julia, like Sage, remembers its history from one session to
    another.
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%julia magics instead.')
    os.system('julia')


def julia_version():
    """
    Return the version of Julia installed.

    EXAMPLES::

        sage: julia_version()    # random; optional - julia
        v"0.6.0"
    """
    return str(julia('VERSION')).strip()
