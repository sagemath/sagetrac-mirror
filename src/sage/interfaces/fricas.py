r"""
Interface to FriCAS

.. TODO::

    - some conversions in ``sage.functions`` are still missing and
      all should be checked and tested

FriCAS is a free GPL-compatible (modified BSD license) general
purpose computer algebra system based on Axiom.  The FriCAS
website can be found at http://fricas.sourceforge.net/.

AUTHORS:

- Mike Hansen (2009-02): Split off the FriCAS interface from
  the Axiom interface.

- Martin Rubey, Bill Page (2016-08): Completely separate from Axiom,
  implement more complete translation from FriCAS to SageMath types.


EXAMPLES::

    sage: fricas('3 * 5')                                                       # optional - fricas
    15
    sage: a = fricas(3) * fricas(5); a                                          # optional - fricas
    15

The type of a is :class:`FriCASElement`, i.e., an element of the
FriCAS interpreter::

    sage: type(a)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: a.parent()                                                            # optional - fricas
    FriCAS

The underlying FriCAS type of a is also available, via the type
method::

    sage: a.typeOf()                                                            # optional - fricas
    PositiveInteger

FriCAS objects are normally displayed using "ASCII art"::

    sage: fricas(2/3)                                                           # optional - fricas
      2
      -
      3
    sage: fricas('x^2 + 3/7')                                                   # optional - fricas
       2   3
      x  + -
           7

Functions defined in FriCAS are available as methods of the :class:`fricas<FriCAS>` object::

    sage: F = fricas.factor('x^5 - y^5'); F                                     # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )
    sage: type(F)                                                               # optional - fricas
    <class 'sage.interfaces.fricas.FriCASElement'>
    sage: F.typeOf()                                                            # optional - fricas
    Factored(Polynomial(Integer))

We can also create a FriCAS polynomial and apply the function
``factor`` from FriCAS.  The notation ``f.factor()`` is consistent
with how the rest of SageMath works::

    sage: f = fricas('x^5 - y^5')                                               # optional - fricas
    sage: f^2                                                                   # optional - fricas
     10      5 5    10
    y   - 2 x y  + x
    sage: f.factor()                                                            # optional - fricas
               4      3    2 2    3     4
    - (y - x)(y  + x y  + x y  + x y + x )

For many FriCAS types, translation to an appropriate SageMath type is
available::

    sage: f.factor().sage()                                                     # optional - fricas
    (y - x) * (y^4 + y^3*x + y^2*x^2 + y*x^3 + x^4)

Control-C interruption works well with the FriCAS interface. For
example, try the following sum but with a much bigger range, and hit
control-C::

    sage:  f = fricas('(x^5 - y^5)^10000')                                      # not tested - fricas
    Interrupting FriCAS...
    ...
    KeyboardInterrupt: Ctrl-c pressed while running FriCAS

Let us demonstrate some features of FriCAS.  FriCAS can guess a
differential equation for the generating function for integer
partitions::

    sage: fricas("guessADE([partition n for n in 0..40], homogeneous==4)")      # optional - fricas
    [
      [
          n
        [x ]f(x):
               2    3 (iv)           2    2 ,              3  ,,,
              x f(x) f    (x) + (20 x f(x) f (x) + 5 x f(x) )f   (x)
    <BLANKLINE>
            +
                    2    2 ,,   2
              - 39 x f(x) f  (x)
    <BLANKLINE>
            +
                   2     ,   2            2 ,            3  ,,         2 ,   4
              (12 x f(x)f (x)  - 15 x f(x) f (x) + 4 f(x) )f  (x) + 6 x f (x)
    <BLANKLINE>
            +
                        ,   3          2 ,   2
              10 x f(x)f (x)  - 16 f(x) f (x)
    <BLANKLINE>
          =
            0
        ,
                         2      3      4
       f(x) = 1 + x + 2 x  + 3 x  + O(x )]
      ]

FriCAS can solve linear ordinary differential equations::

    sage: fricas.set("y", "operator y")                                         # optional - fricas
    sage: fricas.set("deq", "x^3*D(y x, x, 3) + x^2*D(y x, x, 2) - 2*x*D(y x, x) + 2*y x - 2*x^4")  # optional - fricas
    sage: fricas.set("sol", "solve(deq, y, x)"); fricas("sol")                  # optional - fricas
                   5       3       2
                  x  - 10 x  + 20 x  + 4
    [particular = ----------------------,
                           15 x
                 3      2       3       3      2
              2 x  - 3 x  + 1  x  - 1  x  - 3 x  - 1
     basis = [---------------, ------, -------------]]
                     x            x          x

    sage: fricas("sol.particular").sage()                                       # optional - fricas
    1/15*(x^5 - 10*x^3 + 20*x^2 + 4)/x
    sage: fricas("sol.basis").sage()                                            # optional - fricas
    [(2*x^3 - 3*x^2 + 1)/x, (x^3 - 1)/x, (x^3 - 3*x^2 - 1)/x]
    sage: fricas.eval(")clear values y deq sol")                                # optional - fricas
    ''

FriCAS can expand expressions into series::

    sage: x = var('x'); ex = sqrt(cos(x)); a = fricas(ex).series(x=0); a        # optional - fricas
        1  2    1  4    19   6     559   8     29161    10      11
    1 - - x  - -- x  - ---- x  - ------ x  - --------- x   + O(x  )
        4      96      5760      645120      116121600

    sage: a.coefficients()[38].sage()                                           # optional - fricas
    -29472026335337227150423659490832640468979/274214482066329363682430667508979749984665600000000

    sage: ex = sqrt(atan(x)); a = fricas(ex).series(x=0); a                     # optional - fricas
     1      5        9
     -      -        -
     2   1  2    31  2      6
    x  - - x  + --- x  + O(x )
         6      360

    sage: a.coefficient(9/2).sage()                                             # optional - fricas
    31/360

    sage: x = fricas("x::TaylorSeries Fraction Integer")                        # optional - fricas
    sage: y = fricas("y::TaylorSeries Fraction Integer")                        # optional - fricas
    sage: 2*(1+2*x+sqrt(1-4*x)-2*x*y).recip()                                   # optional - fricas
                  2       3     2 2      3       4        4        5
      1 + (x y + x ) + 2 x  + (x y  + 2 x y + 6 x ) + (4 x y + 18 x )
    +
        3 3      4 2       5        6        5 2       6         7
      (x y  + 3 x y  + 13 x y + 57 x ) + (6 x y  + 40 x y + 186 x )
    +
        4 4      5 3       6 2        7         8
      (x y  + 4 x y  + 21 x y  + 130 x y + 622 x )
    +
          6 3       7 2        8          9
      (8 x y  + 66 x y  + 432 x y + 2120 x )
    +
        5 5      6 4       7 3        8 2         9          10
      (x y  + 5 x y  + 30 x y  + 220 x y  + 1466 x y + 7338 x  ) + O(11)

FriCAS does some limits right::

    sage: x = var('x'); ex = x^2*exp(-x)*Ei(x) - x; fricas(ex).limit(x=oo)      # optional - fricas
    1

"""

###########################################################################
#       Copyright (C) 2008 Mike Hansen <mhansen@gmail.com>
#                     2007 Bill Page
#                     2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
###########################################################################

import re
import os
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.interfaces.expect import Expect, ExpectElement, FunctionElement, ExpectFunction
from sage.misc.misc import SAGE_TMP_INTERFACE
from sage.env import DOT_SAGE, LOCAL_IDENTIFIER
from sage.docs.instancedoc import instancedoc
from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.misc.lazy_import import lazy_import
lazy_import('sage.libs.pynac.pynac', ['symbol_table'])
lazy_import('sage.calculus.var', ['var', 'function'])


FRICAS_SINGLE_LINE_START = 3  # where output starts when it fits next to the line number
FRICAS_MULTI_LINE_START = 2   # and when it doesn't
FRICAS_LINE_LENGTH = 80       # length of a line, should match the line length in sage
# the following messages have, unfortunately, no markup.
FRICAS_WHAT_OPERATIONS_STRING = r"Operations whose names satisfy the above pattern\(s\):"
FRICAS_ERROR_IN_LIBRARY_CODE = ">> Error detected within library code:"

# only the last command should be necessary to make the interface
# work, the other are optimizations.  Beware that lisp distinguishes
# between ' and ".
FRICAS_INIT_CODE = (
    ")set functions compile on",
    ")set message autoload off",
    ")set message type off",
    ")set output length " + str(FRICAS_LINE_LENGTH),
    ")lisp (setf |$ioHook|"
    "            (lambda (x &optional args)"
    "              (when (member x '(|startAlgebraOutput| |endOfAlgebraOutput|"
    "                                |startKeyedMsg|      |endOfKeyedMsg|))"
    "               (prin1 x)"
    "               (princ #\\Newline))))")
# code (one-liners!) executed after having set up the prompt
FRICAS_HELPER_CODE = (
    'sageprint(x:InputForm):String == ' +
    '(atom? x => (' +
    'float? x => return float(x)::String;' +
    'integer? x => return integer(x)::String;' +
    'string? x => return concat(["_"", string(x)::String, "_""])$String;' +
    'symbol? x => return string(symbol(x)));' +
    'S: List String := [sageprint y for y in destruct x];' +
    'R: String := new(1 + reduce(_+, [1 + #(s)$String for s in S], 0),' +
    'space()$Character);' +
    'copyInto!(R, "(", 1);' +
    'i := 2;' +
    'for s in S repeat'
    '(copyInto!(R, s, i); i := i + 1 + #(s)$String);' +
    'copyInto!(R, ")", i-1);' +
    'return R)',)

FRICAS_LINENUMBER_OFF_CODE = ")lisp (setf |$IOindex| NIL)"
FRICAS_FIRST_PROMPT = r"\(1\) -> "
FRICAS_LINENUMBER_OFF_PROMPT = r"\(NIL\) -> "


class FriCAS(ExtraTabCompletion, Expect):
    """
    Interface to a FriCAS interpreter.
    """
    def __init__(self, name='fricas', command='fricas -nosman',
                 script_subdirectory=None, logfile=None,
                 server=None, server_tmpdir=None):
        """
        Create an instance of the FriCAS interpreter.

        TESTS::

            sage: fricas == loads(dumps(fricas))                                # optional - fricas
            True

        Check that :trac:`25174` is fixed::

            sage: fricas(I)                                                     # optional - fricas
            %i

            sage: integrate(sin(x)*exp(I*x), x, -pi, 0, algorithm="fricas")     # optional - fricas
            1/2*I*pi

            sage: fricas(I*sin(x)).sage()                                       # optional - fricas
            I*sin(x)

            sage: fricas(I*x).sage()                                            # optional - fricas
            I*x
        """
        eval_using_file_cutoff = 4096 - 5  # magic number from Expect._eval_line (there might be a bug)
        assert max(len(c) for c in FRICAS_INIT_CODE) < eval_using_file_cutoff
        self.__eval_using_file_cutoff = eval_using_file_cutoff
        self._COMMANDS_CACHE = '%s/%s_commandlist_cache.sobj' % (DOT_SAGE, name)
        # we run the init code in _start to avoid spurious output
        Expect.__init__(self,
                        name=name,
                        prompt=FRICAS_FIRST_PROMPT,
                        command=command,
                        script_subdirectory=script_subdirectory,
                        server=server,
                        server_tmpdir=server_tmpdir,
                        restart_on_ctrlc=False,
                        verbose_start=False,
                        init_code=[],
                        logfile=logfile,
                        eval_using_file_cutoff=eval_using_file_cutoff)

    def _start(self):
        """
        Start the FriCAS interpreter and switch off the linenumbers.

        EXAMPLES::

            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
        """
        # setting the prompt properly is necessary for restarting FriCAS
        self._prompt = FRICAS_FIRST_PROMPT
        Expect._start(self)
        for line in FRICAS_INIT_CODE:
            self.eval(line, reformat=False)
        # switching off the line numbers also modified the prompt
        self._prompt = FRICAS_LINENUMBER_OFF_PROMPT
        self.eval(FRICAS_LINENUMBER_OFF_CODE, reformat=False)
        for line in FRICAS_HELPER_CODE:
            self.eval(line, reformat=False)

    def _install_hints(self):
        """
        Hints for installing FriCAS on your computer.

        EXAMPLES::

            sage: print(fricas._install_hints())
            In order...
        """
        return r"""
In order to use the FriCAS interface you need to have FriCAS installed.
You can either run 'sage -i fricas' to install FriCAS as an optional
package within SageMath, or install FriCAS separately, see
http://fricas.sourceforge.net.
"""

    def _quit_string(self):
        """
        Return the string used to quit FriCAS.

        EXAMPLES::

            sage: fricas._quit_string()                                         # optional - fricas
            ')quit'
            sage: a = FriCAS()                                                  # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False
            sage: a._start()                                                    # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            True
            sage: a.quit()                                                      # optional - fricas
            sage: a.is_running()                                                # optional - fricas
            False

        TESTS::

            sage: import psutil                                                 # optional - fricas
            sage: p = fricas.pid(); pr = psutil.Process(p); pr                  # optional - fricas
            <psutil.Process(pid=..., name='FRICASsys') at ...>
            sage: pr.children()                                                 # optional - fricas
            []
        """
        return ')quit'

    def _commands(self):
        """
        Return a list of commands available. This is done by parsing the
        result of the first section of the output of ')what things'.

        EXAMPLES::

            sage: cmds = fricas._commands()                                     # optional - fricas
            sage: len(cmds) > 100                                               # optional - fricas
            True
            sage: '<' in cmds                                                   # optional - fricas
            True
            sage: 'factor' in cmds                                              # optional - fricas
            True
        """
        output = self.eval(")what operations", reformat=False)
        m = re.search(FRICAS_WHAT_OPERATIONS_STRING + r"\n(.*)\n\|startKeyedMsg\|",
                      output, flags=re.DOTALL)
        l = m.groups()[0].split()
        return l

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Return a list of all the commands defined in Fricas and optionally
        (per default) store them to disk.

        EXAMPLES::

            sage: c = fricas._tab_completion(use_disk_cache=False, verbose=False)         # optional - fricas
            sage: len(c) > 100                                                  # optional - fricas
            True
            sage: 'factor' in c                                                 # optional - fricas
            True
            sage: '**' in c                                                     # optional - fricas
            False
            sage: 'upperCase?' in c                                             # optional - fricas
            False
            sage: 'upperCase_q' in c                                            # optional - fricas
            True
            sage: 'upperCase_e' in c                                            # optional - fricas
            True
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(self._COMMANDS_CACHE)
                    return self.__tab_completion
                except IOError:
                    pass
            if verbose:
                print("\nBuilding %s command completion list (this takes" % self)
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % self._COMMANDS_CACHE)
            v = self._commands()

            # process the commands to strip out things which are not
            # valid Python identifiers
            valid = re.compile('[^a-zA-Z0-9_]+')
            names = [x for x in v if valid.search(x) is None]

            # replace trailing ? with _q and trailing ! with _e
            names += [x[:-1] + "_q" for x in v if x.endswith("?")]
            names += [x[:-1] + "_e" for x in v if x.endswith("!")]

            self.__tab_completion = names
            if len(v) > 200:
                # Fricas is actually installed.
                sage.misc.persist.save(v, self._COMMANDS_CACHE)
            return names

    def _read_in_file_command(self, filename):
        """
        Return the FriCAS command to read the file ``filename``.

        INPUT:

        - ``filename``, a string ending in '.input'.

        OUTPUT:

        - a string with the command for reading filename without output.

        TESTS:

        Evaluate a rather long line::

            sage: len(fricas([i for i in range(600)]))                          # optional - fricas, indirect doctest
            600

        """
        if not filename.endswith('.input'):
            raise ValueError("the filename must end with .input")

        return ')read %s )quiet' % filename

    def _local_tmpfile(self):
        """
        Return a local tmpfile ending with ".input" used to buffer long
        command lines sent to FriCAS.

        """
        try:
            return self.__local_tmpfile
        except AttributeError:
            self.__local_tmpfile = os.path.join(SAGE_TMP_INTERFACE, 'tmp' + str(self.pid()) + '.input')
            return self.__local_tmpfile

    def _remote_tmpfile(self):
        """
        Return a remote tmpfile ending with ".input" used to buffer long
        command lines sent to FriCAS.

        """
        try:
            return self.__remote_tmpfile
        except AttributeError:
            self.__remote_tmpfile = self._remote_tmpdir() + "/interface_%s:%s.input" % (LOCAL_IDENTIFIER, self.pid())
            return self.__remote_tmpfile

# what I expect from FriCAS:

# 1.) in set(self, var, value)
#
# no markers:
# there could be some "debugging" output, as in fricas("guessADE([1,1,1,1], debug==true)")
#
# startKeyedMsg: an error happened
#
# 2.) in get(self, var)
#
# |startAlgebraOutput\|...|endOfAlgebraOutput\|
#
# 3.) I also need a routine to send a system command and get its output.

    def _check_errors(self, line, output):
        """
        Check whether output contains an error and, if so, raise it.

        INPUT:

        - ``line``, a string that was sent to FriCAS.

        - ``output``, a string returned by FriCAS

        OUTPUT:

        None

        TESTS::

            sage: fricas.set("x", "[i fo83r i in 0..17]")                       # optional - fricas, indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: An error occurred when FriCAS evaluated '[i fo83r i in 0..17]':
              Line   1: x:=[i fo83r i in 0..17];
                       ...A..........B
              Error  A: Missing mate.
              Error  B: syntax error at top level
              Error  B: Possibly missing a ]
               3 error(s) parsing

            sage: fricas.set("x", "something stupid")                           # optional - fricas, indirect doctest
            Traceback (most recent call last):
            ...
            RuntimeError: An error occurred when FriCAS evaluated 'something stupid':
               There are no library operations named something
                  Use HyperDoc Browse or issue
                                           )what op something
                  to learn if there is any operation containing " something " in its
                  name.
            <BLANKLINE>
               Cannot find a definition or applicable library operation named
                  something with argument type(s)
                                            Variable(stupid)
            <BLANKLINE>
                  Perhaps you should use "@" to indicate the required return type, or
                  "$" to specify which version of the function you need.

        """
        # otherwise there might be a message
        m = re.search(r"\|startKeyedMsg\|\n(.*)\n\|endOfKeyedMsg\|",
                      output, flags=re.DOTALL)
        if m:
            replacements = [('|startKeyedMsg|\n', ''),
                            ('|endOfKeyedMsg|', '')]
            for old, new in replacements:
                output = output.replace(old, new)
            raise RuntimeError("An error occurred when FriCAS evaluated '%s':\n%s" % (line, output))

        # or even an error
        if FRICAS_ERROR_IN_LIBRARY_CODE in output:
            raise RuntimeError("An error occurred when FriCAS evaluated '%s':\n%s" % (line, output))

    def set(self, var, value):
        """
        Set a variable to a value in FriCAS.

        INPUT:

        - ``var``, ``value``: strings, the first representing a valid
          FriCAS variable identifier, the second a FriCAS expression.

        OUTPUT: None

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'

        """
        cmd = '%s%s%s;' % (var, self._assign_symbol(), value)
        output = self.eval(cmd, reformat=False)
        self._check_errors(value, output)

    def get(self, var):
        r"""
        Get the string representation of the value (more precisely, the
        OutputForm) of a variable or expression in FriCAS.

        If FriCAS cannot evaluate `var` an error is raised.

        EXAMPLES::

            sage: fricas.set('xx', '2')                                         # optional - fricas
            sage: fricas.get('xx')                                              # optional - fricas
            '2'
            sage: a = fricas('(1 + sqrt(2))^5')                                 # optional - fricas
            sage: fricas.get(a.name())                                          # optional - fricas
            '    +-+\n29 \\|2  + 41'
            sage: fricas.get('(1 + sqrt(2))^5')                                 # optional - fricas
            '    +-+\n29 \\|2  + 41'
            sage: fricas.new('(1 + sqrt(2))^5')                                 # optional - fricas
                +-+
            29 \|2  + 41
        """
        output = self.eval(str(var), reformat=False)
        # if there is AlgebraOutput we ask no more
        m = re.search(r"\|startAlgebraOutput\|\n(.*)\n\|endOfAlgebraOutput\|",
                      output, flags=re.DOTALL)
        if m:
            lines = m.groups()[0].split("\n")
            if max(len(line) for line in lines) < FRICAS_LINE_LENGTH:
                return "\n".join(line[FRICAS_SINGLE_LINE_START:] for line in lines)
            else:
                return "\n".join(line[FRICAS_MULTI_LINE_START:] for line in lines)

        self._check_errors(var, output)

    def get_string(self, var):
        """
        Return the value of a FriCAS string as a string, without checking
        that it is a string.

        TESTS:

        We test that strings are returned properly::

            sage: r = fricas.get_string('concat([concat(string(i)," ") for i in 0..299])')   # optional - fricas
            sage: r == " ".join([str(i) for i in range(300)]) + ' '                          # optional - fricas
            True

            sage: fricas.get_string('concat([string(1) for i in 1..5])') == "1"*5            # optional - fricas
            True

            sage: fricas.get_string('concat([string(1) for i in 1..10000])') == "1"*10000    # optional - fricas
            True

        A problem with leading space::

            sage: s = "unparse((-1234567890123456789012345678901234567890123456789012345678901234567890*n::EXPR INT)::INFORM)"
            sage: fricas.get_string(s)                                                       # optional - fricas
            '(-1234567890123456789012345678901234567890123456789012345678901234567890)*n'

        Check that :trac:`25628` is fixed::

            sage: var("a b"); f = 1/(1+a*cos(x))                                # optional - fricas
            (a, b)
            sage: lF = integrate(f, x, algorithm="fricas")                      # optional - fricas
            sage: (diff(lF[0], x) - f).simplify_trig()                          # optional - fricas
            0
            sage: (diff(lF[1], x) - f).simplify_trig()                          # optional - fricas
            0
            sage: f = 1/(b*x^2+a); lF = integrate(f, x, algorithm="fricas"); lF # optional - fricas
            [1/2*log((2*a*b*x + (b*x^2 - a)*sqrt(-a*b))/(b*x^2 + a))/sqrt(-a*b),
             arctan(sqrt(a*b)*x/a)/sqrt(a*b)]
            sage: (diff(lF[0], x) - f).simplify_trig()                          # optional - fricas
            0
            sage: (diff(lF[1], x) - f).simplify_trig()                          # optional - fricas
            0

        """
        # strip removes leading and trailing whitespace, after that
        # we can assume that the first and the last character are
        # double quotes
        return self.get(str(var)).replace("\n", "").strip()[1:-1]

    def get_integer(self, var):
        """
        Return the value of a FriCAS integer as an integer, without
        checking that it is an integer.

        TESTS::

            sage: fricas.get_integer('factorial 1111') == factorial(1111)       # optional - fricas
            True

        """
        return int(self.get_unparsed_InputForm(str(var)))

    def get_boolean(self, var):
        """
        Return the value of a FriCAS boolean as a boolean, without checking
        that it is a boolean.

        TESTS::

            sage: fricas.get_boolean('(1=1)::Boolean') == True                  # optional - fricas
            True

            sage: fricas.get_boolean('(1=2)::Boolean') == False                 # optional - fricas
            True
        """
        return self.get(str(var)).replace("\n", "") == "true"

    def get_unparsed_InputForm(self, var):
        """
        Return the unparsed ``InputForm`` as a string.

        .. TODO::

            - catch errors, especially when InputForm is not available:

                - for example when integration returns ``"failed"``

                - ``UnivariatePolynomial``

            - should we provide workarounds, too?

        TESTS::

            sage: fricas.get_unparsed_InputForm('1..3')                         # optional - fricas
            '(1..3)$Segment(PositiveInteger())'

        """
        return self.get_string('unparse((%s)::InputForm)' % str(var))

    def get_InputForm(self, var):
        """
        Return the ``InputForm`` as a string.

        TESTS::

            sage: fricas.get_InputForm('1..3')                                  # optional - fricas
            '(($elt (Segment (PositiveInteger)) SEGMENT) 1 3)'

        """
        return self.get_string('sageprint((%s)::InputForm)' % str(var))

    def _assign_symbol(self):
        """
        Return the symbol used for setting a variable in FriCAS.

        EXAMPLES::

            sage: fricas.set("x", "1");                                         # optional - fricas, indirect doctest
            sage: fricas.get("x")                                               # optional - fricas
            '1'
            sage: fricas.eval(")cl val x")                                      # optional - fricas
            ''
        """
        return ":="

    def _equality_symbol(self):
        """
        Return the equality testing logical symbol in FriCAS.

        EXAMPLES::

            sage: a = fricas(x==6); a                                           # optional - fricas, indirect doctest
            x = 6

        A warning::

            sage: fricas.set("x", 2);                                           # optional - fricas
            sage: a = fricas(x==6); a                                           # optional - fricas
            2 = 6
            sage: fricas.eval(")cl val x")                                      # optional - fricas
            ''
        """
        return "="

    def _true_symbol(self):
        """
        Return the string used for ``True`` in FriCAS.

        EXAMPLES::

            sage: str(fricas("(1=1)@Boolean")) == fricas._true_symbol()         # optional - fricas
            True
        """
        return "true"

    def _false_symbol(self):
        """
        Return the string used for ``False`` in FriCAS.

        EXAMPLES::

            sage: str(fricas("(1~=1)@Boolean")) == fricas._false_symbol()       # optional - fricas
            True
        """
        return "false"

    def _inequality_symbol(self):
        """
        Return the string used for False in FriCAS.

        EXAMPLES::

            sage: fricas(x!=0)                                                  # optional - fricas, indirect doctest
            true
        """
        return '~='

    def _repr_(self):
        """
        EXAMPLES::

            sage: fricas                                                        # indirect doctest
            FriCAS
        """
        return "FriCAS"

    def __reduce__(self):
        """
        EXAMPLES::

            sage: fricas.__reduce__()                                           # optional - fricas
            (<function reduce_load_fricas at 0x...>, ())
            sage: f, args = _                                                   # optional - fricas
            sage: f(*args)                                                      # optional - fricas
            FriCAS
        """
        return reduce_load_fricas, tuple([])

    def eval(self, code, strip=True, synchronize=False, locals=None, allow_use_file=True,
             split_lines="nofile", reformat=True, **kwds):
        """
        Evaluate ``code`` using FriCAS.

        Except ``reformat``, all arguments are passed to
        :meth:`sage.interfaces.expect.Expect.eval`.

        INPUT:

        - ``reformat`` -- bool; remove the output markers when True.

        This can also be used to pass system commands to FriCAS.

        EXAMPLES::

            sage: fricas.set("x", "1783"); fricas("x")                          # optional - fricas
            1783
            sage: fricas.eval(")cl val x");                                     # optional - fricas
            ''
            sage: fricas("x")                                                   # optional - fricas
            x

        """
        output = Expect.eval(self, code, strip=strip,
                             synchronize=synchronize, locals=locals,
                             allow_use_file=allow_use_file, split_lines=split_lines,
                             **kwds)
        # we remove carriage returns (\r) to make parsing easier
        # they are sent depending on how fricas was invoked:
        # on linux, "fricas -nox -noclef" sends "\r\n" and "fricas -nosman" sends "\n"
        output = output.replace('\r', '')
        if reformat:
            replacements = [('|startAlgebraOutput|\n', ''),
                            ('|endOfAlgebraOutput|', ''),
                            ('|startKeyedMsg|\n', ''),
                            ('|endOfKeyedMsg|', '')]
            for old, new in replacements:
                output = output.replace(old, new)

        return output

    def _function_class(self):
        """
        Return the FriCASExpectFunction class.

        EXAMPLES::

            sage: fricas._function_class()                                      # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
            sage: type(fricas.gcd)                                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASExpectFunction'>
        """
        return FriCASExpectFunction

    def _object_class(self):
        """
        EXAMPLES::

            sage: fricas._object_class()                                        # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
            sage: type(fricas(2))                                               # optional - fricas
            <class 'sage.interfaces.fricas.FriCASElement'>
        """
        return FriCASElement

    def _function_element_class(self):
        """
        Return the FriCAS function element class.

        EXAMPLES::

            sage: fricas._function_element_class()                              # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
            sage: type(fricas(2).gcd)                                           # optional - fricas
            <class 'sage.interfaces.fricas.FriCASFunctionElement'>
        """
        return FriCASFunctionElement

    def console(self):
        """
        Spawn a new FriCAS command-line session.

        EXAMPLES::

            sage: fricas.console()                                              # not tested
                             FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
            -----------------------------------------------------------------------------
               Issue )copyright to view copyright notices.
               Issue )summary for a summary of useful system commands.
               Issue )quit to leave AXIOM and return to shell.
            -----------------------------------------------------------------------------
        """
        fricas_console()


@instancedoc
class FriCASElement(ExpectElement):
    """
    Instances of this class represent objects in FriCAS.

    Using the method :meth:`sage` we can translate some of them to
    SageMath objects:

    .. automethod:: _sage_
    """
    def __len__(self):
        """
        Return the length of a list.

        EXAMPLES::

            sage: v = fricas('[x^i for i in 0..5]')    # optional - fricas
            sage: len(v)                               # optional - fricas
            6

        TESTS:

        Streams are not handled yet::

            sage: oh = fricas('[i for i in 1..]')    # optional - fricas
            sage: len(oh)                            # optional - fricas
            Traceback (most recent call last):
            ...
            TypeError: ...
        """
        P = self._check_valid()
        l = P('#(%s)' % self._name)
        return l.sage()

    def __iter__(self):
        """
        Return an iterator over ``self``.

        EXAMPLES::

            sage: L = fricas([4,5,6])   # optional - fricas
            sage: list(L)               # optional - fricas
            [4, 5, 6]

        TESTS:

        Streams are not handled yet::

            sage: oh = fricas('[i for i in 1..]')    # optional - fricas
            sage: next(iter(oh))       # known bug
        """
        for i in range(len(self)):
            yield self[i]

    def __getitem__(self, n):
        """
        We implement the sage conventions here, translating to 0-based iterables.

        We do not check validity, since many objects in FriCAS are
        iterable, in particular Streams

        TESTS::

            sage: fricas("[1,2,3]")[0]                                          # optional - fricas
            1

        Negative indices do work::

            sage: fricas("[1,2,3]")[-1]    # optional - fricas
            3

            sage: fricas("[1,2,3]")[-2]    # optional - fricas
            2

        Invalid indices raise exceptions::

            sage: fricas("[1,2,3]")[3]     # optional - fricas
            Traceback (most recent call last):
            ...
            TypeError: An error occurred when FriCAS evaluated 'elt(...,...)':
            <BLANKLINE>
            >> Error detected within library code:
            index out of range

        And streams are ok too::

            sage: oh = fricas('[i for i in 1..]')    # optional - fricas
            sage: oh[4]                              # optional - fricas
            5
        """
        n = int(n)
        P = self._check_valid()
        if n < -1:
            try:
                l = len(self)
            except TypeError:
                raise
            else:
                n += l
                if not(0 <= n < l):
                    raise IndexError("index out of range")
        # use "elt" instead of "." here because then the error
        # message is clearer
        if n == -1:
            return P.new("elt(%s,last)" % (self._name))
        return P.new("elt(%s,%s)" % (self._name, n + 1))

    def __int__(self):
        """
        TESTS::

            sage: int(fricas(2))                                                # optional - fricas
            2
        """
        return int(self.sage())

    def bool(self):
        """
        Coerce the expression into a boolean.

        EXAMPLES::

            sage: fricas("1=1").bool()                                          # optional - fricas
            True
            sage: fricas("1~=1").bool()                                         # optional - fricas
            False
        """
        P = self._check_valid()
        return P.new(self._name + "::Boolean").sage()

    def __bool__(self):
        """
        Check whether the expression is different from zero.

        EXAMPLES::

            sage: fricas(0).is_zero()                                           # optional - fricas, indirect doctest
            True
        """
        P = self._check_valid()
        return not P.new("zero?(%s)" % self._name).sage()

    __nonzero__ = __bool__

    def __float__(self):
        """
        TESTS::

            sage: float(fricas(2))                                              # optional - fricas
            2.0
        """
        return float(self.sage())

    def _integer_(self, ZZ=None):
        """
        EXAMPLES::

            sage: ZZ(fricas('1'))                                               # optional - fricas
            1
        """
        return ZZ(self.sage())

    def _rational_(self):
        """
        EXAMPLES::

            sage: QQ(fricas('-1/2'))                                            # optional - fricas
            -1/2
        """
        return QQ(self.sage())

    def gen(self, n):
        """
        Return an error, since the n-th generator in FriCAS is not well defined.
        """
        raise NotImplementedError

    def _latex_(self):
        r"""
        EXAMPLES::

            sage: latex(fricas("sin(x+y)/exp(z)*log(1+%e)"))                    # optional - fricas
            {{\log \left( {{e+1}} \right)} \  {\sin \left( {{y+x}} \right)}} \over {{e} ^{z}}

            sage: latex(fricas("matrix([[1,2],[3,4]])"))                        # optional - fricas
            \left[ \begin{array}{cc} 1 & 2 \\ 3 & 4\end{array} \right]

            sage: latex(fricas("integrate(sin(x+1/x),x)"))                      # optional - fricas
            \int ^{\displaystyle x} {{\sin \left( {{{{{ \%O} ^{2}}+1} \over  \%O}} \right)} \  {d \%O}}
        """
        replacements = [(r'\sp ', '^'),
                        (r'\sp{', '^{'),
                        (r'\sb ', '_'),
                        (r'\sb{', '_{')]
        P = self._check_valid()
        s = P.get_string("latex(%s)" % self._name)
        for old, new in replacements:
            s = s.replace(old, new)
        return s

    def _get_sage_type(self, domain):
        """
        INPUT:

        - ``domain``, a FriCAS SExpression

        OUTPUT:

        - a corresponding Sage type

        EXAMPLES::

            sage: m = fricas("dom(1/2)::Any")                                   # optional - fricas
            sage: fricas(0)._get_sage_type(m)                                   # optional - fricas
            Rational Field

        TESTS::

            sage: m = fricas("UP(y, UP(x, AN))::INFORM")                        # optional - fricas
            sage: fricas(0)._get_sage_type(m)                                   # optional - fricas
            Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Algebraic Field
        """
        from sage.rings.all import QQbar, RDF, PolynomialRing
        from sage.rings.fraction_field import FractionField
        from sage.rings.finite_rings.integer_mod_ring import Integers
        from sage.rings.finite_rings.finite_field_constructor import FiniteField
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR

        # first implement domains without arguments
        head = str(domain.car())
        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ
        if head == "String":
            return str
        if head == "Float":
            P = self._check_valid()
            prec = max(P.new("length mantissa(%s)" % self._name).sage(), 53)
            return RealField(prec)
        if head == "DoubleFloat":
            return RDF
        if head == "AlgebraicNumber":
            return QQbar

        # now implement "functorial" types
        if head == "OrderedCompletion" or head == "Complex":
            # this is a workaround, I don't know how translate this
            return SR

        if head == "IntegerMod":
            return Integers(domain[1].integer().sage())

        if head == "PrimeField":
            return FiniteField(domain[1].integer().sage())

        if head == "Fraction":
            return FractionField(self._get_sage_type(domain[1]))

        if head == "Expression":
            return SR

        if head == "Polynomial":
            # this is a workaround, since in sage we always have to specify the variables
            return SR

        if head == "UnivariatePolynomial":
            var = str(domain[1])
            return PolynomialRing(self._get_sage_type(domain[2]), var)

        raise NotImplementedError("the translation of FriCAS type %s to sage is not yet implemented" % domain)

    _WHITESPACE = " "
    _LEFTBRACKET = "("
    _RIGHTBRACKET = ")"
    _STRINGMARKER = '"'
    _ESCAPEMARKER = '_'  # STRINGMARKER must be escaped in strings

    @staticmethod
    def _parse_and_eval(s, start=0):
        """
        Parse and evaluate the string.

        INPUT:

        - ``s`` -- string
        - ``start`` -- integer; specifies where to start parsing

        OUTPUT:

        - a pair ``(L, end)``, where ``L`` is the parsed list and
          ``end`` is the position of the last parsed letter.

        TESTS::

            sage: from sage.interfaces.fricas import FriCASElement
            sage: FriCASElement._parse_and_eval("abc")
            (abc, 2)

            sage: FriCASElement._parse_and_eval("(asin c)")
            (arcsin(c), 7)

            sage: N(FriCASElement._parse_and_eval("(pi)")[0])                   # optional - fricas
            3.14159265358979

            sage: FriCASElement._parse_and_eval('(a "(b c)")')
            Traceback (most recent call last):
            ...
            TypeError: cannot coerce arguments: no canonical coercion from <type 'str'> to Symbolic Ring

        """
        a = start
        while s[a] in FriCASElement._WHITESPACE:
            a += 1

        if s[a] == FriCASElement._LEFTBRACKET:
            return FriCASElement._parse_list(s, start=a)
        elif s[a] == FriCASElement._STRINGMARKER:
            return FriCASElement._parse_string(s, start=a)
        else:
            return FriCASElement._parse_other(s, start=a)

    @staticmethod
    def _parse_list(s, start=0):
        """
        Parse the initial part of a string, assuming that it is a
        whitespace separated list, treating its first element as
        function and the rest as arguments.

        INPUT:

        - ``s`` -- string
        - ``start`` -- integer; specifies the position of the left
          bracket

        TESTS::

            sage: from sage.interfaces.fricas import FriCASElement
            sage: FriCASElement._parse_list("()")
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object is not callable

            sage: FriCASElement._parse_list("(a b c)")
            (a(b, c), 6)

            sage: FriCASElement._parse_list('(bcd)')
            (bcd(), 4)

        """
        a = start
        assert s[a] == FriCASElement._LEFTBRACKET
        a += 1
        # the first element of a list must be a function call
        fun, a = FriCASElement._parse_other(s, start=a, make_fun=True)
        a += 1
        args = []
        while s[a] != FriCASElement._RIGHTBRACKET:
            e, a = FriCASElement._parse_and_eval(s, start=a)
            args.append(e)
            a += 1
        return fun(*args), a

    @staticmethod
    def _parse_other(s, start=0, make_fun=False):
        """
        Parse the initial part of a string, assuming that it is an
        atom, but not a string.

        Symbols and numbers must not contain ``FriCASElement._WHITESPACE`` and
        ``FriCASElement._RIGHTBRACKET``.

        INPUT:

        - ``s`` -- string
        - ``start`` -- integer; specifies where the symbol begins
        - ``make_fun`` -- (default: ``False``) a Boolean; specifying
          whether the atom should be interpreted as a function call

        TESTS::

            sage: from sage.interfaces.fricas import FriCASElement
            sage: FriCASElement._parse_other("abc")
            (abc, 2)
            sage: FriCASElement._parse_other("123 xyz")
            (123, 2)
            sage: FriCASElement._parse_other("abc -1.23", 4)
            (-1.23, 8)

        This function uses the symbol table to translate symbols
        which are not function calls.  At least ``%pi`` is an
        example showing that this may be necessary::

            sage: FriCASElement._parse_other("%pi")
            (pi, 2)

        """
        a = start
        b = len(s)
        while a < b and s[a] not in FriCASElement._WHITESPACE and s[a] != FriCASElement._RIGHTBRACKET:
            a += 1

        e = s[start:a]
        if make_fun:
            try:
                e = symbol_table["fricas"][e]
            except KeyError:
                e = function(e)
        else:
            try:
                e = ZZ(e)
            except TypeError:
                try:
                    e = float(e)
                except ValueError:
                    try:
                        e = symbol_table["fricas"][e]
                    except KeyError:
                        e = var(e.replace("%", "_"))
        return e, a - 1

    @staticmethod
    def _parse_string(s, start=0):
        r"""
        Parse the initial part of a string, assuming that it represents a
        string.

        INPUT:

        - ``s`` -- string
        - ``start`` -- integer; specifies the position of the left
          quote

        TESTS::

            sage: from sage.interfaces.fricas import FriCASElement
            sage: FriCASElement._parse_string('"abc" 123')
            ('abc', 4)

            sage: FriCASElement._parse_string('"" 123')
            ('', 1)

            sage: FriCASElement._parse_string('"____" 123')
            ('__', 5)

            sage: FriCASElement._parse_string('"_a" 123')
            ('a', 3)

            sage: FriCASElement._parse_string('"_" _"" 123')
            ('" "', 6)

            sage: FriCASElement._parse_string('"(b c)"')
            ('(b c)', 6)

        """
        a = start
        assert s[a] == FriCASElement._STRINGMARKER
        b = a + 1
        a = b
        S = []
        while s[b] != FriCASElement._STRINGMARKER:
            if s[b] == FriCASElement._ESCAPEMARKER:
                S.append(s[a:b])
                b += 1
                a = b
            b += 1
        S.append(s[a:b])
        return "".join(S), b

    @staticmethod
    def _sage_expression(fricas_InputForm):
        r"""
        Convert an expression to an element of the Symbolic Ring.

        INPUT:

        - fricas_InputForm, a string, the InputForm of a FriCAS expression.

        TESTS::

            sage: f = fricas('integrate(sin(x^2), x)'); f                       # optional - fricas
                       +---+
                       | 2
            fresnelS(x |--- )
                      \|%pi
            -----------------
                   +---+
                   | 2
                   |---
                  \|%pi
            sage: s = fricas.get_InputForm(f._name); s                          # optional - fricas
            '(/ (fresnelS (* x (^ (/ 2 (pi)) (/ 1 2)))) (^ (/ 2 (pi)) (/ 1 2)))'
            sage: from sage.interfaces.fricas import FriCASElement
            sage: FriCASElement._sage_expression(s)                             # optional - fricas
            1/2*sqrt(2)*sqrt(pi)*fresnel_sin(sqrt(2)*x/sqrt(pi))

        Check that :trac:`22525` is fixed::

            sage: l = [sin, cos, sec, csc, cot, tan, asin, acos, atan, acot, acsc, asec, arcsin, arccos, arctan, arccot, arccsc, arcsec]
            sage: [f(x)._fricas_().sage().subs(x=0.9) for f in l]               # optional - fricas
            [0.783326909627483,
             0.621609968270664,
             1.60872581046605,
             1.27660621345890,
             0.793551147842317,
             1.26015821755034,
             1.11976951499863,
             0.451026811796262,
             0.732815101786507,
             0.837981225008390,
             1.57079632679490 - 0.467145308103262*I,
             0.467145308103262*I,
             1.11976951499863,
             0.451026811796262,
             0.732815101786507,
             0.837981225008390,
             1.57079632679490 - 0.467145308103262*I,
             0.467145308103262*I]
            sage: l = [tanh, sinh, cosh, coth, sech, csch, asinh, acosh, atanh, acoth, asech, acsch, arcsinh, arccosh, arctanh, arccoth, arcsech, arccsch]
            sage: [f(x)._fricas_().sage().subs(x=0.9) for f in l]               # optional - fricas
            [0.716297870199024,
             1.02651672570818,
             1.43308638544877,
             1.39606725303001,
             0.697794641100332,
             0.974168247780004,
             0.808866935652782,
             0.451026811796262*I,
             1.47221948958322,
             1.47221948958322 - 1.57079632679490*I,
             0.467145308103262,
             0.957800449200672,
             0.808866935652782,
             0.451026811796262*I,
             1.47221948958322,
             1.47221948958322 - 1.57079632679490*I,
             0.467145308103262,
             0.957800449200672]

        Check that :trac:`23782` is fixed::

            sage: s = '((3*n^10-25*n^9+50*n^8+62*n^7-229*n^6-25*n^5+320*n^4-12*n^3-144*n^2)/11520)::EXPR INT'
            sage: fricas(s).sage()                                              # optional - fricas
            1/3840*n^10 - 5/2304*n^9 + 5/1152*n^8 + 31/5760*n^7 - 229/11520*n^6 - 5/2304*n^5 + 1/36*n^4 - 1/960*n^3 - 1/80*n^2


        Check that :trac:`25224` is fixed::

            sage: integrate(log(x)/(1-x),x,algorithm='fricas')                  # optional - fricas
            dilog(-x + 1)
            sage: fricas(dilog(-x + 1))                                         # optional - fricas
            dilog(x)
            sage: dilog._fricas_()(1.0)                                         # optional - fricas
            1.6449340668_4822643647_24152
            sage: dilog(1.0)
            1.64493406684823

        Check that :trac:`25987` is fixed::

            sage: integrate(lambert_w(x), x, algorithm="fricas")                # optional - fricas
            (x*lambert_w(x)^2 - x*lambert_w(x) + x)/lambert_w(x)

        Check that :trac:`25838` is fixed::

            sage: F = function('f'); f = SR.var('f')
            sage: FF = fricas(F(f)); FF                                         # optional - fricas
            f(f)
            sage: FF.D(f).sage()                                                # optional - fricas
            diff(f(f), f)
            sage: bool(FF.D(f).integrate(f).sage() == F(f))                     # optional - fricas
            True

        Check that :trac:`25602` is fixed::

            sage: r = fricas.integrate(72000/(1+x^5),x).sage()                  # optional - fricas
            sage: n(r.subs(x=5)-r.subs(x=3))                                    # optional - fricas tol 0.1
            193.020947266210

            sage: var("a"); r = fricas.integrate(72000*a^8/(a^5+x^5),x).sage()  # optional - fricas
            a
            sage: n(r.subs(a=1, x=5)-r.subs(a=1, x=3))                          # optional - fricas tol 0.1
            193.020947266268 - 8.73114913702011e-11*I

        Check conversions of sums and products::

            sage: var("k, m, n")                                                # optional - fricas
            (k, m, n)
            sage: fricas("sum(1/factorial(k), k=1..n)").sage()                  # optional - fricas
            sum(1/factorial(_...), _..., 1, n)
            sage: fricas("eval(sum(x/k, k=1..n), x=k)").sage()                  # optional - fricas
            k*harmonic_number(n)
            sage: fricas("product(1/factorial(k), k=1..n)").sage()              # optional - fricas
            1/product(factorial(_...), _..., 1, n)

            sage: f = fricas.guess([sum(1/k, k,1,n) for n in range(10)])[0]; f  # optional - fricas
             n - 1
              --+       1
              >      -------
              --+    s   + 1
            s   = 0   10
             10

            sage: f.sage()                                                      # optional - fricas
            harmonic_number(n)

            sage: f = fricas.guess([0, 1, 3, 9, 33])[0]; f                      # optional - fricas
                    s  - 1
            n - 1    5
             --+    ++-++
             >       | |    p  + 2
             --+     | |     4
            s  = 0  p  = 0
             5       4

            sage: f.sage()                                                      # optional - fricas
            sum(factorial(_... + 1), _..., 0, n - 1)

        Check that :trac:`26746` is fixed::

            sage: _ = var('x, y, z')
            sage: f = sin(x^2) + y^z
            sage: f.integrate(x, algorithm='fricas')                            # optional - fricas
            1/2*sqrt(2)*sqrt(pi)*(sqrt(2)*x*y^z/sqrt(pi) + fresnel_sin(sqrt(2)*x/sqrt(pi)))

            sage: fricas(fresnel_sin(1))                                        # optional - fricas
            fresnelS(1)
            sage: fricas("fresnelS(1.0)")                                       # optional - fricas
            0.4382591473_9035476607_676

            sage: fricas(fresnel_cos(1))                                        # optional - fricas
            fresnelC(1)
            sage: fricas("fresnelC(1.0)")                                       # optional - fricas
            0.7798934003_7682282947_42

        Check that :trac:`17908` is fixed::

            sage: fricas(abs(x)).sage().subs(x=-1783)                           # optional - fricas
            1783

        Check that :trac:`27310` is fixed::

            sage: fricas.set("F", "operator 'f")                                # optional - fricas
            sage: fricas("eval(D(F(x,y), [x, y], [2, 1]), x=x+y)").sage()       # optional - fricas
            D[0, 0, 1](f)(x + y, y)

        Conversion of hypergeometric functions (:trac:`31298`)::

            sage: a,b,c = var("a b c")
            sage: A = hypergeometric([a, b], [c], x)
            sage: fricas(A).sage() - A                                          # optional - fricas
            0
            sage: fricas(A).D(x).sage() - diff(A, x)                            # optional - fricas
            0
        """
        from sage.libs.pynac.pynac import register_symbol
        from sage.symbolic.constants import e, pi, I
        from sage.calculus.functional import diff
        from sage.functions.log import dilog, lambert_w
        from sage.functions.trig import sin, cos, tan, cot, sec, csc
        from sage.functions.hyperbolic import tanh, sinh, cosh, coth, sech, csch
        from sage.functions.other import abs
        from sage.misc.functional import symbolic_sum, symbolic_prod
        from sage.rings.infinity import infinity
        register_symbol(I, {'fricas': '%i'})
        register_symbol(e, {'fricas': '%e'})
        register_symbol(pi, {'fricas': 'pi'})  # fricas uses both pi and %pi
        register_symbol(lambda: infinity, {'fricas': 'infinity'})
        register_symbol(lambda: infinity, {'fricas': 'plusInfinity'})
        register_symbol(lambda: -infinity, {'fricas': 'minusInfinity'})
        register_symbol(cos, {'fricas': 'cos'})
        register_symbol(sin, {'fricas': 'sin'})
        register_symbol(tan, {'fricas': 'tan'})
        register_symbol(cot, {'fricas': 'cot'})
        register_symbol(sec, {'fricas': 'sec'})
        register_symbol(csc, {'fricas': 'csc'})
        register_symbol(tanh, {'fricas': 'tanh'})
        register_symbol(sinh, {'fricas': 'sinh'})
        register_symbol(cosh, {'fricas': 'cosh'})
        register_symbol(coth, {'fricas': 'coth'})
        register_symbol(sech, {'fricas': 'sech'})
        register_symbol(csch, {'fricas': 'csch'})
        register_symbol(lambda x, y: x + y, {'fricas': '+'})
        register_symbol(lambda x, y: x - y, {'fricas': '-'})
        register_symbol(lambda x, y: x * y, {'fricas': '*'})
        register_symbol(lambda x, y: x / y, {'fricas': '/'})
        register_symbol(lambda x, y: x ** y, {'fricas': '^'})
        register_symbol(lambda f, x: diff(f, x), {'fricas': 'D'})
        register_symbol(lambda x, y: x + y * I, {'fricas': 'complex'})
        register_symbol(lambda x: dilog(1 - x), {'fricas': 'dilog'})
        register_symbol(lambda z: lambert_w(z), {'fricas': 'lambertW'})
        register_symbol(abs, {'fricas': 'abs'})
        # construct occurs in the InputForm of hypergeometricF
        register_symbol(lambda *x: x, {'fricas': 'construct'})
        # the following is a hack to deal with
        # integrate(sin((x^2+1)/x),x)::INFORM giving
        # (integral (sin (/ (+ (^ x 2) 1) x)) (:: x Symbol))
        register_symbol(lambda x, y: x, {'fricas': '::'})

        def _convert_eval(f, a, b):
            # it might be that FriCAS also returns a two-argument
            # eval, where the second argument is a list of equations,
            # in which case this function needs to be adapted
            return f.subs({a: b})

        register_symbol(_convert_eval, {'fricas': 'eval'})

        def _convert_sum(x, y):
            v, seg = y.operands()
            a, b = seg.operands()
            return symbolic_sum(x, v, a, b)

        def _convert_prod(x, y):
            v, seg = y.operands()
            a, b = seg.operands()
            return symbolic_prod(x, v, a, b)

        register_symbol(_convert_sum, {'fricas': 'sum'})
        register_symbol(_convert_prod, {'fricas': 'product'})

        def explicitly_not_implemented(*args):
            raise NotImplementedError("the translation of the FriCAS Expression '%s' to sage is not yet implemented" % args)

        register_symbol(lambda *args: explicitly_not_implemented("rootOfADE"), {'fricas': 'rootOfADE'})
        register_symbol(lambda *args: explicitly_not_implemented("rootOfRec"), {'fricas': 'rootOfRec'})

        rootOf = dict()  # (variable, polynomial)
        rootOf_ev = dict()  # variable -> (complex) algebraic number

        def convert_rootOf(x, y):
            if y in rootOf:
                assert rootOf[y] == x
            else:
                rootOf[y] = x
            return y

        register_symbol(convert_rootOf, {'fricas': 'rootOf'})

        ex, _ = FriCASElement._parse_and_eval(fricas_InputForm)
        # postprocessing of rootOf
        from sage.rings.all import QQbar, PolynomialRing
        while rootOf:
            for var, poly in rootOf.items():
                pvars = poly.variables()
                rvars = [v for v in pvars if v not in rootOf_ev]  # remaining variables
                uvars = [v for v in rvars if v in rootOf]  # variables to evaluate
                if len(uvars) == 1:
                    assert uvars[0] == var, "the only variable in uvars should be %s but is %s" % (var, uvars[0])
                    break
            else:
                assert False, "circular dependency in rootOf expression"
            # substitute known roots
            poly = poly.subs(rootOf_ev)
            evars = [v for v in rvars if v not in rootOf]  # extraneous variables
            assert set(evars) == set(poly.variables()).difference([var])
            del rootOf[var]
            if evars:
                # we just need any root per FriCAS specification
                rootOf_ev[var] = poly.roots(var, multiplicities=False)[0]
            else:
                R = PolynomialRing(QQbar, "x")
                # PolynomialRing does not accept variable names with leading underscores
                poly = R(poly.subs({var: R.gen()}))
                # we just need any root per FriCAS specification
                rootOf_ev[var] = poly.roots(multiplicities=False)[0].radical_expression()

        return ex.subs(rootOf_ev)

    def _sage_(self):
        r"""
        Convert ``self`` to a Sage object.

        EXAMPLES:

        Floats::

            sage: fricas(2.1234).sage()                                         # optional - fricas
            2.12340000000000
            sage: _.parent()                                                    # optional - fricas
            Real Field with 53 bits of precision
            sage: a = RealField(100)(pi)                                        # optional - fricas
            sage: fricas(a).sage()                                              # optional - fricas
            3.1415926535897932384626433833
            sage: _.parent()                                                    # optional - fricas
            Real Field with 100 bits of precision
            sage: fricas(a).sage() == a                                         # optional - fricas
            True
            sage: fricas(2.0).sage()                                            # optional - fricas
            2.00000000000000
            sage: _.parent()                                                    # optional - fricas
            Real Field with 53 bits of precision

        Algebraic numbers::

            sage: a = fricas('(1 + sqrt(2))^5'); a                              # optional - fricas
                +-+
            29 \|2  + 41
            sage: b = a.sage(); b                                               # optional - fricas
            82.0121933088198?
            sage: b.radical_expression()                                        # optional - fricas
            29*sqrt(2) + 41

        Integers modulo n::

            sage: fricas("((42^17)^1783)::IntegerMod(5^(5^5))").sage() == Integers(5^(5^5))((42^17)^1783) # optional - fricas
            True

        Matrices over a prime field::

            sage: fricas("matrix [[1::PF 3, 2],[2, 0]]").sage().parent()        # optional - fricas
            Full MatrixSpace of 2 by 2 dense matrices over Finite Field of size 3

        We can also convert FriCAS's polynomials to Sage polynomials::

            sage: a = fricas(x^2 + 1); a.typeOf()                               # optional - fricas
            Polynomial(Integer)
            sage: a.sage()                                                      # optional - fricas
            x^2 + 1
            sage: _.parent()                                                    # optional - fricas
            Univariate Polynomial Ring in x over Integer Ring
            sage: fricas('x^2 + y^2 + 1/2').sage()                              # optional - fricas
            y^2 + x^2 + 1/2
            sage: _.parent()                                                    # optional - fricas
            Multivariate Polynomial Ring in y, x over Rational Field

            sage: fricas("1$Polynomial Integer").sage()                         # optional - fricas
            1

            sage: fricas("x^2/2").sage()                                        # optional - fricas
            1/2*x^2

            sage: x = polygen(QQ, 'x')
            sage: fricas(x+3).sage()                                            # optional - fricas
            x + 3
            sage: fricas(x+3).domainOf()                                        # optional - fricas
            Polynomial(Integer())

            sage: fricas(matrix([[2,3],[4,x+5]])).diagonal().sage()             # optional - fricas
            (2, x + 5)

            sage: f = fricas("(y^2+3)::UP(y, INT)").sage(); f                   # optional - fricas
            y^2 + 3
            sage: f.parent()                                                    # optional - fricas
            Univariate Polynomial Ring in y over Integer Ring

            sage: fricas("(y^2+sqrt 3)::UP(y, AN)").sage()                      # optional - fricas
            y^2 + 1.732050807568878?

        Rational functions::

            sage: fricas("x^2 + 1/z").sage()                                    # optional - fricas
            x^2 + 1/z

        Expressions::

            sage: fricas(pi).sage()                                             # optional - fricas
            pi

            sage: fricas("sin(x+y)/exp(z)*log(1+%e)").sage()                    # optional - fricas
            e^(-z)*log(e + 1)*sin(x + y)

            sage: fricas("factorial(n)").sage()                                 # optional - fricas
            factorial(n)

            sage: fricas("integrate(sin(x+y), x=0..1)").sage()                  # optional - fricas
            -cos(y + 1) + cos(y)

            sage: fricas("integrate(x*sin(1/x), x=0..1)").sage()                # optional - fricas
            'failed'

            sage: fricas("integrate(sin((x^2+1)/x),x)").sage()                  # optional - fricas
            integral(sin((x^2 + 1)/x), x)

        .. TODO::

            - Converting matrices and lists takes much too long.

        Matrices::

            sage: fricas("matrix [[x^n/2^m for n in 0..5] for m in 0..3]").sage()         # optional - fricas, long time
            [      1       x     x^2     x^3     x^4     x^5]
            [    1/2   1/2*x 1/2*x^2 1/2*x^3 1/2*x^4 1/2*x^5]
            [    1/4   1/4*x 1/4*x^2 1/4*x^3 1/4*x^4 1/4*x^5]
            [    1/8   1/8*x 1/8*x^2 1/8*x^3 1/8*x^4 1/8*x^5]

        Lists::

            sage: fricas("[2^n/x^n for n in 0..5]").sage()                      # optional - fricas, long time
            [1, 2/x, 4/x^2, 8/x^3, 16/x^4, 32/x^5]

            sage: fricas("[matrix [[i for i in 1..n]] for n in 0..5]").sage()   # optional - fricas, long time
            [[], [1], [1 2], [1 2 3], [1 2 3 4], [1 2 3 4 5]]

        Error handling::

            sage: s = fricas.guessPade("[fibonacci i for i in 0..10]"); s       # optional - fricas
                n        x
            [[[x ]- ----------]]
                     2
                    x  + x - 1
            sage: s.sage()                                                      # optional - fricas
            Traceback (most recent call last):
            ...
            NotImplementedError: the translation of the FriCAS Expression 'rootOfADE' to sage is not yet implemented

            sage: s = fricas("series(sqrt(1+x), x=0)"); s                       # optional - fricas
                  1     1  2    1  3    5   4    7   5    21   6    33   7    429   8
              1 + - x - - x  + -- x  - --- x  + --- x  - ---- x  + ---- x  - ----- x
                  2     8      16      128      256      1024      2048      32768
            +
               715   9    2431   10      11
              ----- x  - ------ x   + O(x  )
              65536      262144

            sage: s.sage()                                                      # optional - fricas
            Traceback (most recent call last):
            ...
            NotImplementedError: the translation of the FriCAS object
            <BLANKLINE>
                  1     1  2    1  3    5   4    7   5    21   6    33   7    429   8
              1 + - x - - x  + -- x  - --- x  + --- x  - ---- x  + ---- x  - ----- x
                  2     8      16      128      256      1024      2048      32768
            +
               715   9    2431   10      11
              ----- x  - ------ x   + O(x  )
              65536      262144
            <BLANKLINE>
            to sage is not yet implemented:
            An error occurred when FriCAS evaluated 'unparse(...::InputForm)':
            <BLANKLINE>
               Cannot convert the value from type Any to InputForm .
        """
        from sage.rings.all import PolynomialRing, RDF, I
        from sage.rings.real_mpfr import RealField
        from sage.symbolic.ring import SR
        from sage.matrix.constructor import matrix
        from sage.modules.free_module_element import vector
        from sage.structure.factorization import Factorization
        from sage.misc.sage_eval import sage_eval

        # TODO: perhaps we should translate the type first?
        # TODO: perhaps we should get the InputForm as SExpression?

        # remember: fricas.new gives a FriCASElement

        # the coercion to Any gets rid of the Union domain
        P = self._check_valid()
        domain = P.new("dom((%s)::Any)" % self._name)  # domain is now a fricas SExpression

        # first translate dummy domains such as "failed".  We must
        # not recurse here!
        if P.get_boolean("string?(%s)" % domain._name):
            return P.get_string("string(%s)" % domain._name)

        # now translate domains which cannot be coerced to InputForm,
        # or where we do not need it.
        head = str(domain.car())
        if head == "Record":
            fields = fricas("[string symbol(e.2) for e in rest destruct %s]" % domain._name).sage()
            return {field: self.elt(field).sage() for field in fields}

        if head == "List":
            n = P.get_integer('#(%s)' % self._name)
            return [self.elt(k).sage() for k in range(1, n + 1)]

        if head == "Vector" or head == "DirectProduct":
            n = P.get_integer('#(%s)' % self._name)
            return vector([self.elt(k).sage() for k in range(1, n + 1)])

        if head == "Matrix":
            base_ring = self._get_sage_type(domain[1])
            rows = self.listOfLists().sage()
            return matrix(base_ring, rows)

        if head == "Fraction":
            return self.numer().sage() / self.denom().sage()

        if head == "Complex":
            return self.real().sage() + self.imag().sage() * I

        if head == "Factored":
            l = P.new('[[f.factor, f.exponent] for f in factors(%s)]' % self._name).sage()
            return Factorization([(p, e) for p, e in l])

        if head == "UnivariatePolynomial":
            base_ring = self._get_sage_type(domain[2])
            vars = str(domain[1])
            R = PolynomialRing(base_ring, vars)
            return R([self.coefficient(i).sage()
                      for i in range(ZZ(self.degree()) + 1)])

        # finally translate domains with InputForm
        try:
            unparsed_InputForm = P.get_unparsed_InputForm(self._name)
        except RuntimeError as error:
            raise NotImplementedError("the translation of the FriCAS object\n\n%s\n\nto sage is not yet implemented:\n%s" % (self, error))
        if head == "Boolean":
            return unparsed_InputForm == "true"

        if head in ["Integer", "NonNegativeInteger", "PositiveInteger"]:
            return ZZ(unparsed_InputForm)

        if head == "String":
            return unparsed_InputForm

        if head == "Float":
            # Warning: precision$Float gives the current precision,
            # whereas length(mantissa(self)) gives the precision of
            # self.
            prec = max(P.new("length mantissa(%s)" % self._name).sage(), 53)
            R = RealField(prec)
            x, e, b = unparsed_InputForm.lstrip('float(').rstrip(')').split(',')
            return R(ZZ(x) * ZZ(b)**ZZ(e))

        if head == "DoubleFloat":
            return RDF(unparsed_InputForm)

        if head == "AlgebraicNumber":
            s = unparsed_InputForm[:-len("::AlgebraicNumber()")]
            return sage_eval("QQbar(" + s + ")")

        if head == "IntegerMod" or head == "PrimeField":
            # one might be tempted not to go via InputForm here, but
            # it turns out to be safer to do it.
            n = unparsed_InputForm[len("index("):]
            n = n[:n.find(")")]
            return self._get_sage_type(domain)(n)

        if head == "Polynomial":
            base_ring = self._get_sage_type(domain[1])
            # Polynomial Complex is translated into SR
            if base_ring is SR:
                return FriCASElement._sage_expression(P.get_InputForm(self._name))

            # the following is a bad hack, we should be getting a list here
            vars = P.get_unparsed_InputForm("variables(%s)" % self._name)[1:-1]
            if vars == "":
                return base_ring(unparsed_InputForm)
            else:
                R = PolynomialRing(base_ring, vars)
                return R(unparsed_InputForm)

        if head in ["OrderedCompletion", "OnePointCompletion"]:
            # it would be more correct to get the type parameter
            # (which might not be Expression Integer) and recurse
            return FriCASElement._sage_expression(P.get_InputForm(self._name))

        if head == "Expression" or head == "Pi":
            # we treat Expression Integer and Expression Complex
            # Integer just the same
            return FriCASElement._sage_expression(P.get_InputForm(self._name))

        if head == 'DistributedMultivariatePolynomial':
            base_ring = self._get_sage_type(domain[2])
            vars = domain[1].car()
            R = PolynomialRing(base_ring, vars)
            return R(unparsed_InputForm)

        raise NotImplementedError("the translation of the FriCAS object %s to sage is not yet implemented" % (unparsed_InputForm))


@instancedoc
class FriCASFunctionElement(FunctionElement):
    def __init__(self, object, name):
        """
        Make FriCAS operation names valid python function identifiers.

        TESTS::

            sage: a = fricas('"Hello"')                                         # optional - fricas
            sage: a.upperCase_q                                                 # optional - fricas
            upperCase?
            sage: a.upperCase_e                                                 # optional - fricas
            upperCase!
            sage: a.upperCase_e()                                               # optional - fricas
            "HELLO"

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        FunctionElement.__init__(self, object, name)


@instancedoc
class FriCASExpectFunction(ExpectFunction):
    def __init__(self, parent, name):
        """
        Translate the pythonized function identifier back to a FriCAS
        operation name.

        TESTS::

            sage: fricas.upperCase_q                                            # optional - fricas
            upperCase?
            sage: fricas.upperCase_e                                            # optional - fricas
            upperCase!

        """
        if name.endswith("_q"):
            name = name[:-2] + "?"
        elif name.endswith("_e"):
            name = name[:-2] + "!"
        ExpectFunction.__init__(self, parent, name)


def is_FriCASElement(x):
    """
    Return ``True`` if ``x`` is of type :class:`FriCASElement`.

    EXAMPLES::

        sage: from sage.interfaces.fricas import is_FriCASElement               # optional - fricas
        sage: is_FriCASElement(fricas(2))                                       # optional - fricas
        True
        sage: is_FriCASElement(2)                                               # optional - fricas
        False
    """
    return isinstance(x, FriCASElement)


fricas = FriCAS()


def reduce_load_fricas():
    """
    Return the FriCAS interface object defined in
    :sage.interfaces.fricas.

    EXAMPLES::

        sage: from sage.interfaces.fricas import reduce_load_fricas             # optional - fricas
        sage: reduce_load_fricas()                                              # optional - fricas
        FriCAS
    """
    return fricas


def fricas_console():
    """
    Spawn a new FriCAS command-line session.

    EXAMPLES::

        sage: fricas_console()                                                  # not tested
                         FriCAS (AXIOM fork) Computer Algebra System
                                    Version: FriCAS 1.0.5
                     Timestamp: Thursday February 19, 2009 at 06:57:33
        -----------------------------------------------------------------------------
           Issue )copyright to view copyright notices.
           Issue )summary for a summary of useful system commands.
           Issue )quit to leave AXIOM and return to shell.
        -----------------------------------------------------------------------------
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%fricas magics instead.')
    os.system('fricas -nox')


def __doctest_cleanup():
    """
    EXAMPLES::

        sage: from sage.interfaces.fricas import __doctest_cleanup              # optional - fricas
        sage: a = FriCAS()                                                      # optional - fricas
        sage: two = a(2)                                                        # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        True
        sage: __doctest_cleanup()                                               # optional - fricas
        sage: a.is_running()                                                    # optional - fricas
        False
    """
    import sage.interfaces.quit
    sage.interfaces.quit.expect_quitall()
