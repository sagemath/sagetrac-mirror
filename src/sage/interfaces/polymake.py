r"""
Interface to polymake

"""


#*****************************************************************************
#       Copyright (C) 2017 Simon King <simon.king@uni-jena.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from __future__ import print_function
from __future__ import absolute_import

import os
import re
import sys
import six

from sage.structure.parent import Parent
from .expect import console, Expect, ExpectElement, ExpectFunction, FunctionElement

from sage.env import SAGE_EXTCODE, DOT_SAGE
from sage.misc.misc import get_verbose
from sage.misc.cachefunc import cached_method
from sage.interfaces.tab_completion import ExtraTabCompletion

import pexpect
from random import randrange

from time import sleep
from six.moves import range
from six import reraise as raise_
import warnings

_name_pattern = re.compile('SAGE[0-9]+')

_available_polymake_answers = {
    0: "returns prompt",
    1: "returns continuation prompt",
    2: "requests interactive input",
    3: "kills computation",
    4: "raises error",
    5: "issues warning",
    6: "shows additional information",
    7: "lost connection",
    8: "fails to respond timely"
        }

class PolymakeError(RuntimeError):
    """
    Raised if polymake yields an error message.

    TESTS::

        sage: polymake.eval('print foo;')    # optional polymake
        Traceback (most recent call last):
        ...
        PolymakeError: Unquoted string "foo" may clash with future reserved word...

    """
    pass

def polymake_console(command=''):
    """
    Spawn a new polymake command-line session.

    EXAMPLES::

        sage: from sage.interfaces.polymake import polymake_console
        sage: polymake_console()        # not tested
        Welcome to polymake version ...
        ...
        Ewgenij Gawrilow, Michael Joswig (TU Berlin)
        http://www.polymake.org

        This is free software licensed under GPL; see the source for copying conditions.
        There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

        Press F1 or enter 'help;' for basic instructions.

        Application polytope currently uses following third-party software packages:
        4ti2, bliss, cdd, latte, libnormaliz, lrs, permlib, ppl, sketch, sympol, threejs, tikz, topcom, tosimplex
        For more details:  show_credits;
        polytope >

    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%polymake magics instead.')
    os.system(command or os.getenv('SAGE_POLYMAKE_COMMAND') or 'polymake')

class Polymake(ExtraTabCompletion, Expect):
    r"""
    Interface to the polymake interpreter.

    In order to use this interface, you need to either install the
    optional polymake package for Sage, or install polymake system-wide
    on your computer.

    Type ``polymake.[tab]`` for a list of most functions
    available from your polymake install. Type
    ``polymake.Function?`` for polymake's help about a given ``Function``
    Type ``polymake(...)`` to create a new Magma
    object, and ``polymake.eval(...)`` to run a string using
    polymake and get the result back as a string.

    EXAMPLES::

        sage: p = polymake.rand_sphere(4, 20, seed=5)       # optional - polymake
        sage: p                                             # optional - polymake
        Random spherical polytope of dimension 4; seed=5...
        sage: set_verbose(3)
        sage: p.H_VECTOR                                    # optional - polymake
        used package ppl
          The Parma Polyhedra Library (PPL): A C++ library for convex polyhedra
          and other numerical abstractions.
          http://www.cs.unipr.it/ppl/
        1 16 47 16 1
        sage: set_verbose(0)
        sage: p.F_VECTOR                                    # optional - polymake
        20 101 162 81
        sage: print(p.F_VECTOR._sage_doc_())                # optional - polymake # random
        property_types/Algebraic Types/Vector:
         A type for vectors with entries of type Element.

         You can perform algebraic operations such as addition or scalar multiplication.

         You can create a new Vector by entering its elements, e.g.:
            $v = new Vector<Int>(1,2,3);
         or
            $v = new Vector<Int>([1,2,3]);

    .. automethod:: _eval_line
    """
    def __init__(self, script_subdirectory=None,
                 logfile=None, server=None,server_tmpdir=None,
                 seed=None, command=None):
        """
        TESTS::

            sage: from sage.interfaces.polymake import Polymake
            sage: Polymake()
            Polymake
            sage: Polymake().is_running()
            False

        """
        if command is None:
            command = "env TERM=dumb {}".format(os.getenv('SAGE_POLYMAKE_COMMAND') or 'polymake')
        Expect.__init__(self,
                        name="polymake",
                        command=command,
                        prompt="polytope > ",
                        server=server,
                        server_tmpdir=server_tmpdir,
                        script_subdirectory=script_subdirectory,
                        restart_on_ctrlc=False,
                        logfile=logfile,
                        eval_using_file_cutoff=1024)   # > 1024 causes hangs

        self._seed = seed
        self.__tab_completion = {}

    @cached_method
    def version(self):
        """
        Version of the polymake installation.

        EXAMPLES::

            sage: polymake.version()               # optional - polymake
            '3...'

        TESTS::

            sage: from sage.interfaces.polymake import Polymake
            sage: Polymake(command='foobar').version()
            Traceback (most recent call last):
            ...
            RuntimeError: unable to start polymake because the command 'foobar' failed:
            The command was not found or was not executable: foobar.
            Please install the optional polymake package for sage (but read its SPKG.txt first!)
            or install polymake system-wide

        """
        return self.get('$Polymake::Version')

    # Pickling etc

    def __reduce__(self):
        """
        EXAMPLES::

            sage: loads(dumps(polymake)) is polymake
            True

        """
        return reduce_load_Polymake, tuple([])

    def _object_class(self):
        """
        Return the class by which elements in this interface are implemented.

        TESTS::

            sage: C = polymake('cube(3)')  # indirect doctest   # optional - polymake
            sage: C                                             # optional - polymake
            cube of dimension 3
            sage: type(C)                                       # optional - polymake
            <class 'sage.interfaces.polymake.PolymakeElement'>

        """
        return PolymakeElement

    def _function_element_class(self):
        """
        Return the class by which member functions of this interface are implemented.

        TESTS:

        We use ellipses in the tests, to make it more robust against future
        changes in polymake::

            sage: p = polymake.rand_sphere(4, 20, seed=5)    # optional - polymake
            sage: p.get_schedule                            # optional - polymake  # indirect doctest
            Member function 'get_schedule' of Polymake::polytope::Polytope__Rational object
            sage: p.get_schedule('F_VECTOR')                # optional - polymake  # random
            CONE_DIM : RAYS | INPUT_RAYS
            precondition : BOUNDED ( POINTED : )
            POINTED :
            N_INPUT_RAYS : INPUT_RAYS
            precondition : N_RAYS | N_INPUT_RAYS ( ppl.convex_hull.primal: FACETS, LINEAR_SPAN : RAYS | INPUT_RAYS )
            sensitivity check for FacetPerm
            ppl.convex_hull.primal: FACETS, LINEAR_SPAN : RAYS | INPUT_RAYS
            INPUT_RAYS_IN_FACETS : INPUT_RAYS, FACETS
            sensitivity check for VertexPerm
            RAYS_IN_FACETS, RAYS, LINEALITY_SPACE : INPUT_RAYS_IN_FACETS, INPUT_RAYS
            GRAPH.ADJACENCY : RAYS_IN_FACETS
            DUAL_GRAPH.ADJACENCY : RAYS_IN_FACETS
            N_EDGES : ADJACENCY ( applied to GRAPH )
            N_EDGES : ADJACENCY ( applied to DUAL_GRAPH )
            precondition : POINTED ( LINEALITY_DIM, LINEALITY_SPACE : )
            LINEALITY_DIM, LINEALITY_SPACE :
            COMBINATORIAL_DIM : CONE_DIM, LINEALITY_DIM
            N_RAYS : RAYS
            N_FACETS : FACETS
            precondition : COMBINATORIAL_DIM ( F_VECTOR : N_FACETS, N_RAYS, GRAPH.N_EDGES, DUAL_GRAPH.N_EDGES, COMBINATORIAL_DIM )
            F_VECTOR : N_FACETS, N_RAYS, GRAPH.N_EDGES, DUAL_GRAPH.N_EDGES, COMBINATORIAL_DIM

        """
        return PolymakeFunctionElement

    def function_call(self, function, args=None, kwds=None):
        """
        EXAMPLES::

            sage: polymake.rand_sphere(4, 30, seed=15)           # optional - polymake  # indirect doctest
            Random spherical polytope of dimension 4; seed=15...

        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        s = self._function_call_string(function,
                                       [s.name() for s in args],
                                       ['%s=>%s'%(key,value.name()) for key, value in kwds.items()])
        return self(s)

    def _function_call_string(self, function, args, kwds):
        """
        Returns the string used to make function calls.

        EXAMPLES::

            sage: polymake._function_call_string('cube', ['2','7','3'], ['group=>1']) # optional - polymake
            'cube(2,7,3, group=>1);'
            sage: c = polymake('cube(2,7,3, group=>1)')                 # optional - polymake
            sage: c.VERTICES                                            # optional - polymake
            1 3 3
            1 7 3
            1 3 7
            1 7 7
            sage: c.GROUP                                               # optional - polymake
            full combinatorial group on facets...

        """
        if kwds:
            if args:
                return "%s(%s, %s);"%(function, ",".join(list(args)), ",".join(list(kwds)))
            return "%s(%s);"%(function, ",".join(list(kwds)))
        return "%s(%s);"%(function, ",".join(list(args)))

    def console(self):
        """
        Raise an error, pointing to :meth:`~sage.interfaces.interface.Interface.interact` and :func:`polymake_console`.

        EXAMPLES::

            sage: polymake.console()
            Traceback (most recent call last):
            ...
            NotImplementedError: Please use polymake_console() function or the .interact() method

        """
        raise NotImplementedError("Please use polymake_console() function or the .interact() method")

    # Methods concerning interface communication

    def _install_hints(self):
        """
        TESTS::

            sage: print(polymake._install_hints())
            Please install the optional polymake package for sage (but read its SPKG.txt first!)
            or install polymake system-wide

        """
        return "Please install the optional polymake package for sage  (but read its SPKG.txt first!)"+os.linesep+"or install polymake system-wide"

    def _start(self, alt_message=None):
        """
        Start the polymake interface in the application "polytope".

        NOTE:

        There should be no need to call this explicitly.

        TESTS::

            sage: polymake.application('fan')               # optional - polymake
            sage: 'normal_fan' in dir(polymake)             # optional - polymake
            True
            sage: polymake.quit()                           # optional - polymake
            sage: polymake._start()                         # optional - polymake

        Since 'normal_fan' is not defined in the polymake application 'polytope',
        we now get
        ::

            sage: 'normal_fan' in dir(polymake)             # optional - polymake
            False

        """
        if not self.is_running():
            self._change_prompt("polytope > ")
            Expect._start(self, alt_message=None)
        self.application("polytope")
        self.eval('use Scalar::Util qw(reftype);')
        self.eval('use Scalar::Util qw(blessed);')
        self.eval('use File::Slurp;')

    def _quit_string(self):
        """
        TEST::

            sage: polymake._quit_string()
            'exit;'
        """
        return "exit;"

    def _assign_symbol(self):
        """
        TEST::

            sage: polymake._assign_symbol()
            '='
        """
        return "="

    def _equality_symbol(self):
        """
        TEST::

            sage: polymake._equality_symbol()
            '=='
        """
        return "=="

    def _read_in_file_command(self, filename):
        """
        TEST::

            sage: polymake._read_in_file_command('foobar')
            'eval read_file "foobar";\n'

        Force use of file::

            sage: L = polymake([42] * 400)                      # optional - polymake
            sage: len(L)                                        # optional - polymake
            400

        Just below standard file cutoff of 1024::

            sage: L = polymake([42] * 84)                       # optional - polymake
            sage: len(L)                                        # optional - polymake
            84

        """
        return 'eval read_file "{}";\n'.format(filename)

    def _keyboard_interrupt(self):
        """
        Interrupt a computation with <Ctrl-c>

        TEST:

        For reasons that are not clear to the author, the following test
        is very flaky. Therefore, this test is marked as "not tested".

            sage: c = polymake.cube(15)                         # optional - polymake
            sage: alarm(1)                                      # not tested
            sage: try:                                          # not tested # indirect doctest
            ....:     c.F_VECTOR
            ....: except KeyboardInterrupt:
            ....:     pass
            Interrupting Polymake...
            doctest:warning
            ...
            RuntimeWarning: We ignore that Polymake issues warning during keyboard interrupt
            doctest:warning
            ...
            RuntimeWarning: We ignore that Polymake raises error during keyboard interrupt

        Afterwards, the interface should still be running.  ::

            sage: c.N_FACETS                                    # optional - polymake
            30

        """
        if not self.is_running():
            raise KeyboardInterrupt
        print("Interrupting %s..." % self)
        while True:
            try:
                self._expect.send(chr(3))
            except pexpect.ExceptionPexpect as msg:
                raise pexpect.ExceptionPexpect("THIS IS A BUG -- PLEASE REPORT. This should never happen.\n" + msg)
            sleep(0.1)
            i = self._expect.expect_list(self._prompt, timeout=1)
            if i==0:
                break
            elif i==7:  # EOF
                warnings.warn("Polymake {} during keyboard interrupt".format(_available_polymake_answers[i]), RuntimeWarning)
                self._crash_msg()
                self.quit()
            elif i==8:  # Timeout
                self.quit()
                raise RuntimeError("{} interface is not responding. We closed it".format(self))
            elif i!=3: # Anything but a "computation killed"
                warnings.warn("We ignore that {} {} during keyboard interrupt".format(self, _available_polymake_answers[i]), RuntimeWarning)
        raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)

    def _synchronize(self):
        """
        TEST::

            sage: Q = polymake.cube(4)                          # optional - polymake
            sage: polymake('"ok"')                              # optional - polymake
            ok
            sage: polymake._expect.sendline()                   # optional - polymake
            1

        Now the interface is badly out of sync::

            sage: polymake('"foobar"')                          # optional - polymake
            <repr(<sage.interfaces.polymake.PolymakeElement at ...>) failed:
            PolymakeError: Can't locate object method "description" via package "1"
            (perhaps you forgot to load "1"?)...>
            sage: Q.typeof()                                    # optional - polymake
            ('foobar...', 'Polymake::polytope::Polytope__Rational')
            sage: Q.typeof.clear_cache()                        # optional - polymake

        After synchronisation, things work again as expected::

            sage: polymake._synchronize()                       # optional - polymake
            doctest:warning
            ...
            UserWarning: Polymake seems out of sync:
            The expected output did not appear before reaching the next prompt.
            sage: polymake('"back to normal"')                  # optional - polymake
            back to normal
            sage: Q.typeof()                                    # optional - polymake
            ('Polymake::polytope::Polytope__Rational', 'ARRAY')

        """
        if not self.is_running():
            return
        rnd = randrange(2147483647)
        res = str(rnd+1)
        cmd='print 1+{};'+self._expect.linesep
        self._sendstr(cmd.format(rnd))
        pat = self._expect.expect(self._prompt,timeout=0.5)
        # 0: normal prompt
        # 1: continuation prompt
        # 2: user input expected when requestion "help"
        # 3: what we are looking for when interrupting a computation
        # 4: error
        # 5: warning
        # 6: anything but an error or warning, thus, an information
        # 7: unexpected end of the stream
        # 8: (expected) timeout
        if pat == 8: # timeout
            warnings.warn("{} unexpectedly {} during synchronisation.".format(self, _available_polymake_answers[pat]), RuntimeWarning)
            self.interrupt()
            # ... but we continue, as that probably means we currently are at the end of the buffer
        elif pat == 7: # EOF
            self._crash_msg()
            self.quit()
        elif pat == 0:
            # We got the right prompt, but perhaps in a wrong position in the stream
            # The result of the addition should appear *before* our prompt
            if not res in self._expect.before:
                try:
                    warnings.warn("{} seems out of sync: The expected output did not appear before reaching the next prompt.".format(self))
                    while True:
                        i = self._expect.expect_list(self._prompt, timeout=0.1)
                        if i==8: # This time, we do expect a timeout
                            return
                        elif i>0:
                            raise RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[i]))
                except pexpect.TIMEOUT:
                    warnings.warn("A timeout has occured when synchronising {}.".format(self), RuntimeWarning)
                    self._interrupt()
                except pexpect.EOF:
                    self._crash_msg()
                    self.quit()
            else:
                return
        else:
            raise RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))

    def _next_var_name(self):
        r"""
        Returns the next unused variable name.

        TEST::

            sage: print(polymake._next_var_name())
            SAGE...

        """
        if len(self._available_vars) != 0:
            return self._available_vars.pop(0)
        try:
            self.__seq += 1
        except AttributeError:
            self.__seq = 0
        return r'SAGE%s'%self.__seq

    def clear(self, var):
        """
        Clear the variable named var.

        NOTE:

        This is implicitly done when deleting an element in the interface.

        TESTS::

            sage: c = polymake.cube(15)                 # optional - polymake
            sage: polymake._available_vars = []         # optional - polymake
            sage: old = c._name                         # optional - polymake
            sage: del c                                 # optional - polymake  # indirect doctest
            sage: len(polymake._available_vars)         # optional - polymake
            1
            sage: polymake._next_var_name() in old      # optional - polymake
            True

        """
        self._available_vars.append(_name_pattern.search(var).group())

    def _create(self, value, name=None):
        """
        Assign a value to a name in the polymake interface.

        INPUT:

        - ``value``, a string: Polymake command (or value) whose result
          is to be assigned to a variable
        - ``name``, optional string: If given, the new variable has this name.
          Otherwise, the name is automatically generated.

        RETURN:

        The command by which the assigned value can now be retrieved.

        NOTE:

        In order to overcome problems with the perl programming language,
        we store *all* data as arrays. If the given value is an array
        of length different from one, then the new variable contains that
        array. Otherwise, the new variable is an array of length one whose
        only entry is the given value, which has to be a scalar (which
        also includes Perl references). In other words, perl hashes
        are not suitable.

        EXAMPLES::

            sage: polymake._create("('foo', 'bar')", name="my_array")   # optional - polymake
            '@my_array'
            sage: print(polymake.eval('print join(", ", @my_array);'))  # optional - polymake
            foo, bar
            sage: polymake._create("'foobar'", name="my_string")        # optional - polymake
            '$my_string[0]'
            sage: print(polymake.eval('print $my_string[0];'))          # optional - polymake
            foobar

        """
        name = self._next_var_name() if name is None else name
        self.set(name, value)
        # If value is a list, then @name is now equal to that list.
        # Otherwise, value is obtained by $name[0]. So, we modify
        # the name returned by _create so that it can be used to
        # access the wrapped value.
        if self.eval('print scalar @{};'.format(name)).strip() == '1':
            return '$'+name+'[0]'
        return '@'+name

    def set(self, var, value):
        """
        Set the variable var to the given value.

        Eventually, ``var`` is a reference to ``value``.

        .. WARNING::

            This method, although it doesn't start with an underscore, is
            an internal method and not part of the interface. So, please do
            not try to call it explicitly. Instead, use the polymake interface
            as shown in the examples.

        REMARK:

        Polymake's user language is Perl. In Perl, if one wants to assign
        the return value of a function to a variable, the syntax to do so
        depends on the type of the return value. While this is fine in
        compiled code, it seems quite awkward in user interaction.

        To make this polymake pexpect interface a bit more user friendly,
        we treat *all* variables as arrays. A scalar value (most typically
        a reference) is thus interpreted as the only item in an array of
        length one. It is, of course, possible to use the interface without
        knowing these details.

        EXAMPLES::

            sage: c = polymake('cube(3)')                       # optional - polymake # indirect doctest
            sage: d = polymake.cube(3)                          # optional - polymake

        Equality is, for "big" objects such as polytopes, comparison by
        identity::

            sage: c == d                                        # optional - polymake
            False

        However, the list of vertices is equal::

            sage: c.VERTICES == d.VERTICES                      # optional - polymake
            True

        TESTS:

        The following shows how polymake variables are wrapped in our interface.
        It should, however, **never** be needed to do the following
        *explicitly*::

            sage: polymake.set('myvar', 'cube(3)')              # optional - polymake
            sage: polymake.get('$myvar[0]')                     # optional - polymake
            'Polymake::polytope::Polytope__Rational=ARRAY(...)'

        The following tests against :trac:`22658`::

            sage: P = polymake.new_object("Polytope", FACETS=[[12, -2, -3, -5, -8, -13, -21, -34, -55],   # optional - polymake
            ....:  [0, 1, 0, 0, 0, 0, 0, 0, 0],
            ....:  [0, 0, 0, 0, 0, 0, 0, 0, 1],
            ....:  [0, 0, 0, 0, 0, 0, 0, 1, 0],
            ....:  [0, 0, 0, 0, 0, 0, 1, 0, 0],
            ....:  [0, 0, 0, 0, 0, 1, 0, 0, 0],
            ....:  [0, 0, 0, 0, 1, 0, 0, 0, 0],
            ....:  [0, 0, 0, 1, 0, 0, 0, 0, 0],
            ....:  [0, 0, 1, 0, 0, 0, 0, 0, 0]])
            sage: P.VERTICES                        # optional - polymake
            1 6 0 0 0 0 0 0 0
            1 0 4 0 0 0 0 0 0
            1 0 0 0 0 0 0 0 0
            1 0 0 12/5 0 0 0 0 0
            1 0 0 0 0 0 0 0 12/55
            1 0 0 0 0 0 0 6/17 0
            1 0 0 0 0 0 4/7 0 0
            1 0 0 0 0 12/13 0 0 0
            1 0 0 0 3/2 0 0 0 0
            sage: P.F_VECTOR                        # optional - polymake
            9 36 84 126 126 84 36 9

        """
        if isinstance(value, six.string_types):
            value = value.strip().rstrip(';').strip()
        cmd = '@%s%s(%s);'%(var,self._assign_symbol(), value)
        self.eval(cmd)

    def get(self, cmd):
        """
        Return the string representation of an object in the polymake interface.

        EXAMPLES::

            sage: polymake.get('cube(3)')                     # optional - polymake
            'Polymake::polytope::Polytope__Rational=ARRAY(...)'

        Note that the above string representation is what polymake provides.
        In our interface, we use what polymake calls a "description"::

            sage: polymake('cube(3)')                         # optional - polymake
            cube of dimension 3


        """
        return self.eval("print {};".format(cmd)).strip()

    def help(self, topic, pager=True):
        """
        Displays polymake's help on a given topic, as a string.

        INPUT:

        - ``topic``, a string
        - ``pager``, optional bool, default ``True``: When True, display help, otherwise return as a string.

        EXAMPLES::

            sage: print(polymake.help('Polytope', pager=False))         # optional - polymake # random
            objects/Polytope:
             Not necessarily bounded or unbounded polyhedron.
             Nonetheless, the name "Polytope" is used for two reasons:
             Firstly, combinatorially we always deal with polytopes; see the description of VERTICES_IN_FACETS for details.
             The second reason is historical.
             We use homogeneous coordinates, which is why Polytope is derived from Cone.
             Note that a pointed polyhedron is projectively equivalent to a polytope.
             Scalar is the numeric data type used for the coordinates.

        In some cases, polymake expects user interaction to choose from
        different available help topics. In these cases, a warning is given,
        and the available help topics are displayed resp. printed, without
        user interaction::

            sage: polymake.help('TRIANGULATION')                        # optional - polymake # random
            doctest:warning
            ...
            UserWarning: Polymake expects user interaction. We abort and return the options that Polymake provides.
            There are 5 help topics matching 'TRIANGULATION':
            1: objects/Visualization/Visual::Polytope/methods/TRIANGULATION
            2: objects/Visualization/Visual::PointConfiguration/methods/TRIANGULATION
            3: objects/Cone/properties/Triangulation and volume/TRIANGULATION
            4: objects/PointConfiguration/properties/Triangulation and volume/TRIANGULATION
            5: objects/Polytope/properties/Triangulation and volume/TRIANGULATION

        If an unkown help topic is requested, a :class:`PolymakeError` results::

            sage: polymake.help('Triangulation')                        # optional - polymake
            Traceback (most recent call last):
            ...
            PolymakeError: unknown help topic 'Triangulation'

        """
        H = self.eval('help("{}");\n'.format(topic))
        if pager:
            from IPython.core.page import page
            page(H, start = 0)
        else:
            return H

    def _eval_line(self, line, allow_use_file=True, wait_for_prompt=True, restart_if_needed=True, **kwds):
        r"""
        Evaluate a command.

        INPUT:

        - ``line``, a command (string) to be evaluated
        - ``allow_use_file`` (optional bool, default ``True``), whether or not
          to use a file if the line is very long.
        - ``wait_for_prompt`` (optional, default ``True``), whether or not
          to wait before polymake returns a prompt. If it is a string, it is considered
          as alternative prompt to be waited for.
        - ``restart_if_needed`` (optional bool, default ``True``), whether or
          not to restart polymake in case something goes wrong
        - further optional arguments (e.g., timeout) that will be passed to
          :meth:`pexpect.pty_spawn.spawn.expect`. Note that they are ignored
          if the line is too long and thus is evaluated via a file. So,
          if a timeout is defined, it should be accompanied by ``allow_use_file=False``.

        Different reaction types of polymake, including warnings, comments,
        errors, request for user interaction, and yielding a continuation prompt,
        are taken into account.

        Usually, this method is indirectly called via :meth:`~sage.interfaces.expect.Expect.eval`.

        EXAMPLES::

            sage: p = polymake.cube(3)              # optional - polymake  # indirect doctest

        Here we see that remarks printed by polymake are displayed if
        the verbosity is positive::

            sage: set_verbose(1)
            sage: p.N_LATTICE_POINTS                # optional - polymake
            used package latte
              LattE (Lattice point Enumeration) is a computer software dedicated to the
              problems of counting lattice points and integration inside convex polytopes.
              Copyright by Matthias Koeppe, Jesus A. De Loera and others.
              http://www.math.ucdavis.edu/~latte/
            27
            sage: set_verbose(0)

        If polymake raises an error, the polymake *interface* raises
        a :class:`PolymakeError`::

            sage: polymake.eval('FOOBAR(3);')       # optional - polymake
            Traceback (most recent call last):
            ...
            PolymakeError: Undefined subroutine &Polymake::User::FOOBAR called...

        If a command is incomplete, then polymake returns a continuation
        prompt. In that case, we raise an error::

            sage: polymake.eval('print 3')          # optional - polymake
            Traceback (most recent call last):
            ...
            SyntaxError: Incomplete polymake command 'print 3'
            sage: polymake.eval('print 3;')         # optional - polymake
            '3'

        However, if the command contains line breaks but eventually is complete,
        no error is raised::

            sage: print(polymake.eval('$tmp="abc";\nprint $tmp;'))  # optional - polymake
            abc

        When requesting help, polymake sometimes expect the user to choose
        from a list. In that situation, we abort with a warning, and show
        the list from which the user can chose; we could demonstrate this using
        the :meth:`help` method, but here we use an explicit code evaluation::

            sage: print(polymake.eval('help "TRIANGULATION";'))     # optional - polymake # random
            doctest:warning
            ...
            UserWarning: Polymake expects user interaction. We abort and return
            the options that Polymake provides.
            There are 5 help topics matching 'TRIANGULATION':
            1: objects/Cone/properties/Triangulation and volume/TRIANGULATION
            2: objects/Polytope/properties/Triangulation and volume/TRIANGULATION
            3: objects/Visualization/Visual::PointConfiguration/methods/TRIANGULATION
            4: objects/Visualization/Visual::Polytope/methods/TRIANGULATION
            5: objects/PointConfiguration/properties/Triangulation and volume/TRIANGULATION

        By default, we just wait until polymake returns a result. However,
        it is possible to explicitly set a timeout. The following usually does
        work in an interactive session and often in doc tests, too. However,
        sometimes it hangs, and therefore we remove it from the tests, for now::

            sage: c = polymake.cube(15)             # optional - polymake
            sage: polymake.eval('print {}->F_VECTOR;'.format(c.name()), timeout=1) # optional - polymake # not tested
            Traceback (most recent call last):
            ...
            RuntimeError: Polymake fails to respond timely

        We verify that after the timeout, polymake is still able to give answers::

            sage: c                                 # optional - polymake
            cube of dimension 15
            sage: c.N_VERTICES                      # optional - polymake
            32768

        Note, however, that the recovery after a timeout is not perfect.
        It may happen that in some situation the interface collapses and
        thus polymake would automatically be restarted, thereby losing all
        data that have been computed before.

        """
        line = line.strip()
        if allow_use_file and wait_for_prompt and self._eval_using_file_cutoff and len(line) > self._eval_using_file_cutoff:
            return self._eval_line_using_file(line)
        try:
            if not self.is_running():
                self._start()
            E = self._expect
            try:
                if len(line) >= 4096:
                    raise RuntimeError("Sending more than 4096 characters with %s on a line may cause a hang and you're sending %s characters"%(self, len(line)))
                E.sendline(line)
                if not wait_for_prompt:
                    return ''

            except OSError as msg:
                if restart_if_needed:
                    # The subprocess most likely crashed.
                    # If it's really still alive, we fall through
                    # and raise RuntimeError.
                    if sys.platform.startswith('sunos'):
                        # On (Open)Solaris, we might need to wait a
                        # while because the process might not die
                        # immediately. See Trac #14371.
                        for t in [0.5, 1.0, 2.0]:
                            if E.isalive():
                                time.sleep(t)
                            else:
                                break
                    if not E.isalive():
                        try:
                            self._synchronize()
                        except (TypeError, RuntimeError):
                            pass
                        return self._eval_line(line,allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False, **kwds)
                raise_(RuntimeError, "%s\nError evaluating %s in %s"%(msg, line, self), sys.exc_info()[2])

            p_warnings = []
            p_errors = []
            have_warning = False
            have_error = False
            have_log = False
            if len(line)>0:
                first = True
                while True:
                    try:
                        if isinstance(wait_for_prompt, six.string_types):
                            pat = E.expect(wait_for_prompt, **kwds)
                        else:
                            pat = E.expect_list(self._prompt, **kwds)
                    except pexpect.EOF as msg:
                        try:
                            if self.is_local():
                                tmp_to_use = self._local_tmpfile()
                            else:
                                tmp_to_use = self._remote_tmpfile()
                            if self._read_in_file_command(tmp_to_use) in line:
                                raise pexpect.EOF(msg)
                        except NotImplementedError:
                            pass
                        if self._quit_string() in line:
                            # we expect to get an EOF if we're quitting.
                            return ''
                        elif restart_if_needed: # the subprocess might have crashed
                            try:
                                self._synchronize()
                                return self._eval_line(line,allow_use_file=allow_use_file, wait_for_prompt=wait_for_prompt, restart_if_needed=False, **kwds)
                            except (TypeError, RuntimeError):
                                pass
                        raise RuntimeError("%s\n%s crashed executing %s"%(msg,self, line))
                    if self._terminal_echo:
                        out = E.before
                    else:
                        out = E.before.rstrip('\n\r')
                    if self._terminal_echo and first:
                        i = out.find("\n")
                        j = out.rfind("\r")
                        out = out[i+1:j].replace('\r\n','\n')
                    else:
                        out = out.strip().replace('\r\n','\n')
                    first = False
                    if have_error:
                        p_errors.append(out)
                        have_error = False
                        out = ""
                    elif have_warning:
                        p_warnings.append(out)
                        have_warning = False
                        out = ""
                    elif have_log:
                        if get_verbose() > 0:
                            print(out)
                        have_log = False
                        out = ""
                    # 0: normal prompt
                    # 1: continuation prompt
                    # 2: user input expected when requestion "help"
                    # 3: what we are looking for when interrupting a computation
                    # 4: error
                    # 5: warning
                    # 6: anything but an error or warning, thus, an information
                    # 7: unexpected end of the stream
                    # 8: (expected) timeout
                    if pat == 0:
                        have_log = False
                        have_error = False
                        have_warning = False
                        if E.buffer:
                            if not E.buffer.strip():
                                E.send(chr(3))
                                sleep(0.1)
                                pat = E.expect_list(self._prompt)
                                if E.buffer or pat:
                                    raise RuntimeError("Couldn't return to prompt after command '{}'".format(line))
                        break
                    elif pat == 1: # unexpected continuation prompt
                        # Return to normal prompt
                        i = pat
                        E.send(chr(3))
                        sleep(0.1)
                        i = E.expect_list(self._prompt)
                        assert i==0, "Command '{}': Couldn't return to normal prompt after polymake {}. Instead, polymake {}".format(line,_available_polymake_answers[pat],_available_polymake_answers[i])
                        raise SyntaxError("Incomplete polymake command '{}'".format(line))
                    elif pat == 2: # request for user interaction
                        # Return to normal prompt
                        warnings.warn("{} expects user interaction. We abort and return the options that {} provides.".format(self,self))
                        i = pat
                        while i:
                            self._expect.send(chr(3))
                            sleep(0.1)
                            i = self._expect.expect(self._prompt, timeout=0.1)
                        # User interaction is expected to happen when requesting help
                        if line.startswith('help'):
                            out = os.linesep.join(out.split(os.linesep)[:-1])
                            break
                        else:
                            RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))
                    elif pat == 3: # killed by signal
                        i = pat
                        while pat != 0:
                            E.send(chr(3))
                            sleep(0.1)
                            i = E.expect_list(self._prompt)
                        RuntimeError("Polymake unexpectedly {}".format(_available_polymake_answers[pat]))
                    elif pat == 4: # polymake error
                        have_error = True
                    elif pat == 5: # polymake warning
                        have_warning = True
                    elif pat == 6: # apparently polymake prints a comment
                        have_log = True
                    elif pat == 7: # we have reached the end of the buffer
                        warnings.warn("Polymake unexpectedly {}".format(_available_polymake_answers[pat]), RuntimeWarning)
                        E.buffer = E.before + E.after + E.buffer
                        break
                    else: # timeout or some other problem
                        # Polymake would still continue with the computation. Thus, we send an interrupt
                        E.send(chr(3))
                        sleep(0.1)
                        while E.expect_list(self._prompt, timeout=0.1):
                            # ... and since a single Ctrl-c just interrupts *one* of polymake's
                            # rule chains, we repeat until polymake is running out of rules.
                            E.send(chr(3))
                            sleep(0.1)
                        raise RuntimeError("Polymake {}".format(_available_polymake_answers[pat]))
            else:
                out = ''
        except KeyboardInterrupt:
            self._keyboard_interrupt()
            raise KeyboardInterrupt("Ctrl-c pressed while running %s"%self)
        for w in p_warnings:
            warnings.warn(w, RuntimeWarning)
        for e in p_errors:
            raise PolymakeError(e)
        return out

    def _tab_completion(self):
        """
        Returns a list of polymake function names.

        NOTE:

        - The list of functions depends on the current application. The
          result is cached, of course separately for each application.
        - It is generally not the case that all the returned function names
          can actually successfully be called.

        TESTS::

            sage: polymake.application('fan')                   # optional - polymake
            sage: 'normal_fan' in dir(polymake)                 # optional - polymake  # indirect doctest
            True
            sage: polymake.quit()                               # optional - polymake
            sage: polymake._start()                             # optional - polymake

        Since 'normal_fan' is not defined in the polymake application 'polytope',
        we now get
        ::

            sage: 'normal_fan' in dir(polymake)                 # optional - polymake
            False

        """
        if not self.is_running():
            self._start()
        try:
            return self.__tab_completion[self._application]
        except KeyError:
            pass
        s = self.eval("apropos '';").split(self._expect.linesep)
        out = []
        for name in s:
            if name.startswith("/function"):
                out.append(name.split("/")[-1])
        self.__tab_completion[self._application] = sorted(out)
        return self.__tab_completion[self._application]

    # Polymake specific methods

    def application(self, app):
        """
        Change to a given polymake application.

        INPUT:

        - ``app``, a string, one of "common", "fulton", "group", "matroid", "topaz",
          "fan", "graph", "ideal", "polytope", "tropical"

        EXAMPLES:

        We expose a computation that uses both the 'polytope' and the 'fan'
        application of polymake. Let us start by defining a polytope `q` in
        terms of inequalities. Polymake knows to compute the f- and h-vector
        and finds that the polytope is very ample::

            sage: q = polymake.new_object("Polytope", INEQUALITIES=[[5,-4,0,1],[-3,0,-4,1],[-2,1,0,0],[-4,4,4,-1],[0,0,1,0],[8,0,0,-1],[1,0,-1,0],[3,-1,0,0]]) # optional - polymake
            sage: q.H_VECTOR                    # optional - polymake
            1 5 5 1
            sage: q.F_VECTOR                    # optional - polymake
            8 14 8
            sage: q.VERY_AMPLE                  # optional - polymake
            1

        In the application 'fan', polymake can now compute the normal fan
        of `q` and its (primitive) rays::

            sage: polymake.application('fan')   # optional - polymake
            sage: g = q.normal_fan()            # optional - polymake
            sage: g.RAYS                        # optional - polymake
            -1 0 1/4
            0 -1 1/4
            1 0 0
            1 1 -1/4
            0 1 0
            0 0 -1
            0 -1 0
            -1 0 0
            sage: g.RAYS.primitive()            # optional - polymake
            -4 0 1
            0 -4 1
            1 0 0
            4 4 -1
            0 1 0
            0 0 -1
            0 -1 0
            -1 0 0

        Note that the list of functions available by tab completion depends
        on the application.

        TESTS:

        Since 'tubing_of_graph' is not defined in the polymake application 'polytope'
        but only in 'tropical', the following shows the effect of changing
        the application. ::

            sage: polymake.application('polytope')                   # optional - polymake
            sage: 'tubing_of_graph' in dir(polymake)                 # optional - polymake
            False
            sage: polymake.application('tropical')                   # optional - polymake
            sage: 'tubing_of_graph' in dir(polymake)                 # optional - polymake
            True
            sage: polymake.application('polytope')                   # optional - polymake
            sage: 'tubing_of_graph' in dir(polymake)                 # optional - polymake
            False

        For completeness, we show what happens when asking for an application
        that doesn't exist::

            sage: polymake.application('killerapp')                  # optional - polymake
            Traceback (most recent call last):
            ...
            ValueError: Unknown polymake application 'killerapp'

        Of course, a different error results when we send an explicit
        command in polymake to change to an unknown application::

            sage: polymake.eval('application "killerapp";')         # optional - polymake
            Traceback (most recent call last):
            ...
            PolymakeError: Unknown application killerapp

        """
        if not self.is_running():
            self._start()
        if app not in ["common", "fulton", "group", "matroid", "topaz", "fan", "graph", "ideal", "polytope", "tropical"]:
            raise ValueError("Unknown polymake application '{}'".format(app))
        self._application = app
        patterns = ["{} > ".format(app),            # 0: normal prompt
                    "{} \([0-9]+\)> ".format(app),  # 1: continuation prompt
                    "Please choose ".format(app),   # 2: user input expected when requesting "help"
                    "killed by signal",             # 3: what we are looking for when interrupting a computation
                    "polymake: +ERROR: +",          # 4: error
                    "polymake: +WARNING: +",        # 5: warning
                    "polymake: +",                  # 6: anything but an error or warning, thus, an information
                    pexpect.EOF,                    # 7: unexpected end of the stream
                    pexpect.TIMEOUT]                # 8: timeout
        self._change_prompt(self._expect.compile_pattern_list(patterns))
        self._sendstr('application "{}";{}'.format(app, self._expect.linesep))
        pat = self._expect.expect_list(self._prompt)
        if pat:
            raise RuntimeError("When changing the application, polymake unexpectedly {}".format(_available_polymake_answers[pat]))

    def new_object(self, name, *args, **kwds):
        """
        Return a new instance of a given polymake type, with given positional or named arguments.

        INPUT:

        - ``name`` of a polymake class (potentially templated), as string.
        - further positional or named arguments, to be passed to the constructor.

        EXAMPLES::

            sage: q = polymake.new_object("Polytope<Rational>", INEQUALITIES=[[4,-4,0,1],[-4,0,-4,1],[-2,1,0,0],[-4,4,4,-1],[0,0,1,0],[8,0,0,-1]]) # optional - polymake
            sage: q.N_VERTICES                  # optional - polymake
            4
            sage: q.BOUNDED                     # optional - polymake
            1
            sage: q.VERTICES                    # optional - polymake
            1 2 0 4
            1 3 0 8
            1 2 1 8
            1 3 1 8
            sage: q.full_typename()             # optional - polymake
            'Polytope<Rational>'

        """
        try:
            f = self.__new[name]
        except AttributeError:
            self.__new = {}
            f = self.__new[name] = self._function_class()(self, "new %s"%name)
        except KeyError:
            f = self.__new[name] = self._function_class()(self, "new %s"%name)
        return f(*args, **kwds)

polymake = Polymake()

def reduce_load_Polymake():
    """
    Returns the polymake interface object defined in :mod:`sage.interfaces.polymake`.

    EXAMPLES::

        sage: from sage.interfaces.polymake import reduce_load_Polymake
        sage: reduce_load_Polymake()
        Polymake
    """
    return polymake

########################################
## Elements

from warnings import warn

class PolymakeElement(ExtraTabCompletion, ExpectElement):
    """
    Elements in the polymake interface.

    EXAMPLES:

    We support all "big" polymake types, Perl arrays of length
    different from one, and Perl scalars::

        sage: p = polymake.rand_sphere(4, 20, seed=5)            # optional - polymake
        sage: p.typename()                                      # optional - polymake
        'Polytope'
        sage: p                                                 # optional - polymake
        Random spherical polytope of dimension 4; seed=5...

    Now, one can work with that element in Python syntax, for example::

        sage: p.VERTICES[2][2]                                  # optional - polymake
        -3319173990813887/4503599627370496

    """
    def _repr_(self):
        """
        String representation of polymake elements.

        EXAMPLES:

        In the case of a "big" object, if polymake provides a description
        of the object that is not longer than single line, it is used for
        printing::

            sage: p = polymake.rand_sphere(3, 12, seed=15)           # optional - polymake
            sage: p                                                 # optional - polymake
            Random spherical polytope of dimension 3; seed=15...
            sage: c = polymake.cube(4)                              # optional - polymake
            sage: c                                                 # optional - polymake
            cube of dimension 4

        We use the print representation of scalars to display scalars::

            sage: p.N_VERTICES                                      # optional - polymake
            12

        The items of a Perl arrays are shown separated by commas::

            sage: p.get_member('list_properties')                   # optional - polymake  # random
            POINTS, CONE_AMBIENT_DIM, BOUNDED, FEASIBLE, N_POINTS, POINTED,
            CONE_DIM, FULL_DIM, LINEALITY_DIM, LINEALITY_SPACE,
            COMBINATORIAL_DIM, AFFINE_HULL, VERTICES, N_VERTICES

        We chose to print rule chains explicitly, so that the user doesn't
        need to know how to list the rules using polymake commands::

            sage: r = p.get_schedule('"H_VECTOR"')                  # optional - polymake
            sage: r                                                 # optional - polymake  # random
            precondition : N_RAYS | N_INPUT_RAYS ( ppl.convex_hull.primal: FACETS, LINEAR_SPAN : RAYS | INPUT_RAYS )
            sensitivity check for FacetPerm
            ppl.convex_hull.primal: FACETS, LINEAR_SPAN : RAYS | INPUT_RAYS
            RAYS_IN_FACETS : RAYS, FACETS
            SIMPLICIAL : COMBINATORIAL_DIM, RAYS_IN_FACETS
            N_FACETS : FACETS
            precondition : COMBINATORIAL_DIM ( F_VECTOR : N_FACETS, N_RAYS, COMBINATORIAL_DIM )
            F_VECTOR : N_FACETS, N_RAYS, COMBINATORIAL_DIM
            precondition : SIMPLICIAL ( H_VECTOR : F_VECTOR )
            H_VECTOR : F_VECTOR
            sage: r.typeof()                                        # optional - polymake
            ('Polymake::Core::Scheduler::RuleChain', 'ARRAY')

        Similarly, polymake matrices and vectors are explicitly listed::

            sage: c.VERTICES.typename()                         # optional - polymake
            'Matrix'
            sage: c.VERTICES[0].typename()                      # optional - polymake
            'Vector'
            sage: c.VERTICES                                    # optional - polymake # random
            1 -1 -1 -1 -1
            1 1 -1 -1 -1
            1 -1 1 -1 -1
            1 1 1 -1 -1
            1 -1 -1 1 -1
            1 1 -1 1 -1
            1 -1 1 1 -1
            1 1 1 1 -1
            1 -1 -1 -1 1
            1 1 -1 -1 1
            1 -1 1 -1 1
            1 1 1 -1 1
            1 -1 -1 1 1
            1 1 -1 1 1
            1 -1 1 1 1
            1 1 1 1 1
            sage: c.VERTICES[0]                                 # optional - polymake
            1 -1 -1 -1 -1

        For other types, we simply use the print representation offered
        by polymake::

            sage: p.TWO_FACE_SIZES.typename()                   # optional - polymake
            'Map'
            sage: p.TWO_FACE_SIZES                              # optional - polymake
            {(3 20)}

    """
        T1, T2 = self.typeof()
        P = self._check_valid()
        name = self._name
        if T1:
            Temp = self.typename()
            if Temp:
                T1 = Temp
        if T1 in ['Matrix', 'Vector']:
            out = P.get(name).strip()
        elif 'RuleChain' in T1:
            out = os.linesep.join(P.get('join("##",{}->list)'.format(name)).split('##'))
        else:
            try:
                out = P.get('{}->description'.format(name)).strip()
            except PolymakeError:
                out = ''
            if os.linesep in out:
                out = ''
        if not out:
            if "Polytope" == T1:
                out = "{}[{}]".format(P.get("{}->type->full_name".format(name)) or "PolymakeElement", _name_pattern.search(name).group())
            elif T1=='' and T2=='ARRAY':
                out = P.eval('print join(", ", @{});'.format(name)).strip()
            elif T1=='' and T2=='HASH':
                out = P.get('%{}'.format(name)).strip()
            elif self._name[0] == '@':
                out = P.eval('print join(", ", {});'.format(name)).strip()
            else:
                out = P.get(name).strip()
        return out

    def __cmp__(self, other):
        """
        Comparison of polymake elements.

        EXAMPLES:

        The default for comparing equality for polytopes is *identity*::

            sage: p1 = polymake.rand_sphere(3, 12, seed=15)          # optional - polymake
            sage: p2 = polymake.rand_sphere(3, 12, seed=15)          # optional - polymake
            sage: p1 == p2                                          # optional - polymake
            False

        However, other data types are compared by equality, not identity::

            sage: p1.VERTICES == p2.VERTICES                        # optional - polymake
            True

        A computation applied to a polytope can change the available
        properties, and thus we have
        ::

            sage: p1.get_member('list_properties') == p2.get_member('list_properties')  # optional - polymake
            True
            sage: p1.F_VECTOR                                       # optional - polymake
            12 30 20
            sage: p1.get_member('list_properties') == p2.get_member('list_properties')  # optional - polymake
            False

        """
        P = self._check_valid()
        if P.eval("print %s %s %s;"%(self.name(), P._equality_symbol(), other.name())).strip() == P._true_symbol():
            return 0
        if P.eval("print %s %s %s;"%(self.name(), P._lessthan_symbol(), other.name())).strip() == P._true_symbol():
            return -1
        if P.eval("print %s %s %s;"%(self.name(), P._greaterthan_symbol(), other.name())).strip() == P._true_symbol():
            return 1
        return -2 # that's supposed to be an error value.

    def bool(self):
        """
        Return whether this polymake element is equal to ``True``.

        EXAMPLES::

            sage: from sage.interfaces.polymake import polymake
            sage: polymake(0).bool()                # optional polymake
            False
            sage: polymake(1).bool()                # optional polymake
            True

        """
        P = self._check_valid()
        t = P._true_symbol()
        cmd = '%s %s %s;'%(self._name, P._equality_symbol(), t)
        return P.get(cmd) == t

    def known_properties(self):
        """
        List the names of properties that have been computed so far on this element.

        NOTE:

        This is in many cases equivalent to use polymake's ``list_properties``,
        which returns a blank separated string representation of the list of properties.
        However, on some elements, ``list_properties`` would simply result in
        an error.

        EXAMPLES::

            sage: c = polymake.cube(4)                      # optional - polymake
            sage: c.known_properties()                      # optional - polymake
            ['AFFINE_HULL',
             'BOUNDED',
             'CONE_AMBIENT_DIM',
             'CONE_DIM',
            ...
             'VERTICES_IN_FACETS']
            sage: c.list_properties()                       # optional - polymake
            CONE_AMBIENT_DIM, CONE_DIM, FACETS, AFFINE_HULL, VERTICES_IN_FACETS,
            BOUNDED...

        A computation can change the list of known properties::

            sage: c.F_VECTOR                                # optional - polymake
            16 32 24 8
            sage: c.known_properties()                      # optional - polymake
            ['AFFINE_HULL',
             'BOUNDED',
             'COMBINATORIAL_DIM',
             'CONE_AMBIENT_DIM',
             'CONE_DIM',
             'DUAL_H_VECTOR',
             'FACETS',
             'FAR_FACE',
             'FEASIBLE',
             'FULL_DIM',
             'F_VECTOR',
             'GRAPH',
             'LINEALITY_DIM',
             'LINEALITY_SPACE',
             'N_FACETS',
             'N_VERTICES',
             'POINTED',
             'SIMPLE',
             'SIMPLICIAL',
             'VERTICES',
             'VERTICES_IN_FACETS']

        """
        P = self._check_valid()
        try:
            return sorted(P.get('join(", ", {}->list_properties)'.format(self._name)).split(', '))
        except PolymakeError:
            return []

    @cached_method
    def _member_list(self):
        """
        The list of properties that polymake knows to compute for this element.

        The resulting list is used for tab completion.

        TESTS::

            sage: c = polymake.cube(4)                          # optional - polymake
            sage: c._member_list()                              # optional - polymake
            ['AFFINE_HULL',
             'ALTSHULER_DET',
             'BALANCE',
             'BALANCED',
            ...
             'WEAKLY_CENTERED',
             'ZONOTOPE_INPUT_POINTS']

        """
        ### return the members of a "big" object.
        P = self._check_valid()
        try:
            P.eval('$SAGETMP = typeof {+'+self._name+'};')
        except (TypeError, PolymakeError):  # this happens for a perl type that isn't a Polymake type
            return []
        cmd = 'print join(", ", sorted_uniq(sort { $a cmp $b } map { keys %{$_->properties} }$SAGETMP, @{$SAGETMP->super}));'
        try:
            out = P.eval(cmd).split(', ')
        except PolymakeError as msg:
            return []
        return sorted(out)

    def typename(self):
        """
        The name of the underlying base type of this element in polymake.

        EXAMPLES::

            sage: c = polymake.cube(4)          # optional - polymake
            sage: c.typename()                  # optional - polymake
            'Polytope'
            sage: c.VERTICES.typename()         # optional - polymake
            'Matrix'

        """
        P = self._check_valid()
        try:
            return P.eval('print {}->type->name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def full_typename(self):
        """
        The name of the specialised type of this element.

        EXAMPLES::

            sage: c = polymake.cube(4)          # optional - polymake
            sage: c.full_typename()             # optional - polymake
            'Polytope<Rational>'
            sage: c.VERTICES.full_typename()    # optional - polymake
            'Matrix<Rational, NonSymmetric>'

        """
        P = self._check_valid()
        try:
            return P.eval('print {}->type->full_name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def qualified_typename(self):
        """
        The qualified name of the type of this element.

        EXAMPLES::

            sage: c = polymake.cube(4)              # optional - polymake
            sage: c.qualified_typename()            # optional - polymake
            'polytope::Polytope<Rational>'
            sage: c.VERTICES.qualified_typename()   # optional - polymake
            'common::Matrix<Rational, NonSymmetric>'

        """
        P = self._check_valid()
        try:
            return P.eval('print {}->type->qualified_name;'.format(self._name))
        except PolymakeError as msg:
            return ''

    def _tab_completion(self):
        """
        Return a list of available function and property names.

        NOTE:

        This currently returns the names of functions defined in the current
        application, regardless whether they can be applied to this element
        or not, together with the list of properties of this element that
        polymake knows how to compute. It does not contain the list of available
        member functions of this element. This may change in future versions
        of polymake.

        EXAMPLES::

            sage: c = polymake.cube(4)              # optional - polymake
            sage: c._tab_completion()               # optional - polymake
            ['AFFINE_HULL',
             'ALTSHULER_DET',
             'BALANCE',
             'BALANCED',
             'BOUNDARY_LATTICE_POINTS',
             ...
             'zero_vector',
             'zonotope',
             'zonotope_tiling_lattice',
             'zonotope_vertices_fukuda']

        """
        return sorted(self._member_list()+self.parent()._tab_completion())

    def __getattr__(self, attrname):
        """
        Return a property of this element, or a polymake function with this
        element as first argument, or a member function of this element.

        NOTE:

        If the attribute name is known as the name of a property, it is
        interpreted as such. Otherwise, if it is known as a function in
        the currenty application, the function is returned with this
        element inserted as first argument, and potential further arguments,
        when called. Otherwise, it is assumed that it is a member function
        of this element, and treated as such. Note that member functions
        are currently invisible in tab completion, thus, the user has
        to know the name of the member function.

        EXAMPLES:

        A property::

            sage: c = polymake.cube(3)                  # optional - polymake
            sage: c.H_VECTOR                            # optional - polymake
            1 5 5 1
            sage: c.N_VERTICES                          # optional - polymake
            8
            sage: d = polymake.cross(3)                 # optional - polymake
            sage: d.N_VERTICES                          # optional - polymake
            6

        A function::

            sage: c.minkowski_sum_fukuda                # optional - polymake
            minkowski_sum_fukuda (bound to Polymake::polytope::Polytope__Rational object)
            sage: s = c.minkowski_sum_fukuda(d)         # optional - polymake
            sage: s.N_VERTICES                          # optional - polymake
            24
            sage: s                                     # optional - polymake
            Polytope<Rational>[SAGE...]

        A member function::

            sage: c = polymake.cube(2)                          # optional - polymake
            sage: V = polymake.new_object('Vector', [1,0,0])    # optional - polymake
            sage: V                                             # optional - polymake
            1 0 0
            sage: c.contains                                    # optional - polymake
            Member function 'contains' of Polymake::polytope::Polytope__Rational object
            sage: c.contains(V)                                 # optional - polymake
            1

        """
        P = self._check_valid()
        if attrname[:1] == "_":
            raise AttributeError
        if attrname not in P._tab_completion():
            if attrname in self._member_list():
                try:
                    return P('{}->{}'.format(self._name, attrname))
                except (TypeError, PolymakeError):
                    raise AttributeError
            else:
                try:
                    return P._function_element_class()(self, '{}->{}'.format(self._name, attrname), memberfunction=True)
                except (TypeError, PolymakeError):
                    raise AttributeError
        return P._function_element_class()(self, attrname, memberfunction=False)

    def get_member_function(self, attrname):
        """
        Request a member function of this element.

        NOTE:

        It is not checked whether a member function with the given name
        exists.

        EXAMPLES::

            sage: c = polymake.cube(2)                                  # optional - polymake
            sage: c.contains                                            # optional - polymake
            Member function 'contains' of Polymake::polytope::Polytope__Rational object
            sage: V = polymake.new_object('Vector', [1,0,0])            # optional - polymake
            sage: V                                                     # optional - polymake
            1 0 0
            sage: c.contains(V)                                         # optional - polymake
            1

        Whether a member function of the given name actually exists for that
        object will only be clear when calling it::

            sage: c.get_member_function('foo')                          # optional - polymake
            Member function 'foo' of Polymake::polytope::Polytope__Rational object
            sage: c.get_member_function('foo')()                        # optional - polymake
            Traceback (most recent call last):
            ...
            TypeError: Can't locate object method "foo" via package "Polymake::polytope::Polytope__Rational" at input line 1.

        """
        P = self._check_valid()
        return P._function_element_class()(self, '{}->{}'.format(self._name, attrname), memberfunction=True)

    def get_member(self, attrname):
        """
        Get a member/property of this element.

        NOTE:

        Normally, it should be possible to just access the property
        in the usual Python syntax for attribute access. However, if
        that fails, one can request the member explicitly.

        EXAMPLES::

            sage: p = polymake.rand_sphere(4, 20, seed=5)    # optional - polymake

        Normally, a property would be accessed as follows::

            sage: p.F_VECTOR                                # optional - polymake
            20 101 162 81

        However, explicit access is possible as well::

            sage: p.get_member('F_VECTOR')                  # optional - polymake
            20 101 162 81

        In some cases, the explicit access works better::

            sage: p.type                                    # optional - polymake
            Member function 'type' of Polymake::polytope::Polytope__Rational object
            sage: p.get_member('type')                      # optional - polymake
            Polytope<Rational>[SAGE...]
            sage: p.get_member('type').get_member('name')   # optional - polymake
            Polytope

        Note that in the last example calling the erroneously constructed
        member function ``type`` still works::

            sage: p.type()                                  # optional - polymake
            Polytope<Rational>[SAGE...]

        """
        P = self._check_valid()
        return P('%s->%s'%(self.name(), attrname))

    def __getitem__(self, key):
        """
        Indexing and slicing.

        Slicing returns a Python list.

        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 12, seed=15)  # optional - polymake
            sage: p.VERTICES[3]                             # optional - polymake
            1 -6157731020575175/18014398509481984 4184896164481703/4503599627370496 -2527292586301447/18014398509481984
            sage: p.list_properties()[2]                    # optional - polymake
            BOUNDED

        Slicing::

            sage: p.F_VECTOR[:]                             # optional - polymake
            [12, 30, 20]
            sage: p.F_VECTOR[0:1]                           # optional - polymake
            [12]
            sage: p.F_VECTOR[0:3:2]                         # optional - polymake
            [12, 20]
        """
        P = self._check_valid()
        if isinstance(key, slice):
            indices = key.indices(len(self))
            return [ self[i] for i in range(*indices) ]
        _, T = self.typeof()
        if self._name.startswith('@'):
            return P('${}[{}]'.format(self._name[1:], key))
        if T=='ARRAY':
            return P('{}[{}]'.format(self._name, key))
        if T=='HASH':
            try:
                if key.parent() is self.parent():
                    key = key._name
                else:
                    key = str(key)
            except AttributeError:
                key = str(key)
            return P(name+"{"+key+"}")
        raise NotImplementedError("Cannot get items from Perl type {}".format(T))

    def __iter__(self):
        """
        Return an iterator for ``self``.

        OUTPUT: iterator

        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 12, seed=15)  # optional - polymake
            sage: [ x for x in p.VERTICES[3] ]              # optional - polymake
            [1, -6157731020575175/18014398509481984, 4184896164481703/4503599627370496, -2527292586301447/18014398509481984]
        """
        for i in range(len(self)):
            yield self[i]

    def __len__(self):
        """
        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 12, seed=15)           # optional - polymake
            sage: len(p.FACETS)                                     # optional - polymake
            20
            sage: len(p.list_properties())                          # optional - polymake
            13

        """
        P = self._check_valid()
        T1,T2 = self.typeof()
        name = self._name
        if T2=='ARRAY':
            return int(P.eval('print scalar @{+%s};'%name))
        if T2=='HASH':
            return int(P.eval('print scalar keys %{+%s};'%name))
        if T1:
            raise TypeError("Don't know how to compute the length of {} object".format(T1))
        return int(P.eval('print scalar {};'.format(name)))

    @cached_method
    def typeof(self):
        """
        Returns the type of a polymake "big" object, and its underlying Perl type.

        NOTE:

        This is mainly for internal use.

        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 13, seed=12)               # optional - polymake
            sage: p.typeof()                                            # optional - polymake
            ('Polymake::polytope::Polytope__Rational', 'ARRAY')
            sage: p.VERTICES.typeof()                                   # optional - polymake
            ('Polymake::common::Matrix_A_Rational_I_NonSymmetric_Z', 'ARRAY')
            sage: p.get_schedule("F_VECTOR").typeof()                   # optional - polymake
            ('Polymake::Core::Scheduler::RuleChain', 'ARRAY')

        On "small" objects, it just returns empty strings::

            sage: p.N_VERTICES.typeof()                                 # optional - polymake
            ('', '')
            sage: p.list_properties().typeof()                          # optional - polymake
            ('', '')
        """
        P = self._check_valid()
        name = self._name
        return P.eval('print ref({});'.format(name)), P.eval('print reftype({});'.format(name))

    def _sage_(self):
        """
        Convert self to a Sage object.

        EXAMPLES::

            sage: a = polymake(1/2); a    # optional - polymake
            1/2
            sage: a.sage()                # optional - polymake
            1/2
            sage: _.parent()              # optional - polymake
            Rational Field

        Quadratic extensions::

            sage: K.<sqrt5> = QuadraticField(5)
            sage: polymake(K(0)).sage()   # optional - polymake
            0
            sage: _.parent()              # optional - polymake
            Rational Field
            sage: polymake(sqrt5).sage()   # optional - polymake
            a
            sage: polymake(-sqrt5).sage()   # optional - polymake
            -a
            sage: polymake(1/3-1/2*sqrt5).sage()   # optional - polymake
            -1/2*a + 1/3
            sage: polymake(-1+sqrt5).sage()   # optional - polymake
            a - 1

        Vectors::

            sage: PP = polymake.cube(3)   # optional - polymake
            sage: PP.F_VECTOR.sage()      # optional - polymake
            (8, 12, 6)
            sage: _.parent()              # optional - polymake
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        Matrices::

            sage: polymake.unit_matrix(2).sage()   # optional - polymake
            [1 0]
            [0 1]
            sage: _.parent()              # optional - polymake
            Full MatrixSpace of 2 by 2 dense matrices over Integer Ring

        """
        T1, T2 = self.typeof()
        P = self._check_valid()
        if T1:
            Temp = self.typename()
            if Temp:
                T1 = Temp
        if T1 == 'QuadraticExtension':
            # We can't seem to access a, b, r by method calls, so let's parse.
            from re import match
            m = match(r'(-?[0-9/]+)[+]?((-?[0-9/]+)r([0-9/]+))?', repr(self))
            if m is None:
                raise NotImplementedError("Cannot parse QuadraticExtension element: {}".format(self))
            a, b, r = m.group(1), m.group(3), m.group(4)
            from sage.rings.rational_field import QQ
            if r is None:
                # Prints like a rational, so we can't know the extension. Coerce to rational.
                return QQ(a)
            else:
                from sage.rings.number_field.number_field import QuadraticField
                K = QuadraticField(r)
                return QQ(a) + QQ(b) * K.gen()
        elif T1 == 'Vector' or T1 == 'SparseVector':
            from sage.modules.free_module_element import vector
            return vector([x.sage() for x in self])
        elif T1 == 'Matrix' or T1 == 'SparseMatrix':
            from sage.matrix.constructor import matrix
            return matrix([x.sage() for x in self])
        else:
            return super(PolymakeElement, self)._sage_()

    def _sage_doc_(self):
        """
        EXAMPLES::

            sage: c = polymake.cube(3)                  # optional - polymake
            sage: print(c._sage_doc_())                 # optional - polymake # random
            objects/Polytope:
             Not necessarily bounded or unbounded polyhedron.
             Nonetheless, the name "Polytope" is used for two reasons:
             Firstly, combinatorially we always deal with polytopes; see the description of VERTICES_IN_FACETS for details.
             The second reason is historical.
             We use homogeneous coordinates, which is why Polytope is derived from Cone.
             Note that a pointed polyhedron is projectively equivalent to a polytope.
             Scalar is the numeric data type used for the coordinates.
            <BLANKLINE>
            objects/Polytope/specializations/Polytope<Rational>:
             A rational polyhedron realized in Q^d
            sage: print(c.FACETS._sage_doc_())          # optional - polymake # random
            property_types/Algebraic Types/SparseMatrix:
             A SparseMatrix is a two-dimensional associative array with row and column indices as keys; elements equal to the default value (ElementType(), which is 0 for most numerical types) are not stored, but implicitly encoded by the gaps in the key set. Each row and column is organized as an AVL-tree.
            <BLANKLINE>
             Use dense to convert this into its dense form.
            <BLANKLINE>
             You can create a new SparseMatrix by entering its entries row by row, as a list of SparseVectors e.g.:
                $A = new SparseMatrix<Int>(<< '.');
                (5) (1 1)
                (5) (4 2)
                (5)
                (5) (0 3) (1 -1)
                .

        """
        P = self._check_valid()
        # according to Julian Pfeifle, the only case in which the fully qualified
        # typename would not provide the doc.
        Tname = self.typename()
        Tqname = self.qualified_typename()
        Tfname = self.full_typename()
        if Tname == 'Polytope':
            try:
                doc = P.eval('help "Polytope";')
            except PolymakeError:
                doc = ''
        else:
            try:
                doc = P.eval('help "{}";'.format(Tname))
            except PolymakeError:
                doc = ''
            try:
                doc2 = P.eval('help "{}";'.format(Tqname))
            except PolymakeError:
                doc2 = ''
            if doc:
                if doc2:
                    doc = doc+os.linesep+doc2
            else:
                doc = doc2
        try:
            doc3 = P.eval('help "{}";'.format(Tfname))
        except PolymakeError:
            doc3 = ''
        if doc:
            if doc3:
                doc = doc+os.linesep+doc3
        else:
            doc = doc3
        if doc:
            return doc
        return "Undocumented polymake type '{}'".format(self.full_typename())

class PolymakeFunctionElement(FunctionElement):
    """
    A callable (function or member function) bound to a polymake element.

    EXAMPLES::

        sage: c = polymake.cube(2)                          # optional - polymake
        sage: V = polymake.new_object('Vector', [1,0,0])    # optional - polymake
        sage: V                                             # optional - polymake
        1 0 0
        sage: c.contains                                    # optional - polymake
        Member function 'contains' of Polymake::polytope::Polytope__Rational object
        sage: c.contains(V)                                 # optional - polymake
        1

    """
    def __init__(self, obj, name, memberfunction=False):
        """
        INPUT:

        - Polymake object that this function is bound to
        - name (string): It actually says how to call this function in polymake.
          So, if it is a member function, it will look like `"$SAGE123[0]->func_name"`.
        - ``memberfunction`` (bool, default False): Whether this is a member function
          or a plain function applied with this element as first argument.

        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 13, seed=12)   # optional - polymake
            sage: p.minkowski_sum_fukuda                    # optional - polymake
            minkowski_sum_fukuda (bound to Polymake::polytope::Polytope__Rational object)
            sage: p.get_schedule                            # optional - polymake
            Member function 'get_schedule' of Polymake::polytope::Polytope__Rational object

        """
        self._obj = obj
        self._name = name
        self._is_memberfunc = memberfunction
    def _repr_(self):
        """
        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 13, seed=12)  # optional - polymake
            sage: p.minkowski_sum_fukuda                    # optional - polymake
            minkowski_sum_fukuda (bound to Polymake::polytope::Polytope__Rational object)
            sage: p.contains                                # optional - polymake
            Member function 'contains' of Polymake::polytope::Polytope__Rational object

        """
        if self._is_memberfunc:
            return "Member function '{}' of {} object".format(self._name.split("->")[-1], self._obj.typeof()[0])
        return "{} (bound to {} object)".format(self._name, self._obj.typeof()[0])

    def __call__(self, *args, **kwds):
        """
        EXAMPLES:

        We consider both member functions of an element and global functions
        bound to an element::

            sage: p = polymake.rand_sphere(3, 13, seed=12)      # optional - polymake
            sage: p.get_schedule('VERTICES')                    # optional - polymake  # random
            sensitivity check for VertexPerm
            cdd.convex_hull.canon: POINTED, RAYS, LINEALITY_SPACE : INPUT_RAYS
            sage: p.minkowski_sum_fukuda(p).F_VECTOR            # optional - polymake
            13 33 22

        """
        if self._is_memberfunc:
            return self._obj._check_valid().function_call(self._name, list(args), kwds)
        return self._obj._check_valid().function_call(self._name, [self._obj] + list(args), kwds)

    def _sage_doc_(self):
        """
        Return documentation of this function.

        NOTE:

        For unclear reasons, accessing documentation with `?` sometimes
        does not include the return value of this method.

        EXAMPLES::

            sage: p = polymake.rand_sphere(3, 13, seed=12)           # optional - polymake
            sage: print(p.get_schedule._sage_doc_())                 # optional - polymake # random
            objects/Core::Object/methods/get_schedule:
            get_schedule(request;  ... ) -> Core::RuleChain
            <BLANKLINE>
             Compose an optimal chain of production rules providing all requested properties.
             The returned RuleChain object can be applied to the original object as well as to any other object
             with the same initial set of properties.  If no feasible rule chain exists, `undef' is returned.
            <BLANKLINE>
             To watch the rule scheduler at work, e.g. to see announcements about tried preconditions,
             you may temporarily increase the verbosity levels $Verbose::rules and $Verbose::scheduler.
            <BLANKLINE>
            Arguments:
              String request : name of a property with optional alternatives or a property path in dotted notation.
                Several requests may be listed.
            <BLANKLINE>
            Returns Core::RuleChain
            sage: print(p.minkowski_sum_fukuda._sage_doc_())        # optional - polymake # random
            functions/Producing a polytope from polytopes/minkowski_sum_fukuda:
            minkowski_sum_fukuda(summands) -> Polytope<Scalar>
            <BLANKLINE>
             Computes the (VERTICES of the) Minkowski sum of a list of polytopes using the algorithm by Fukuda described in
                   Komei Fukuda, From the zonotope construction to the Minkowski addition of convex polytopes, J. Symbolic Comput., 38(4):1261-1272, 2004.
            <BLANKLINE>
            Arguments:
              Array<Polytope<Scalar>> summands
            <BLANKLINE>
            Returns Polytope<Scalar>
            <BLANKLINE>
            Example:
                > $p = minkowski_sum_fukuda([cube(2),simplex(2),cross(2)]);
                > print $p->VERTICES;
                1 -2 -1
                1 -1 -2
                1 3 -1
                1 3 1
                1 2 -2
                1 -2 2
                1 -1 3
                1 1 3

        """
        P = self._obj._check_valid()
        return P.help(self._name.split("->")[-1], pager=False)
