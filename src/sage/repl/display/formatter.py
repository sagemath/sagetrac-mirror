# -*- coding: utf-8 -*-
r"""
IPython Displayhook Formatters

The classes in this module can be used as IPython displayhook
formatters. It has two main features, by default the displayhook
contains a new facility for displaying lists of matrices in an easier
to read format::

    sage: [identity_matrix(i) for i in range(2,5)]
    [
                    [1 0 0 0]
           [1 0 0]  [0 1 0 0]
    [1 0]  [0 1 0]  [0 0 1 0]
    [0 1], [0 0 1], [0 0 0 1]
    ]

This facility uses :meth:`_repr_` (and a simple string) to try do a nice read
format (see :meth:`sage.structure.parent.Parent._repr_option` for details).

With this displayhook there exists an other way for displaying object and more
generally, all sage expression as an ASCII art object::

    sage: from sage.repl.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: shell.run_cell('%display ascii_art')
    sage: shell.run_cell('integral(x^2/pi^x, x)')
     / 2    2                      \  -x*log(pi)
    -\x *log (pi) + 2*x*log(pi) + 2/*e
    ---------------------------------------------
                         3
                      log (pi)
    sage: shell.run_cell("i = var('i')")
    sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
        10      9      8      7      6      5      4      3      2
    10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
    sage: shell.run_cell('StandardTableaux(4).list()')
    [
    [                                                                  1  4
    [                 1  3  4    1  2  4    1  2  3    1  3    1  2    2
    [   1  2  3  4,   2      ,   3      ,   4      ,   2  4,   3  4,   3   ,
    <BLANKLINE>
                       1 ]
       1  3    1  2    2 ]
       2       3       3 ]
       4   ,   4   ,   4 ]
    sage: shell.run_cell('%display default')
    sage: shell.quit()

This other facility uses a simple
:class:`~sage.typeset.ascii_art.AsciiArt` object (see and
:meth:`sage.structure.sage_object.SageObject._ascii_art_`).  """

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from IPython.core.formatters import DisplayFormatter, PlainTextFormatter
from IPython.utils.py3compat import str_to_unicode, unicode_to_str

from sage.structure.sage_object import SageObject
from sage.repl.display.pretty_print import SagePrettyPrinter



class SageDisplayFormatter(DisplayFormatter):

    def __init__(self, *args, **kwds):
        """
        This is where the Sage rich objects are translated to IPython

        INPUT/OUTPUT:

        See the IPython documentation.

        EXAMPLES:

        This is part of how Sage works with the IPython output
        system. It cannot be used in doctests::

            sage: from sage.repl.display.formatter import SageDisplayFormatter
            sage: fmt = SageDisplayFormatter()
            Traceback (most recent call last):
            ...
            RuntimeError: check failed: current backend is invalid
        """
        super(SageDisplayFormatter, self).__init__(*args, **kwds)
        from sage.repl.rich_output.display_manager import get_display_manager
        self.dm = get_display_manager()
        from sage.repl.rich_output.backend_ipython import BackendIPython
        self.dm.check_backend_class(BackendIPython)

    def format(self, obj, include=None, exclude=None):
        r"""
        Use the Sage rich output instead of IPython

        INPUT/OUTPUT:

        See the IPython documentation.

        EXAMPLES::

            sage: [identity_matrix(i) for i in range(3,7)]
            [
                                             [1 0 0 0 0 0]
                                [1 0 0 0 0]  [0 1 0 0 0 0]
                     [1 0 0 0]  [0 1 0 0 0]  [0 0 1 0 0 0]
            [1 0 0]  [0 1 0 0]  [0 0 1 0 0]  [0 0 0 1 0 0]
            [0 1 0]  [0 0 1 0]  [0 0 0 1 0]  [0 0 0 0 1 0]
            [0 0 1], [0 0 0 1], [0 0 0 0 1], [0 0 0 0 0 1]
            ]
            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.run_cell('%display ascii_art')   # indirect doctest
            sage: shell.run_cell("i = var('i')")
            sage: shell.run_cell('sum(i*x^i, i, 0, 10)')
                10      9      8      7      6      5      4      3      2
            10*x   + 9*x  + 8*x  + 7*x  + 6*x  + 5*x  + 4*x  + 3*x  + 2*x  + x
            sage: shell.run_cell('%display default')
            sage: shell.quit()

        TESTS::

            sage: import os
            sage: from sage.env import SAGE_EXTCODE
            sage: example_png = os.path.join(SAGE_EXTCODE, 'doctest', 'rich_output', 'example.png')
            sage: from sage.repl.rich_output.backend_ipython import BackendIPython
            sage: backend = BackendIPython()
            sage: shell = get_test_shell()
            sage: backend.install(shell=shell)
            sage: shell.run_cell('get_ipython().display_formatter')
            <sage.repl.display.formatter.SageDisplayFormatter object at 0x...>
            sage: shell.run_cell('from IPython.display import Image')
            sage: shell.run_cell('ipython_image = Image("{0}")'.format(example_png))
            sage: shell.run_cell('ipython_image')
            <IPython.core.display.Image object>
            sage: shell.run_cell('get_ipython().display_formatter.format(ipython_image)')
            ({u'image/png': '\x89PNG...',
              u'text/plain': u'<IPython.core.display.Image object>'},
            {})

        Test that ``__repr__`` is only called once when generating text output::

            sage: class Repper(object):
            ....:    def __repr__(self):
            ....:        print('__repr__ called')
            ....:        return 'I am repper'
            sage: Repper()
            __repr__ called
            I am repper
        """
        # First, use Sage rich output if there is any
        PLAIN_TEXT = u'text/plain'
        sage_format, sage_metadata = self.dm.displayhook(obj)
        assert PLAIN_TEXT in sage_format, 'plain text is always present'
        if sage_format.keys() != [PLAIN_TEXT]:
            return sage_format, sage_metadata
        # Second, try IPython widgets (obj._ipython_display_ and type registry)
        if self.ipython_display_formatter(obj):
            return {}, {}
        # Finally, try IPython rich representation (obj._repr_foo_ methods)
        if exclude is not None:
            exclude = list(exclude) + [PLAIN_TEXT]
        else:
            exclude = [PLAIN_TEXT]
        ipy_format, ipy_metadata = super(SageDisplayFormatter, self).format(
            obj, include=include, exclude=exclude)
        ipy_format.update(sage_format)
        ipy_metadata.update(sage_metadata)
        return ipy_format, ipy_metadata


class SagePlainTextFormatter(PlainTextFormatter):

    def __init__(self, *args, **kwds):
        r"""
        Improved plain text IPython formatter.

        In particular, it correctly print lists of matrices or other
        objects (see
        :meth:`sage.structure.parent.Parent._repr_option`).

        .. warning::

            This IPython formatter is NOT used. You could use it to
            enable Sage formatting in IPython, but Sage uses its own
            rich output system that is more flexible and supports
            different backends.

        INPUT/OUTPUT:

        See the IPython documentation.

        EXAMPLES::

            sage: from sage.repl.interpreter import get_test_shell
            sage: shell = get_test_shell()
            sage: shell.display_formatter.formatters['text/plain']
            <IPython.core.formatters.PlainTextFormatter object at 0x...>
            sage: shell.quit()
        """
        super(SagePlainTextFormatter, self).__init__(*args, **kwds)

    def __call__(self, obj):
        """
        Compute the pretty representation of the object.

        Adapted from ``IPython.core.formatters.PlainTextPrettyPrint``.

        INPUT:

        - ``obj`` -- anything.

        OUTPUT:

        String. The plain text representation.

        EXAMPLES::

            sage: from sage.repl.display.formatter import SagePlainTextFormatter
            sage: fmt = SagePlainTextFormatter()
            sage: fmt
            <sage.repl.display.formatter.SagePlainTextFormatter object at 0x...>
            sage: fmt(2)
            ---- calling ipython formatter ----
            '2'
            sage: a = identity_matrix(ZZ, 2)
            sage: fmt([a, a])
            ---- calling ipython formatter ----
            '[\n[1 0]  [1 0]\n[0 1], [0 1]\n]'
        """
        from sage.doctest import DOCTEST_MODE
        if DOCTEST_MODE:
            # Just to show that this is never executed in any other doctests in the Sage library
            print('---- calling ipython formatter ----')
        from six import StringIO
        stream = StringIO()
        printer = SagePrettyPrinter(
            stream, self.max_width, unicode_to_str(self.newline))
        printer.pretty(obj)
        printer.flush()
        return stream.getvalue()

