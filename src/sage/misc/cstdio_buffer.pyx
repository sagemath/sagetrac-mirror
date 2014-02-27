"""
Line Buffered Output

Running a program non-interactively (for example, through a shell
pipe) defaults to block buffered input/output. This is optimal for
throughput, but has the disadvantage that output only becomes visible
once the buffer is full. This module allows one to run the current
session in line buffer mode, just as if it were an interactive session
but without a virtual terminal. This lets us control Sage from another
session without a virtual terminal (which is comparatively slow and
prone to undesirable side effects like dependency on terminal
settings, terminal echo, ...)

EXAMPLES::

    sage: from sage.misc.cstdio_buffer import linebuffer_output
    sage: linebuffer_output()
    sage: 1+1
    2

TESTS::

    sage: from sage.misc.interpreter import get_test_shell
    sage: shell = get_test_shell()
    sage: linebuffer_output()
    sage: shell.run_cell('integral(x,x)')
    1/2*x^2
"""

from __future__ import print_function

cimport libc.stdio as std

import sys


cdef cstdio_line_buffer_mode():
    std.setvbuf(std.stdout, NULL, std._IOLBF, 0)
    std.setvbuf(std.stderr, NULL, std._IOLBF, 0)


class OutputStreamWrapper(object):

    def __init__(self, stream):
        """
        Wrapper for a Python stream to make it line buffered.

        INPUT:

        - ``stream`` -- output stream to wrap.

        EXAMPLES::

            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: import sys
            sage: OutputStreamWrapper(sys.stdout)
            <sage.misc.cstdio_buffer.OutputStreamWrapper object at 0x...>
        """
        self.__stream = stream

    @staticmethod
    def install():
        """
        Install the wrapper for Python stdout and stderr.

        EXAMPLES::

            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: OutputStreamWrapper.install()
        """
        if not isinstance(sys.stdout, OutputStreamWrapper):
            sys.stdout = OutputStreamWrapper(sys.stdout)
        if not isinstance(sys.stderr, OutputStreamWrapper):
            sys.stderr = OutputStreamWrapper(sys.stderr)

    def write(self, data):
        """
        Wrap the stream write method

        INPUT:

        - ``data`` -- string. To write to the stream

        OUTPUT:

        Integer. The number of characters written.

        EXAMPLES::

            sage: import sys
            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: stream = OutputStreamWrapper(sys.stdout)
            sage: stream.write('hello')
            hello
        """
        n = self.__stream.write(data)
        self.__stream.flush()
        return n

    def writelines(self, data):
        """
        Wrap the stream writelines method

        INPUT:

        - ``data`` -- list of strings. Lines to write to the stream.

        EXAMPLES::

            sage: import sys
            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: stream = OutputStreamWrapper(sys.stdout)
            sage: stream.writelines(['hello'])
            hello
        """
        self.__stream.writelines(data)
        self.__stream.flush()

    def isatty(self):
        """
        Wrap the stream isatty method

        OUTPUT:

        This method always returns ``True`` to fake a tty.

        EXAMPLES::

            sage: import sys
            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: stream = OutputStreamWrapper(sys.stdout)
            sage: stream.isatty()
            True
        """
        return True

    def flush(self):
        """
        Wrap the stream flush method

        EXAMPLES::

            sage: import sys
            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: stream = OutputStreamWrapper(sys.stdout)
            sage: stream.flush()
        """
        self.__stream.flush()

    def fileno(self):
        """
        Wrap the stream fileno method

        OUTPUT:

        Integer. The file descriptor of the stream.

        EXAMPLES::

            sage: import sys
            sage: from sage.misc.cstdio_buffer import OutputStreamWrapper
            sage: stream = OutputStreamWrapper(sys.stdout)
            sage: stream.fileno()
            1
        """
        return self.__stream.fileno()


def linebuffer_output():
    """
    Switch stdout/stderr to line buffer mode

    It is safe to call this function repeatedly.

    EXAMPLES::

        sage: from sage.misc.cstdio_buffer import linebuffer_output
        sage: linebuffer_output()
        sage: 1+1
        2
    """
    cstdio_line_buffer_mode()
    OutputStreamWrapper.install()

