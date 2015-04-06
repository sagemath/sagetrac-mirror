# -*- encoding: utf-8 -*-
"""
GNU Compiler Suite

These are mostly for testing purposes.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.interfaces.cmdline import CompilerTool
from sage.interfaces.cmdline.posix import test_nonzero_output


def test_compile_c(tool):
    """
    Test that a tool compiles C code

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.gcc import test_compile_c
        sage: test_compile_c(cat)
    """
    src = """
    int main(int argc, char **argv) 
    {
        return 0;
    }
    """
    binary = tool(src)
    assert len(binary) > 0, \
        '{0} ran without error but did not produce output'.format(tool)


gcc = CompilerTool(
    'gcc',
    ['gcc', '-o', 'output', 'source.c'],
    'Part of GNU GCC',
    test_compile_c,
    filename=['source.c', 'output'],
    rpm='gcc',
)
