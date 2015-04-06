# -*- encoding: utf-8 -*-
"""
Standard Commandline Tools

These are mostly for testing purposes, the Python standard library
provides generally better functionality than calling external posix
command line tools.
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sage.interfaces.cmdline import Tool, ToolAdapter


def test_nonzero_output(tool):
    """
    Test that a tool returns anything

    INPUT:

    - ``tool`` -- :class:`sage.interfaces.cmdline.tool.Tool`
      instance. The tool to test.

    EXAMPLES::

        sage: from sage.interfaces.cmdline.posix import cat
        sage: from sage.interfaces.cmdline.posix import test_nonzero_output
        sage: test_nonzero_output(cat)
    """
    output = tool('this is some input')
    assert len(output) > 0, \
        '{0} ran without error but did not produce output'.format(tool)


cat = Tool(
    'cat',
    ['cat'],
    'Part of GNU coreutils',
    test_nonzero_output,
    rpm='coreutils',
)


cp = ToolAdapter(
    'cp',
    ['cp', '{input}', '{output}'],
    'Part of GNU coreutils',
    test_nonzero_output,
    rpm='coreutils',
)
