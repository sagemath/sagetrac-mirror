"""
    trim doctest flags
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    transforms for code-blocks.

    modified from sphinx.transforms.post_transforms.code

    :copyright: Copyright 2007-2019 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import sys
from typing import NamedTuple

from docutils import nodes
from pygments.lexers import PythonConsoleLexer, guess_lexer

from sphinx import addnodes
from sphinx.ext import doctest
from sphinx.transforms import SphinxTransform

class TrimDoctestFlagsTransformSage(SphinxTransform):
    """
    Trim doctest flags like ``# doctest: +FLAG`` from sage code-blocks.

    The only difference from Sphinx' TrimDoctestFlagsTransform is that
    here lines starting with 'sage:' rather than '>>>' are transformed.
    """
    default_priority = 401

    def apply(self, **kwargs):
        # type: (Any) -> None
        if not self.config.trim_doctest_flags:
            return

        for node in self.document.traverse(nodes.literal_block):
            if self.is_pyconsole(node):
                source = node.rawsource
                source = doctest.blankline_re.sub('', source)
                source = doctest.doctestopt_re.sub('', source)
                node.rawsource = source
                node[:] = [nodes.Text(source)]

    @staticmethod
    def is_pyconsole(node):
        # type: (nodes.literal_block) -> bool
        if node.rawsource != node.astext():
            return False  # skip parsed-literal node

        language = node.get('language')
        if language in ('pycon', 'pycon3'):
            return True
        elif language in ('py', 'py3', 'python', 'python3', 'default'):
            return node.rawsource.startswith('sage:')
        elif language == 'guess':
            try:
                lexer = guess_lexer(node.rawsource)
                return isinstance(lexer, PythonConsoleLexer)
            except Exception:
                pass

        return False


def setup(app):
    app.add_post_transform(TrimDoctestFlagsTransformSage)
