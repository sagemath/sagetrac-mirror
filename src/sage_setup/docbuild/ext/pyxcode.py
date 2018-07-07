"""
Reimplements and/or extends some source code analysis utilities from Sphinx
such that they provide a modicum of support for Cython modules.

These are used to implement :mod:`sage_setup.docbuild.ext.sage_viewcode`.
"""

import os
import sys

import six

from sphinx import highlighting
from sphinx.errors import PycodeError
from sphinx.pycode import ModuleAnalyzer, Parser

from pygments.lexers import CythonLexer
from pygments.token import Operator

from sage.misc.sageinspect import sage_getfile, sage_getsource


CYTHON_EXTS = ('.pyx', '.pxd', '.pxi')


class _CythonLexer(CythonLexer):
    """Patched CythonLexer to support the matmul @ operator."""

    tokens = {'root': CythonLexer.tokens['root'][:]}
    tokens['root'].insert(tokens['root'].index('keywords'),
                          ('@', Operator))

# patch Sphinx to use our patched Cython lexer
highlighting.lexers['cython'] = _CythonLexer(stripnl=False)


def sage_get_full_modname(modname, attribute):
    """
    Same as Sphinx's `sphinx.util.get_full_modname` but supports method
    descriptors on built-in types.

    The name of the original function is something of a misnomer.  It should
    really be called something like ``get_real_modname``.  The purpose of this
    function is, given the name of a module and a member of that module, to
    return the name of the actual module in which that member was defined (e.g.
    in the case of classes that were imported from other modules, etc.

    EXAMPLES:

    Demonstrate that this function works on built-in method descriptors on
    built-in types (whereas the original function did not)::

        sage: from sphinx.util import get_full_modname
        sage: from sage_setup.docbuild.ext.pyxcode import (
        ....:     sage_get_full_modname)
        sage: get_full_modname('sage.structure.sage_object',
        ....:                  'SageObject.parent') is None
        True
        sage: sage_get_full_modname('sage.structure.sage_object',
        ....:                       'SageObject.parent')
        'sage.structure.sage_object'
    """

    # Original source from https://github.com/sphinx-doc/sphinx/blob/v1.7.5/sphinx/util/__init__.py
    # See sage_viewcode.py in the same directory as this module for the license
    # of the original code.
    if modname is None:
        # Prevents a TypeError: if the last getattr() call will return None
        # then it's better to return it directly
        return None
    __import__(modname)
    module = sys.modules[modname]

    # Allow an attribute to have multiple parts and incidentially allow
    # repeated .s in the attribute.
    value = module
    for attr in attribute.split('.'):
        if attr:
            value = getattr(value, attr)

    mod = getattr(value, '__module__', None)
    if mod is None and hasattr(value, '__objclass__'):
        # There's no reason this should ever be missing that I can think of
        # but use getattr just in case...
        mod = getattr(value.__objclass__, '__module__', None)

    return mod


def sage_get_module_source(modname):
    """
    Like `sphinx.util.get_module_source` but uses the `sage.misc.sage_inspect`
    module and supports Cython modules.

    EXAMPLES::

        sage: from sage_setup.docbuild.ext.pyxcode import (
        ....:     sage_get_module_source)
        sage: sage_get_module_source('sage.structure.sage_object')
        ('file', '...sage/structure/sage_object.pyx')
    """

    if modname not in sys.modules:
        try:
            __import__(modname)
        except Exception as exc:
            raise PycodeError('error importing {!r}'.format(modname), exc)

    mod = sys.modules[modname]

    try:
        filename = sage_getfile(mod)
    except Exception as exc:
        raise PycodeError('error getting filename for {!r}'.format(modname),
                exc)

    if filename and os.path.exists(filename):
        # Just return the source filename so we don't have to load its contents
        # into memory straight away
        return ('file', filename)

    try:
        source = sage_getsource(mod)
    except Exception as exc:
        raise PycodeError('error getting source for {!r}'.format(modname),
                exc)

    if not isinstance(source, six.text_type):
        source = source.decode('utf-8')

    return ('string', source)


class SageModuleAnalyzer(ModuleAnalyzer):
    """
    Subclass `sphinx.pycode.ModuleAnalyzer` to support loading and analyzing
    Cython modules.

    Currently does not fill the ``.tagorder`` or ``.attr_docs`` attibutes, only
    the ``.tags`` attribute needed to support the ``sage.ext.viewcode``
    extension.
    """

    # Ensure that SageModuleAnalyzer has its own cache independent from the
    # base class's cache
    cache = {}

    @classmethod
    def for_module(cls, modname):
        """
        Instantiate a :class:`SageModuleAnalyzer` given the name of a module.

        This fails in the base class for Cython modules.

        EXAMPLES::

            sage: from sage_setup.docbuild.ext.pyxcode import (
            ....:     SageModuleAnalyzer)
            sage: from sphinx.pycode import ModuleAnalyzer
            sage: ModuleAnalyzer.for_module('sage.structure.sage_object')
            Traceback (most recent call last):
            ...
            PycodeError: source is not a .py file: ...
            sage: analyzer = SageModuleAnalyzer.for_module(
            ....:     'sage.structure.sage_object')
            sage: analyzer.parse()
            sage: sorted(analyzer.find_tags().items())
            [(u'SageObject', ('class', ..., ...)),
             ...,
             (u'SageObject.save', ('def', ..., ...))]
        """
        # Reimplement ModuleAnalyzer.for_module go through
        # sage_get_module_source
        cache_key = ('module', modname)
        entry = cls.cache.get(cache_key)
        if entry is None:
            try:
                type_, source = sage_get_module_source(modname)
                if type_ == 'string':
                    entry = cls.for_string(source, modname)
                else:
                    entry = cls.for_file(source, modname)
            except PycodeError as exc:
                entry = exc

            cls.cache[cache_key] = entry

        if isinstance(entry, PycodeError):
            raise entry

        return entry

    def parse(self):
        """
        If the module is a Cython module, parse using the more limited
        :class:`SageParser`, otherwise fall back on the default.
        """
        ext = os.path.splitext(self.srcname)[1]
        if ext in CYTHON_EXTS:
            # Try our best with the hacked Parser
            try:
                parser = SageParser(self.code, self.encoding)
                parser.parse()
                self.tags = parser.definitions
            except Exception as exc:
                raise PycodeError('parsing {!r} failed: {!r}'.format(
                    self.srcname, exc))
        else:
            # A standard Python module, presumably; try the default
            # implementation
            super(SageModuleAnalyzer, self).parse()


class SageParser(Parser):
    """
    Simplified extension to `sphinx.pycode.Parser` which skips the full parsing
    into an AST and attribute comment loading, and just returns the list of
    "tags"--class and function definitions in a module.

    EXAMPLES::

        sage: from sage_setup.docbuild.ext.pyxcode import (
        ....:     SageParser)
        sage: from sage.misc.sageinspect import sage_getsource
        sage: from sphinx.pycode import Parser
        sage: import sage.structure.sage_object
        sage: import six
        sage: code = sage_getsource(sage.structure.sage_object)
        sage: if not isinstance(code, six.text_type):
        ....:     code = code.decode('utf-8')  # Parser expects unicode
        sage: parser = Parser(code, 'utf-8')
        sage: parser.parse()  # Tries to parse Cython and fails...
        Traceback (most recent call last):
        ...
        SyntaxError: invalid syntax
        sage: parser = SageParser(code, 'utf-8')
        sage: parser.parse()
        sage: sorted(parser.definitions.items())
        [(u'SageObject', ('class', ..., ...)),
         ...,
         (u'SageObject.save', ('def', ..., ...))]
    """

    def parse_comments(self):
        """
        No-op--would otherwise raise a ``SyntaxError`` on most Cython sources.
        """
