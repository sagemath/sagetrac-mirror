# -*- coding: utf-8 -*-
"""
Sage autodoc extension based on sphinx.ext.autodoc from Sphinx

From sphinx.ext.autodoc:

    Automatically insert docstrings for functions, classes or whole modules into
    the doctree, thus avoiding duplication between docstrings and documentation
    for those who like elaborate docstrings.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.

The upstream original can be found at
https://github.com/sphinx-doc/sphinx/blob/master/sphinx/ext/autodoc/__init__.py

AUTHORS:

- ? (before 2011): initial version derived from sphinx.ext.autodoc

- Johan S. R. Nielsen: support for _sage_argspec_

- Simon King (2011-04): use sageinspect; include public cython attributes
  only in the documentation if they have a docstring

- Kwankyu Lee (2018-12-26): rebase on the latest sphinx.ext.autodoc

"""

import inspect
import re
import sys

from docutils.statemachine import ViewList

import sphinx
from sphinx.ext.autodoc import mock, ObjectMember
from sphinx.ext.autodoc.importer import import_object, get_object_members, get_module_members
from sphinx.locale import _, __
from sphinx.pycode import ModuleAnalyzer
from sphinx.errors import PycodeError
from sphinx.util import logging
from sphinx.util.docstrings import prepare_docstring
from sphinx.util.inspect import isdescriptor, \
    safe_getattr, object_description, is_builtin_class_method, \
    isenumattribute, isclassmethod, isstaticmethod, getdoc

from sage.misc.sageinspect import (sage_getdoc_original,
                                   sage_getargspec, isclassinstance,
                                   sage_formatargspec,
                                   is_function_or_cython_function)
from sage.misc.lazy_import import LazyImport

# This is used to filter objects of classes that inherit from
# ClasscallMetaclass. See the class AttributeDocumenter below.
from sage.misc.classcall_metaclass import ClasscallMetaclass


logger = logging.getLogger(__name__)


# This type isn't exposed directly in any modules, but can be found
# here in most Python versions
MethodDescriptorType = type(type.__subclasses__)


#: extended signature RE: with explicit module name separated by ::
py_ext_sig_re = re.compile(
    r'''^ ([\w.]+::)?            # explicit module name
          ([\w.]+\.)?            # module and/or class name(s)
          (\w+)  \s*             # thing name
          (?: \((.*)\)           # optional: arguments
           (?:\s* -> \s* (.*))?  #           return annotation
          )? $                   # and nothing more
          ''', re.VERBOSE)


def identity(x):
    # type: (Any) -> Any
    return x


ALL = object()
INSTANCEATTR = object()


def members_option(arg):
    # type: (Any) -> Union[object, List[unicode]]
    """Used to convert the :members: option to auto directives."""
    if arg is None:
        return ALL
    return [x.strip() for x in arg.split(',')]


def members_set_option(arg):
    # type: (Any) -> Union[object, Set[unicode]]
    """Used to convert the :members: option to auto directives."""
    if arg is None:
        return ALL
    return set(x.strip() for x in arg.split(','))


SUPPRESS = object()


def annotation_option(arg):
    # type: (Any) -> Any
    if arg is None:
        # suppress showing the representation of the object
        return SUPPRESS
    else:
        return arg


def bool_option(arg):
    # type: (Any) -> bool
    """Used to convert flag options to auto directives.  (Instead of
    directives.flag(), which returns None).
    """
    return True


def formatargspec(function, args, varargs=None, varkw=None, defaults=None,
                  kwonlyargs=(), kwonlydefaults={}, annotations={}):
    """
    Sphinx's version of formatargspec is deprecated, so use Sage's instead.
    """
    return sage_formatargspec(args, varargs=varargs, varkw=varkw, defaults=defaults,
                              kwonlyargs=kwonlyargs, kwonlydefaults=kwonlydefaults,
                              annotations=annotations)


class AutodocReporter(object):
    """
    A reporter replacement that assigns the correct source name
    and line number to a system message, as recorded in a ViewList.
    """
    def __init__(self, viewlist, reporter):
        # type: (ViewList, Reporter) -> None
        self.viewlist = viewlist
        self.reporter = reporter

    def __getattr__(self, name):
        # type: (unicode) -> Any
        return getattr(self.reporter, name)

    def system_message(self, level, message, *children, **kwargs):
        # type: (int, unicode, Any, Any) -> nodes.system_message
        if 'line' in kwargs and 'source' not in kwargs:
            try:
                source, line = self.viewlist.items[kwargs['line']]
            except IndexError:
                pass
            else:
                kwargs['source'] = source
                kwargs['line'] = line
        return self.reporter.system_message(level, message,
                                            *children, **kwargs)

    def debug(self, *args, **kwargs):
        # type: (Any, Any) -> nodes.system_message
        if self.reporter.debug_flag:
            return self.system_message(0, *args, **kwargs)

    def info(self, *args, **kwargs):
        # type: (Any, Any) -> nodes.system_message
        return self.system_message(1, *args, **kwargs)

    def warning(self, *args, **kwargs):
        # type: (Any, Any) -> nodes.system_message
        return self.system_message(2, *args, **kwargs)

    def error(self, *args, **kwargs):
        # type: (Any, Any) -> nodes.system_message
        return self.system_message(3, *args, **kwargs)

    def severe(self, *args, **kwargs):
        # type: (Any, Any) -> nodes.system_message
        return self.system_message(4, *args, **kwargs)


# Some useful event listener factories for autodoc-process-docstring.

def cut_lines(pre, post=0, what=None):
    # type: (int, int, unicode) -> Callable
    """Return a listener that removes the first *pre* and last *post*
    lines of every docstring.  If *what* is a sequence of strings,
    only docstrings of a type in *what* will be processed.

    Use like this (e.g. in the ``setup()`` function of :file:`conf.py`)::

       from sphinx.ext.autodoc import cut_lines
       app.connect('autodoc-process-docstring', cut_lines(4, what=['module']))

    This can (and should) be used in place of :confval:`automodule_skip_lines`.
    """
    def process(app, what_, name, obj, options, lines):
        # type: (Sphinx, unicode, unicode, Any, Any, List[unicode]) -> None
        if what and what_ not in what:
            return
        del lines[:pre]
        if post:
            # remove one trailing blank line.
            if lines and not lines[-1]:
                lines.pop(-1)
            del lines[-post:]
        # make sure there is a blank line at the end
        if lines and lines[-1]:
            lines.append('')
    return process


def between(marker, what=None, keepempty=False, exclude=False):
    # type: (unicode, Sequence[unicode], bool, bool) -> Callable
    """Return a listener that either keeps, or if *exclude* is True excludes,
    lines between lines that match the *marker* regular expression.  If no line
    matches, the resulting docstring would be empty, so no change will be made
    unless *keepempty* is true.

    If *what* is a sequence of strings, only docstrings of a type in *what* will
    be processed.
    """
    marker_re = re.compile(marker)

    def process(app, what_, name, obj, options, lines):
        # type: (Sphinx, unicode, unicode, Any, Any, List[unicode]) -> None
        if what and what_ not in what:
            return
        deleted = 0
        delete = not exclude
        orig_lines = lines[:]
        for i, line in enumerate(orig_lines):
            if delete:
                lines.pop(i - deleted)
                deleted += 1
            if marker_re.match(line):
                delete = not delete
                if delete:
                    lines.pop(i - deleted)
                    deleted += 1
        if not lines and not keepempty:
            lines[:] = orig_lines
        # make sure there is a blank line at the end
        if lines and lines[-1]:
            lines.append('')
    return process




class Documenter(object):
    """
    A Documenter knows how to autodocument a single object type.  When
    registered with the AutoDirective, it will be used to document objects
    of that type when needed by autodoc.

    Its *objtype* attribute selects what auto directive it is assigned to
    (the directive name is 'auto' + objtype), and what directive it generates
    by default, though that can be overridden by an attribute called
    *directivetype*.

    A Documenter has an *option_spec* that works like a docutils directive's;
    in fact, it will be used to parse an auto directive's options that matches
    the documenter.
    """
    #: name by which the directive is called (auto...) and the default
    #: generated directive name
    objtype = 'object'
    #: indentation by which to indent the directive content
    content_indent = u'   '
    #: priority if multiple documenters return True from can_document_member
    priority = 0
    #: order if autodoc_member_order is set to 'groupwise'
    member_order = 0
    #: true if the generated content may contain titles
    titles_allowed = False

    option_spec = {'noindex': bool_option}  # type: Dict[unicode, Callable]

    def get_attr(self, obj, name, *defargs):
        # type: (Any, unicode, Any) -> Any
        """getattr() override for types such as Zope interfaces."""
        return autodoc_attrgetter(self.env.app, obj, name, *defargs)

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        """Called to see if a member can be documented by this documenter."""
        raise NotImplementedError('must be implemented in subclasses')

    def __init__(self, directive, name, indent=u''):
        # type: (DocumenterBridge, unicode, unicode) -> None
        self.directive = directive
        self.env = directive.env    # type: BuildEnvironment
        self.options = directive.genopt
        self.name = name
        self.indent = indent
        # the module and object path within the module, and the fully
        # qualified name (all set after resolve_name succeeds)
        self.modname = None         # type: str
        self.module = None          # type: ModuleType
        self.objpath = None         # type: List[unicode]
        self.fullname = None        # type: unicode
        # extra signature items (arguments and return annotation,
        # also set after resolve_name succeeds)
        self.args = None            # type: unicode
        self.retann = None          # type: unicode
        # the object to document (set after import_object succeeds)
        self.object = None          # type: Any
        self.object_name = None     # type: unicode
        # the parent/owner of the object to document
        self.parent = None          # type: Any
        # the module analyzer to get at attribute docs, or None
        self.analyzer = None        # type: Any

    @property
    def documenters(self):
        # type: () -> Dict[unicode, Type[Documenter]]
        """Returns registered Documenter classes"""
        return get_documenters(self.env.app)

    def add_line(self, line, source, *lineno):
        # type: (unicode, unicode, int) -> None
        """Append one line of generated reST to the output."""
        self.directive.result.append(self.indent + line, source, *lineno)

    def resolve_name(self, modname, parents, path, base):
        # type: (str, Any, str, Any) -> Tuple[str, List[unicode]]
        """Resolve the module and name of the object to document given by the
        arguments and the current module/class.

        Must return a pair of the module name and a chain of attributes; for
        example, it would return ``('zipfile', ['ZipFile', 'open'])`` for the
        ``zipfile.ZipFile.open`` method.
        """
        raise NotImplementedError('must be implemented in subclasses')

    def parse_name(self):
        # type: () -> bool
        """Determine what module to import and what attribute to document.

        Returns True and sets *self.modname*, *self.objpath*, *self.fullname*,
        *self.args* and *self.retann* if parsing and resolving was successful.
        """
        # first, parse the definition -- auto directives for classes and
        # functions can contain a signature which is then used instead of
        # an autogenerated one
        try:
            explicit_modname, path, base, args, retann = \
                py_ext_sig_re.match(self.name).groups()  # type: ignore
        except AttributeError:
            logger.warning('invalid signature for auto%s (%r)' % (self.objtype, self.name))
            return False

        # support explicit module and class name separation via ::
        if explicit_modname is not None:
            modname = explicit_modname[:-2]
            parents = path and path.rstrip('.').split('.') or []
        else:
            modname = None
            parents = []

        self.modname, self.objpath = self.resolve_name(modname, parents, path, base)

        if not self.modname:
            return False

        self.args = args
        self.retann = retann
        self.fullname = (self.modname or '') + \
                        (self.objpath and '.' + '.'.join(self.objpath) or '')
        return True

    def import_object(self):
        # type: () -> bool
        """Import the object given by *self.modname* and *self.objpath* and set
        it as *self.object*.

        Returns True if successful, False if an error occurred.
        """
        with mock(self.env.config.autodoc_mock_imports):
            try:
                ret = import_object(self.modname, self.objpath, self.objtype,
                                    attrgetter=self.get_attr,
                                    warningiserror=self.env.config.autodoc_warningiserror)
                self.module, self.parent, self.object_name, self.object = ret
                return True
            except ImportError as exc:
                logger.warning(exc.args[0], type='autodoc', subtype='import_object')
                self.env.note_reread()
                return False

    def get_real_modname(self):
        # type: () -> str
        """Get the real module name of an object to document.

        It can differ from the name of the module through which the object was
        imported.
        """
        return self.get_attr(self.object, '__module__', None) or self.modname

    def check_module(self):
        # type: () -> bool
        """Check if *self.object* is really defined in the module given by
        *self.modname*.
        """
        if self.options.imported_members:
            return True

        modname = self.get_attr(self.object, '__module__', None)
        if modname and modname != self.modname:
            return False
        return True

    def format_args(self):
        # type: () -> unicode
        """Format the argument signature of *self.object*.

        Should return None if the object does not have a signature.
        """
        return None

    def format_name(self):
        # type: () -> unicode
        """Format the name of *self.object*.

        This normally should be something that can be parsed by the generated
        directive, but doesn't need to be (Sphinx will display it unparsed
        then).
        """
        # normally the name doesn't contain the module (except for module
        # directives of course)
        return '.'.join(self.objpath) or self.modname

    def format_signature(self):
        # type: () -> unicode
        """Format the signature (arguments and return annotation) of the object.

        Let the user process it via the ``autodoc-process-signature`` event.
        """
        if self.args is not None:
            # signature given explicitly
            args = "(%s)" % self.args  # type: unicode
        else:
            # try to introspect the signature
            try:
                args = self.format_args()
            except Exception as err:
                logger.warning('error while formatting arguments for %s: %s' %
                               (self.fullname, err))
                args = None

        retann = self.retann

        result = self.env.app.emit_firstresult(
            'autodoc-process-signature', self.objtype, self.fullname,
            self.object, self.options, args, retann)
        if result:
            args, retann = result

        if args is not None:
            return args + (retann and (' -> %s' % retann) or '')
        else:
            return ''

    def add_directive_header(self, sig):
        # type: (unicode) -> None
        """Add the directive header and options to the generated content."""
        domain = getattr(self, 'domain', 'py')
        directive = getattr(self, 'directivetype', self.objtype)
        name = self.format_name()
        sourcename = self.get_sourcename()
        self.add_line(u'.. %s:%s:: %s%s' % (domain, directive, name, sig),
                      sourcename)
        if self.options.noindex:
            self.add_line(u'   :noindex:', sourcename)
        if self.objpath:
            # Be explicit about the module, this is necessary since .. class::
            # etc. don't support a prepended module name
            self.add_line(u'   :module: %s' % self.modname, sourcename)

    def get_doc(self, encoding=None, ignore=1):
        # type: (unicode, int) -> List[List[unicode]]
        """Decode and return lines of the docstring(s) for the object."""
        docstring = sage_getdoc_original(self.object)
        if docstring is None and self.env.config.autodoc_inherit_docstrings:
            docstring = getdoc(self.object)
        # make sure we have Unicode docstrings, then sanitize and split
        # into lines
        if isinstance(docstring, str):
            return [prepare_docstring(docstring)]
        # ... else it is something strange, let's ignore it
        return []

    def process_doc(self, docstrings):
        # type: (List[List[unicode]]) -> Iterator[unicode]
        """Let the user process the docstrings before adding them."""
        for docstringlines in docstrings:
            if self.env.app:
                # let extensions preprocess docstrings
                self.env.app.emit('autodoc-process-docstring',
                                  self.objtype, self.fullname, self.object,
                                  self.options, docstringlines)
            for line in docstringlines:
                yield line

    def get_sourcename(self):
        # type: () -> unicode
        if self.analyzer:
            # prevent encoding errors when the file name is non-ASCII
            if not isinstance(self.analyzer.srcname, str):
                filename = str(self.analyzer.srcname,
                               sys.getfilesystemencoding(), 'replace')
            else:
                filename = self.analyzer.srcname
            return u'%s:docstring of %s' % (filename, self.fullname)
        return u'docstring of %s' % self.fullname

    def add_content(self, more_content, no_docstring=False):
        # type: (Any, bool) -> None
        """Add content from docstrings, attribute documentation and user."""
        # set sourcename and add content from attribute documentation
        sourcename = self.get_sourcename()
        if self.analyzer:
            attr_docs = self.analyzer.find_attr_docs()
            if self.objpath:
                key = ('.'.join(self.objpath[:-1]), self.objpath[-1])
                if key in attr_docs:
                    no_docstring = True
                    docstrings = [attr_docs[key]]
                    for i, line in enumerate(self.process_doc(docstrings)):
                        self.add_line(line, sourcename, i)

        # add content from docstrings
        if not no_docstring:
            encoding = self.analyzer
            docstrings = self.get_doc(encoding)
            if not docstrings:
                # append at least a dummy docstring, so that the event
                # autodoc-process-docstring is fired and can add some
                # content if desired
                docstrings.append([])
            for i, line in enumerate(self.process_doc(docstrings)):
                self.add_line(line, sourcename, i)

        # add additional content (e.g. from document), if present
        if more_content:
            for line, src in zip(more_content.data, more_content.items):
                self.add_line(line, src[0], src[1])

    def get_object_members(self, want_all):
        # type: (bool) -> Tuple[bool, List[Tuple[unicode, Any]]]
        """Return `(members_check_module, members)` where `members` is a
        list of `(membername, member)` pairs of the members of *self.object*.

        If *want_all* is True, return all members.  Else, only return those
        members given by *self.options.members* (which may also be none).
        """
        members = get_object_members(self.object, self.objpath, self.get_attr, self.analyzer)
        if not want_all:
            if not self.options.members:
                return False, []
            # specific members given
            selected = []
            for name in self.options.members:
                if name in members:
                    selected.append((name, members[name].value))
                else:
                    logger.warning('missing attribute %s in object %s' %
                                   (name, self.fullname))
            return False, sorted(selected)
        elif self.options.inherited_members:
            return False, sorted((m.name, m.value) for m in members.values())
        else:
            return False, sorted((m.name, m.value) for m in members.values()
                                 if m.directly_defined)

    def filter_members(self, members, want_all):
        # type: (List[Tuple[unicode, Any]], bool) -> List[Tuple[unicode, Any, bool]]
        """Filter the given member list.

        Members are skipped if

        - they are private (except if given explicitly or the private-members
          option is set)
        - they are special methods (except if given explicitly or the
          special-members option is set)
        - they are undocumented (except if the undoc-members option is set)

        The user can override the skipping decision by connecting to the
        ``autodoc-skip-member`` event.
        """
        ret = []

        # search for members in source code too
        namespace = '.'.join(self.objpath)  # will be empty for modules

        if self.analyzer:
            attr_docs = self.analyzer.find_attr_docs()
        else:
            attr_docs = {}

        # process members and determine which to skip
        for (membername, member) in members:
            # Immediately skip lazy imports to avoid deprecation
            # messages (#17455).
            if isinstance(member, LazyImport):
                continue

            # if isattr is True, the member is documented as an attribute
            isattr = False

            doc = self.get_attr(member, '__doc__', None)
            if doc is None and self.env.config.autodoc_inherit_docstrings:
                doc = getdoc(member)

            # if the member __doc__ is the same as self's __doc__, it's just
            # inherited and therefore not the member's doc
            cls = self.get_attr(member, '__class__', None)
            if cls:
                cls_doc = self.get_attr(cls, '__doc__', None)
                if cls_doc == doc:
                    doc = None
            has_doc = bool(doc)

            keep = False
            if want_all and membername.startswith('__') and \
                    membername.endswith('__') and len(membername) > 4:
                # special __methods__
                if self.options.special_members is ALL and \
                        membername != '__doc__':
                    keep = has_doc or self.options.undoc_members
                elif self.options.special_members and \
                    self.options.special_members is not ALL and \
                        membername in self.options.special_members:
                    keep = has_doc or self.options.undoc_members
            elif (namespace, membername) in attr_docs:
                if want_all and membername.startswith('_'):
                    # ignore members whose name starts with _ by default
                    keep = self.options.private_members
                else:
                    # keep documented attributes
                    keep = True
                isattr = True
            elif want_all and membername.startswith('_'):
                # ignore members whose name starts with _ by default
                keep = self.options.private_members and \
                    (has_doc or self.options.undoc_members)
            else:
                # ignore undocumented members if :undoc-members: is not given
                keep = has_doc or self.options.undoc_members

            # give the user a chance to decide whether this member
            # should be skipped
            if self.env.app:
                # let extensions preprocess docstrings
                try:
                    skip_user = self.env.app.emit_firstresult(
                        'autodoc-skip-member', self.objtype, membername, member,
                        not keep, self.options)
                    if skip_user is not None:
                        keep = not skip_user
                except Exception as exc:
                    logger.warning(__('autodoc: failed to determine %r to be documented.'
                                      'the following exception was raised:\n%s'),
                                   member, exc)
                    keep = False

            if keep:
                ret.append((membername, member, isattr))

        return ret

    def document_members(self, all_members=False):
        # type: (bool) -> None
        """Generate reST for member documentation.

        If *all_members* is True, do all members, else those given by
        *self.options.members*.
        """
        # set current namespace for finding members
        self.env.temp_data['autodoc:module'] = self.modname
        if self.objpath:
            self.env.temp_data['autodoc:class'] = self.objpath[0]

        want_all = all_members or self.options.inherited_members or \
            self.options.members is ALL
        # find out which members are documentable
        members_check_module, members = self.get_object_members(want_all)

        # remove members given by exclude-members
        if self.options.exclude_members:
            members = [(membername, member) for (membername, member) in members
                       if membername not in self.options.exclude_members]

        # document non-skipped members
        memberdocumenters = []  # type: List[Tuple[Documenter, bool]]
        for (mname, member, isattr) in self.filter_members(members, want_all):
            classes = [cls for cls in self.documenters.values()
                       if cls.can_document_member(member, mname, isattr, self)]
            if not classes:
                # don't know how to document this member
                continue
            # prefer the documenter with the highest priority
            classes.sort(key=lambda cls: cls.priority)
            # give explicitly separated module name, so that members
            # of inner classes can be documented
            full_mname = self.modname + '::' + \
                '.'.join(self.objpath + [mname])
            documenter = classes[-1](self.directive, full_mname, self.indent)
            memberdocumenters.append((documenter, isattr))
        member_order = self.options.member_order or \
            self.env.config.autodoc_member_order
        if member_order == 'groupwise':
            # sort by group; relies on stable sort to keep items in the
            # same group sorted alphabetically
            memberdocumenters.sort(key=lambda e: e[0].member_order)
        elif member_order == 'bysource' and self.analyzer:
            # sort by source order, by virtue of the module analyzer
            tagorder = self.analyzer.tagorder

            def keyfunc(entry):
                # type: (Tuple[Documenter, bool]) -> int
                fullname = entry[0].name.split('::')[1]
                return tagorder.get(fullname, len(tagorder))
            memberdocumenters.sort(key=keyfunc)

        for documenter, isattr in memberdocumenters:
            documenter.generate(
                all_members=True, real_modname=self.real_modname,
                check_module=members_check_module and not isattr)

        # reset current objects
        self.env.temp_data['autodoc:module'] = None
        self.env.temp_data['autodoc:class'] = None

    def generate(self, more_content=None, real_modname=None,
                 check_module=False, all_members=False):
        # type: (Any, str, bool, bool) -> None
        """Generate reST for the object given by *self.name*, and possibly for
        its members.

        If *more_content* is given, include that content. If *real_modname* is
        given, use that module name to find attribute docs. If *check_module* is
        True, only generate if the object is defined in the module name it is
        imported from. If *all_members* is True, document all members.
        """
        if not self.parse_name():
            # need a module to import
            logger.warning(
                'don\'t know which module to import for autodocumenting '
                '%r (try placing a "module" or "currentmodule" directive '
                'in the document, or giving an explicit module name)' %
                self.name)
            return

        # now, import the module and get object to document
        if not self.import_object():
            return

        # If there is no real module defined, figure out which to use.
        # The real module is used in the module analyzer to look up the module
        # where the attribute documentation would actually be found in.
        # This is used for situations where you have a module that collects the
        # functions and classes of internal submodules.
        self.real_modname = real_modname or self.get_real_modname()  # type: str

        # try to also get a source code analyzer for attribute docs
        try:
            self.analyzer = ModuleAnalyzer.for_module(self.real_modname)
            # parse right now, to get PycodeErrors on parsing (results will
            # be cached anyway)
            self.analyzer.find_attr_docs()
        except PycodeError as err:
            logger.debug('[autodoc] module analyzer failed: %s', err)
            # no source file -- e.g. for builtin and C modules
            self.analyzer = None
            # at least add the module as a dependency
            name = self.module.__name__ if hasattr(self.module, '__name__') else None
            fname = self.module.__file__ if hasattr(self.module, '__file__') else None
            if name != self.real_modname:
                # !!! SageMath Specific for make fast-rebuild-clean !!!
                # try not to record a dependency to a .pyc file but to the corresponding .py files instead.
                try:
                    fname = ModuleAnalyzer.for_module(name).srcname
                except PycodeError:
                    pass
            if fname is not None:
                self.directive.filename_set.add(fname)
        else:
            self.directive.filename_set.add(self.analyzer.srcname)

        # check __module__ of object (for members not given explicitly)
        if check_module:
            if not self.check_module():
                return

        sourcename = self.get_sourcename()

        # make sure that the result starts with an empty line.  This is
        # necessary for some situations where another directive preprocesses
        # reST and no starting newline is present
        self.add_line('', sourcename)

        # format the object's signature, if any
        sig = self.format_signature()

        # generate the directive header and options, if applicable
        self.add_directive_header(sig)
        self.add_line('', sourcename)

        # e.g. the module directive doesn't have content
        self.indent += self.content_indent

        # add all content (from docstrings, attribute docs etc.)
        self.add_content(more_content)

        # document members, if possible
        self.document_members(all_members)


class ModuleDocumenter(Documenter):
    """
    Specialized Documenter subclass for modules.
    """
    objtype = 'module'
    content_indent = u''
    titles_allowed = True

    option_spec = {
        'members': members_option, 'undoc-members': bool_option,
        'noindex': bool_option, 'inherited-members': bool_option,
        'show-inheritance': bool_option, 'synopsis': identity,
        'platform': identity, 'deprecated': bool_option,
        'member-order': identity, 'exclude-members': members_set_option,
        'private-members': bool_option, 'special-members': members_option,
        'imported-members': bool_option,
    }  # type: Dict[unicode, Callable]

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        # don't document submodules automatically
        return False

    def resolve_name(self, modname, parents, path, base):
        # type: (str, Any, str, Any) -> Tuple[str, List[unicode]]
        if modname is not None:
            logger.warning('"::" in automodule name doesn\'t make sense')
        return (path or '') + base, []

    def parse_name(self):
        # type: () -> bool
        ret = Documenter.parse_name(self)
        if self.args or self.retann:
            logger.warning('signature arguments or return annotation '
                           'given for automodule %s' % self.fullname)
        return ret

    def add_directive_header(self, sig):
        # type: (unicode) -> None
        Documenter.add_directive_header(self, sig)

        sourcename = self.get_sourcename()

        # add some module-specific options
        if self.options.synopsis:
            self.add_line(
                u'   :synopsis: ' + self.options.synopsis, sourcename)
        if self.options.platform:
            self.add_line(
                u'   :platform: ' + self.options.platform, sourcename)
        if self.options.deprecated:
            self.add_line(u'   :deprecated:', sourcename)

    def get_module_members(self):
        """Get members of target module."""
        if self.analyzer:
            attr_docs = self.analyzer.attr_docs
        else:
            attr_docs = {}

        members = {}  # type: Dict[str, ObjectMember]
        for name in dir(self.object):
            try:
                value = safe_getattr(self.object, name, None)
                docstring = attr_docs.get(('', name), [])
                members[name] = ObjectMember(name, value, docstring="\n".join(docstring))
            except AttributeError:
                continue

        # annotation only member (ex. attr: int)
        try:
            for name in inspect.getannotations(self.object):
                if name not in members:
                    docstring = attr_docs.get(('', name), [])
                    members[name] = ObjectMember(name, INSTANCEATTR,
                                                 docstring="\n".join(docstring))
        except AttributeError:
            pass

        return members

    def get_object_members(self, want_all):
        # type: (bool) -> Tuple[bool, List[Tuple[unicode, object]]]
        members = self.get_module_members()
        if want_all:
            if not hasattr(self.object, '__all__'):
                # for implicit module members, check __module__ to avoid
                # documenting imported objects
                return True, list(members.values())
            else:
                memberlist = self.object.__all__
                # Sometimes __all__ is broken...
                if not isinstance(memberlist, (list, tuple)) or not \
                   all(isinstance(entry, str) for entry in memberlist):
                    logger.warning(
                        '__all__ should be a list of strings, not %r '
                        '(in module %s) -- ignoring __all__' %
                        (memberlist, self.fullname))
                    # fall back to all members
                    return True, list(members.values())
        else:
            memberlist = self.options.members or []
        ret = []
        for mname in memberlist:
            if mname in members:
                ret.append(members[mname])
            else:
                logger.warning(
                    'missing attribute mentioned in :members: or __all__: '
                    'module %s, attribute %s' %
                    (safe_getattr(self.object, '__name__', '???'), mname))
        return False, ret


class ModuleLevelDocumenter(Documenter):
    """
    Specialized Documenter subclass for objects on module level (functions,
    classes, data/constants).
    """
    def resolve_name(self, modname, parents, path, base):
        # type: (str, Any, str, Any) -> Tuple[str, List[unicode]]
        if modname is None:
            if path:
                modname = path.rstrip('.')
            else:
                # if documenting a toplevel object without explicit module,
                # it can be contained in another auto directive ...
                modname = self.env.temp_data.get('autodoc:module')
                # ... or in the scope of a module directive
                if not modname:
                    modname = self.env.ref_context.get('py:module')
                # ... else, it stays None, which means invalid
        return modname, parents + [base]


class ClassLevelDocumenter(Documenter):
    """
    Specialized Documenter subclass for objects on class level (methods,
    attributes).
    """
    def resolve_name(self, modname, parents, path, base):
        # type: (str, Any, str, Any) -> Tuple[str, List[unicode]]
        if modname is None:
            if path:
                mod_cls = path.rstrip('.')
            else:
                mod_cls = None
                # if documenting a class-level object without path,
                # there must be a current class, either from a parent
                # auto directive ...
                mod_cls = self.env.temp_data.get('autodoc:class')
                # ... or from a class directive
                if mod_cls is None:
                    mod_cls = self.env.ref_context.get('py:class')
                # ... if still None, there's no way to know
                if mod_cls is None:
                    return None, []
            modname, _, cls = mod_cls.rpartition('.')  # type: ignore
            parents = [cls]
            # if the module name is still missing, get it like above
            if not modname:
                modname = self.env.temp_data.get('autodoc:module')
            if not modname:
                modname = self.env.ref_context.get('py:module')
            # ... else, it stays None, which means invalid
        return modname, parents + [base]


class DocstringSignatureMixin(object):
    """
    Mixin for FunctionDocumenter and MethodDocumenter to provide the
    feature of reading the signature from the docstring.
    """

    def _find_signature(self, encoding=None):
        # type: (unicode) -> Tuple[str, str]
        docstrings = self.get_doc(encoding)
        self._new_docstrings = docstrings[:]
        result = None
        for i, doclines in enumerate(docstrings):
            # no lines in docstring, no match
            if not doclines:
                continue
            # match first line of docstring against signature RE
            match = py_ext_sig_re.match(doclines[0])  # type: ignore
            if not match:
                continue
            exmod, path, base, args, retann = match.groups()
            # the base name must match ours
            valid_names = [self.objpath[-1]]  # type: ignore
            if isinstance(self, ClassDocumenter):
                valid_names.append('__init__')
                if hasattr(self.object, '__mro__'):
                    valid_names.extend(cls.__name__ for cls in self.object.__mro__)
            if base not in valid_names:
                continue
            # re-prepare docstring to ignore more leading indentation
            self._new_docstrings[i] = prepare_docstring('\n'.join(doclines[1:]))
            result = args, retann
            # don't look any further
            break
        return result

    def get_doc(self, encoding=None, ignore=1):
        # type: (unicode, int) -> List[List[unicode]]
        lines = getattr(self, '_new_docstrings', None)
        if lines is not None:
            return lines
        return Documenter.get_doc(self, encoding, ignore)  # type: ignore

    def format_signature(self):
        # type: () -> unicode
        if self.args is None and self.env.config.autodoc_docstring_signature:  # type: ignore
            # only act if a signature is not explicitly given already, and if
            # the feature is enabled
            result = self._find_signature()
            if result is not None:
                self.args, self.retann = result
        return Documenter.format_signature(self)


class DocstringStripSignatureMixin(DocstringSignatureMixin):
    """
    Mixin for AttributeDocumenter to provide the
    feature of stripping any function signature from the docstring.
    """
    def format_signature(self):
        # type: () -> unicode
        if self.args is None and self.env.config.autodoc_docstring_signature:
            # only act if a signature is not explicitly given already, and if
            # the feature is enabled
            result = self._find_signature()
            if result is not None:
                # Discarding _args is a only difference with
                # DocstringSignatureMixin.format_signature.
                # Documenter.format_signature use self.args value to format.
                _args, self.retann = result
        return Documenter.format_signature(self)


class FunctionDocumenter(DocstringSignatureMixin, ModuleLevelDocumenter):
    """
    Specialized Documenter subclass for functions.
    """
    objtype = 'function'
    member_order = 30

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        # It can be documented if it is a genuine function.
        # Often, a class instance has the same documentation as its class,
        # and then we typically want to document the class and not the instance.
        # However, there is an exception: CachedFunction(f) returns a class instance,
        # whose doc string coincides with that of f and is thus different from
        # that of the class CachedFunction. In that situation, we want that f is documented.
        # This is part of trac #9976.
        return (is_function_or_cython_function(member) or inspect.isbuiltin(member)
                or (isclassinstance(member)
                    and sage_getdoc_original(member) != sage_getdoc_original(member.__class__)))

    # Trac #9976: This function has been rewritten to support the
    # _sage_argspec_ attribute which makes it possible to get argument
    # specification of decorated callables in documentation correct.
    # See e.g. sage.misc.decorators.sage_wraps
    def args_on_obj(self, obj):
        if hasattr(obj, "_sage_argspec_"):
            return obj._sage_argspec_()
        if inspect.isbuiltin(obj) or \
               inspect.ismethoddescriptor(obj):
            # cannot introspect arguments of a C function or method
            # unless a function to do so is supplied
            if self.env.config.autodoc_builtin_argspec:
                argspec = self.env.config.autodoc_builtin_argspec(obj)
                return argspec
            else:
                return None
        argspec = sage_getargspec(obj)
        if isclassinstance(obj) or inspect.isclass(obj):
            # if a class should be documented as function, we try
            # to use the constructor signature as function
            # signature without the first argument.
            if argspec is not None and argspec[0]:
                del argspec[0][0]
        return argspec

    def format_args(self):
        # type: () -> unicode
        argspec = self.args_on_obj(self.object)
        if argspec is None:
            return None
        args = formatargspec(self.object, *argspec)
        # escape backslashes for reST
        args = args.replace('\\', '\\\\')
        return args

    def document_members(self, all_members=False):
        # type: (bool) -> None
        pass


class ClassDocumenter(DocstringSignatureMixin, ModuleLevelDocumenter):
    """
    Specialized Documenter subclass for classes.
    """
    objtype = 'class'
    member_order = 20
    option_spec = {
        'members': members_option, 'undoc-members': bool_option,
        'noindex': bool_option, 'inherited-members': bool_option,
        'show-inheritance': bool_option, 'member-order': identity,
        'exclude-members': members_set_option,
        'private-members': bool_option, 'special-members': members_option,
    }

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        return isinstance(member, type)

    def import_object(self):
        # type: () -> Any
        ret = ModuleLevelDocumenter.import_object(self)
        # Objective: if the class is documented under another name, document it
        # as data/attribute
        #
        # Notes from trac #27692:
        #
        # For Python 3, we make use of self.object.__qualname__.
        #
        # Notes from trac #7448:
        #
        # The original goal of this was that if some class is aliased, the
        # alias is generated as a link rather than duplicated. For example in
        #
        #     class A:
        #         pass
        #     B = A
        #
        # Then B is an alias of A, and should be generated as such.
        #
        # The way it was solved is to compare the name under which the
        # current class is found and the actual name of the class (stored in
        # the attribute __name__):
        #
        #   if hasattr(self.object, '__name__'):
        #       self.doc_as_attr = (self.objpath[-1] != self.object.__name__)
        #   else:
        #       self.doc_as_attr = True
        #
        # Now, to work around a pickling bug of nested class in Python,
        # by using the metaclass NestedMetaclass, we change the attribute
        # __name__ of the nested class. For example, in
        #
        #     class A(metaclass=NestedMetaclass):
        #        class B(object):
        #            pass
        #
        # the class B get its name changed to 'A.B'. Such dots '.' in names
        # are not supposed to occur in normal python name. I use it to check
        # if the class is a nested one and to compare its __name__ with its
        # path.
        #
        # The original implementation as well as the new one here doesn't work
        # if a class is aliased from a different place under the same name. For
        # example, in the following,
        #
        #     class A:
        #         pass
        #     class Container:
        #         A = A
        #
        # The nested copy Container.A is also documented. Actually, it seems
        # that there is no way to solve this by introspection. I'll submbit
        # this problem on sphinx trac.
        #
        # References: trac #5986, file sage/misc/nested_class.py
        if ret:
            module = getattr(self.object, '__module__', False)
            name = getattr(self.object, '__name__', False)
            qualname = getattr(self.object, '__qualname__', name)
            if qualname and module:
                # walk the standard attribute lookup path for this object
                qualname_parts = qualname.split('.')
                cls = getattr(sys.modules[module], qualname_parts[0], None)
                for part in qualname_parts[1:]:
                    if cls is None:
                        break
                    cls = getattr(cls, part, None)
                self.doc_as_attr = (self.objpath != qualname_parts and
                                    self.object is cls)
            else:
                self.doc_as_attr = True
        return ret

    def format_args(self):
        # type: () -> unicode
        # for classes, the relevant signature is the __init__ method's
        initmeth = self.get_attr(self.object, '__init__', None)
        # classes without __init__ method, default __init__ or
        # __init__ written in C?
        if initmeth is None or \
                is_builtin_class_method(self.object, '__init__') or \
                not(inspect.ismethod(initmeth) or is_function_or_cython_function(initmeth)):
            return None
        try:
            argspec = sage_getargspec(initmeth)
            if argspec[0] and argspec[0][0] in ('cls', 'self'):
                del argspec[0][0]
            return formatargspec(initmeth, *argspec)
        except TypeError:
            # still not possible: happens e.g. for old-style classes
            # with __init__ in C
            return None

    def format_signature(self):
        # type: () -> unicode
        if self.doc_as_attr:
            return ''

        return DocstringSignatureMixin.format_signature(self)

    def add_directive_header(self, sig):
        # type: (unicode) -> None
        if self.doc_as_attr:
            self.directivetype = 'attribute'
        Documenter.add_directive_header(self, sig)

        # add inheritance info, if wanted
        if not self.doc_as_attr and self.options.show_inheritance:
            sourcename = self.get_sourcename()
            self.add_line(u'', sourcename)
            if hasattr(self.object, '__bases__') and len(self.object.__bases__):
                bases = [b.__module__ in ('__builtin__', 'builtins') and
                         u':class:`%s`' % b.__name__ or
                         u':class:`%s.%s`' % (b.__module__, b.__name__)
                         for b in self.object.__bases__]
                self.add_line(u'   ' + _(u'Bases: %s') % ', '.join(bases),
                              sourcename)

    def get_doc(self, encoding=None, ignore=1):
        # type: (unicode, int) -> List[List[unicode]]
        lines = getattr(self, '_new_docstrings', None)
        if lines is not None:
            return lines

        content = self.env.config.autoclass_content

        docstrings = []
        attrdocstring = sage_getdoc_original(self.object)
        if attrdocstring:
            docstrings.append(attrdocstring)

        # for classes, what the "docstring" is can be controlled via a
        # config value; the default is only the class docstring
        if content in ('both', 'init'):
            initdocstring = sage_getdoc_original(
                self.get_attr(self.object, '__init__', None))
            # for new-style classes, no __init__ means default __init__
            if (initdocstring is not None and
                (initdocstring == object.__init__.__doc__ or  # for pypy
                 initdocstring.strip() == object.__init__.__doc__)):  # for !pypy
                initdocstring = None
            if not initdocstring:
                # try __new__
                initdocstring = self.get_attr(
                    self.get_attr(self.object, '__new__', None), '__doc__')
                # for new-style classes, no __new__ means default __new__
                if (initdocstring is not None and
                    (initdocstring == object.__new__.__doc__ or  # for pypy
                     initdocstring.strip() == object.__new__.__doc__)):  # for !pypy
                    initdocstring = None
            if initdocstring:
                if content == 'init':
                    docstrings = [initdocstring]
                else:
                    docstrings.append(initdocstring)
        doc = []
        for docstring in docstrings:
            if isinstance(docstring, str):
                doc.append(prepare_docstring(docstring))
        return doc

    def add_content(self, more_content, no_docstring=False):
        # type: (Any, bool) -> None
        if self.doc_as_attr:
            # We cannot rely on __qualname__ yet for Python 2, because of a
            # Cython bug: https://github.com/cython/cython/issues/2772. See
            # trac #27002.
            classname = safe_getattr(self.object, '__qualname__', None)
            if not classname:
                classname = safe_getattr(self.object, '__name__', None)
            if classname:
                module = safe_getattr(self.object, '__module__', None)
                parentmodule = safe_getattr(self.parent, '__module__', None)
                if module and module != parentmodule:
                    classname = str(module) + u'.' + str(classname)
                content = ViewList(
                    [_('alias of :class:`%s`') % classname], source='')
                ModuleLevelDocumenter.add_content(self, content,
                                                  no_docstring=True)
        else:
            ModuleLevelDocumenter.add_content(self, more_content)

    def document_members(self, all_members=False):
        # type: (bool) -> None
        if self.doc_as_attr:
            return
        ModuleLevelDocumenter.document_members(self, all_members)

    def generate(self, more_content=None, real_modname=None,
                 check_module=False, all_members=False):
        # Do not pass real_modname and use the name from the __module__
        # attribute of the class.
        # If a class gets imported into the module real_modname
        # the analyzer won't find the source of the class, if
        # it looks in real_modname.
        return super(ClassDocumenter, self).generate(more_content=more_content,
                                                     check_module=check_module,
                                                     all_members=all_members)

class ExceptionDocumenter(ClassDocumenter):
    """
    Specialized ClassDocumenter subclass for exceptions.
    """
    objtype = 'exception'
    member_order = 10

    # needs a higher priority than ClassDocumenter
    priority = 10

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        return isinstance(member, type) and issubclass(member, BaseException)


class DataDocumenter(ModuleLevelDocumenter):
    """
    Specialized Documenter subclass for data items.
    """
    objtype = 'data'
    member_order = 40
    priority = -10
    option_spec = dict(ModuleLevelDocumenter.option_spec)
    option_spec["annotation"] = annotation_option

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        return isinstance(parent, ModuleDocumenter) and isattr

    def add_directive_header(self, sig):
        # type: (unicode) -> None
        ModuleLevelDocumenter.add_directive_header(self, sig)
        sourcename = self.get_sourcename()
        if not self.options.annotation:
            try:
                objrepr = object_description(self.object)
            except ValueError:
                pass
            else:
                self.add_line(u'   :annotation: = ' + objrepr, sourcename)
        elif self.options.annotation is SUPPRESS:
            pass
        else:
            self.add_line(u'   :annotation: %s' % self.options.annotation,
                          sourcename)

    def document_members(self, all_members=False):
        # type: (bool) -> None
        pass


class MethodDocumenter(DocstringSignatureMixin, ClassLevelDocumenter):
    """
    Specialized Documenter subclass for methods (normal, static and class).
    """
    objtype = 'method'
    member_order = 50
    priority = 1  # must be more than FunctionDocumenter

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        return inspect.isroutine(member) and not \
                isinstance(parent, ModuleDocumenter)

    def import_object(self):
        # type: () -> Any
        ret = ClassLevelDocumenter.import_object(self)
        if not ret:
            return ret

        # to distinguish classmethod/staticmethod
        obj = self.parent.__dict__.get(self.object_name)
        if obj is None:
            obj = self.object

        if isclassmethod(obj):
            self.directivetype = 'classmethod'
            # document class and static members before ordinary ones
            self.member_order = self.member_order - 1
        elif isstaticmethod(obj, cls=self.parent, name=self.object_name):
            self.directivetype = 'staticmethod'
            # document class and static members before ordinary ones
            self.member_order = self.member_order - 1
        else:
            self.directivetype = 'method'
        return ret

    # Trac #9976: This function has been rewritten to support the
    # _sage_argspec_ attribute which makes it possible to get argument
    # specification of decorated callables in documentation correct.
    # See e.g. sage.misc.decorators.sage_wraps.
    #
    # Note, however, that sage.misc.sageinspect.sage_getargspec already
    # uses a method _sage_argspec_, that only works on objects, not on classes, though.
    def args_on_obj(self, obj):
        if hasattr(obj, "_sage_argspec_"):
            argspec = obj._sage_argspec_()
        elif inspect.isbuiltin(obj) or inspect.ismethoddescriptor(obj):
            # can never get arguments of a C function or method unless
            # a function to do so is supplied
            if self.env.config.autodoc_builtin_argspec:
                argspec = self.env.config.autodoc_builtin_argspec(obj)
            else:
                return None
        else:
            # The check above misses ordinary Python methods in Cython
            # files.
            argspec = sage_getargspec(obj)
        if argspec is not None and argspec[0] and argspec[0][0] in ('cls', 'self'):
            del argspec[0][0]
        return argspec

    def format_args(self):
        # type: () -> unicode
        argspec = self.args_on_obj(self.object)
        if argspec is None:
            return None
        args = formatargspec(self.object, *argspec)
        # escape backslashes for reST
        args = args.replace('\\', '\\\\')
        return args

    def document_members(self, all_members=False):
        # type: (bool) -> None
        pass


class AttributeDocumenter(DocstringStripSignatureMixin, ClassLevelDocumenter):  # type: ignore
    """
    Specialized Documenter subclass for attributes.
    """
    objtype = 'attribute'
    member_order = 60
    option_spec = dict(ModuleLevelDocumenter.option_spec)
    option_spec["annotation"] = annotation_option

    # must be higher than the MethodDocumenter, else it will recognize
    # some non-data descriptors as methods
    priority = 10

    @staticmethod
    def is_function_or_method(obj):
        return is_function_or_cython_function(obj) or inspect.isbuiltin(obj) or inspect.ismethod(obj)

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        non_attr_types = (type, MethodDescriptorType)
        isdatadesc = isdescriptor(member) and not \
                inspect.isroutine(member) and not \
                isinstance(member, non_attr_types) and not \
                type(member).__name__ == "instancemethod"

        isattribute = isdatadesc or (not isinstance(parent, ModuleDocumenter) and isattr)

        # Trac #26522: This condition is here just to pass objects of classes
        # that inherit ClasscallMetaclass as attributes rather than method
        # descriptors.
        isattribute = isattribute or isinstance(type(member), ClasscallMetaclass)

        return isattribute

        # We ignore the obscure case supported in the following return
        # statement. The additional check opens a door for attributes without
        # docstrings to appear in the Sage documentation, and more seriously
        # effectively prevents certain attributes to get properly documented.
        # See trac #28698.

        # That last condition addresses an obscure case of C-defined
        # methods using a deprecated type in Python 3, that is not otherwise
        # exported anywhere by Python
        return isattribute or (not isinstance(parent, ModuleDocumenter) and
                              not inspect.isroutine(member) and
                              not isinstance(member, type))

    def document_members(self, all_members=False):
        # type: (bool) -> None
        pass

    def import_object(self):
        # type: () -> Any
        ret = ClassLevelDocumenter.import_object(self)
        if isenumattribute(self.object):
            self.object = self.object.value
        if isdescriptor(self.object) and \
                not self.is_function_or_method(self.object):
            self._datadescriptor = True
        else:
            # if it's not a data descriptor
            self._datadescriptor = False
        return ret

    def get_real_modname(self):
        # type: () -> str
        return self.get_attr(self.parent or self.object, '__module__', None) \
            or self.modname

    def add_directive_header(self, sig):
        # type: (unicode) -> None
        ClassLevelDocumenter.add_directive_header(self, sig)
        sourcename = self.get_sourcename()
        if not self.options.annotation:
            if not self._datadescriptor:
                try:
                    objrepr = object_description(self.object)
                except ValueError:
                    pass
                else:
                    self.add_line(u'   :annotation: = ' + objrepr, sourcename)
        elif self.options.annotation is SUPPRESS:
            pass
        else:
            self.add_line(u'   :annotation: %s' % self.options.annotation,
                          sourcename)

    def add_content(self, more_content, no_docstring=False):
        # type: (Any, bool) -> None
        if not self._datadescriptor:
            # if it's not a data descriptor, its docstring is very probably the
            # wrong thing to display
            no_docstring = True
        ClassLevelDocumenter.add_content(self, more_content, no_docstring)


class InstanceAttributeDocumenter(AttributeDocumenter):
    """
    Specialized Documenter subclass for attributes that cannot be imported
    because they are instance attributes (e.g. assigned in __init__).
    """
    objtype = 'instanceattribute'
    directivetype = 'attribute'
    member_order = 60

    # must be higher than AttributeDocumenter
    priority = 11

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        # type: (Any, unicode, bool, Any) -> bool
        """This documents only INSTANCEATTR members."""
        return isattr and (member is INSTANCEATTR)

    def import_object(self):
        # type: () -> bool
        """Never import anything."""
        # disguise as an attribute
        self.objtype = 'attribute'
        self._datadescriptor = False
        return True

    def add_content(self, more_content, no_docstring=False):
        # type: (Any, bool) -> None
        """Never try to get a docstring from the object."""
        AttributeDocumenter.add_content(self, more_content, no_docstring=True)


def get_documenters(app):
    # type: (Sphinx) -> Dict[unicode, Type[Documenter]]
    """Returns registered Documenter classes"""
    return app.registry.documenters


def autodoc_attrgetter(app, obj, name, *defargs):
    # type: (Sphinx, Any, unicode, Any) -> Any
    """Alternative getattr() for types"""
    for typ, func in app.registry.autodoc_attrgettrs.items():
        if isinstance(obj, typ):
            return func(obj, name, *defargs)

    return safe_getattr(obj, name, *defargs)


def setup(app):
    # type: (Sphinx) -> Dict[unicode, Any]
    app.add_autodocumenter(ModuleDocumenter)
    app.add_autodocumenter(ClassDocumenter)
    app.add_autodocumenter(ExceptionDocumenter)
    app.add_autodocumenter(DataDocumenter)
    app.add_autodocumenter(FunctionDocumenter)
    app.add_autodocumenter(MethodDocumenter)
    app.add_autodocumenter(AttributeDocumenter)
    app.add_autodocumenter(InstanceAttributeDocumenter)

    app.add_config_value('autoclass_content', 'class', True)
    app.add_config_value('autodoc_member_order', 'alphabetic', True)
    app.add_config_value('autodoc_default_flags', [], True) # deprecated since Sphinx 1.8
    app.add_config_value('autodoc_default_options', {}, True)
    app.add_config_value('autodoc_docstring_signature', False, True)
    app.add_config_value('autodoc_mock_imports', [], True)
    app.add_config_value('autodoc_warningiserror', True, True)
    app.add_config_value('autodoc_inherit_docstrings', True, True)
    app.add_config_value('autodoc_builtin_argspec', None, True)

    app.add_event('autodoc-process-docstring')
    app.add_event('autodoc-process-signature')
    app.add_event('autodoc-skip-member')

    return {'version': sphinx.__display_version__, 'parallel_read_safe': True}


class testcls:
    """test doc string"""

    def __getattr__(self, x):
        # type: (Any) -> Any
        return x

    def __setattr__(self, x, y):
        # type: (Any, Any) -> None
        """Attr setter."""
