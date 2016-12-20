"""Miscellaneous utilities for the documentation builder."""

from __future__ import absolute_import
from __future__ import print_function

import logging
import math
import optparse
import os
import sys

from six.moves import range

from . import get_documents, get_formats


logger = logging.getLogger('docbuild')


# File utilities

def delete_empty_directories(root_dir):
    """
    Delete all empty directories found under ``root_dir``

    INPUT:

    - ``root_dir`` -- string. A valid directory name.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir, topdown=False):
        if not dirnames + filenames:
            logger.warning('Deleting empty directory {0}'.format(dirpath))
            os.rmdir(dirpath)


# Help formatting utilities for the docbuild main script

class IndentedHelpFormatter2(optparse.IndentedHelpFormatter, object):
    """
    Custom help formatter class for optparse's OptionParser.
    """
    def format_description(self, description):
        """
        Returns a formatted description, preserving any original
        explicit new line characters.
        """
        if description:
            lines_in = description.split('\n')
            lines_out = [self._format_text(line) for line in lines_in]
            return "\n".join(lines_out) + "\n"
        else:
            return ""

    def format_heading(self, heading):
        """
        Returns a formatted heading using the superclass' formatter.
        If the heading is 'options', up to case, the function converts
        it to ALL CAPS.  This allows us to match the heading 'OPTIONS' with
        the same token in the builder's usage message.
        """
        if heading.lower() == 'options':
            heading = "OPTIONS"
        return super(IndentedHelpFormatter2, self).format_heading(heading)


def help_usage(s=u"", compact=False):
    """
    Appends and returns a brief usage message for the Sage
    documentation builder.  If 'compact' is False, the function adds a
    final newline character.
    """
    s += "sage --docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)"
    if not compact:
        s += "\n"
    return s


def help_description(s=u"", compact=False):
    """
    Appends and returns a brief description of the Sage documentation
    builder.  If 'compact' is False, the function adds a final newline
    character.
    """
    s += "Build or return information about Sage documentation.\n\n"
    s += "    DOCUMENT    name of the document to build\n"
    s += "    FORMAT      document output format\n"
    s += "    COMMAND     document-specific command\n\n"
    s += "Note that DOCUMENT may have the form 'file=/path/to/FILE',\n"
    s += "which builds the documentation for the specified file.\n\n"
    s += "A DOCUMENT and either a FORMAT or a COMMAND are required,\n"
    s += "unless a list of one or more of these is requested."
    if not compact:
        s += "\n"
    return s


def help_documents(s=u""):
    """
    Appends and returns a tabular list of documents, including a
    shortcut 'all' for all documents, available to the Sage
    documentation builder.
    """
    docs = get_documents(default_lang='en')
    s += "DOCUMENTs:\n"
    s += format_columns(docs + ['all  (!)'])
    s += "(!) Builds everything.\n\n"
    if 'reference' in docs:
        s+= "Other valid document names take the form 'reference/DIR', where\n"
        s+= "DIR is a subdirectory of SAGE_DOC_SRC/en/reference/.\n"
        s+= "This builds just the specified part of the reference manual.\n"
    s += "DOCUMENT may also have the form 'file=/path/to/FILE', which builds\n"
    s += "the documentation for the specified file.\n"
    return s


def help_formats(s=u""):
    """
    Appends and returns a tabular list of output formats available to
    the Sage documentation builder.
    """
    s += "FORMATs:\n"
    s += format_columns(get_formats())
    return s


def help_commands(name='all', s=u""):
    """
    Appends and returns a tabular list of commands, if any, the Sage
    documentation builder can run on the indicated document.  The
    default is to list all commands for all documents.
    """
    # To do: Generate the lists dynamically, using class attributes,
    # as with the Builders above.
    command_dict = { 'reference' : [
        'print_included_modules',   'print_modified_modules       (*)',
        'print_unincluded_modules', 'print_newly_included_modules (*)',
        ] }
    for doc in command_dict:
        if name == 'all' or doc == name:
            s += "COMMANDs for the DOCUMENT '" + doc + "':\n"
            s += format_columns(command_dict[doc])
            s += "(*) Since the last build.\n"
    return s


def help_examples(s=u""):
    """
    Appends and returns some usage examples for the Sage documentation
    builder.
    """
    s += "Examples:\n"
    s += "    sage --docbuild -FDC all\n"
    s += "    sage --docbuild constructions pdf\n"
    s += "    sage --docbuild reference html -jv3\n"
    s += "    sage --docbuild --mathjax tutorial html\n"
    s += "    sage --docbuild reference print_unincluded_modules\n"
    s += "    sage --docbuild developer -j html --sphinx-opts -q,-aE --verbose 2"
    return s


def help_wrapper(option, opt_str, value, parser):
    """
    A helper wrapper for command-line options to the Sage
    documentation builder that print lists, such as document names,
    formats, and document-specific commands.
    """
    if option.dest == 'commands':
        print(help_commands(value), end="")
    if option.dest == 'documents':
        print(help_documents(), end="")
    if option.dest == 'formats':
        print(help_formats(), end="")
    setattr(parser.values, 'printed_list', 1)


def help_message_short(option=None, opt_str=None, value=None, parser=None,
                       error=False):
    """
    Prints a help message for the Sage documentation builder.  The
    message includes command-line usage and a list of options.  The
    message is printed only on the first call.  If error is True
    during this call, the message is printed only if the user hasn't
    requested a list (e.g., documents, formats, commands).
    """
    if not hasattr(parser.values, 'printed_help'):
        if error == True:
            if not hasattr(parser.values, 'printed_list'):
                parser.print_help()
        else:
            parser.print_help()
        setattr(parser.values, 'printed_help', 1)


def help_message_long(option, opt_str, value, parser):
    """
    Prints an extended help message for the Sage documentation builder
    and exits.
    """
    help_funcs = [ help_usage, help_description, help_documents,
                   help_formats, help_commands, parser.format_option_help,
                   help_examples ]
    for f in help_funcs:
        print(f())
    sys.exit(0)


def format_columns(lst, align=u'<', cols=None, indent=4, pad=3, width=80):
    """
    Utility function that formats a list as a simple table and returns
    a Unicode string representation.  The number of columns is
    computed from the other options, unless it's passed as a keyword
    argument.  For help on Python's string formatter, see

    http://docs.python.org/library/string.html#format-string-syntax
    """
    # Can we generalize this (efficiently) to other / multiple inputs
    # and generators?
    size = max(map(len, lst)) + pad
    if cols is None:
        cols = math.trunc((width - indent) / size)
    s = u" " * indent
    for i in range(len(lst)):
        if i != 0 and i % cols == 0:
            s += u"\n" + u" " * indent
        s += u"{0:{1}{2}}".format(lst[i], align, size)
    s += u"\n"
    return s


# Monkey-patches into Sphinx

class IntersphinxCache(object):
    """
    Replace sphinx.ext.intersphinx.fetch_inventory by an in-memory
    cached version.
    """
    def __init__(self):
        import sphinx.ext.intersphinx
        self.inventories = {}
        self.real_fetch_inventory = sphinx.ext.intersphinx.fetch_inventory
        sphinx.ext.intersphinx.fetch_inventory = self.fetch_inventory

    def fetch_inventory(self, app, uri, inv):
        """
        Return the result of ``sphinx.ext.intersphinx.fetch_inventory()``
        from a cache if possible. Otherwise, call
        ``sphinx.ext.intersphinx.fetch_inventory()`` and cache the result.
        """
        t = (uri, inv)
        try:
            return self.inventories[t]
        except KeyError:
            i = self.real_fetch_inventory(app, uri, inv)
            self.inventories[t] = i
            return i
