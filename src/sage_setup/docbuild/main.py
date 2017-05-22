"""
The documentation builder

It is the main entry point for building the documentation, and is responsible
for figuring out what to build and with which options. The actual documentation
build for each individual document is then done in a subprocess call to sphinx,
see :func:`output_formatter`.

* The builder can be configured in build_options.py
* The sphinx subprocesses are configured in conf.py
"""

from __future__ import absolute_import

import logging
import optparse
import os
import sys

from six.moves import range

from sage.env import SAGE_DOC_SRC

from . import build_options, get_builder
from .utils import (help_usage, IndentedHelpFormatter2, help_description,
                    help_message_short, help_message_long, help_wrapper,
                    IntersphinxCache, delete_empty_directories)


logger = logging.getLogger('docbuild')


def excepthook(*exc_info):
    """
    Print docbuild error and hint how to solve it
    """
    logger.error('Error building the documentation.', exc_info=exc_info)
    if build_options.INCREMENTAL_BUILD:
        logger.error('''
Note: incremental documentation builds sometimes cause spurious
error messages. To be certain that these are real errors, run
"make doc-clean" first and try again.''')


def setup_parser():
    """
    Sets up and returns a command-line OptionParser instance for the
    Sage documentation builder.
    """
    # Documentation: http://docs.python.org/library/optparse.html
    parser = optparse.OptionParser(add_help_option=False,
                                   usage=help_usage(compact=True),
                                   formatter=IndentedHelpFormatter2(),
                                   description=help_description(compact=True))

    # Standard options. Note: We use explicit option.dest names
    # to avoid ambiguity.
    standard = optparse.OptionGroup(parser, "Standard")
    standard.add_option("-h", "--help",
                        action="callback", callback=help_message_short,
                        help="show a help message and exit")
    standard.add_option("-H", "--help-all",
                        action="callback", callback=help_message_long,
                        help="show an extended help message and exit")
    standard.add_option("-D", "--documents", dest="documents",
                        action="callback", callback=help_wrapper,
                        help="list all available DOCUMENTs")
    standard.add_option("-F", "--formats", dest="formats",
                        action="callback", callback=help_wrapper,
                        help="list all output FORMATs")
    standard.add_option("-C", "--commands", dest="commands",
                        type="string", metavar="DOC",
                        action="callback", callback=help_wrapper,
                        help="list all COMMANDs for DOCUMENT DOC; use 'all' "
                             "to list all")

    standard.add_option("-i", "--inherited", dest="inherited",
                        default=False, action="store_true",
                        help="include inherited members in reference manual; "
                             "may be slow, may fail for PDF output")
    standard.add_option("-u", "--underscore", dest="underscore",
                        default=False, action="store_true",
                        help="include variables prefixed with '_' in "
                             "reference manual; may be slow, may fail for "
                             "PDF output")

    standard.add_option("-j", "--mathjax", "--jsmath", dest="mathjax",
                        action="store_true",
                        help="render math using MathJax; FORMATs: html, "
                             "json, pickle, web")
    standard.add_option("--no-plot", dest="no_plot",
                        action="store_true",
                        help="do not include graphics auto-generated using "
                             "the '.. plot' markup")
    standard.add_option("--include-tests-blocks", dest="skip_tests", default=True,
                        action="store_false",
                        help="include TESTS blocks in the reference manual")
    standard.add_option("--no-pdf-links", dest="no_pdf_links",
                        action="store_true",
                        help="do not include PDF links in DOCUMENT "
                             "'website'; FORMATs: html, json, pickle, web")
    standard.add_option("--warn-links", dest="warn_links",
                        default=False, action="store_true",
                        help="issue a warning whenever a link is not "
                             "properly resolved; equivalent to "
                             "'--sphinx-opts -n' (sphinx option: nitpicky)")
    standard.add_option("--check-nested", dest="check_nested",
                        action="store_true",
                        help="check picklability of nested classes in "
                             "DOCUMENT 'reference'")
    standard.add_option("-N", "--no-colors", dest="color", default=True,
                        action="store_false",
                        help="do not color output; does not affect "
                             "children")
    standard.add_option("-q", "--quiet", dest="verbose",
                        action="store_const", const=0,
                        help="work quietly; same as --verbose=0")
    standard.add_option("-v", "--verbose", dest="verbose",
                        type="int", default=1, metavar="LEVEL",
                        action="store",
                        help="report progress at LEVEL=0 (quiet), 1 (normal), "
                             "2 (info), or 3 (debug); does not "
                             "affect children")
    standard.add_option("-o", "--output", dest="output_dir", default=None,
                        metavar="DIR", action="store",
                        help="if DOCUMENT is a single file ('file=...'), "
                             "write output to this directory")
    parser.add_option_group(standard)

    # Advanced options.
    advanced = optparse.OptionGroup(parser, "Advanced",
                                    "Use these options with care.")
    advanced.add_option("-S", "--sphinx-opts", dest="sphinx_opts",
                        type="string", metavar="OPTS",
                        action="store",
                        help="pass comma-separated OPTS to sphinx-build")
    advanced.add_option("-U", "--update-mtimes", dest="update_mtimes",
                        default=False, action="store_true",
                        help="before building reference manual, update "
                             "modification times for auto-generated ReST "
                             "files")
    advanced.add_option("-k", "--keep-going", dest="keep_going",
                        default=False, action="store_true",
                        help="Do not abort on errors but continue as "
                             "much as possible after an error")
    parser.add_option_group(advanced)

    return parser


def setup_logger(name='docbuild', verbose=1, color=True):
    """
    Sets up the docbuild logger for use by the command-line application.
    The optional arguments set the logger's level and message formats.
    """

    logger = logging.getLogger(name)

    # Set up colors. Adapted from sphinx.cmdline.
    import sphinx.util.console as c
    if not color or not sys.stdout.isatty() or not c.color_terminal():
        c.nocolor()

    # Available colors: black, darkgray, (dark)red, dark(green),
    # brown, yellow, (dark)blue, purple, fuchsia, turquoise, teal,
    # lightgray, white.  Available styles: reset, bold, faint,
    # standout, underline, blink.

    # Set up log record formats.
    format_std = "%(message)s"
    formatter = logging.Formatter(format_std)

    # format_debug = "%(module)s #%(lineno)s %(funcName)s() %(message)s"
    fields = ['%(module)s', '#%(lineno)s', '%(funcName)s()', '%(message)s']
    colors = ['darkblue', 'darkred', 'brown', 'reset']
    styles = ['reset', 'reset', 'reset', 'reset']
    format_debug = ""
    for i in range(len(fields)):
        format_debug += c.colorize(styles[i], c.colorize(colors[i],
                                                         fields[i]))
        if i != len(fields):
            format_debug += " "

    # Note: There's also Handler.setLevel().  The argument is the
    # lowest severity message that the respective logger or handler
    # will pass on.  The default levels are DEBUG, INFO, WARNING,
    # ERROR, and CRITICAL.  We use "WARNING" for normal verbosity and
    # "ERROR" for quiet operation.  It's possible to define custom
    # levels.  See the documentation for details.
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    if verbose == 1:
        logger.setLevel(logging.WARNING)
    if verbose == 2:
        logger.setLevel(logging.INFO)
    if verbose == 3:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(format_debug)

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)


def main():
    # Parse the command-line.
    parser = setup_parser()
    options, args = parser.parse_args()

    # Get the name and type (target format) of the document we are
    # trying to build.
    try:
        name, type = args
    except ValueError:
        help_message_short(parser=parser, error=True)
        sys.exit(1)

    # Set up module-wide logging.
    setup_logger('docbuild', options.verbose, options.color)

    sys.excepthook = excepthook

    # Process selected options.
    #
    # MathJax: this check usually has no practical effect, since
    # SAGE_DOC_MATHJAX is set to "True" by the script sage-env.
    # To disable MathJax, set SAGE_DOC_MATHJAX to "no" or "False".
    if options.mathjax or (os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'no'
                           and os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'False'):
        os.environ['SAGE_DOC_MATHJAX'] = 'True'

    if options.check_nested:
        os.environ['SAGE_CHECK_NESTED'] = 'True'

    if options.inherited:
        build_options.INHERITED = True

    if options.underscore:
        os.environ['SAGE_DOC_UNDERSCORE'] = "True"
        build_options.UNDERSCORE = True

    if options.update_mtimes:
        build_options.UPDATE_MTIMES = True

    if options.sphinx_opts:
        build_options.ALLSPHINXOPTS += options.sphinx_opts.replace(',', ' ') + " "
    if options.no_pdf_links:
        build_options.WEBSITESPHINXOPTS = " -A hide_pdf_links=1 "
    if options.warn_links:
        build_options.ALLSPHINXOPTS += "-n "
    if options.no_plot:
        os.environ['SAGE_SKIP_PLOT_DIRECTIVE'] = 'yes'
    if options.skip_tests:
        os.environ['SAGE_SKIP_TESTS_BLOCKS'] = 'True'

    if options.output_dir:
        build_options.OUTPUT_DIR = options.output_dir

    build_options.ABORT_ON_ERROR = not options.keep_going

    # Delete empty directories. This is needed in particular for empty
    # directories due to "git checkout" which never deletes empty
    # directories it leaves behind. See Trac #20010.
    delete_empty_directories(SAGE_DOC_SRC)

    # Set up Intersphinx cache
    C = IntersphinxCache()

    builder = getattr(get_builder(name), type)
    builder()
