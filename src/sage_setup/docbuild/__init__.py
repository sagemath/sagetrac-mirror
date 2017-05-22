from __future__ import absolute_import
from __future__ import print_function
from six.moves import range

import logging, optparse, os, shutil, subprocess, sys, re

import sphinx.cmdline
import sphinx.util.console
import sphinx.ext.intersphinx

from .build_options import (LANGUAGES, SPHINXOPTS, PAPER, OMIT,
     PAPEROPTS, ALLSPHINXOPTS, NUM_THREADS, WEBSITESPHINXOPTS,
     INCREMENTAL_BUILD, ABORT_ON_ERROR)


def get_builder(name):
    """
    Returns an appropriate *Builder object for the document ``name``.
    DocBuilder and its subclasses do all the real work in building the
    documentation.
    """
    if name == 'all':
        return AllBuilder()
    elif name.endswith('reference'):
        return ReferenceBuilder(name)
    elif 'reference' in name and os.path.exists(os.path.join(SAGE_DOC_SRC, 'en', name)):
        return ReferenceSubBuilder(name)
    elif name.endswith('website'):
        return WebsiteBuilder(name)
    elif name.startswith('file='):
        path = name[5:]
        if path.endswith('.sage') or path.endswith('.pyx'):
            raise NotImplementedError('Building documentation for a single file only works for Python files.')
        return SingleFileBuilder(path)
    elif name in get_documents() or name in AllBuilder().get_all_documents():
        return DocBuilder(name)
    else:
        print("'%s' is not a recognized document. Type 'sage --docbuild -D' for a list"%name)
        print("of documents, or 'sage --docbuild --help' for more help.")
        sys.exit(1)


def get_documents():
    """
    Returns a list of document names the Sage documentation builder
    will accept as command-line arguments.
    """
    all_b = AllBuilder()
    docs = all_b.get_all_documents()
    docs = [(d[3:] if d[0:3] == 'en/' else d) for d in docs]
    return docs


def get_formats():
    """
    Returns a list of output formats the Sage documentation builder
    will accept on the command-line.
    """
    tut_b = DocBuilder('en/tutorial')
    formats = tut_b._output_formats()
    formats.remove('html')
    return ['html', 'pdf'] + formats
