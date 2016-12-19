from __future__ import absolute_import
from __future__ import print_function

import logging
import os
import re
import shutil
import sys
import time
import types

from functools import partial, update_wrapper

from sage.env import SAGE_DOC, SAGE_DOC_SRC
from sage.misc.cachefunc import cached_method
from sage.misc.misc import sage_makedirs

from six.moves import cPickle

from sphinx.environment import BuildEnvironment

from . import AllBuilder, build_many, get_builder
from .docbuilder import DocBuilder
from .. import build_options as opts


__all__ = ['ReferenceBuilder', 'ReferenceSubBuilder']


logger = logging.getLogger('docbuild')


class ReferenceBuilder(AllBuilder):
    """
    This class builds the reference manual.  It uses DocBuilder to
    build the top-level page and ReferenceSubBuilder for each
    sub-component.
    """
    def __init__(self, name, lang='en'):
        """
        Records the reference manual's name, in case it's not
        identical to 'reference'.
        """
        super(ReferenceBuilder, self).__init__()
        doc = name.split(os.path.sep)

        if doc[0] in opts.LANGUAGES:
            lang = doc[0]
            doc.pop(0)

        self.name = doc[0]
        self.lang = lang

    def _output_dir(self, type, lang='en'):
        """
        Returns the directory where the output of type type is stored.
        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceBuilder
            sage: b = ReferenceBuilder('reference')
            sage: b._output_dir('html')
            '.../html/en/reference'
        """
        d = os.path.join(SAGE_DOC, type, lang, self.name)
        sage_makedirs(d)
        return d

    @staticmethod
    def _build_ref_doc(args):
        """Target for parallel reference doc builds."""

        doc = args[0]
        lang = args[1]
        format = args[2]
        kwds = args[3]
        args = args[4:]
        if format == 'inventory':
            # you must not use the inventory to build the inventory
            kwds['use_multidoc_inventory'] = False
        getattr(ReferenceSubBuilder(doc, lang), format)(*args, **kwds)

    def _wrapper(self, format, *args, **kwds):
        """
        Builds reference manuals.  For each language, it builds the
        top-level document and its components.
        """
        for lang in opts.LANGUAGES:
            refdir = os.path.join(SAGE_DOC_SRC, lang, self.name)
            if not os.path.exists(refdir):
                continue
            output_dir = self._output_dir(format, lang)
            L = [(doc, lang, format, kwds) + args for doc in self.get_all_documents(refdir)]
            build_many(self._build_ref_doc, L)
            # The html refman must be build at the end to ensure correct
            # merging of indexes and inventories.
            # Sphinx is run here in the current process (not in a
            # subprocess) and the IntersphinxCache gets populated to be
            # used for the second pass of the reference manual and for
            # the other documents.
            getattr(DocBuilder(self.name, lang), format)(*args, **kwds)

            # PDF: we need to build master index file which lists all
            # of the PDF file.  So we create an html file, based on
            # the file index.html from the "website" target.
            if format == 'pdf':
                # First build the website page.  (This only takes a
                # few seconds.)
                getattr(get_builder('website'), 'html')()
                # Copy the relevant pieces of
                # SAGE_DOC/html/en/website/_static to output_dir.
                # (Don't copy all of _static to save some space: we
                # don't need all of the MathJax stuff, and in
                # particular we don't need the fonts.)
                website_dir = os.path.join(SAGE_DOC, 'html',
                                           'en', 'website')
                static_files = ['COPYING.txt', 'basic.css', 'blank.gif',
                         'default.css', 'doctools.js', 'favicon.ico',
                         'file.png', 'jquery.js', 'minus.png',
                         'pdf.png', 'plus.png', 'pygments.css',
                         'sage.css', 'sageicon.png',
                         'logo_sagemath.svg', 'logo_sagemath_black.svg'
                         'searchtools.js', 'sidebar.js', 'underscore.js']
                sage_makedirs(os.path.join(output_dir, '_static'))
                for f in static_files:
                    try:
                        shutil.copyfile(os.path.join(website_dir, '_static', f),
                                        os.path.join(output_dir, '_static', f))
                    except IOError: # original file does not exist
                        pass
                # Now modify website's index.html page and write it
                # to output_dir.
                f = open(os.path.join(website_dir, 'index.html'))
                html = f.read().replace('Documentation', 'Reference')
                f.close()
                html_output_dir = os.path.dirname(website_dir)
                html = html.replace('http://www.sagemath.org',
                                    os.path.join(html_output_dir, 'index.html'))
                # From index.html, we want the preamble and the tail.
                html_end_preamble = html.find('<h1>Sage Reference')
                html_bottom = html.rfind('</table>') + len('</table>')
                # For the content, we modify doc/en/reference/index.rst,
                # which has two parts: the body and the table of contents.
                f = open(os.path.join(SAGE_DOC_SRC, lang, 'reference', 'index.rst'))
                rst = f.read()
                f.close()
                # Replace rst links with html links.  There are two forms:
                #
                #   `blah`__    followed by __ LINK
                #
                #   :doc:`blah <module/index>`
                #
                # Change the first form to
                #
                #   <a href="LINK">blah</a>
                #
                # Change the second form to
                #
                #   <a href="module/module.pdf">blah <img src="_static/pdf.png" /></a>
                #
                rst = re.sub('`([^`]*)`__\.\n\n__ (.*)',
                                  r'<a href="\2">\1</a>.', rst)
                rst = re.sub(r':doc:`([^<]*?)\s+<(.*)/index>`',
                             r'<a href="\2/\2.pdf">\1 <img src="_static/pdf.png" /></a>',
                             rst)
                # Get rid of todolist and miscellaneous rst markup.
                rst = rst.replace('.. toctree::', '')
                rst = rst.replace(':maxdepth: 2', '')
                rst = rst.replace('todolist', '')
                start = rst.find('=\n') + 1
                end = rst.find('Table of Contents')
                # Body: add paragraph <p> markup.
                rst_body = rst[start:end]
                rst_body = rst_body.replace('\n\n', '</p>\n<p>')
                start = rst.find('Table of Contents') + 2*len('Table of Contents') + 1
                # Don't include the indices.
                end = rst.find('Indices and Tables')
                # TOC: change * to <li>, change rst headers to html headers.
                rst_toc = rst[start:end]
                rst_toc = rst_toc.replace('*', '<li>')
                rst_toc = re.sub('\n([A-Z][a-zA-Z, ]*)\n-*\n',
                             '</ul>\n\n\n<h2>\\1</h2>\n\n<ul>\n', rst_toc)
                # Now write the file.
                new_index = open(os.path.join(output_dir, 'index.html'), 'w')
                new_index.write(html[:html_end_preamble])
                new_index.write('<h1>' + rst[:rst.find('\n')] +
                                ' (PDF version)'+ '</h1>')
                new_index.write(rst_body)
                new_index.write('<h2>Table of Contents</h2>\n\n<ul>')
                new_index.write(rst_toc)
                new_index.write('</ul>\n\n')
                new_index.write(html[html_bottom:])
                new_index.close()
                logger.warning('''
PDF documents have been created in subdirectories of

  %s

Alternatively, you can open

  %s

for a webpage listing all of the documents.''' % (output_dir,
                                                 os.path.join(output_dir,
                                                              'index.html')))

    def get_all_documents(self, refdir):
        """
        Returns a list of all reference manual components to build.
        We add a component name if it's a subdirectory of the manual's
        directory and contains a file named 'index.rst'.

        We return the largest component (most subdirectory entries)
        first since they will take the longest to build.

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceBuilder
            sage: b = ReferenceBuilder('reference')
            sage: refdir = os.path.join(os.environ['SAGE_DOC_SRC'], 'en', b.name)
            sage: sorted(b.get_all_documents(refdir))
            ['reference/algebras',
             'reference/arithgroup',
             ...,
             'reference/tensor_free_modules']
        """
        documents = []

        for doc in os.listdir(refdir):
            directory = os.path.join(refdir, doc)
            if os.path.exists(os.path.join(directory, 'index.rst')):
                n = len(os.listdir(directory))
                documents.append((-n, os.path.join(self.name, doc)))

        return [ doc[1] for doc in sorted(documents) ]


class ReferenceSubBuilder(DocBuilder):
    """
    This class builds sub-components of the reference manual.  It is
    resposible for making sure the auto generated ReST files for the
    Sage library are up to date.

    When building any output, we must first go through and check
    to see if we need to update any of the autogenerated ReST
    files.  There are two cases where this would happen:

    1. A new module gets added to one of the toctrees.

    2. The actual module gets updated and possibly contains a new
       title.
    """
    def __init__(self, *args, **kwds):
        super(ReferenceSubBuilder, self).__init__(*args, **kwds)
        self._wrap_builder_helpers()

    def _wrap_builder_helpers(self):
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                f = partial(self._wrapper, attr)
                f.is_output_format = True
                update_wrapper(f, getattr(self, attr))
                setattr(self, attr, f)

    def _wrapper(self, build_type, *args, **kwds):
        """
        This is the wrapper around the builder_helper methods that
        goes through and makes sure things are up to date.
        """
        # Delete the auto-generated .rst files, if the inherited
        # and/or underscored members options have changed.
        inherit_prev = self.get_cache().get('option_inherited')
        underscore_prev = self.get_cache().get('option_underscore')
        if (inherit_prev is None or inherit_prev != opts.INHERITED or
            underscore_prev is None or underscore_prev != opts.UNDERSCORE):
            logger.info("Detected change(s) in inherited and/or underscored members option(s).")
            self.clean_auto()
            self.get_cache.clear_cache()

        # After "sage -clone", refresh the .rst file mtimes in
        # environment.pickle.
        if opts.UPDATE_MTIMES:
            logger.info("Checking for .rst file mtimes to update...")
            self.update_mtimes()

        #Update the .rst files for modified Python modules
        logger.info("Updating .rst files with modified modules...")
        for module_name in self.get_modified_modules():
            self.write_auto_rest_file(module_name.replace(os.path.sep, '.'))

        #Write the .rst files for newly included modules
        logger.info("Writing .rst files for newly-included modules...")
        for module_name in self.get_newly_included_modules(save=True):
            self.write_auto_rest_file(module_name)

        #Copy over the custom .rst files from _sage
        _sage = os.path.join(self.dir, '_sage')
        if os.path.exists(_sage):
            logger.info("Copying over custom .rst files from %s ...", _sage)
            shutil.copytree(_sage, os.path.join(self.dir, 'sage'))

        getattr(DocBuilder, build_type)(self, *args, **kwds)

    def cache_filename(self):
        """
        Returns the filename where the pickle of the dictionary of
        already generated ReST files is stored.
        """
        return os.path.join(self._doctrees_dir(), 'reference.pickle')

    @cached_method
    def get_cache(self):
        """
        Retrieve the cache of already generated ReST files.  If it
        doesn't exist, then we just return an empty dictionary.  If it
        is corrupted, return an empty dictionary.
        """
        filename = self.cache_filename()
        if not os.path.exists(filename):
            return {}
        file = open(self.cache_filename(), 'rb')
        try:
            cache = cPickle.load(file)
        except Exception:
            logger.debug("Cache file '%s' is corrupted; ignoring it..."% filename)
            cache = {}
        else:
            logger.debug("Loaded .rst file cache: %s", filename)
        finally:
            file.close()
        return cache

    def save_cache(self):
        """
        Save the cache of already generated ReST files.
        """
        cache = self.get_cache()

        cache['option_inherited'] = opts.INHERITED
        cache['option_underscore'] = opts.UNDERSCORE

        file = open(self.cache_filename(), 'wb')
        cPickle.dump(cache, file)
        file.close()
        logger.debug("Saved .rst file cache: %s", self.cache_filename())

    def get_sphinx_environment(self):
        """
        Returns the Sphinx environment for this project.
        """
        class Foo(object):
            pass
        config = Foo()
        config.values = []

        env_pickle = os.path.join(self._doctrees_dir(), 'environment.pickle')
        try:
            env = BuildEnvironment.frompickle(self.dir, config, env_pickle)
            logger.debug("Opened Sphinx environment: %s", env_pickle)
            return env
        except IOError as err:
            logger.debug("Failed to open Sphinx environment: %s", err)
            pass

    def update_mtimes(self):
        """
        Updates the modification times for ReST files in the Sphinx
        environment for this project.
        """
        env = self.get_sphinx_environment()
        if env is not None:
            for doc in env.all_docs:
                env.all_docs[doc] = time.time()
            logger.info("Updated %d .rst file mtimes", len(env.all_docs))
            # This is the only place we need to save (as opposed to
            # load) Sphinx's pickle, so we do it right here.
            env_pickle = os.path.join(self._doctrees_dir(),
                                      'environment.pickle')

            # When cloning a new branch (see
            # SAGE_LOCAL/bin/sage-clone), we hard link the doc output.
            # To avoid making unlinked, potentially inconsistent
            # copies of the environment, we *don't* use
            # env.topickle(env_pickle), which first writes a temporary
            # file.  We adapt sphinx.environment's
            # BuildEnvironment.topickle:

            # remove unpicklable attributes
            env.set_warnfunc(None)
            del env.config.values
            picklefile = open(env_pickle, 'wb')
            # remove potentially pickling-problematic values from config
            for key, val in vars(env.config).items():
                if key.startswith('_') or isinstance(val, (types.ModuleType,
                                                           types.FunctionType,
                                                           type)):
                    del env.config[key]
            try:
                cPickle.dump(env, picklefile, cPickle.HIGHEST_PROTOCOL)
            finally:
                picklefile.close()

            logger.debug("Saved Sphinx environment: %s", env_pickle)

    def get_modified_modules(self):
        """
        Returns an iterator for all the modules that have been modified
        since the documentation was last built.
        """
        env = self.get_sphinx_environment()
        if env is None:
            logger.debug("Stopped check for modified modules.")
            return
        try:
            added, changed, removed = env.get_outdated_files(False)
            logger.info("Sphinx found %d modified modules", len(changed))
        except OSError as err:
            logger.debug("Sphinx failed to determine modified modules: %s", err)
            self.clean_auto()
            return
        for name in changed:
            # Only pay attention to files in a directory sage/... In
            # particular, don't treat a file like 'sagetex.rst' in
            # doc/en/reference/misc as an autogenerated file: see
            # #14199.
            if name.startswith('sage' + os.sep):
                yield name

    def print_modified_modules(self):
        """
        Prints a list of all the modules that have been modified since
        the documentation was last built.
        """
        for module_name in self.get_modified_modules():
            print(module_name)

    def get_all_rst_files(self, exclude_sage=True):
        """
        Returns an iterator for all rst files which are not
        autogenerated.
        """
        for directory, subdirs, files in os.walk(self.dir):
            if exclude_sage and directory.startswith(os.path.join(self.dir, 'sage')):
                continue
            for filename in files:
                if not filename.endswith('.rst'):
                    continue
                yield os.path.join(directory, filename)

    def get_all_included_modules(self):
        """
        Returns an iterator for all modules which are included in the
        reference manual.
        """
        for filename in self.get_all_rst_files():
            for module in self.get_modules(filename):
                yield module

    def get_newly_included_modules(self, save=False):
        """
        Returns an iterator for all modules that appear in the
        toctrees that don't appear in the cache.
        """
        cache = self.get_cache()
        new_modules = 0
        for module in self.get_all_included_modules():
            if module not in cache:
                cache[module] = True
                new_modules += 1
                yield module
        logger.info("Found %d newly included modules", new_modules)
        if save:
            self.save_cache()

    def print_newly_included_modules(self):
        """
        Prints all of the modules that appear in the toctrees that
        don't appear in the cache.
        """
        for module_name in self.get_newly_included_modules():
            print(module_name)

    def get_modules(self, filename):
        """
        Given a filename for a ReST file, return an iterator for
        all of the autogenerated ReST files that it includes.
        """
        #Create the regular expression used to detect an autogenerated file
        auto_re = re.compile('^\s*(..\/)*(sage(nb)?\/[\w\/]*)\s*$')

        #Read the lines
        f = open(filename)
        lines = f.readlines()
        f.close()

        for line in lines:
            match = auto_re.match(line)
            if match:
                yield match.group(2).replace(os.path.sep, '.')

    def get_module_docstring_title(self, module_name):
        """
        Returns the title of the module from its docstring.
        """
        #Try to import the module
        try:
            __import__(module_name)
        except ImportError as err:
            logger.error("Warning: Could not import %s %s", module_name, err)
            return "UNABLE TO IMPORT MODULE"
        module = sys.modules[module_name]

        #Get the docstring
        doc = module.__doc__
        if doc is None:
            doc = module.doc if hasattr(module, 'doc') else ""

        #Extract the title
        i = doc.find('\n')
        if i != -1:
            return doc[i+1:].lstrip().splitlines()[0]
        else:
            return doc

    def auto_rest_filename(self, module_name):
        """
        Returns the name of the file associated to a given module

        EXAMPLES::

            sage: from sage_setup.docbuild import ReferenceSubBuilder
            sage: ReferenceSubBuilder("reference").auto_rest_filename("sage.combinat.partition")
            '.../doc/en/reference/sage/combinat/partition.rst'
        """
        return self.dir + os.path.sep + module_name.replace('.',os.path.sep) + '.rst'

    def write_auto_rest_file(self, module_name):
        """
        Writes the autogenerated ReST file for module_name.
        """
        if not module_name.startswith('sage'):
            return
        filename = self.auto_rest_filename(module_name)
        sage_makedirs(os.path.dirname(filename))

        outfile = open(filename, 'w')

        title = self.get_module_docstring_title(module_name)

        if title == '':
            logger.error("Warning: Missing title for %s", module_name)
            title = "MISSING TITLE"

        # Don't doctest the autogenerated file.
        outfile.write(".. nodoctest\n\n")
        # Now write the actual content.
        outfile.write(".. _%s:\n\n"%(module_name.replace(".__init__","")))
        outfile.write(title + '\n')
        outfile.write('='*len(title) + "\n\n")
        outfile.write('.. This file has been autogenerated.\n\n')

        inherited = ':inherited-members:' if opts.INHERITED else ''

        automodule = '''
.. automodule:: %s
   :members:
   :undoc-members:
   :show-inheritance:
   %s

'''
        outfile.write(automodule % (module_name, inherited))

        outfile.close()

    def clean_auto(self):
        """
        Remove the cache file for the autogenerated files as well as
        the files themselves.
        """
        if os.path.exists(self.cache_filename()):
            os.unlink(self.cache_filename())
            logger.debug("Deleted .rst cache file: %s", self.cache_filename())

        try:
            shutil.rmtree(os.path.join(self.dir, 'sage'))
            logger.debug("Deleted auto-generated .rst files in: %s",
                         os.path.join(self.dir, 'sage'))
        except OSError:
            pass

    def get_unincluded_modules(self):
        """
        Returns an iterator for all the modules in the Sage library
        which are not included in the reference manual.
        """
        #Make a dictionary of the included modules
        included_modules = {}
        for module_name in self.get_all_included_modules():
            included_modules[module_name] = True

        base_path = os.path.join(SAGE_SRC, 'sage')
        for directory, subdirs, files in os.walk(base_path):
            for filename in files:
                if not (filename.endswith('.py') or
                        filename.endswith('.pyx')):
                    continue

                path = os.path.join(directory, filename)

                #Create the module name
                module_name = path[len(base_path):].replace(os.path.sep, '.')
                module_name = 'sage' + module_name
                module_name = module_name[:-4] if module_name.endswith('pyx') else module_name[:-3]

                #Exclude some ones  -- we don't want init the manual
                if module_name.endswith('__init__') or module_name.endswith('all'):
                    continue

                if module_name not in included_modules:
                    yield module_name

    def print_unincluded_modules(self):
        """
        Prints all of the modules which are not included in the Sage
        reference manual.
        """
        for module_name in self.get_unincluded_modules():
            print(module_name)

    def print_included_modules(self):
        """
        Prints all of the modules that are included in the Sage reference
        manual.
        """
        for module_name in self.get_all_included_modules():
            print(module_name)
