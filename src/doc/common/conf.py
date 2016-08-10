import sys, os, sphinx
from sage.env import SAGE_DOC_SRC, SAGE_DOC, SAGE_SRC
from datetime import date


# If your extensions are in another directory, add it here.
sys.path.append(os.path.join(SAGE_SRC, "sage_setup", "docbuild", "ext"))

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['inventory_builder', 'multidocs',
              'sage_autodoc',  'sphinx.ext.graphviz',
              'sphinx.ext.inheritance_diagram', 'sphinx.ext.todo',
              'sphinx.ext.extlinks', 'matplotlib.sphinxext.plot_directive']

# This code is executed before each ".. PLOT::" directive in the Sphinx
# documentation. It defines a 'sphinx_plot' function that displays a Sage object
# through mathplotlib, so that it will be displayed in the HTML doc
plot_html_show_source_link = False
plot_pre_code = """
def sphinx_plot(plot):
    import matplotlib.image as mpimg
    from sage.misc.temporary_file import tmp_filename
    import matplotlib.pyplot as plt
    if os.environ.get('SAGE_SKIP_PLOT_DIRECTIVE', 'no') != 'yes':
        fn = tmp_filename(ext=".png")
        plot.plot().save(fn)
        img = mpimg.imread(fn)
        plt.imshow(img)
        plt.margins(0)
        plt.axis("off")
        plt.tight_layout(pad=0)

from sage.all_cmdline import *
"""

plot_html_show_formats = False

# We do *not* fully initialize intersphinx since we call it by hand
# in find_sage_dangling_links.
#, 'sphinx.ext.intersphinx']


# Add any paths that contain templates here, relative to this directory.
templates_path = [os.path.join(SAGE_DOC_SRC, 'common', 'templates'), 'templates']

# The suffix of source filenames.
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u""
copyright = u"2005--{}, The Sage Development Team".format(date.today().year)

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
#
# The short X.Y version.
from sage.version import version
release = version

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of glob-style patterns that should be excluded when looking for
# source files. [1] They are matched against the source file names
# relative to the source directory, using slashes as directory
# separators on all platforms.
exclude_patterns = ['.build']

# The reST default role (used for this markup: `text`) to use for all documents.
default_role = 'math'

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.  NOTE:
# This overrides a HTML theme's corresponding setting (see below).
pygments_style = 'sphinx'

# GraphViz includes dot, neato, twopi, circo, fdp.
graphviz_dot = 'dot'
inheritance_graph_attrs = { 'rankdir' : 'BT' }
inheritance_node_attrs = { 'height' : 0.5, 'fontsize' : 12, 'shape' : 'oval' }
inheritance_edge_attrs = {}

# Extension configuration
# -----------------------

# include the todos
todo_include_todos = True


# Cross-links to other project's online documentation.
intersphinx_mapping = {
    'python': ('https://docs.python.org/',
                os.path.join(SAGE_DOC_SRC, "common", "python.inv"))}

def set_intersphinx_mappings(app):
    """
    Add precompiled inventory (the objects.inv)
    """
    refpath = os.path.join(SAGE_DOC, "html", "en", "reference")
    invpath = os.path.join(SAGE_DOC, "inventory", "en", "reference")
    if app.config.multidoc_first_pass == 1 or \
            not (os.path.exists(refpath) and os.path.exists(invpath)):
        app.config.intersphinx_mapping = {}
        return

    app.config.intersphinx_mapping = intersphinx_mapping

    # Add master intersphinx mapping
    dst = os.path.join(invpath, 'objects.inv')
    app.config.intersphinx_mapping['sagemath'] = (refpath, dst)

    # Add intersphinx mapping for subdirectories
    # We intentionally do not name these such that these get higher
    # priority in case of conflicts
    for directory in os.listdir(os.path.join(invpath)):
        if os.path.isdir(os.path.join(invpath, directory)):
            src = os.path.join(refpath, directory)
            dst = os.path.join(invpath, directory, 'objects.inv')
            app.config.intersphinx_mapping[src] = dst


pythonversion = sys.version.split(' ')[0]
# Python and Sage trac ticket shortcuts. For example, :trac:`7549` .

# Sage trac ticket shortcuts. For example, :trac:`7549` .
extlinks = {
    'python': ('https://docs.python.org/release/'+pythonversion+'/%s', ''),
    'trac': ('https://trac.sagemath.org/%s', 'trac ticket #'),
    'wikipedia': ('https://en.wikipedia.org/wiki/%s', 'Wikipedia article '),
    'arxiv': ('http://arxiv.org/abs/%s', 'Arxiv '),
    'oeis': ('https://oeis.org/%s', 'OEIS sequence '),
    'doi': ('https://dx.doi.org/%s', 'doi:'),
    'pari': ('http://pari.math.u-bordeaux.fr/dochtml/html/4.html#%s', 'pari:')
    'mathscinet': ('http://www.ams.org/mathscinet-getitem?mr=%s', 'MathSciNet ')
    }

# By default document are not master.
multidocs_is_master = True

# Options for HTML output
# -----------------------

# HTML theme (e.g., 'default', 'sphinxdoc').  We use a custom Sage
# theme to set a Pygments style, stylesheet, and insert MathJax macros. See
# the directory doc/common/themes/sage/ for files comprising the custom Sage
# theme.
html_theme = 'sage'

# Theme options are theme-specific and customize the look and feel of
# a theme further.  For a list of options available for each theme,
# see the documentation.
html_theme_options = {}

# Add any paths that contain custom themes here, relative to this directory.
html_theme_path = [os.path.join(SAGE_DOC_SRC, 'common', 'themes')]

# HTML style sheet NOTE: This overrides a HTML theme's corresponding
# setting.
#html_style = 'default.css'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
#html_title = None

# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The name of an image file (within the static path) to place at the top of
# the sidebar.
#html_logo = 'sagelogo-word.ico'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = 'favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [os.path.join(SAGE_DOC_SRC, 'common', 'static'), 'static']

# We use MathJax to build the documentation unless the environment
# variable SAGE_DOC_MATHJAX is set to "no" or "False".  (Note that if
# the user does not set this variable, then the script sage-env sets
# it to "True".)

if (os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'no'
            and os.environ.get('SAGE_DOC_MATHJAX', 'no') != 'False'):

    extensions.append('sphinx.ext.mathjax')
    mathjax_path = 'MathJax.js?config=TeX-AMS_HTML-full,../mathjax_sage.js'

    from sage.misc.latex_macros import sage_mathjax_macros
    html_theme_options['mathjax_macros'] = sage_mathjax_macros()

    from pkg_resources import Requirement, working_set
    sagenb_path = working_set.find(Requirement.parse('sagenb')).location
    mathjax_relative = os.path.join('sagenb','data','mathjax')

    # It would be really nice if sphinx would copy the entire mathjax directory,
    # (so we could have a _static/mathjax directory), rather than the contents of the directory

    mathjax_static = os.path.join(sagenb_path, mathjax_relative)
    html_static_path.append(mathjax_static)
    exclude_patterns += ['**/'+os.path.join(mathjax_relative, i)
                         for i in ('docs', 'README*', 'test',
                                   'unpacked', 'LICENSE')]
else:
     extensions.append('sphinx.ext.pngmath')


# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_use_modindex = True

# A list of prefixes that are ignored for sorting the Python module index ( if
# this is set to ['foo.'], then foo.bar is shown under B, not F). Works only
# for the HTML builder currently.
modindex_common_prefix = ['sage.']

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
html_split_index = True

# If true, the reST sources are included in the HTML build as _sources/<name>.
#html_copy_source = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = ''

# Output file base name for HTML help builder.
#htmlhelp_basename = ''


# Options for LaTeX output
# ------------------------
# See http://sphinx-doc.org/config.html#confval-latex_elements
latex_elements = {}

# The paper size ('letterpaper' or 'a4paper').
#latex_elements['papersize'] = 'letterpaper'

# The font size ('10pt', '11pt' or '12pt').
#latex_elements['pointsize'] = '10pt'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = []

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = 'sagelogo-word.png'

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
latex_elements['preamble'] = r"""
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{textcomp}
\usepackage{mathrsfs}

% Only declare unicode characters when compiling with pdftex; E.g. japanese 
% tutorial does not use pdftex
\ifPDFTeX
    \DeclareUnicodeCharacter{01CE}{\capitalcaron a}
    \DeclareUnicodeCharacter{0428}{cyrillic Sha}
    \DeclareUnicodeCharacter{250C}{+}
    \DeclareUnicodeCharacter{2510}{+}
    \DeclareUnicodeCharacter{2514}{+}
    \DeclareUnicodeCharacter{2518}{+}
    \DeclareUnicodeCharacter{253C}{+}

    \DeclareUnicodeCharacter{03B1}{\ensuremath{\alpha}}
    \DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
    \DeclareUnicodeCharacter{03B3}{\ensuremath{\gamma}}
    \DeclareUnicodeCharacter{0393}{\ensuremath{\Gamma}}
    \DeclareUnicodeCharacter{03B4}{\ensuremath{\delta}}
    \DeclareUnicodeCharacter{0394}{\ensuremath{\Delta}}
    \DeclareUnicodeCharacter{03B5}{\ensuremath{\varepsilon}}
    \DeclareUnicodeCharacter{03B6}{\ensuremath{\zeta}}
    \DeclareUnicodeCharacter{03B7}{\ensuremath{\eta}}
    \DeclareUnicodeCharacter{03B8}{\ensuremath{\vartheta}}
    \DeclareUnicodeCharacter{0398}{\ensuremath{\Theta}}
    \DeclareUnicodeCharacter{03BA}{\ensuremath{\kappa}}
    \DeclareUnicodeCharacter{03BB}{\ensuremath{\lambda}}
    \DeclareUnicodeCharacter{039B}{\ensuremath{\Lambda}}
    \DeclareUnicodeCharacter{00B5}{\ensuremath{\mu}}      % micron sign
    \DeclareUnicodeCharacter{03BC}{\ensuremath{\mu}}
    \DeclareUnicodeCharacter{03BD}{\ensuremath{\nu}}
    \DeclareUnicodeCharacter{03BE}{\ensuremath{\xi}}
    \DeclareUnicodeCharacter{039E}{\ensuremath{\Xi}}
    \DeclareUnicodeCharacter{03B9}{\ensuremath{\iota}}
    \DeclareUnicodeCharacter{03C0}{\ensuremath{\pi}}
    \DeclareUnicodeCharacter{03A0}{\ensuremath{\Pi}}
    \DeclareUnicodeCharacter{03C1}{\ensuremath{\rho}}
    \DeclareUnicodeCharacter{03C3}{\ensuremath{\sigma}}
    \DeclareUnicodeCharacter{03A3}{\ensuremath{\Sigma}}
    \DeclareUnicodeCharacter{03C4}{\ensuremath{\tau}}
    \DeclareUnicodeCharacter{03C6}{\ensuremath{\varphi}}
    \DeclareUnicodeCharacter{03A6}{\ensuremath{\Phi}}
    \DeclareUnicodeCharacter{03C7}{\ensuremath{\chi}}
    \DeclareUnicodeCharacter{03C8}{\ensuremath{\psi}}
    \DeclareUnicodeCharacter{03A8}{\ensuremath{\Psi}}
    \DeclareUnicodeCharacter{03C9}{\ensuremath{\omega}}
    \DeclareUnicodeCharacter{03A9}{\ensuremath{\Omega}}
    \DeclareUnicodeCharacter{03C5}{\ensuremath{\upsilon}}
    \DeclareUnicodeCharacter{03A5}{\ensuremath{\Upsilon}}
    \DeclareUnicodeCharacter{2113}{\ell}
    
    \DeclareUnicodeCharacter{221A}{\ensuremath{\sqrt{}}}
    \DeclareUnicodeCharacter{2264}{\leq}
    \DeclareUnicodeCharacter{2265}{\geq}
    \DeclareUnicodeCharacter{221E}{\infty}
    \DeclareUnicodeCharacter{2211}{\sum}
    \DeclareUnicodeCharacter{2208}{\in}
    \DeclareUnicodeCharacter{2209}{\notin}
    \DeclareUnicodeCharacter{2202}{\partial}
    \DeclareUnicodeCharacter{222B}{\ensuremath{\int}}
    \DeclareUnicodeCharacter{2148}{\id}
    \DeclareUnicodeCharacter{2248}{\approx}
    \DeclareUnicodeCharacter{2260}{\neq}
    \DeclareUnicodeCharacter{00B1}{\pm}
    \DeclareUnicodeCharacter{2A02}{\otimes}
    \DeclareUnicodeCharacter{2A01}{\oplus}
    \DeclareUnicodeCharacter{00BD}{\nicefrac{1}{2}}
    \DeclareUnicodeCharacter{00D7}{\times}
    \DeclareUnicodeCharacter{00B7}{\cdot}
    \DeclareUnicodeCharacter{230A}{\lfloor}
    \DeclareUnicodeCharacter{230B}{\rfloor}
    \DeclareUnicodeCharacter{2308}{\lceil}
    \DeclareUnicodeCharacter{2309}{\rceil}
    \DeclareUnicodeCharacter{22C5}{\ensuremath{\cdot}}
    
    \newcommand{\sageMexSymbol}[1]
    {{\fontencoding{OMX}\fontfamily{cmex}\selectfont\raisebox{0.75em}{\symbol{#1}}}}
    \DeclareUnicodeCharacter{239B}{\sageMexSymbol{"30}} % parenlefttp
    \DeclareUnicodeCharacter{239C}{\sageMexSymbol{"42}} % parenleftex
    \DeclareUnicodeCharacter{239D}{\sageMexSymbol{"40}} % parenleftbt
    \DeclareUnicodeCharacter{239E}{\sageMexSymbol{"31}} % parenrighttp
    \DeclareUnicodeCharacter{239F}{\sageMexSymbol{"43}} % parenrightex
    \DeclareUnicodeCharacter{23A0}{\sageMexSymbol{"41}} % parenrightbt
    \DeclareUnicodeCharacter{23A1}{\sageMexSymbol{"32}} % bracketlefttp
    \DeclareUnicodeCharacter{23A2}{\sageMexSymbol{"36}} % bracketleftex
    \DeclareUnicodeCharacter{23A3}{\sageMexSymbol{"34}} % bracketleftbt
    \DeclareUnicodeCharacter{23A4}{\sageMexSymbol{"33}} % bracketrighttp
    \DeclareUnicodeCharacter{23A5}{\sageMexSymbol{"37}} % bracketrightex
    \DeclareUnicodeCharacter{23A6}{\sageMexSymbol{"35}} % bracketrightbt
    
    \DeclareUnicodeCharacter{23A7}{\sageMexSymbol{"38}} % curly brace left top
    \DeclareUnicodeCharacter{23A8}{\sageMexSymbol{"3C}} % curly brace left middle
    \DeclareUnicodeCharacter{23A9}{\sageMexSymbol{"3A}} % curly brace left bottom
    \DeclareUnicodeCharacter{23AA}{\sageMexSymbol{"3E}} % curly brace extension
    \DeclareUnicodeCharacter{23AB}{\sageMexSymbol{"39}} % curly brace right top
    \DeclareUnicodeCharacter{23AC}{\sageMexSymbol{"3D}} % curly brace right middle
    \DeclareUnicodeCharacter{23AD}{\sageMexSymbol{"3B}} % curly brace right bottom
    \DeclareUnicodeCharacter{23B0}{\{} % 2-line curly brace left top half  (not in cmex)
    \DeclareUnicodeCharacter{23B1}{\}} % 2-line curly brace right top half (not in cmex)
    
    \DeclareUnicodeCharacter{2320}{\ensuremath{\int}} % top half integral
    \DeclareUnicodeCharacter{2321}{\ensuremath{\int}} % bottom half integral
    \DeclareUnicodeCharacter{23AE}{\ensuremath{\|}} % integral extenison
    
    \DeclareUnicodeCharacter{2571}{/}   % Box drawings light diagonal upper right to lower left
\fi

\let\textLaTeX\LaTeX
\renewcommand*{\LaTeX}{\hbox{\textLaTeX}}
"""

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_use_modindex = True

#####################################################
# add LaTeX macros for Sage

from sage.misc.latex_macros import sage_latex_macros

try:
    pngmath_latex_preamble  # check whether this is already defined
except NameError:
    pngmath_latex_preamble = ""

for macro in sage_latex_macros():
    # used when building latex and pdf versions
    latex_elements['preamble'] += macro + '\n'
    # used when building html version
    pngmath_latex_preamble += macro + '\n'

#####################################################

def process_docstring_aliases(app, what, name, obj, options, docstringlines):
    """
    Change the docstrings for aliases to point to the original object.
    """
    basename = name.rpartition('.')[2]
    if hasattr(obj, '__name__') and obj.__name__ != basename:
        docstringlines[:] = ['See :obj:`%s`.' % name]

def process_directives(app, what, name, obj, options, docstringlines):
    """
    Remove 'nodetex' and other directives from the first line of any
    docstring where they appear.
    """
    if len(docstringlines) == 0:
        return
    first_line = docstringlines[0]
    directives = [ d.lower() for d in first_line.split(',') ]
    if 'nodetex' in directives:
        docstringlines.pop(0)

def process_docstring_cython(app, what, name, obj, options, docstringlines):
    """
    Remove Cython's filename and location embedding.
    """
    if len(docstringlines) <= 1:
        return

    first_line = docstringlines[0]
    if first_line.startswith('File:') and '(starting at' in first_line:
        #Remove the first two lines
        docstringlines.pop(0)
        docstringlines.pop(0)

def process_docstring_module_title(app, what, name, obj, options, docstringlines):
    """
    Removes the first line from the beginning of the module's docstring.  This
    corresponds to the title of the module's documentation page.
    """
    if what != "module":
        return

    #Remove any additional blank lines at the beginning
    title_removed = False
    while len(docstringlines) > 1 and not title_removed:
        if docstringlines[0].strip() != "":
            title_removed = True
        docstringlines.pop(0)

    #Remove any additional blank lines at the beginning
    while len(docstringlines) > 1:
        if docstringlines[0].strip() == "":
            docstringlines.pop(0)
        else:
            break

skip_picklability_check_modules = [
    #'sage.misc.nested_class_test', # for test only
    'sage.misc.latex',
    'sage.misc.explain_pickle',
    '__builtin__',
]

def check_nested_class_picklability(app, what, name, obj, skip, options):
    """
    Print a warning if pickling is broken for nested classes.
    """
    if hasattr(obj, '__dict__') and hasattr(obj, '__module__'):
        # Check picklability of nested classes.  Adapted from
        # sage.misc.nested_class.modify_for_nested_pickle.
        module = sys.modules[obj.__module__]
        for (nm, v) in obj.__dict__.iteritems():
            if (isinstance(v, type) and
                v.__name__ == nm and
                v.__module__ == module.__name__ and
                getattr(module, nm, None) is not v and
                v.__module__ not in skip_picklability_check_modules):
                # OK, probably this is an *unpicklable* nested class.
                app.warn('Pickling of nested class %r is probably broken. '
                         'Please set __metaclass__ of the parent class to '
                         'sage.misc.nested_class.NestedClassMetaclass.' % (
                        v.__module__ + '.' + name + '.' + nm))

def skip_member(app, what, name, obj, skip, options):
    """
    To suppress Sphinx warnings / errors, we

    - Don't include [aliases of] builtins.

    - Don't include the docstring for any nested class which has been
      inserted into its module by
      :class:`sage.misc.NestedClassMetaclass` only for pickling.  The
      class will be properly documented inside its surrounding class.

    - Don't include
      sagenb.notebook.twist.userchild_download_worksheets.zip.

    - Optionally, check whether pickling is broken for nested classes.

    - Optionally, include objects whose name begins with an underscore
      ('_'), i.e., "private" or "hidden" attributes, methods, etc.

    Otherwise, we abide by Sphinx's decision.  Note: The object
    ``obj`` is excluded (included) if this handler returns True
    (False).
    """
    if 'SAGE_CHECK_NESTED' in os.environ:
        check_nested_class_picklability(app, what, name, obj, skip, options)

    if getattr(obj, '__module__', None) == '__builtin__':
        return True

    objname = getattr(obj, "__name__", None)
    if objname is not None:
        if objname.find('.') != -1 and objname.split('.')[-1] != name:
            return True

    if name.find("userchild_download_worksheets.zip") != -1:
        return True

    if 'SAGE_DOC_UNDERSCORE' in os.environ:
        if name.split('.')[-1].startswith('_'):
            return False

    return skip

def process_dollars(app, what, name, obj, options, docstringlines):
    r"""
    Replace dollar signs with backticks.
    See sage.misc.sagedoc.process_dollars for more information
    """
    if len(docstringlines) > 0 and name.find("process_dollars") == -1:
        from sage.misc.sagedoc import process_dollars as sagedoc_dollars
        s = sagedoc_dollars("\n".join(docstringlines))
        lines = s.split("\n")
        for i in range(len(lines)):
            docstringlines[i] = lines[i]

def process_inherited(app, what, name, obj, options, docstringlines):
    """
    If we're including inherited members, omit their docstrings.
    """
    if not options.get('inherited-members'):
        return

    if what in ['class', 'data', 'exception', 'function', 'module']:
        return

    name = name.split('.')[-1]

    if what == 'method' and hasattr(obj, 'im_class'):
        if name in obj.im_class.__dict__.keys():
            return

    if what == 'attribute' and hasattr(obj, '__objclass__'):
        if name in obj.__objclass__.__dict__.keys():
            return

    for i in xrange(len(docstringlines)):
        docstringlines.pop()

dangling_debug = False

def debug_inf(app, message):
    if dangling_debug: app.info(message)

def call_intersphinx(app, env, node, contnode):
    """
    Call intersphinx and make links between Sage manuals relative.

    TESTS:

    Check that the link from the thematic tutorials to the reference
    manual is relative, see :trac:`20118`::

        sage: from sage.env import SAGE_DOC
        sage: thematic_index = os.path.join(SAGE_DOC, "html", "en", "thematic_tutorials", "index.html")
        sage: for line in open(thematic_index).readlines():
        ....:     if "padics" in line:
        ....:         sys.stdout.write(line)
        <li><a class="reference external" href="../reference/padics/sage/rings/padics/tutorial.html#sage-rings-padics-tutorial" title="(in Sage Reference Manual: p-Adics ...)"><span>Introduction to the -adics</span></a></li>
    """
    debug_inf(app, "???? Trying intersphinx for %s"%node['reftarget'])
    builder = app.builder
    res =  sphinx.ext.intersphinx.missing_reference(
        app, env, node, contnode)
    if res:
        # Replace absolute links to $SAGE_DOC by relative links: this
        # allows to copy the whole documentation tree somewhere else
        # without breaking links, see Trac #20118.
        if res['refuri'].startswith(SAGE_DOC):
            here = os.path.dirname(os.path.join(builder.outdir,
                                                node['refdoc']))
            res['refuri'] = os.path.relpath(res['refuri'], here)
            debug_inf(app, "++++ Found at %s"%res['refuri'])
    else:
        debug_inf(app, "---- Intersphinx: %s not Found"%node['reftarget'])
    return res

def find_sage_dangling_links(app, env, node, contnode):
    """
    Try to find dangling link in local module imports or all.py.
    """
    debug_inf(app, "==================== find_sage_dangling_links ")

    reftype = node['reftype']
    reftarget  = node['reftarget']
    try:
        doc = node['refdoc']
    except KeyError:
        debug_inf(app, "-- no refdoc in node %s"%node)
        return None

    debug_inf(app, "Searching %s from %s"%(reftarget, doc))

    # Workaround: in Python's doc 'object', 'list', ... are documented as a
    # function rather than a class
    if reftarget in base_class_as_func and reftype == 'class':
        node['reftype'] = 'func'

    res = call_intersphinx(app, env, node, contnode)
    if res:
        debug_inf(app, "++ DONE %s"%(res['refuri']))
        return res

    if node.get('refdomain') != 'py': # not a python file
       return None

    try:
        module = node['py:module']
        cls    = node['py:class']
    except KeyError:
        debug_inf(app, "-- no module or class for :%s:%s"%(reftype, reftarget))
        return None

    basename = reftarget.split(".")[0]
    try:
        target_module = getattr(sys.modules['sage.all'], basename).__module__
    except AttributeError:
        debug_inf(app, "-- %s not found in sage.all"%(basename))
        return None
    if target_module is None:
        target_module = ""
        debug_inf(app, "?? found in None !!!")

    debug_inf(app, "++ found %s using sage.all in %s"%(basename, target_module))

    newtarget = target_module+'.'+reftarget
    node['reftarget'] = newtarget

    # adapted  from sphinx/domains/python.py
    builder = app.builder
    searchmode = node.hasattr('refspecific') and 1 or 0
    matches =  builder.env.domains['py'].find_obj(
        builder.env, module, cls, newtarget, reftype, searchmode)
    if not matches:
        debug_inf(app, "?? no matching doc for %s"%newtarget)
        return call_intersphinx(app, env, node, contnode)
    elif len(matches) > 1:
        env.warn(target_module,
                 'more than one target found for cross-reference '
                 '%r: %s' % (newtarget,
                             ', '.join(match[0] for match in matches)),
                 node.line)
    name, obj = matches[0]
    debug_inf(app, "++ match = %s %s"%(name, obj))

    from docutils import nodes
    newnode = nodes.reference('', '', internal=True)
    if name == target_module:
        newnode['refid'] = name
    else:
        newnode['refuri'] = builder.get_relative_uri(node['refdoc'], obj[0])
        newnode['refuri'] += '#' + name
        debug_inf(app, "++ DONE at URI %s"%(newnode['refuri']))
    newnode['reftitle'] = name
    newnode.append(contnode)
    return newnode

# lists of basic Python class which are documented as functions
base_class_as_func = [
    'bool', 'complex', 'dict', 'file', 'float',
    'frozenset', 'int', 'list', 'long', 'object',
    'set', 'slice', 'str', 'tuple', 'type', 'unicode', 'xrange']

# Nit picky option configuration: Put here broken links we want to ignore. For
# link to the Python documentation several links where broken because there
# where class listed as functions. Expand the list 'base_class_as_func'
# above instead of marking the link as broken.
nitpick_ignore = [
    ('py:class', 'twisted.web2.resource.Resource'),
    ('py:class', 'twisted.web2.resource.PostableResource')]

def nitpick_patch_config(app):
    """
    Patch the default config for nitpicky

    Calling path_config ensure that nitpicky is not considered as a Sphinx
    environment variable but rather as a Sage environment variable. As a
    consequence, changing it doesn't force the recompilation of the entire
    documentation.
    """
    app.config.values['nitpicky'] = (False, 'sage')
    app.config.values['nitpick_ignore'] = ([], 'sage')

def skip_TESTS_block(app, what, name, obj, options, docstringlines):
    """
    Skip blocks labeled "TESTS:".

    See sage.misc.sagedoc.skip_TESTS_block for more information.
    """
    from sage.misc.sagedoc import skip_TESTS_block as sagedoc_skip_TESTS
    if not docstringlines:
        # No docstring, so don't do anything. See Trac #19932.
        return
    s = sagedoc_skip_TESTS("\n".join(docstringlines))
    lines = s.split("\n")
    for i in range(len(lines)):
        docstringlines[i] = lines[i]
    while len(docstringlines) > len(lines):
        del docstringlines[len(lines)]

from sage.misc.sageinspect import sage_getargspec
autodoc_builtin_argspec = sage_getargspec

def setup(app):
    app.connect('autodoc-process-docstring', process_docstring_cython)
    app.connect('autodoc-process-docstring', process_directives)
    app.connect('autodoc-process-docstring', process_docstring_module_title)
    app.connect('autodoc-process-docstring', process_dollars)
    app.connect('autodoc-process-docstring', process_inherited)
    if os.environ.get('SAGE_SKIP_TESTS_BLOCKS', False):
        app.connect('autodoc-process-docstring', skip_TESTS_block)
    app.connect('autodoc-skip-member', skip_member)

    # When building the standard docs, app.srcdir is set to SAGE_DOC_SRC +
    # 'LANGUAGE/DOCNAME', but when doing introspection, app.srcdir is
    # set to a temporary directory.  We don't want to use intersphinx,
    # etc., when doing introspection.
    if app.srcdir.startswith(SAGE_DOC_SRC):
        app.add_config_value('intersphinx_mapping', {}, False)
        app.add_config_value('intersphinx_cache_limit', 5, False)
        # We do *not* fully initialize intersphinx since we call it by hand
        # in find_sage_dangling_links.
        #   app.connect('missing-reference', missing_reference)
        app.connect('missing-reference', find_sage_dangling_links)
        import sphinx.ext.intersphinx
        app.connect('builder-inited', set_intersphinx_mappings)
        app.connect('builder-inited', sphinx.ext.intersphinx.load_mappings)
        app.connect('builder-inited', nitpick_patch_config)
