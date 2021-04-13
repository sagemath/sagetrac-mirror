# -*- coding: utf-8 -*-
r"""
Interfaces to R

This is the reference to the Sagemath R interface, usable from any
Sage program.

The %r interface creating an R cell in the sage
notebook is decribed in the Notebook manual.

The %R and %%R interface creating an R line or an R cell in the
Jupyter notebook are briefly decribed at the end of this page. This
documentation will be expanded and placed in the Jupyter notebook
manual when this manual exists.

The following examples try to follow "An Introduction to R" which can
be found at http://cran.r-project.org/doc/manuals/R-intro.html .

EXAMPLES:

Simple manipulations; numbers and vectors

The simplest data structure in R is the numeric vector which
consists of an ordered collection of numbers.  To create a
vector named $x$ using the R interface in Sage, you pass the
R interpreter object a list or tuple of numbers::

    sage: x = r([10.4,5.6,3.1,6.4,21.7]); x  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7

You can invert elements of a vector x in R by using the
invert operator or by doing 1/x::

    sage: ~x  # optional - rpy2
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295
    sage: 1/x  # optional - rpy2
    [1] 0.09615385 0.17857143 0.32258065 0.15625000 0.04608295

The following assignment creates a vector $y$ with 11 entries which
consists of two copies of $x$ with a 0 in between::

    sage: y = r([x,0,x]); y  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7  0.0 10.4  5.6  3.1  6.4 21.7

Vector Arithmetic

The following command generates a new vector $v$ of length 11 constructed
by adding together (element by element) $2x$ repeated 2.2 times, $y$
repeated just once, and 1 repeated 11 times::

    sage: v = 2*x+y+1; v  # optional - rpy2
    [1] 32.2 17.8 10.3 20.2 66.1 21.8 22.6 12.8 16.9 50.8 43.5

One can compute the sum of the elements of an R vector in the following
two ways::

    sage: sum(x)  # optional - rpy2
    [1] 47.2
    sage: x.sum()  # optional - rpy2
    [1] 47.2

One can calculate the sample variance of a list of numbers::

    sage: ((x-x.mean())^2/(x.length()-1)).sum()  # optional - rpy2
    [1] 53.853
    sage: x.var()  # optional - rpy2
    [1] 53.853

    sage: x.sort()  # optional - rpy2
    [1] 3.1  5.6  6.4 10.4 21.7
    sage: x.min()  # optional - rpy2
    [1] 3.1
    sage: x.max()  # optional - rpy2
    [1] 21.7
    sage: x  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7

    sage: r(-17).sqrt()  # optional - rpy2
    [1] NaN
    sage: r('-17+0i').sqrt()  # optional - rpy2
    [1] 0+4.123106i

Generating an arithmetic sequence::

    sage: r('1:10')  # optional - rpy2
    [1] 1  2  3  4  5  6  7  8  9 10

Because ``from`` is a keyword in Python, it can't be used
as a keyword argument.  Instead, ``from_`` can be passed, and
R will recognize it as the correct thing::

    sage: r.seq(length=10, from_=-1, by=.2)  # optional - rpy2
    [1] -1.0 -0.8 -0.6 -0.4 -0.2  0.0  0.2  0.4  0.6  0.8

    sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
    sage: x.rep(2)  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(times=2)  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7 10.4  5.6  3.1  6.4 21.7
    sage: x.rep(each=2)  # optional - rpy2
    [1] 10.4 10.4  5.6  5.6  3.1  3.1  6.4  6.4 21.7 21.7

Missing Values::

    sage: na = r('NA')  # optional - rpy2
    sage: z = r([1,2,3,na])  # optional - rpy2
    sage: z  # optional - rpy2
    [1]  1  2  3 NA
    sage: ind = r.is_na(z)  # optional - rpy2
    sage: ind  # optional - rpy2
    [1] FALSE FALSE FALSE  TRUE
    sage: zero = r(0)  # optional - rpy2
    sage: zero / zero  # optional - rpy2
    [1] NaN
    sage: inf = r('Inf')  # optional - rpy2
    sage: inf-inf  # optional - rpy2
    [1] NaN
    sage: r.is_na(inf)  # optional - rpy2
    [1] FALSE
    sage: r.is_na(inf-inf)  # optional - rpy2
    [1] TRUE
    sage: r.is_na(zero/zero)  # optional - rpy2
    [1] TRUE
    sage: r.is_na(na)  # optional - rpy2
    [1] TRUE
    sage: r.is_nan(inf-inf)  # optional - rpy2
    [1] TRUE
    sage: r.is_nan(zero/zero)  # optional - rpy2
    [1] TRUE
    sage: r.is_nan(na)  # optional - rpy2
    [1] FALSE


Character Vectors::

    sage: labs = r.paste('c("X","Y")', '1:10', sep='""'); labs  # optional - rpy2
    [1] "X1"  "Y2"  "X3"  "Y4"  "X5"  "Y6"  "X7"  "Y8"  "X9"  "Y10"


Index vectors; selecting and modifying subsets of a data set::

    sage: na = r('NA')  # optional - rpy2
    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: x['!is.na(self)']  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7

    sage: x = r([10.4,5.6,3.1,6.4,21.7,na]); x  # optional - rpy2
    [1] 10.4  5.6  3.1  6.4 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']  # optional - rpy2
    [1] 11.4  6.6  4.1  7.4 22.7
    sage: x = r([10.4,-2,3.1,-0.5,21.7,na]); x  # optional - rpy2
    [1] 10.4 -2.0  3.1 -0.5 21.7   NA
    sage: (x+1)['(!is.na(self)) & self>0']  # optional - rpy2
    [1] 11.4  4.1  0.5 22.7

Distributions::

    sage: r.options(width="60")  # optional - rpy2
    $width
    [1] 80

    sage: rr = r.dnorm(r.seq(-3,3,0.1))  # optional - rpy2
    sage: rr  # optional - rpy2
     [1] 0.004431848 0.005952532 0.007915452 0.010420935
     [5] 0.013582969 0.017528300 0.022394530 0.028327038
     [9] 0.035474593 0.043983596 0.053990967 0.065615815
    [13] 0.078950158 0.094049077 0.110920835 0.129517596
    [17] 0.149727466 0.171368592 0.194186055 0.217852177
    [21] 0.241970725 0.266085250 0.289691553 0.312253933
    [25] 0.333224603 0.352065327 0.368270140 0.381387815
    [29] 0.391042694 0.396952547 0.398942280 0.396952547
    [33] 0.391042694 0.381387815 0.368270140 0.352065327
    [37] 0.333224603 0.312253933 0.289691553 0.266085250
    [41] 0.241970725 0.217852177 0.194186055 0.171368592
    [45] 0.149727466 0.129517596 0.110920835 0.094049077
    [49] 0.078950158 0.065615815 0.053990967 0.043983596
    [53] 0.035474593 0.028327038 0.022394530 0.017528300
    [57] 0.013582969 0.010420935 0.007915452 0.005952532
    [61] 0.004431848

Convert R Data Structures to Python/Sage::

    sage: rr = r.dnorm(r.seq(-3,3,0.1))  # optional - rpy2
    sage: sum(rr._sage_())  # optional - rpy2
    9.9772125168981...

Or you get a dictionary to be able to access all the information::

    sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))  # optional - rpy2
    sage: rs  # optional - rpy2
       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
      1.000   1.750   3.000   2.875   4.000   5.000
      sage: d = rs._sage_()  # optional - rpy2
      sage: d['DATA']  # optional - rpy2
      [1, 1.75, 3, 2.875, 4, 5]
      sage: d['_Names']  # optional - rpy2
      ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.']
      sage: d['_r_class']  # optional - rpy2
      ['summaryDefault', 'table']

It is also possible to access the plotting capabilities of R
through Sage.  For more information see the documentation of
r.plot() or r.png().

THE JUPYTER NOTEBOOK INTERFACE (work in progress).

The %r interface described in the Sage notebook manual is not useful
in the Jupyter notebook : it creates a inferior R interpreter which
cannot be escaped.

The RPy2 library allows the creation of an R cell in the Jupyter
notebook analogous to the %r escape in command line or %r cell in a
Sage notebook.

The interface is loaded by a cell containing the sole code:

"%load_ext rpy2.ipython"

After execution of this code, the %R and %%R magics are available:

- %R allows the execution of a single line of R code. Data exchange is
   possible via the -i and -o options. Do "%R?" in a standalone cell
   to get the documentation.

- %%R allows the execution in R of the whole text of a cell, with
    similar options (do "%%R?" in a standalone cell for
    documentation).

A few important points must be noted:

- The R interpreter launched by this interface IS (currently)
  DIFFERENT from the R interpreter used br other r... functions.

- Data exchanged via the -i and -o options have a format DIFFERENT
  from the format used by the r... functions (RPy2 mostly uses arrays,
  and bugs the user to use the pandas Python package).

- R graphics are (beautifully) displayed in output cells, but are not
  directly importable. You have to save them as .png, .pdf or .svg
  files and import them in Sage for further use.

In its current incarnation, this interface is mostly useful to
statisticians needing Sage for a few symbolic computations but mostly
using R for applied work.

AUTHORS:

- Mike Hansen (2007-11-01)
- William Stein (2008-04-19)
- Harald Schilly (2008-03-20)
- Mike Hansen (2008-04-19)
- Emmanuel Charpentier (2015-12-12, RPy2 interface)
"""

##########################################################################
#
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#                     2007 Mike Hansen   <mhansen@gmail.com>
#                     2008 Harald Schilly <harald.schilly@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

from .interface import Interface, InterfaceElement, InterfaceFunction, InterfaceFunctionElement
from sage.env import DOT_SAGE
import re
from sage.structure.element import parent
from sage.interfaces.tab_completion import ExtraTabCompletion
from sage.docs.instancedoc import instancedoc

# see the _lazy_init for some reasoning behind the lazy imports
from sage.misc.lazy_import import lazy_import
lazy_import("rpy2", "robjects")
lazy_import("rpy2.robjects", "packages", "rpy2_packages")
lazy_import("rpy2.robjects.conversion", "localconverter")

# for help page fetching
lazy_import("rpy2.robjects.help", "Package")
lazy_import("rpy2", "rinterface")

COMMANDS_CACHE = '%s/r_commandlist.sobj'%DOT_SAGE

#there is a mirror network, but lets take #1 for now
RRepositoryURL = "http://cran.r-project.org/"
RFilteredPackages = ['.GlobalEnv']

# crosscheck with https://svn.r-project.org/R/trunk/src/main/names.c
# but package:base should cover this. i think.
RBaseCommands = ['c', "NULL", "NA", "True", "False", "Inf", "NaN"]

def _setup_r_to_sage_converter():
    """
    Set up a the converter used to convert from rpy2's
    representation of R objects to the one sage expects.

    EXAMPLES::

    Test

    Simple numeric values are represented as vectors in R. So `1` would actually
    be an array of length 1. We convert all vectors of length 1 to simple values,
    whether or not they "originally" were simple values or not:

        sage: r([42]).sage()  # optional - rpy2
        42

        sage: r(42).sage()  # optional - rpy2
        42

        sage: r('c("foo")').sage()  # optional - rpy2
        'foo'

    Arrays of length greater than one are treated normally:

        sage: r([42, 43]).sage()  # optional - rpy2
        [42, 43]

    We also convert all numeric values to integers if that is possible without
    loss of precision:

        sage: type(r([1.0]).sage()) == int  # optional - rpy2
        True

        sage: r([1.0, 42.5]).sage()  # optional - rpy2
        [1, 42.5]

    Matrices are converted to sage matrices:

        sage: r('matrix(c(2,4,3,1,5,7), nrow=2, ncol=3)').sage()  # optional - rpy2
        [2 3 5]
        [4 1 7]

    More complex r structures are represented by dictionaries:

        sage: r.summary(1).sage()  # optional - rpy2
        {'DATA': [1, 1, 1, 1, 1, 1],
         '_Names': ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.'],
         '_r_class': ['summaryDefault', 'table']}

        sage: r.options(width="60").sage()  # optional - rpy2
        {'DATA': {'width': 60}, '_Names': 'width'}

    The conversion can handle "not a number", infintiy, imaginary values and
    missing values:

        sage: r(-17).sqrt().sage()  # optional - rpy2
        nan
        sage: r('-17+0i').sqrt().sage()  # optional - rpy2
        4.123105625617661j
        sage: r('NA').sage()  # optional - rpy2
        NA
        sage: inf = r('Inf'); inf.sage()  # optional - rpy2
        inf


    Character Vectors are represented by regular python arrays:

        sage: labs = r.paste('c("X","Y")', '1:10', sep='""'); labs.sage()  # optional - rpy2
        ['X1', 'Y2', 'X3', 'Y4', 'X5', 'Y6', 'X7', 'Y8', 'X9', 'Y10']
    """
    from rpy2.rinterface import SexpVector, ListSexpVector, FloatSexpVector
    from rpy2.robjects.conversion import Converter

    # convert rpy2's representation of r objects to the one sage expects (as defined by the old
    # expect interface)
    cv = Converter('r to sage converter')

    # support rpy version 2 and 3
    try:
        rpy2py = cv.rpy2py
    except AttributeError:
        rpy2py = cv.ri2py

    # fallback
    rpy2py.register(object, lambda obj: obj)

    def float_to_int_if_possible(f):
        # First, round the float to at most 15 significant places.
        # This is what R does by default when using `dput`. It prevents
        # platform-specific fluctuations.
        f = float('%.15g' % f)
        # Preserve the behaviour of the old r parser, e.g. return 1 instead of 1.0
        float_or_int = int(f) if isinstance(f, int) or f.is_integer() else f
        return float_or_int
    rpy2py.register(float, float_to_int_if_possible)

    def list_to_singleton_if_possible(l):
        if len(l) == 1:
            return l[0]
        else:
            return l

    def _vector(vec):
        attrs = vec.list_attrs()
        # Recursive calls have to be made explicitly
        # https://bitbucket.org/rpy2/rpy2/issues/363/custom-converters-are-not-applied
        data = list_to_singleton_if_possible([ rpy2py(val) for val in vec ])
        rclass = list(vec.do_slot('class')) if 'class' in attrs else vec.rclass

        if 'names' in attrs:
            # separate names and values, call rpy2py recursively to convert elements
            names = list_to_singleton_if_possible(list(vec.do_slot('names')))
            return {
                'DATA': data,
                '_Names': names,
                '_r_class': rclass,
            }
        else:
            # if no names are present, convert to a normal list or a single value
            return data
    rpy2py.register(SexpVector, _vector)

    def _matrix(mat):
        if 'dim' in mat.list_attrs():
            try:
                from sage.matrix.constructor import matrix
                dimensions = mat.do_slot("dim")
                if len(dimensions) != 2:
                    raise NotImplementedError("Higher-dimension matrices are currently not supported")
                (nrow, ncol) = dimensions
                # Since R does it the other way round, we assign transposed and
                # then transpose the matrix :)
                m = matrix(ncol, nrow, [rpy2py(i) for i in mat])
                return m.transpose()
            except TypeError:
                pass
        else:
            return _vector(mat)
    rpy2py.register(FloatSexpVector, _matrix)

    def _list_vector(vec):
        # we have a R list (vector of arbitrary elements)
        attrs = vec.list_attrs()
        names = vec.do_slot('names')
        values = [ rpy2py(val) for val in vec ]
        rclass = list(vec.do_slot('class')) if 'class' in attrs else vec.rclass
        data = zip(names, values)
        return {
            'DATA': dict(data),
            '_Names': rpy2py(names),
            # We don't give the rclass here because the old expect interface
            # didn't do that either and we want to maintain compatibility.
        };
    rpy2py.register(ListSexpVector, _list_vector)

    return cv

class R(ExtraTabCompletion, Interface):
    def __init__(self,
                 maxread=None,
                 logfile=None,
                 init_list_length=1024,
                 seed=None):
        """
        An interface to the R interpreter.

        R is a comprehensive collection of methods for statistics,
        modelling, bioinformatics, data analysis and much more.
        For more details, see http://www.r-project.org/about.html

        Resources:

        * http://r-project.org/ provides more information about R.
        * http://rseek.org/ R's own search engine.

        EXAMPLES::

             sage: r.summary(r.c(1,2,3,111,2,3,2,3,2,5,4))  # optional - rpy2
             Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
             1.00    2.00    3.00   12.55    3.50  111.00

        TESTS::

            sage: r == loads(dumps(r))  # optional - rpy2
            True
        """

        Interface.__init__(
                self,
                name = 'r', # The capitalized version of this is used for printing.
        )
        self._seed = seed
        self._initialized = False # done lazily

    def _lazy_init(self):
        """
        Initialize the R interpreter.

        This will set the initial options and implicitly (through
        lazy_import) import rpy2 if it is not already imported.

        Importing rpy2 takes something in the order of hundreds of milliseconds.
        It also takes tens of megabytes of RAM. Since an instance of R is
        assigned to the global variable `r` at sage startup, it is important to
        be as lazy as possible here.
        For some discussion, see https://bitbucket.org/rpy2/rpy2/issues/490.

        Also, importing rpy2 too early (e.g. before numpy) can cause issues with
        the blas implementation that is used.
        For details, see https://bitbucket.org/rpy2/rpy2/issues/491.

        TESTS::

        Initialization happens on eval:

             sage: my_r = R()  # optional - rpy2
             sage: my_r._initialized  # optional - rpy2
             False
             sage: my_r(42) # indirect doctest  # optional - rpy2
             [1] 42
             sage: my_r._initialized  # optional - rpy2
             True

        And on package import:

             sage: my_r = R()  # optional - rpy2
             sage: my_r._initialized  # optional - rpy2
             False
             sage: my_r.library('grid')  # optional - rpy2
             sage: my_r._initialized  # optional - rpy2
             True

        And when fetching help pages:

             sage: my_r = R()  # optional - rpy2
             sage: my_r._initialized  # optional - rpy2
             False
             sage: _ = my_r.help('c')  # optional - rpy2
             sage: my_r._initialized  # optional - rpy2
             True
        """
        if not self._initialized:
            # Set this to True *before* the call to start, since that will call eval() which will in turn call this function.
            # Setting this to True early prevents infinite recursion.
            self._initialized = True
            self._r_to_sage_converter = _setup_r_to_sage_converter()
            self._start()

    def _coerce_impl(self, x, use_special=True):
        """
        TESTS:

        Check conversion of Booleans (:trac:`28705`)::

            sage: repr(r(True)) == r._true_symbol()  # indirect doctest  # optional - rpy2
            True
        """
        # We overwrite _coerce_impl here because r._true_symbol() and
        # r._false_symbol() are output strings that start with "[1] " and thus
        # cannot be used as input
        if isinstance(x, bool):
            return self('TRUE' if x else 'FALSE')
        return super(R, self)._coerce_impl(x, use_special=use_special)

    def set_seed(self, seed=None):
        """
        Set the seed for R interpreter.

        The seed should be an integer.

        EXAMPLES::

            sage: r = R()  # optional - rpy2
            sage: r.set_seed(1)  # optional - rpy2
            1
            sage: r.sample("1:10", 5) # random  # optional - rpy2
            [1] 3 4 5 7 2
        """
        if seed is None:
            seed = self.rand_seed()
        self.eval('set.seed(%d)' % seed)
        self._seed = seed
        return seed

    def _start(self):
        """
        Start up the R interpreter and sets the initial prompt and options.

        This is called the first time the R interface is actually used.

        EXAMPLES::

            sage: r = R()  # optional - rpy2
            sage: r._start()  # optional - rpy2
        """
        # pager needed to replace help view from less to printout
        # option device= is for plotting, is set to x11, NULL would be better?
        self.eval('options(pager="cat",device="png")')
        self.eval('options(repos="%s")'%RRepositoryURL)
        self.eval('options(CRAN="%s")'%RRepositoryURL)

        # don't abort on errors, just raise them!
        # necessary for non-interactive execution
        self.eval('options(error = expression(NULL))')

        # set random seed
        self.set_seed(self._seed)

    def png(self, *args, **kwds):
        """
        Creates an R PNG device.

        This should primarily be used to save an R graphic to a custom file.  Note
        that when using this in the notebook, one must plot in the same cell that
        one creates the device.  See r.plot() documentation for more information
        about plotting via R in Sage.

        These examples won't work on the many platforms where R still gets
        built without graphics support.

        EXAMPLES::

            sage: filename = tmp_filename() + '.png'  # optional - rpy2
            sage: r.png(filename='"%s"'%filename)             # optional -- rgraphics  # optional - rpy2
            NULL
            sage: x = r([1,2,3])  # optional - rpy2
            sage: y = r([4,5,6])  # optional - rpy2
            sage: r.plot(x,y) # This saves to filename, but is not viewable from command line; optional -- rgraphics  # optional - rpy2
            null device
                      1
            sage: import os; os.unlink(filename) # We remove the file for doctesting; optional -- rgraphics  # optional - rpy2

        We want to make sure that we actually can view R graphics, which happens
        differently on different platforms::

            sage: s = r.eval('capabilities("png")') # Should be on Linux and Solaris  # optional - rpy2
            sage: t = r.eval('capabilities("aqua")') # Should be on all supported Mac versions  # optional - rpy2
            sage: "TRUE" in s+t                      # optional -- rgraphics  # optional - rpy2
            True
        """
        #Check to see if R has PNG support
        s = self.eval('capabilities("png")')
        t = self.eval('capabilities("aqua")')
        if "TRUE" not in s+t:
            raise RuntimeError("R was not compiled with PNG support")

        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE:
            self.setwd('"%s"'%os.path.abspath('.'))
        return RFunction(self, 'png')(*args, **kwds)

    def convert_r_list(self, l):
        r"""
        Converts an R list to a Python list.

        EXAMPLES::

            sage: s = 'c(".GlobalEnv", "package:stats", "package:graphics", "package:grDevices", \n"package:utils", "package:datasets", "package:methods", "Autoloads", \n"package:base")'  # optional - rpy2
            sage: r.convert_r_list(s)  # optional - rpy2
            ['.GlobalEnv',
             'package:stats',
             'package:graphics',
             'package:grDevices',
             'package:utils',
             'package:datasets',
             'package:methods',
             'Autoloads',
             'package:base']
        """
        # This function is only kept for legacy reasons. It was used internally
        # in the old expect based interface and for some reason was made part
        # of the public api.
        return self(l).sage()

    def install_packages(self, package_name):
        """
        Install an R package into Sage's R installation.

        EXAMPLES::

            sage: r.install_packages('aaMI')       # not tested  # optional - rpy2
            ...
            R is free software and comes with ABSOLUTELY NO WARRANTY.
            You are welcome to redistribute it under certain conditions.
            Type 'license()' or 'licence()' for distribution details.
            ...
            Please restart Sage in order to use 'aaMI'.
        """
        cmd = """options(repos="%s"); install.packages("%s")"""%(RRepositoryURL, package_name)
        os.system("time echo '%s' | R --vanilla"%cmd)
        print("Please restart Sage in order to use '%s'." % package_name)

    def _repr_(self):
        """
        Return string representation of this R interface.

        EXAMPLES::

            sage: r                 # indirect doctest  # optional - rpy2
            R Interpreter
        """
        return 'R Interpreter'

    def __reduce__(self):
        """
        Used in serializing an R interface.

        EXAMPLES::

            sage: rlr, t = r.__reduce__()  # optional - rpy2
            sage: rlr(*t)  # optional - rpy2
            R Interpreter
        """
        return reduce_load_R, tuple([])

    def __getattr__(self, attrname):
        """
        Called when you get an attribute of the R interface.  This
        manufactures an R function, which is a Python function that
        can then be called with various inputs.

        EXAMPLES::

            sage: c = r.c; c  # optional - rpy2
            c
            sage: type(c)  # optional - rpy2
            <class 'sage.interfaces.r.RFunction'>
        """
        try:
            # First try to get a regular python attribute. This makes it
            # possible to still use attributes like _r_to_sage_converter
            # internally.
            self.__getattribute__(attrname)
        except AttributeError:
            # if there is no such attribute, get the r attribute
            if attrname[:1] == "_":
                raise AttributeError("Attribute {} is not allowed to start with an underscore.".format(attrname))
            return RFunction(self, attrname)


    def _read_in_file_command(self, filename):
        r"""
        Return the R command (as a string) to read in a file named
        filename into the R interpreter.

        EXAMPLES::

            sage: r._read_in_file_command('file.txt')  # optional - rpy2
            'file=file("file.txt",open="r")\nsource(file)'
        """
        return 'file=file("%s",open="r")\nsource(file)'%filename

    def read(self, filename):
        r"""
        Read filename into the R interpreter by calling R's source function on a
        read-only file connection.

        EXAMPLES::

            sage: filename = tmp_filename()  # optional - rpy2
            sage: f = open(filename, 'w')  # optional - rpy2
            sage: _ = f.write('a <- 2+2\n')  # optional - rpy2
            sage: f.close()  # optional - rpy2
            sage: r.read(filename)  # optional - rpy2
            sage: r.get('a')  # optional - rpy2
            '[1] 4'
        """
        self.eval( self._read_in_file_command(filename) )

    def _install_hints(self):
        """
        EXAMPLES::

            sage: print(r._install_hints())  # optional - rpy2
            R is currently installed with Sage.
        """
        return "R is currently installed with Sage.\n"

    def _source(self, s):
        """
        Returns the source code of an R function as a string.

        INPUT:

        - s -- the name of the function as a string

        EXAMPLES::

            sage: print(r._source("c"))  # optional - rpy2
            function (...)  .Primitive("c")
        """
        if s[-2:] == "()":
            s = s[-2:]
        return self.eval('%s'%s)

    def source(self, s):
        """
        Display the R source (if possible) about the function named s.

        INPUT:

        - s -- a string representing the function whose source code you want to see

        OUTPUT: string -- source code

        EXAMPLES::

            sage: print(r.source("c"))  # optional - rpy2
            function (...)  .Primitive("c")
        """
        return self._source(s)

    def version(self):
        """
        Return the version of R currently running.

        OUTPUT: tuple of ints; string

        EXAMPLES::

            sage: r.version() # not tested  # optional - rpy2
            ((3, 0, 1), 'R version 3.0.1 (2013-05-16)')
            sage: rint, rstr = r.version()  # optional - rpy2
            sage: rint[0] >= 3  # optional - rpy2
            True
            sage: rstr.startswith('R version')  # optional - rpy2
            True
        """
        major_re = re.compile(r'^major\s*(\d.*?)$', re.M)
        minor_re = re.compile(r'^minor\s*(\d.*?)$', re.M)
        version_string_re = re.compile(r'^version.string\s*(R.*?)$', re.M)

        s = self.eval('version')

        major = int( major_re.findall(s)[0].strip() )
        minor = tuple(int(i) for i in minor_re.findall(s)[0].strip().split(".") )
        version_string = version_string_re.findall(s)[0].strip()

        return ( (major,) + minor, version_string )

    def library(self, library_name):
        """
        Load the library library_name into the R interpreter.

        This function raises an ImportError if the given library
        is not known.

        INPUT:

        - library_name -- string

        EXAMPLES::

            sage: r.library('grid')  # optional - rpy2
            sage: 'grid' in r.eval('(.packages())')  # optional - rpy2
            True
            sage: r.library('foobar')  # optional - rpy2
            Traceback (most recent call last):
            ...
            ImportError: ...
        """
        self._lazy_init()
        if rpy2_packages.isinstalled(library_name):
            rpy2_packages.importr(library_name)
        else:
            raise ImportError("R library {} not installed".format(library_name))

        try:
            # We need to rebuild keywords!
            del self.__tab_completion
        except AttributeError:
            pass
        self._tab_completion(verbose=False, use_disk_cache=False)

    require = library #overwrites require

    def available_packages(self):
        """
        Returns a list of all available R package names.

        This list is not necessarily sorted.

        OUTPUT: list of strings

        .. note::

            This requires an internet connection. The CRAN server is
            that is checked is defined at the top of sage/interfaces/r.py.

        EXAMPLES::

            sage: ap = r.available_packages()   # optional - internet  # optional - rpy2
            sage: len(ap) > 20                  # optional - internet  # optional - rpy2
            True
        """
        p = self.new('available.packages("%s/src/contrib")'%RRepositoryURL)
        s = str(p).splitlines()[1:]
        v = [x.split()[0].strip("'") for x in s]
        return v
        #The following was more structural, but breaks on my machine.  (stein)
        #p = p._sage_()
        #s = p['_Dim'][0]
        #l = [[p['DATA'][i],p['DATA'][s+1+i]] for i in range(0,s)]
        #return l

    def _object_class(self):
        """
        Return the underlying class of elements of the R interface object.

        OUTPUT: a Python class

        EXAMPLES::

            sage: r._object_class()  # optional - rpy2
            <class 'sage.interfaces.r.RElement'>
        """
        return RElement

    def _true_symbol(self):
        """
        Return the symbol that represents True in R.

        OUTPUT: string

        EXAMPLES::

            sage: r._true_symbol()  # optional - rpy2
            '[1] TRUE'
        """
        # return the string rep of truth, i.e., what the system outputs
        # when you type 1==1.
        return "[1] TRUE"

    def _false_symbol(self):
        """
        Return the symbol that represents True in R.

        OUTPUT: string

        EXAMPLES::

            sage: r._false_symbol()  # optional - rpy2
            '[1] FALSE'
        """
        # return the string rep of false, i.e., what the system outputs
        # when you type 1==2.
        return "[1] FALSE"

    def _equality_symbol(self):
        """
        EXAMPLES::

            sage: r._equality_symbol()  # optional - rpy2
            '=='
        """
        # return the symbol for checking equality, e.g., == or eq.
        return "=="

    # A replacement for rpy2's help.pages that only considers loaded packages
    # (as R's help function does by default). Hopefully upstream will support
    # this in the future: https://bitbucket.org/rpy2/rpy2/issues/498
    def _loaded_package_pages(self, topic):
        # for some reason `except` doesn't work with lazy import, so import this here
        from rpy2.robjects.help import HelpNotFoundError
        self._lazy_init()
        res = list()

        for name in rinterface.baseenv['loadedNamespaces']():
            pack = Package(name)
            try:
                page = pack.fetch(topic)
                res.append(page)
            except HelpNotFoundError:
                pass

        return tuple(res)

    def help(self, command):
        """
        Returns help string for a given command.

        INPUT:
        - command -- a string

        OUTPUT: HelpExpression -- a subclass of string whose __repr__ method is __str__, so it prints nicely

        EXAMPLES::

            sage: r.help('c')  # optional - rpy2
            title
            -----
            <BLANKLINE>
            Combine Values into a Vector or List
            <BLANKLINE>
            name
            ----
            <BLANKLINE>
            c
            ...
        """
        # This is looking for the topic in all existing namespaces.
        # Theoretically, there may be multiple options. We ignore that.
        pages_for_topic = self._loaded_package_pages(command)
        if len(pages_for_topic) == 0:
            raise ValueError("There is no help page for the given topic")

        s = pages_for_topic[0].to_docstring()
        return HelpExpression(s)

    def _assign_symbol(self):
        """
        Return the symbol used in R for assignment, which is ' <- '.

        OUTPUT: string

        EXAMPLES::

            sage: r._assign_symbol()  # optional - rpy2
            ' <- '
        """
        return " <- "

    def _left_list_delim(self):
        """
        Return the left delimiter for lists in R, which is 'c('

        OUTPUT: string

        EXAMPLES::

            sage: r._left_list_delim()  # optional - rpy2
            'c('
        """
        return "c("

    def _right_list_delim(self):
        """
        Return the right delimiter for lists in R, which is 'c('

        OUTPUT: string

        EXAMPLES::

            sage: r._right_list_delim()  # optional - rpy2
            ')'
        """
        return ")"

    def console(self):
        """
        Runs the R console as a separate new R process.

        EXAMPLES::

            sage: r.console()                    # not tested  # optional - rpy2
                R version 2.6.1 (2007-11-26)
                Copyright (C) 2007 The R Foundation for Statistical Computing
                ISBN 3-900051-07-0
                ...
        """
        r_console()

    def function_call(self, function, args=None, kwds=None):
        """
        Return the result of calling an R function, with given args and keyword args.

        OUTPUT: RElement -- an object in R

        EXAMPLES::

            sage: r.function_call('length', args=[ [1,2,3] ])  # optional - rpy2
            [1] 3
        """
        args, kwds = self._convert_args_kwds(args, kwds)
        self._check_valid_function_name(function)
        return self.new("%s(%s)"%(function, ",".join([s.name() for s in args] +
                                                     [self._sage_to_r_name(key)+'='+kwds[key].name() for key in kwds ] )))

    def call(self, function_name, *args, **kwds):
        r"""
        This is an alias for :meth:`function_call`.

        EXAMPLES::

            sage: r.call('length', [1,2,3])  # optional - rpy2
            [1] 3
        """
        return self.function_call(function_name, args=args, kwds=kwds)

    def _an_element_impl(self):
        """
        Returns an element belonging to the R interpreter.  This is used
        behind the scenes when doing things like comparisons, etc.

        OUTPUT: RElement -- an R element.

        EXAMPLES::

            sage: r._an_element_impl()  # optional - rpy2
            [1] 0
            sage: type(_)  # optional - rpy2
            <class 'sage.interfaces.r.RElement'>
        """
        return self(0)

    def set(self, var, value):
        """
        Set the variable var in R to what the string value evaluates to in R.

        INPUT:

        - var -- a string
        - value -- a string

        EXAMPLES::

            sage: r.set('a', '2 + 3')  # optional - rpy2
            sage: r.get('a')  # optional - rpy2
            '[1] 5'

        """
        cmd = '%s <- %s'%(var,value)
        out = self.eval(cmd)

    def get(self, var):
        """
        Returns the string representation of the variable var.

        INPUT:

        - var -- a string

        OUTPUT: string

        EXAMPLES::

            sage: r.set('a', 2)  # optional - rpy2
            sage: r.get('a')  # optional - rpy2
            '[1] 2'
        """
        return self.eval('%s'%var)

    def na(self):
        """
        Returns the NA in R.

        OUTPUT: RElement -- an element of R

        EXAMPLES::

            sage: r.na()  # optional - rpy2
            [1] NA
        """
        return self('NA')

    def completions(self, s):
        """
        Return all commands names that complete the command starting with the
        string s.   This is like typing s[Ctrl-T] in the R interpreter.

        INPUT:

        - s -- string

        OUTPUT: list -- a list of strings

        EXAMPLES::

            sage: dummy = r._tab_completion(use_disk_cache=False)    #clean doctest  # optional - rpy2
            sage: 'testInheritedMethods' in r.completions('tes')  # optional - rpy2
            True
        """
        return [name for name in self._tab_completion() if name[:len(s)] == s]

    def _commands(self):
        """
        Return list of all commands defined in R.

        OUTPUT: list -- a sorted list of strings

        EXAMPLES::

            sage: l = r._commands()  # optional - rpy2
            sage: 'AIC' in l  # optional - rpy2
            True
            sage: len(l) > 200  # optional - rpy2
            True
        """
        v = RBaseCommands

        ll = self('search()')._sage_() # loaded libs

        for lib in ll:
            if lib in RFilteredPackages:
                continue

            if lib.find("package:") != 0:
                continue #only packages

            raw = self('objects("%s")'%lib)._sage_()

            #TODO are there others? many of them are shortcuts or
            #should be done on another level, like selections in lists
            #instead of calling obj.[[( fun-args) or other crazy stuff like that

            #TODO further filtering, check if strings are now
            #really functions, in R: exists(s, mode = "function"))
            # (apply to vector with sapply(vec,func))

            #filter only python compatible identifiers
            valid = re.compile('[^a-zA-Z0-9_]+')
            raw = [x for x in raw if valid.search(x) is None]
            v += raw

        v.sort()
        return v

    def _tab_completion(self, verbose=True, use_disk_cache=True):
        """
        Return list of all R functions.

        INPUT:

        - verbose -- bool (default: True); if True, display debugging information
        - use_disk_cache -- bool (default: True); if True, use the disk cache of
          tab completions to save time.

        OUTPUT: list -- list of string

        EXAMPLES::

            sage: t = r._tab_completion(verbose=False)  # optional - rpy2
            sage: len(t) > 200  # optional - rpy2
            True
        """
        try:
            return self.__tab_completion
        except AttributeError:
            import sage.misc.persist
            if use_disk_cache:
                try:
                    self.__tab_completion = sage.misc.persist.load(COMMANDS_CACHE)
                    return self.__tab_completion
                except IOError:
                    pass
            if verbose and use_disk_cache:
                print("\nBuilding R command completion list (this takes")
                print("a few seconds only the first time you do it).")
                print("To force rebuild later, delete %s." % COMMANDS_CACHE)
            v = self._commands()
            self.__tab_completion = v
            if len(v) > 200 and use_disk_cache:
                sage.misc.persist.save(v, COMMANDS_CACHE)
            return v

    def plot(self, *args, **kwds):
        """
        The R plot function.  Type r.help('plot') for much more extensive
        documentation about this function.  See also below for a brief
        introduction to more plotting with R.

        If one simply wants to view an R graphic, using this function is
        is sufficient (because it calls dev.off() to turn off the device).

        However, if one wants to save the graphic to a specific file, it
        should be used as in the example below to write the output.

        EXAMPLES:

        This example saves a plot to the standard R output, usually
        a filename like ``Rplot001.png`` - from the command line, in
        the current directory, and in the cell directory in the notebook::

            sage: d=r.setwd('"%s"'%SAGE_TMP)    # for doctesting only; ignore if you are trying this  # optional - rpy2
            sage: r.plot("1:10")                # optional -- rgraphics  # optional - rpy2
            null device
                      1

        To save to a specific file name, one should use :meth:`png` to set
        the output device to that file.  If this is done in the notebook, it
        must be done in the same cell as the plot itself::

            sage: filename = tmp_filename() + '.png'  # optional - rpy2
            sage: r.png(filename='"%s"'%filename) # Note the double quotes in single quotes!; optional -- rgraphics  # optional - rpy2
            NULL
            sage: x = r([1,2,3])  # optional - rpy2
            sage: y = r([4,5,6])  # optional - rpy2
            sage: r.plot(x,y)         # optional -- rgraphics  # optional - rpy2
            null device
                      1
            sage: import os; os.unlink(filename) # For doctesting, we remove the file; optional -- rgraphics  # optional - rpy2

        Please note that for more extensive use of R's plotting
        capabilities (such as the lattices package), it is advisable
        to either use an interactive plotting device or to use the
        notebook.  The following examples are not tested, because they
        differ depending on operating system::

            sage: r.X11() # not tested - opens interactive device on systems with X11 support  # optional - rpy2
            sage: r.quartz() # not tested - opens interactive device on OSX  # optional - rpy2
            sage: r.hist("rnorm(100)") # not tested - makes a plot  # optional - rpy2
            sage: r.library("lattice") # not tested - loads R lattice plotting package  # optional - rpy2
            sage: r.histogram(x = "~ wt | cyl", data="mtcars") # not tested - makes a lattice plot  # optional - rpy2
            sage: r.dev_off() # not tested, turns off the interactive viewer  # optional - rpy2

        In the notebook, one can use r.png() to open the device, but
        would need to use the following since R lattice graphics do
        not automatically print away from the command line::

            sage: filename = tmp_filename() + '.png' # Not needed in notebook, used for doctesting  # optional - rpy2
            sage: r.png(filename='"%s"'%filename) # filename not needed in notebook, used for doctesting; optional -- rgraphics  # optional - rpy2
            NULL
            sage: r.library("lattice")  # optional - rpy2
            sage: r("print(histogram(~wt | cyl, data=mtcars))") # plot should appear; optional -- rgraphics  # optional - rpy2
            sage: import os; os.unlink(filename) # We remove the file for doctesting, not needed in notebook; optional -- rgraphics  # optional - rpy2
        """
        # We have to define this to override the plot function defined in the
        # superclass.
        from sage.server.support import EMBEDDED_MODE
        if EMBEDDED_MODE:
            self.setwd('"%s"'%os.path.abspath('.'))
        RFunction(self, 'plot')(*args, **kwds)
        return RFunction(self, 'dev.off')()

    def eval(self, code, *args, **kwds):
        """
        Evaluates a command inside the R interpreter and returns the output
        as a string.

        EXAMPLES::

            sage: r.eval('1+1')  # optional - rpy2
            '[1] 2'
        """
        self._lazy_init()
        return str(robjects.r(code)).rstrip()


    def _r_to_sage_name(self, s):
        """
        Returns a Sage/Python identifier from an R one.  This involves
        replacing periods with underscores, <- with __, and prepending
        _ in front of Python keywords.

        INPUT:

        - s -- a string

        OUTPUT: a string

        EXAMPLES::

            sage: f = r._r_to_sage_name  # optional - rpy2
            sage: f('t.test')  # optional - rpy2
            't_test'
            sage: f('attr<-')  # optional - rpy2
            'attr__'
            sage: f('parent.env<-')  # optional - rpy2
            'parent_env__'
            sage: f('class')  # optional - rpy2
            'class_'
        """
        from keyword import iskeyword
        s = s.replace('.', '_')
        s = s.replace('<-', '__')
        if iskeyword(s):
            s += '_'
        return s

    def _sage_to_r_name(self, s):
        r"""
        The reverse of :meth:`_r_to_sage_name`.  See the docs for that method.

        EXAMPLES::

            sage: f = r._sage_to_r_name  # optional - rpy2
            sage: f('t_test')  # optional - rpy2
            't.test'
            sage: f('attr__')  # optional - rpy2
            'attr<-'
            sage: f('parent_env__')  # optional - rpy2
            'parent.env<-'
            sage: r._r_to_sage_name(f('parent_env__'))  # optional - rpy2
            'parent_env__'
            sage: f('class_')  # optional - rpy2
            'class'
        """
        if len(s) > 1 and s[-2:] == "__":
            s = s[:-2] + '<-'
        if len(s) > 0 and s[-1] == '_':
            s = s[:-1]
        s = s.replace('_', '.')
        return s

    def __getitem__(self, s):
        """
        Returns the RFunction with name s.

        INPUT:

        - s -- a string
        OUTPUT: RFunction -- the R function that in R has name s

        EXAMPLES::

            sage: r['as.data.frame']  # optional - rpy2
            as.data.frame
            sage: r['print']  # optional - rpy2
            print
        """
        return RFunction(self, s, r_name=True)

    def chdir(self, dir):
        """
        Changes the working directory to ``dir``

        INPUT:

        - ``dir`` -- the directory to change to.

        EXAMPLES::

            sage: import tempfile  # optional - rpy2
            sage: tmpdir = tempfile.mkdtemp()  # optional - rpy2
            sage: r.chdir(tmpdir)  # optional - rpy2

        Check that ``tmpdir`` and ``r.getwd()`` refer to the same
        directory.  We need to use ``realpath()`` in case ``$TMPDIR``
        (by default ``/tmp``) is a symbolic link (see :trac:`10264`).

        ::

            sage: os.path.realpath(tmpdir) == sageobj(r.getwd())  # known bug (trac #9970)  # optional - rpy2
            True
        """
        self.execute('setwd(%r)' % dir)


@instancedoc
class RElement(ExtraTabCompletion, InterfaceElement):

    def _tab_completion(self):
        """
        Return a list of all methods of this object.

        .. note::

            Currently returns all R commands.

        EXAMPLES::

            sage: a = r([1,2,3])  # optional - rpy2
            sage: t = a._tab_completion()  # optional - rpy2
            sage: len(t) > 200  # optional - rpy2
            True
        """
        # TODO: rewrite it, just take methods(class=class(self))
        return self.parent()._tab_completion()

    def tilde(self, x):
        """
        The tilde regression operator in R.

        EXAMPLES::

            sage: x = r([1,2,3,4,5])  # optional - rpy2
            sage: y = r([3,5,7,9,11])  # optional - rpy2
            sage: a = r.lm( y.tilde(x) ) # lm( y ~ x )  # optional - rpy2
            sage: d = a._sage_()  # optional - rpy2
            sage: d['DATA']['coefficients']['DATA'][1]  # optional - rpy2
            2
        """
        par = self.parent()
        rx = par(x)
        return par.new("%s ~ %s" % (self.name(), rx.name()))

    stat_model = tilde

    def is_string(self):
        """
        Tell whether this element is a string.

        EXAMPLES::

            sage: r('"abc"').is_string()  # optional - rpy2
            True
            sage: r([1,2,3]).is_string()  # optional - rpy2
            False

        """
        return isinstance(self.sage(), str)

    def __len__(self):
        """
        Return the length of this object.

        OUTPUT: integer

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: len(x)  # optional - rpy2
            5
        """
        return self.parent()('length(%s)'%self.name()).sage()

    def __getattr__(self, attrname):
        """
        Return attribute of this object, which is an R function with this object
        as the first input.

        INPUT:

        - attrname -- string

        OUTPUT: RFunctionElement

        EXAMPLES::

            sage: x = r([1,2,3])  # optional - rpy2
            sage: length = x.length  # optional - rpy2
            sage: type(length)  # optional - rpy2
            <class 'sage.interfaces.r.RFunctionElement'>
            sage: length()  # optional - rpy2
            [1] 3
        """
        try:
            # First try to get a regular python attribute. This makes it
            # possible to still use attributes like _r_to_sage_converter
            # internally.
            self.__getattribute__(attrname)
        except AttributeError:
            self._check_valid()
            if attrname[:1] == "_":
                raise AttributeError("Attribute {} is not allowed to start with an underscore.".format(attrname))
            return RFunctionElement(self, attrname)

    def __getitem__(self, n):
        """
        Return element(s) of self.

        INPUT:

        - n -- an integer, a tuple, a string that makes sense to R, or an RElement

        OUTPUT: RElement

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x[0]  # optional - rpy2
            numeric(0)
            sage: x[1]  # optional - rpy2
            [1] 10.4
            sage: x[-1]  # optional - rpy2
            [1] 5.6  3.1  6.4 21.7
            sage: x[-2]  # optional - rpy2
            [1] 10.4  3.1  6.4 21.7
            sage: x[-3]  # optional - rpy2
            [1] 10.4  5.6  6.4 21.7
            sage: x['c(2,3)']  # optional - rpy2
            [1]  5.6 3.1
            sage: key = r.c(2,3)  # optional - rpy2
            sage: x[key]  # optional - rpy2
            [1]  5.6 3.1
            sage: m = r.array('1:3',r.c(2,4,2))  # optional - rpy2
            sage: m  # optional - rpy2
            , , 1
                 [,1] [,2] [,3] [,4]
            [1,]    1    3    2    1
            [2,]    2    1    3    2
            , , 2
                 [,1] [,2] [,3] [,4]
            [1,]    3    2    1    3
            [2,]    1    3    2    1
            sage: m[1,2,2]  # optional - rpy2
            [1] 2
            sage: m[1,r.c(1,2),1]  # optional - rpy2
            [1] 1 3
        """
        P = self._check_valid()
        if isinstance(n, str):
            n = n.replace('self', self._name)
            return P.new('%s[%s]'%(self._name, n))
        elif parent(n) is P:  # the key is RElement itself
            return P.new('%s[%s]'%(self._name, n.name()))
        elif not isinstance(n,tuple):
            return P.new('%s[%s]'%(self._name, n))
        else:
            L = []
            for i in range(len(n)):
                if parent(n[i]) is P:
                    L.append(n[i].name())
                else:
                    L.append(str(n[i]))
            return P.new('%s[%s]'%(self._name, ','.join(L)))

    def __bool__(self):
        """
        Implements bool(self).

        .. note::

            bool(self) will only return True if self == 0 contains a FALSE in its representation.

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: bool(x)  # optional - rpy2
            True
            sage: y = r([0,0,0,0])  # optional - rpy2
            sage: bool(y)  # optional - rpy2
            False
            sage: bool(r(0))  # optional - rpy2
            False
            sage: bool(r(1))  # optional - rpy2
            True
        """
        return "FALSE" in repr(self==0)

    __nonzero__ = __bool__

    def _comparison(self, other, symbol):
        """
        Used to implement comparison of two objects.

        INPUT:

        - other -- RElement
        - symbol -- string

        OUTPUT: RElement -- output is an R element; not a bool!

        TESTS::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x._comparison(10.4, "==")  # optional - rpy2
            [1] TRUE FALSE FALSE FALSE FALSE
        """
        P = self.parent()
        other = P(other)
        return P('%s %s %s'%(self.name(), symbol, other.name()))

    def __eq__(self, other):
        """
        Equality testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x == 10.4  # optional - rpy2
            [1] TRUE FALSE FALSE FALSE FALSE
        """
        return self._comparison(other, "==")

    def __lt__(self, other):
        """
        Less than testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x < 7  # optional - rpy2
            [1] FALSE  TRUE  TRUE  TRUE FALSE
        """
        return self._comparison(other, "<")

    def __gt__(self, other):
        """
        Greater than testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES:

        Notice that comparison is term by term and returns an R element. ::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x > 8  # optional - rpy2
            [1] TRUE FALSE FALSE FALSE  TRUE
        """
        return self._comparison(other, ">")

    def __le__(self, other):
        """
        Less than or equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x <= 10.4  # optional - rpy2
            [1] TRUE  TRUE  TRUE  TRUE FALSE
        """
        return self._comparison(other, "<=")

    def __ge__(self, other):
        """
        Greater than or equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x >= 10.4  # optional - rpy2
            [1] TRUE FALSE FALSE FALSE  TRUE
        """
        return self._comparison(other, ">=")

    def __ne__(self, other):
        """
        Not equal testing term by term.

        INPUT:

        - other -- RElement

        OUTPUT: RElement -- an R element; not a bool!

        EXAMPLES::

            sage: x = r([10.4,5.6,3.1,6.4,21.7])  # optional - rpy2
            sage: x != 10.4  # optional - rpy2
            [1] FALSE  TRUE  TRUE  TRUE  TRUE

        """
        return self._comparison(other, "!=")

    def dot_product(self, other):
        """
        Implements the notation self . other.

        INPUT:

        - self, other -- R elements

        OUTPUT: R element

        EXAMPLES::

            sage: c = r.c(1,2,3,4)  # optional - rpy2
            sage: c.dot_product(c.t())  # optional - rpy2
                 [,1] [,2] [,3] [,4]
            [1,]    1    2    3    4
            [2,]    2    4    6    8
            [3,]    3    6    9   12
            [4,]    4    8   12   16

            sage: v = r([3,-1,8])  # optional - rpy2
            sage: v.dot_product(v)  # optional - rpy2
                 [,1]
            [1,]   74
        """
        P = self._check_valid()
        Q = P(other)
        # the R operator is %*% for matrix multiplication
        return P('%s %%*%% %s'%(self.name(), Q.name()))

    def _sage_(self):
        r"""
        Returns Sage representation of the R object.

        R objects are basic C structures, of different kind, that can
        be stacked together.  This is similar to Python lists with
        variable objects, including lists of lists.  If R lists have
        names, they are translated to a Python dictionary, with anonymous
        list entries called ``#{number}``.

        OUTPUT: object -- Python object

        EXAMPLES::

            sage: rs = r.summary(r.c(1,4,3,4,3,2,5,1))  # optional - rpy2
            sage: d = rs._sage_()  # optional - rpy2
            sage: sorted(d.items())  # optional - rpy2
            [('DATA', [1, 1.75, 3, 2.875, 4, 5]),
             ('_Names', ['Min.', '1st Qu.', 'Median', 'Mean', '3rd Qu.', 'Max.']),
             ('_r_class', ['summaryDefault', 'table'])]
        """
        self._check_valid()
        P = self.parent()

        with localconverter(P._r_to_sage_converter) as cv:
            parsed = robjects.r(self.name())
            return parsed


    def _latex_(self):
        r"""
        Return LaTeX representation of this R object.

        This calls the ``latex`` command in R.

        OUTPUT: a latex expression (basically a string)

        EXAMPLES::

            sage: latex(r(2))  # optional - Hmisc (R package)  # optional - rpy2
            2
        """
        from sage.misc.latex import LatexExpr
        self._check_valid()
        P = self.parent()
        # latex is in Hmisc, this is currently not part of Sage's R!!!
        try:
            P.library('Hmisc')
        except ImportError:
            raise RuntimeError("The R package 'Hmisc' is required for R to LaTeX conversion, but it is not available.")
        return LatexExpr(P.eval('latex(%s, file="");' % self.name()))


@instancedoc
class RFunctionElement(InterfaceFunctionElement):
    def __reduce__(self):
        """
        EXAMPLES::

            sage: a = r([1,2,3])  # optional - rpy2
            sage: a.mean  # optional - rpy2
            mean
            sage: dumps(a.mean)  # optional - rpy2
            Traceback (most recent call last):
            ...
            NotImplementedError: pickling of R element methods is not yet supported
        """
        raise NotImplementedError("pickling of R element methods is not yet supported")

    def _instancedoc_(self):
        """
        Returns the help for self as a string.

        EXAMPLES::

            sage: a = r([1,2,3])  # optional - rpy2
            sage: length = a.length  # optional - rpy2
            sage: print(length.__doc__)  # optional - rpy2
            title
            -----
            <BLANKLINE>
            Length of an Object
            <BLANKLINE>
            name
            ----
            <BLANKLINE>
            length
            ...
        """
        M = self._obj.parent()
        return M.help(self._name)

    def _sage_src_(self):
        """
        Returns the source code of self.

        EXAMPLES::

            sage: a = r([1,2,3])  # optional - rpy2
            sage: length = a.length  # optional - rpy2
            sage: print(length._sage_src_())  # optional - rpy2
            function (x)  .Primitive("length")
        """
        M = self._obj.parent()
        return M.source(self._name)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: a = r([1,2,3])  # optional - rpy2
            sage: length = a.length  # optional - rpy2
            sage: length()  # optional - rpy2
            [1] 3
        """
        return self._obj.parent().function_call(self._name, args=[self._obj] + list(args), kwds=kwds)


@instancedoc
class RFunction(InterfaceFunction):
    def __init__(self, parent, name, r_name=None):
        """
        A Function in the R interface.

        INPUT:

        - parent -- the R interface
        - name -- the name of the function for Python
        - r_name -- the name of the function in R itself (which can have dots in it)

        EXAMPLES::

            sage: length = r.length  # optional - rpy2
            sage: type(length)  # optional - rpy2
            <class 'sage.interfaces.r.RFunction'>
            sage: loads(dumps(length))  # optional - rpy2
            length
        """
        self._parent = parent
        if r_name:
            self._name = name
        else:
            self._name = parent._sage_to_r_name(name)

    def __eq__(self, other):
        """
        EXAMPLES::

            sage: r.mean == loads(dumps(r.mean))  # optional - rpy2
            True
            sage: r.mean == r.lr  # optional - rpy2
            False
        """
        return (isinstance(other, RFunction) and
            self._name == other._name)

    def __ne__(self, other):
        """
        EXAMPLES::

            sage: r.mean != loads(dumps(r.mean))  # optional - rpy2
            False
            sage: r.mean != r.lr  # optional - rpy2
            True
        """
        return not (self == other)

    def _instancedoc_(self):
        """
        Returns the help for self.

        EXAMPLES::

            sage: length = r.length  # optional - rpy2
            sage: print(length.__doc__)  # optional - rpy2
            title
            -----
            <BLANKLINE>
            Length of an Object
            <BLANKLINE>
            name
            ----
            <BLANKLINE>
            length
            ...
        """
        M = self._parent
        return M.help(self._name)

    def _sage_src_(self):
        """
        Returns the source of self.

        EXAMPLES::

            sage: length = r.length  # optional - rpy2
            sage: print(length._sage_src_())  # optional - rpy2
            function (x)  .Primitive("length")

        """
        M = self._parent
        return M.source(self._name)

    def __call__(self, *args, **kwds):
        """
        EXAMPLES::

            sage: length = r.length  # optional - rpy2
            sage: length([1,2,3])  # optional - rpy2
            [1] 3
        """
        return self._parent.function_call(self._name, args=list(args), kwds=kwds)

def is_RElement(x):
    """
    Return True if x is an element in an R interface.

    INPUT:

    - x -- object

    OUTPUT: bool

    EXAMPLES::

        sage: from sage.interfaces.r import is_RElement  # optional - rpy2
        sage: is_RElement(2)  # optional - rpy2
        False
        sage: is_RElement(r(2))  # optional - rpy2
        True
    """
    return isinstance(x, RElement)

# An instance of R
r = R()

def reduce_load_R():
    """
    Used for reconstructing a copy of the R interpreter from a pickle.

    EXAMPLES::

        sage: from sage.interfaces.r import reduce_load_R  # optional - rpy2
        sage: reduce_load_R()  # optional - rpy2
        R Interpreter
    """
    return r

import os
def r_console():
    """
    Spawn a new R command-line session.

    EXAMPLES::

        sage: r.console()                    # not tested  # optional - rpy2
            R version 2.6.1 (2007-11-26)
            Copyright (C) 2007 The R Foundation for Statistical Computing
            ISBN 3-900051-07-0
            ...
    """
    from sage.repl.rich_output.display_manager import get_display_manager
    if not get_display_manager().is_in_terminal():
        raise RuntimeError('Can use the console only in the terminal. Try %%r magics instead.')
    # This will only spawn local processes
    os.system('R --vanilla')

def r_version():
    """
    Return the R version.

    EXAMPLES::

        sage: r_version() # not tested  # optional - rpy2
        ((3, 0, 1), 'R version 3.0.1 (2013-05-16)')
        sage: rint, rstr = r_version()  # optional - rpy2
        sage: rint[0] >= 3  # optional - rpy2
        True
        sage: rstr.startswith('R version')  # optional - rpy2
        True
    """
    return r.version()

class HelpExpression(str):
    """
    Used to improve printing of output of r.help.
    """
    def __repr__(self):
        """
        Return string representation of self.

        OUTPUT: string

        EXAMPLES::

            sage: a = sage.interfaces.r.HelpExpression("This\nis\nR!")  # optional - rpy2
            sage: type(a)  # optional - rpy2
            <class 'sage.interfaces.r.HelpExpression'>
            sage: a  # optional - rpy2
            This
            is
            R!
        """
        return str(self)

