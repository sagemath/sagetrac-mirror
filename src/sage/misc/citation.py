r"""
Citation management.

Facilities to track (internal and external) components used by specific computations.  This modules provides two different ways, to access citations: A eval like command and a context manager.

EXAMPLES:

The context manager ``citations`` is set up by surounding computations
by ``with citatitons():``.  The compuations that we want to track may
contain more than one statement.  They will be evaluated in as if the
same code was provided without the citations statement.

::

    sage: with citations():
    ...       var('y')
    ...       integrate(y^2, y, 0, 1)
    y
    1/3
    The computation used the following components.
    Access them as a list by calling latest_citations().
        ginac, Maxima

As the statement tells us, we can access the citation items as a list
of :class:`~sage.misc.citation_items.citation_item.CitationItem`'s.
Currently, this list has no further purpose, but in the future it will
allow for automatically generating Bibtex code.

::

    sage: latest_citations()
    [ginac, Maxima]
    sage: latest_citations().print_latex_citation()
    Print BibTeX code by calling method print_bibtex().
    \cite{software-ginac, software-maxima}

For certain purposes, it is desirable to only record citations of some
parts of the computation.  This can be achieved by passing a list,
which will be extended by citation items for invoked components.

::

    sage: record = CitationRecord()
    sage: with citations(record):
    ...       K = QuadraticField(-3, 'a')
    sage: record
    [GAP, GMP, MPFI, MPFR, NTL]
    sage: record.print_latex_citation()
    Print BibTeX code by calling method print_bibtex().
    \cite{software-gap, software-gmp, software-mpfi, software-mpfr, software-ntl}
    sage: with citations(record):
    ...       integrate(y^2, y, 0, 1)
    1/3
    sage: record
    [GAP, GMP, MPFI, MPFR, NTL, ginac, Maxima]

As an alternative way to access citations, we provide an eval like statement.

::

    sage: eval_citations("integrate(y^2, y, 0, 10)")
    [ginac, Maxima]

AUTHORS:

- Martin Westerholt-Raum (martin@raum-brothers.eu)
"""

###########################################################################
#       Copyright (C) 2014 Martin Raum <martin@raum-brothers.eu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

from contextlib import contextmanager
from sage.misc.all import tmp_filename
from sage.misc.citation_items.all import CitationRecord
import os, re, sys


_latest_citations = CitationRecord()

def latest_citations():
    r"""
    A list of citations the latest compuation.  See :meth:`citations`.

    EXAMPLE::

        sage: a = var('a')
        sage: with citations():
        ...       h = ((a+1)^2).expand()
        The computation used the following components.
        Access them as a list by calling latest_citations().
            ginac, GMP
        sage: latest_citations()
        [ginac, GMP]
    """
    return _latest_citations

@contextmanager
def citations(record = None):
    r"""
    Print or append to ``record`` a list of the systems used in
    running block of statements.  Note that the results can sometimes
    include systems that did not actually contribute to the
    computation. Due to caching , it could miss some dependencies as
    well.

    INPUT:

    - ``record`` -- A list or ``None``.

    EXAMPLES::

        sage: R.<x,y,z> = QQ[]
        sage: I = R.ideal(x^2+y^2, z^2+y)
        sage: with citations():
        ...       I.primary_decomposition()
        [Ideal (z^2 + y, x^2 + y^2) of Multivariate Polynomial Ring in x, y, z over Rational Field]
        The computation used the following components.
        Access them as a list by calling latest_citations().
            Macaulay2, Singular
        sage: latest_citations()
        [Macaulay2, Singular]
    """
    import warnings
    import cProfile, pstats

    from sage.misc.citation_items.all import citation_items


    cprofiler = cProfile.Profile()

    cprofiler.enable()

    yield

    cprofiler.disable()
    cprofiler_calls = map(_cprofile_stat_to_function_string,
                          pstats.Stats(cprofiler).stats.keys())


    #Remove trivial functions
    bad_res = map(re.compile,
                  [r'is_.*Element'])
    calls = [c for c in cprofiler_calls
             if all(r.match(c) is None for r in bad_res)]


    #Check to see which citations appear in the profiled run
    called_items = [item for item in citation_items
                    if any(r.match(c) is not None for r in item.re() for c in calls)]

    if record is None:
        import sage
        sage.misc.citation._latest_citations = CitationRecord(called_items)

        print "The computation used the following components."
        print "Access them as a list by calling latest_citations()."
        print "    " + ", ".join(map(repr, called_items))
    else:
        record.extend(called_items)

def eval_citations(cmd, locals = None):
    r"""
    Returns a list of the systems used in running the command
    ``cmd``.  Note that the results can sometimes include systems
    that did not actually contribute to the computation. Due
    to caching , it could miss some dependencies as well.

    INPUT:

    - ``cmd`` - a string to run

    - ``locals`` - a dictionary of locals.

    OUPUT:

    A list of citation items

    EXAMPLES::

        sage: s = eval_citations( "integrate(x^2, x)" ); #priming coercion model
        sage: eval_citations('integrate(x^2, x)')
        [ginac, Maxima]
    """
    import inspect
    from sage.misc.all import sage_eval

    if locals is None:
        locals = inspect.stack()[1][0].f_globals

    record = list()
    with citations(record):
        sage_eval(cmd, locals=locals)
    return record

def get_systems(cmd):
    """
    DEPRECATED:

    Use :meth:`eval_citations` instead.

    Returns a list of the systems used in running the command
    cmd.  Note that the results can sometimes include systems
    that did not actually contribute to the computation. Due
    to caching, it could miss some dependencies as well.

    INPUT:

    - ``cmd`` - a string to run

    EXAMPLES::

        sage: from sage.misc.citation import get_systems
        sage: s = get_systems('integrate(x^2, x)')
        doctest:...DeprecationWarning: call eval_citations instead of get_systems
        See http://trac.sagemath.org/16777 for details.
    """
    from sage.misc.superseded import deprecation
    deprecation(16777, 'call eval_citations instead of get_systems')

    import inspect

    return eval_citations(cmd, locals=inspect.stack()[1][0].f_globals)

def _cprofile_stat_to_function_string(stat_key):
    r"""
    Parse the profiling data collected by ``cProfile`` and ``pstat``.

    INPUT:

    - ``stat_key`` -- A triple ``(string, number, string)``, which
                      occurs as a key in the dictionary ``pstats.Stats.stats``.

    OUTPUT:

    A string.

    EXAMPLES::

        sage: from sage.misc.citation import _cprofile_stat_to_function_string

    Modules accessed as site packages are reconstructed using file
    path information.

    ::

        sage: _cprofile_stat_to_function_string( ("/home/.../sage/local/lib/python2.7/site-packages/sage/rings/number_field/number_field.py", 1178, "_element_constructor_") )
        'sage.rings.number_field.number_field._element_constructor_'

    Builtin methods are reconstructed from the object clauses, if possible.

    ::

        sage: _cprofile_stat_to_function_string( ('~', 0, "<method 'base' of 'sage.structure.category_object.CategoryObject' objects>") )
        'sage.structure.category_object.CategoryObject.base'
        sage: _cprofile_stat_to_function_string( ('~', 0, "<posix.WIFEXITED>") )
        'posix.WIFEXITED'
    """
    if stat_key[0] == '~':
        if stat_key[2].startswith("<method "):
            object_start = stat_key[2].find("of '") + len("of '")
            object_end = stat_key[2].find("'", object_start)
            module_part = stat_key[2][object_start:object_end]

            function_start = stat_key[2].find("method '") + len("method '")
            function_end = stat_key[2].find("'", function_start)
            function_part = stat_key[2][function_start:function_end]
        else:
            module_part = None
            function_part = stat_key[2][1:-1]
    else:
        module_part = (stat_key[0].split("site-packages/")[-1]
                       .split('.')[0]
                       .replace("/", "."))
        function_part = stat_key[2]

    if module_part:
        return module_part + "." + function_part
    else:
        return function_part
