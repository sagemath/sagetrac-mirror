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
        GMP, ginac, Maxima

As the statement tells us, we can access the citation items as a list
of :class:`~sage.misc.citation_items.citation_item.CitationItem`'s.
Currently, this list has no further purpose, but in the future it will
allow for automatically generating Bibtex code.

::

    sage: latest_citations()
    [GMP, ginac, Maxima]
    sage: latest_citations().print_latex_citation()
    Print BibTeX code by calling method print_bibtex().
    \cite{software-gmp, software-ginac, software-maxima}

For certain purposes, it is desirable to only record citations of some
parts of the computation.  This can be achieved by passing a list,
which will be extended by citation items for invoked components.

::

    sage: record = CitationRecord()
    sage: with citations(record):
    ...       K = QuadraticField(-3, 'a')
    sage: record
    [GMP, MPFR, MPFI, NTL]
    sage: record.print_latex_citation()
    Print BibTeX code by calling method print_bibtex().
    \cite{software-gmp, software-mpfr, software-mpfi, software-ntl}
    sage: with citations(record):
    ...       integrate(y^2, y, 0, 1)
    1/3
    sage: record
    [GMP, MPFR, MPFI, NTL, ginac, Maxima]

As an alternative way to access citations, we provide an eval like statement.

::

    sage: eval_citations("integrate(y^2, y, 0, 10)")
    [GMP, ginac, Maxima]

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
from sage.misc.cite import citation_enable, citation_disable
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
            GMP, ginac
        sage: latest_citations()
        [GMP, ginac]
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
            GMP, Singular
        sage: latest_citations()
        [GMP, Singular]

    ::

        sage: with citations():
        ...       d = dict()
        The computation did not use any registered components.
    """
    import warnings

    if record is None:
        _record = CitationRecord()
    else:
        _record = record

    citation_enable(_record)
    yield
    citation_disable(_record)

    if record is None:
        import sage
        sage.misc.citation._latest_citations = _record

        if len(_record) == 0:
            print "The computation did not use any registered components."
        else:
            print "The computation used the following components."
            print "Access them as a list by calling latest_citations()."
            print "    " + ", ".join(map(repr, _record))

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
        [GMP, ginac, Maxima]
    """
    import inspect
    from sage.misc.all import sage_eval

    if locals is None:
        locals = inspect.stack()[1][0].f_globals

    record = CitationRecord()
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
