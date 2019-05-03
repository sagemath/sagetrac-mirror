"""
LaTeX printing support for sagetex

In order to support latex formatting specifically for sagetex, an object 
should define a special method ``_sagetex_(self)`` that returns a string, 
which will be typeset in a mathematical mode (the exact mode depends on 
circumstances). If not set, this should default to ``_latex_(self)``.


AUTHORS:

- Aram Dermenjian: original implementation
"""

#*****************************************************************************
#       Copyright (C) 2019 Aram Dermenjian <aram.dermenjian@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#
#
#  Code is based off that found in ./latex.py from William Stein and Joel Mohler
#*****************************************************************************



from sage.misc.latex import LatexExpr
from sage.misc.latex import LatexCall
from sage.misc.latex import Latex



def has_sagetex_attr(x):
    """
    Return ``True`` if ``x`` has a ``_sagetex_`` attribute, except if ``x``
    is a ``type``, in which case return ``False``.

    INPUT:

    - ``x`` -- a Sage object

    EXAMPLES::

        sage: from sage.misc.sagetex import has_sagetex_attr
        sage: has_sagetex_attr(identity_matrix(3))
        False
        sage: has_sagetex_attr("abc")  # strings have no _latex_ method
        False

    Types inherit the ``_sagetex_`` method of the class to which they refer,
    but calling it is broken::

        sage: T = type(identity_matrix(3)); T
        <type 'sage.matrix.matrix_integer_dense.Matrix_integer_dense'>
        sage: hasattr(T, '_sagetex_')
        False
        sage: T._sagetex_()
        Traceback (most recent call last):
        ...
        AttributeError: type object 'sage.matrix.matrix_integer_dense.Matrix_integer_de' has no attribute '_sagetex_'
        sage: has_sagetex_attr(T)
        False
    """
    return hasattr(x, '_sagetex_') and not isinstance(x, type)


class SagetexCall:
    r"""
    Typeset Sage objects via a ``__call__`` method to this class,
    typically by calling those objects' ``_sagetex_`` methods.  The
    class :class:`Sagetex` inherits from this. If a ``_sagetex_`` 
    method is not defined, then it attempts to call ``_latex_``
    using :class:`LatexCall`.

    EXAMPLES::

        sage: from sage.misc.sagetex import SagetexCall
        sage: SagetexCall()(ZZ)
        \Bold{Z}
        sage: SagetexCall().__call__(ZZ)
        \Bold{Z}

    This returns an instance of the class :class:`LatexExpr`::

        sage: type(SagetexCall()(ZZ))
        <class 'sage.misc.latex.LatexExpr'>
    """
    def __call__(self, x, combine_all=False):
        r"""
        Return a :class:`LatexExpr` built out of the argument ``x``.

        INPUT:

        - ``x`` -- a Sage object

        - ``combine_all`` -- boolean (Default: ``False``) If ``combine_all``
          is ``True`` and the input is a tuple, then it does not return a
          tuple and instead returns a string with all the elements separated by
          a single space.

        OUTPUT:

        A :class:`LatexExpr` built from ``x``

        EXAMPLES::

            sage: from sage.misc.sagetex import sagetex
            sage: sagetex(Integer(3))  # indirect doctest
            3
            sage: sagetex(1==0)
            \mathrm{False}
            sage: print(sagetex([x,2]))
            \left[x, 2\right]

        """
        if has_sagetex_attr(x):
            return LatexExpr(x._sagetex_())

        return LatexCall().__call__(x,combine_all)


# SagetexCall before Latex in order for __call__ to be forced from SagetexCall and not Latex
class Sagetex(SagetexCall, Latex):
    r"""nodetex
    Enter, e.g.,

    Use ``sagetex(...)`` to typeset a Sage object.

    .. WARNING::

       You must have dvipng (or dvips and convert) installed
       on your operating system, or this command won't work.

    EXAMPLES::

        sage: from sage.misc.sagetex import sagetex
        sage: sagetex(x^20 + 1)
        x^{20} + 1
        sage: sagetex(FiniteField(25,'a'))
        \Bold{F}_{5^{2}}
        sage: sagetex("hello")
        \text{\texttt{hello}}
        sage: sagetex(Tableau([[1,2]]))
        {\def\lr##1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$##1$}}}
        \raisebox{-.6ex}{$\begin{array}[b]{*{2}c}\cline{1-2}
        \lr{1}&\lr{2}\\\cline{1-2}
        \end{array}$}
        }
    """
    pass

# Note: This assignment is meant to replicate the Latex implementation

sagetex = Sagetex()

# Ensure that sagetex appear in the sphinx doc as a function
# so that the link :func:`sagetex` is correctly set up.
sagetex.__doc__  = Sagetex.__call__.__doc__
#########################################
