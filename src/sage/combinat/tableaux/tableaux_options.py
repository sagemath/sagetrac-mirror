r"""
GlobalOptions for (skew) tableau.

Constrols display options for tableau, like using English
vs. French notation.

AUTHORS:

- Josh Swanson (2015): tableaux refactoring/cleanup
"""
#*****************************************************************************
#       Copyright (C) 2015 Josh Swanson
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.structure.global_options import GlobalOptions

# GlobalOptions input docstrings are not doctested, so force doctesting using
#   dummy methods.
def _doc():
    r"""
    Sets the global options for elements of the :class:`SkewTableau` class.
    By default, they are displayed as a list, `\LaTeX`ed as a Young
    diagram, and English convention is used.
    """
    pass

# TODO: fix failing doctests
def _end_doc():
    r"""
    .. NOTE::

        Changing the ``convention`` for tableaux also changes the
        ``convention`` for partitions.

    If no parameters are set, then the function returns a copy of the
    options dictionary.

    EXAMPLES::

        sage: T = Tableau([[1,2,3],[4,5]])
        sage: T
        [[1, 2, 3], [4, 5]]
        sage: Tableaux.global_options(display="array")
        sage: T
          1  2  3
          4  5
        sage: Tableaux.global_options(convention="french")
        sage: T
          4  5
          1  2  3

    Changing the ``convention`` for tableaux also changes the ``convention``
    for partitions and vice versa::

        sage: P = Partition([3,3,1])
        sage: print P.ferrers_diagram()
        *
        ***
        ***
        sage: Partitions.global_options(convention="english")
        sage: print P.ferrers_diagram()
        ***
        ***
        *
        sage: T
          1  2  3
          4  5

    The ASCII art can also be changed::

        sage: t = Tableau([[1,2,3],[4,5]])
        sage: ascii_art(t)
          1  2  3
          4  5
        sage: Tableaux.global_options(ascii_art="table")
        sage: ascii_art(t)
        +---+---+
        | 4 | 5 |
        +---+---+---+
        | 1 | 2 | 3 |
        +---+---+---+
        sage: Tableaux.global_options(ascii_art="compact")
        sage: ascii_art(t)
        |4|5|
        |1|2|3|
        sage: Tableaux.global_options.reset()
    """
    pass

TableauOptions=GlobalOptions(name='skew tableaux',
    doc        = _doc.__doc__,
    end_doc    = _end_doc.__doc__,
    display    = dict(default        = "list",
                      description    = 'Controls the way in which tableaux are printed',
                      values         = dict(list    = 'print tableaux as lists',
                                            diagram = 'display as Young diagram (similar to :meth:`~sage.combinat.tableau.Tableau.pp()`',
                                            compact = 'minimal length string representation'),
                      alias          = dict(array           = "diagram",
                                            ferrers_diagram = "diagram",
                                            young_diagram   = "diagram"),
                      case_sensitive = False),
    ascii_art  = dict(default        = "repr",
                      description    = 'Controls the ascii art output for tableaux',
                      values         = dict(repr    = 'display using the diagram string representation',
                                            table   = 'display as a table',
                                            compact = 'minimal length ascii art'),
                      case_sensitive = False),
    latex      = dict(default        = "diagram",
                      description    = 'Controls the way in which tableaux are latexed',
                      values         = dict(list    = 'as a list',
                                            diagram = 'as a Young diagram'),
                      alias          = dict(array           = "diagram",
                                            ferrers_diagram = "diagram",
                                            young_diagram   = "diagram"),
                      case_sensitive = False),
    convention = dict(default        = "English",
                      description    = 'Sets the convention used for displaying tableaux and partitions',
                      values         = dict(English = 'use the English convention',
                                            French  = 'use the French convention'),
                      case_sensitive = False),
    notation   = dict(alt_name       = "convention")
)
