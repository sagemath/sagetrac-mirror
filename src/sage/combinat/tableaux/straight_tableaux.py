r"""
Tableax Parent class and subclasses.

See Tableau and its subclasses for the corresponding Element classes.
A better name for might be ``StraightTableaux`` instead of ``Tableaux``,
though the latter is used for historical reasons.

AUTHORS:

- Mike Hansen (2007): initial version
- Jason Bandlow (2011): updated to use Parent/Element model, and many
  minor fixes
- Andrew Mathas (2012-13): completed the transition to the parent/element model
  begun by Jason Bandlow
- Travis Scrimshaw (11-22-2012): Added tuple options, changed ``*katabolism*``
  to ``*catabolism*``. Cleaned up documentation.
- Josh Swanson (and others) (2015): tableaux refactoring/cleanup

CLASSES:

.. TODO:: List classes.
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>
#                     2011 Jason Bandlow <jbandlow@gmail.com>
#                     2015 Josh Swanson
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

import six

from sage.combinat.tableaux.skew_tableaux import SkewTableaux
from sage.combinat.tableaux.straight_tableau import Tableau

class Tableaux(SkewTableaux):
    Element = Tableau

    def _repr_(self):
        r"""
        Return the representation string.

        OUTPUT:

        A string.
        """
        return "Straight Tableaux"

