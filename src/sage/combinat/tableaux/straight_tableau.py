r"""
Tableau Element class and subclasses.

See Talbeaux and its subclasses for the corresponding Parent classes.
A better name for might be ``StraightTableau`` instead of ``Tableau``,
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

from sage.combinat.tableaux.skew_tableau import SkewTableau

class StraightTableau(SkewTableau):
    _generic_parent = parent_class('StraightTableaux')

    def __init__(self, parent, l=[], dictionary=[], check=False):
        SkewTableau.__init__(self, parent, l, dictionary, check)

class ExampleTableau(BadShapeTableau):
    def __init__(self):
        BadShapeTableau.__init__(self, d = {(0,1): 3, (-3, 5): 17, (1,1): "a"})
