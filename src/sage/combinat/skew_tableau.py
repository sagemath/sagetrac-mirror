r"""
Skew Tableaux

AUTHORS:

- Mike Hansen: Initial version
- Travis Scrimshaw, Arthur Lubovsky (2013-02-11):
  Factored out ``CombinatorialClass``
"""
#*****************************************************************************
#       Copyright (C) 2007 Mike Hansen <mhansen@gmail.com>,
#       Copyright (C) 2013 Travis Scrimshaw <tscrim at ucdavis.edu>
#       Copyright (C) 2013 Arthur Lubovsky
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

from sage.misc.lazy_import import lazy_import
lazy_import('sage.combinat.tableaux.skew_tableau',
            'SkewTableauFactory',
            'SkewTableau',
             deprecation=(18013, 'Importing SkewTableau from here is deprecated. If you need to use it, please import it directly at sage.combinat.tableaux.skew_tableau.SkewTableauFactory'))
lazy_import('sage.combinat.tableaux.skew_tableaux',
            'SkewTableaux',
             deprecation=18013)
lazy_import('sage.combinat.tableaux.skew_tableaux',
            'SemistandardSkewTableauxFactory',
            'SemistandardSkewTableaux',
             deprecation=(18013, 'Importing SemistandardSkewTableaux from here is deprecated. If you need to use it, please import it directly at sage.combinat.tableaux.skew_tableaux.SemistandardSkewTableauxFactory'))
lazy_import('sage.combinat.tableaux.skew_tableaux',
            'StandardSkewTableauxFactory',
            'StandardSkewTableaux',
             deprecation=(18013, 'Importing StandardSkewTableaux from here is deprecated. If you need to use it, please import it directly at sage.combinat.tableaux.skew_tableaux.StandardSkewTableauxFactory'))
