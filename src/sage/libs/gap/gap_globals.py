"""Common globals defined by GAP."""

###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################

# TODO: As mentioned in sage.libs.gap.gap_functions, gappy should provide an
# interface for extending the list of common GAP global variables.
from sage.misc.superseded import deprecation
from gappy.gap_globals import *
deprecation(31297,
    'this module has been subsumed by the gappy package and will be removed')
