"""
Hyperelliptic curves of genus 3 over the rationals
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_g3_generic, hyperelliptic_rational_field

class HyperellipticCurve_g3_rational_field(
    hyperelliptic_g3_generic.HyperellipticCurve_g3_generic,
    hyperelliptic_rational_field.HyperellipticCurve_rational_field):
    pass
