"""
Hyperelliptic curves of genus 3 over a finite field
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_finite_field, hyperelliptic_g3_generic

class HyperellipticCurve_g3_finite_field(
    hyperelliptic_g3_generic.HyperellipticCurve_g3_generic,
    hyperelliptic_finite_field.HyperellipticCurve_finite_field):
    pass

