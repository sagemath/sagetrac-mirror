"""
Hyperelliptic curves of genus 3 over a p-adic field
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_g3_generic, hyperelliptic_padic_field

class HyperellipticCurve_g3_padic_field(
    hyperelliptic_g3_generic.HyperellipticCurve_g3_generic,
    hyperelliptic_padic_field.HyperellipticCurve_padic_field):
    pass
