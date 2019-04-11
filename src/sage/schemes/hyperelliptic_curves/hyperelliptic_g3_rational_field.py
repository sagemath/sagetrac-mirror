"""
Hyperelliptic curves of genus 3 over the rationals
"""
from __future__ import absolute_import

#*****************************************************************************
#  Copyright (C) 2017 Sorina Ionica <sorina.ionica@gmail.com>
#                2017 Elisa Lorenzo Garcia <elisa.lorenzo@gmail.com>
#                2017 Anna Somoza <anna.somoza@upc.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_g3_generic, hyperelliptic_rational_field

class HyperellipticCurve_g3_rational_field(
    hyperelliptic_g3_generic.HyperellipticCurve_g3_generic,
    hyperelliptic_rational_field.HyperellipticCurve_rational_field):
    pass
