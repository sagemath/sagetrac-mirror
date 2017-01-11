"""
Hyperelliptic curves of genus 3 over a general ring
"""
from __future__ import absolute_import
#*****************************************************************************
#  Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from . import hyperelliptic_generic
#from . import jacobian_g3(?)
from . import invariants
#from .g3_reconstruction import fromshiodas as reconstruct


class HyperellipticCurve_g3_generic(hyperelliptic_generic.HyperellipticCurve_generic):
    def shioda_invariants(self):
        r"""
        Return the Shioda invariants `(J2, J3, J4, J5, J6, J7, J8, J9)` of Shioda, [Sh]_.

        .. SEEALSO::

            :meth:`sage.schemes.hyperelliptic_curves.invariants`

        EXAMPLES::

            sage: R.<x> = QQ[]
            sage: f = x^7 - x^4 + 3

        TESTS::

        """
        f, h = self.hyperelliptic_polynomials()
	assert h == 0, 'Argument must be a simplified model of genus 3.'
        return invariants.shioda_invariants(f)