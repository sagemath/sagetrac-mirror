r"""
Set of homomorphisms
"""

#*****************************************************************************
#       Copyright (C) 2014 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.schemes.generic.homset import SchemeHomset_points


class SchemeHomset_points_product_projective_spaces_ring(SchemeHomset_points):
    """
    Set of points of a product of projective spaces.

    INPUT:

    See :class:`~sage.schemes.generic.homset.SchemeHomset_generic`.

    EXAMPLES::

        sage: from sage.schemes.product_projective.homset import \
        ....:     SchemeHomset_points_product_projective_spaces_ring
        sage: P1xP1 = ProductProjectiveSpaces([1, 1], QQ)
        sage: SchemeHomset_points_product_projective_spaces_ring(Spec(QQ), P1xP1)
        Set of rational points of Product of projective
        spaces P^1 x P^1 over Rational Field
    """

