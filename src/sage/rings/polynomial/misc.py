# ****************************************************************************
#       Copyright (C) 2021, Vincent Delecroix <20100.delecroix@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

# TODO: move this at the level of category
def inverse_of_unit(elt):
    """
    Generic inversion of univariate or multivariate element of a graded algebra.

    TESTS::

        sage: Integers(1)['x'](0).inverse_of_unit()
        0
    """
    P = elt.parent()
    if not elt.is_unit():
        raise ArithmeticError(f"{elt} is not a unit in {elt.parent()}")

    u = elt.constant_coefficient()
    ui = u.inverse_of_unit()
    n = - ui * (elt - u)  # nilpotent
    nn = n
    res = P.one()
    while nn:
        res += nn
        nn *= n
    return ui * res
