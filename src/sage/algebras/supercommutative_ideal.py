r"""
Ideals of supercommutative algebras

A supercommutative algebra is a super algebra satisfying

.. MATH::

	x*y = (-1)^{|x|\cdot |y|}y*x

EXAMPLES:

	The classic example of a supercommutative algebra is the exterior
	algebra, where `|x|=1` for all `x`. For example, we can construct
	the exterior algebra in four variables::

	sage: E.<x,y,z,w> = ExteriorAlgebra(QQ)
	sage: y*x
	-x*y

AUTHORS:
- Trevor K. Karn (2021-07-21): initial version

# ****************************************************************************
# Copyright (C) 2021 Trevor K. Karn <karnx018 at umn.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
# https://www.gnu.org/licenses/
# ****************************************************************************

"""

class Ideal_sc(Ideal_nc):
    r"""
    Ideals of supercommutative algebras
    """

    def __init__(self, ring, gens, coerce=True, side='twosided'):

        # check that ring is supercommutative
        if ring not in SupercommutativeAlgebras
            raise ValueError(f'{ring} must be in the category of SupercommutativeAlgebras')

        Ideal_nc.__init__(self,ring,gens,coerce,side)

    def groebner_basis(self, degbound = None):

        raise NotImplementedError

    def __contains__(self,x):

        raise NotImplementedError

    def reduce(self,x):

        raise NotImplementedError


#     # Q: how to  make this different from the following:

#sage: A.<x,y,z> = FreeAlgebra(QQ,3)
#sage: G = A.g_algebra({x*y:-y*x})
