# -*- coding: utf-8 -*-
r"""

Some characteristic formal power series

"""
# *****************************************************************************
#  Copyright (C) 2015 Jean-Baptiste Priez <jbp at kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
# ******************************************************************************
from sage.combinat.species2.formal_power_series import FPS
from sage.rings.infinity import Infinity
from sage.rings.integer import Integer


class ZeroFPS(FPS):

    def coefficient(self, n):
        return Integer(0)

    def _repr_(self):
        return "0"

    def _valuation_(self):
        return Infinity


class OneFPS(FPS):

    def coefficient(self, n):
        return Integer(1) if n == 0 else Integer(0)

    def _repr_(self):
        return "1"

    def _valuation_(self):
        return 0


class SingletonsFPS(FPS):

    def coefficient(self, n):
        return Integer(1) if n == 1 else Integer(0)

    def _repr_(self):
        return "x"

    def _valuation_(self):
        return 1


class SetsFPS(FPS):

    def coefficient(self, n):
        return Integer(1)

    def _repr_(self):
        return "1/(1-x)"

    def _valuation_(self):
        return 0