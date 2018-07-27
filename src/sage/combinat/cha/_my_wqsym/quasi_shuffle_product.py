# -*- coding: utf-8 -*-
r"""
Shuffle product of ordered partitions sets

AUTHOR:

- Jean-Baptiste Priez
"""
#*****************************************************************************
#       Copyright (C) 2012 Jean-Baptiste Priez <jbp@kerios.fr>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.combinat import CombinatorialClass
from sage.combinat.integer_vector import IntegerVectors
from sage.sets.set import Set


class QuasiShuffleProduct(CombinatorialClass):

    def __init__(self, p1, p2):
        self._p1 = p1
        self._p2 = p2

    def __repr__(self):
        return "Quasi shuffle product of %s and %s" % (
            repr(self._p1),
            repr(self._p2)
        )

    def __contains__(self, x):
        # TODO:: test the instance
        # FIXME:: implement
        return False

    def cardinality(self):
        from sage.rings.arith import binomial
        # C(p,q) = sum_k b(q, k) x b(p+q-k, q)
        p = len(self._p1)
        q = len(self._p2)
        acc = 0
        # sum k <- 0 to p
        for k in range(p + 1):
            acc += binomial(q, k) * binomial(p + q - k, q)
        return acc

    def __iter__(self):
        if len(self._p1) == 0:
            yield self._p2
        elif len(self._p2) == 0:
            yield self._p1
        else:
            # {a}S::{b}S' = {a} (S::{b}S') + {a,b} (S::S') + {b} ({a}S::S')
            # # {a} (S::{b}S')
            a = [self._p1[0]]
            for osp in QuasiShuffleProduct(self._p1[1:], self._p2):
                yield a + osp
            # # {a,b} (S::S')
            b = [self._p1[0] + self._p2[0]]
            for osp in QuasiShuffleProduct(self._p1[1:], self._p2[1:]):
                yield b + osp
            # # {b} ({a}S::S')
            c = [self._p2[0]]
            for osp in QuasiShuffleProduct(self._p1, self._p2[1:]):
                yield c + osp


class Shifted_QuasiShuffleProduct(QuasiShuffleProduct):

    def __init__(self, p1, p2):
        self.__p2 = p2
        maxi = sum(len(x) for x in p1)
                   # max(reduce(lambda x, y: x.add(y), p1,frozenset([0])))
        p2s = [Set([i + maxi for i in s]) for s in p2]
        QuasiShuffleProduct.__init__(self, p1, p2s)

    def __repr__(self):
        return "Shifted Quasi shuffle product of %s and %s" % (
            repr(self._p1),
            repr(self.__p2)
        )
