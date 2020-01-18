# -*- coding: utf-8 -*-
"""
Jacobi motives

EXAMPLES::

    sage: from sage.modular.jacobi_motives import JacobiMotive
    sage: M = JacobiMotive((1/5,)*5); M
    Jacobi Motive for (1/5),(1/5),(1/5),(1/5),(1/5)

    sage: M = JacobiMotive((1/2,)*2); M
    Jacobi Motive for (1/2),(1/2)

    sage: M = JacobiMotive((2/3,2/3),(1/3,)); M
    Jacobi Motive for (2/3),(2/3) - (1/3)

REFERENCES:

- Mark Watkins, Jacobi sums and Grössencharacters, 2018, p. 111-122.

- http://pmb.univ-fcomte.fr/2018/PMB_Watkins.pdf

- https://magma.maths.usyd.edu.au/magma/handbook/text/1544

- https://magma.maths.usyd.edu.au/magma/handbook/text/1545

- https://magma.maths.usyd.edu.au/magma/handbook/text/1546

"""
# ****************************************************************************
#       Copyright (C) 2020     Frédéric Chapoton
#                              Kiran S. Kedlaya <kskedl@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import annotations

from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.arith.functions import lcm
from sage.arith.misc import euler_phi
from sage.rings.number_field.number_field import CyclotomicField


class JacobiMotive(object):
    def __init__(self, positive, negative=None):
        r"""
        INPUT:

        - ``positive`` and ``negative`` -- lists of rational numbers

        The ``negative`` argument is optional.

        The rational numbers in both arguments are considered modulo `\ZZ`.
        """
        if negative is None:
            negative = tuple()
        posi = (QQ(z) for z in positive)
        nega = (QQ(z) for z in negative)
        posi = [z - z.floor() for z in posi if z not in ZZ]
        nega = [z - z.floor() for z in nega if z not in ZZ]
        clean_nega = []
        for z in nega:
            if z in posi:
                posi.remove(z)
            else:
                clean_nega.append(z)
        self._posi = tuple(sorted(z for z in posi))
        self._nega = tuple(sorted(z for z in clean_nega))
        if (ZZ.sum(self._posi) - ZZ.sum(self._nega)).denominator() != 1:
            raise ValueError('sum of input is not an integer')
        self._m = lcm(z.denominator() for z in self._posi + self._nega)

    def __repr__(self) -> str:
        """
        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/5,)*5); M
            Jacobi Motive for (1/5),(1/5),(1/5),(1/5),(1/5)
        """
        text_posi = ",".join("({})".format(z) for z in self._posi)
        text_nega = ",".join(" - ({})".format(z) for z in self._nega)
        return "Jacobi Motive for " + text_posi + text_nega

    def scale(self, u) -> JacobiMotive:
        """
        Return another Jacobi motive by scaling ``self``.

        INPUT:

        - `u` -- integer, assumed to be invertible modulo `m`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/13,3/13,9/13))
            sage: M.scale(2)
            Jacobi Motive for (2/13),(5/13),(6/13)
            sage: M.scale(-1)
            Jacobi Motive for (4/13),(10/13),(12/13)
        """
        new_posi = tuple(u * z for z in self._posi)
        new_nega = tuple(u * z for z in self._nega)
        return JacobiMotive(new_posi, new_nega)

    def is_totally_real(self) -> bool:
        """
        Return whether ``self`` is totally real.

        This means invariant by scaling by `-1`.

        .. SEEALSO:: :meth:`is_cm`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/13,3/13,9/13))
            sage: M.is_totally_real()
            False
            sage: M = JacobiMotive((1/3,2/3))
            sage: M.is_totally_real()
            True
        """
        return self == self.scale(-1)

    def is_cm(self) -> bool:
        """
        Return whether ``self`` is cm (complex-multiplication).

        This means not invariant by scaling by `-1`.

        .. SEEALSO:: :meth:`is_totally_real`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/13,3/13,9/13))
            sage: M.is_totally_real()
            False
            sage: M = JacobiMotive((1/3,2/3))
            sage: M.is_totally_real()
            True
        """
        return not self.is_totally_real()

    def degree(self):
        """
        Return the degree.

        This is the sum of the Hodge numbers.

        TODO : CHECK : is this correct ?
        or should we use the degree of the field of definition ?

        .. SEEALSO::

            :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/13,3/13,9/13))
            sage: M.degree()
            12
        """
        return euler_phi(self._m)

    def weight(self):
        """
        Return the weight of ``self``.

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/5,)*5)
            sage: M.weight()
            5
            sage: M = JacobiMotive((1/3,2/3))
            sage: M.weight()
            2
            sage: M = JacobiMotive((2/3,2/3),(1/3,))
            sage: M.weight()
            1
        """
        resu = sum(1 for e in self._posi if e)
        resu -= sum(1 for e in self._nega if e)
        return ZZ(resu)

    def hodge_numbers(self):
        """
        Return the Hodge numbers.

        .. SEEALSO::

            :meth:`degree`

        TODO : how can one compute that ?
        """
        m = self._m
        for u in m.coprime_integers(m + 1):
            iter_posi = (u * z for z in self._posi)
            local_hodge_weight = ZZ.sum(y - y.floor() for y in iter_posi)
            iter_posi = (u * z for z in self._nega)
            local_hodge_weight -= ZZ.sum(y - y.floor() for y in iter_posi)
        return local_hodge_weight

    def effective_weight(self):
        """
        Return the effective weight of ``self``.

        This is the width of the Hodge structure.

        .. SEEALSO::

            :meth:`hodge_numbers`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/5,)*5)
            sage: M.effective_weight()  # not tested
            3
        """
        pass

    def __eq__(self, other) -> bool:
        """
        Return whether ``self`` is equal to ``other``.

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M1 = JacobiMotive((1/13,3/13,9/13))
            sage: M2 = JacobiMotive((1/13,3/13,9/13,1/2),(1/2,))
            sage: M3 = JacobiMotive((1/2,) * 2)
            sage: M1 == M2
            True
            sage: M1 == M3
            False
        """
        return self._posi == other._posi and self._nega == other._nega

    def __ne__(self, other) -> bool:
        """
        Return whether ``self`` is not equal to ``other``.

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M1 = JacobiMotive((1/13,3/13,9/13))
            sage: M2 = JacobiMotive((1/13,3/13,9/13,1/2),(1/2,))
            sage: M3 = JacobiMotive((1/2,) * 2)
            sage: M1 != M2
            False
            sage: M1 != M3
            True
        """
        return not (self == other)

    def __hash__(self) -> int:
        """
        Return the hash of ``self``.

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M1 = JacobiMotive((1/13,3/13,9/13))
            sage: M2 = JacobiMotive((1/13,3/13,9/13,1/2),(1/2,))
            sage: hash(M1) == hash(M2)
            True
        """
        return hash((self._posi, self._nega))

    def field_of_definition(self):
        """
        Return the field of definition of ``self``.

        See :meth:`invariant scalings`

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive([2/3,2/3],[1/3])
            sage: M.field_of_definition()
            Cyclotomic Field of order 3 and degree 2

            sage: E = JacobiMotive([1/7,2/7,4/7])
            sage: E.field_of_definition()
            Number Field in z0 with defining polynomial x^2 - x + 2 with ...
        """
        m = self._m
        L = self.invariant_scalings()
        K = CyclotomicField(m, 'z')
        z = K.gen()
        G = K.galois_group()
        power_list = z.powers(m)
        convert = ((sigma, power_list.index(sigma(z))) for sigma in G)
        G_gens = [g for g, i in convert if i in L]
        return G.subgroup(G_gens).fixed_field()[0]

    def invariant_scalings(self) -> list:
        """
        Return the list of scalings that leave the motive invariant.

        EXAMPLES::

            sage: from sage.modular.jacobi_motives import JacobiMotive
            sage: M = JacobiMotive((1/13,3/13,9/13))
            sage: M.invariant_scalings()
            [1, 3, 9]
            sage: M = JacobiMotive([1/7,2/7,4/7])
            sage: M.invariant_scalings()
            [1, 2, 4]
        """
        m = self._m
        return [u for u in m.coprime_integers(m) if self.scale(u) == self]

    def euler_factor(self, p):
        """
        Return the Euler factor at the prime `p`.

        TODO, using either complex or p-adic Gauss sums
        """
        pass

    def __mul__(self, other) -> JacobiMotive:
        """
        Return the tensor product of ``self`` and ``other``.
        """
        pass

    def __div__(self, other) -> JacobiMotive:
        """
        Return the tensor quotient of ``self`` and ``other``, if possible.
        """
        pass
