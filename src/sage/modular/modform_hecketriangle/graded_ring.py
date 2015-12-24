r"""
Graded rings of modular forms for Hecke triangle groups

AUTHORS:

- Jonas Jermann (2013): initial version

"""

#*****************************************************************************
#       Copyright (C) 2013-2014 Jonas Jermann <jjermann2@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.all import ZZ, QQ, infinity

from sage.rings.ring import CommutativeAlgebra
from sage.categories.all import CommutativeAlgebras
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.cachefunc import cached_method

from hecke_triangle_groups import HeckeTriangleGroup
from abstract_ring import FormsRing_abstract
from analytic_type import AT


def canonical_parameters(group, base_ring, red_hom, n=None, frac=False):
    r"""
    Return a canonical version of the parameters.

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import canonical_parameters
        sage: canonical_parameters(4, ZZ, 1)
        (Hecke triangle group for n = 4, Integer Ring, True, 4, False)
        sage: canonical_parameters(infinity, RR, 0)
        (Hecke triangle group for n = +Infinity,
         Real Field with 53 bits of precision,
         False,
         +Infinity,
         False)
    """

    if not (n is None):
        group = n

    if (group == infinity):
        group = HeckeTriangleGroup(infinity)
    else:
        try:
            group = HeckeTriangleGroup(ZZ(group))
        except TypeError:
            group = HeckeTriangleGroup(group.n())

    red_hom = bool(red_hom)
    n = group.n()

    if (n == infinity):
        frac = bool(frac)
    else:
        frac = False

    return (group, base_ring, red_hom, n, frac)


class QuasiMeromorphicModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi meromorphic modular forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiMeromorphicModularFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(4, ZZ, 1)
            sage: QuasiMeromorphicModularFormsRing(4, ZZ, 1) == QuasiMeromorphicModularFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) quasi meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiMeromorphicModularFormsRing
            sage: MR = QuasiMeromorphicModularFormsRing(4, ZZ, 1)
            sage: MR
            QuasiMeromorphicModularFormsRing(n=4) over Integer Ring
            sage: MR.analytic_type()
            quasi meromorphic modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring

            sage: QuasiMeromorphicModularFormsRing(n=infinity)
            QuasiMeromorphicModularFormsRing(n=+Infinity) over Integer Ring
        """

        analytic_type = AT(["quasi", "mero"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class QuasiWeakModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi weakly holomorphic modular forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiWeakModularFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(5, CC, 0)
            sage: QuasiWeakModularFormsRing(5, CC, 0) == QuasiWeakModularFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) quasi weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiWeakModularFormsRing
            sage: MR = QuasiWeakModularFormsRing(5, CC, 0)
            sage: MR
            QuasiWeakModularFormsRing(n=5) over Complex Field with 53 bits of precision
            sage: MR.analytic_type()
            quasi weakly holomorphic modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision
        """

        analytic_type = AT(["quasi", "weak"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class QuasiModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiModularFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(6, ZZ, True)
            sage: QuasiModularFormsRing(6, ZZ, True) == QuasiModularFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) quasi modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiModularFormsRing
            sage: MR = QuasiModularFormsRing(6, ZZ, True)
            sage: MR
            QuasiModularFormsRing(n=6) over Integer Ring
            sage: MR.analytic_type()
            quasi modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        analytic_type = AT(["quasi", "holo"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class QuasiCuspFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) quasi cusp forms
    for the given group and base ring.
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, QuasiCuspFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(7, ZZ, 1)
            sage: QuasiCuspFormsRing(7, ZZ, 1) == QuasiCuspFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) quasi cusp forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) quasi cusp forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import QuasiCuspFormsRing
            sage: MR = QuasiCuspFormsRing(7, ZZ, 1)
            sage: MR
            QuasiCuspFormsRing(n=7) over Integer Ring
            sage: MR.analytic_type()
            quasi cuspidal
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        analytic_type = AT(["quasi", "cusp"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class MeromorphicModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) meromorphic modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, MeromorphicModularFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(4, ZZ, 1)
            sage: MeromorphicModularFormsRing(4, ZZ, 1) == MeromorphicModularFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) meromorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import MeromorphicModularFormsRing
            sage: MR = MeromorphicModularFormsRing(4, ZZ, 1)
            sage: MR
            MeromorphicModularFormsRing(n=4) over Integer Ring
            sage: MR.analytic_type()
            meromorphic modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        analytic_type = AT(["mero"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class WeakModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) weakly holomorphic modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, WeakModularFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(5, ZZ, 0)
            sage: WeakModularFormsRing(5, ZZ, 0) == WeakModularFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) weakly holomorphic modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import WeakModularFormsRing
            sage: MR = WeakModularFormsRing(5, ZZ, 0)
            sage: MR
            WeakModularFormsRing(n=5) over Integer Ring
            sage: MR.analytic_type()
            weakly holomorphic modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        analytic_type = AT(["weak"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class ModularFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) modular forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: ModularFormsRing(3, ZZ, 0) == ModularFormsRing()
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) modular forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            (for instance ``theta``) in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) modular forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import ModularFormsRing
            sage: MR = ModularFormsRing()
            sage: MR
            ModularFormsRing(n=3) over Integer Ring
            sage: MR.analytic_type()
            modular
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Integer Ring
        """

        analytic_type = AT(["holo"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()

class CuspFormsRing(FormsRing_abstract, CommutativeAlgebra, UniqueRepresentation):
    r"""
    Graded ring of (Hecke) cusp forms
    for the given group and base ring
    """

    @staticmethod
    def __classcall__(cls, group = HeckeTriangleGroup(3), base_ring = ZZ, red_hom = False, n=None, frac=False):
        r"""
        Return a (cached) instance with canonical parameters.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import (canonical_parameters, CuspFormsRing)
            sage: (group, base_ring, red_hom, n, frac) = canonical_parameters(5, CC, True)
            sage: CuspFormsRing(5, CC, True) == CuspFormsRing(group, base_ring, red_hom, n, frac)
            True
        """

        (group, base_ring, red_hom, n, frac) = canonical_parameters(group, base_ring, red_hom, n, frac)
        return super(FormsRing_abstract,cls).__classcall__(cls, group=group, base_ring=base_ring, red_hom=red_hom, n=n, frac=frac)

    def __init__(self, group, base_ring, red_hom, n, frac):
        r"""
        Return the graded ring of (Hecke) cusp forms
        for the given ``group`` and ``base_ring``.

        INPUT:

        - ``group``      -- The Hecke triangle group (default: ``HeckeTriangleGroup(3)``)

        - ``base_ring``  -- The base_ring (default: ``ZZ``).

        - ``red_hom``    -- If True then results of binary operations are considered
                            homogeneous whenever it makes sense (default: False).
                            This is mainly used by the spaces of homogeneous elements.

        - ``frac``       -- If True then include forms with fractional orders at the cusp ``-1``
                            in case of ``n=infinity`` (default: ``True``).
                            If ``n!=infinity`` this parameter is ignored.

        OUTPUT:

        The corresponding graded ring of (Hecke) cusp forms
        for the given ``group`` and ``base_ring``.

        EXAMPLES::

            sage: from sage.modular.modform_hecketriangle.graded_ring import CuspFormsRing
            sage: MR = CuspFormsRing(5, CC, True)
            sage: MR
            CuspFormsRing(n=5) over Complex Field with 53 bits of precision
            sage: MR.analytic_type()
            cuspidal
            sage: MR.category()
            Category of commutative algebras over Fraction Field of Univariate Polynomial Ring in d over Complex Field with 53 bits of precision

            sage: CuspFormsRing(n=infinity, base_ring=CC, red_hom=True)
            CuspFormsRing(n=+Infinity) over Complex Field with 53 bits of precision
        """

        analytic_type = AT(["cusp"])
        if (n==infinity and frac):
            analytic_type = analytic_type.extend_by("frac")

        FormsRing_abstract.__init__(self, group=group, base_ring=base_ring, red_hom=red_hom, n=n, analytic_type=analytic_type)
        CommutativeAlgebra.__init__(self, base_ring=self.coeff_ring(), category=CommutativeAlgebras(self.coeff_ring()))
        self._post_init()


def ThetaQuasiMeromorphicModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of quasi meromorphic modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaQuasiMeromorphicModularFormsRing
        sage: MR = ThetaQuasiMeromorphicModularFormsRing()
        sage: MR
        ThetaQuasiMeromorphicModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta quasi meromorphic modular
    """

    return QuasiMeromorphicModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaQuasiWeakModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of quasi weak modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaQuasiWeakModularFormsRing
        sage: MR = ThetaQuasiWeakModularFormsRing()
        sage: MR
        ThetaQuasiWeakModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta quasi weakly holomorphic modular
    """

    return QuasiWeakModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaQuasiModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of quasi modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaQuasiModularFormsRing
        sage: MR = ThetaQuasiModularFormsRing()
        sage: MR
        ThetaQuasiModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta quasi modular
    """

    return QuasiModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaQuasiCuspFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of quasi cusp forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaQuasiCuspFormsRing
        sage: MR = ThetaQuasiCuspFormsRing()
        sage: MR
        ThetaQuasiCuspFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta quasi cuspidal
    """

    return QuasiCuspFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaMeromorphicModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of meromorphic modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaMeromorphicModularFormsRing
        sage: MR = ThetaMeromorphicModularFormsRing()
        sage: MR
        ThetaMeromorphicModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta meromorphic modular
    """

    return MeromorphicModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaWeakModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of weakly holomorphic modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaWeakModularFormsRing
        sage: MR = ThetaWeakModularFormsRing()
        sage: MR
        ThetaWeakModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta weakly holomorphic modular
    """

    return WeakModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaModularFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of modular forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaModularFormsRing
        sage: MR = ThetaModularFormsRing()
        sage: MR
        ThetaModularFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta modular
    """

    return ModularFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)

def ThetaCuspFormsRing(base_ring = ZZ, red_hom = False):
    r"""
    Graded ring of cusp forms for the Theta group and the given base ring (with fractional orders at the cusp ``-1``)

    EXAMPLES::

        sage: from sage.modular.modform_hecketriangle.graded_ring import ThetaCuspFormsRing
        sage: MR = ThetaCuspFormsRing()
        sage: MR
        ThetaCuspFormsRing() over Integer Ring
        sage: MR.analytic_type()
        theta cuspidal
    """

    return CuspFormsRing(group = HeckeTriangleGroup(infinity), base_ring = base_ring, red_hom = red_hom, n=None, frac=True)
