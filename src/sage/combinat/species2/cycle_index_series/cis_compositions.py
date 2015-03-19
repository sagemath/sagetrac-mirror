# Série indicatrice de cycles des compositions
from sage.combinat.partition import Partitions
from sage.combinat.species2.cycle_index_series import CIS
from sage.rings.arith import multinomial


class CISCompo(CIS):

    def Frobenius_characteristic(self, n):
        """
        TEST::

            sage: ZComp = CISCompo()
            sage: list(ZComp.generating_series().coefficients(6))
            [1, 1, 3, 13, 75, 541, 4683]
            sage: list(ZComp.type_generating_series().coefficients(6))
            [1, 1, 2, 4, 8, 16, 32]

        """
        h = self._sym_.h()
        return sum(multinomial(lambd.to_exp()) * h(lambd)
                   for lambd in Partitions(n))

    def _repr_(self):
        return "Z_C"
