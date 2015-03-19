from sage.combinat.partition import Partitions
from sage.combinat.species2.cycle_index_series import CIS


class CISPerm(CIS):
    def Frobenius_characteristic(self, n):
        """
        TEST::
        
            sage: ZPerm = CISPerm()
            sage: list(ZPerm.generating_series().coefficients(6))
            [1, 1, 2, 6, 24, 120, 720]
            sage: print list(ZPerm.type_generating_series().coefficients(19))
            [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77, 101, 135, 176, 231, 297, 385, 490]

        """
        p = self._sym_.p()
        return p.sum_of_monomials(Partitions(n))

    def _repr_(self):
        return "Z_P"