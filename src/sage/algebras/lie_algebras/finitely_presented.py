"""
Finitely Presented Lie Algebras

AUTHORS:

- Travis Scrimshaw (2013-05-03): Initial version
"""

#*****************************************************************************
#  Copyright (C) 2013 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty
#    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
#  See the GNU General Public License for more details; the full text
#  is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.algebras.lie_algebras.lie_algebra import FinitelyGeneratedLieAlgebra
from sage.algebras.lie_algebras.ideal import LieAlgebraIdeal
from sage.algebras.lie_algebras.quotient import QuotientLieAlgebra, QuotientLieAlgebraElement

class FinitelyPresentedLieAlgebra(QuotientLieAlgebra):
    r"""
    A finitely presented Lie algebra.

    INPUT:

    - ``F`` -- the defining free Lie algebra
    - ``relations`` -- the defining relations (as elements in ``F``)
    - ``names`` -- (optional) the variables names; if not provided then
      default to the variable names of ``F``
    - ``category`` -- (optional) the category

    REFERENCES:

    .. [GK] *An algorithm for analysis of the structure of finitely presented
       Lie algebras*. Vladimir Gerdt and Vladimir Kornyak. Discrete
       Mathematics and Theoretical Computer Science. **1** (1997). 217-228.
       http://www.emis.de/journals/DMTCS/volumes/abstracts/pdfpapers/dm010113.pdf
    """
    def __init__(self, F, relations, names=None, category=None):
        """
        Initialize ``self``.
        """
        self._rels = relations
        I = LieAlgebraIdeal(F, relations, False)
        if names is None:
            names = F.variable_names()
        QuotientLieAlgebra.__init__(self, F, I, names, category)

    def relations(self):
        """
        Return the relations of ``self``.
        """
        return self._rels

    def _construct_UEA(self):
        """
        Construct the universal enveloping algebra of ``self``.
        """
        UEA = self._free.UEA()
        I = UEA.ideal([UEA(x) for x in self._rels])
        return UEA.quotient(I)

    class Element(QuotientLieAlgebraElement):
        def _repr_(self):
            """
            Return a string representation of ``self``.
            """
            return repr(self.value)

        def _reduce_(self):
            """
            Reduce the element modulo the defining relations of this
            Lie algebra.

            (Note that this has nothing to do with pickling.)

            TESTS::
            """
            I = self.parent().defining_ideal()
            self.value = I.reduce(self.value)

        def lift(self):
            """
            Return the lift of ``self`` to the universal enveloping algbra.
            """
            UEA = self.parent().universal_enveloping_algebra()
            return UEA(self.value.lift())

