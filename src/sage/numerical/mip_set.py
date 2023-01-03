"""
Sets of feasible solutions of mixed integer linear programs
"""

from sage.numerical.mip import MixedIntegerLinearProgram
from sage.structure.parent import Parent, Set_generic
from sage.sets.set import Set_base, Set_boolean_operators, Set_add_sub_operators


class Set_mip(Set_generic, Set_base, Set_boolean_operators, Set_add_sub_operators):
    r"""
    Set of (projections of) feasible solutions of a mixed integer linear program

    INPUT:

    - ``mip`` - a :class:`~sage.numerical.mip.MixedIntegerLinearProgram`

    - ``variables`` - an instance of :class:`~sage.numerical.mip.MIPVariable` (or one of its
      elements), or lists of them; same as
      :meth:`~sage.numerical.mip.MixedIntegerLinearProgram.get_values`.
    """

    def __init__(self, mip: MixedIntegerLinearProgram,
                 variables, convert=None, category=None):
        self._mip = mip
        self._variables = variables
        self._convert = convert

        # TODO: Refine category to Sets().Enumerated() etc. based on variable types and bounds

        Parent.__init__(self, category=category)

    def mip(self):
        return self._mip

    def __iter__(self):
        r"""
        Enumerate solutions in the order of the objective function (starting with an
        optimal solution), refined to some well-ordering (lex when bounded)
        """
        raise NotImplementedError

    def __contains__(self, x):
        """
        Return ``True`` if ``x`` is in ``self``.
        """
        # Add constraints for x and check feasibility of the projected out variables
        raise NotImplementedError
