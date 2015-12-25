r"""
Kernel Subgroups

The kernel subgroup of a homomorphism.

AUTHORS:

- Travis Scrimshaw (12-2015): initial version
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrimsh at umn.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.cachefunc import cached_method
from sage.categories.groups import Groups
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

class KernelSubgroup(Parent, UniqueRepresentation):
    r"""
    The kernel (normal) subgroup.

    Let `\phi : G \to H` be a group homomorphism. Then the kernel
    `K = \{\phi(g) = 1 | g \in G\}` is a normal subgroup of `G`.
    """
    def __init__(self, morphism):
        """
        Initialize ``self``.
        """
        self._morphism = morphism
        cat = Groups().Subobjects()
        if morphism.domain() in Groups().Finite():
            cat = cat.Finite()
        Parent.__init__(self, category=cat)

    def _repr_(self):
        """
        Return a string representation of ``self``.
        """
        return "Kernel subgroup defined by {}".format(self._morphism)

    def morphism(self):
        """
        Return the defining morphism of ``self``.
        """
        return self._morphism

    @cached_method
    def ambient(self):
        """
        Return the ambient group of ``self``.
        """
        return self._morphism.domain()

    def _an_element_(self):
        """
        Return an element of ``self``.
        """
        return self.element_class(self, self.ambient().one())

    def lift(self, x):
        """
        Lift ``x`` to the ambient group of ``self``.
        """
        return x.value

    def retract(self, x):
        """
        Convert ``x`` to an element of ``self``.
        """
        return self._element_constructor_(x)

    def _element_constructor_(self, x):
        """
        Construct an element of ``self`` from ``x``.
        """
        if self._morphism(x) != self._morphism.codomain().one():
            raise ValueError("{} is not in the kernel of {}".format(x, self._morphism))
        return self.element_class(self, x)

    def __iter__(self):
        """
        Iterate through ``self``.
        """
        for g in self.ambient():
            try:
                yield self(g)
            except ValueError:
                pass

    class Element(ElementWrapper):
        def _mul_(self, other):
            """
            Multiply ``self`` and ``other``.
            """
            return type(self)(self.parent(), self.value * other.value)

        def __invert__(self):
            """
            Return the inverse of ``self``.
            """
            return type(self)(self.parent(), ~self.value)

