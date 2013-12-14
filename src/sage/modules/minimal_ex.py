r"""
Minimal Example for Nicolas
"""

from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation
from sage.misc.bindable_class import BindableClass
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.realizations import Realizations, Category_realization_of_parent
from sage.combinat.free_module import CombinatorialFreeModule

class BaseObject(Parent, UniqueRepresentation):
    """
    TESTS::

        sage: from sage.modules.minimal_ex import BaseObject
        sage: B = BaseObject(QQ)
        sage: TestSuite(B).run() # This breaks
        sage: R = B.a_realization()
        sage: TestSuite(B).run() # Now it does not
        sage: TestSuite(R).run()
    """
    def __init__(self, R):
        Parent.__init__(self, base=R, category=ModulesWithBasis(R).WithRealizations())

    def _repr_(self):
        return "A (breaking) base object"

    def a_realization(self):
        return self.Foo()

    class Foo(CombinatorialFreeModule, BindableClass):
        def __init__(self, base):
            self._basis_name = "Foo"
            CombinatorialFreeModule.__init__(self, base.base_ring(), ['a','b','c'],
                                             prefix='F',
                                             category=Bases(base))

    class Bar(CombinatorialFreeModule, BindableClass):
        def __init__(self, base):
            self._basis_name = "Bar"
            CombinatorialFreeModule.__init__(self, base.base_ring(), ['a','b','c'],
                                             prefix='B',
                                             category=Bases(base))

class Bases(Category_realization_of_parent):
    def __init__(self, base):
        Category_realization_of_parent.__init__(self, base)

    def _repr_(self):
        return "Category of bases of {}".format(self.base())

    def super_categories(self):
        return [ModulesWithBasis(self.base().base_ring()), Realizations(self.base())]

    class ParentMethods:
        def _repr_(self):
            return "{} in the {} basis".format(self.realization_of(), self._basis_name)

