import pytest

from sage.categories.examples.finite_semigroups import LeftRegularBand
from .test_categories import CategoryTests

class TestFiniteSemigroup(CategoryTests):

    test_elements = LeftRegularBand(alphabet = ('a','b','c')).some_elements()

    @pytest.fixture(params=test_elements)
    def set_element(self, request):
        # This will invoke test_elements_eq_reflexive (and all other test methods in this class that use set_element test fixture)
        # with ['a', 'b', 'c', 'ab', 'ac', 'ba', 'bc', 'ca', 'cb', 'abc']
        return request.param
