import pytest

class SetsTests:
    """
    Generic tests for the category of sets.
    """

    @pytest.fixture
    def set_elements(self, category_instance):
        return category_instance.some_elements()

    def test_element_eq_reflexive(self, set_elements):
        """
        Test that equality is reflexive.
        """
        for element in set_elements:
            assert element == element


    def test_cardinality_return_type(self, category_instance):
        """
        Test that the method :meth:`cardinality` has the correct return type.
        """
        try:
            cardinality = category_instance.cardinality()
        except (AttributeError,NotImplementedError):
            return

        from sage.structure.element import parent
        from sage.rings.infinity import Infinity
        from sage.rings.integer_ring import ZZ
        assert cardinality is Infinity or parent(cardinality) is ZZ
