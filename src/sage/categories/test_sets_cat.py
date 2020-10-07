
class GenericTests:
    def test_element_eq_reflexive(self, set_element):
        """
        Generic test on the equality of elements.

        Test that ``==`` is reflexive.
        """
        assert set_element == set_element

    # TODO: Add other generic test methods here