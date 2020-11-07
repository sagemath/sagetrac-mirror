
class SetsTests:
    """
    Tests for the category of sets.
    """

    def test_element_eq_reflexive(self, set_element):
        """
        Test that equality is reflexive.
        """
        assert set_element == set_element
