class SemigroupTests:
    """
    Generic tests for semigroups.
    """

    def test_associativity(self, set_elements, max_runs):
        """
        Test associativity of multiplication of elements.
        """
        from sage.misc.misc import some_tuples

        for x, y, z in some_tuples(set_elements, 3, max_runs):
            assert (x * y) * z == x * (y * z)
