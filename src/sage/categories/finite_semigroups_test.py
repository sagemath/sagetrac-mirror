class TestFiniteSemigroup():
    """
    Tests for finite semigroups.
    """

    @staticmethod
    def category_instances():
        from sage.categories.examples.finite_semigroups import LeftRegularBand
        return [
            LeftRegularBand(alphabet = ('a','b')), 
            LeftRegularBand(alphabet = ('a','b','c')),
            LeftRegularBand(alphabet = ('a','b','c', 'd'))
        ]

    def test_associativity(self, set_elements, max_runs):
        """
        Test associativity of multiplication of elements.
        """
        from sage.misc.misc import some_tuples
        for x, y, z in some_tuples(set_elements, 3, max_runs):
            assert (x * y) * z == x * (y * z)
