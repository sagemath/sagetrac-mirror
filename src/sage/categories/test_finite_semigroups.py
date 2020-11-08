class TestFiniteSemigroup():
    """
    Tests for finite semigroups.
    """

    @staticmethod
    def category_instances():
        from sage.categories.examples.finite_semigroups import LeftRegularBand
        return [
            LeftRegularBand(alphabet = ('a','b')), 
            LeftRegularBand(alphabet = ('a','b','c'))
        ]
