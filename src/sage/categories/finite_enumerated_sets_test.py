class FiniteEnumeratedSetsTests:
    """
    Generic tests for the category of finite enumerated sets.
    """

    def test_enumerated_set_iter_cardinality(self, category_instance, max_runs):
        """
        Checks that the methods :meth:`.cardinality` and
        :meth:`.__iter__` are consistent. Also checks that
        :meth:`.cardinality` returns an ``Integer``.

        For efficiency reasons, those tests are not run if
        :meth:`.cardinality` is
        :meth:`._cardinality_from_iterator`, or if ``category_instance`` is too
        big.
        """
        if category_instance.cardinality != category_instance._cardinality_from_iterator:
            card = category_instance.cardinality()
            if card <= max_runs:
                assert card == category_instance._cardinality_from_iterator()
