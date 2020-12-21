import pytest

class EnumeratedSetsTests:
    """
    Generic tests for the category of enumerated sets.
    """

    def test_enumerated_set_contains(self, category_instance, max_runs):
        """
        Test that the methods :meth:`.__contains__` and :meth:`.__iter__` are consistent.
        """
        i = 0
        for w in category_instance:
            assert w in category_instance

            i += 1
            if i > max_runs:
                return

    def test_enumerated_set_iter_list(self, category_instance, max_runs):
        """
        Test that the methods :meth:`.list` and :meth:`.__iter__` are consistent.

        .. NOTE::

            This test does nothing if the cardinality of the set
            is larger than the max_runs argument.
        """
        if category_instance.list != category_instance._list_default:
            # TODO: if self._cardinality is self._cardinality_from_iterator
            # we could make sure to stop the counting at
            # self.max_test_enumerated_set_loop
            if category_instance.cardinality() > max_runs:
                print("Enumerated set too big; skipping test; increase max_runs")
                return
            ls = category_instance.list()
            i = 0
            for obj in category_instance:
                assert obj == ls[i]
                i += 1
            assert i == len(ls)
