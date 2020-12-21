import pytest
from sage.rings.infinity import infinity

class InfiniteEnumeratedSetsTests:
    """
    Generic tests for the category of infinite enumerated sets.
    """    
    def test_enumerated_set_iter_cardinality(self, category_instance):
        """
        Test that the methods :meth:`.cardinality` and
        :meth:`.__iter__` are consistent.

        For infinite enumerated sets:

        * :meth:`.cardinality` is supposed to return `infinity`

        * :meth:`.list` is supposed to raise a ``NotImplementedError``.
        """
        assert category_instance.cardinality() == infinity
        with pytest.raises(NotImplementedError):
            category_instance.list
