import pytest

class SetsTests:
    """
    Generic tests for the category of sets.
    """
    @pytest.fixture
    def set_elements(self, category_instance):
        return category_instance.some_elements()

    def test_an_element(self, category_instance):
        """
        Test that self.an_element() is in self
        """
        assert category_instance.an_element() in category_instance

    def test_an_element_idempotent(self, category_instance):
        """
        Test that element construction is not idempotent
        """
        an_element = category_instance.an_element()
        if category_instance.is_parent_of(an_element):
            assert category_instance(an_element) == an_element
        else: # Allows self(an_element) to fails for facade parent.
            try:
                rebuilt_element = category_instance(an_element)
            except NotImplementedError:
                print(f"The set {category_instance} doesn't seems to implement __call__; skipping test of construction idempotency")
            else:
                assert rebuilt_element == an_element

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

    def test_construction(self, category_instance):
        """
        Test that the construction returned by self really yields self.

        :meth:`construction` either returns None or a pair ``(F,O)``,
        and if it returns the latter, then it is supposed that ``F(O)==self`.
        """
        FO = category_instance.construction()
        if FO is None:
            return
        assert FO[0](FO[1]) == category_instance

    def test_element_eq_reflexive(self, set_elements):
        """
        Test that equality is reflexive.
        """
        for element in set_elements:
            assert element == element

    def test_elements_eq_symmetric(self, set_elements, max_runs):
        """
        Test that equality is symmetric.
        """
        S = list(set_elements) + [None, 0]
        from sage.misc.misc import some_tuples
        for x, y in some_tuples(S, 2, max_runs):
            assert (x==y) == (y==x), f"non symmetric equality: {print_compare(x, y)} but {print_compare(y, x)}"

    def test_elements_eq_transitive(self, set_elements, max_runs):
        """
        Test that equality is transitive.
        """
        S = set_elements
        if (len(S)+2)**3 <= max_runs:
            S = list(S) + [None, 0]
        else:
            from random import sample
            from sage.rings.integer import Integer
            S = sample(S, Integer(max_runs).nth_root(3, truncate_mode=1)[0] - 2) + [None, 0]

        for x in S:
            for y in S:
                if not x == y:
                    continue
                for z in S:
                    if not y == z:
                        continue

                    assert x == z, f"non transitive equality:\n {print_compare(x, y)} and {print_compare(y, z)} but {print_compare(x, z)}"

    def test_elements_neq(self, set_elements, max_runs):
        """
        Test that ``==`` and ``!=`` are consistent.
        """
        S = list(set_elements) + [None, 0]

        from sage.misc.misc import some_tuples
        for x,y in some_tuples(S, 2, max_runs):
            assert (x == y) != (x != y), f"__eq__ and __ne__ inconsistency:\n {x} == {y} returns {x == y}  but  {x} != {y} returns {x != y}"
