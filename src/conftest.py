# pyright: strict

"""Configuration and fixtures for pytest.

This file configures pytest and provides some global fixtures.
See https://docs.pytest.org/en/latest/index.html for more details.

At the moment, Sage is not yet using any tests based on pytest.
"""

from __future__ import annotations
from typing import Any
import pytest
from _pytest.python import Metafunc, PyCollector
import inspect

import sage.all # TODO: Remove once categories can be used without previous sage.all import
from sage.categories.sets_cat import Sets
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.infinite_enumerated_sets import InfiniteEnumeratedSets
from sage.categories.sets_cat_test import SetsTests
from sage.categories.enumerated_sets_test import EnumeratedSetsTests
from sage.categories.finite_enumerated_sets_test import FiniteEnumeratedSetsTests
from sage.categories.infinite_enumerated_sets_test import InfiniteEnumeratedSetsTests

# Dictionary relating the category to its test class
categories_with_tests = {
    Sets: SetsTests,
    EnumeratedSets: EnumeratedSetsTests,
    FiniteEnumeratedSets: FiniteEnumeratedSetsTests,
    InfiniteEnumeratedSets: InfiniteEnumeratedSetsTests
}


def pytest_pycollect_makeitem(collector: PyCollector, name: str, obj: type):
    if inspect.isclass(obj) and "category_instances" in dir(obj):
        # Enrich test class by functions from the corresponding category test class
        for category, category_test_class in categories_with_tests.items():
            if category() in obj.category_instances()[0].category().all_super_categories():
                methods = [method for method in dir(category_test_class)
                           if not method.startswith('__') and callable(getattr(category_test_class, method))]

                for method in methods:
                    setattr(obj, method, getattr(category_test_class, method))

        return pytest.Class.from_parent(collector, name=name, obj=obj)
    else:
        return None


def pytest_generate_tests(metafunc: Metafunc):
    # Add support for the "category_instance" parametrization
    categor_instance_fixure_name = "category_instance"
    if categor_instance_fixure_name in metafunc.fixturenames:
        params = []
        for category_instance in metafunc.cls.category_instances():
            params.append(pytest.param(category_instance, id=str(category_instance)))
        metafunc.parametrize(categor_instance_fixure_name, params)


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    # TODO: Remove this workaround as soon as sage objects can be used without previous sage.all import
    import sage.all  # type: ignore # implicitly used below by calling locals()
    doctest_namespace.update(**locals())


@pytest.fixture
def max_runs():
    return 20
