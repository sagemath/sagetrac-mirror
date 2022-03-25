# pyright: strict

"""Configuration and fixtures for pytest.

This file configures pytest and provides some global fixtures.
See https://docs.pytest.org/en/latest/index.html for more details.

At the moment, Sage is not yet using any tests based on pytest.
"""

from __future__ import annotations
from typing import Any, List
import pytest


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: List[pytest.Item]
):
    skip_as_false_positive = pytest.mark.skip(
        reason="Skipping this because its not a pytest test but an ordinary method that happens to start with 'test_'"
    )
    for item in items:
        if item.name in [
            "test_relation_maxima",
            "test_b2_local",
            "test_b2_global",
            "test_a1a3_local",
            "test_a1a3_global",
            "test_rst_global",
            "test_comparison",
            "test_signed_infinity",
            "test_pickle",
        ]:
            item.add_marker(skip_as_false_positive)


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    import sage.all  # type: ignore # implicitly used below by calling locals()
    doctest_namespace.update(**locals())
