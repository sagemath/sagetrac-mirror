# pyright: strict

"""Configuration and fixtures for pytest.

This file configures pytest and provides some global fixtures.
See https://docs.pytest.org/en/latest/index.html for more details.

At the moment, Sage is not yet using any tests based on pytest.
"""

from __future__ import annotations

import sys
from importlib.abc import MetaPathFinder, SourceLoader
from importlib.machinery import ModuleSpec
from importlib.util import spec_from_loader
from pathlib import Path
from typing import Any, List, Sequence

import pytest


def pytest_collection_modifyitems(
    session: pytest.Session, config: pytest.Config, items: List[pytest.Item]
):
    """
    This hook is called after collection has been performed, and can be used to
    modify the list of items that will be run.

    See `pytest documentation <https://docs.pytest.org/en/latest/reference/reference.html#std-hook-pytest_collection_modifyitems>`_.
    """
    skip_as_false_positive = pytest.mark.skip(
        reason="Skipping this because its not a pytest test but an ordinary"
        + "method that happens to start with 'test_'"
    )
    for item in items:
        # Add a mark to all tests that should be skipped
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


def pytest_collect_file(parent: pytest.Collector, file_path: Path):
    if file_path.suffix == ".sage":
        return pytest.Module.from_parent(parent, path=file_path)  # type: ignore # pytest typing is inclomplete


class SageSourceFinder(MetaPathFinder):
    """A MetaPathFinder that finds modules in Sage files."""

    def find_spec(
        self,
        fullname: str,
        path: Sequence[str | bytes] | None = None,
        target: Any | None = None,
    ) -> ModuleSpec | None:
        if not path:
            path = sys.path
        for directory in path:
            if not isinstance(directory, str):
                directory = directory.decode()

            module_name = fullname.split(".")[-1]
            sage_file = Path(directory, f"{module_name}.sage")
            if sage_file.is_file():
                return spec_from_loader(
                    module_name,
                    SageSourceLoader(sage_file),
                    origin=sage_file.as_uri(),
                )


class SageSourceLoader(SourceLoader):
    """
    A SourceLoader that can load modules from a Sage file.
    """

    def __init__(self, path: Path):
        self.path = path

    def get_data(self, path: str | bytes) -> bytes:
        from sage.repl.preparse import preparse_file

        with open(path) as file:
            # We need to import the sage library first
            preamble = "from sage.all_cmdline import *\n"
            content = preparse_file(file.read())
            return (preamble + content).encode()

    def get_filename(self, fullname: str) -> str:
        return self.path.absolute().as_posix()


sys.meta_path.append(SageSourceFinder())


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    import sage.all  # type: ignore # implicitly used below by calling locals()
    doctest_namespace.update(**locals())
