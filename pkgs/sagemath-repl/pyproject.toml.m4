[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
                    ')]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-repl"
description = "Sage: Open Source Mathematics Software: System and software environment"
license = {text = "GNU General Public License (GPL) v2 or later"}
authors = [{name = "The Sage Developers", email = "sage-support@googlegroups.com"}]
classifiers = [
    "Development Status :: 6 - Mature",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
    "Operating System :: POSIX",
    "Operating System :: MacOS :: MacOS X",
    "Programming Language :: Python :: 3 :: Only",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: Implementation :: CPython",
    "Topic :: Scientific/Engineering :: Mathematics",
]
urls = {Homepage = "https://www.sagemath.org"}
requires-python = ">=3.8, <3.12"
dependencies = [esyscmd(`sage-get-system-packages install-requires-toml sagemath_objects')]
dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[tool.setuptools]
py-modules = [
    "sage.all__sagemath_repl",
    "sage.misc.all__sagemath_repl",
    "sage.misc.banner",
    "sage.misc.sagedoc",
    "sage.misc.sage_input",
    "sage.misc.sage_eval",
]
packages = [
    "sage.doctest",
    "sage.repl",
    "sage.repl.display",
    "sage.repl.ipython_kernel",
    "sage.repl.rich_output",
]
script-files = [
    # Other scripts that should be in the path also for OS packaging of sage:
    "bin/sage-eval",
    # Included because it is useful for doctesting/coverage testing user scripts too:
    "bin/sage-runtests",
    "bin/sage-fixdoctests",
    "bin/sage-coverage",
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    "bin/sage-cachegrind",
    "bin/sage-callgrind",
    "bin/sage-massif",
    "bin/sage-omega",
    "bin/sage-valgrind",
    "bin/sage-cleaner",
    # Uncategorized scripts in alphabetical order
    "bin/sage-inline-fortran",
    "bin/sage-ipynb2rst",
    "bin/sage-ipython",
    "bin/sage-notebook",
    "bin/sage-preparse",
    "bin/sage-run",
    "bin/sage-run-cython",
    "bin/sage-startuptime.py",
]
include-package-data = false

[tool.setuptools.package-data]
"sage.doctest" = ["tests/*"]
"sage.repl.rich_output" = ["example*"]

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
