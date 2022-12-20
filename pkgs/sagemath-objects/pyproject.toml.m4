[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        sage_setup     \
        sagemath_environment \
        cython         \
        gmpy2          \
        cysignals      \
                    ')]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-objects"
description = "Sage: Open Source Mathematics Software: Sage objects, elements, parents, categories, coercion, metaclasses"
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
dependencies = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        gmpy2                         \
        cysignals                     \
    ')]
dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
# Currently we do not use the sage doctester to test sagemath-objects,
# so we do not list sagemath-repl here.
test = []

[tool.setuptools]
include-package-data = false

[tool.setuptools.package-data]
"sage.cpython" = [
    "pyx_visit.h",
    "string_impl.h",
    "cython_metaclass.h",
    "python_debug.h",
]
"sage.rings" = ["integer_fake.h"]

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
