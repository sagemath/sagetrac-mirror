[build-system]
# Minimum requirements for the build system to execute.
requires = [
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools_wheel     \
                    ')]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-environment"
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
dependencies = ["esyscmd(`sage-get-system-packages install-requires | sed \"2,\\$s/^/    /;\"')dnl"]
dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
# sage.env can optionally use sage_conf
conf = ["esyscmd(`sage-get-system-packages install-requires sage_conf')"]
# For "sage --docbuild"
docbuild = ["esyscmd(`sage-get-system-packages install-requires sage_docbuild')"]
# For "sage", "sage -t", ...
sage = ["esyscmd(`sage-get-system-packages install-requires sagelib')"]
# For "sage --cython"
cython = ["esyscmd(`sage-get-system-packages install-requires cython')"]
# For "sage --pytest"
pytest = ["esyscmd(`sage-get-system-packages install-requires pytest')"]
# For "sage --rst2ipynb"
rst2ipynb = ["esyscmd(`sage-get-system-packages install-requires rst2ipynb')"]
# For "sage --sws2rst"
sws2rst = ["esyscmd(`sage-get-system-packages install-requires sage_sws2rst')"]

[tool.setuptools]
py-modules = [
    "sage.all__sagemath_environment",
    "sage.env",
    "sage.version",
    "sage.misc.all__sagemath_environment",
    "sage.misc.package",
    "sage.misc.package_dir",
    "sage.misc.temporary_file",
    "sage.misc.viewer",
]
packages = ["sage.features"]
script-files = [
    # The sage script
    "bin/sage",
    # Auxiliary scripts for setting up the environment
    "bin/sage-env",
    "bin/sage-num-threads.py",
    "bin/sage-venv-config",
    "bin/sage-version.sh",
    # Auxiliary script for invoking Python in the Sage environment
    "bin/sage-python",
    # Not included:
    # - bin/sage-env-config                  -- installed by sage_conf
    # - bin/sage-env-config.in               -- not to be installed
    # - bin/sage-run, bin/sage-runtests, ... -- installed by sagemath-repl
    # - bin/sage-ipython  -- uses sage.repl, so installed by sagemath-repl
]
include-package-data = false

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
