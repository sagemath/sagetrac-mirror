[build-system]
# Minimum requirements for the build system to execute.
requires = [
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional sage_spkg dependencies:
    # https://github.com/pypa/pip/issues/6144
     esyscmd(`sage-get-system-packages install-requires-toml \
        sage_conf      \
        setuptools     \
        wheel          \
        sage_setup     \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        numpy          \
        pkgconfig      \
        pplpy          \
        memory_allocator \
                    ')]
build-backend = "setuptools.build_meta"

[project]
name = "sagemath-standard"
description = "Sage: Open Source Mathematics Software: Standard Python Library"
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
        sage_conf                     \
        six                           \
    ')
dnl From build/pkgs/sagelib/dependencies
    esyscmd(`sage-get-system-packages install-requires-toml \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        lrcalc_python  \
        memory_allocator \
        numpy          \
        pkgconfig      \
        pplpy          \
        primecountpy   \
        requests       \
    ')
dnl From Makefile.in: SAGERUNTIME
    esyscmd(`sage-get-system-packages install-requires-toml \
        ipython        \
        pexpect        \
    ')
dnl From Makefile.in: DOC_DEPENDENCIES
    esyscmd(`sage-get-system-packages install-requires-toml \
        sphinx         \
        networkx       \
        scipy          \
        sympy          \
        matplotlib     \
        pillow         \
        mpmath         \
        ipykernel      \
        jupyter_client \
        ipywidgets     \
    ')
dnl Other Python packages that are standard spkg, used in doctests
    esyscmd(`sage-get-system-packages install-requires-toml \
        fpylll         \
    ')
dnl pycryptosat  # Sage distribution installs it as part of cryptominisat. According to its README on https://pypi.org/project/pycryptosat/: "The pycryptosat python package compiles while compiling CryptoMiniSat. It cannot be compiled on its own, it must be compiled at the same time as CryptoMiniSat."
dnl Packages with important upper version bounds
    esyscmd(`sage-get-system-packages install-requires-toml \
        ptyprocess     \
    ')
]

dynamic = ["version"]

[project.readme]
file = "README.rst"
content-type = "text/x-rst"

[project.optional-dependencies]
r = [esyscmd(`sage-get-system-packages install-requires-toml rpy2')]

[tool.setuptools]
script-files = [
    # The sage script
    "bin/sage",
    # Other scripts that should be in the path also for OS packaging of sage:
    "bin/sage-eval",
    # Included because it is useful for doctesting/coverage testing user scripts too:
    "bin/sage-runtests",
    "bin/sage-fixdoctests",
    "bin/sage-coverage",
    # The following is deprecated but might still be used in user package install scripts
    "bin/sage-cython",
    # Helper scripts invoked by sage script
    # (they would actually belong to something like libexec)
    "bin/sage-cachegrind",
    "bin/sage-callgrind",
    "bin/sage-massif",
    "bin/sage-omega",
    "bin/sage-valgrind",
    "bin/sage-venv-config",
    "bin/sage-version.sh",
    "bin/sage-cleaner",
    # Only makes sense in sage-the-distribution. TODO: Move to another installation script.
    "bin/sage-list-packages",
    # Uncategorized scripts in alphabetical order
    "bin/math-readline",
    "bin/sage-env",
    # sage-env-config -- installed by sage_conf
    # sage-env-config.in -- not to be installed
    "bin/sage-grep",
    "bin/sage-grepdoc",
    "bin/sage-inline-fortran",
    "bin/sage-ipynb2rst",
    "bin/sage-ipython",
    "bin/sage-notebook",
    "bin/sage-num-threads.py",
    "bin/sage-preparse",
    "bin/sage-python",
    "bin/sage-rebase.bat",
    "bin/sage-rebase.sh",
    "bin/sage-rebaseall.bat",
    "bin/sage-rebaseall.sh",
    "bin/sage-run",
    "bin/sage-run-cython",
    "bin/sage-startuptime.py",
    "bin/sage-update-version",
]
license-files = ["LICENSE.txt"]
include-package-data = false

[tool.setuptools.package-data]
"sage.libs.gap" = ["sage.gaprc"]
"sage.interfaces" = ["sage-maxima.lisp"]
"sage.doctest" = ["tests/*"]
"sage.repl.rich_output" = ["example*"]
sage = [
    "ext_data/*",
    "ext_data/kenzo/*",
    "ext_data/singular/*",
    "ext_data/singular/function_field/*",
    "ext_data/images/*",
    "ext_data/doctest/*",
    "ext_data/doctest/invalid/*",
    "ext_data/gap/*",
    "ext_data/gap/joyner/*",
    "ext_data/mwrank/*",
    "ext_data/notebook-ipython/*",
    "ext_data/nbconvert/*",
    "ext_data/graphs/*",
    "ext_data/pari/*",
    "ext_data/pari/dokchitser/*",
    "ext_data/pari/buzzard/*",
    "ext_data/pari/simon/*",
    "ext_data/magma/*",
    "ext_data/magma/latex/*",
    "ext_data/magma/sage/*",
    "ext_data/valgrind/*",
    "ext_data/threejs/*",
]

[tool.setuptools.dynamic]
version = {file = ["VERSION.txt"]}
