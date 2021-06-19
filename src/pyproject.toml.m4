[build-system]
# Minimum requirements for the build system to execute.
requires = [
    # Some version of sage-conf is required.
    # Note that PEP517/518 have no notion of optional build dependencies:
    # https://github.com/pypa/pip/issues/6144
    'sage-conf',
    esyscmd(`sage-get-system-packages install-requires-toml \
        setuptools     \
        wheel          \
        cypari         \
        cysignals      \
        cython         \
        gmpy2          \
        jinja2         \
        jupyter_core   \
        numpy          \
        pkgconfig      \
        pplpy          \
                    ')]
build-backend = "setuptools.build_meta"
