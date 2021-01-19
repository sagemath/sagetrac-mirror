[build-system]
# Minimum requirements for the build system to execute.
requires = [esyscmd(`sage-get-system-packages install-requires-toml \
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
