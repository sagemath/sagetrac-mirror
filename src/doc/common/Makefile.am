
all-local: build_options.py

CLEANFILES = build_options.pyc

EXTRA_DIST = \
    builder.py \
    build_options.py.in \
    conf.py \
    custom-sphinx-build.py \
    __init__.py \
    inventory_builder.py \
    multidocs.py \
    python.inv \
    sage_autodoc.py \
    themes \
    update-python-inv.sh

dist-hook:
	-find themes -name *~ -delete
	-find themes -name .swp -delete
