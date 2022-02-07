sage_legacy_install_cleaner
===========================

Running `setup.py install` on this pseudo-distribution package removes
modules of the Sage library that are no longer present in the source
tree from the installation in `site-packages`.

Sage-the-distribution runs it after installing the Sage library using
the legacy distutils-based installation method (when `configure
--enable-editable` is not in use).

`sage_legacy_install_cleaner` is not intended for distribution on PyPI
and will also not be needed in downstream packaging of SageMath.
