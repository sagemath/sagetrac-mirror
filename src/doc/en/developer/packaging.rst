.. _chapter-packaging:

==========================
Packaging Third-Party Code
==========================

One of the mottoes of the Sage project is to not reinvent the wheel: If
an algorithm is already implemented in a well-tested library then
consider incorporating that library into Sage. The current list of
available packages are the subdirectories of ``SAGE_ROOT/build/pkgs/``.
The management of packages is done through a bash script located in
``SAGE_ROOT/build/bin/sage-spkg``. This script is typically invoked by
giving the command::

    [user@localhost]$ sage -i <options> <package name>...

options can be:

- f: install a package even if the same version is already installed
- s: do not delete temporary build directory
- c: after installing, run the test suite for the spkg. This should
  override the settings of ``SAGE_CHECK`` and ``SAGE_CHECK_PACKAGES``.
- d: only download the package

Not all packages are built by default, they are divided into standard,
optional and experimental ones. Standard packages are built by default
and have stringent quality requirements. Optional packages are subject
to the same requirements (although these are tested not as much,
so in practice they might break more often than standard packages).
For experimental packages, the bar is much lower: even if there are some
problems, the package can still be accepted.

The section :ref:`section-directory-structure` describes the structure
of each individual package in ``SAGE_ROOT/build/pkgs``. In section
:ref:`section-manual-build` we see how you can install and test a new
spkg that you or someone else wrote. Finally,
:ref:`section-inclusion-procedure` explains how to submit a new package
for inclusion in the Sage source code.


.. _section-directory-structure:

Directory Structure
===================

Third-party packages in Sage consists of two parts:

#. The tarball as it is distributed by the third party, or as close as
   possible. Valid reasons for modifying the tarball are deleting
   unnecessary files to keep the download size manageable or
   regenerating auto-generated files if necessary. But the actual code
   must be unmodified. See also :ref:`section-spkg-src`.

#. The build scripts and associated files are in a subdirectory
   ``SAGE_ROOT/build/pkgs/package``, where you replace ``package``
   with a lower-case version of the upstream project name.

As an example, let us consider a hypothetical FoO project. They
(upstream) distribute a tarball ``foo-1.3.tar.gz`` (that will be
automatically placed in ``SAGE_ROOT/upstream`` during the installation
process). To package it in Sage, we create a subdirectory containing the
following::

    SAGE_ROOT/build/pkgs/foo
    |-- patches
    |   |-- bar.patch
    |   `-- baz.patch
    |-- checksums.ini
    |-- package-version.txt
    |-- spkg-check
    |-- spkg-install
    |-- spkg-src
    `-- SPKG.txt

When installing Sage these files are used to patch the tarball and to
start the build and install process of the package.

We discuss the individual files in the following.


.. _section-spkg-install:

Install Script
--------------

The ``spkg-install`` file is a shell script installing the package,
with ``PACKAGE_NAME`` replaced by the the package name. In the best
case, the upstream project can simply be installed by the usual
configure / make / make install steps. In that case, the build script
would simply consist of::

    #!/usr/bin/env bash

    cd src

    ./configure --prefix="$SAGE_LOCAL" --libdir="$SAGE_LOCAL/lib"
    if [ $? -ne 0 ]; then
        echo >&2 "Error configuring PACKAGE_NAME."
        exit 1
    fi

    $MAKE
    if [ $? -ne 0 ]; then
        echo >&2 "Error building PACKAGE_NAME."
        exit 1
    fi

    $MAKE -j1 install
    if [ $? -ne 0 ]; then
        echo >&2 "Error installing PACKAGE_NAME."
        exit 1
    fi


Note that the top-level directory inside the tarball is renamed to
``src`` before calling the ``spkg-install`` script, so you can just use
``cd src`` instead of ``cd foo-1.3``.

If there is any meaningful documentation included but not installed by
``make install``, then you can add something like the following to
install it::

    if [ "$SAGE_SPKG_INSTALL_DOCS" = yes ] ; then
        $MAKE doc
        if [ $? -ne 0 ]; then
            echo >&2 "Error building PACKAGE_NAME docs."
            exit 1
        fi
        mkdir -p "$SAGE_LOCAL/share/doc/PACKAGE_NAME"
        cp -R doc/* "$SAGE_ROOT/local/share/doc/PACKAGE_NAME"
    fi




.. _section-spkg-check:

Self-Tests
----------

The ``spkg-check`` file is an optional, but highly recommended, script
to run self-tests of the package. It is run after building and
installing if the ``SAGE_CHECK`` environment variable is set, see the
Sage installation guide. Ideally, upstream has some sort of tests suite
that can be run with the standard ``make check`` target. In that case,
the ``spkg-check`` script would simply contain::

    #!/usr/bin/env bash

    cd src
    $MAKE check


.. _section-spkg-versioning:

Package Versioning
------------------

The ``package-version.txt`` file containts just the version. So if
upstream is ``foo-1.3.tar.gz`` then the package version file would only
contain ``1.3``.

If the upstream package is taken from some revision other than a stable
version, you should use the date at which the revision is made, e.g. the
Singular package ``20090818`` is made with the revision as of
2009-08-18.

If you made any changes to the upstream tarball (see
:ref:`section-directory-structure` for allowable changes) then you
should append a ``.p1`` to the version. If you make further changes,
increase the patch level as necessary. So the different versions would
be ``1.3``, ``1.3.p1``, ``1.3.p2``, ...


.. _section-spkg-SPKG-txt:

The SPKG.txt File
-----------------

The ``SPKG.txt`` file should follow this pattern::

     = PACKAGE_NAME =

     == Description ==

     What does the package do?

     == License ==

     What is the license? If non-standard, is it GPLv3+ compatible?

     == SPKG Maintainers ==

     * Mary Smith
     * Bill Jones
     * Leonhard Euler

     == Upstream Contact ==

     Provide information for upstream contact.

     == Dependencies ==

     Put a bulleted list of dependencies here:

     * python
     * readline

     == Special Update/Build Instructions ==

     If the tarball was modified by hand and not via a spkg-src
     script, describe what was changed.


with ``PACKAGE_NAME`` replaced by the the package name. Legacy
``SPKG.txt`` files have an additional changelog section, but this
information is now kept in the git repository.


.. _section-spkg-patching:

Patching Sources
----------------

Actual changes to the source code must be via patches, which should be
placed in the ``patches`` directory. GNU patch is distributed with
Sage, so you can rely on it being available. Patches must include
documentation in their header (before the first diff hunk), so a
typical patch file should look like this::

    Add autodoc_builtin_argspec config option

    Following the title line you can add a multi-line description of
    what the patch does, where you got it from if you did not write it
    yourself, if they are platform specific, if they should be pushed
    upstream, etc...

    diff -dru Sphinx-1.2.2/sphinx/ext/autodoc.py.orig Sphinx-1.2.2/sphinx/ext/autodoc.py
    --- Sphinx-1.2.2/sphinx/ext/autodoc.py.orig  2014-03-02 20:38:09.000000000 +1300
    +++ Sphinx-1.2.2/sphinx/ext/autodoc.py  2014-10-19 23:02:09.000000000 +1300
    @@ -1452,6 +1462,7 @@

         app.add_config_value('autoclass_content', 'class', True)
         app.add_config_value('autodoc_member_order', 'alphabetic', True)
    +    app.add_config_value('autodoc_builtin_argspec', None, True)
         app.add_config_value('autodoc_default_flags', [], True)
         app.add_config_value('autodoc_docstring_signature', True, True)
         app.add_event('autodoc-process-docstring')

Patches to files in ``src/`` need to be applied in ``spkg-install``,
that is, if there are any patches then your ``spkg-install`` script
should contain a section like this::

    for patch in ../patches/*.patch; do
        [ -r "$patch" ] || continue  # Skip non-existing or non-readable patches
        patch -p1 <"$patch"
        if [ $? -ne 0 ]; then
            echo >&2 "Error applying '$patch'"
            exit 1
        fi
    done

which applies the patches to the sources.


.. _section-spkg-src:

Modified Tarballs
-----------------

The ``spkg-src`` file is optional and only to document how the upstream
tarball was changed. Ideally it is not modified, then there would be no
``spkg-src`` file present either.

However, if you really must modify the upstream tarball then it is
recommended that you write a script, called ``spkg-src``, that makes the
changes. This not only serves as documentation but also makes it easier
to apply the same modifications to future versions.


Checksums
---------

The ``checksums.ini`` file contains checksums of the upstream tarball.
It is autogenerated, so you just have to place the upstream tarball in
the ``SAGE_ROOT/upstream/`` directory and run::

    [user@localhost]$ sage -sh sage-fix-pkg-checksums


.. _section-manual-build:

Manual package build and installation
=====================================

At this stage you have a new tarball that is not yet distributed with
Sage (``foo-1.3.tar.gz`` in the example of section
:ref:`section-directory-structure`). Now you need to manually place it
in the ``SAGE_ROOT/upstream/`` directory. Then you can run the
installation via::

    [user@localhost]$ sage -i package_name

or::

    [user@localhost]$ sage -f package_name

to force a reinstallation. If your package contains a ``spkg-check``
script (see :ref:`section-spkg-check`) it can be run with::

    [user@localhost]$ sage -i -c package_name

or::

    [user@localhost]$ sage -f -c package_name

If all went fine, open a ticket, put a link to the original tarball in
the ticket and upload a branch with the code under
``SAGE_ROOT/build/pkgs``.


.. _section-inclusion-procedure:

Inclusion Procedure for New and Updated Packages
================================================

Packages that are not part of Sage will first become optional or
experimental (the latter if they will not build on all supported
systems). After they have been in optional for some time without
problems they can be proposed to be included as standard packages in
Sage.

To propose a package for optional/experimental inclusion please open a
trac ticket with the respective ``Component:`` field set to either
``packages:experimental`` or ``packages:optional``. The associated code
requirements are described in the following sections.

After the ticket was reviewed and included, optional packages stay in
that status for at least a year, after which they can be proposed to be
included as standard packages in Sage. For this a trac ticket is opened
with the ``Component:`` field set to ``packages:standard``. Note that
the script in ``SAGE_ROOT/build/make/deps`` is called when building Sage so
please include the build command for your standard package there. Then
make a proposal in the Google Group ``sage-devel``.

Upgrading packages to new upstream versions or with additional patches
includes opening a ticket in the respective category too, as described
above.

License Information
-------------------

If you are patching a standard Sage spkg, then you should make sure that
the license information for that package is up-to-date, both in its
``SPKG.txt`` file and in the file ``SAGE_ROOT/COPYING.txt``.  For
example, if you are producing an spkg which upgrades the vanilla source
to a new version, check whether the license changed between versions.

Prerequisites for New Standard Packages
---------------------------------------

For a package to become part of Sage's standard distribution, it
must meet the following requirements:

- **License**. For standard packages, the license must be compatible
  with the GNU General Public License, version 3. The Free Software
  Foundation maintains a long list of `licenses and comments about
  them <http://www.gnu.org/licenses/license-list.html>`_.

- **Build Support**. The code must build on all the `fully supported
  platforms
  <http://wiki.sagemath.org/SupportedPlatforms#Fully_supported>`_.

  A standard package should also work on all the platforms where Sage
  is `expected to work
  <http://wiki.sagemath.org/SupportedPlatforms#Expected_to_work>`_ and
  on which Sage `almost works
  <http://wiki.sagemath.org/SupportedPlatforms#Almost_works>`_ but
  since we don't fully support these platforms and often lack the
  resources to test on them, you are not expected to confirm your
  packages works on those platforms.

- **Quality**. The code should be "better" than any other available
  code (that passes the two above criteria), and the authors need to
  justify this. The comparison should be made to both Python and other
  software. Criteria in passing the quality test include:

  - Speed

  - Documentation

  - Usability

  - Absence of memory leaks

  - Maintainable

  - Portability

  - Reasonable build time, size, dependencies

- **Previously an optional package**. A new standard package must have
  spent some time as an optional package. Or have a good reason why
  this is not possible.

- **Refereeing**. The code must be refereed, as discussed in
  :ref:`chapter-sage-trac`.


