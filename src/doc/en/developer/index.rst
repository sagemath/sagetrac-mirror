.. _developers-guide:

======================================
Welcome to the Sage Developer's Guide!
======================================

Everybody who uses SageMath is encouraged to contribute something back to it at
some point.

Note that the SageMath project is a major integrating force in the mathematical
software landscape.
Contributing to SageMath does not necessarily mean to work on the Sage library.
Follow the ticket numbers to see how you can help; for example, by surveying software,
building or extending interfaces, or contributing to distribution packaging.

.. list-table::
   :widths: 20 40 40
   :header-rows: 1

   * -
     - **In the SageMath distribution**
     - Not in the SageMath distribution
   * - Unused math software
     - `SymEngine <../reference/spkg/symengine.html>`__,
       `cocoalib <../reference/spkg/cocoalib.html>`__,
       `isl <../reference/spkg/isl.html>`__,
       `CyLP <../reference/spkg/cylp.html>`__, ...
       :trac:`#33773 <33773>`
     - `sagemath.org links <https://www.sagemath.org/links.html>`__,
       `Open-source CAS <https://wiki.sagemath.org/OSCAS>`__,
       `SoftwareToIntegrate <https://wiki.sagemath.org/devel/SoftwareToIntegrate>`__,
       `SoftwareSurveys <https://wiki.sagemath.org/devel/SoftwareSurveys>`__ :trac:`#33725 <33725>`,
       :trac:`Python optimization packages <26511>` :trac:`#26511 <26511>`,
       `polyhedral geometry packages <https://wiki.sagemath.org/OptiPolyGeom>`__,
       `swMATH.org <https://swmath.org/>`__,
       math databases :trac:`#30914 <30914>`
   * - Upstream non-Python math software
     - `Packages <../reference/spkg/>`__:
       `GMP <../reference/spkg/gmp.html>`__,
       `FLINT <../reference/spkg/flint.html>`__ :trac:`#31408 <31408>`,
       `arb <../reference/spkg/arb.html>`__,
       `Singular <../reference/spkg/singular.html>`__,
       `GAP <../reference/spkg/gap.html>`__,
       `PARI <../reference/spkg/pari.html>`__,
       `Maxima <../reference/spkg/maxima.html>`__,
       `R <../reference/spkg/r.html>`__,
       `FPLLL <../reference/spkg/fplll.html>`__,
       `LinBox <../reference/spkg/linbox.html>`__,
       `GSL <../reference/spkg/gsl.html>`__,
       `polymake <../reference/spkg/polymake.html>`__,
       ...
     - `Free software that Sage interfaces to <https://wiki.sagemath.org/Interfaces#Interfaces_to_other_software_in_SageMath>`__:
       Macaulay2, Octave, ...;
       `Non-free software that Sage interfaces to <https://wiki.sagemath.org/Interfaces#Interfaces_to_other_software_in_SageMath>`__:
       Magma, Maple, Mathematica, CPLEX, Gurobi, SCIP, ...
   * - Upstream Python math software
     - `Packages <../reference/spkg/>`__:
       `NumPy <../reference/spkg/numpy.html>`__,
       `SciPy <../reference/spkg/scipy.html>`__,
       `SymPy <../reference/spkg/sympy.html>`__,
       `CVXOPT <../reference/spkg/cvxopt.html>`__,
       `NetworkX <../reference/spkg/networkx.html>`__, ...
     -
   * - Upstream distributions
     -
     - `Distributions providing system packages <https://wiki.sagemath.org/ReleaseTours/sage-9.1#Portability_improvements.2C_increased_use_of_system_packages>`__: homebrew, conda-forge, archlinux, ... :trac:`#27330 <27330>`
   * - Absorbed libraries
     - :trac:`GiNaC/Pynac <32386>` :trac:`#33401 <33401>`,
       :trac:`giacpy <29171>`,
       :trac:`sage_brial <30332>`
     -
   * - `SageMath library <../reference/index.html>`__
     - `Source <https://github.com/sagemath/sage/tree/develop/src/sage>`__:
       `sage.algebras.\* <https://github.com/sagemath/sage/tree/develop/src/sage/algebras>`__,
       `sage.categories.\* <https://github.com/sagemath/sage/tree/develop/src/sage/categories>`__,
       `sage.combinat.\* <https://github.com/sagemath/sage/tree/develop/src/sage/combinat>`__,
       `sage.crypto.\* <https://github.com/sagemath/sage/tree/develop/src/sage/crypto>`__,
       `sage.geometry.\* <https://github.com/sagemath/sage/tree/develop/src/sage/geometry>`__,
       `sage.graphs.\* <https://github.com/sagemath/sage/tree/develop/src/sage/graphs>`__,
       `sage.groups.\* <https://github.com/sagemath/sage/tree/develop/src/sage/groups>`__,
       `sage.manifolds.\* <https://github.com/sagemath/sage/tree/develop/src/sage/manifolds>`__,
       `sage.modular.\* <https://github.com/sagemath/sage/tree/develop/src/sage/modular>`__,
       `sage.rings.\* <https://github.com/sagemath/sage/tree/develop/src/sage/rings>`__,
       `sage.schemes.\* <https://github.com/sagemath/sage/tree/develop/src/sage/schemes>`__,
       `sage.symbolic.\* <https://github.com/sagemath/sage/tree/develop/src/sage/symbolic>`__,
       `sage.tensor.\* <https://github.com/sagemath/sage/tree/develop/src/sage/tensor>`__,
       `sage.topology.\* <https://github.com/sagemath/sage/tree/develop/src/sage/topology>`__,
       ...
     -
   * - Absorbed SageMath extensions
     - `Sage-Combinat <https://wiki.sagemath.org/combinat>`__,
       `SageManifolds <https://sagemanifolds.obspm.fr/>`__
     -
   * - Downstream software
     - `SageTeX <../reference/spkg/sagetex>`__,
       :trac:`In-distribution Sage user packages <31164>`:
       `admcycles <../reference/spkg/admcycles>`__,
       `sage-flatsurf <../reference/spkg/sage_flatsurf>`__,
       ...
     - `External Sage user packages <https://wiki.sagemath.org/SageMathExternalPackages>`__ :trac:`#31164 <31164>`,
       `GAP-homalg <https://wiki.sagemath.org/Interfaces#Interfaces_to_SageMath_in_other_software>`__
   * - User interfaces
     - Jupyter, JupyterLab :trac:`#30399 <30399>`
     - `Emacs sage-shell-mode <https://wiki.sagemath.org/Emacs>`__,
       `Interfaces <https://wiki.sagemath.org/Interfaces#Interfaces_to_SageMath_in_other_software>`__,
       :trac:`IDEs <30500>` :trac:`#30500 <30500>`
   * - Downstream distributions
     - SageMath distribution :trac:`#33774 <33774>`,
       :trac:`Sage Docker images <wiki/Distribution#Dockerimages>`,
       modularized distributions on PyPI :trac:`#29705 <29705>`,
       wheels :trac:`#31251 <31251>`
     - :trac:`Distributions carrying SageMath <wiki/Distribution#Packagemanagers>`:
       conda-forge, archlinux, debian/ubuntu, fedora,
       `sage-on-gentoo <https://github.com/cschwan/sage-on-gentoo>`__,
       nix, voidlinux, ... :trac:`#33775 <33775>`;
       :trac:`Third-party Docker images <wiki/Distribution#Dockerimages>`
   * - Downstream deployments
     - :trac:`Gitpod <wiki/ReleaseTours/sage-9.6#SagedevelopmentinthecloudwithGitpod>` :trac:`#33113 <33113>`
     - `Running SageMath in the cloud <../installation/>`__:
       `CoCalc <https://cocalc.com/?utm_source=trac.sagemath.org>`__,
       `SageCell <https://sagecell.sagemath.org/>`__,
       `Binder <https://github.com/sagemath/sage-binder-env>`__
   * - Downstream of downstream
     - Devcontainers :trac:`#33671 <33671>`, :trac:`#34363 <34363>`
     -

(In software development, "upstream" refers to a direction toward the original authors or
maintainers of software (:wikipedia:`Upstream_(software_development)`).  The SageMath library
depends on packages such as Singular, so from the viewpoint of the SageMath library, these
packages are "upstream packages".  Distributions such as Arch Linux, Debian, or conda-forge have
a double role: When we build Sage from source, we can make use of packages that are provided by
these distributions, which makes Sage downstream from them.  But the distributions also provide
a SageMath package; in this role, we are upstream.)

But even just within the SageMath library, you could:

* Add examples to the documentation
* Find bugs or typos
* Fix a bug
* Implement a new function
* Contribute a useful tutorial for a mathematical topic
* Translate an existing document to a new language
* Create a new class, create a fast new C library, etc.

This document tells you what you need to know to do all the above, from
reporting bugs to modifying and extending Sage and its documentation.  We also
discuss how to share your new and modified code with other Sage users around the
globe.

Here are brief overviews of each part; for more details, see the extended table
of contents below.  No matter where you start, good luck and welcome to Sage
development!

- **Trac server:** all changes go through `the Sage Trac server
  <https://trac.sagemath.org>`_ at some point. It contains bug reports, upgrade
  requests, changes in progress, and those already part of Sage
  today. :ref:`Click here <chapter-sage-trac>` for more information.

  Importantly, you will need to :ref:`create a trac account
  <section-trac-account>` in order to contribute.

- **Source code:** You need your own copy of Sage's source code to change it.
  `Go there <http://doc.sagemath.org/html/en/installation/source.html>`_ to get it
  and for instructions to build it.

  If you have never worked on software before, pay close attention to the
  `prerequisites to compile
  <http://doc.sagemath.org/html/en/installation/source.html#prerequisites>`_ on your
  system.

  As an easy way to get started, you can run and edit Sage's code and contribute
  your changes using `Gitpod <https://www.gitpod.io>`_,
  a free online development environment based on VS Code.
  It will launch a pre-made workspace with all dependencies and tools installed 
  so that you can start contributing straight away.
  Start by `going to Gitpod <https://gitpod.io/#https://github.com/sagemath/sage>`_,
  and read :ref:`our Gitpod guidelines <section-gitpod>` to learn more.

- **Conventions:** read our :ref:`conventions and guidelines
  <section-writing-code-for-sage>` for code and documentation.

  For everything related to manuals, tutorials, and languages, :ref:`click here
  <chapter-sage_manuals>`.

- **Git (revision control):** To share changes with the Sage community, you will
  need to learn about revision control; we use the software Git for this
  purpose.

  - :ref:`How to install it? <section-git-install>`
  - :ref:`How to configure it for use with Trac? <section-git-setup-name>`
  - :ref:`Here is <chapter-walkthrough>` an overview of our development flow.

Git and Trac for Sage development
=================================

First Steps with Git
--------------------

Sage uses git for version control.

.. toctree::
   :maxdepth: 3

   git_setup
   walk_through

.. _section-git-tricks-and-tips:

Using Git with Trac
-------------------

To contribute back your changes to Sage source code to the project,
you will need a ticket on the
`Sage trac server <http://trac.sagemath.org>`_.

.. toctree::
   :maxdepth: 2

   trac
   manual_git
   git_background
   advanced_git
   workflows
   git_trac


.. _section-writing-code-for-sage:

Writing Code for Sage
=====================

.. toctree::
   :maxdepth: 3

   workspace
   coding_basics
   reviewer_checklist

Running Sage's tests
--------------------

.. toctree::
   :maxdepth: 3

   doctesting

Testing on multiple platforms
-----------------------------

.. toctree::
   :maxdepth: 3

   portability_testing

Additional development and testing tools
----------------------------------------

.. toctree::
   :maxdepth: 3

   tools

Contributing to Manuals and Tutorials
-------------------------------------

.. toctree::
   :maxdepth: 3

   sage_manuals

Sage Coding Details
-------------------

.. toctree::
   :maxdepth: 3

   coding_in_python
   coding_in_cython
   coding_in_other

Packaging the Sage Library
--------------------------

.. toctree::
   :maxdepth: 3

   packaging_sage_library

Packaging Third-Party Code
--------------------------

.. toctree::
   :maxdepth: 3

   packaging


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License <http://creativecommons.org/licenses/by-sa/3.0/>`_.
