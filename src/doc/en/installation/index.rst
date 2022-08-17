.. _installation-guide:

=======================================
Welcome to the Sage Installation Guide!
=======================================

If you are reading this manual at https://doc.sagemath.org/, note that
it was built at the time the most recent stable release of SageMath
was made.

More up-to-date information and details regarding supported platforms
may have become available afterwards and can be found in the section
"Availability and installation help" of the
`release tour <https://wiki.sagemath.org/ReleaseTours>`_ for each
SageMath release.

**Where would you like to run SageMath?** Pick one of the following sections.

macOS
=====

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

  - **No development:**

    - Install the `binary build of SageMath <https://github.com/3-manifolds/Sage_macOS/releases>`_
      from the 3-manifolds project.  It is a signed and notarized app, which
      works for macOS 10.12 and newer. It is completely self-contained and
      provides the standard Sage distribution together with many optional
      packages. Additional optional Python packages can be installed with the
      ``%pip`` magic command and will go into your ``~/.sage`` directory.

    - Alternatively, install SageMath from the `conda-forge
      <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

Windows
=======

- **Recommended: Use Visual Studio Code with SageMath devcontainers**

  - Install `VS Code <https://code.visualstudio.com/>`_

  - In VS Code, if the `Visual Studio Code Remote - Containers
    <https://code.visualstudio.com/docs/remote/containers>`_ extension
    is not already installed, click the "Extension" icon on the left
    (or press :kbd:`Ctrl` + :kbd:`Shift` + :kbd:`X`) to open a list of
    extensions. Search for "Remote - Containers" and install it.

    The extension will walk you through the installation of Docker
    Desktop and activating Windows Subsystem for Linux (WSL 2).

  - In VS Code, `clone
    <https://code.visualstudio.com/docs/editor/versioncontrol#_cloning-a-repository>`_
    the `SageMath git repository
    <https://github.com/sagemath/sage>`_. For example, use the VS Code
    Command Palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`), enter the
    command "Git: Clone", hit :kbd:`Enter`, paste the URL https://github.com/sagemath/sage
    and hit :kbd:`Enter` again.

  - In VS Code, find the folder ``.devcontainer``.

  - **Do you want to do SageMath development?**

    - **Yes, development:**

      Copy the file ``devcontainer.json`` from the subfolder ``portability-ubuntu-jammy-standard``
      to the folder ``.devcontainer``.

    - **No development:**

      Copy the file ``devcontainer.json`` from the subfolder ``downstream-archlinux-latest``
      to the folder ``.devcontainer``.

  - In VS Code, use the command palette (:kbd:`Ctrl` + :kbd:`Shift` + :kbd:`P`),
    enter the command "Remote-Containers: Reopen Folder in Container" , and hit :kbd:`Enter`.

    By clicking on "show log", you can see what it does; it will take a while.

  - To use Sage in a terminal, `open a new terminal in VS Code
    <https://code.visualstudio.com/docs/terminal/basics>`_ and type ``./sage`` (development)
    or ``sage`` (no development).

  - To use Sage in a Jupyter notebook, click the "Extension" icon on the left (or press
    :kbd:`Ctrl` + :kbd:`Shift` + :kbd:`X`) to open a list of
    extensions. Install the "Jupyter" extension. In the command
    palette, enter "Create: New Jupyter Notebook", and hit
    :kbd:`Enter`. Click "Select Kernel" on the right (or press :kbd:`Ctrl` +
    :kbd:`Alt` + :kbd:`Enter`), select SageMath, and hit :kbd:`Enter`.

  We recommend this installation method in particular to users who are not already
  familiar with Windows Subsystem for Linux. It is also suitable for machines that
  cannot run WSL 2 (because of hardware constraints or policy set by system
  administrators) but can run Docker.

- **Alternative: Manual installation using Windows Subsystem for Linux**

  Enable Windows Subsystem for Linux (WSL) by following the
  `official WSL setup guide
  <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`_. Be
  sure to do the steps to install WSL2 and set it as default.
  Then go to the Microsoft Store and install Ubuntu (or another
  Linux distribution). Start Ubuntu from the start menu.

  - **Do you want to do SageMath development?**

    - **Yes, development:**

      Then follow the instructions for development on Linux below.

    - **No development:**

      On the Linux running on WSL, you always have root access, so you
      can use any of the installation methods described below for
      Linux.

  - To use Sage in a Jupyter notebook, it is convenient to use Visual Studio Code.
    See :ref:`sec-launching-vscode-wsl`.

Linux
=====

- **Do you want to do SageMath development?**

  - **Yes, development:**

    Obtain the SageMath sources via ``git`` as described in `The Sage
    Developer's Guide
    <https://doc.sagemath.org/html/en/developer/walk_through.html#chapter-walkthrough>`_.

    - Then build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

    - Alternatively, follow the instructions in section
      :ref:`sec-installation-conda-develop`;
      these describe an experimental method that gets all required
      packages, including Python packages, from conda-forge.

  - No development: **Do you have root access (sudo)?**

    - **Yes, root access:** Then the easiest way to install SageMath is
      through a Linux distribution that provides it as a package.  Most
      major Linux distributions have up-to-date versions of SageMath,
      see `repology.org: sagemath
      <https://repology.org/project/sagemath/versions>`_ for an
      overview.  See :ref:`sec-GNU-Linux` for additional information.

      If you are on an older version of your distribution and a recent
      version of SageMath is only available on a newer version of the
      distribution, consider upgrading your distribution.
      In particular, do not install a version of Sage older than 9.2.

    - **No root access, or on an older distribution** Install SageMath from
      the `conda-forge <https://conda-forge.org/>`_ project, as described in section
      :ref:`sec-installation-conda`.

    - Alternatively, build SageMath from source as described in section
      :ref:`sec-installation-from-sources`.

In the cloud
============

- `CoCalc <https://cocalc.com/>`_: an online service that provides SageMath and
  many other tools.

- On any system that allows you to bring your own Docker images to run in
  a container:  Use the `Docker image sagemath/sagemath <https://hub.docker.com/r/sagemath/sagemath/>`_.

- `Sage Cell Server <https://sagecell.sagemath.org/>`_: an online service for
  elementary SageMath computations.


More information
================

.. toctree::
   :maxdepth: 2

   linux
   binary
   conda
   source
   launching
   troubles

This work is licensed under a `Creative Commons Attribution-Share Alike
3.0 License`__.

__ http://creativecommons.org/licenses/by-sa/3.0/
