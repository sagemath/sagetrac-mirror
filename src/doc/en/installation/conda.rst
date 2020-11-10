.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed from `conda-forge <https://conda-forge.org>`_ on Linux
and macOS running x86_64 that most current desktops and laptops use and also ``aarch64`` on Linux.

To install SageMath, install `Miniforge <https://github.com/conda-forge/miniforge#miniforge3>`_ and
then type in the following commands in a terminal:

* Create a new environment containing SageMath: ``conda create -n sage sage python=X``, where
  ``X`` is version of Python, e.g. ``3.7``
* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

If you are using Miniconda/Anaconda, then do the following before the above steps.
* Add the conda-forge channel: ``conda config --add channels conda-forge``
* Change channel priority to strict: ``conda config --set channel_priority strict``

.. note:: The dependency resolution process of conda is slow; a faster experimental resolver
   we tried with success is `mamba <https://github.com/QuantStack/mamba>`_. It can be used as follows.

   .. code-block:: shell

      conda install mamba         # installs mamba
      mamba create -n sage sage   # replaces "conda create..."


