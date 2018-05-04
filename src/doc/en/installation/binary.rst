.. _sec-installation-from-binaries:

Install from Pre-built Binaries
===============================

Linux and OS X
--------------

Installation from a pre-built binary tarball should in the long run
be the easiest and fastest way to install Sage. This is not
necessarily the case right now. Note that Sage is itself a
programming environment, so building it from source guarantees you
maximum flexibility in the long run. Nonetheless, we provide
pre-built binaries.

Assumptions: You have a computer with at least 2 GB of free
disk space and the operating system is Linux (32-bit or 64-bit) or
OS X (10.4 or later). An SSL library (OpenSSL recommended) is
installed in your operating system.

Highly Recommended: It is highly recommended that you have LaTeX
installed.  If you want to view animations, you should install either
ImageMagick or ffmpeg.  ImageMagick or dvipng is also used for
displaying some LaTeX output in the Sage notebook.

Download the latest binary tarball from
http://www.sagemath.org/download.html. For example, it might be
called ``sage-x.y.z-x86_64-Linux.tgz``. Unpack it on your computer
in a directory which you have permission to read and write::

    tar zxvf sage-x.y.z-x86_64-Linux.tgz

You can move the resulting directory ``SageMath``
anywhere and still run ``./sage`` from it, as long as the full path name
has **no spaces** in it. The first time you run Sage, you will see a
message like

::

   Rewriting paths for your new installation directory
   ===================================================

   This might take a few minutes but only has to be done once.

   patching ...  (long list of files)

At this point, you can no longer move your Sage installation and
expect Sage to function.

You can also copy the file ``sage`` from
the ``SageMath`` directory and put it anywhere, e.g., ``/usr/local/bin/``, but
then you have to edit the ``#SAGE_ROOT=/path/to/sage-version`` line
at the top of the copied file ``/usr/local/bin/sage``
(you should not edit the original ``sage`` executable).
The variable ``SAGE_ROOT`` should point to the directory
``SageMath`` of the extracted Sage tarball.
As long as ``/usr/local/bin`` is in your
``$PATH``, you can then type ``sage`` from the command line to run
Sage.  Another approach is to create a symbolic link, say
``/usr/local/bin/sage``, pointing to ``SageMath/sage`` ::

    ln -s /path/to/SageMath/sage /usr/local/bin/sage

With this approach, there is no need to edit ``/usr/local/bin/sage``,
the ``SAGE_ROOT`` path will be discovered automatically thanks to the
symbolic link.
When you want to install a new version of Sage, just delete the old
link and create a new one.

We currently distribute ``.dmg`` files for recent versions of OS X, but
we would like to make Sage more of a native application. Work for that
is ongoing, and help is always welcome.


Microsoft Windows
-----------------

The best way to install Sage on Windows is to install
`VirtualBox for Windows <http://www.virtualbox.org/wiki/Downloads>`_
and then download and install the VirtualBox distribution of Sage. See
`this URL <http://www.sagemath.org/download-windows.html>`_ for
further instructions on installing Sage on Windows. For more information
on how to set up and run the Sage Virtual Machine, please read these
`additional instructions <http://wiki.sagemath.org/SageAppliance>`_.
