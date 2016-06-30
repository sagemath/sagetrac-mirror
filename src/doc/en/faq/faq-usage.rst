.. _chapter-faq-usage:

===============
FAQ: Using Sage
===============


How do I get started?
"""""""""""""""""""""

You can try out Sage without downloading anything:

* **SageMathCloud™:** Go to http://cloud.sagemath.org and set up a free
  account.

  If you log in, you will gain access to the latest version of Sage and to 
  many other programs.

* **Sage cell:** A "one-off" version of Sage, available for doing one
  computation at a time. http://sagecell.sagemath.org/

To download a **pre-built binary** Sage distribution, visit
http://sagemath.org/download.html and click on the link for the binary for your
operating system.

The **source code** of Sage is also available for you to download and use. Go to
http://www.sagemath.org/download-source.html to download the tar archive for any
release of Sage.

The Sage notebook runs within a web browser. To start the notebook,
issue the following command in a terminal, if ``sage`` is in your ``PATH`` ::

    sage -notebook

You can also run it from the command line of sage::

    sage: notebook() # not tested

Where can I find more information about using SageMathCloud™?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Visit the Frequently Asked Questions page for SageMathCloud™ at https://github.com/sagemath/cloud/wiki/FAQ .

What are the prerequisites for installing a copy of Sage on my computer?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Most of the dependencies of Sage are shipped with Sage itself. In most
cases, you can download a pre-built binary and use that without
installing any dependencies. If you use Windows, you will need to
install
`VirtualBox <http://www.virtualbox.org>`_, which can be downloaded
from the page http://www.virtualbox.org/wiki/Downloads. After
installing VirtualBox, you need to download a VirtualBox distribution
of Sage available at
http://www.sagemath.org/download-windows.html. Ensure you follow the
instructions at that page. Now you can start the Sage virtual machine
using the VirtualBox software, wait for the virtual machine to boot
up, then type ``notebook`` at the prompt.

You can get the complete source for Sage to compile it on your own
Linux or Mac OS X system. Sage lives in an isolated directory and does
not interfere with your surrounding system. It ships together with
everything necessary to develop Sage, the source code, all its
dependencies and the complete changelog. On Linux systems like
Debian/Ubuntu, you may have to install the ``build essential``
package and the ``m4`` macro processor. Your system
needs to have a working C compiler if you want to compile Sage
from source. On
Debian/Ubuntu, you can install these prerequisites as follows::

    sudo apt-get install build-essential m4

If you have a multi-core system, you can opt for a parallel build of
Sage. The command ::

    export MAKE='make -j8'

will enable 8 threads for parts of the build that support
parallelism. Change the number 8 as appropriate to suit the number of
cores on your system.


How to get Sage's Python to recognize my system's Tcl/Tk install?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It may be that you have Tcl/Tk installed and that your system's Python
recognizes it but Sage's Python does not. To fix that, install the
tcl/tk development library. On Ubuntu, this is the command ::

    sudo apt-get install tk8.5-dev

or something along that line. Next, reinstall Sage's Python::

    sage -f python2

This will pick up the tcl/tk library automatically. After successfully
reinstalling Sage's Python, from within the Sage command line interface,
issue these commands::

    import _tkinter
    import Tkinter

If they do not raise an ``ImportError`` then it worked.


How do I import Sage into a Python script?
""""""""""""""""""""""""""""""""""""""""""

You can import Sage as a library in a Python script. One caveat is
that you need to run that Python script using the version of Python
that is bundled with Sage; currently Python 2.6.x. To import Sage, put
the following in your Python script::

    from sage.all import *

When you want to run your script, you need to invoke Sage with the
option ``-python`` which would run your script using the version of
Python that comes with Sage. For example, if Sage is in your ``PATH``
variable then you can do this::

    sage -python /path/to/my/script.py

Another way is to write a Sage script and run that script using Sage
itself. A Sage script has the file extension ``.sage`` and is more or
less a Python script but uses Sage-specific functions and
commands. You can then run that Sage script like so::

    sage /path/to/my/script.sage

This will take care of loading the necessary environment variables and
default imports for you.

How can I reload a Python script in a Sage session?
"""""""""""""""""""""""""""""""""""""""""""""""""""

You can load a Python script in a Sage session with the command **load**. For example, we could use Sage to import a file called simple.py with::

    load("simple.py")

and repeat this command every time that we change the file simple.py. However, if we type::

    attach("simple.py")

every change applied to the file simple.py will be automatically updated in Sage.

Can I use Sage with Python 3.x?
"""""""""""""""""""""""""""""""

Currently, no. Sage depends on the
`SciPy <http://www.scipy.org>`_
stack of Python numerical and scientific packages. As of 2010, SciPy
still uses Python 2.x. So until SciPy is ported to run with Python
3.x and
`Cython <http://www.cython.org>`_
supports Python 3.x, Sage will continue to use Python 2.x.


I'm seeing an error about "Permission denied" on a file called "sage-flags.txt".
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When sage is built from source, it keeps track of what special
instructions your CPU supports (such as SSE2) and records these. This
is so that if you try running the code on a different machine, which
does not support these extra instructions, you get a sensible error
message instead of a segfault or illegal instruction. Since this
should be stored with Sage itself (as opposed to a user's ``.sage``
directory), it has to be created by someone with the appropriate
permissions. So if you are seeing something like this ::

    Traceback (most recent call last):
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 174, in <module>
        t, R = install_moved()
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 18, in install_moved
        write_flags_file()
      File "/usr/local/sage-4.0.2/local/bin/sage-location", line 82, in write_flags_file
        open(flags_file,'w').write(get_flags_info())
    IOError: [Errno 13] Permission denied:
      '/usr/local/sage-4.0.2/local/lib/sage-flags.txt'

it probably means that you compiled/installed Sage as one user, but
have not run it to let it generate the ``sage-flags.txt`` file. Just
run Sage one time as whatever user installed it and this problem
should go away. This would also be easy to fix by having Sage run once
as part of the install process; see
`trac ticket #6375 <http://trac.sagemath.org/sage_trac/ticket/6375>`_
for this fix.


I downloaded a Sage binary and it crashes on startup with "Illegal instruction". What can I do?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

One way to fix this is to build Sage entirely from source. Another
option is to fix your Sage installation by rebuilding MPIR and ATLAS
by typing the following from the ``SAGE_ROOT`` of your Sage
installation directory and wait about 15 to 20 minutes::

    rm spkg/installed/mpir* spkg/installed/atlas*
    make

It is possible that the binaries have been built for a newer
architecture than what you have. Nobody has yet figured out how to
build Sage in such a way that MPIR and ATLAS work on all
hardware. This will eventually get fixed. Any help is appreciated.


I used Debian/Ubuntu to install Sage 3.0.5 and that version is giving lots of errors. What can I do?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The version of Sage, i.e. Sage version 3.0.5, that is available
through ``apt-get`` in Debian and Ubuntu is very old. No one has yet
found time to update the Debian/Ubuntu version of Sage. Any help is
greatly appreciated. You should download the latest version of Sage
from the
`download page <http://www.sagemath.org/download.html>`_.
If you would like to help with updating the Debian/Ubuntu version of
Sage, please email the
`sage-devel <http://groups.google.com/group/sage-devel>`_
mailing list.


Should I use the official version or development version?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You are encouraged to use the latest official version of
Sage. Development versions are frequently announced on the
`sage-devel <http://groups.google.com/group/sage-devel>`_
and
`sage-release <http://groups.google.com/group/sage-release>`_
mailing lists. An easy way of helping out with Sage development is to
download the latest development release, compile it on your system,
run all doctests, and report any compilation errors or doctest
failures.


Is Sage difficult to learn?
"""""""""""""""""""""""""""

Basic features of Sage should be as easy to learn as learning the
basics of Python. Numerous tutorials are available online to help you
learn Sage. To get the most out of Sage, you are encouraged to learn
some features of the Python programming language. Here is an
incomplete list of resources on Python. Further resources can be found
by a web search.

* `Building Skills in Python <http://homepage.mac.com/s_lott/books/python.html>`_
  by Steven F. Lott
* `Dive into Python <http://www.diveintopython.net>`_ by Mark Pilgrim
* `How to Think Like a Computer Scientist <http://www.openbookproject.net/thinkCSpy>`_
  by Jeffrey Elkner, Allen B. Downey, and Chris Meyers
* `Official Python Tutorial <http://docs.python.org/tutorial>`_
* `Python <http://www.python.org>`_ home page and the
  `Python standard documentation <http://docs.python.org>`_


Can I do X in Sage?
"""""""""""""""""""

You are encouraged to use Sage's tab autocompletion. Just type a few
characters, hit the tab key, and see if the command you want appears
in the list of tab autocompletion. If you have a command called
``mycmd``, then type ``mycmd.`` and hit the tab key to get a list of
functionalities that are supported by that command. To read the
documentation of ``mycmd``, type ``mycmd?`` and press the enter key to
read the documentation for that command. Similarly, type ``mycmd??``
and hit the enter key to get the source code of that command. You are
also encouraged to search through the source code and documentation of
the Sage library. To search through the source code of the Sage
library, use the command ``search_src("<search-keyword>")`` where you
should replace ``<search-keyword>`` with the key words you are looking
for. Also, you can search through the documentation of the Sage
library using the command ``search_doc("<search-keyword>")``.


What exactly does Sage do when I type "0.6**2"?
"""""""""""""""""""""""""""""""""""""""""""""""

When you type "0.6**2" in Python, it returns something like
0.35999999999999999. But when you do the same in Sage it returns
0.360000000000000. To understand why Python behaves as it does, see
the
`Python Tutorial <http://docs.python.org/tutorial/floatingpoint.html>`_,
especially the chapter "Floating Point Arithmetic: Issues and
Limitations". What Sage does is "preparse" the input and transforms it
like this::

    sage: preparse("0.6**2")
    "RealNumber('0.6')**Integer(2)"

So what is *actually* run is::

    RealNumber('0.6')**Integer(2)

The Sage developers (in fact, Carl Witty) decided that Sage floating
point numbers should by default print only the known correct decimal
digits, when possible, thus skirting the problem that Python has. This
decision has its pros and cons. Note that ``RealNumber`` and
``Integer``  are Sage specific, so you would not be able to just type
the above into Python and expect it to work without first an import
statement such as::

    from sage.all import RealNumber, Integer, preparse


Why is Sage's command history different from Magma's?
"""""""""""""""""""""""""""""""""""""""""""""""""""""

Using Sage, you are missing a feature of the Magma command line
interface. In Magma, if you enter a line found in history using up
arrow key and then press down arrow key, then the next line in history
is fetched. This feature allows you to fetch as many successive lines
in history as you like. However, Sage does not have a similar
feature. The
`IPython <http://ipython.scipy.org>`_
command prompt uses the readline library (via pyreadline), which
evidently does not support this feature. Magma has its own custom
"readline-like" library, which does support this feature. (Since so
many people have requested this feature, if anybody can figure out how
to implement it, then such an implementation would certainly be
welcome!)


I have type issues using SciPy, cvxopt or NumPy from Sage.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You are using SciPy or cvxopt or NumPy from Sage and you get type
errors, e.g. ::

    TypeError: function not supported for these types, and can't coerce safely to supported types.

When you type in numbers into Sage, the pre-processor converts them to
a base ring, which you can see by doing::

    sage: preparse("stats.uniform(0,15).ppf([0.5,0.7])")
    "stats.uniform(Integer(0),Integer(15)).ppf([RealNumber('0.5'),RealNumber('0.7')])"

Unfortunately, NumPy support of these advanced Sage types like
``Integer`` or ``RealNumber`` is not yet at 100%. As a solution,
redefine ``RealNumber`` and/or ``Integer`` to change the behavior of
the Sage preparser, so decimal literals are floats instead of Sage
arbitrary precision real numbers, and integer literals are Python
ints. For example::

    sage: RealNumber = float; Integer = int
    sage: from scipy import stats
    sage: stats.ttest_ind(list([1,2,3,4,5]),list([2,3,4,5,.6]))
    Ttest_indResult(statistic=0.076752955645333687, pvalue=0.94070490247380478)
    sage: stats.uniform(0,15).ppf([0.5,0.7])
    array([  7.5,  10.5])

Alternatively, be explicit about data types, e.g. ::

    sage: from scipy import stats
    sage: stats.uniform(int(0),int(15)).ppf([float(0.5),float(0.7)])
    array([  7.5,  10.5])

As a third alternative, use the raw suffix::

    sage: from scipy import stats
    sage: stats.uniform(0r,15r).ppf([0.5r,0.7r])
    array([  7.5,  10.5])

You can also disable the preprocessor in your code via
``preparse(False)``. You can may start IPython alone from the command
line ``sage -ipython`` which does not pre-load anything
Sage-specific. Or switching the Notebook language to "Python".


How do I save an object so I don't have to compute it each time I open a worksheet?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The ``save`` and ``load`` commands will save and load an object,
respectively. In the notebook, the ``DATA`` variable is the location
of the data storage area of the worksheet. To save the object
``my_stuff`` in a worksheet, you could do ::

    save(my_stuff, DATA + "my_stuff")

and to reload it, you would just do ::

    my_stuff = load(DATA + "my_stuff")


Does Sage contain a function similar to Mathematica's ToCharacterCode[]?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You might want to convert ASCII characters such as "Big Mac" to ASCII
numerals for further processing. In Sage and Python, you can use ``ord``,
e.g. ::

    sage: map(ord, "abcde")
    [97, 98, 99, 100, 101]
    sage: map(ord, "Big Mac")
    [66, 105, 103, 32, 77, 97, 99]


Can I make Sage automatically execute commands on startup?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Yes, just make a file ``$HOME/.sage/init.sage`` and it will be
executed any time you start Sage. This assumes that the Sage
environment variable ``DOT_SAGE`` points to the hidden directory
``$HOME/.sage``, which by default is the case.


My Sage upgrade failed with missing gmp symbols on OSX 10.4. What can I do?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Moving a Sage install on Mac OS X 10.4 and then upgrading anything
that is linked against NTL leads to link errors due to missing gmp
symbols. The problem is the link mode with which the dynamic NTL is
created. There is have a fix, but it still being verified that it
really fixes the issue. Everything that is linked against NTL needs to
be recompiled, i.e. singular and cremona at the moment. To add to the
confusion: This is not an issue on Mac OS X 10.5. A fix for this issue
went into Sage 2.8.15, so please report if you see this with a more
current Sage release.


When I compile Sage my computer beeps and shuts down or hangs.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Compiling Sage is quite taxing on the CPU. The above behavior usually
indicates that your computer has overheated. In many cases this can be
fixed by cleaning the CPU fan and assuring proper ventilation of the
system. Please ask your system administrator or a professional to do
this in case you have never done this. Such hardware maintenance, if
not performed by a skilled professional, you can potentially damage
your system.

For Linux users, if you suspect that the compilation fails because of
a resource issue, a fix might be to edit your ``/etc/inittab`` so that
Linux boots into run level 3. The file ``/etc/inittab`` usually
contains something similar to the following snippet::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    id:5:initdefault:

which directs your Linux distribution to boot into a graphical login
screen. Comment out the line ``id:5:initdefault:`` and add the line
``id:3:initdefault:``, so that you now have something like::

    #   0 - halt (Do NOT set initdefault to this)
    #   1 - Single user mode
    #   2 - Multiuser, without NFS (The same as 3, if you do not have
    #   networking)
    #   3 - Full multiuser mode
    #   4 - unused
    #   5 - X11
    #   6 - reboot (Do NOT set initdefault to this)
    #
    # id:5:initdefault:
    id:3:initdefault:

Now if you reboot your system, you will be greeted with a text based
login screen. This allows you to log into your system with a text
based session from within a virtual terminal. A text based session
usually does not consume as much system resources as would be the case
with a graphical session. Then build your Sage source distribution
from within your text based session. You need to make sure that you
can first restore your graphical session, before you attempt to log
into a text based session.


When I run doctests on Mac OS X I see the messages with "malloc", but in the end Sage reports that everything went fine.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The "malloc" messages you refer to might be something such as the
following::

    sage -t  src/sage/libs/pari/gen.pyx
    python(4563) malloc: *** vm_allocate(size=4096000000) failed (error code=3)
    python(4563) malloc: *** error: can't allocate region
    python(4563) malloc: *** set a breakpoint in szone_error to debug

The issue above is not a doctest failure. It is an error message
printed by the system and it is exactly what one expects to see. In
that particular doctest, we try to allocate a very large list in Pari
that does not fit into physical memory (it is at least 100GB in
size). So Mac OS X tells you that it could not allocate a chunk of
memory roughly 4 GB in size, which is expected, if you are using Sage
on a 32-bit version of OS X and you have compiled Sage in 32-bit bit
mode or your binary Sage distribution is 32-bit.


Sage 2.9 and higher fails compiling ATLAS on Linux. How can I fix this?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The most likely cause is enabled power management. Disabling it should
fix the problem. Depending on your flavor of distribution, this might
either be possible with some nice GUI tool or not. On the command line
do the following as root for each CPU you have::

    /usr/bin/cpufreq-selector -g performance -c #number CPU

On Ubuntu, try disabling "Power Manager" via ::

    System --> Preferences --> Sessions

under the "Startup Programs" or using ``cpufreq-set`` via the command
line.

When I start Sage, SELinux complains that "/path/to/libpari-gmp.so.2" requires text-relocation. How can I fix it?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The problem can be fixed by running the following command::

    chcon -t textrel_shlib_t /path/to/libpari-gmp.so.2


Upgrading Sage went fine, but now the banner still shows the old version. How can I fix this?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The banner is stored and not computed at every new start of Sage. If
it has not been updated, this should not prevent Sage to run
correctly. Type ``banner()`` in a Sage session to check the real
version. If you want the correct banner, you need to build Sage again
by typing ``make build`` in a terminal.


How do I run sage in daemon mode, i.e. as a service?
""""""""""""""""""""""""""""""""""""""""""""""""""""

There are several possibilities. Use ``screen``, ``nohup`` or ``disown``.


The show command for plotting 3-D objects does not work.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The default live 3-D plotting for Sage 6.4+ uses
`Jmol/JSmol <http://jmol.sourceforge.net>`_
for viewing. From the command line the Jmol Java application is used,
and for in browser viewing either pure javascript or a Java applet
is used.  By default in browsers pure javascript is used to avoid
the problems with some browsers that do not support java applet
plugins (namely Chrome).  On each browser worksheet there is a
checkbox which must be checked before a 3-D plot is generated if
the user wants to use the Java applet (the applet is a little faster
with complex plots).

The most likely reason for a malfunction is that you do not have
a Java Run Time Environment (JRE) installed or you have one older than
version 1.7.  If things work from the command line another possibility
is that your browser does not have the proper plugin to support Java
applets (at present, 2014, plugins do not work with most versions of
Chrome).  Make sure you have installed either the IcedTea browser
plugin (for linux see your package manager), see:
`IcedTea <http://icedtea.classpath.org/wiki/IcedTea-Web>`_,
or the Oracle Java plugin see:
`Java <https://java.com/en/download/help/index_installing.xml>`_.

If you are using a Sage server over the web and even javascript rendering
does not work, you may have a problem with your browser's javascript
engine or have it turned off.

May I use Sage tools in a commercial environment?
"""""""""""""""""""""""""""""""""""""""""""""""""

Yes! Absolutely! Basically the *only* constraint is that if you make
changes to Sage itself and redistribute this changed version of Sage
publicly, then you must make these changes available to us so that we
can put them into the standard version of Sage (if we
want). Otherwise, you are free to use as many copies of Sage as you
want completely for free to make money, etc. without paying any
license fees at all.


I want to write some Cython code that uses finite field arithmetic but "cimport sage.rings.finite_field_givaro" fails. What can I do?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You need to give hints to Sage so that it uses C++ (both Givaro and
NTL are C++ libraries), and it also needs the GMP and STDC C++
libraries. Here is a small example::

    # These comments are hints to Sage/Pyrex about the compiler and
    # libraries needed for the Givaro library:
    #
    #clang c++
    #clib givaro gmpxx gmp m stdc++
    cimport sage.rings.finite_field_givaro
    # Construct a finite field of order 11.
    cdef sage.rings.finite_field_givaro.FiniteField_givaro K
    K = sage.rings.finite_field_givaro.FiniteField_givaro(11)
    print "K is a", type(K)
    print "K cardinality =", K.cardinality()
    # Construct two values in the field:
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement x
    cdef sage.rings.finite_field_givaro.FiniteField_givaroElement y
    x = K(3)
    y = K(6)
    print "x is a", type(x)
    print "x =", x
    print "y =", y
    print "x has multiplicative order =", x.multiplicative_order()
    print "y has multiplicative order =", y.multiplicative_order()
    print "x*y =", x*y
    # Show that x behaves like a finite field element:
    for i in range(1, x.multiplicative_order() + 1):
        print i, x**i
    assert x*(1/x) == K.one()

To find out more, type ::

    sage.rings.finite_field_givaro.FiniteField_givaro.

at the Sage prompt and hit tab, then use ``??`` to get more
information on each function. For example::

    sage.rings.finite_field_givaro.FiniteField_givaro.one??

tells you more about the multiplicative unit element in the finite
field.


I'm getting weird build failures on Mac OS X. How do I fix this?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Search the build log (install.log) to see if you are getting the
following log message::

    fork: Resource temporarily unavailable.

If so, try the following. Create (or edit) ``/etc/launchd.conf`` and
include the following::

    limit maxproc 512 2048

then reboot. See
`this page <http://www.macosxhints.com/article.php?story=20050709233920660>`_
for more details.

How do I plot the cube root (or other odd roots) for negative input?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This is one of the most frequently asked questions.  There are several
methods mentioned in the plot documentation, but this one is easiest::

    sage: plot(sign(x)*abs(x)^(1/3),-1,1)
    Graphics object consisting of 1 graphics primitive

The *reason* this is necessary is that Sage returns complex numbers
for odd roots of negative numbers when numerically approximated, which
is a `standard convention <http://en.wikipedia.org/wiki/Cube_root#Complex_numbers>`_.

    sage: N((-1)^(1/3))
    0.500000000000000 + 0.866025403784439*I

How do I use the bitwise XOR operator in Sage?
""""""""""""""""""""""""""""""""""""""""""""""

The exclusive or operator in Sage is ``^^``. This also works for
the inplace operator ``^^=``::

   sage: 3^^2
   1
   sage: a = 2
   sage: a ^^= 8
   sage: a
   10

If you define two variables and then evaluate as follows::

    sage: a = 5; b = 8
    sage: a.__xor__(b), 13
    (13, 13)

You can also do ::

    sage: (5).__xor__(8)
    13

The parentheses are necessary so that Sage does not think you have a
real number. There are several ways to define a function::

    sage: xor = lambda x, y: x.__xor__(y)
    sage: xor(3, 8)
    11

Another option, which sneaks around the Sage
preparser, is ::

    sage: def xor(a, b):
    ....:     return eval("%s^%s" % (a, b))
    sage: xor(3, 8)
    11

You can also turn off the Sage preparser with ``preparser(False)``,
then ``^`` will work just like in Python. You can later turn on the
preparser with ``preparser(True)``. That only works in command line
Sage. In a notebook, switch to Python mode.


When I try to use LaTeX in the notebook, it says it cannot find fullpage.sty.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The general---but perhaps not very helpful---answer is that you need
to install ``fullpage.sty`` into a directory searched by TeX. On
Ubuntu (and probably many other Linux distributions), you should
install the ``texlive-latex-extra`` package. If that is not available,
try installing the ``tetex-extra package``. If you are using Mac OS X,
you will have to use whatever TeX distribution you use to get
``fullpage.sty`` (if you use MacTeX, it is likely already
installed). If you are using the VirtualBox image on Windows, you will
need to log into the VirtualBox image and install
``texlive-latex-extra`` there.


With objects a and b and a function f, I accidentally typed f(a) = b instead of f(a) == b. This returned a TypeError (as expected), but also deleted the object a. Why?
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

It is because of how functions are defined in Sage with the
``f(x) = expr`` notation using the preparser. Also notice that if you
make this mistake inside of an ``if`` statement, you will get a
``SyntaxError`` before anything else goes wrong. So in this case,
there is no problem.


How do I use a different browser with the Sage notebook?
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You will need to do this from the command line.  Just run a command like this.

* Linux (assuming you have Sage in ``/usr/bin``)::

    env SAGE_BROWSER=opera /usr/bin/sage -notebook

* Mac (assuming you are in the directory of your downloaded Sage)::

    SAGE_BROWSER='open -a Firefox' ./sage -notebook
    SAGE_BROWSER='open -a Google\ Chrome' ./sage -notebook
