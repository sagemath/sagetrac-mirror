########################################################################
##
## Customize the Extensions processed by Cython
##
########################################################################

import os
import sys
import time
import json
from distutils import log
from distutils.cmd import Command
from distutils.errors import (DistutilsModuleError,
                              DistutilsOptionError)

from Cython.Build.Dependencies import default_create_extension
from sage_setup.util import stable_uniq, have_module
from sage_setup.find import find_extra_files
from sage_setup.library_order import library_order

from sage.env import (cython_aliases, sage_include_directories)

# Do not put all, but only the most common libraries and their headers
# (that are likely to change on an upgrade) here:
# [At least at the moment. Make sure the headers aren't copied with "-p",
# or explicitly touch them in the respective spkg's spkg-install.]
lib_headers = dict()

# Set by build/bin/sage-build-env-config. Empty if the system package is used.
gmp_prefix = os.environ.get("SAGE_GMP_PREFIX", "")
if gmp_prefix:
    lib_headers["gmp"] = [os.path.join(gmp_prefix, 'include', 'gmp.h')]  # cf. #8664, #9896
    lib_headers["gmpxx"] = [os.path.join(gmp_prefix, 'include', 'gmpxx.h')]
ntl_prefix = os.environ.get("SAGE_NTL_PREFIX", "")
if ntl_prefix:
    lib_headers["ntl"] = [os.path.join(ntl_prefix, 'include', 'NTL', 'config.h')]

# Manually add -fno-strict-aliasing, which is needed to compile Cython
# and disappears from the default flags if the user has set CFLAGS.
#
# Add -DCYTHON_CLINE_IN_TRACEBACK=1 which causes the .c line number to
# always appear in exception tracebacks (by default, this is a runtime
# setting in Cython which causes some overhead every time an exception
# is raised).
extra_compile_args = ["-fno-strict-aliasing", "-DCYTHON_CLINE_IN_TRACEBACK=1"]
extra_link_args = [ ]

DEVEL = False
if DEVEL:
    extra_compile_args.append('-ggdb')

import subprocess
# Work around GCC-4.8 bug which miscompiles some sig_on() statements:
# * http://trac.sagemath.org/sage_trac/ticket/14460
# * http://trac.sagemath.org/sage_trac/ticket/20226
# * http://gcc.gnu.org/bugzilla/show_bug.cgi?id=56982
if subprocess.call("""$CC --version | grep -i 'gcc.* 4[.]8' >/dev/null """, shell=True) == 0:
    extra_compile_args.append('-fno-tree-copyrename')

# Search for dependencies in the source tree and add to the list of include directories
include_dirs = sage_include_directories(use_sources=True)

# Look for libraries only in what is configured already through distutils
# and environment variables
library_dirs = []

class sage_build_cython(Command):
    name = 'build_cython'
    description = "compile Cython extensions into C/C++ extensions"

    user_options = [
        # TODO: Temporarily disabled since the value for this option is
        # hard-coded; change as part of work on #21525
        #('build-dir=', 'd',
        # "directory for compiled C/C++ sources and header files"),
        ('profile', 'p',
         "enable Cython profiling support"),
        ('parallel=', 'j',
         "run cythonize in parallel with N processes"),
        ('force=', 'f',
         "force files to be cythonized even if the are not changed")
    ]

    boolean_options = ['debug', 'profile', 'force']

    def initialize_options(self):
        self.extensions = None
        self.build_base = None
        self.build_dir = None

        # Always have Cython produce debugging info by default, unless
        # SAGE_DEBUG=no explicitly
        self.debug = True
        self.profile = None
        self.parallel = None
        self.force = None

        self.cython_directives = None
        self.compile_time_env = None

        self.build_lib = None
        self.cythonized_files = None

    def finalize_options(self):
        self.extensions = self.distribution.ext_modules

        # Let Cython generate its files in the "cythonized"
        # subdirectory of the build_base directory.
        self.set_undefined_options('build', ('build_base', 'build_base'))
        self.build_dir = os.path.join(self.build_base, "cythonized")

        # Inherit some options from the 'build_ext' command if possible
        # (this in turn implies inheritance from the 'build' command)
        inherit_opts = [('build_lib', 'build_lib'),
                        ('debug', 'debug'),
                        ('force', 'force')]

        # Python 3.5 now has a parallel option as well
        if sys.version_info[:2] >= (3, 5):
            inherit_opts.append(('parallel', 'parallel'))

        self.set_undefined_options('build_ext', *inherit_opts)

        # Always produce debugging output unless SAGE_DEBUG=no is given
        # explicitly
        self.debug = os.environ.get('SAGE_DEBUG', None) != 'no'

        if self.debug:
            log.info('Enabling Cython debugging support')

        if self.profile is None:
            self.profile = os.environ.get('SAGE_PROFILE') == 'yes'

        if self.profile:
            log.info('Enabling Cython profiling support')

        if self.parallel is None:
            self.parallel = os.environ.get('SAGE_NUM_THREADS', '0')

        try:
            self.parallel = int(self.parallel)
        except ValueError:
            raise DistutilsOptionError("parallel should be an integer")

        try:
            import Cython
        except ImportError:
            raise DistutilsModuleError(
                "Cython must be installed and importable in order to run "
                "the cythonize command")

        # Cython compiler directives
        self.cython_directives = dict(
            auto_pickle=False,
            autotestdict=False,
            cdivision=True,
            embedsignature=True,
            fast_getattr=True,
            language_level="3str",
            preliminary_late_includes_cy28=True,
            profile=self.profile,
        )
        self.compile_time_env = dict(
            PY_VERSION_HEX=sys.hexversion,
            PY_MAJOR_VERSION=sys.version_info[0],
        )

        # We check the Cython version and some relevant configuration
        # options from the earlier build to see if we need to force a
        # recythonization. If the version or options have changed, we
        # must recythonize all files.
        self._version_file = os.path.join(self.build_dir, '.cython_version')
        self._version_stamp = json.dumps({
            'version': Cython.__version__,
            'debug': self.debug,
            'directives': self.cython_directives,
            'compile_time_env': self.compile_time_env,
        }, sort_keys=True)

        # Read an already written version file if it exists and compare to the
        # current version stamp
        try:
            if open(self._version_file).read() == self._version_stamp:
                force = False
            else:
                # version_file exists but its contents are not what we
                # want => recythonize all Cython code.
                force = True
                # In case this cythonization is interrupted, we end up
                # in an inconsistent state with C code generated by
                # different Cython versions or with different options.
                # To ensure that this inconsistent state will be fixed,
                # we remove the version_file now to force a
                # recythonization the next time we build Sage.
                os.unlink(self._version_file)
        except IOError:
            # Most likely, the version_file does not exist
            # => (re)cythonize all Cython code.
            force = True

        # If the --force flag was given at the command line, always force;
        # otherwise use what we determined from reading the version file
        if self.force is None:
            self.force = force

    def get_cythonized_package_files(self):
        """
        Return a list of files found in the Sage sources and/or Cythonize
        output directory that should be installed with Python packages (a la
        ``package_files``).
        """

        if self.cythonized_files is not None:
            return self.cythonized_files

        self.cythonized_files = list(find_extra_files(
            ".", ["sage"], self.build_dir, []).items())

        return self.cythonized_files

    def run(self):
        """
        Call ``cythonize()`` to replace the ``ext_modules`` with the
        extensions containing Cython-generated C code.
        """
        from Cython.Build import cythonize
        import Cython.Compiler.Options

        Cython.Compiler.Options.embed_pos_in_docstring = True

        log.info("Updating Cython code....")
        t = time.time()
        extensions = cythonize(
            self.extensions,
            nthreads=self.parallel,
            build_dir=self.build_dir,
            force=self.force,
            aliases=cython_aliases(),
            compiler_directives=self.cython_directives,
            compile_time_env=self.compile_time_env,
            create_extension=self.create_extension,
            # Debugging
            gdb_debug=self.debug,
            output_dir=os.path.join(self.build_lib, "sage"),
            # Disable Cython caching, which is currently too broken to
            # use reliably: http://trac.sagemath.org/ticket/17851
            cache=False,
            )

        # Filter out extensions with skip_build=True
        extensions = [ext for ext in extensions if not getattr(ext, "skip_build", False)]

        # We use [:] to change the list in-place because the same list
        # object is pointed to from different places.
        self.extensions[:] = extensions

        log.info("Finished Cythonizing, time: %.2f seconds." % (time.time() - t))

        with open(self._version_file, 'w') as f:
            f.write(self._version_stamp)

        # Finally, copy relevant cythonized files from build/cythonized
        # tree into the build-lib tree
        for (dst_dir, src_files) in self.get_cythonized_package_files():
            dst = os.path.join(self.build_lib, dst_dir)
            self.mkpath(dst)
            for src in src_files:
                self.copy_file(src, dst, preserve_mode=False)

    def create_extension(self, template, kwds):
        """
        Create a distutils Extension given data from Cython.

        This adjust the ``kwds`` in the following ways:

        - Make everything depend on *this* setup.py file

        - Add dependencies on header files for certain libraries

        - Sort the libraries according to the library order

        - Add some default compile/link args and directories

        - Choose C99 standard for C code and C++11 for C++ code

        - Drop -std=c99 and similar from C++ extensions

        - Ensure that each flag, library, ... is listed at most once
        """
        lang = kwds.get('language', 'c')
        cplusplus = (lang == "c++")

        # Libraries: sort them
        libs = kwds.get('libraries', [])
        kwds['libraries'] = sorted(set(libs),
                key=lambda lib: library_order.get(lib, 0))

        # Dependencies: add setup.py and lib_headers
        depends = kwds.get('depends', []) + [self.distribution.script_name]
        for lib, headers in lib_headers.items():
            if lib in libs:
                depends += headers
        kwds['depends'] = depends  # These are sorted and uniq'ed by Cython

        # Process extra_compile_args
        cflags = []
        have_std_flag = False
        for flag in kwds.get('extra_compile_args', []):
            if flag.startswith("-std="):
                if cplusplus and "++" not in flag:
                    continue  # Skip -std=c99 and similar for C++
                have_std_flag = True
            cflags.append(flag)
        if not have_std_flag:  # See Trac #23919
            if sys.platform == 'cygwin':
                # Cygwin (particularly newlib, Cygwin's libc) has some bugs
                # with strict ANSI C/C++ in some headers; using the GNU
                # extensions typically fares better:
                # https://trac.sagemath.org/ticket/24192
                if cplusplus:
                    cflags.append("-std=gnu++11")
                else:
                    cflags.append("-std=gnu99")
            else:
                if cplusplus:
                    cflags.append("-std=c++11")
                else:
                    cflags.append("-std=c99")
        cflags = extra_compile_args + cflags
        kwds['extra_compile_args'] = stable_uniq(cflags)

        # Process extra_link_args
        ldflags = kwds.get('extra_link_args', []) + extra_link_args
        kwds['extra_link_args'] = stable_uniq(ldflags)

        # Process library_dirs
        lib_dirs = kwds.get('library_dirs', []) + library_dirs
        kwds['library_dirs'] = stable_uniq(lib_dirs)

        # Process include_dirs
        inc_dirs = kwds.get('include_dirs', []) + include_dirs + [self.build_dir]
        kwds['include_dirs'] = stable_uniq(inc_dirs)

        return default_create_extension(template, kwds)
