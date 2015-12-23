"""
Test Binary Files
"""

import os
import subprocess
from glob import glob

from sage.env import SAGE_LOCAL
from sage.tests.local.test_base import (
    LocalFileTest,
    FileMagicException, TestException
)



class ELFBinaryFile(LocalFileTest):

    SAGE_LIBRARIES = tuple(
        os.path.basename(lib).split('.')[0] for lib in
        glob(os.path.join(SAGE_LOCAL, 'lib/lib*.so')) + 
        glob(os.path.join(SAGE_LOCAL, 'lib64/lib*.so'))
    )

    def __init__(self, filename, head):
        r"""
        ELF Binary File

        This is a wrapper around a file whose magic matches that of an
        ELF binary.

        INPUT:

        - ``filename`` -- string. The filename.

        EXAMPLES::

            sage: from sage.tests.local.binary_audit import ELFBinaryFile
            sage: patch = ELFBinaryFile.load(os.path.join(SAGE_LOCAL, 'bin', 'patch'))
            sage: patch
            <sage.tests.local.binary_audit.ELFBinaryFile object at 0x...>
            sage: patch.is_dynamic
            True
            sage: patch.rpath_string
            '/.../local/lib:/.../local/lib64'
            sage: patch.needed   # random output
            ['libattr.so.1', 'libc.so.6']
        """
        if not head.startswith('\x7fELF'):
            raise FileMagicException('not an ELF binary')
        super(ELFBinaryFile, self).__init__(filename, head)
        self.needed = []
        self.rpath_string = None
        self.is_dynamic = False
        objdump = subprocess.check_output(['objdump', '-p', self.filename])
        for line in objdump.splitlines():
            if line.startswith('Dynamic Section:'):
                self.is_dynamic = True
            if line.startswith('  RPATH'):
                self.rpath_string = line[7:].strip()
            if line.startswith('  NEEDED'):
                self.needed.append(line[8:].strip())
        return None

    @classmethod
    def is_sage_library(cls, lib):
        r"""
        Test whether lib is a library shipped with sage

        INPUT:

        - ``lib`` -- string. The library name from the ELF needed
          section, for example ``'libmpir.so.16'``.

        OUTPUT:

        Boolean. Whether the library is supposed to be part of Sage.

        EXAMPLES::

            sage: from sage.tests.local.binary_audit import ELFBinaryFile
            sage: ELFBinaryFile.is_sage_library('libmpir.so.2.3.4')
            True
            sage: ELFBinaryFile.is_sage_library('libaudit.so')
            False
        """
        return any(lib.startswith(name) for name in cls.SAGE_LIBRARIES)
        
    def perform_tests(self):
        r"""
        Perform tests

        RAISES:

        :class:`TestException` if a test fails.

        EXAMPLES::

            sage: from sage.tests.local.binary_audit import ELFBinaryFile
            sage: patch = ELFBinaryFile.load(os.path.join(SAGE_LOCAL, 'bin', 'patch'))
            sage: patch.perform_tests()
        """
        if not self.is_dynamic:
            return   # only dynamic binaries need to have their rpath set
        if not any(self.is_sage_library(lib) for lib in self.needed):
            return   # No dependencies that would require rpath
        rpath = self.rpath()
        local_lib = os.path.abspath(os.path.join(SAGE_LOCAL, 'lib'))
        if local_lib not in rpath:
            raise TestException('{0} not in rpath: {1}'.format(local_lib, rpath))
        local_lib64 = os.path.abspath(os.path.join(SAGE_LOCAL, 'lib64'))
        if local_lib64 not in rpath:
            raise TestException('{0} not in rpath: {1}'.format(local_lib64, rpath))
    
    def rpath(self):
        r"""
        Return the rpath components

        OUTPUT:

        Tuple of strings.

        EXAMPLES::

            sage: from sage.tests.local.binary_audit import ELFBinaryFile
            sage: patch = ELFBinaryFile.load(os.path.join(SAGE_LOCAL, 'bin', 'patch'))
            sage: patch.rpath()
            ('/.../sage/local/lib', '/.../sage/local/lib64')
        """
        if self.rpath_string is None:
            return tuple()
        else:
            return tuple(
                os.path.abspath(component)
                for component in self.rpath_string.split(':')
            )


