"""
Tar file support
"""

#*****************************************************************************
#       Copyright (C) 2016 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from __future__ import print_function

import os
import copy
import tarfile
import stat
import subprocess
import time

from io import BytesIO

from sage_bootstrap.uncompress.filter_os_files import filter_os_files


class SageBaseTarFile(tarfile.TarFile):
    """
    Sage as tarfile.TarFile, but applies a reasonable umask (0022) to the
    permissions of all extracted files and directories.

    Previously this applied the user's current umask per the default behavior
    of the ``tar`` utility, but this did not provide sufficiently reliable
    behavior in all cases, such as when the user's umask is not strict enough.

    This also sets the modified timestamps on all extracted files to the same
    time, not the timestamps stored in the tarball.

    See http://trac.sagemath.org/ticket/20218#comment:16 and
    https://trac.sagemath.org/ticket/24567 for more background.
    """

    umask = 0o022

    def __init__(self, *args, **kwargs):
        # Unfortunately the only way to get the current umask is to set it
        # and then restore it
        super(SageBaseTarFile, self).__init__(*args, **kwargs)

        # When extracting files, this will be set to the time of the first
        # file extracted, so that all files can be set to the same mtime
        self._extracted_mtime = None

    @property
    def names(self):
        """
        List of filenames in the archive.

        Filters out names of OS-related files that shouldn't be in the
        archive (.DS_Store, etc.)
        """

        return filter_os_files(self.getnames())

    def chmod(self, tarinfo, target):
        """Apply ``self.umask`` instead of the permissions in the TarInfo."""
        tarinfo = copy.copy(tarinfo)
        tarinfo.mode &= ~self.umask
        tarinfo.mode &= ~(stat.S_ISUID | stat.S_ISGID)
        return super(SageBaseTarFile, self).chmod(tarinfo, target)

    def utime(self, tarinfo, target):
        """Override to keep the extraction time as the file's timestamp."""

        tarinfo.mtime = self._extracted_mtime
        return super(SageBaseTarFile, self).utime(tarinfo, target)

    def extractall(self, path='.', members=None):
        """
        Same as tarfile.TarFile.extractall but allows filenames for
        the members argument (like zipfile.ZipFile).
        """

        self._extracted_mtime = None

        if members:
            name_to_member = dict([member.name, member] for member in self.getmembers())
            members = [m if isinstance(m, tarfile.TarInfo)
                       else name_to_member[m]
                       for m in members]
        return super(SageBaseTarFile, self).extractall(path=path,
                                                       members=members)

    def extractbytes(self, member):
        """
        Return the contents of the specified archive member as bytes.

        If the member does not exist, returns None.
        """

        if member in self.getnames():
            reader = self.extractfile(member)
            return reader.read()

    def _extract_member(self, tarinfo, targetpath):
        """
        Override to ensure that our custom umask is applied over the entire
        directory tree, even for directories that are not explicitly listed in
        the tarball.
        """

        if self._extracted_mtime is None:
            self._extracted_mtime = time.time()

        old_umask = os.umask(self.umask)
        try:
            super(SageBaseTarFile, self)._extract_member(tarinfo, targetpath)
        finally:
            os.umask(old_umask)


class SageTarFile(SageBaseTarFile):
    """
    A wrapper around SageBaseTarFile such that SageTarFile(filename) is
    essentially equivalent to TarFile.open(filename) which is more
    flexible than the basic TarFile.__init__
    """
    def __new__(cls, filename):
        return SageBaseTarFile.open(filename)

    @staticmethod
    def can_read(filename):
        """
        Given an archive filename, returns True if this class can read and
        process the archive format of that file.
        """
        return tarfile.is_tarfile(filename)


class SageTarXZFile(SageBaseTarFile):
    """
    A ``.tar.xz`` file which is uncompressed in memory.
    """
    def __new__(cls, filename):
        # Read uncompressed data through a pipe
        proc = subprocess.Popen(["xz", "-d", "-c", filename], stdout=subprocess.PIPE)
        data, _ = proc.communicate()
        return SageBaseTarFile(mode="r", fileobj=BytesIO(data))

    @staticmethod
    def can_read(filename):
        """
        Given an archive filename, returns True if this class can read and
        process the archive format of that file.
        """
        devnull = open(os.devnull, 'w')
        try:
            subprocess.check_call(["xz", "-l", filename], stdout=devnull, stderr=devnull)
        except Exception:
            return False
        return True
