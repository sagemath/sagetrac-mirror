# -*- encoding: utf-8 -*-
r"""
Output Buffer

This is the fundamental unit of rich output, a single immutable buffer
(either in-memory or as a file). Rich output always consists of one or
more buffers. Ideally, the Sage library always uses the buffer object
as an in-memory buffer. But you can also ask it for a filename, and it
will save the data to a file if necessary. Either way, the buffer
object presents the same interface for getting the content of an
in-memory buffer or a temporary file. So any rich output backends do
not need to know where the buffer content is actually stored.

EXAMPLES::

    sage: from sage.repl.rich_output.buffer import OutputBuffer
    sage: buf = OutputBuffer('this is the buffer content');  buf
    buffer containing 26 bytes
    sage: buf.get()
    'this is the buffer content'
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


import os
from sage.structure.sage_object import SageObject


class OutputBuffer(SageObject):

    def __init__(self, data):
        """
        Data stored either in memory or as a file

        This class is an abstraction for "files", in that they can
        either be defined by a bytes array (Python 3) or string
        (Python 2) or by a file (see :meth:`from_file`).

        INPUT:

        - ``data`` -- bytes. The data that is stored in the buffer.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: buf = OutputBuffer('this is the buffer content');  buf
            buffer containing 26 bytes

            sage: buf2 = OutputBuffer(buf);  buf2
            buffer containing 26 bytes

            sage: buf.get()
            'this is the buffer content'
            sage: buf.filename(ext='.txt')
            '/....txt'
        """
        if isinstance(data, OutputBuffer):
            self._filename = data._filename
            self._data = data._data
        else:
            self._filename = None
            self._data = bytes(data)

    @classmethod
    def from_file(cls, filename):
        """
        Construct buffer from data in file.

        .. warning::

            The buffer assumes that the file content remains the same
            during the lifetime of the Sage session. To communicate
            this to the user, the file permissions will be changed to
            read only.

        INPUT:

        - ``filename`` -- string. The filename under which the data is
          stored.

        OUTPUT:
        
        String containing the buffer data.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: name = sage.misc.temporary_file.tmp_filename()
            sage: with open(name, 'w') as f:
            ....:    f.write('file content')
            sage: buf = OutputBuffer.from_file(name);  buf
            buffer containing 12 bytes

            sage: buf.filename() == name
            True
            sage: buf.get()
            'file content'
        """
        buf = cls.__new__(cls)
        buf._filename = filename
        buf._data = None
        buf._chmod_readonly(buf._filename)
        return buf

    @classmethod
    def _chmod_readonly(cls, filename):
        """
        Make file readonly

        INPUT:

        - ``filename`` -- string. Name of an already-existing file.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: tmp = sage.misc.temporary_file.tmp_filename()
            sage: with open(tmp, 'w') as f:
            ....:    f.write('file content')
            sage: OutputBuffer._chmod_readonly(tmp)
            sage: import os, stat
            sage: stat.S_IMODE(os.stat(tmp).st_mode) & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
            0
        """
        import os
        from sage.env import SAGE_EXTCODE
        filename = os.path.abspath(filename)
        if filename.startswith(os.path.abspath(SAGE_EXTCODE)):
            # Do not change permissions on the sample rich output
            # files, as it will cause trouble when upgrading Sage
            return
        import stat
        mode = os.stat(filename).st_mode
        mode = stat.S_IMODE(mode) & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
        os.chmod(filename, mode)
        
    def _repr_(self):
        """
        Return a string representation

        OUTPUT:

        String

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: OutputBuffer('test1234')
            buffer containing 8 bytes
        """
        return 'buffer containing {0} bytes'.format(len(self.get()))
        
    def get(self):
        """
        Return the buffer content
        
        OUTPUT:

        Bytes. A string in Python 2.x.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: OutputBuffer('test1234').get()
            'test1234'
        """
        if self._data is None:
            with open(self._filename) as f:
                self._data = f.read()
        return self._data

    def filename(self, ext=None):
        """
        Return the filename.

        INPUT:

        - ``ext`` -- string. The file extension.

        OUTPUT:

        Name of a file, most likely a temporary file. If ``ext`` is
        specified, the filename will have that extension.

        You must not modify the returned file. Its permissions are set
        to readonly to help with that.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: buf = OutputBuffer('test')
            sage: buf.filename()           # random output
            '/home/user/.sage/temp/hostname/26085/tmp_RNSfAc'

            sage: os.path.isfile(buf.filename())
            True
            sage: buf.filename(ext='txt')  # random output
            '/home/user/.sage/temp/hostname/26085/tmp_Rjjp4V.txt'
            sage: buf.filename(ext='txt').endswith('.txt')
            True
        """
        if ext is None:
            ext = ''
        elif not ext.startswith('.'):
            ext = '.' + ext

        if self._filename is None or not self._filename.endswith(ext):
            from sage.misc.temporary_file import tmp_filename
            output = tmp_filename(ext=ext)
        else:
            output = self._filename

        if self._filename is None:
            assert self._data is not None
            with open(output, 'wb') as f:
                f.write(self._data)
            self._filename = output
        elif self._filename != output:
            try:
                os.link(self._filename, output)
            except (OSError, AttributeError):
                import shutil
                shutil.copy2(self._filename, output)

        self._chmod_readonly(output)
        return output
        
    def save_as(self, filename):
        """
        Save a copy of the buffer content.

        You may edit the returned file, unlike the file returned by
        :meth:`filename`.

        INPUT:

        - ``filename`` -- string. The file name to save under.

        EXAMPLES::

            sage: from sage.repl.rich_output.buffer import OutputBuffer
            sage: buf = OutputBuffer('test')
            sage: buf.filename(ext='txt')  # random output
            sage: tmp = tmp_dir()
            sage: filename = os.path.join(tmp, 'foo.txt')
            sage: buf.save_as(filename)
            sage: with open(filename, 'r') as f:
            ....:     f.read()
            'test'
        """
        with open(filename, 'w') as f:
            f.write(self.get())
