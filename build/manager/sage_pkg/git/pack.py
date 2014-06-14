"""
Packed Git Objects

Taken from the Dulwich Git client and all write functionality
removed.

# EXAMPLES::

#     >>> import os
#     >>> from sage_pkg.git.pack import load_pack_index, Pack
#     >>> p = load_pack_index(os.path.join(config.path.dot_git, 'objects', 'pack', 'pack-1d714b665805cc32d140495c84abb23a8ec22c34.idx'))
#     >>> p = load_pack_index(os.path.join(config.path.dot_git, 'objects', 'pack', 'pack-40ca8d1eac39ffa0fc6a071ae834e59e1a49ef9f.idx'))
#     >>> # p.object_index('b16df56703138b0e6bed3e0bb48ce0fa276cb575')
#     >>> p = Pack(os.path.join(config.path.dot_git, 'objects', 'pack', 'pack-40ca8d1eac39ffa0fc6a071ae834e59e1a49ef9f'))
#     >>> p.get_raw('b16df56703138b0e6bed3e0bb48ce0fa276cb575')
"""


# pack.py -- For dealing with packed git objects.
# Copyright (C) 2007 James Westby <jw+debian@jameswestby.net>
# Copyright (C) 2008-2013 Jelmer Vernooij <jelmer@samba.org>
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; version 2
# of the License or (at your option) a later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA  02110-1301, USA.

"""Classes for dealing with packed git objects.

A pack is a compact representation of a bunch of objects, stored
using deltas where possible.

They have two parts, the pack file, which stores the data, and an index
that tells you where the data is.

To find an object you look in all of the index files 'til you find a
match for the object name. You then use the pointer got from this as
a pointer in to the corresponding packfile.
"""


import binascii

from collections import defaultdict

from io import BytesIO
from collections import (
    deque,
    )

from itertools import chain, imap, izip

import mmap
from hashlib import sha1
import os
from os import (
    SEEK_CUR,
    SEEK_END,
    )
import struct
from struct import unpack_from
import warnings
import zlib

from sage_pkg.git.util import sha_to_binary


OFS_DELTA = 6
REF_DELTA = 7

DELTA_TYPES = (OFS_DELTA, REF_DELTA)

DEFAULT_PACK_DELTA_WINDOW_SIZE = 10


class ChecksumMismatch(Exception):
    pass

def take_msb_bytes(read, crc32=None):
    """Read bytes marked with most significant bit.

    :param read: Read function
    """
    ret = []
    while len(ret) == 0 or ret[-1] & 0x80:
        b = read(1)
        if crc32 is not None:
            crc32 = binascii.crc32(b, crc32)
        ret.append(ord(b))
    return ret, crc32


class UnpackedObject(object):
    """Class encapsulating an object unpacked from a pack file.

    These objects should only be created from within unpack_object. Most
    members start out as empty and are filled in at various points by
    read_zlib_chunks, unpack_object, DeltaChainIterator, etc.

    End users of this object should take care that the function they're getting
    this object from is guaranteed to set the members they need.
    """

    __slots__ = [
      'offset',         # Offset in its pack.
      '_sha',           # Cached binary SHA.
      'obj_type_num',   # Type of this object.
      'obj_chunks',     # Decompressed and delta-resolved chunks.
      'pack_type_num',  # Type of this object in the pack (may be a delta).
      'delta_base',     # Delta base offset or SHA.
      'comp_chunks',    # Compressed object chunks.
      'decomp_chunks',  # Decompressed object chunks.
      'decomp_len',     # Decompressed length of this object.
      'crc32',          # CRC32.
      ]

    # TODO(dborowitz): read_zlib_chunks and unpack_object could very well be
    # methods of this object.
    def __init__(self, pack_type_num, delta_base, decomp_len, crc32):
        self.offset = None
        self._sha = None
        self.pack_type_num = pack_type_num
        self.delta_base = delta_base
        self.comp_chunks = None
        self.decomp_chunks = []
        self.decomp_len = decomp_len
        self.crc32 = crc32

        if pack_type_num in DELTA_TYPES:
            self.obj_type_num = None
            self.obj_chunks = None
        else:
            self.obj_type_num = pack_type_num
            self.obj_chunks = self.decomp_chunks
            self.delta_base = delta_base

    # Only provided for backwards compatibility with code that expects either
    # chunks or a delta tuple.
    def _obj(self):
        """Return the decompressed chunks, or (delta base, delta chunks)."""
        if self.pack_type_num in DELTA_TYPES:
            return (self.delta_base, self.decomp_chunks)
        else:
            return self.decomp_chunks

    def __eq__(self, other):
        if not isinstance(other, UnpackedObject):
            return False
        for slot in self.__slots__:
            if getattr(self, slot) != getattr(other, slot):
                return False
        return True

    def __ne__(self, other):
        return not (self == other)

    def __repr__(self):
        data = ['%s=%r' % (s, getattr(self, s)) for s in self.__slots__]
        return '%s(%s)' % (self.__class__.__name__, ', '.join(data))


_ZLIB_BUFSIZE = 4096


def read_zlib_chunks(read_some, unpacked, include_comp=False,
                     buffer_size=_ZLIB_BUFSIZE):
    """Read zlib data from a buffer.

    This function requires that the buffer have additional data following the
    compressed data, which is guaranteed to be the case for git pack files.

    :param read_some: Read function that returns at least one byte, but may
        return less than the requested size.
    :param unpacked: An UnpackedObject to write result data to. If its crc32
        attr is not None, the CRC32 of the compressed bytes will be computed
        using this starting CRC32.
        After this function, will have the following attrs set:
        * comp_chunks    (if include_comp is True)
        * decomp_chunks
        * decomp_len
        * crc32
    :param include_comp: If True, include compressed data in the result.
    :param buffer_size: Size of the read buffer.
    :return: Leftover unused data from the decompression.
    :raise zlib.error: if a decompression error occurred.
    """
    if unpacked.decomp_len <= -1:
        raise ValueError('non-negative zlib data stream size expected')
    decomp_obj = zlib.decompressobj()

    comp_chunks = []
    decomp_chunks = unpacked.decomp_chunks
    decomp_len = 0
    crc32 = unpacked.crc32

    while True:
        add = read_some(buffer_size)
        if not add:
            raise zlib.error('EOF before end of zlib stream')
        comp_chunks.append(add)
        decomp = decomp_obj.decompress(add)
        decomp_len += len(decomp)
        decomp_chunks.append(decomp)
        unused = decomp_obj.unused_data
        if unused:
            left = len(unused)
            if crc32 is not None:
                crc32 = binascii.crc32(add[:-left], crc32)
            if include_comp:
                comp_chunks[-1] = add[:-left]
            break
        elif crc32 is not None:
            crc32 = binascii.crc32(add, crc32)
    if crc32 is not None:
        crc32 &= 0xffffffff

    if decomp_len != unpacked.decomp_len:
        raise zlib.error('decompressed data does not match expected size')

    unpacked.crc32 = crc32
    if include_comp:
        unpacked.comp_chunks = comp_chunks
    return unused


def iter_sha1(iter):
    """Return the hexdigest of the SHA1 over a set of names.

    :param iter: Iterator over string objects
    :return: 40-byte hex sha1 digest
    """
    sha = sha1()
    for name in iter:
        sha.update(name)
    return sha.hexdigest()


def load_pack_index(path):
    """Load an index file by path.

    :param filename: Path to the index file
    :return: A PackIndex loaded from the given path
    """
    with open(path, 'rb') as f:
        return load_pack_index_file(path, f)


def load_pack_index_file(path, f):
    """Load an index file from a file-like object.

    :param path: Path for the index file
    :param f: File-like object
    :return: A PackIndex loaded from the given file
    """
    contents, size = _load_file_contents(f)
    if contents[:4] != '\377tOc':
        raise ValueError('old format or invalid pack index')
    version = struct.unpack('>L', contents[4:8])[0]
    if version != 2:
        raise ValueError('only pack version 2 is supported')        
    return PackIndex2(path, f, contents, size)


def _load_file_contents(f, size=None):
    fileno = getattr(f, 'fileno', None)
    assert fileno is not None  # must be actual file
    fd = f.fileno()
    if size is None:
        size = os.fstat(fd).st_size
    contents = mmap.mmap(fd, size, access=mmap.ACCESS_READ)
    return contents, size


def bisect_find_sha(start, end, sha, unpack_name):
    """Find a SHA in a data blob with sorted SHAs.

    :param start: Start index of range to search
    :param end: End index of range to search
    :param sha: Sha to find
    :param unpack_name: Callback to retrieve SHA by index
    :return: Index of the SHA, or None if it wasn't found
    """
    assert start <= end
    while start <= end:
        i = (start + end) // 2
        file_sha = unpack_name(i)
        x = cmp(file_sha, sha)
        if x < 0:
            start = i + 1
        elif x > 0:
            end = i - 1
        else:
            return i
    return None


class PackIndex(object):
    """An index in to a packfile.

    Given a sha id of an object a pack index can tell you the location in the
    packfile of that object if it has it.
    """

    def __eq__(self, other):
        if not isinstance(other, PackIndex):
            return False

        for (name1, _, _), (name2, _, _) in izip(self.iterentries(),
                                                 other.iterentries()):
            if name1 != name2:
                return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __len__(self):
        """Return the number of entries in this pack index."""
        raise NotImplementedError(self.__len__)

    def get_pack_checksum(self):
        """Return the SHA1 checksum stored for the corresponding packfile.

        :return: 20-byte binary digest
        """
        raise NotImplementedError(self.get_pack_checksum)

    def object_index(self, sha):
        """Return the index in to the corresponding packfile for the object.

        Given the name of an object it will return the offset that object
        lives at within the corresponding pack file. If the pack file doesn't
        have the object then None will be returned.
        """
        if len(sha) == 40:
            sha = sha_to_binary(sha)
        return self._object_index(sha)

    def _object_index(self, sha):
        """See object_index.

        :param sha: A *binary* SHA string. (20 characters long)_
        """
        raise NotImplementedError(self._object_index)



class FilePackIndex(PackIndex):
    """Pack index that is based on a file.

    To do the loop it opens the file, and indexes first 256 4 byte groups
    with the first byte of the sha id. The value in the four byte group indexed
    is the end of the group that shares the same starting byte. Subtract one
    from the starting byte and index again to find the start of the group.
    The values are sorted by sha id within the group, so do the math to find
    the start and end offset and then bisect in to find if the value is present.
    """

    def __init__(self, filename, file, contents, size):
        """Create a pack index object.

        Provide it with the name of the index file to consider, and it will map
        it whenever required.
        """
        self._filename = filename
        self._file = file
        self._contents = contents
        self._size = size

    def __eq__(self, other):
        # Quick optimization:
        if (isinstance(other, FilePackIndex) and
            self._fan_out_table != other._fan_out_table):
            return False

        return super(FilePackIndex, self).__eq__(other)

    def close(self):
        self._file.close()
        if getattr(self._contents, "close", None) is not None:
            self._contents.close()

    def __len__(self):
        """Return the number of entries in this pack index."""
        return self._fan_out_table[-1]

    def _unpack_entry(self, i):
        """Unpack the i-th entry in the index file.

        :return: Tuple with object name (SHA), offset in pack file and CRC32
            checksum (if known).
        """
        raise NotImplementedError(self._unpack_entry)

    def _unpack_name(self, i):
        """Unpack the i-th name from the index file."""
        raise NotImplementedError(self._unpack_name)

    def _unpack_offset(self, i):
        """Unpack the i-th object offset from the index file."""
        raise NotImplementedError(self._unpack_offset)

    def _unpack_crc32_checksum(self, i):
        """Unpack the crc32 checksum for the i-th object from the index file."""
        raise NotImplementedError(self._unpack_crc32_checksum)

    def _read_fan_out_table(self, start_offset):
        ret = []
        for i in range(0x100):
            fanout_entry = self._contents[start_offset+i*4:start_offset+(i+1)*4]
            ret.append(struct.unpack('>L', fanout_entry)[0])
        return ret

    def calculate_checksum(self):
        """Calculate the SHA1 checksum over this pack index.

        :return: This is a 20-byte binary digest
        """
        return sha1(self._contents[:-20]).digest()

    def get_pack_checksum(self):
        """Return the SHA1 checksum stored for the corresponding packfile.

        :return: 20-byte binary digest
        """
        return str(self._contents[-40:-20])

    def _object_index(self, sha):
        """
        See object_index.

        :param sha: A *binary* SHA string. (20 characters long)_
        """
        assert len(sha) == 20
        idx = ord(sha[0])
        if idx == 0:
            start = 0
        else:
            start = self._fan_out_table[idx-1]
        end = self._fan_out_table[idx]
        i = bisect_find_sha(start, end, sha, self._unpack_name)
        if i is None:
            raise KeyError(sha)
        return self._unpack_offset(i)


class PackIndex2(FilePackIndex):
    """Version 2 Pack Index file."""

    def __init__(self, filename, file, contents, size):
        super(PackIndex2, self).__init__(filename, file, contents, size)
        if self._contents[:4] != '\377tOc':
            raise AssertionError('Not a v2 pack index file')
        (self.version, ) = unpack_from('>L', self._contents, 4)
        if self.version != 2:
            raise AssertionError('Version was %d' % self.version)
        self._fan_out_table = self._read_fan_out_table(8)
        self._name_table_offset = 8 + 0x100 * 4
        self._crc32_table_offset = self._name_table_offset + 20 * len(self)
        self._pack_offset_table_offset = (self._crc32_table_offset +
                                          4 * len(self))
        self._pack_offset_largetable_offset = (self._pack_offset_table_offset +
                                          4 * len(self))

    def _unpack_entry(self, i):
        return (self._unpack_name(i), self._unpack_offset(i),
                self._unpack_crc32_checksum(i))

    def _unpack_name(self, i):
        offset = self._name_table_offset + i * 20
        return self._contents[offset:offset+20]

    def _unpack_offset(self, i):
        offset = self._pack_offset_table_offset + i * 4
        offset = unpack_from('>L', self._contents, offset)[0]
        if offset & (2**31):
            offset = self._pack_offset_largetable_offset + (offset&(2**31-1)) * 8
            offset = unpack_from('>Q', self._contents, offset)[0]
        return offset



def read_pack_header(read):
    """Read the header of a pack file.

    :param read: Read function
    :return: Tuple of (pack version, number of objects). If no data is available
        to read, returns (None, None).
    """
    header = read(12)
    if not header:
        return None, None
    if header[:4] != 'PACK':
        raise AssertionError('Invalid pack header %r' % header)
    (version,) = unpack_from('>L', header, 4)
    if version not in (2, 3):
        raise AssertionError('Version was %d' % version)
    (num_objects,) = unpack_from('>L', header, 8)
    return (version, num_objects)


def chunks_length(chunks):
    return sum(imap(len, chunks))


def unpack_object(read_all, read_some=None, compute_crc32=False,
                  include_comp=False, zlib_bufsize=_ZLIB_BUFSIZE):
    """Unpack a Git object.

    :param read_all: Read function that blocks until the number of requested
        bytes are read.
    :param read_some: Read function that returns at least one byte, but may not
        return the number of bytes requested.
    :param compute_crc32: If True, compute the CRC32 of the compressed data. If
        False, the returned CRC32 will be None.
    :param include_comp: If True, include compressed data in the result.
    :param zlib_bufsize: An optional buffer size for zlib operations.
    :return: A tuple of (unpacked, unused), where unused is the unused data
        leftover from decompression, and unpacked in an UnpackedObject with
        the following attrs set:

        * obj_chunks     (for non-delta types)
        * pack_type_num
        * delta_base     (for delta types)
        * comp_chunks    (if include_comp is True)
        * decomp_chunks
        * decomp_len
        * crc32          (if compute_crc32 is True)
    """
    if read_some is None:
        read_some = read_all
    if compute_crc32:
        crc32 = 0
    else:
        crc32 = None

    bytes, crc32 = take_msb_bytes(read_all, crc32=crc32)
    type_num = (bytes[0] >> 4) & 0x07
    size = bytes[0] & 0x0f
    for i, byte in enumerate(bytes[1:]):
        size += (byte & 0x7f) << ((i * 7) + 4)

    raw_base = len(bytes)
    if type_num == OFS_DELTA:
        bytes, crc32 = take_msb_bytes(read_all, crc32=crc32)
        raw_base += len(bytes)
        if bytes[-1] & 0x80:
            raise AssertionError
        delta_base_offset = bytes[0] & 0x7f
        for byte in bytes[1:]:
            delta_base_offset += 1
            delta_base_offset <<= 7
            delta_base_offset += (byte & 0x7f)
        delta_base = delta_base_offset
    elif type_num == REF_DELTA:
        delta_base = read_all(20)
        if compute_crc32:
            crc32 = binascii.crc32(delta_base, crc32)
        raw_base += 20
    else:
        delta_base = None

    unpacked = UnpackedObject(type_num, delta_base, size, crc32)
    unused = read_zlib_chunks(read_some, unpacked, buffer_size=zlib_bufsize,
                              include_comp=include_comp)
    return unpacked, unused



class PackStreamReader(object):
    """Class to read a pack stream.

    The pack is read from a ReceivableProtocol using read() or recv() as
    appropriate.
    """

    def __init__(self, read_all, read_some=None, zlib_bufsize=_ZLIB_BUFSIZE):
        self.read_all = read_all
        if read_some is None:
            self.read_some = read_all
        else:
            self.read_some = read_some
        self.sha = sha1()
        self._offset = 0
        self._rbuf = BytesIO()
        # trailer is a deque to avoid memory allocation on small reads
        self._trailer = deque()
        self._zlib_bufsize = zlib_bufsize

    def _read(self, read, size):
        """Read up to size bytes using the given callback.

        As a side effect, update the verifier's hash (excluding the last 20
        bytes read).

        :param read: The read callback to read from.
        :param size: The maximum number of bytes to read; the particular
            behavior is callback-specific.
        """
        data = read(size)

        # maintain a trailer of the last 20 bytes we've read
        n = len(data)
        self._offset += n
        tn = len(self._trailer)
        if n >= 20:
            to_pop = tn
            to_add = 20
        else:
            to_pop = max(n + tn - 20, 0)
            to_add = n
        for _ in xrange(to_pop):
            self.sha.update(self._trailer.popleft())
        self._trailer.extend(data[-to_add:])

        # hash everything but the trailer
        self.sha.update(data[:-to_add])
        return data

    def _buf_len(self):
        buf = self._rbuf
        start = buf.tell()
        buf.seek(0, SEEK_END)
        end = buf.tell()
        buf.seek(start)
        return end - start

    @property
    def offset(self):
        return self._offset - self._buf_len()

    def read(self, size):
        """Read, blocking until size bytes are read."""
        buf_len = self._buf_len()
        if buf_len >= size:
            return self._rbuf.read(size)
        buf_data = self._rbuf.read()
        self._rbuf = BytesIO()
        return buf_data + self._read(self.read_all, size - buf_len)

    def recv(self, size):
        """Read up to size bytes, blocking until one byte is read."""
        buf_len = self._buf_len()
        if buf_len:
            data = self._rbuf.read(size)
            if size >= buf_len:
                self._rbuf = BytesIO()
            return data
        return self._read(self.read_some, size)

    def __len__(self):
        return self._num_objects

    def read_objects(self, compute_crc32=False):
        """Read the objects in this pack file.

        :param compute_crc32: If True, compute the CRC32 of the compressed
            data. If False, the returned CRC32 will be None.
        :return: Iterator over UnpackedObjects with the following members set:
            offset
            obj_type_num
            obj_chunks (for non-delta types)
            delta_base (for delta types)
            decomp_chunks
            decomp_len
            crc32 (if compute_crc32 is True)
        :raise ChecksumMismatch: if the checksum of the pack contents does not
            match the checksum in the pack trailer.
        :raise zlib.error: if an error occurred during zlib decompression.
        :raise IOError: if an error occurred writing to the output file.
        """
        pack_version, self._num_objects = read_pack_header(self.read)
        if pack_version is None:
            return

        for i in xrange(self._num_objects):
            offset = self.offset
            unpacked, unused = unpack_object(
              self.read, read_some=self.recv, compute_crc32=compute_crc32,
              zlib_bufsize=self._zlib_bufsize)
            unpacked.offset = offset

            # prepend any unused data to current read buffer
            buf = BytesIO()
            buf.write(unused)
            buf.write(self._rbuf.read())
            buf.seek(0)
            self._rbuf = buf

            yield unpacked

        if self._buf_len() < 20:
            # If the read buffer is full, then the last read() got the whole
            # trailer off the wire. If not, it means there is still some of the
            # trailer to read. We need to read() all 20 bytes; N come from the
            # read buffer and (20 - N) come from the wire.
            self.read(20)

        pack_sha = ''.join(self._trailer)
        if pack_sha != self.sha.digest():
            raise ChecksumMismatch()




class PackData(object):
    """The data contained in a packfile.

    Pack files can be accessed both sequentially for exploding a pack, and
    directly with the help of an index to retrieve a specific object.

    The objects within are either complete or a delta aginst another.

    The header is variable length. If the MSB of each byte is set then it
    indicates that the subsequent byte is still part of the header.
    For the first byte the next MS bits are the type, which tells you the type
    of object, and whether it is a delta. The LS byte is the lowest bits of the
    size. For each subsequent byte the LS 7 bits are the next MS bits of the
    size, i.e. the last byte of the header contains the MS bits of the size.

    For the complete objects the data is stored as zlib deflated data.
    The size in the header is the uncompressed object size, so to uncompress
    you need to just keep feeding data to zlib until you get an object back,
    or it errors on bad data. This is done here by just giving the complete
    buffer from the start of the deflated object on. This is bad, but until I
    get mmap sorted out it will have to do.

    Currently there are no integrity checks done. Also no attempt is made to
    try and detect the delta case, or a request for an object at the wrong
    position.  It will all just throw a zlib or KeyError.
    """

    def __init__(self, filename, file=None, size=None):
        """
        Create a PackData object representing the pack in the given filename.

        The file must exist and stay readable until the object is disposed of. It
        must also stay the same size. It will be mapped whenever needed.

        Currently there is a restriction on the size of the pack as the python
        mmap implementation is flawed.
        """
        self._filename = filename
        self._size = size
        self._header_size = 12
        if file is None:
            self._file = open(self._filename, 'rb')
        else:
            self._file = file
        (version, self._num_objects) = read_pack_header(self._file.read)
        self.pack = None

    @property
    def filename(self):
        return os.path.basename(self._filename)

    @classmethod
    def from_file(cls, file, size):
        return cls(str(file), file=file, size=size)

    @classmethod
    def from_path(cls, path):
        return cls(filename=path)

    def close(self):
        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def _get_size(self):
        if self._size is not None:
            return self._size
        self._size = os.path.getsize(self._filename)
        if self._size < self._header_size:
            errmsg = ('%s is too small for a packfile (%d < %d)' %
                      (self._filename, self._size, self._header_size))
            raise AssertionError(errmsg)
        return self._size

    def __len__(self):
        """Returns the number of objects in this pack."""
        return self._num_objects

    def calculate_checksum(self):
        """Calculate the checksum for this pack.

        :return: 20-byte binary SHA1 digest
        """
        return compute_file_sha(self._file, end_ofs=-20).digest()

    def get_ref(self, sha):
        """Get the object for a ref SHA, only looking in this pack."""
        # TODO: cache these results
        if self.pack is None:
            raise KeyError(sha)
        offset = self.pack.index.object_index(sha)
        type, obj = self.get_object_at(offset)
        return offset, type, obj

    def resolve_object(self, offset, type, obj, get_ref=None):
        """Resolve an object, possibly resolving deltas when necessary.

        :return: Tuple with object type and contents.
        """
        if type not in DELTA_TYPES:
            return type, obj

        if get_ref is None:
            get_ref = self.get_ref
        if type == OFS_DELTA:
            (delta_offset, delta) = obj
            # TODO: clean up asserts and replace with nicer error messages
            assert isinstance(offset, int) or isinstance(offset, long)
            assert isinstance(delta_offset, int) or isinstance(offset, long)
            base_offset = offset-delta_offset
            type, base_obj = self.get_object_at(base_offset)
            assert isinstance(type, int)
        elif type == REF_DELTA:
            (basename, delta) = obj
            assert isinstance(basename, str) and len(basename) == 20
            base_offset, type, base_obj = get_ref(basename)
            assert isinstance(type, int)
        type, base_chunks = self.resolve_object(base_offset, type, base_obj)
        chunks = apply_delta(base_chunks, delta)
        return type, chunks

    def get_stored_checksum(self):
        """Return the expected checksum stored in this pack."""
        self._file.seek(-20, SEEK_END)
        return self._file.read(20)

    def get_object_at(self, offset):
        """Given an offset in to the packfile return the object that is there.

        Using the associated index the location of an object can be looked up,
        and then the packfile can be asked directly for that object using this
        function.
        """
        assert isinstance(offset, long) or isinstance(offset, int),\
                'offset was %r' % offset
        assert offset >= self._header_size
        self._file.seek(offset)
        unpacked, _ = unpack_object(self._file.read)
        return (unpacked.pack_type_num, unpacked._obj())



def apply_delta(src_buf, delta):
    """Based on the similar function in git's patch-delta.c.

    :param src_buf: Source buffer
    :param delta: Delta instructions
    """
    if not isinstance(src_buf, str):
        src_buf = ''.join(src_buf)
    if not isinstance(delta, str):
        delta = ''.join(delta)
    out = []
    index = 0
    delta_length = len(delta)
    def get_delta_header_size(delta, index):
        size = 0
        i = 0
        while delta:
            cmd = ord(delta[index])
            index += 1
            size |= (cmd & ~0x80) << i
            i += 7
            if not cmd & 0x80:
                break
        return size, index
    src_size, index = get_delta_header_size(delta, index)
    dest_size, index = get_delta_header_size(delta, index)
    assert src_size == len(src_buf), '%d vs %d' % (src_size, len(src_buf))
    while index < delta_length:
        cmd = ord(delta[index])
        index += 1
        if cmd & 0x80:
            cp_off = 0
            for i in range(4):
                if cmd & (1 << i):
                    x = ord(delta[index])
                    index += 1
                    cp_off |= x << (i * 8)
            cp_size = 0
            for i in range(3):
                if cmd & (1 << (4+i)):
                    x = ord(delta[index])
                    index += 1
                    cp_size |= x << (i * 8)
            if cp_size == 0:
                cp_size = 0x10000
            if (cp_off + cp_size < cp_size or
                cp_off + cp_size > src_size or
                cp_size > dest_size):
                break
            out.append(src_buf[cp_off:cp_off+cp_size])
        elif cmd != 0:
            out.append(delta[index:index+cmd])
            index += cmd
        else:
            raise ApplyDeltaError('Invalid opcode 0')

    if index != delta_length:
        raise ApplyDeltaError('delta not empty: %r' % delta[index:])

    if dest_size != chunks_length(out):
        raise ApplyDeltaError('dest size incorrect')

    return out


class Pack(object):
    """A Git pack object."""

    def __init__(self, basename, resolve_ext_ref=None):
        self._basename = basename
        self._data = None
        self._idx = None
        self._idx_path = self._basename + '.idx'
        self._data_path = self._basename + '.pack'
        self._data_load = lambda: PackData(self._data_path)
        self._idx_load = lambda: load_pack_index(self._idx_path)
        self.resolve_ext_ref = resolve_ext_ref

    @classmethod
    def from_lazy_objects(self, data_fn, idx_fn):
        """Create a new pack object from callables to load pack data and
        index objects."""
        ret = Pack('')
        ret._data_load = data_fn
        ret._idx_load = idx_fn
        return ret

    @classmethod
    def from_objects(self, data, idx):
        """Create a new pack object from pack data and index objects."""
        ret = Pack('')
        ret._data_load = lambda: data
        ret._idx_load = lambda: idx
        return ret

    def name(self):
        """The SHA over the SHAs of the objects in this pack."""
        return self.index.objects_sha1()

    @property
    def data(self):
        """The pack data object being used."""
        if self._data is None:
            self._data = self._data_load()
            self._data.pack = self
            self.check_length_and_checksum()
        return self._data

    @property
    def index(self):
        """The index being used.

        :note: This may be an in-memory index
        """
        if self._idx is None:
            self._idx = self._idx_load()
        return self._idx

    def close(self):
        if self._data is not None:
            self._data.close()
        if self._idx is not None:
            self._idx.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __eq__(self, other):
        return isinstance(self, type(other)) and self.index == other.index

    def __len__(self):
        """Number of entries in this pack."""
        return len(self.index)

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, self._basename)

    def check_length_and_checksum(self):
        """Sanity check the length and checksum of the pack index and data."""
        assert len(self.index) == len(self.data)
        idx_stored_checksum = self.index.get_pack_checksum()
        data_stored_checksum = self.data.get_stored_checksum()
        if idx_stored_checksum != data_stored_checksum:
            raise ChecksumMismatch()

    def check(self):
        """Check the integrity of this pack.

        :raise ChecksumMismatch: if a checksum for the index or data is wrong
        """
        self.index.check()
        self.data.check()
        for obj in self.iterobjects():
            obj.check()
        # TODO: object connectivity checks

    def get_stored_checksum(self):
        return self.data.get_stored_checksum()

    def get_raw(self, sha1):
        offset = self.index.object_index(sha1)
        obj_type, obj = self.data.get_object_at(offset)
        type_num, chunks = self.data.resolve_object(offset, obj_type, obj)
        return type_num, ''.join(chunks)



