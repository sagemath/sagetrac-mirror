"""
Git Blobs
"""


import os
import re
import zlib
import glob

from sage_pkg.logger import logger
from sage_pkg.git.pack import Pack
from sage_pkg.config import config


BLOB_RE = re.compile(r'(?P<type>[a-z]*) (?P<size>[0-9]*)\0(?P<content>.*)', flags=re.DOTALL)

BLOB_COMMIT_TREE_RE = re.compile(r'tree (?P<sha1>[a-f0-9]{40,40})')


def Blob_from_file(filename, sha1):
    """
    Load the unpacked git object with given sha1 hash
    """
    with open(filename, 'rb') as f:
        blob = zlib.decompress(f.read())
    match = BLOB_RE.match(blob)
    if not match:
        raise ValueError('file is not a git object')
    if blob.startswith('commit'):
        return BlobCommit(match.group('content'), sha1)
    if blob.startswith('tree'):
        return BlobTree(match.group('content'), sha1)
    if blob.startswith('blob'):
        return BlobFile(match.group('content'), sha1)
    raise ValueError('unsupported blob: ' + repr(blob))


def Blob_from_pack(sha1):
    """
    Load object from pack
    
    EXAMPLES::

        >>> git.get('b16df56703138b0e6bed3e0bb48ce0fa276cb575')   # doctest: +ELLIPSIS
        BlobTree:
             40000 dcb048ff5ef59140c1fcd3bc67d80913bd821f36 atlas
             ...
             40000 19b12481f7f4267c59d11eab2c615fb0b9b6509c zn_poly
    """
    #P = Pack(os.path.join(config.path.dot_git, 'objects', 'pack', 'pack-40ca8d1eac39ffa0fc6a071ae834e59e1a49ef9f'))
    #return P.get_raw(sha1)

    pack_glob = os.path.join(config.path.dot_git, 'objects', 'pack', 'pack-*.pack')
    for pack_name in glob.glob(pack_glob):
        basename = pack_name[:-5]   # strip off '.pack'
        with Pack(basename) as P:
            try:
                type_num, content = P.get_raw(sha1)
            except KeyError:
                continue
        if type_num == 1:
            return BlobCommit(content, sha1)
        if type_num == 2:
            return BlobTree(content, sha1)
        if type_num == 3:
            return BlobFile(content, sha1)
        if type_num == 4:
            raise NotImplementedError('git tag object not implemented')
        raise ValueError('unsupported blob type_num: ' + str(type_num))
    raise ValueError('not in any pack file: ' + sha1)



class BlobABC(object):
    
    def __init__(self, content, sha1):
        self._content = content
        self.sha1 = sha1

    def __repr__(self):
        s = [self.__class__.__name__ + ':']
        for line in self._content.splitlines():
            if len(line.strip()) == 0:
                continue
            s.append('    ' + line.rstrip())
        return '\n'.join(s)


class BlobCommit(BlobABC):

    @property
    def tree(self):
        """
        Return the tree sha1 that is being committed
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> type(commit)
            <class 'sage_pkg.git.blob.BlobCommit'>
            >>> commit.tree    # doctest: +SKIP
            'f1795efaddd86273a715baf26a80c908a132a035'
            >>> len(commit.tree) == 40 and int(commit.tree, 16) >= 0
            True
        """
        match = BLOB_COMMIT_TREE_RE.match(self._content)
        return match.group('sha1')

    
class BlobTree(BlobABC):

    MODE_NORMAL = '100644'
    MODE_EXEC = '100755'
    MODE_DIR = '40000'

    def ls(self):
        """
        Iterate over the tree content
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> type(tree)
            <class 'sage_pkg.git.blob.BlobTree'>
            >>> for item in tree.ls():   # doctest: +ELLIPSIS
            ...     print(item)
            ('100644', '.gitignore', ...
        """
        pos = 0
        while pos < len(self._content):
            pos_space = self._content.find(' ', pos + 1)
            mode = self._content[pos:pos_space]
            pos_zero = self._content.find('\0', pos_space + 1)
            name = self._content[pos_space + 1:pos_zero]
            sha1_binary = self._content[pos_zero + 1:pos_zero + 21]
            sha1 = ''.join('{0:02x}'.format(ord(c)) for c in sha1_binary)
            yield (mode, name, sha1)
            pos = pos_zero + 21

    def __repr__(self):
        """
        Iterate over the tree content
        
        Override the generic repr since the tree object contains binary sha1's.

        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> git.get(commit.tree)   # doctest: +ELLIPSIS
            BlobTree:
                100644 ...
        """
        s = [self.__class__.__name__ + ':']
        for mode, name, sha1 in self.ls():
            s.append('    {0:>6} {1} {2}'.format(mode, sha1, name))
        return '\n'.join(s)

    def ls_dirs(self):
        """
        Iterater over the sub-trees
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> type(tree)
            <class 'sage_pkg.git.blob.BlobTree'>
            >>> for item in tree.ls_dirs():   # doctest: +ELLIPSIS
            ...     print(item)
            ('40000', 'build', ...
        """
        for dirent in self.ls():
            if dirent[0] == self.MODE_DIR:
                yield dirent

    def get(self, filename):
        """
        Return the sha1 of a file in the tree

        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> sha1 = tree.get('.gitignore')
            >>> sha1  # doctests: +SKIP
            'a9b6be08742e31b728f04f6f89c4b93f28ac4b92'
            >>> len(sha1) == 40 and int(sha1, 16) >= 0
            True
        """
        for mode, name, sha1 in self.ls():
            if name == filename:
                return sha1
        raise ValueError('file is not in the tree')

    def _compare_files(self, dirname, verbose=False):
        """
        Helper to just compare the files
        """


class BlobFile(BlobABC):
    pass
    

