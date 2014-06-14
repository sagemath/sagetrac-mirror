"""
Minimal Git Interface

This module implements a minimal git interface. We only use it to find
out the SHA1 of tree objects of the current branch.

Limitations:

* Readonly!

* Only new git repo format (e.g. pack index v2)
"""


import os
import re
import zlib

from sage_pkg.logger import logger
from sage_pkg.git.git_ignore import is_ignored

from sage_pkg.git.blob import Blob_from_file, Blob_from_pack, BlobTree, BlobCommit, BlobFile
from sage_pkg.git.util import full_split_repo_path


_blob_cache = dict()


class GitRepository(object):

    def __init__(self, dot_git):
        if not os.path.isdir(dot_git):
            raise ValueError('directory does not exist')
        self.dot_git = dot_git

    def __repr__(self):
        """
        Return a string representation

        EXAMPLES::

            >>> git    # doctest: +ELLIPSIS
            Git repository at .../.git
        """
        return 'Git repository at {0}'.format(self.dot_git)

    def get(self, sha1):
        """
        Load the blob with given sha1

        INPUT:

        - ``sha1`` -- 40-digit hex number as string. The sha1 of the
          object to load. Must be an unpacked blob.
        """
        global _blob_cache
        try:
            return _blob_cache[sha1]
        except KeyError:
            pass
        logger.debug('Loading blob with sha1 = %s', sha1)
        filename = os.path.join(self.dot_git, 'objects', sha1[0:2], sha1[2:])
        if os.path.isfile(filename):
            blob = Blob_from_file(filename)
        else:
            blob = Blob_from_pack(sha1)
        _blob_cache[sha1] = blob
        return blob

    def get_symbolic_ref(self, symbolic_ref):
        """
        Load the blob with given symbolic ref

        INPUT:

        - ``symbolic_ref`` -- string. A symbolic ref like ``'HEAD'``,
          ``'FETCH_HEAD'``, ....

        EXAMPLES::

            >>> git.get_symbolic_ref('HEAD')     # doctest: +ELLIPSIS
            BlobCommit:
                tree ...
        """
        logger.debug('Looking up symbolic ref = %s', symbolic_ref)
        filename = os.path.join(self.dot_git, symbolic_ref)
        with open(filename, 'r') as f:
            ref = f.read()
        before, after = ref.split('ref: ')
        if len(before) > 0:
            raise ValueError('not a symbolic ref')
        return self.get_branch(after.rstrip())

    def get_branch(self, branch):
        """
        Load the blob that the branch head points to

        INPUT:

        - ``branch`` -- a full branch name, starting with ``'refs/'``.
        """
        logger.debug('Looking up branch = %s', branch)
        if not branch.startswith('refs/'):
            raise ValueError('branch name must be full')
        filename = os.path.join(*([self.dot_git] + branch.split(os.path.sep)))
        if os.path.exists(filename):
            # Ref is unpacked
            with open(filename, 'r') as f:
                sha1 = f.read().rstrip()
            return self.get(sha1)
        # Otherwise, ref is packed
        packed_filename = os.path.join(self.dot_git, 'packed-refs')
        with open(packed_filename, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                if line.rstrip().endswith(branch):
                    sha1 = line[0:40]
                    return sha1
        raise ValueError('no such branch')
        
    def get_tree(self, dirname):
        """
        Return the tree object in the HEAD of the given directory

        INPUT:

        - ``dirname`` -- string. Directory name relative to the repository root.

        EXAMPLES::

            >>> subtree = git.get_tree('build/manager')
            >>> type(subtree)
            <class 'sage_pkg.git.blob.BlobTree'>
            >>> sha1 = subtree.get('sage-pkg')
            >>> len(sha1) == 40 and int(sha1, 16) >= 0
            True
        """
        commit = self.get_symbolic_ref('HEAD')
        tree = self.get(commit.tree)
        for name in full_split_repo_path(dirname):
            sha1 = tree.get(name)
            tree = self.get(sha1)
        return tree

    def is_clean_dir(self, dirname, verbose=False):
        """
        Check that the directory is clean (identical to the repo state)

        INPUT:

        - ``dirname`` -- string. The name of the directory to compare to.

        - ``verbose`` -- boolean.

        OUTPUT:

        Boolean.

        EXAMPLES::

            >>> git.is_clean_dir(config.path.packages, verbose=True)
            True
        """
        tree = self.get_tree(dirname)
        return self._is_clean_dir(dirname, tree, verbose)

    def _is_clean_dir(self, dirname, tree, verbose):
        dirent = set(os.listdir(dirname))
        for mode, name, sha1 in tree.ls():
            if name not in dirent:
                if verbose: print('deleted: ' + name)
                return False
            if mode == BlobTree.MODE_DIR:
                sub_dirname = os.path.join(dirname, name)
                if not os.path.isdir(sub_dirname):
                    if verbose: print('missing directory: ' + sub_dirname)
                    return False
                sub_tree = self.get(sha1)
                clean = self._is_clean_dir(sub_dirname, sub_tree, verbose)
            else:
                filename = os.path.join(dirname, name)
                clean = self._is_clean_file(filename, sha1, verbose)
            if not clean:
                return False
            dirent.remove(name)
        if any(not is_ignored(name) for name in dirent) > 0:
            if verbose: print('untracked files: ' + ', '.join(dirent))
            return False
        return True
        
    def _is_clean_file(self, filename, sha1, verbose):
        if not os.path.isfile(filename):
            if verbose: print('missing file: ' + filename)
            return False
        blob = self.get(sha1)
        with open(filename, 'r') as f:
            if f.read() != blob._content:
                if verbose: print('file modified: ' + filename)
                return False
        return True
        



