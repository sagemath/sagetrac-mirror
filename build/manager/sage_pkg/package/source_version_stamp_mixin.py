"""
Mixin for richer dependencies

This mixin lets the version stamp (the only version that is used for
dependency computatons) not only depend on the package directory, but
also no other directories of the repository.
"""
 
from sage_pkg.config import config
from sage_pkg.utils import random_sha1, cached_property
from sage_pkg.git.repo import GitRepository
from sage_pkg.package.dependency_accumulator import CommutingSha1Accumulator


class SourceVersionStampMixin(object):

    def _directory_version_stamp(self, directory):
        git = GitRepository(config.path.dot_git)
        if git.is_clean_dir(directory):
            tree = git.get_tree(directory)
            return tree.sha1
        else:
            return random_sha1()

    @cached_property
    def version_stamp(self):
        """
        Unique machine-readable version.

        Some builder types allow you to specify a list of source
        directories as ``config.source``. Using this mixin, their git
        tree sha1 will be added to the package directory sha1. Hence
        it depends on both directory contents.

        OUTPUT:

        String (ASCII SHA1).
        """
        acc = CommutingSha1Accumulator()
        acc += self._version_stamp
        for source_dir in self._config.source:
            acc += self._directory_version_stamp(source_dir)
        return str(acc)
