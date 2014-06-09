"""
Minimal Git Interface

This module implements a minimal (and read-only) git interface. We
only use it to find out the SHA1 of tree objects.
"""

import os


class GitRepository(object):

    def __init__(self, dot_git):
        if not os.path.isdir(dot_git):
            raise ValueError('directory does not exist')
        self.dot_git = dot_git
