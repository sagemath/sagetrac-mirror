# -*- coding: utf-8 -*-
"""
Access the List of Sage Download Mirrors
"""

#*****************************************************************************
#       Copyright (C) 2015 Volker Braun <vbraun.name@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os
import contextlib
import logging
log = logging.getLogger()

from sage_bootstrap.compat import urllib, urlparse
from sage_bootstrap.env import SAGE_DISTFILES

from fcntl import flock, LOCK_SH, LOCK_EX
from errno import ENOLCK


def try_lock(fd, operation):
    """
    Try flock() but ignore ``ENOLCK`` errors, which could happen if the
    file system does not support locking.
    """
    try:
        flock(fd, operation)
    except IOError as e:
        if e.errno != ENOLCK:
            raise


class MirrorList(object):
    URL = 'http://www.sagemath.org/mirror_list'
    MAXAGE = 24*60*60   # seconds

    def __init__(self):
        self.filename = os.path.join(SAGE_DISTFILES, 'mirror_list')
        self.mirrors = None

        try:
            self.mirrorfile = open(self.filename, 'r+t')
        except IOError:
            self.mirrorfile = open(self.filename, 'w+t')

        with self.mirrorfile:
            self.mirrorfd = self.mirrorfile.fileno()
            try_lock(self.mirrorfd, LOCK_SH)  # shared (read) lock
            if self._must_refresh():
                try_lock(self.mirrorfd, LOCK_EX)  # exclusive (write) lock
                # Maybe the mirror list file was updated by a different
                # process while we waited for the lock?  Check again.
                if self._must_refresh():
                    self._refresh()
            if self.mirrors is None:
                self.mirrors = self._load()

    def _load(self, mirror_list=None):
        """
        Load and return `mirror_list` (defaults to the one on disk) as
        a list of strings
        """
        if mirror_list is None:
            try:
                self.mirrorfile.seek(0)
                mirror_list = self.mirrorfile.read()
            except IOError:
                log.critical('Failed to load the cached mirror list')
                return []
        import ast
        return ast.literal_eval(mirror_list)

    def _save(self):
        """
        Save the mirror list for (short-term) future  use.
        """
        self.mirrorfile.seek(0)
        self.mirrorfile.write(repr(self.mirrors))
        self.mirrorfile.truncate()
        self.mirrorfile.flush()

    def _port_of_mirror(self, mirror):
        if mirror.startswith('http://'):
            return 80
        if mirror.startswith('https://'):
            return 443
        if mirror.startswith('ftp://'):
            return 21

    def _rank_mirrors(self):
        """
        Sort the mirrors by speed, fastest being first

        This method is used by the YUM fastestmirror plugin
        """
        timed_mirrors = []
        import time, socket
        log.info('Searching fastest mirror')
        timeout = 1
        for mirror in self.mirrors:
            if not mirror.startswith('http'):
                log.debug('we currently can only handle http, got %s', mirror)
                continue
            port = self._port_of_mirror(mirror)
            mirror_hostname = urlparse.urlsplit(mirror).netloc
            time_before = time.time()
            try:
                sock = socket.create_connection((mirror_hostname, port), timeout)
                sock.close()
            except (IOError, socket.error, socket.timeout) as err:
                log.warning(str(err).strip() + ': ' + mirror)
                continue
            result = time.time() - time_before
            result_ms = int(1000 * result)
            log.info(str(result_ms).rjust(5) + 'ms: ' + mirror)
            timed_mirrors.append((result, mirror))
        if len(timed_mirrors) == 0:
            # We cannot reach any mirror directly, most likely firewall issue
           if 'http_proxy' not in os.environ:
               log.error('Could not reach any mirror directly and no proxy set')
               raise RuntimeError('no internet connection')
           log.info('Cannot time mirrors via proxy, using default order')
        else:
            timed_mirrors.sort()
            self.mirrors = [m[1] for m in timed_mirrors]
        log.info('Fastest mirror: ' + self.fastest)

    def _age(self):
        """
        Return the age of the cached mirror list in seconds
        """
        import time
        mtime = os.fstat(self.mirrorfd).st_mtime
        now = time.mktime(time.localtime())
        return now - mtime

    def _must_refresh(self):
        """
        Return whether we must download the mirror list.

        If and only if this method returns ``False`` is it admissible
        to use the cached mirror list.
        """
        if os.fstat(self.mirrorfd).st_size == 0:
            return True
        return self._age() > self.MAXAGE

    def _refresh(self):
        """
        Download and rank the mirror list.
        """
        log.info('Downloading the Sage mirror list')
        try:
            with contextlib.closing(urllib.urlopen(self.URL)) as f:
                mirror_list = f.read().decode("ascii")
        except IOError:
            log.critical('Downloading the mirror list failed')
        else:
            self.mirrors = self._load(mirror_list)
            self._rank_mirrors()
            self._save()

    def __iter__(self):
        """
        Iterate through the list of mirrors.

        This is the main entry point into the mirror list. Every
        script should just use this function to try mirrors in order
        of preference. This will not just yield the official mirrors,
        but also urls for packages that are currently being tested.
        """
        try:
            yield os.environ['SAGE_SERVER']
        except KeyError:
            pass
        for mirror in self.mirrors:
            yield mirror
        # If all else fails: Try the packages we host ourselves
        yield 'http://sagepad.org/'

    @property
    def fastest(self):
        return next(iter(self))
