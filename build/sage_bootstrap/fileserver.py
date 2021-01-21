# -*- coding: utf-8 -*-
"""
Interface to the Sage fileserver
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



import os
import subprocess



class FileServer(object):

    def __init__(self):
        self.user = 'sagemath'
        self.hostname = 'fileserver.sagemath.org'

    def upstream_directory(self, package):
        """
        Return the directory where the tarball resides on the server
        """
        return os.path.join(
            '/', 'data', 'files', 'spkg', 'upstream', package.name,
        )

    def upload(self, package):
        """
        Upload the current tarball of package
        """
        if not package.tarball.is_distributable():
            raise ValueError('Tarball of {} is marked as not distributable'.format(package))
        subprocess.check_call([
            'ssh', 'sagemath@fileserver.sagemath.org',
            'mkdir -p {0} && touch {0}/index.html'.format(self.upstream_directory(package))
        ])
        subprocess.check_call([
            'rsync', '-av', '--checksum', '-e', 'ssh -l sagemath',
            package.tarball.upstream_fqn,
            'fileserver.sagemath.org:{0}'.format(self.upstream_directory(package))
        ])

    def publish(self):
        """
        Publish the files
        """
        subprocess.check_call([
            'ssh', 'sagemath@fileserver.sagemath.org', './publish-files.sh'
        ])
