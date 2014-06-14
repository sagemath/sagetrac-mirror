"""
Mixing for Sage environment variables
"""

import os
from sage_pkg.config import config


class SageEnvironmentMixin(object):

    def get_environment(self):
        """
        Set up the build-time environment for Sage packages
        """
        env = super(SageEnvironmentMixin, self).get_environment()
        env['SAGE_ROOT'] = SAGE_ROOT = config.path.root
        env['SAGE_LOCAL'] = SAGE_LOCAL = config.path.install
        env['SAGE_SHARE'] = os.path.join(SAGE_LOCAL, 'share')
        env['SAGE_EXTCODE'] = os.path.join(SAGE_LOCAL, 'share', 'sage', 'ext')
        env['PATH'] = ':'.join([
            os.path.join(SAGE_ROOT, 'src', 'bin'),
            os.path.join(SAGE_LOCAL, 'bin'),
            env['PATH'],
        ])
        env['PYTHONPATH'] = SAGE_LOCAL
        return env

