"""
Make the required directories
"""

import os
from sage_pkg.config import config


def make_required_directories():
    required = [
        config.path.root,
        config.path.packages,
        config.path.download_cache,
        config.path.build,
        config.path.install,
        config.path.install_metadata,
    ]
    for path in required:
        if os.path.exists(path):
            if os.path.isdir(path):
                continue
            raise ValueError('required path name exists but is not directory: %s'.format(path))
        else:
            os.makedirs(path)
