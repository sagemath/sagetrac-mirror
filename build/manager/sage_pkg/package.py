"""
Package Data Reader
"""

import os
import yaml

def all_packages(package_dir, git):
    for name in os.listdir(package_dir):
        if not os.path.isdir(name):
            continue
        app_yaml = os.path.join(package_dir, name, 'app.yaml')
        if not os.path.isfile(name):
            continue
        yield _load(app_yaml, git)


def load_package(pkg_dir, git):
    app_yaml = os.path.join(package_dir, name, 'app.yaml')
    pkg_name = os.path.basename(pkg_dir)
    return _load(app_yaml, git)


def _load(app_yaml, git):
    config = yaml.load_all(app_yaml).next()
    return config
    
