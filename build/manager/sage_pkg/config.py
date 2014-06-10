import os


config = None

def set_configuration(*configuration_files):
    from sage_pkg.config_yaml import SagePkgConfig
    global config
    assert config is None # only allow setting the configuration once
    config = SagePkgConfig(*configuration_files)


