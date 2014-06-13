import os


config = None

def set_configuration(*configuration_files):
    from sage_pkg.manager_config import ManagerConfig
    global config
    assert config is None # only allow setting the configuration once
    config = ManagerConfig(*configuration_files)


