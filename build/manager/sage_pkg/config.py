import os


config = None

def set_configuration(configuration_file):
    from sage_pkg.manager_config import ManagerConfig
    global config
    assert config is None # only allow setting the configuration once
    base_dir = os.path.dirname(configuration_file)
    config = ManagerConfig(configuration_file, base_dir=base_dir)


