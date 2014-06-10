import os

_HOME = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

class _ConfigDefault(object):

    HOME_DIR = _HOME

    APP_DIR = os.path.join(_HOME, 'sage_pkg')

    TEST_DIR = os.path.join(_HOME, 'sage_pkg', 'test')

    CONFIG_TEST_YAML = os.path.join(_HOME, 'test.yaml')



config = _ConfigDefault()


from sage_pkg.config import set_configuration
set_configuration(config.CONFIG_TEST_YAML)
