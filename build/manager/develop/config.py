import os

_HOME = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

class _ConfigDefault(object):

    HOME_DIR = _HOME

    APP_DIR = os.path.join(_HOME, 'sage_pkg')

    TEST_DIR = os.path.join(_HOME, 'sage_pkg', 'test')




config = _ConfigDefault()


