
import sys
from subprocess import check_call, CalledProcessError



def run_doctests(doctest_dirs, filename=None):
    from develop.tester import DocTester
    dt = DocTester()
    if filename is None:
        dt.add_dir(*doctest_dirs)
    else:
        dt.add_file(filename)
    dt.run()


def run_unittests(unittest_dirs, filename=None):
    from develop.tester import UnitTester
    ut = UnitTester()
    if filename is None:
        ut.add_dir(*unittest_dirs)
    else:
        ut.add_file(filename)
    ut.run()


def run_all(filename=None):
    from subprocess import check_call
    file_arg = [] if filename is None else [filename]
    check_call([sys.argv[0], '--doc']  + file_arg)
    check_call([sys.argv[0], '--unit'] + file_arg)
