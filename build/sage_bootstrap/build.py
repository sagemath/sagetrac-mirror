from __future__ import absolute_import

def build(targets):
    from sage_conf import SAGE_ROOT
    from os import system
    command = 'cd {} && make {}'.format(SAGE_ROOT, " ".join(targets))
    print("Running {}".format(command))
    retval = system(command)
    if retval != 0:
        raise RuntimeError("Failed with status {}".format(retval))
