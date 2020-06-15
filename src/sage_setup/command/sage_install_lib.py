from setuptools.command.install_lib import install_lib

class sage_install_lib(install_lib):

    def get_exclusions(self):

        exclusions = install_lib.get_exclusions(self)
        print("EXCLUSIONS: ", exclusions)
        # FIXME: Exclude files outside of the current distribution
        return exclusions

    pass
