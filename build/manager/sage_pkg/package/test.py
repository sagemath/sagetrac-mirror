"""
Builder for tests
"""


from .base import PackageBase

class TestPackage(PackageBase):

    def download(self):
        print('{self.name}: downloading'.format(self=self))

    def unpack(self):
        print('{self.name}: unpacking'.format(self=self))

    def prepare(self):
        print('{self.name}: preparing'.format(self=self))

    def configure(self):
        print('{self.name}: configuring'.format(self=self))

    def build(self):
        print('{self.name}: building'.format(self=self))

    def check(self):
        print('{self.name}: checking'.format(self=self))

    def install(self):
        print('{self.name}: installing'.format(self=self))


