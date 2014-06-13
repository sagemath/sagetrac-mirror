"""
Mixin to download from the Sage mirror network
"""

from sage_pkg.config import config


class SageMirrorMixin(object):

    @property
    def tarball(self):
        return self._config.source.tarball

    def tarball_url_iter(self):
        try:
            mirrors = config.mirrors
        except AttributeError:
            raise ValueError('you need to define at least one mirror')
        for url in mirrors:
            url = url.rstrip('/')
            yield '/'.join(url, self.name, self.tarball)


    def download(self):
        pass
        
