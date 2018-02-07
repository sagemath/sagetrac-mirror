import os
import sys
import tokenize
import types

import importlib.machinery

from importlib.abc import PathEntryFinder
from importlib.machinery import ExtensionFileLoader, FileFinder
from importlib.util import spec_from_file_location


class RestrictedFileFinder(FileFinder):
    """
    Like FileFinder, but its path hook is restricted to a given directory and
    its subdirectories.
    """

    @classmethod
    def path_hook(cls, basepath, *loader_details):
        orig_hook = super().path_hook(*loader_details)

        if basepath is None:
            return orig_hook

        basepath = os.path.abspath(basepath)

        def path_hook_for_RestrictedFileFinder(path):
            path = os.path.abspath(path)
            if os.path.commonpath([path, basepath]) != basepath:
                raise ImportError(
                    'only directories under {} are supported'.format(
                        basepath), path=path)

            return orig_hook(path)

        return path_hook_for_RestrictedFileFinder


@PathEntryFinder.register
class MetaFileFinder:
    """
    A 'middleware', if you will, between the PathFinder sys.meta_path hook,
    and sys.path_hooks hooks--particularly FileFinder.

    The hook returned by FileFinder.path_hook is rather 'promiscuous' in that
    it will handle *any* directory.  So if one wants to insert another
    FileFinder.path_hook into sys.path_hooks, that will totally take over
    importing for any directory, and previous path hooks will be ignored.

    This class provides its own sys.path_hooks hook as follows: If inserted
    on sys.path_hooks (it should be inserted early so that it can supersede
    anything else).  Its find_spec method then calls each hook on
    sys.path_hooks after itself and, for each hook that can handle the given
    sys.path entry, it calls the hook to create a finder, and calls that
    finder's find_spec.  So each sys.path_hooks entry is tried until a spec is
    found or all finders are exhausted.
    """

    class hook:
        """
        Use this little internal class rather than a function with a closure
        or a classmethod or anything like that so that it's easier to
        identify our hook and skip over it while processing sys.path_hooks.
        """

        def __init__(self, basepath=None):
            self.basepath = os.path.abspath(basepath)

        def __call__(self, path):
            if not os.path.isdir(path):
                raise ImportError('only directories are supported', path=path)
            elif not self.handles(path):
                raise ImportError(
                    'only directories under {} are supported'.format(
                        self.basepath), path=path)

            return MetaFileFinder(path)

        def handles(self, path):
            """
            Return whether this hook will handle the given path, depending on
            what its basepath is.
            """

            path = os.path.abspath(path)

            return (self.basepath is None or
                    os.path.commonpath([self.basepath, path]) == self.basepath)

    def __init__(self, path):
        self.path = path
        self._finder_cache = {}

    def __repr__(self):
        return '{}({!r})'.format(self.__class__.__name__, self.path)

    def find_spec(self, fullname, target=None):
        if not sys.path_hooks:
            return None

        last = len(sys.path_hooks) - 1

        for idx, hook in enumerate(sys.path_hooks):
            if isinstance(hook, self.__class__.hook):
                continue

            finder = None
            try:
                if hook in self._finder_cache:
                    finder = self._finder_cache[hook]
                    if finder is None:
                        # We've tried this finder before and got an ImportError
                        continue
            except TypeError:
                # The hook is unhashable
                pass

            if finder is None:
                try:
                    finder = hook(self.path)
                except ImportError:
                    pass

            try:
                self._finder_cache[hook] = finder
            except TypeError:
                # The hook is unhashable for some reason so we don't bother
                # caching it
                pass

            if finder is not None:
                spec = finder.find_spec(fullname, target)
                if (spec is not None and
                        (spec.loader is not None or idx == last)):
                    # If no __init__.<suffix> was found by any Finder,
                    # we may be importing a namespace package (which
                    # FileFinder.find_spec returns in this case).  But we
                    # only want to return the namespace ModuleSpec if we've
                    # exhausted every other finder first.
                    return spec

        # Module spec not found through any of the finders
        return None

    def invalidate_caches(self):
        for finder in self._finder_cache.values():
            finder.invalidate_caches()

    @classmethod
    def install(cls, basepath=None):
        """
        Install the MetaFileFinder in the front sys.path_hooks, so that
        it can support any existing sys.path_hooks and any that might
        be appended later.

        If given, only support paths under and including basepath.  In this
        case it's not necessary to invalidate the entire
        sys.path_importer_cache, but only any existing entries under basepath.
        """

        if basepath is not None:
            basepath = os.path.abspath(basepath)

        hook = cls.hook(basepath)
        sys.path_hooks.insert(0, hook)
        if basepath is None:
            sys.path_importer_cache.clear()
        else:
            for path in list(sys.path_importer_cache):
                if hook.handles(path):
                    del sys.path_importer_cache[path]


class CythonModule(types.ModuleType):
    """
    Custom module type for modules compiled from Cython sources, which
    allows the original Cython source (if available) to be accessed
    through its __cython_source__ attribute.
    """

    @property
    def __cython_source__(self):
        if not (hasattr(self, '__loader__') and
                isinstance(self.__loader__, CompiledPyxFileLoader)):
            return None

        return self.__loader__.get_source(self.__name__)


class CompiledPyxFileLoader(ExtensionFileLoader):
    """
    Thin wrapper around ExtensionFileLoader which, given the path to a
    .pyx file, loads its associated extension module, assuming the extension
    module is already compiled and under the same path.

    In other words, this simply replaces the .pyx in the path with the
    appropriate extension module suffix.  It also adds a __cython_source__
    attribute to loaded modules, containing the Cython source code of the
    module.
    """

    SUFFIXES = ('.pyx',)

    def __init__(self, name, path):
        pyx_path = path
        basename, _ = os.path.splitext(os.path.basename(pyx_path))
        dirname = os.path.dirname(pyx_path)
        for suffix in importlib.machinery.EXTENSION_SUFFIXES:
            path = os.path.join(dirname, basename + suffix)
            if os.path.isfile(path):
                break
        else:
            raise ImportError(
                'could not find the compiled extension module associated '
                'with the Cython module {}'.format(pyx_path), path=pyx_path)

        super().__init__(name, path)
        self.pyx_path = pyx_path

    def create_module(self, spec):
        """This hook allows us to customize the module class."""

        # Note: This is passed the ModuleSpec for the .pyx file, but
        # the parent class's create_module really expects the spec for the
        # actual extension module, which we thus create here, and ignore
        # the passed in spec
        spec = spec_from_file_location(self.name, self.path)
        module = super().create_module(spec)
        module.__class__ = CythonModule
        return module

    def get_source(self, fullname):
        """Retrieve the Cython source from the .pyx file."""

        if not os.path.isfile(self.pyx_path):
            # The file should not have disappeared.  Maybe issue a warning
            # here?
            return None

        with tokenize.open(self.pyx_path) as f:
            return f.read()

    @classmethod
    def install(cls, basepath=None):
        """
        Install a FileFinder hook to load .pyx files with the
        CompiledPyxFileLoader.


        Place it early on the path so that it overrides the default FileFinder
        hook for extension modules.
        """

        hook = RestrictedFileFinder.path_hook(basepath, (cls, cls.SUFFIXES))
        sys.path_hooks.insert(0, hook)


# Install the hook for .pyx files first, so that the MetaFileFinder hook
# precedes it (this isn't strictly necessary but is still more logical)
CompiledPyxFileLoader.install(os.path.dirname(__file__))

# Install the MetaFileFinder for the sage package only, for now.
# However, this may also be useful for packages that use Sage so
# consider extending this in the future...
MetaFileFinder.install(os.path.dirname(__file__))
