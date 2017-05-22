from __future__ import absolute_import

import importlib
import logging


logger = logging.getLogger('docbuild')


def get_builder(name):
    """
    Returns an appropriate *Builder object for the document ``name``.
    DocBuilder and its subclasses do all the real work in building the
    documentation.

    This uses the Builder.match methods, called in order of Builder.priority to
    determine the Builder to use for a document.  If any of the match() methods
    raise an exception, it will be passed over to check if any other of the
    Builders match.  If no Builders are found the first exception encountered
    is raised, or a default exception if no exceptions were raised.
    """

    exceptions = []

    for builder_cls in sorted(get_builders(), key=lambda b: -b.priority):
        logger.debug("Matching document '%s' to builder '%s'" %
                     (name, builder_cls.__name__))
        try:
            builder = builder_cls.match(name)
        except Exception as e:
            exceptions.append((builder_cls, e))
            continue

        if builder is not None:
            logger.debug("Match of '%s' to builder '%s' succeeded" %
                         (name, builder_cls.__name__))
            break

    if builder is not None:
        for cls, exception in exceptions:
            logger.warning(
                'Exception raised by %s when trying to match it to '
                'document "%s": %s' % (cls.__name__, name, exception))

        return builder
    else:
        if exceptions:
            raise exceptions[0][1]
        else:
            raise NotImplementedError(
                "'%s' is not a recognized document. Type "
                "'sage --docbuild -D' for a list of documents, or "
                "'sage --docbuild --help' for more help." % name)


_all_builders = None
def get_builders():
    """
    Returns a list of all `Builder` classes.

    Currently all modules to search for `Builder` classes are submodules of
    this package, and are hard-coded below.  The `Builder` classes to return
    must be listed in the ``__all__`` of that module as well.
    """

    from .builders import Builder, AllBuilder

    global _all_builders

    if _all_builders is not None:
        return _all_builders

    # Just hard-coded for now
    submodules = [
        'docbuilder',
        'reference',
        'singlefile',
        'website'
    ]

    builders = [AllBuilder]

    for submodule in submodules:
        mod = importlib.import_module(__package__ + '.builders.' + submodule)
        for name in mod.__all__:
            obj = getattr(mod, name)
            if isinstance(obj, type) and issubclass(obj, Builder):
                builders.append(obj)

    _all_builders = builders

    return builders


def get_documents(default_lang=None):
    """
    Returns a list of document names the Sage documentation builder
    will accept as command-line arguments.
    """

    from .builders import AllBuilder

    all_b = AllBuilder()
    docs = all_b.get_all_documents(default_lang=default_lang)
    return docs


def get_formats():
    """
    Returns a list of output formats the Sage documentation builder
    will accept on the command-line.
    """

    from .builders.docbuilder import DocBuilder

    formats = DocBuilder.output_formats()

    # Place html and pdf as the first two formats, as they are the most common
    formats.remove('html')
    formats.remove('pdf')
    return ['html', 'pdf'] + formats
