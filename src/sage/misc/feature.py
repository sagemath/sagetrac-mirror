# -*- coding: utf-8 -*-
r"""
Testing for Features of the Environment at Runtime

A computation can require a certain package to be installed in the runtime
environment. Abstractly such a package describes a :class`Feature` which can be
tested for at runtime. It can be of various kinds, most prominently a
:class:`Executable` in the PATH or an additional package for some installed
system such as a :class:`GapPackage`.

AUTHORS:

- Julian Rüth (2016-04-07): Initial version

EXAMPLES:

There are special features already defined for commonly used package types such
as packages in GAP::

    sage: from sage.misc.feature import GapPackage, Executable
    sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
    True

To test for the existence of a binary, one can rely on a :class:`Executable`
feature::

    sage: Executable(name="sh", executable="sh").is_present()
    FeatureTestResult('sh', True)

Features try to produce a helpful error message if a feature is not present::

    sage: Executable(name="random", executable="randomOochoz6x", spkg="random", url="http://rand.om").require()
    Traceback (most recent call last):
    ...
    FeatureNotPresentError: random is not available.
    No `randomOochoz6x` found on PATH.
    To install random you can try to run `sage -i random`.
    Further installation instructions might be available at http://rand.om.

"""
from sage.misc.cachefunc import cached_method

class FeatureNotPresentError(RuntimeError):
    r"""
    Indicates a missing feature.

    EXAMPLES::

        sage: from sage.misc.feature import Feature, FeatureTestResult
        sage: class Missing(Feature):
        ....:     def is_present(self): return FeatureTestResult(self, False)

        sage: Missing(name="missing").require() # indirect doctest
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: missing is not available.

    """
    def __init__(self, feature, explanation=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.misc.feature import Feature, FeatureNotPresentError, FeatureTestResult
            sage: class Missing(Feature):
            ....:     def is_present(self): return FeatureTestResult(self, False)

            sage: try:
            ....:     Missing(name="missing").require() # indirect doctest
            ....: except FeatureNotPresentError: pass

        """
        self.feature = feature
        self.explanation = explanation
        self.resolution = resolution

    def __str__(self):
        r"""
        An error message.

        EXAMPLES::

            sage: from sage.misc.feature import GapPackage
            sage: GapPackage("gapZuHoh8Uu").require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package gapZuHoh8Uu is not available.
            `TestPackageAvailability("gapZuHoh8Uu")` evaluated to `fail` in GAP.

        """
        return "\n".join(filter(None,(
            "%s is not available."%self.feature.name,
            self.explanation,
            self.resolution
            )))

class FeatureTestResult(object):
    r"""
    The result of a :method:`Feature.is_present` call.

    Behaves like a boolean with some extra data which may explain why a feature
    is not present and how this may be resolved.

    EXAMPLES::

        sage: from sage.misc.feature import GapPackage
        sage: presence = GapPackage("NOT_A_PACKAGE").is_present(); presence # indirect doctest
        FeatureTestResult('GAP package NOT_A_PACKAGE', False)
        sage: bool(presence)
        False

    Explanatory messages might be available as ``explanation`` and
    ``resolution``::

        sage: presence.explanation
        '`TestPackageAvailability("NOT_A_PACKAGE")` evaluated to `fail` in GAP.'
        sage: presence.resolution is None
        True

    """
    def __init__(self, feature, is_present, explanation=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.misc.feature import Executable, FeatureTestResult
            sage: isinstance(Executable(name="sh", executable="sh").is_present(), FeatureTestResult)
            True

        """
        self.feature = feature
        self.is_present = is_present
        self.explanation = explanation
        self.resolution = resolution

    def __nonzero__(self):
        r"""
        Whether the tested :class:`Feature` is present.

        TESTS::

            sage: from sage.misc.feature import Feature, FeatureTestResult
            sage: bool(FeatureTestResult(Feature("SomePresentFeature"), True)) # indirect doctest
            True
            sage: bool(FeatureTestResult(Feature("SomeMissingFeature"), False))
            False

        """
        return bool(self.is_present)

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.misc.feature import Feature, FeatureTestResult
            sage: FeatureTestResult(Feature("SomePresentFeature"), True) # indirect doctest
            FeatureTestResult('SomePresentFeature', True)

        """
        return "FeatureTestResult(%s, %s)"%(repr(self.feature.name), repr(self.is_present))

class Feature(object):
    r"""
    A feature of the runtime environment

    Overwrite :meth:`is_present` to add feature checks.

    EXAMPLES::

        sage: from sage.misc.feature import GapPackage
        sage: GapPackage("grape", spkg="gap_packages") # indirect doctest
        Feature("GAP package grape")

    """
    VERBOSE_LEVEL = 50

    def __init__(self, name, spkg = None):
        r"""
        TESTS::

            sage: from sage.misc.feature import GapPackage, Feature
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), Feature) # indirect doctest
            True

        """
        self.name = name;
        self.spkg = spkg

    def is_present(self, explain=False):
        r"""
        Whether the feature is present.

        OUTPUT:

        A :class:`FeatureTestResult` which can be used as a boolean and
        contains additional information about the feature test.

        EXAMPLES::

            sage: from sage.misc.feature import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
            True

        """
        return FeatureTestResult(self, True)

    def require(self):
        r"""
        Raise a :class:`FeatureNotPresentError` if the feature is not present.

        EXAMPLES::

            sage: from sage.misc.feature import GapPackage
            sage: GapPackage("ve1EeThu").require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package ve1EeThu is not available.
            `TestPackageAvailability("ve1EeThu")` evaluated to `fail` in GAP.

        """
        presence = self.is_present()
        if not presence:
            raise FeatureNotPresentError(self, presence.explanation, presence.resolution)

    def __repr__(self):
        r"""
        A printable representation of this object.

        EXAMPLES::

            sage: from sage.misc.feature import GapPackage
            sage: GapPackage("grape") # indirect doctest
            Feature("GAP package grape")

        """
        return 'Feature("%s")'%(self.name,)

class GapPackage(Feature):
    r"""
    A feature describing the presence of a GAP package.

    EXMAPLES::

        sage: from sage.misc.feature import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")
        Feature("GAP package grape")

    """
    def __init__(self, package, spkg=None):
        Feature.__init__(self, "GAP package %s"%package)
        self.package = package
        self.spkg = spkg

    @cached_method
    def is_present(self):
        r"""
        Whether the package is available in GAP.

        This does not check whether this package is functional.

        EXAMPLES::

            sage: from sage.misc.feature import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present() # optional: gap_packages
            True

        """
        from sage.libs.gap.libgap import libgap
        command = 'TestPackageAvailability("%s")'%self.package
        presence = libgap.eval(command)
        if presence:
            return FeatureTestResult(self, True, explanation = "`%s` evaluated to `%s` in GAP."%(command, presence))
        else:
            resolution = "You can try to install this package with sage -i %s."%(self.spkg,) if self.spkg else None
            return FeatureTestResult(self, False, explanation = "`%s` evaluated to `%s` in GAP."%(command, presence), resolution = resolution)

class Executable(Feature):
    r"""
    A feature describing the presence of an executable in the PATH.

    .. NOTE::

        Overwrite :meth:`is_functional` if you also want to check if the executable
        shows proper behaviour.
        Calls to :meth:`is_present` are cached. You might want to cache the
        :class:`Executable` object to prevent unnecessary calls to the
        excutable.

    EXAMPLES::

        sage: from sage.misc.feature import Executable
        sage: Executable(name="sh", executable="sh").is_present()
        FeatureTestResult('sh', True)

    """
    def __init__(self, name, executable, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.misc.feature import Executable
            sage: isinstance(Executable(name="sh", executable="sh"), Executable)
            True

        """
        Feature.__init__(self, name=name)
        self.executable = executable
        self.spkg = spkg
        self.url = url

    @cached_method
    def is_present(self):
        r"""
        Test whether the executable is on the current PATH and functional.

        EXAMPLES::

            sage: from sage.misc.feature import Executable
            sage: Executable(name="sh", executable="sh").is_present()
            FeatureTestResult('sh', True)

        """
        from distutils.spawn import find_executable
        if find_executable(self.executable) is None:
            return FeatureTestResult(self, False, "No `%s` found on PATH."%(self.executable,), self.resolution())
        return self.is_functional()

    def is_functional(self):
        r"""
        Run as part of :meth:`is_present` to decide whether the executable on
        PATH is functional.

        EXAMPLES:

        Returns ``True`` unless explicitly overwritten::

            sage: from sage.misc.feature import Executable
            sage: Executable(name="sh", executable="sh").is_functional()
            FeatureTestResult('sh', True)

        """
        return FeatureTestResult(self, True)

    def resolution(self):
        r"""
        EXAMPLES::

            sage: from sage.misc.feature import Executable
            sage: Executable(name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp").resolution()
            'To install CSDP you can try to run `sage -i csdp`.\nFurther installation instructions might be available at http://github.org/dimpase/csdp.'

        """
        return "\n".join(filter(None,[
            "To install %s you can try to run `sage -i %s`."%(self.name, self.spkg) if self.spkg else "",
            "Further installation instructions might be available at %s."%(self.url,) if self.url else ""]))
