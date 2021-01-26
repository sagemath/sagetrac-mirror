# -*- coding: utf-8 -*-
r"""
Base classes for features
"""

from sage.misc.lazy_string import lazy_string

class TrivialClasscallMetaClass(type):
    """
    A trivial version of :class:`ClasscallMetaclass` without Cython dependencies.
    """
    def __call__(cls, *args, **kwds):
        r"""
        This method implements ``cls(<some arguments>)``.
        """
        if hasattr(cls, '__classcall__'):
            return cls.__classcall__(cls, *args, **kwds)
        else:
            return type.__call__(cls, *args, **kwds)

_trivial_unique_representation_cache = dict()

class TrivialUniqueRepresentation(metaclass=TrivialClasscallMetaClass):
    r"""
    A trivial version of :class:`UniqueRepresentation` without Cython dependencies.
    """

    @staticmethod
    def __classcall__(cls, *args, **options):
        r"""
        Construct a new object of this class or reuse an existing one.
        """
        key = (cls, tuple(args), frozenset(options.items()))
        cached = _trivial_unique_representation_cache.get(key, None)
        if cached is None:
            cached = _trivial_unique_representation_cache[key] = type.__call__(cls, *args, **options)
        return cached

class Feature(TrivialUniqueRepresentation):
    r"""
    A feature of the runtime environment

    Overwrite :meth:`_is_present` to add feature checks.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: GapPackage("grape", spkg="gap_packages")  # indirect doctest
        Feature('GAP package grape')

    For efficiency, features are unique::

        sage: GapPackage("grape") is GapPackage("grape")
        True
    """
    def __init__(self, name, spkg=None, url=None):
        r"""
        TESTS::

            sage: from sage.features import Feature
            sage: from sage.features.gap import GapPackage
            sage: isinstance(GapPackage("grape", spkg="gap_packages"), Feature)  # indirect doctest
            True
        """
        self.name = name
        self.spkg = spkg
        self.url = url
        self._cache_is_present = None
        self._cache_resolution = None

    def is_present(self):
        r"""
        Return whether the feature is present.

        OUTPUT:

        A :class:`FeatureTestResult` which can be used as a boolean and
        contains additional information about the feature test.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape", spkg="gap_packages").is_present()  # optional: gap_packages
            FeatureTestResult('GAP package grape', True)
            sage: GapPackage("NOT_A_PACKAGE", spkg="gap_packages").is_present()
            FeatureTestResult('GAP package NOT_A_PACKAGE', False)

        The result is cached::

            sage: from sage.features import Feature
            sage: class TestFeature(Feature):
            ....:     def _is_present(self):
            ....:         print("checking presence")
            ....:         return True
            sage: TestFeature("test").is_present()
            checking presence
            FeatureTestResult('test', True)
            sage: TestFeature("test").is_present()
            FeatureTestResult('test', True)
            sage: TestFeature("other").is_present()
            checking presence
            FeatureTestResult('other', True)
            sage: TestFeature("other").is_present()
            FeatureTestResult('other', True)
        """
        # We do not use @cached_method here because we wish to use
        # Feature early in the build system of sagelib.
        if self._cache_is_present is None:
            res = self._is_present()
            if not isinstance(res, FeatureTestResult):
                res = FeatureTestResult(self, res)
            self._cache_is_present = res
        return self._cache_is_present

    def _is_present(self):
        r"""
        Override this in a derived class to implement the feature check.

        This should return either an instance of
        :class:`FeatureTestResult` or a boolean.
        """
        raise NotImplementedError("_is_present not implemented for feature {!r}".format(self.name))

    def require(self):
        r"""
        Raise a :class:`FeatureNotPresentError` if the feature is not present.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("ve1EeThu").require()
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package ve1EeThu is not available.
            `TestPackageAvailability("ve1EeThu")` evaluated to `fail` in GAP.
        """
        presence = self.is_present()
        if not presence:
            raise FeatureNotPresentError(self, presence.reason, presence.resolution)

    def __repr__(self):
        r"""
        Return a printable representation of this object.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("grape")  # indirect doctest
            Feature('GAP package grape')
        """
        return 'Feature({name!r})'.format(name=self.name)

    def resolution(self):
        r"""
        Return a suggestion on how to make :meth:`is_present` pass if it did not
        pass.

        OUTPUT:

        A string, a lazy string, or ``None``.  The default implementation always
        returns a lazy string.

        EXAMPLES::

            sage: from sage.features import Executable
            sage: Executable(name="CSDP", spkg="csdp", executable="theta", url="http://github.org/dimpase/csdp").resolution()  # optional - sage_spkg
            l'...To install CSDP...you can try to run...sage -i csdp...Further installation instructions might be available at http://github.org/dimpase/csdp.'
        """
        def find_resolution():
            if self._cache_resolution is not None:
                return self._cache_resolution
            lines = []
            if self.spkg:
                from .pkg_systems import package_systems
                for ps in package_systems():
                    lines.append(ps.spkg_installation_hint(self.spkg, feature=self.name))
            if self.url:
                lines.append("Further installation instructions might be available at {url}.".format(url=self.url))
            self._cache_resolution = "\n".join(lines)
            return self._cache_resolution

        return lazy_string(find_resolution)

class FeatureNotPresentError(RuntimeError):
    r"""
    A missing feature error.

    EXAMPLES::

        sage: from sage.features import Feature, FeatureTestResult
        sage: class Missing(Feature):
        ....:     def _is_present(self):
        ....:         return False

        sage: Missing(name="missing").require()
        Traceback (most recent call last):
        ...
        FeatureNotPresentError: missing is not available.
    """
    def __init__(self, feature, reason=None, resolution=None):
        self.feature = feature
        self.reason = reason
        self.resolution = resolution or feature.resolution()

    def __str__(self):
        r"""
        Return the error message.

        EXAMPLES::

            sage: from sage.features.gap import GapPackage
            sage: GapPackage("gapZuHoh8Uu").require()  # indirect doctest
            Traceback (most recent call last):
            ...
            FeatureNotPresentError: GAP package gapZuHoh8Uu is not available.
            `TestPackageAvailability("gapZuHoh8Uu")` evaluated to `fail` in GAP.
        """
        lines = ["{feature} is not available.".format(feature=self.feature.name)]
        if self.reason:
            lines.append(self.reason)
        if self.resolution:
            lines.append(str(self.resolution))
        return "\n".join(lines)


class FeatureTestResult(object):
    r"""
    The result of a :meth:`Feature.is_present` call.

    Behaves like a boolean with some extra data which may explain why a feature
    is not present and how this may be resolved.

    EXAMPLES::

        sage: from sage.features.gap import GapPackage
        sage: presence = GapPackage("NOT_A_PACKAGE").is_present(); presence  # indirect doctest
        FeatureTestResult('GAP package NOT_A_PACKAGE', False)
        sage: bool(presence)
        False

    Explanatory messages might be available as ``reason`` and
    ``resolution``::

        sage: presence.reason
        '`TestPackageAvailability("NOT_A_PACKAGE")` evaluated to `fail` in GAP.'
        sage: bool(presence.resolution)
        False

    If a feature is not present, ``resolution`` defaults to
    ``feature.resolution()`` if this is defined. If you do not want to use this
    default you need explicitly set ``resolution`` to a string::

        sage: from sage.features import FeatureTestResult
        sage: package = GapPackage("NOT_A_PACKAGE", spkg="no_package")
        sage: str(FeatureTestResult(package, True).resolution)  # optional - sage_spkg
        '...To install GAP package NOT_A_PACKAGE...you can try to run...sage -i no_package...'
        sage: str(FeatureTestResult(package, False).resolution) # optional - sage_spkg
        '...To install GAP package NOT_A_PACKAGE...you can try to run...sage -i no_package...'
        sage: FeatureTestResult(package, False, resolution="rtm").resolution
        'rtm'
    """
    def __init__(self, feature, is_present, reason=None, resolution=None):
        r"""
        TESTS::

            sage: from sage.features import Executable, FeatureTestResult
            sage: isinstance(Executable(name="sh", executable="sh").is_present(), FeatureTestResult)
            True
        """
        self.feature = feature
        self.is_present = is_present
        self.reason = reason
        self.resolution = resolution or feature.resolution()

    def __bool__(self):
        r"""
        Whether the tested :class:`Feature` is present.

        TESTS::

            sage: from sage.features import Feature, FeatureTestResult
            sage: bool(FeatureTestResult(Feature("SomePresentFeature"), True))  # indirect doctest
            True
            sage: bool(FeatureTestResult(Feature("SomeMissingFeature"), False))
            False
        """
        return bool(self.is_present)

    __nonzero__ = __bool__

    def __repr__(self):
        r"""
        TESTS::

            sage: from sage.features import Feature, FeatureTestResult
            sage: FeatureTestResult(Feature("SomePresentFeature"), True)  # indirect doctest
            FeatureTestResult('SomePresentFeature', True)
        """
        return "FeatureTestResult({feature!r}, {is_present!r})".format(feature=self.feature.name, is_present=self.is_present)
