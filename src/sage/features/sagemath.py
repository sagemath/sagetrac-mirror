r"""
Check for SageMath Python modules
"""
from . import PythonModule
from .join_feature import JoinFeature


class sage__combinat(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.combinat``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__combinat
        sage: sage__combinat().is_present()  # optional - sage.combinat
        FeatureTestResult('sage.combinat', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__combinat
            sage: isinstance(sage__combinat(), sage__combinat)
            True
        """
        # sage.combinat will be a namespace package.
        # Testing whether sage.combinat itself can be imported is meaningless.
        # Hence, we test a Python module within the package.
        JoinFeature.__init__(self, 'sage.combinat',
                             [PythonModule('sage.combinat.combinations')])


class sage__geometry__polyhedron(PythonModule):

    def __init__(self):
        PythonModule.__init__(self, 'sage.geometry.polyhedron')


class sage__graphs(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.graphs``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__graphs
        sage: sage__graphs().is_present()  # optional - sage.graphs
        FeatureTestResult('sage.graphs', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__graphs
            sage: isinstance(sage__graphs(), sage__graphs)
            True
        """
        JoinFeature.__init__(self, 'sage.graphs',
                             [PythonModule('sage.graphs.graph')])


class sage__plot(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.plot``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__plot
        sage: sage__plot().is_present()  # optional - sage.plot
        FeatureTestResult('sage.plot', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__plot
            sage: isinstance(sage__plot(), sage__plot)
            True
        """
        JoinFeature.__init__(self, 'sage.plot',
                             [PythonModule('sage.plot.plot')])


class sage__rings__number_field(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.rings.number_field``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__number_field
        sage: sage__rings__number_field().is_present()  # optional - sage.rings.number_field
        FeatureTestResult('sage.rings.number_field', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__number_field
            sage: isinstance(sage__rings__number_field(), sage__rings__number_field)
            True
        """
        JoinFeature.__init__(self, 'sage.rings.number_field',
                             [PythonModule('sage.rings.number_field.number_field_element')])


class sage__rings__real_double(PythonModule):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.rings.real_double``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__rings__real_double
        sage: sage__rings__real_double().is_present()  # optional - sage.rings.real_double
        FeatureTestResult('sage.rings.real_double', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__rings__real_double
            sage: isinstance(sage__rings__real_double(), sage__rings__real_double)
            True
        """
        PythonModule.__init__(self, 'sage.rings.real_double')


class sage__symbolic(JoinFeature):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``sage.symbolic``.

    EXAMPLES::

        sage: from sage.features.sagemath import sage__symbolic
        sage: sage__symbolic().is_present()  # optional - sage.symbolic
        FeatureTestResult('sage.symbolic', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.sagemath import sage__symbolic
            sage: isinstance(sage__symbolic(), sage__symbolic)
            True
        """
        JoinFeature.__init__(self, 'sage.symbolic',
                             [PythonModule('sage.symbolic.expression')],
                             spkg="sagemath_symbolics")


def sage_features():
    """
    Return features corresponding to parts of the Sage library.

    These tags are named after Python packages/modules (e.g., :mod:`~sage.symbolic`),
    not distribution packages (``sagemath-symbolics``).

    This design is motivated by a separation of concerns: The author of a module that depends
    on some functionality provided by a Python module usually already knows the
    name of the Python module, so we do not want to force the author to also
    know about the distribution package that provides the Python module.

    Instead, we associate distribution packages to Python modules in
    :mod:`sage.features.sagemath` via the ``spkg`` parameter of :class:`Feature`.

    EXAMPLES::

        sage: from sage.features.sagemath import sage_features
        sage: list(sage_features())  # random
        [Feature('sage.graphs'),
         Feature('sage.plot'),
         Feature('sage.rings.number_field'),
         Feature('sage.rings.real_double')]
    """
    for feature in [sage__combinat(),
                    sage__geometry__polyhedron(),
                    sage__graphs(),
                    sage__plot(),
                    sage__rings__number_field(),
                    sage__rings__real_double(),
                    sage__symbolic()]:
        if feature.is_present():
            yield feature
