"""
Constructor Object for Axioms

The purpose of this module is two-fold. First, it defines an object
that is available in the Sage global namespace as ``axioms``, allowing
you to easily construct axioms:: 

    sage: axioms.Finite()
    Finite

The other purpose is to define the sort order of the built-in axioms::

    sage: axioms.Finite() < axioms.Unital()
    True

To add new axioms, you would add them to
:meth:`AxiomFactory._init_axioms` at the desired place in the sort
order.
"""


class AxiomLazyImporter(object):

    def __init__(self, module_name, axiom_name, sort_key):
        self._module_name = module_name
        self._axiom_name = axiom_name
        self._sort_key = sort_key

    def __cmp__(self, other):
        return cmp(self._sort_key, other._sort_key)

    @property
    def __name__(self):
        return self._axiom_name

    def __repr__(self):
        return 'Lazy importer for the axiom ' + self._axiom_name

    def __get__(self, instance, owner):
        module = __import__(self._module_name, fromlist=[self._axiom_name])
        cls = module.__dict__[self._axiom_name]
        cls._sort_key = self._sort_key
        assert cls.__module__ == self._module_name
        assert cls.__name__ == self._axiom_name
        instance._add_class(cls)
        return cls


class AxiomFactory(object):
    """
    Library of all axioms
    """

    _instance = None

    def __init__(self):
        """
        Python constructor
    
        TESTS::

            sage: from sage.categories.axioms.factory import AxiomFactory
            sage: AxiomFactory()
            Traceback (most recent call last):
            ...
            AssertionError
        """
        assert self._instance is None   # ensure this is a singleton
        AxiomFactory._instance = self
        self._counter = 0
        self._lazy_axiom_sort_key = dict()

    def _init_axioms(self):
        """
        Initialize the list of built-in axioms.

        This is where the sort order of the built-in axioms is
        determined.
        """
        self._add_lazy('sage.categories.axioms.testing', 'Flying')
        self._add_lazy('sage.categories.axioms.testing', 'Blue')
        self._add_lazy('sage.categories.axioms.remaining', 'Facade') 
        self._add_lazy('sage.categories.axioms.finite', 'Finite') 
        self._add_lazy('sage.categories.axioms.finite', 'Infinite') 
        self._add_lazy('sage.categories.axioms.finite', 'FiniteDimensional')
        self._add_lazy('sage.categories.axioms.remaining', 'Connected')
        self._add_lazy('sage.categories.axioms.remaining', 'WithBasis')
        self._add_lazy('sage.categories.axioms.remaining', 'Irreducible')
        self._add_lazy('sage.categories.axioms.commutative', 'Commutative') 
        self._add_lazy('sage.categories.axioms.remaining', 'Associative')
        self._add_lazy('sage.categories.axioms.inverse', 'Inverse') 
        self._add_lazy('sage.categories.axioms.unital', 'Unital')
        self._add_lazy('sage.categories.axioms.remaining', 'Division')
        self._add_lazy('sage.categories.axioms.remaining', 'NoZeroDivisors')
        self._add_lazy('sage.categories.axioms.commutative', 'AdditiveCommutative') 
        self._add_lazy('sage.categories.axioms.remaining', 'AdditiveAssociative')
        self._add_lazy('sage.categories.axioms.inverse', 'AdditiveInverse')
        self._add_lazy('sage.categories.axioms.unital', 'AdditiveUnital')
        
    def _set_axiom(self, name, axiom, overwrite=False):
        """
        Add another axiom to the list
        """
        already_set = (name in AxiomFactory.__dict__.keys())
        if overwrite != already_set:
            raise ValueError('overwrite argument mismatch')
        setattr(AxiomFactory, name, axiom)

    def __repr__(self):
        return('Axiom factory object')

    def _add_lazy(self, module_name, axiom_name):
        """
        Lazily add axiom
        """
        sort_key = self._counter
        importer = AxiomLazyImporter(module_name, axiom_name, sort_key)
        self._set_axiom(axiom_name, importer)
        assert (module_name, axiom_name) not in self._lazy_axiom_sort_key
        self._lazy_axiom_sort_key[(module_name, axiom_name)] = sort_key
        self._counter += 1
        
    def _add_instance(self, axiom):
        """
        Add instance of axiom
        """
        axiom_type = type(axiom)
        self._add_class(axiom_type)

    def _add_class(self, axiom_class):
        """
        Add axiom class
        """
        module_name = axiom_class.__module__
        axiom_name = axiom_class.__name__
        try:
            sort_key = self._lazy_axiom_sort_key[(module_name, axiom_name)]
            overwrite = True
        except KeyError:
            sort_key = self._counter
            axiom_class._sort_key = sort_key
            self._counter += 1
            overwrite = False
        self._set_axiom(axiom_name, axiom_class, overwrite=overwrite)

    def register(self, new_axiom_class):
        from axiom import Axiom
        if not issubclass(new_axiom_class, Axiom):
            raise ValueError('argument must be a subclass of Axiom')
        if hasattr(new_axiom_class, '_sort_key'):
            raise ValueError('already registered')
        self._add_class(new_axiom_class)

    def _list_axioms(self):
        """
        Return the (sorted) current list of axiom holders.
        """
        from axiom import Axiom
        def is_axiom(ax):
            return isinstance(ax, AxiomLazyImporter) or \
                (isinstance(ax, type) and issubclass(ax, Axiom))
        return sorted([axiom for axiom in AxiomFactory.__dict__.values() if is_axiom(axiom)])

    def explain(self):
        print("Sorted list of known axioms:")
        for axiom in self._list_axioms():
            instantiated = ' ' if isinstance(axiom, AxiomLazyImporter) else '*'
            print('  {0} {1}'.format(instantiated, axiom.__name__))
        
    def _test_axioms(self):
        from sage.categories.axioms.axiom import Axiom
        for axiom_class in self._list_axioms():
            axiom = axiom_class()
            assert isinstance(axiom, Axiom)
        
    def deprecated_with_name(self, name):
        """
        TODO: remove this method
        """
        cls = getattr(self, name)
        return cls()

    def deprecated_all_names(self):
        """
        TODO: remove this method
        """
        return tuple(axiom.__name__ for axiom in self._list_axioms())


axioms = AxiomFactory()
axioms._init_axioms()

