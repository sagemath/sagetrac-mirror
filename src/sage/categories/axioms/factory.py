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


class AxiomHolder(object):

    def __init__(self, axiom, sort_key):
        self._axiom = axiom
        self._sort_key = sort_key

    def __cmp__(self, other):
        return cmp(self._sort_key, other._sort_key)

    def __call__(self):
        return self._axiom

    def __str__(self):
        return str(self._axiom)

    def __repr__(self):
        return 'Method returning the axiom ' + str(self)


class AxiomLazyImporter(AxiomHolder):

    def __init__(self, module_name, axiom_name, sort_key):
        self._module_name = module_name
        super(AxiomLazyImporter, self).__init__(axiom_name, sort_key)

    def __call__(self):
        axiom_name = self._axiom
        module = __import__(self._module_name, fromlist=[axiom_name])
        cls = module.__dict__[axiom_name]
        instance = cls()
        return instance


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
        
    def _set_axiom(self, name, axiom_holder, overwrite=False):
        """
        Add another axiom to the list
        """
        if overwrite != hasattr(self, name):
            raise ValueError('overwrite argument mismatch')
        setattr(self, name, axiom_holder)

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
        return sort_key
        
    def _add_instance(self, axiom):
        """
        Add instance of axiom

        Note that this method is called by every axiom constructor, in
        particular during the resolution of the lazily loaded axioms
        from :meth:`_add_lazy`.
        """
        axiom_type = type(axiom)
        module_name = axiom_type.__module__
        axiom_name = axiom_type.__name__
        try:
            sort_key = self._lazy_axiom_sort_key[(module_name, axiom_name)]
            overwrite = True
        except KeyError:
            sort_key = self._counter
            self._counter += 1
            overwrite = False
        holder = AxiomHolder(axiom, sort_key)
        self._set_axiom(axiom_name, holder, overwrite=overwrite)
        return sort_key

    def _list_axioms(self):
        """
        Return the (sorted) current list of axiom holders.
        """
        return sorted([axiom for axiom in self.__dict__.values()
                       if isinstance(axiom, AxiomHolder)])

    def explain(self):
        print("Sorted list of known axioms:")
        for axiom in self._list_axioms():
            instantiated = ' ' if isinstance(axiom, AxiomLazyImporter) else '*'
            print('  {0} {1}'.format(instantiated, str(axiom)))
        
    def _test_axioms(self):
        from sage.categories.axioms.axiom import Axiom
        for axiom_class in self._list_axioms():
            axiom = axiom_class()
            assert isinstance(axiom, Axiom)
        
    def deprecated_with_name(self, name):
        """
        TODO: remove this method
        """
        for axiom_holder in self._list_axioms():
            if str(axiom_holder) == name:
                return axiom_holder()
        raise ValueError('no such axiom: '+str(name))

    def deprecated_all_names(self):
        """
        TODO: remove this method
        """
        return map(str, self._list_axioms())


axioms = AxiomFactory()
axioms._init_axioms()
