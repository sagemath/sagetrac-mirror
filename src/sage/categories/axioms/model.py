"""
Category/Axiom Model

The :class:`CategoryAxiomModel` is the central factory object for
categories with axioms. This is where relations are used to deduce the
simplest presentation of a category with axioms.

You should always use the ``category_axiom_model`` instance
"""


from sage.structure.sage_object import SageObject
from sage.misc.classcall_metaclass import typecall


DEBUG = True


class ClassHolder(object):

    def __init__(self, category_class):
        self.cls = category_class
        self.base = category_class.__base_category_class__
        self.axioms = tuple(sorted([ax() for ax in category_class.__axiom_classes__]))

    def __repr__(self):
        return '{0} implements {1} with {2}'.format(
            self.cls.__name__, self.base.__name__, self.axioms)


class BaseHolder(object):

    def __init__(self, base_category_class, args, kwds, base_category=None):
        self.cls = base_category_class
        self.args = args
        self.kwds = kwds
        self.base = base_category

    def __repr__(self):
        return 'Base {0} = {1}(*{2}, **{3})'.format(
            self.base, self.cls.__name__, self.args, self.kwds)


class CategoryHolder(object):

    def __init__(self, base_category, axioms, category_instance=None):
        if DEBUG:
            assert isinstance(axioms, tuple)
            assert len(set(axioms)) == len(axioms)
        self.base = base_category
        self.axioms = axioms
        self.category = category_instance

    def __repr__(self):
        return '{0} = {1} with ({2})'.format(self.category, self.base, self.axioms)
    


class CAMException(Exception):
    pass


class ClassNotFound(CAMException):
    pass


class InstanceNotFound(CAMException):
    pass


class CategoryAxiomModel(SageObject):
    
    def __init__(self):
        # The category classes that have been declared
        self._classes = []
        self._bases = []
        self._categories = []

    def _repr_(self):
        s = '=== Category Axiom Model ===\n'
        s += 'registered category classes:\n'
        for cls in self._classes:
            s += '    ' + repr(cls) + '\n'
        s += 'registered base category instances:\n'
        for base in self._bases:
            s += '    ' + repr(base) + '\n'
        s += 'registered instances of categories with axioms:\n'
        for cat in self._categories:
            s += '    ' + repr(cat) + '\n'
        s += '============================'            
        return s

    def find_class(self, base_category_class, axioms):
        for cls in self._classes:
            if cls.base is base_category_class and cls.axioms == axioms:
                return cls.cls
        return None
                
    def find_base(self, base_category_class, args, kwds):
        for base in self._bases:
            if base.cls is base_category_class and base.args == args and base.kwds == kwds:
                return base.base
        return None

    def find_category(self, base_category, axioms):
        for category in self._categories:
            if category.base is base_category and category.axioms == axioms:
                return category.category
        return None

    def register_category_class(self, category_class):
        if DEBUG:
            # ensure that a category is registered only once
            base_category_class = category_class.__base_category_class__
            axioms = tuple(sorted([ax() for ax in category_class.__axiom_classes__]))
            assert not self.find_class(base_category_class, axioms)
        #print("CAM category class {0} = {1} with {2}".format(
        #    category_class.__name__, category_class.__base_category_class__.__name__, 
        #    [axiom.__name__ for axiom in category_class.__axiom_classes__]))
        holder = ClassHolder(category_class)
        self._classes.append(holder)
        
    def register_base_category(self, base_category, base_category_class, args, kwds):
        if DEBUG:
            assert not self.find_base(base_category_class, args, kwds)
        # print("CAM base category instance of {0}".format(base_category_class))
        holder = BaseHolder(base_category_class, args, kwds, base_category)
        self._bases.append(holder)

    def register_category_instance(self, category):
        if DEBUG:
            assert not self.find_category(category._base_category, category.axioms())
        # print("CAM category instance of {0}".format(category.__class__))
        holder = CategoryHolder(category._base_category, category.axioms(), category)
        self._categories.append(holder)

    def add_relation(self, category_1, category_2, description=None):
        print("CAM relation {0} = {1}".format(category_1, category_2))
        # todo

    def applicable_axioms(self, category):
        """
        Return all axioms that can be added to the category.
        
        OUTPUT:

        A sorted duplicate-free tuple containing the axioms that can
        be added, but have not yet been added.
        """
        base_category_class = base_category.__base_category_class__
        axioms = set()
        for cls in self._classes:
            if cls.base is base_category_class:
                axioms.update(cls.axioms)
        axioms.difference(category.axioms())
        return tuple(sorted(axioms))

    def construct_base(self, base_category_class, args, kwds):
        """
        Construct a base category (without axioms)
        """
        base = self.find_base(base_category_class, args, kwds)
        if base is None:
            base = typecall(base_category_class, *args, **kwds)
            self.register_base_category(base, base_category_class, args, kwds)
        return base

    def construct_classcall(self, category_class, args, kwds, base_category=None):
        """
        Implementation of the classcall constructor of categories
        
        INPUT:

        - ``base_category_class`` -- the class of the base category

        - ``axiom_list`` -- a list/tuple/iterable of axioms
        
        - ``*args`` -- arguments to intstantiate the category

        - ``**kwds`` -- optional keyword arguments to instatiate the category
        
        OUTPUT:

        A category that is equal or equivalent to the specified
        category-with-axioms. Returns ``None`` if no such class has
        been registered.
        """
        if DEBUG:
            from sage.categories.category import Category
            assert issubclass(category_class, Category)
        base_category_class = category_class.__base_category_class__
        if base_category_class is category_class:
            # We are constructing the base class
            assert base_category is None  # you can't know yet
            return self.construct_base(base_category_class, args, kwds)
        if base_category is None:
            base_category = self.find_base(base_category_class, args, kwds)
        if base_category is None:
            raise ValueError('base category {0}({1}, {2}) not registered'.format(
                base_category_class.__name__, args, kwds))
        axioms = tuple(sorted(ax() for ax in category_class.__axiom_classes__))
        return self.construct_with_axioms(base_category, axioms)

    def construct_with_axioms(self, base_category, axioms):
        """
        Implementation of the Category.with_axioms factory method of categories
        
        INPUT:

        - ``category`` -- the base category

        - ``axiom_list`` -- a list/tuple/iterable of axioms
        
        OUTPUT:

        A category that is equal or equivalent to the specified
        category-with-axioms.
        """
        if DEBUG: 
            from sage.categories.category import Category
            assert isinstance(base_category, Category)
            assert isinstance(axioms, tuple)
            from sage.categories.axioms.axiom import Axiom
            assert all(isinstance(ax, Axiom) for ax in axioms)
            assert len(set(axioms)) == len(axioms)

        # TODO: bring base_category, axioms into normal form here

        category = self.find_category(base_category, axioms)
        if category is not None:
            return category
        base_category_class = base_category.__base_category_class__
        cls = self.find_class(base_category_class, axioms)
        if cls is not None:
            instance = typecall(cls, base_category=base_category)
            self.register_category_instance(instance)
            return instance
        print('not found', base_category_class, axioms)

    def construct_with_one_axiom(self, base_category, axiom):
        """
        Construct class for base category + single axiom
        
        If necessary, a dynamic class will be generated.

        EXAMPLES::

            sage: from sage.categories.category_with_axioms import Test2
            sage: t = Test2()
            sage: # category_axiom_model.construct_with_one_axiom(t, axioms.Associative())
            sage: t.with_axioms(axioms.Associative())
            <class 'sage.categories.metaclass.dynamic.Test2.Associative'>
        """
        base_category_class = base_category.__base_category_class__
        axioms = (axiom,)
        cls = self.find_class(base_category_class, axioms)
        if cls is None:
            cls = base_category_class.__metaclass__.__new_dynamic__(base_category, axioms)
        return cls

        
    

category_axiom_model = CategoryAxiomModel()
