"""
Categories with Axioms



"""

import importlib
import re
from sage.misc.cachefunc import cached_method, cached_function
from sage.misc.lazy_attribute import lazy_class_attribute
from sage.misc.lazy_import import LazyImport
from sage.misc.misc import call_method, uniq
from sage.categories.axioms.factory import axioms
from sage.categories.category import Category
from sage.categories.category_singleton import Category_singleton
from sage.categories.category_types import Category_over_base_ring
from sage.structure.dynamic_class import dynamic_class, DynamicMetaclass
from sage.misc.classcall_metaclass import typecall


#################################################################3
#  Category/Axiom Model

from sage.structure.sage_object import SageObject


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
        # The categories that have been declared as classes
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
            assert isinstance(base_category, Category)
            assert isinstance(axioms, tuple)
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

            sage: from sage.categories.category_with_axioms import Test2, category_axiom_model
            sage: t = Test2()
            sage: # category_axiom_model.construct_with_one_axiom(t, axioms.Associative())
            sage: t.with_axioms(axioms.Associative())
            <class 'sage.categories.category_with_axioms.dynamic.Test2.Associative'>
        """
        base_category_class = base_category.__base_category_class__
        axioms = (axiom,)
        cls = self.find_class(base_category_class, axioms)
        if cls is None:
            cls = base_category_class.__metaclass__.__new_dynamic__(base_category, axioms)
        return cls
        
    


category_axiom_model = CategoryAxiomModel()


#################################################################3
#  Metaclass 


from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.structure.dynamic_class import DynamicClasscallMetaclass
from sage.categories.axioms.axiom import Axiom



class CategoryWithAxiomsPlaceholder(object):
    """
    Auxiliary class to hold a list of axioms for the CategoryMetaclass
    """
    def __init__(self, axiom_classes):
        self.axiom_classes = axiom_classes

def subcategory_with_axioms(*axiom_classes):
    return CategoryWithAxiomsPlaceholder(axiom_classes)


class InnerCategory(object):
    """
    A marker for inner categories.
    
    This class is only useful as in inner class in a category. There,
    the category metaclass detects it while assembling the class,
    turns it into a sub-category defined by axioms.
    """
    def __init__(self, *args, **kwds):
        raise TypeError('only inner classes of categories may be of this type')


def separate_axiom_classes(class_list):
    """
    Helper to split a list of classes into categories and axioms

    Also, enforce that axioms are at the end of ``class_list``.

    EXAMPLES::

        sage: from sage.categories.category_with_axioms import separate_axiom_classes
        sage: class_list = [Category, int, axioms.Finite, axioms.Commutative]
        sage: cats, axioms = separate_axiom_classes(class_list)
        sage: cats
        [<class 'sage.categories.category.Category'>, <type 'int'>]
        sage: sorted(axioms)
        [<class 'sage.categories.axioms.commutative.Commutative'>, 
         <class 'sage.categories.axioms.finite.Finite'>]
    """
    remainder = []
    axioms = []
    for cls in class_list:
        if issubclass(cls, Axiom):
            axioms.append(cls)
        else:
            if len(axioms) > 0:
                raise TypeError('Axioms must be at the end of {0}'.format(class_list))
            remainder.append(cls)        
    return remainder, axioms


def replace_InnerCategory(classes_list, replacement):
    """
    Helper to replace :class:`InnerCategory` in a list
    
    EXAMPLES::

        sage: from sage.categories.category_with_axioms import replace_InnerCategory, InnerCategory
        sage: class_list = ['A', InnerCategory, 'B']
        sage: replace_InnerCategory(class_list, ['Replaced'])
        ('A', 'Replaced', 'B')
    """
    result = []
    for cls in classes_list:
        if cls is InnerCategory:
            result.extend(replacement)
        else:
            result.append(cls)
    return tuple(result)


def subcategory_with_axiom_factory(name, *axiom_classes):
    def self_with_axioms(self):
        """
        Return the category with added axioms.

        Calls ``self.with_axioms(...)``
        """
        axiom_list = [axiom() for axiom in axiom_classes]
        return self.with_axioms(*axiom_list)
    self_with_axioms.__name__ = name
    return self_with_axioms


def static_class_property():
    def static_class(self):
        """
        The category class without the dynamically generated ``Foo_with_category``
        """
        cls = self.__class__
        if isinstance(cls, DynamicMetaclass):
            return cls.__base__
        else:
            return cls
    static_class.__name__ = name
    return property(static_class)


def add_generated_attributes(category, base_category, axiom_classes):
    """
    Add the autogenerated attributes to the category class
    """
    category.__class_without_category__ = category
    category.__base_category_class__ = base_category
    category.__axiom_classes__ = axiom_classes


def read_off_axioms(bases, dictionary):
    """
    Figure out the axioms for a class
    
    INPUT:

    - ``bases`` -- tuple of base classes.

    - ``dictionary`` -- the class dictionary.

    OUTPUT:

    A list of axioms. A ``TypeError`` is raised if both inheritance
    and the ``__axiom_classes__`` attribute is used to specify the
    axioms.

    EXAMPLES::

        sage: from sage.categories.category_with_axioms import read_off_axioms
        sage: read_off_axioms((Category,), 
        ....:                 dict(__axiom_classes__=frozenset([axioms.Finite])))
        frozenset([<class 'sage.categories.axioms.finite.Finite'>])

        sage: read_off_axioms((Category, axioms.Finite), dict())
        frozenset([<class 'sage.categories.axioms.finite.Finite'>])

        sage: read_off_axioms((Category, axioms.Finite), 
        ....:                 dict(__axiom_classes__=frozenset([axioms.Finite])))
        Traceback (most recent call last):
        ...
        TypeError: must not both inherit from Axiom and set the __axiom_classes__ attribute
    """
    base_remainder, base_axioms = separate_axiom_classes(bases)
    attr_axioms = list(dictionary.get('__axiom_classes__', []))
    if len(attr_axioms) > 0 and len(base_axioms) > 0:
        raise TypeError('must not both inherit from Axiom and '
                        'set the __axiom_classes__ attribute')
    return frozenset(base_axioms + attr_axioms)


class CategoryMetaclass(ClasscallMetaclass): 
    """
    Metaclass to define subcategories by inner classes
    
    This is the solution to the following chicken-and-egg problem:
    While the inner class is constructed, the containing class does
    not yet exist.
    """
    
    @classmethod
    def __new_dynamic__(cls, base_category, axioms):
        """
        Construct a dynamic category-with-axiom class

        EXAMPLES::

            sage: from sage.categories.category_with_axioms import Test1, CategoryMetaclass, CategoryWithAxioms
            sage: cls = CategoryWithAxioms.__metaclass__.__new_dynamic__(Test1(), [axioms.WithBasis()])
            sage: cls().explain()
            Category sage.categories.category_with_axioms.dynamic.Test1.WithBasis_with_category
            Set by meta class:
                __class_without_category__ = dynamic.Test1.WithBasis
                __base_category_class__ = Test1
                __axiom_classes__ = frozenset([WithBasis])
            Set by the category axiom manager:
                _base_category = Category of test1
            MRO = (dynamic.Test1.WithBasis_with_category, dynamic.Test1.WithBasis, CategoryWithAxioms, NewCategory, Category, UniqueRepresentation, CachedRepresentation, WithEqualityById, SageObject, dynamic.Test1.WithBasis.subcategory_class, object)
        """
        # print('Meta {0} {1} {2}'.format(cls, base_category, axioms))
        base_category_class = base_category.__base_category_class__
        axiom_classes = frozenset(ax.__class__ for ax in axioms)
        name = 'dynamic.{0}.{1}'.format(
            base_category_class.__name__, 
            ''.join(map(str, axioms)))
        bases = (CategoryWithAxioms,)
        dictionary = {
            '__base_category_class__': base_category_class,
            '__axiom_classes__': axiom_classes,
        }
        category = super(CategoryMetaclass, cls).__new__(
            cls, name, bases, dictionary)
        add_generated_attributes(category, base_category_class, axiom_classes)
        category_axiom_model.register_category_class(category)
        return category

    def __new__(cls, name, bases, dictionary):
        """
        Construct a new Category subclass

        .. NOTE::

            The metacategory works on classes only, it never sees
            instances of the categories. So in the implementation of
            the metaclass, note that ``category`` really always refers
            to the class (type) of the category. However, in the
            interest of breverity, we do not add ``_class`` to every
            variable name.
        """
        # don't do it twice for the top-most dynamic class
        # TODO: we need some proper hook there
        #       shoud generate the dynamic _with_category in the metaclass 
        #       as well and only set its attributes in Category.__init__
        if issubclass(cls, DynamicMetaclass):
            return super(CategoryMetaclass, cls).__new__(
                cls, name, bases, dictionary)

        # print('Meta {0} {1} {2} {3}'.format(cls, name, bases, dictionary))
        base_category = dictionary.pop('__base_category_class__', None)
        axioms = read_off_axioms(bases, dictionary)
        
        # Scan through the dictionary of the class we are going to construct
        inner_categories = []
        for key, attr in dictionary.items():
            if isinstance(attr, CategoryWithAxiomsPlaceholder):
                # Turn "SomeName = subcategory_with_axioms(axioms.Foo)" 
                # into a factory method for category-with-Foo 
                dictionary[key] = subcategory_with_axiom_factory(key, *attr.axiom_classes)
            elif isinstance(attr, type) and issubclass(attr, InnerCategory):
                # Turn inner class into a factory method for new category-with-axioms
                attr_remainder, attr_axioms = separate_axiom_classes(attr.__bases__)
                inner_categories.append([key, attr])
                if attr.__bases__[0] is not InnerCategory:
                    raise TypeError('InnerCategory must be the first parent class for '
                                    'for subcategory {0}'.format(key))
                dictionary[key] = subcategory_with_axiom_factory(key, *attr_axioms)
            elif isinstance(attr, type) and issubclass(attr, Category):
                raise TypeError('subcategory {0} defined by inner class must derive '
                                'from InnerCategory and axioms only'.format(key))
            else:
                pass

        # Generate the category class
        category = super(CategoryMetaclass, cls).__new__(
            cls, name, bases, dictionary)
        if base_category is None:
            base_category = category
        add_generated_attributes(category, base_category, axioms)
        category_axiom_model.register_category_class(category)

        # Generate the subcategory classes
        for key, subcategory in inner_categories:
            # Base classes for the category defined by the inner class
            sub_bases = replace_InnerCategory(subcategory.__bases__, [CategoryWithAxioms])

            sub_dictionary = dict(subcategory.__dict__)
            if '__base_category_class__' in sub_dictionary:
                raise ValueError('InnerCategory must not set __base_category_class__, '
                                 'it wil automatically be set to the outermost class')
            sub_dictionary['__base_category_class__'] = base_category

            sub_axioms = read_off_axioms(sub_bases, sub_dictionary)
            if axioms.intersection(sub_axioms) != set():
                raise ValueError('subcategory from inner class {0} adds axiom that is '
                                     'already in the containing class'.format(key))
            sub_dictionary['__axiom_classes__'] = axioms.union(sub_axioms)

            # strip out the axioms from the oop base classes
            sub_bases, _ = separate_axiom_classes(sub_bases)
            sub_bases = tuple(sub_bases)

            # Generate the subcategory class with a recursive call to __new__
            subcategory = cls.__new__(cls, name + '.' + key, sub_bases, sub_dictionary)

        return category

    def __init__(cls, name, bases, dictionary):
        #print '-----------------------------------'
        #print "Initializing class", name, cls
        #print cls
        #print bases
        #print dictionary
        super(CategoryMetaclass, cls).__init__(name, bases, dictionary)
        #print '__init__', name, bases, dictionary
        #cls._base_category


class DynamicCategoryMetaclass(DynamicMetaclass, CategoryMetaclass):
    pass



#################################################################3
#  New category class

DEBUG = True

class NewCategory(Category):
    r"""
    This class is only for demonstration purposes, it should
    be merged with Category.
    """
    __metaclass__ = CategoryMetaclass

    # The class of the base category, auto-generated by the CategoryMetaclass
    # __base_category_class__ = NewCategory

    # Frozenset of axiom classes, auto-generated by the CategoryMetaclass
    # __axiom_classes__ = frozenset()
    
    @staticmethod
    def __classcall__(cls, *args, **kwds):
        """
        Intercept all classclass and redirect to the category axiom manager

        EXAMPLES::

            sage: FiniteGroups()
            Category of finite groups
            sage: ModulesWithBasis(ZZ)
            Category of modules with basis over Integer Ring
            sage: AlgebrasWithBasis(QQ)
            Category of algebras with basis over Rational Field

        This is relevant when e.g. ``Foos(**)`` does some non trivial
        transformations::

            sage: Modules(QQ) is VectorSpaces(QQ)
            True
            sage: type(Modules(QQ))
            <class 'sage.categories.vector_spaces.VectorSpaces_with_category'>

            sage: ModulesWithBasis(QQ) is VectorSpaces(QQ).WithBasis()
            True
            sage: type(ModulesWithBasis(QQ))
            <class 'sage.categories.vector_spaces.VectorSpaces.WithBasis_with_category'>
        """
        #print('{2}__classcall__({0}, {1})__'.format(args, kwds, cls))
        return category_axiom_model.construct_classcall(cls, args, kwds)
        # return super(NewCategory, cls).__classcall__(cls, *args, **kwds)


    def __init__(self, base_category=None):
        """
        Initializes this category.

        EXAMPLES::

            sage: class SemiprimitiveRings(Category):
            ....:     def super_categories(self):
            ....:         return [Rings()]
            ....:
            ....:     class ParentMethods:
            ....:         def jacobson_radical(self):
            ....:             return self.ideal(0)
            ....:
            sage: C = SemiprimitiveRings()
            sage: C
            Category of semiprimitive rings
            sage: C.__class__
            <class '__main__.SemiprimitiveRings_with_category'>

        TESTS::

            sage: C = Sets.Finite(); C
            Category of finite sets
            sage: type(C)
            <class 'sage.categories.finite_sets.FiniteSets_with_category'>
            sage: type(C).__base__.__base__
            <class 'sage.categories.category_with_axiom.CategoryWithAxiom_singleton'>

            sage: TestSuite(C).run()
        """
        cls = self.__class__
        self.__class__ = dynamic_class('{0}_with_category'.format(cls.__name__),
                                       (cls, self.subcategory_class),
                                       cache=False, reduction=None,
                                       doccls=cls, 
                                       metaclass=DynamicCategoryMetaclass)
        self._base_category = self if base_category is None else base_category
        if DEBUG:
            self._init_debug()

    def _init_debug(self):
        """
        Perform additional checks on instances
        """
        assert issubclass(self.__base_category_class__, Category)
        assert isinstance(self._base_category, self.__base_category_class__)
        assert isinstance(self.__axiom_classes__, frozenset)
        assert self.__class__.__base__ is self.__class_without_category__

    def _test_category_with_axioms(self, **options):
        r"""
        Run generic tests on this category with axioms.

        .. SEEALSO:: :class:`TestSuite`.

        This check that an axiom category of a
        :class:`Category_singleton` is a singleton category, and
        similarwise for :class`Category_over_base_ring`.

        EXAMPLES::

            sage: Sets().Finite()._test_category_with_axiom()
            sage: Modules(ZZ).FiniteDimensional()._test_category_with_axiom()
        """
        tester = self._tester(**options)
        base = self.base_category()
        if isinstance(base, Category_singleton):
            tester.assertIsInstance(self, NewCategory_singleton)
        if isinstance(base, Category_over_base_ring):
            tester.assertIsInstance(self, NewCategory_over_base_ring)

    @staticmethod
    def _repr_object_names_static(category, axiom_list):
        r"""
        INPUT:

        - ``base_category`` -- a category
        - ``axioms`` -- a list or iterable of strings

        EXAMPLES::

            sage: from sage.categories.category_with_axiom import CategoryWithAxiom
            sage: CategoryWithAxiom._repr_object_names_static(Semigroups(), [axioms.Flying(), axioms.Blue()])
            'flying blue semigroups'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ), 
            ....:     [axioms.Flying(), axioms.WithBasis(), axioms.Blue()])
            'flying blue algebras with basis over Rational Field'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ), [axioms.WithBasis()])
            'algebras with basis over Rational Field'
            sage: CategoryWithAxiom._repr_object_names_static(
            ....:     Sets().Finite().Subquotients(), [axioms.Finite()])
            'subquotients of finite sets'
            sage: CategoryWithAxiom._repr_object_names_static(Monoids(), [axioms.Unital()])
            'monoids'
            sage: CategoryWithAxiom._repr_object_names_static(Algebras(QQ['x']['y']), 
            ....:     [axioms.Flying(), axioms.WithBasis(), axioms.Blue()])
            'flying blue algebras with basis over Univariate Polynomial Ring in y over Univariate Polynomial Ring in x over Rational Field'

        If the axioms is a set or frozen set, then they are first
        sorted using :func:`canonicalize_axioms`::

            sage: CategoryWithAxiom._repr_object_names_static(Semigroups(), 
            ....:     set([axioms.Finite(), axioms.Commutative(), axioms.Facade()]))
            'facade finite commutative semigroups'

        .. SEEALSO:: :meth:`_repr_object_names`

        .. NOTE:: The logic here is shared between :meth:`_repr_object_names`
            and :meth:`.category.JoinCategory._repr_object_names`
        """
        base_category = category._without_axioms(named=True)
        if isinstance(base_category, NewCategory): # Smelly runtime type checking
            result = super(NewCategory, base_category)._repr_object_names()
        else:
            result = base_category._repr_object_names()
        for axiom in reversed(sorted(set(axiom_list))):
            # TODO: find a more generic way to handle the special cases below
            if axiom in base_category.axioms():
                # If the base category already has this axiom, we
                # need not repeat it here. See the example with
                # Sets().Finite().Subquotients() or Monoids()
                continue
            if axiom == axioms.WithBasis():
                result = result.replace(" over ", " with basis over ", 1)
            elif axiom == axioms.Connected() and "graded " in result:
                result = result.replace("graded ", "graded connected ", 1)
            else:
                result = ' '.join([str(ax).lower() for ax in self.axioms()]) + " " + result
        return result

    def _repr_object_names(self):
        r"""
        The names of the objects of this category, as used by `_repr_`

        .. SEEALSO:: :meth:`Category._repr_object_names`

        EXAMPLES::

            sage: FiniteSets()._repr_object_names()
            'finite sets'
            sage: AlgebrasWithBasis(QQ).FiniteDimensional()._repr_object_names()
            'finite dimensional algebras with basis over Rational Field'
            sage: Monoids()._repr_object_names()
            'monoids'
            sage: Semigroups().Unital().Finite()._repr_object_names()
            'finite monoids'
            sage: Algebras(QQ).Commutative()._repr_object_names()
            'commutative algebras over Rational Field'

        .. NOTE::

            This is implemented by taking _repr_object_names from
            self._without_axioms(named=True), and adding the names
            of the relevant axioms in appropriate order.
        """
        return NewCategory._repr_object_names_static(self, self.axioms())

    def base_category(self):
        r"""
        Return the base category of ``self``.

        EXAMPLES::

            sage: C = Sets.Finite(); C
            Category of finite sets
            sage: C.base_category()
            Category of sets
            sage: C._without_axioms()
            Category of sets

        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects, CategoryWithAxiom
            sage: C = TestObjects().Commutative().Facade()
            sage: assert isinstance(C, CategoryWithAxiom)
            sage: C._without_axioms()
            Category of test objects
        """
        return self._base_category

    def __reduce__(self):
        r"""
        Implement the pickle protocol.

        This overides the implementation in
        :meth:`UniqueRepresentation.__reduce__` in order to not
        exposes the implementation detail that, for example, the
        category of magmas which distribute over an associative
        additive magma is implemented as
        ``DistributiveMagmasAndAdditiveMagmas.AdditiveAssociative.AdditiveCommutative``
        and not
        ``DistributiveMagmasAndAdditiveMagmas.AdditiveCommutative.AdditiveAssociative``::
        """
        return (call_method, (self.base_category(), "_with_axioms", self.defining_axioms()))

    def with_axioms(self, *axioms_list):
        base_category = self._base_category
        axioms = tuple(uniq(self.axioms() + tuple(axioms_list)))
        if len(axioms) == 0:
            return base_category
        if len(axioms) == 1:
            return category_axiom_model.construct_with_one_axiom(base_category, axioms[0])
        category = category_axiom_model.construct_with_axioms(base_category, axioms)
        if isinstance(category, tuple):
            return Category.join(category)
        else:
            return category
        

    @cached_method
    def axioms(self):
        r"""
        Return the axioms known to be satisfied by all the objects of ``self``.

        .. SEEALSO:: :meth:`Category.axioms`
        """
        return tuple(sorted(axiom() for axiom in self.__axiom_classes__))

    def explain(self):
        """
        Print information about the category

        #EXAMPLES::

        #    sage: from sage.categories.category_with_axioms import Test1
        #    sage: Test1().explain()
        #    sage: Test1().with_axioms(axioms.Finite(), axioms.Commutative()).explain()
        """
        cls = self.__class__
        print('Category {0}.{1}'.format(cls.__module__, cls.__name__))
        print('Set by meta class:')
        print('    __class_without_category__ = {0}'.format(
            self.__class_without_category__.__name__))
        print('    __base_category_class__ = {0}'.format(
            self.__base_category_class__.__name__))
        print('    __axiom_classes__ = frozenset([{0}])'.format(
            ', '.join(map(str, self.axioms()))))
        print('Set by the category axiom manager:')
        print('    _base_category = {0}'.format(self._base_category))
        print('MRO = ({0})'.format(
            ', '.join([parent.__name__ for parent in cls.mro()])))




class CategoryWithAxioms(NewCategory):

    @cached_method
    def super_categories(self):
        """
        Return a list of the (immediate) super categories of
        ``self``, as per :meth:`Category.super_categories`.

        This implements the property that if ``As`` is a subcategory
        of ``Bs``, then the intersection of As with ``FiniteSets()``
        is a subcategory of ``As`` and of the intersection of ``Bs``
        with ``FiniteSets()``.

        EXAMPLES::

            sage: FiniteSets().super_categories()
            [Category of sets]

            sage: FiniteSemigroups().super_categories()
            [Category of semigroups, Category of finite enumerated sets]

        EXAMPLES:

        A finite magma is both a magma and a finite set::

            sage: Magmas().Finite().super_categories()
            [Category of magmas, Category of finite sets]

        TESTS::

            sage: from sage.categories.category_with_axiom import TestObjects
            sage: C = TestObjects().FiniteDimensional().Unital().Commutative().Finite()
            sage: sorted(C.super_categories(), key=str)
            [Category of finite commutative test objects,
             Category of finite dimensional commutative unital test objects,
             Category of finite finite dimensional test objects]
        """
        return []


##############################################################################
# Utilities and tests tools

# manual construction

class Test1(NewCategory):
    """
    EXAMPLES::

        sage: from sage.categories.category_with_axioms import Test1
        sage: Test1().__class__.__base__
        <class 'sage.categories.category_with_axioms.Test1'>
        sage: Test1().FiniteCommutative1().__class__.__base__
        <class 'sage.categories.category_with_axioms.FiniteCommutativeTest1'>
    """

    def super_categories(self):
        return []

    Infinite1 = subcategory_with_axioms(axioms.Infinite)  # parsed by the metaclass
        
    FiniteCommutative1 = subcategory_with_axioms(axioms.Finite, axioms.Commutative)

    

class FiniteCommutativeTest1(CategoryWithAxioms):
    
    __base_category_class__ = Test1

    __axiom_classes__ = frozenset([axioms.Commutative, axioms.Finite])

    def is_finite_commutative_1(self):
        """ Example of extra code """
        pass


# Inner class syntactic sugar

class Test2(NewCategory):
    """
    EXAMPLES::
    
        sage: from sage.categories.category_with_axioms import Test2
        sage: Test2().FiniteCommutative2()
        Category of finite commutative2

        sage: Test2().FiniteCommutative2().explain()
        Category sage.categories.category_with_axioms.Test2.FiniteCommutative2_with_category
        Set by meta class:
            __class_without_category__ = Test2.FiniteCommutative2
            __base_category_class__ = Test2
            __axiom_classes__ = frozenset([Finite, Commutative])
        Set by the category axiom manager:
            _base_category = Category of test2
        MRO = (Test2.FiniteCommutative2_with_category, Test2.FiniteCommutative2, 
               CategoryWithAxioms, NewCategory, Category, UniqueRepresentation, 
               CachedRepresentation, WithEqualityById, SageObject,  
               Test2.FiniteCommutative2.subcategory_class, object)
    """

    def super_categories(self):
        return []

    Infinite2 = axioms.Infinite
        
    class Finite2(InnerCategory, axioms.Finite):
        def is_finite_2(self):
            """ Example of extra code """
            return True

        class Associative2(InnerCategory, axioms.Associative):
            def is_finite_associative_2(self):
                """ Example of extra code """
                return True

    class FiniteCommutative2(InnerCategory, axioms.Finite, axioms.Commutative):
        def is_finite_commutative_2(self):
            """ Example of extra code """
            return True


    
# relation Test2.Finite is FiniteTest2

#category_axiom_model.register(FiniteTest2)



# class Blahs(Category):

#     def super_categories(self):
#         """
#         TESTS::

#              sage: from sage.categories.category_with_axioms import Blahs
#              sage: Blahs().super_categories()
#              [Category of sets]
#             sage: TestSuite(Blahs()).run()
#         """
#         from sage.categories.sets_cat import Sets
#         return [Sets()]

#     class SubcategoryMethods:
#         FiniteDimensional = with_axioms(axioms.FiniteDimensional())
#         Commutative       = with_axioms(axioms.Commutative())
#         Unital            = with_axioms(axioms.Unital())
#         Connected         = with_axioms(axioms.Connected())
#         Flying            = with_axioms(axioms.Flying())
#         Blue              = with_axioms(axioms.Blue())

#     class FiniteDimensional(NewCategory):
#         pass
#     class Commutative(NewCategory):
#         pass
#     class Connected(NewCategory):
#         pass
#     class Unital(NewCategory):
#         class Blue(NewCategory):
#             pass
#     class Flying(NewCategory):
#         def extra_super_categories(self):
#             """
#             This illustrates a way to have an axiom imply another one.

#             Here, we want ``Flying`` to imply ``Unital``, and to put
#             the class for the category of unital flying blahs in
#             ``Blahs.Flying`` rather than ``Blahs.Unital.Flying``.

#             TESTS::

#                 sage: from sage.categories.category_with_axiom import Blahs, TestObjects, Bars
#                 sage: Blahs().Flying().extra_super_categories()
#                 [Category of unital blahs]
#                 sage: Blahs().Flying()
#                 Category of flying unital blahs
#             """
#             return [Blahs().Unital()]
