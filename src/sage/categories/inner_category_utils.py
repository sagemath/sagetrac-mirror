"""
Utilities to Define Subcategories from Inner (Nested) Classes

These are defined in a separate file to avoid circular imports.
"""


class CategoryWithAxiomsPlaceholder(object):
    """
    Auxiliary class to hold a list of axioms for the CategoryMetaclass
    """
    def __init__(self, axiom_classes):
        self.axiom_classes = axiom_classes

    def __repr__(self):
        return 'object should have been replaced by the category metaclass'


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

