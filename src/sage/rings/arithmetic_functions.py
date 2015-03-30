from sage.categories.rings import Rings
from sage.structure.parent import Parent
from sage.rings.all import ZZ
from sage.structure.element import get_coercion_model


def ArithmeticFunction(f, repr_function=None, repr_latex=None, check=True):
    return ArithmeticFunctions().element_class(ArithmeticFunctions(), f, repr_function, repr_latex, check, convert=True)


class ArithmeticFunctions(Parent):
    r"""
    Arithmetic function ring with multiplication given by convolution.
    """

    from arithmetic_function import ArithmeticFunctionElement
    Element = ArithmeticFunctionElement

    def __init__(self, codomain=ZZ):
        r"""
        Arithmetic function ring with multiplication given by convolution.
        The functions are from N to the given codomain.
        """

        if not (codomain in Rings()):
            raise ValueError("{} is not in {}".format(codomain, Rings()))
        self.codomain = codomain
        self.reset_variables()

        Parent.__init__(self, category=Rings())

    def _repr_(self):
        r"""
        Return the string representation of ``self``.
        """

        return "Ring of functions from N to {}".format(self.codomain)

    def _element_constructor_(self, el):
        r"""
        Return ``el`` coerced/converted into this ring.
        """

        if isinstance(el, self.element_class):
            return self.element_class(self, el.function, el.repr_function, el.repr_latex, convert=False)
        elif el in self.codomain:
            return self.constant_function(el)
        elif hasattr(el, "__call__"):
            # If a function is given as an argument we construct a corresponding element (without any representations)
            return self.element_class(self, el)
        else:
            raise ValueError("Couldn't convert the given element el = {}".format(el))

    def _coerce_map_from_(self, S):
        r"""
        Return whether or not there exists a coercion from ``S`` to ``self``.
        """

        if self.codomain.has_coerce_map_from(S):
            return True
        if isinstance(S, self.__class__) and self.codomain.has_coerce_map_from(S.codomain):
            return True

    def _an_element_(self):
        r"""
        Return an element of ``self``.
        """

        return self(self.one())

    def extend_codomain(self, codomain):
        cm = get_coercion_model()
        new_codomain = cm.common_parent(self.codomain, codomain)

        return self.__class__(codomain=new_codomain)

    def reset_variables(self, variable_list=None):
        if variable_list:
            self._variable_list = variable_list
        else:
            self._variable_list = ["n","m","k","l","a","b","c","d","x","y","z","w"]

    def get_variable(self):
        if len(self._variable_list) > 0:
            return self._variable_list.pop(0)
        else:
            raise Exception("No more variable names available!")

    def constant_function(self, el):
        el = (self.codomain)(el)

        def f(n):
            return el
        def repr_f(var, get_var):
            return "{}".format(el)
        def repr_latex(var, get_var):
            return "{}".format(latex(el))

        return self.element_class(self, f, repr_f, repr_latex, convert=False)

    def trivial_function(self, el):
        el = (self.codomain)(el)
        zero = (self.codomain)(0)

        def f(n):
            return el if n==ZZ(1) else zero
        def repr_f(var, get_var):
            return "({} if {}==1 else 0)".format(el, var)
        def repr_latex(var, get_var):
            return "el_{{{}=1}}".format(var)

        return self.element_class(self, f, repr_f, repr_latex, convert=False)

    def one(self):
        return self.trivial_function(self.codomain.one())

    def example(self):
        from sage.rings.arith import divisors

        def f(n):
            return len(divisors(n))
        def repr_f(var, get_var):
            return "sigma({}, 0)".format(var)
            #return "len(divisors({}))".format(var)
        def repr_latex(var, get_var):
            return "\\sigma_0\\left({}\\right)".format(var)
            #var2=get_var()
            #return "\\left(\\sum_{{\\left.{}\\middle|{}\\right.}}1\\right)".format(var2,var)

        return self.element_class(self, f, repr_f, repr_latex, convert=False)
