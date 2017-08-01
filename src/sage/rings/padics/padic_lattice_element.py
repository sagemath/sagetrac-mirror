import _weakref as weakref

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity

from sage.rings.padics.generic_nodes import pAdicRingBaseGeneric
from sage.rings.padics.padic_generic_element import pAdicGenericElement
from sage.rings.padics.lattice_precision import pRational


class pAdicLatticeElement(pAdicGenericElement):
    def __init__(self, parent, x, prec=None, dx=[], dx_mode='linear_combinaison', valuation=None, check=True):
        self._parent = parent
        pAdicGenericElement.__init__(self, parent)
        self._precision = parent.precision()
        if check:
            if isinstance(x, pAdicGenericElement):
                if parent.prime() != x.parent().prime():
                    raise TypeError("conversion between different p-adic rings/fields not supported")
                if prec is None:
                    prec = x.precision_absolute()
                else:
                    prec = min(prec, x.precision_absolute())
            if isinstance(parent, pAdicRingBaseGeneric):
                x = ZZ(x)
            else:
                x = QQ(x)
            cap = parent.precision_cap()
            if prec is None or prec > cap:
                prec = cap
        self._precision.new_element(self, dx, bigoh=prec, dx_mode=dx_mode)
        if isinstance(x, pRational):
            self._value = x.reduce(prec)
        else:
            self._value = pRational(parent.prime(), QQ(x)).reduce(prec)
        self._approx_one = pRational(parent.prime(), 1)
        self._approx_minusone = pRational(parent.prime(), -1)

    def __hash__(self):
        return id(self)

    def _is_base_elt(self, p):
        return p == self._parent.prime()

    def approximation(self):
        prec = self.precision_absolute()
        app = self._value.reduce(prec)
        return app.value()

    def value(self):
        return self._value.value()

    def precision_lattice(self):
        return self._precision

    def precision_absolute(self):
        return self._precision.precision_absolute(self)

    def valuation(self, secure=False):
        p = self._parent.prime()
        val = self._value.valuation()
        prec = self.precision_absolute()
        if val < prec: 
            return val
        elif secure:
            raise PrecisionError("Not enough precision")
        else:
            return prec

    def precision_relative(self, secure=False):
        return self.precision_absolute() - self.valuation(secure=secure)

    def _cmp_(self, other):
        if (self-other).is_zero():
            return 0
        else:
            return cmp(self.lift(), other.lift())

    def _add_(self, other):
        x = self._value + other._value
        dx = [  [self, self._approx_one], 
               [other, self._approx_one] ]
        return self.__class__(self._parent, x, self._parent.precision_cap(), dx=dx, check=False)

    def _sub_(self, other):
        x = self._value - other._value
        dx = [  [self, self._approx_one], 
               [other, self._approx_minusone] ]
        return self.__class__(self._parent, x, self._parent.precision_cap(), dx=dx, check=False)

    def _mul_(self, other):
        x_self = self._value
        x_other = other._value
        x = x_self * x_other
        dx = [  [self, x_other],
               [other, x_self ] ]
        return self.__class__(self._parent, x, self._parent.precision_cap(), dx=dx, check=False)

    def _div_(self, other):
        if other.is_zero():
            raise PrecisionError("cannot divide by something indistinguishable from zero")
        p = self._parent.prime()
        cap = self._parent.precision_cap()
        x_self = self._value
        x_other = other._value
        x = x_self / x_other
        # dx = (1/other)*dself - (self/other^2)*dother
        dx = [  [self, self._approx_one/x_other],
               [other, -x_self/(x_other*x_other)] ]
        return self.__class__(self._parent.fraction_field(), x, self._parent.precision_cap(), dx=dx, check=False)

    def __invert__(self):
        if self.is_zero():
            raise PrecisionError("cannot invert something indistinguishable from zero")
        p = self._parent.prime()
        cap = self._parent.precision_cap()
        x_self = self._value
        x = self._approx_one / x_self
        # dx = -(1/self^2)*dself
        dx = [  [self, self._approx_minusone/(x_self*x_self)] ]
        return self.__class__(self._parent.fraction_field(), x, self._parent.precision_cap(), dx=dx, check=False)

    def add_bigoh(self, prec):
        x = self._value
        dx = [ [self, self._approx_one ] ]
        return self.__class__(self._parent, x, prec, dx=dx, check=False)

    def lift_to_precision(self, prec=None, infer_precision=False):
        from warnings import warn
        warn("use lift_to_precision with extreme caution in the framework of lattice precision")
        parent = self._parent
        if prec is None:
            prec = parent.precision_cap()
        if infer_precision:
            lift = self.__copy__()
            parent.precision().lift_to_precision(lift, prec)
        else:
            if prec < self.precision_absolute():
                prec = self.precision_absolute()
            lift = self.__class__(parent, self._value, prec, check=False)
        return lift

    def _is_exact_zero(self):
        return False

    def _is_inexact_zero(self):
        absprec = self.precision_absolute()
        return self._value.valuation() >= absprec

    def is_zero(self, prec=None):
        absprec = self.precision_absolute()
        if prec is None:
            prec = absprec
        else:
            prec = min(absprec, prec)
        return self._value.valuation() >= prec

    def value(self):
        return self._value

    def lift(self):
        return self._value.value()

    def rational_reconstruction(self):
        return QQ(self._value.value())

    def __rshift__(self, n):
        return self << (-n)

    def __lshift__(self, n):
        p = self._parent.prime()
        powp = pRational(p, ZZ(1), n)
        x = self._value * powp
        dx = [ [self, powp] ]
        return self.__class__(self._parent, x, self._parent.precision_cap(), dx=dx, check=False)

    def unit_part(self):
        v = self.valuation(secure=True)
        return self >> v

    def val_unit(self):
        v = self.valuation(secure=True)
        return v, self >> v

    def __copy__(self, parent=None):
        if parent is None:
            parent = self._parent
        else:
            if parent.precision() is not self._parent.precision():
                raise TypeError("parent must share the same precision object")
            if isinstance(parent, pAdicRingBaseGeneric) and self.valuation() < 0:
                raise ValueError("element of negative valuation cannot be convert to the integer ring")
        dx = [ [ self, self._approx_one ] ]
        return self.__class__(parent, self._value, self._parent.precision_cap(), dx=dx, check=False)

    def list(self, lift_mode='simple', start_val=None):
        # TODO: implement other lift modes
        p = self._parent.prime()
        prec = self.precision_absolute()
        return self._value.list(prec)
