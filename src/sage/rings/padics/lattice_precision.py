import _weakref as weakref
from sage.misc.misc import walltime

from sage.structure.sage_object import SageObject
from sage.structure.unique_representation import UniqueRepresentation

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.rings.infinity import Infinity


# Class pRational
#################

# It is a temporary class implementing rational numbers
# as approximations of p-adic numbers

class pRational:
    def __init__(self, p, x, exponent=0, valuation=None):
        self.p = p
        if x in ZZ:
            self.x = ZZ(x)
        else:
            self.x = x
        self.exponent = exponent
        self._valuation = valuation

    def __repr__(self):
        if self.exponent == 0:
            return str(self.x)
        else:
            return "%s^%s * %s" % (self.p, self.exponent, self.x)

    def reduce(self, prec):
        if prec is Infinity:
            return self
        x = self.x
        exp = self.exponent
        if x.parent() is ZZ:
            if prec > exp:
                x = x % (self.p ** (prec-exp))
            else:
                x = 0
        elif x.parent() is QQ:
            num = x.numerator()
            denom = x.denominator()
            valdenom = denom.valuation(self.p)
            denom /= self.p ** valdenom
            exp -= valdenom
            modulo = self.p ** (prec - exp)
            # probably we should use Newton iteration instead 
            # (but it is actually slower for now - Python implementation)
            _, inv, _ = denom.xgcd(modulo)
            x = (num*inv) % modulo
        if self.x == 0:
            val = Infinity
        else:
            val = self._valuation
        return self.__class__(self.p, x, exp, valuation=val)

    def normalize(self):
        if self.x == 0:
            self.exponent = 0
        else:
            val = self.valuation()
            exp = self.exponent
            self.x /= self.p ** (val-exp)
            if self.x in ZZ:
                self.x = ZZ(self.x)
            self.exponent = val

    def valuation(self):
        if self._valuation is None:
            valx = self.x.valuation(self.p)
            self._valuation = self.exponent + valx
        return self._valuation

    def is_p_power(self):
        self.normalize()
        return self.x == 1

    def is_zero(self):
        return self.x == 0

    def __add__(self, other):
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if self._valuation is None or other._valuation is None:
            val = None
        elif self._valuation < other._valuation:
            val = self._valuation
        elif self._valuation > other._valuation:
            val = other._valuation
        else:
            val = None
        if sexp < oexp:
            return self.__class__(p, self.x + other.x * p**(oexp-sexp), sexp, valuation=val)
        else:
            return self.__class__(p, self.x * p**(sexp-oexp) + other.x, oexp, valuation=val)

    def __sub__(self, other):
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if self._valuation is None or other._valuation is None:
            val = None
        elif self._valuation < other._valuation:
            val = self._valuation
        elif self._valuation > other._valuation:
            val = other._valuation
        else:
            val = None
        if sexp < oexp:
            return self.__class__(p, self.x - other.x * p**(oexp-sexp), sexp, valuation=val)
        else:
            return self.__class__(p, self.x * p**(sexp-oexp) - other.x, oexp, valuation=val)

    def __neg__(self):
        return self.__class__(self.p, -self.x, self.exponent, valuation=self._valuation)

    def __mul__(self, other):
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation + other._valuation
        return self.__class__(self.p, self.x * other.x, self.exponent + other.exponent, valuation=val)

    def __div__(self, other):
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = self._valuation - other._valuation
        return self.__class__(self.p, self.x / other.x, self.exponent - other.exponent, valuation=val)

    def __lshift__(self, n):
        if self._valuation is None:
            val = None
        else:
            val = self._valuation + n
        return self.__class__(self.p, self.x, self.exponent + n, valuation=val)

    def __rshift__(self, n):
        return self << (-n)

    def unit_part(self):
        if self.is_zero():
            raise ValueError("the unit part of zero is not defined")
        p = self.p
        val = self.valuation()
        x = self.x / (p ** (val-self.exponent))
        return self.__class__(p, x, 0, valuation=0)

    def xgcd(self,other):
        p = self.p
        sexp = self.exponent
        oexp = other.exponent
        if sexp < oexp:
            a = ZZ(self.x)
            b = ZZ(other.x * (p ** (oexp-sexp)))
            exp = sexp
        else:
            a = ZZ(self.x * (p ** (sexp-oexp)))
            b = ZZ(other.x)
            exp = oexp
        d, u, v = a.xgcd(b)
        if self._valuation is None or other._valuation is None:
            val = None
        else:
            val = min(self._valuation, other._valuation)
        d = self.__class__(p, d, exp, valuation=val)
        u = self.__class__(p, u)
        v = self.__class__(p, v)
        return d, u, v

    def value(self):
        return (self.p ** self.exponent) * self.x

    def list(self, prec):
        if self.x not in ZZ:
            raise NotImplementedError
        val = self.valuation()
        p = self.p
        x = ZZ(self.x * p**(self.exponent - val))
        l = [ ]; i = val
        while i < prec:
            x, digit = x.quo_rem(p)
            l.append(digit)
            i += 1
        return l


# Class PrecisionLattice
########################

def list_of_padics(elements):
    from sage.rings.padics.padic_lattice_element import pAdicLatticeElement
    if isinstance(elements, pAdicLatticeElement):
        return [ weakref.ref(elements) ]
    if not isinstance(elements, list):
        elements = list(elements)
    ans = [ ]
    for x in elements:
        ans += list_of_padics(x)
    return ans

def format_history(tme, status, timings):
    status = ''.join(status)
    if timings:
        if tme < 0.000001:
            s = "   ---   "
        elif tme >= 10:
            s = "  >= 10s "
        else:
            s = "%.6fs" % tme
        return s + "  " + status
    else:
        return status

class PrecisionLattice(UniqueRepresentation, SageObject):
    # Internal variables:
    #  . self._cap
    #    a cap for the (working) precision
    #    meaning that the precision lattice always contains p^(self._cap)
    #  . self._elements
    #    list of weak references of elements in this parent
    #  . self._lattice
    #    an upper triangular matrix over ZZ representing the 
    #    lattice of precision
    #    (its columns are indexed by self._elements)
    #  . self._absolute_precisions
    def __init__(self, p, label):
        self._p = p
        self._label = label
        self._elements = [ ]
        self._lattice = { }
        self._absolute_precisions = { }
        self._marked_for_deletion = [ ]
        self._approx_zero = pRational(p, ZZ(0))
        # History
        self._history_init = None
        self._history = None

    def label(self):
        return self._label

    def _repr_(self):
        n = len(self._elements)
        if self._label is None:
            if n > 1:
                return "Precision Lattice on %s objects" % len(self._elements)
            else:
                return "Precision Lattice on %s object" % len(self._elements)
        else:
            if n > 1:
                return "Precision Lattice on %s objects (label: %s)" % (len(self._elements), self._label)
            else:
                return "Precision Lattice on %s object (label: %s)" % (len(self._elements), self._label)

    def prime(self):
        return self._p

    def reduce(self, index=0, partial=False):
        n = len(self._elements)
        if index >= n-1:
            return
        if partial:
            # Partial reduction
            # Cost: O(m^2) with m = n-index
            tme = walltime()
            diffval = (n-index) * [0]
            for j in range(n-1, index, -1):
                col = self._lattice[self._elements[j]]
                prec = col[j].valuation() - diffval[j-index]
                for i in range(index,j):
                    col[i] = col[i].reduce(prec)
                    col[i].normalize()  # seems to be faster then
                    dval = col[i].valuation() - prec
                    if dval < diffval[i-index]:
                        diffval[i-index] = dval
            # We update history
            if self._history is not None:
                self._history.append(('partial reduce', index, walltime(tme)))
        else:
            # Full Hermite reduction
            # Cost: O(m^3) with m = n-index
            tme = walltime()
            for j in range(index+1, n):
                # In what follows, we assume that col[j] is a power of p
                col = self._lattice[self._elements[j]]
                valpivot = col[j].valuation()
                for i in range(index, j):
                    reduced = col[i].reduce(valpivot)
                    scalar = (col[i] - reduced) >> valpivot
                    if scalar.is_zero(): continue
                    col[i] = reduced
                    col[i].normalize()
                    for j2 in range(j+1, n):
                        col2 = self._lattice[self._elements[j2]]
                        col2[i] -= scalar*col2[i]
                        col2[i].normalize()
            # We update history
            if self._history is not None:
                self._history.append(('full reduce', index, walltime(tme)))

    def new_element(self, x, dx, bigoh, dx_mode='linear_combinaison'):
        # First we delete some elements marked for deletion
        if self._marked_for_deletion:
            self.del_elements(thresold=50)

        # Then we add the new element
        tme = walltime()
        p = self._p
        n = len(self._elements)
        x_ref = weakref.ref(x, self.mark_for_deletion)
        self._elements.append(x_ref)
        col = n * [self._approx_zero]
        if dx_mode == 'linear_combinaison':
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                c = self._lattice[ref]
                for i in range(len(c)):
                    col[i] += scalar * c[i]
        elif dx_mode == 'values':
            for elt,scalar in dx:
                ref = weakref.ref(elt)
                if not isinstance(scalar, pRational):
                    scalar = pRational(p, scalar)
                i = len(self._lattice[ref]) - 1
                col[i] = scalar
        else:
            raise ValueError("dx_mode must be either 'linear_combinaison' or 'values'")
        for i in range(n):
            col[i] = col[i].reduce(bigoh)
        col.append(pRational(p, ZZ(1), bigoh))
        self._lattice[x_ref] = col

        # We compute the absolute precision of the new element and cache it
        self._absolute_precisions[x_ref] = min([ c.valuation() for c in col ])

        # We update history
        if self._history is not None:
            self._history.append(('add', None, walltime(tme)))

    def _set_precision(self, x, column={}):
        p = self._p
        x_ref = weakref.ref(x)
        index = len(self._lattice[x_ref]) - 1
        n = len(self._elements)
        col = n * [self._approx_zero]
        self._lattice[x_ref] = col
        self._absolute_precisions[x_ref] = min([ c.valuation() for c in col ])

    def mark_for_deletion(self, ref):
        tme = walltime()
        try:
            index = len(self._lattice[ref]) - 1
        except IndexError:
            return
        self._marked_for_deletion.append(index)
        if self._history is not None:
            self._history.append(('mark', index, walltime(tme)))

    def del_elements(self, thresold=None):
        p = self._p
        n = len(self._elements)

        self._marked_for_deletion.sort(reverse=True)
        count = 0
        for index in self._marked_for_deletion:
            if thresold is not None and index < n - thresold: break
            n -= 1; count += 1

            tme = walltime()
            del self._lattice[self._elements[index]]
            del self._elements[index]

            # Now, we echelonize
            for i in range(index,n):
                col = self._lattice[self._elements[i]]
                vali = col[i].valuation()
                valj = col[i+1].valuation()
                d, u, v = col[i].xgcd(col[i+1])
                up, vp = col[i+1]/d, col[i]/d
                col[i] = d
                del col[i+1]
                for j in range(i+1,n):
                    col = self._lattice[self._elements[j]]
                    col[i], col[i+1] = u*col[i] + v*col[i+1], up*col[i] - vp*col[i+1]

            # We update history
            if self._history is not None:
                self._history.append(('del', index, walltime(tme)))

            # And we reduce a bit
            # (we do not perform a complete reduction because it is costly)
            self.reduce(index, partial=True)

        del self._marked_for_deletion[:count]

    def lift_to_precision(self, x, prec):
        ref = weakref.ref(x)
        col = self._lattice[ref]
        n = len(self._elements)

        rows_by_val = { }
        for i in range(len(col)):
            v = col[i].valuation()
            if v >= prec: continue
            if rows_by_val.has_key(v):
                rows_by_val[v].append(i)
            else:
                rows_by_val[v] = [i]
        vals = rows_by_val.keys()
        vals.sort()
        vals.append(prec)

        for t in range(len(vals)-1):
            v, w = vals[t], vals[t+1]
            rows = rows_by_val[v]
            piv = max(rows)
            alpha = col[piv].unit_part()
            for i in rows:
                if i == piv: continue
                # We clear the entry on the i-th line
                beta = col[i].unit_part()
                for j in range(piv,n):
                    col_cur = self._lattice[self._elements[j]]
                    col_cur[i] = alpha*col_cur[i] - beta*col_cur[piv]
            # We rescale the piv-th line
            for j in range(piv,n):
                col_cur = self._lattice[self._elements[j]]
                col_cur[piv] = col_cur[piv] << (w-v)
            # Now the entry on the piv-th line has valuation w
            # We update the dictionary accordingly
            if w < prec:
                rows_by_val[w].append(piv)

    def precision_absolute(self, x):
        ref = weakref.ref(x)
        return self._absolute_precisions[ref]

    def precision_lattice(self, elements=None, echelon=True):
        if elements is None:
            elements = self._elements
        else:
            elements = list_of_padics(elements)
        n = len(self._elements)
        rows = [ ]; val = 0
        for ref in elements:
            col = self._lattice[ref]
            row = [ x.value() for x in col ]
            valcol = min([ x.valuation() for x in col ])
            if valcol < val: val = valcol
            row += (n-len(row)) * [ZZ(0)]
            rows.append(row)
        from sage.matrix.constructor import matrix
        M = matrix(rows).transpose()
        if val < 0:
            M *= self._p ** (-val)
        if echelon:
            M = M.change_ring(ZZ)
            M.echelonize()
            n = len(elements)
            M = M.submatrix(0,0,n,n)
        if val < 0:
            M *= self._p ** val
        return M

    def number_of_diffused_digits(self, elements=None):
        M = self.precision_lattice(elements)
        n = M.nrows()
        p = self._p
        diffused = 0
        for j in range(n):
            val = minval = M[j,j].valuation(p)
            for i in range(j):
                v = M[i,j].valuation(p)
                if v < minval: minval = v
            diffused += val - minval
        return diffused

    def number_of_tracked_elements(self):
        return len(self._elements)

    def tracked_elements(self, values=True):
        if values:
            return [ ref() for ref in self._elements ]
        else:
            return list(self._elements)


    # History

    def history_enable(self):
        if self._history is None:
            self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
            self._history = [ ]

    def history_disable(self):
        self._history = self._history_init = None

    def history_clear(self):
        if self._history is None:
            raise ValueError("History is not tracked")
        self._history_init = ( len(self._elements), list(self._marked_for_deletion) )
        self._history = [ ]

    def history(self, compact=True, include_reduce=False, timings=True, output_type='asciiart'):
        if self._history is None:
            raise ValueError("History is not tracked")
        total_time = 0
        if output_type == 'asciiart':
            # Legend:
            #  o : tracked element
            #  ~ : element marked for deletion
            #  r : partial reduction
            #  R : full Hermite reduction
            (n, mark) = self._history_init
            status = n*['o']
            for index in mark:
                status[index] = '~'
            hist = [ ]; oldevent = ''
            for (event, index, tme) in self._history:
                if not include_reduce and (event == 'partial reduce' or event == 'full reduce'):
                    continue
                if not compact or oldevent != event:
                    hist.append(format_history(total_time, status, timings))
                    total_time = 0
                    oldevent = event
                total_time += tme
                if event == 'add':
                    status.append('o')
                elif event == 'mark':
                    status[index] = '~'
                elif event == 'del':
                    del status[index]
                elif event == 'partial reduce':
                    status_red = status[:index] + (len(status) - index)*['r']
                    hist.append(format_history(0, status_red, timings))
                elif event == 'full reduce':
                    status_red = status[:index] + (len(status) - index)*['R']
                    hist.append(format_history(0, status_red, timings))
                else:
                    raise RuntimeError("event not known in the history")
            hist.append(format_history(total_time, status, timings))
            return '\n'.join(hist)
        else:
            raise NotImplementedError

    def timings(self, action=None):
        if self._history is None:
            raise ValueError("History is not tracked")
        tme_by_event = { }
        for (event, _, tme) in self._history:
            if tme_by_event.has_key(event):
                tme_by_event[event] += tme
            else:
                tme_by_event[event] = tme
        if action is None:
            return tme_by_event
        if tme_by_event.has_key(action):
            return tme_by_event[action]
        else:
            return 0 
