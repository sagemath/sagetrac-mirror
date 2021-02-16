"""Converters between Sage types and `~gappy.gapobj.GapObj types."""


import copyreg
import itertools

from sage.combinat.permutation import Permutation
from sage.groups.perm_gps.permgroup_element cimport PermutationGroupElement
from sage.libs.gmp.mpz cimport mpz_set_si
from sage.libs.gmp.mpq cimport mpq_numref, mpq_denref
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.all import ZZ, QQ, RDF
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.structure.coerce cimport coercion_model as cm
from sage.structure.sage_object cimport SageObject
from .libgap import libgap

from gappy.core cimport Gap
from gappy.gapobj cimport (GapObj, GapInteger, GapFloat, GapIntegerMod,
        GapFiniteField, GapCyclotomic, GapRational, GapRing, GapBoolean,
        GapString, GapList, GapPermutation, GapRecord)


# TODO: It might be good to implement a "fast lane" for some built-in
# Sage types like Integer (e.g. currently Sage Integers are first
# converted to Python ints, and then to GAP Integers, whereas it would
# be much faster to convert Sage Integers directly to GAP Integers since
# they are both basically mpz_t limbs under the hood).
@libgap.convert_from(SageObject)
def sageobject_to_gapobj(Gap gap, SageObject obj):
    r"""
    gappy converter for converting generic `.SageObject`\s to their
    corresponding `~gappy.gapobj.GapObj` if any.

    This implements the libgap conversion functions already documented for
    `.SageObject`\s: `.SageObject._libgap_` and `.SageObject._libgap_init_`.
    """

    # NOTE: In the default implementation of SageObject._libgap_ it defers
    # to _libgap_init_, so we just need to try calling _libgap_
    try:
        ret = obj._libgap_(gap)
    except TypeError:
        # Support older code where _libgap_ does not take an argument
        ret = obj._libgap_()
    if isinstance(ret, gap.supported_builtins):
        return gap(ret)
    elif isinstance(ret, GapObj):
        return ret
    else:
        raise RuntimeError(
            f'{type(obj).__name__}._libgap_ returned something that cannot '
            f'be converted to a GAP object: {ret!r} ({type(ret).__name__})')


# Pickling for GapObjs using a Sage-specific protocol which tries to pickle
# the objects' Sage representations where possible.
def _reduce_gapobj(obj):
   """
   Attempt to pickle GAP objects from gappy.

   This is inspired in part by ``sage.interfaces.interface.Interface._reduce``,
   though for a fallback we use ``str(self)`` instead of ``repr(self)``, since
   the former is equivalent in the libgap interface to the latter in the
   pexpect interface.

   TESTS:

   This workaround was motivated in particular by this example from the
   permutation groups implementation::

       sage: CC = libgap.eval('ConjugacyClass(SymmetricGroup([ 1 .. 5 ]), (1,2)(3,4))')
       sage: CC.sage()
       Traceback (most recent call last):
       ...
       NotImplementedError: cannot construct equivalent Sage object
       sage: libgap.eval(str(CC))
       (1,2)(3,4)^G
       sage: loads(dumps(CC))
       (1,2)(3,4)^G
   """

   if obj.is_string():
       elem = repr(obj.sage())
   try:
       elem = obj.sage()
   except NotImplementedError:
       elem = str(obj)

   return (_construct_gapobj, (elem,))


def _construct_gapobj(reduced):
   """
   Currently just used for unpickling; equivalent to calling ``libgap(elem)``
   to convert a Sage object to a `~gappy.gapobj.GapObj` where possible.
   """

   if isinstance(reduced, str):
       return libgap.eval(reduced)

   return libgap(reduced)


for cls in itertools.chain([GapObj], GapObj.__subclasses__()):
    copyreg.pickle(cls, _reduce_gapobj, _construct_gapobj)

del cls


@GapObj.convert_to('sage')
def gapobj_to_sage(obj):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapObj`.

    EXAMPLES::

        sage: libgap(1).sage()
        1
        sage: type(_)
        <type 'sage.rings.integer.Integer'>

        sage: libgap(3/7).sage()
        3/7
        sage: type(_)
        <type 'sage.rings.rational.Rational'>

        sage: libgap.eval('5 + 7*E(3)').sage()
        7*zeta3 + 5

        sage: libgap(Infinity).sage()
        +Infinity
        sage: libgap(-Infinity).sage()
        -Infinity

        sage: libgap(True).sage()
        True
        sage: libgap(False).sage()
        False
        sage: type(_)
        <... 'bool'>

        sage: libgap('this is a string').sage()
        'this is a string'
        sage: type(_)
        <... 'str'>

        sage: x = libgap.Integers.Indeterminate("x")

        sage: p = x^2 - 2*x + 3
        sage: p.sage()
        x^2 - 2*x + 3
        sage: p.sage().parent()
        Univariate Polynomial Ring in x over Integer Ring

        sage: p = x^-2 + 3*x
        sage: p.sage()
        x^-2 + 3*x
        sage: p.sage().parent()
        Univariate Laurent Polynomial Ring in x over Integer Ring

        sage: p = (3 * x^2 + x) / (x^2 - 2)
        sage: p.sage()
        (3*x^2 + x)/(x^2 - 2)
        sage: p.sage().parent()
        Fraction Field of Univariate Polynomial Ring in x over Integer Ring

    TESTS:

    Check :trac:`30496`::

        sage: x = libgap.Integers.Indeterminate("x")

        sage: p = x^2 - 2*x
        sage: p.sage()
        x^2 - 2*x
    """

    if obj.IsInfinity():
        from sage.rings.infinity import Infinity
        return Infinity

    elif obj.IsNegInfinity():
        from sage.rings.infinity import Infinity
        return -Infinity

    elif obj.IsUnivariateRationalFunction():
        var = obj.IndeterminateOfUnivariateRationalFunction().String()
        var = var.sage()
        num, den, val = obj.CoefficientsOfUnivariateRationalFunction()
        num = num.sage()
        den = den.sage()
        val = val.sage()
        base_ring = cm.common_parent(*(num + den))

        if obj.IsUnivariatePolynomial():
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(base_ring, var)
            x = R.gen()
            return x**val * R(num)

        elif obj.IsLaurentPolynomial():
            from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
            R = LaurentPolynomialRing(base_ring, var)
            x = R.gen()
            return x**val * R(num)

        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            R = PolynomialRing(base_ring, var)
            x = R.gen()
            return x**val * R(num) / R(den)

    elif obj.IsList():
        # May be a list-like collection of some other type of GapElements
        # that we can convert
        return [item.sage() for item in obj.AsList()]

    raise NotImplementedError('cannot construct equivalent Sage object')


@GapInteger.convert_to('sage')
def gapinteger_to_sage(obj, ring=None, z=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapInteger`

    - ``ring`` -- Integer ring or ``None`` (default). If not
      specified, a the default Sage integer ring is used.

    - ``z`` -- An existing Sage Integer on which to set a value.

    OUTPUT:

    A Sage integer

    EXAMPLES::

        sage: libgap([ 1, 3, 4 ]).sage()
        [1, 3, 4]
        sage: all( x in ZZ for x in _ )
        True

        sage: libgap(132).sage(ring=IntegerModRing(13))
        2
        sage: parent(_)
        Ring of integers modulo 13

    TESTS::

        sage: large = libgap.eval('2^130');  large
        1361129467683753853853498429727072845824
        sage: large.sage()
        1361129467683753853853498429727072845824

        sage: huge = libgap.eval('10^9999');  huge     # gap abbreviates very long ints
        <integer 100...000 (10000 digits)>
        sage: huge.sage().ndigits()
        10000
    """
    cdef Integer zz
    cdef GapInteger gz

    if ring is None:
        ring = ZZ

    if z is None:
        zz = Integer.__new__(Integer)
    else:
        zz = z

    gz = <GapInteger>obj

    if gz.is_C_int():
        mpz_set_si(zz.value, gz.to_C_int())
    else:
        gz.to_mpz(zz.value)

    if ring is ZZ:
        return zz
    else:
        return ring(zz)


@GapFloat.convert_to('sage')
def gapfloat_to_sage(obj, ring=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapFloat`

    - ``ring`` -- a floating point field or ``None`` (default). If not
      specified, the default Sage ``RDF`` is used.

    OUTPUT:

    A Sage double precision floating point number

    EXAMPLES::

        sage: a = libgap.eval("Float(3.25)").sage()
        sage: a
        3.25
        sage: parent(a)
        Real Double Field
    """
    if ring is None:
        ring = RDF
    return ring(float(obj))


@GapIntegerMod.convert_to('sage')
def gapintegermod_to_sage(obj, ring=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapIntegerMod`

    INPUT:

    - ``ring`` -- Sage integer mod ring or ``None`` (default). If not
      specified, a suitable integer mod ringa is used automatically.

    OUTPUT:

    A Sage integer modulo another integer.

    EXAMPLES::

        sage: n = libgap.eval('One(ZmodnZ(123)) * 13')
        sage: n.sage()
        13
        sage: parent(_)
        Ring of integers modulo 123
    """
    if ring is None:
        # ring = obj.DefaultRing().sage()
        characteristic = obj.Characteristic().sage()
        ring = ZZ.quotient_ring(characteristic)
    return obj.lift().sage(ring=ring)


@GapFiniteField.convert_to('sage')
def gapfinitefield_to_sage(obj, ring=None, var='a'):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapFiniteField`.

    INPUT:

    - ``ring`` -- a Sage finite field or ``None`` (default). The field to
      return ``self`` in. If not specified, a suitable finite field will be
      constructed.

    OUTPUT:

    An Sage finite field element. The isomorphism is chosen such that the Gap
    ``PrimitiveRoot()`` maps to the Sage
    :meth:`~sage.rings.finite_rings.finite_field_prime_modn.multiplicative_generator`.

    EXAMPLES::

        sage: n = libgap.eval('Z(25)^2')
        sage: n.sage()
        a + 3
        sage: parent(_)
        Finite Field in a of size 5^2

        sage: n.sage(ring=GF(5))
        Traceback (most recent call last):
        ...
        ValueError: the given ring is incompatible ...

    TESTS::

        sage: n = libgap.eval('Z(2^4)^2 + Z(2^4)^1 + Z(2^4)^0')
        sage: n
        Z(2^2)^2
        sage: n.sage()
        a + 1
        sage: parent(_)
        Finite Field in a of size 2^2
        sage: n.sage(ring=ZZ)
        Traceback (most recent call last):
        ...
        ValueError: the given ring is incompatible ...
        sage: n.sage(ring=CC)
        Traceback (most recent call last):
        ...
        ValueError: the given ring is incompatible ...
        sage: n.sage(ring=GF(5))
        Traceback (most recent call last):
        ...
        ValueError: the given ring is incompatible ...
        sage: n.sage(ring=GF(2^3))
        Traceback (most recent call last):
        ...
        ValueError: the given ring is incompatible ...
        sage: n.sage(ring=GF(2^2, 'a'))
        a + 1
        sage: n.sage(ring=GF(2^4, 'a'))
        a^2 + a + 1
        sage: n.sage(ring=GF(2^8, 'a'))
        a^7 + a^6 + a^4 + a^2 + a + 1

    Check that :trac:`23153` is fixed::

        sage: n = libgap.eval('Z(2^4)^2 + Z(2^4)^1 + Z(2^4)^0')
        sage: n.sage(ring=GF(2^4, 'a'))
        a^2 + a + 1
    """
    deg = obj.DegreeFFE().sage()
    char = obj.Characteristic().sage()
    if ring is None:
        from sage.rings.finite_rings.finite_field_constructor import GF
        ring = GF(char**deg, name=var)
    elif not (ring.is_field() and ring.is_finite() and \
              ring.characteristic() == char and ring.degree() % deg == 0):
        raise ValueError(f'the given ring is incompatible (must be a '
                         f'finite field of characteristic {char} and degree '
                         f'divisible by {deg})')

    if obj.IsOne():
        return ring.one()
    if deg == 1 and char == ring.characteristic():
        return ring(obj.lift().sage())
    else:
        gap_field = obj.parent(ring)
        exp = obj.LogFFE(gap_field.PrimitiveRoot())
        return ring.multiplicative_generator() ** exp.sage()


@GapCyclotomic.convert_to('sage')
def gapcyclotomic_to_sage(obj, ring=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapCyclotomic`.

    INPUT:

    - ``ring`` -- a Sage cyclotomic field or ``None`` (default). If not
      specified, a suitable minimal cyclotomic field will be constructed.

    OUTPUT:

    A Sage cyclotomic field element.

    EXAMPLES::

        sage: n = libgap.eval('E(3)')
        sage: n.sage()
        zeta3
        sage: parent(_)
        Cyclotomic Field of order 3 and degree 2

        sage: n.sage(ring=CyclotomicField(6))
        zeta6 - 1

        sage: libgap.E(3).sage(ring=CyclotomicField(3))
        zeta3
        sage: libgap.E(3).sage(ring=CyclotomicField(6))
        zeta6 - 1

    TESTS:

    Check that :trac:`15204` is fixed::

        sage: libgap.E(3).sage(ring=UniversalCyclotomicField())
        E(3)
        sage: libgap.E(3).sage(ring=CC)
        -0.500000000000000 + 0.866025403784439*I
    """
    if ring is None:
        conductor = obj.Conductor()
        from sage.rings.number_field.number_field import CyclotomicField
        ring = CyclotomicField(conductor.sage())
    else:
        try:
            conductor = ring._n()
        except AttributeError:
            from sage.rings.number_field.number_field import CyclotomicField
            conductor = obj.Conductor()
            cf = CyclotomicField(conductor.sage())
            return ring(cf(obj.CoeffsCyc(conductor).sage()))
    coeff = obj.CoeffsCyc(conductor).sage()
    return ring(coeff)


@GapRational.convert_to('sage')
def gaprational_to_sage(obj, ring=None, q=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapRational`.

    INPUT:

    - ``ring`` -- the Sage rational ring or ``None`` (default). If not
      specified, the rational ring is used automatically.

    - ``q`` -- An existing Sage Rational on which to set a value.

    OUTPUT:

    A Sage rational number.

    EXAMPLES::

        sage: r = libgap(123/456);  r
        41/152
        sage: type(_)
        <class 'gappy.gapobj.GapRational'>
        sage: r.sage()
        41/152
        sage: type(_)
        <type 'sage.rings.rational.Rational'>
    """

    cdef Rational qq
    cdef GapInteger num, den

    if ring is None:
        ring = QQ

    num = <GapInteger>(obj.NumeratorRat())
    den = <GapInteger>(obj.DenominatorRat())

    if ring is QQ:
        if q is None:
            qq = Rational.__new__(Rational)
        else:
            qq = <Rational>q

        if num.is_C_int():
            mpz_set_si(mpq_numref(qq.value), num.to_C_int())
        else:
            num.to_mpz(mpq_numref(qq.value))
        if den.is_C_int():
            mpz_set_si(mpq_denref(qq.value), den.to_C_int())
        else:
            den.to_mpz(mpq_denref(qq.value))

        return qq
    else:
        return num.sage(ring=ring) / den.sage(ring=ring)


@GapRing.convert_to('sage')
def gapring_to_sage(obj, var='a'):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapRing`.

    INPUT:

    - ``var`` -- the variable name to label the generator of finite fields

    OUTPUT:

    A Sage ring.

    EXAMPLES::

        sage: libgap.eval('Integers').sage()
        Integer Ring

        sage: libgap.eval('Rationals').sage()
        Rational Field

        sage: libgap.eval('ZmodnZ(15)').sage()
        Ring of integers modulo 15

        sage: libgap.GF(3,2).sage(var='A')
        Finite Field in A of size 3^2

        sage: libgap.CyclotomicField(6).sage()
        Cyclotomic Field of order 3 and degree 2

        sage: libgap(QQ['x','y']).sage()
        Multivariate Polynomial Ring in x, y over Rational Field
    """

    if obj.IsField():
        if obj.IsRationals():
            return ZZ.fraction_field()
        if obj.IsCyclotomicField():
            conductor = obj.Conductor()
            from sage.rings.number_field.number_field import CyclotomicField
            return CyclotomicField(conductor.sage())
        if obj.IsFinite():
            size = obj.Size().sage()
            from sage.rings.finite_rings.finite_field_constructor import GF
            return GF(size, name=var)
    else:
        if obj.IsIntegers():
            return ZZ
        if obj.IsFinite():
            characteristic = obj.Characteristic().sage()
            return ZZ.quotient_ring(characteristic)
        if obj.IsPolynomialRing():
            base_ring = obj.CoefficientsRing().sage()
            vars = [x.String().sage()
                    for x in obj.IndeterminatesOfPolynomialRing()]
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            return PolynomialRing(base_ring, vars)

    raise NotImplementedError('cannot convert GAP ring to Sage')


@GapBoolean.convert_to('sage')
def gapboolean_to_sage(obj):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapBoolean`

    OUTPUT:

    A Python boolean if the values is either true or false. GAP
    booleans can have the third value ``Fail``, in which case a
    ``ValueError`` is raised.

    EXAMPLES::

        sage: b = libgap.eval('true');  b
        true
        sage: type(_)
        <class 'gappy.gapobj.GapBoolean'>
        sage: b.sage()
        True
        sage: type(_)
        <... 'bool'>

        sage: libgap.eval('fail')
        fail
        sage: _.sage()
        Traceback (most recent call last):
        ...
        ValueError: the GAP boolean value "fail" cannot be represented in Sage
    """
    if repr(obj) == 'fail':
        raise ValueError(
            'the GAP boolean value "fail" cannot be represented in Sage')
    return obj.__bool__()


@GapString.convert_to('sage')
def gapstring_to_string(obj):
    """
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapString`
    """

    return obj.__str__()


@GapList.convert_to('sage')
def gaplist_to_sage(obj, **kwds):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapList`

    OUTPUT:

    A Python list.

    EXAMPLES::

        sage: libgap([ 1, 3, 4 ]).sage()
        [1, 3, 4]
        sage: all( x in ZZ for x in _ )
        True
    """
    return [x.sage(**kwds) for x in obj]


# For compatibility with the old libgap interface, GapLists also have a
# conversion named simply 'matrix' which converts them to a Sage matrix.
@GapList.convert_to('matrix')
def gaplist_to_matrix(self, ring=None):
    """
    Return the list as a matrix.

    GAP does not have a special matrix data type, they are just lists of lists.
    This function converts a GAP list of lists to a Sage matrix.

    OUTPUT:

    A Sage matrix.

    EXAMPLES::

        sage: F = libgap.GF(4)
        sage: a = F.PrimitiveElement()
        sage: m = libgap([[a,a^0],[0*a,a^2]]); m
        [ [ Z(2^2), Z(2)^0 ],
          [ 0*Z(2), Z(2^2)^2 ] ]
        sage: m.IsMatrix()
        true
        sage: matrix(m)
        [    a     1]
        [    0 a + 1]
        sage: matrix(GF(4,'B'), m)
        [    B     1]
        [    0 B + 1]

        sage: M = libgap.eval('SL(2,GF(5))').GeneratorsOfGroup()[1]
        sage: type(M)
        <type 'gappy.gapobj.GapList'>
        sage: M[0][0]
        Z(5)^2
        sage: M.IsMatrix()
        true
        sage: M.matrix()
        [4 1]
        [4 0]
    """

    if not self.IsMatrix():
        raise ValueError('not a GAP matrix')

    if not self.IsRectangularTable():
        raise ValueError('not a rectangular list of lists')

    n, m = self.DimensionsMat()
    entries = self.Flat()

    if ring is None:
        ring = entries.DefaultRing().sage()

    MS = MatrixSpace(ring, n, m)
    return MS([x.sage(ring=ring) for x in entries])


@GapPermutation.convert_to('sage')
def gappermutation_to_sage(obj, parent=None):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapPermutation`

    If the permutation group is given as parent, this method is *much* faster.

    EXAMPLES::

        sage: perm_gap = libgap.eval('(1,5,2)(4,3,8)');  perm_gap
        (1,5,2)(3,8,4)
        sage: perm_gap.sage()
        [5, 1, 8, 3, 2, 6, 7, 4]
        sage: type(_)
        <class 'sage.combinat.permutation.StandardPermutations_all_with_category.element_class'>
        sage: perm_gap.sage(PermutationGroup([(1,2),(1,2,3,4,5,6,7,8)]))
        (1,5,2)(3,8,4)
        sage: type(_)
        <type 'sage.groups.perm_gps.permgroup_element.PermutationGroupElement'>
    """
    cdef PermutationGroupElement one_c

    lst = obj.ListPerm()

    if parent is None:
        return Permutation(lst.sage(), check_input=False)
    else:
        return parent.one()._generate_new_GAP(lst)


@GapRecord.convert_to('sage')
def sage(obj):
    r"""
    Return the Sage equivalent of the :class:`~gappy.gapobj.GapRecord`

    EXAMPLES::

        sage: libgap.eval('rec(a:=1, b:=2)').sage()
        {'a': 1, 'b': 2}
        sage: all( isinstance(key,str) and val in ZZ for key,val in _.items() )
        True

        sage: rec = libgap.eval('rec(a:=123, b:=456, Sym3:=SymmetricGroup(3))')
        sage: rec.sage()
        {'Sym3': NotImplementedError('cannot construct equivalent Sage object'...),
         'a': 123,
         'b': 456}
    """
    result = {}
    for key, val in obj:
        try:
            val = val.sage()
        except Exception as ex:
            val = ex
        result[str(key)] = val
    return result
