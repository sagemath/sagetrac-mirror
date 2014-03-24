r"""
Databases for translation surfaces.

This file contain different databases (with different implementations) for
algorithms related to translation surfaces:

- structure of Strebel differentials for quadratic differentials in order to
  differentiate
- database of separatrix and cylinder diagrams up to isomorphism
- database of volume of connected components of Abelian strata
"""
import os
import sage.misc.misc
from sage.rings.integer import Integer
FLAT_DB_HOME='%s/data/flat_surfaces'%sage.misc.misc.SAGE_ROOT

def line_count(filename):
    r"""
    Returns the number of lines in the file whose name is filename.
    """
    f = open(filename)
    lines = 0
    read_f = f.read # loop optimization

    buf = read_f(0X100000)
    while buf:
        lines += buf.count('\n')
        buf = read_f(0X100000)  # 1024 x 1024

    f.close()

    return Integer(lines)


class GenericRepertoryDatabase:
    r"""
    Database that consists of a list of files in a repertory.
    """
    default_name = "generic_db"
    default_path = sage.misc.misc.SAGE_TMP

    def __init__(self, path=None, name=None, read_only=True):
        r"""
        INPUT:

        - ``path`` - string (default is SAGE_TMP) - path to the database. If the
          repertory does not exists, it is created.

        -  ``name`` - string (default is 'generic_db') - name of the repertory
           that will contain files of the database.

        - ``read_only`` - boolean
        """
        if path is None:
            path = self.default_path

        if name is None:
            name = self.default_name

        if not os.path.isdir(FLAT_DB_HOME):
            try:
                os.mkdir(FLAT_DB_HOME)
            except OSError:
                raise ValueError, "not able to create the database"

        full_path = os.path.join(path,name)

        if not os.path.isdir(full_path):
            try:
                os.mkdir(full_path)
            except OSError:
                raise ValueError, "not able to create the database"

        self.path = os.path.abspath(full_path)
        self.read_only = read_only

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from sage.databases.flat_surfaces import CylinderDiagrams
            sage: from sage.misc.misc import SAGE_TMP
            sage: import os

            sage: rep1 = os.path.join(SAGE_TMP, "rep1")
            sage: rep2 = os.path.join(SAGE_TMP, "rep2")

            sage: C1 = CylinderDiagrams(rep1)
            sage: C2 = CylinderDiagrams(rep1)
            sage: C3 = CylinderDiagrams(rep2)
            sage: C1 == C1
            True
            sage: C1 == C2
            True
            sage: C1 == C3
            False
        """
        return isinstance(other, CylinderDiagrams) and other.path == self.path

    def __ne__(self, other):
        r"""
        TESTS::

            sage: from sage.databases.flat_surfaces import CylinderDiagrams
            sage: from sage.misc.misc import SAGE_TMP
            sage: import os

            sage: rep1 = os.path.join(SAGE_TMP, "A")
            sage: rep2 = os.path.join(SAGE_TMP, "B")


            sage: C1 = CylinderDiagrams(rep1)
            sage: C2 = CylinderDiagrams(rep1)
            sage: C3 = CylinderDiagrams(rep2)
            sage: C1 != C1
            False
            sage: C1 != C2
            False
            sage: C1 != C3
            True
        """
        return not self.__eq__(other)

    def clean(self):
        r"""
        Clean the database.

        EXAMPLES::

            sage: from sage.databases.flat_surfaces import CylinderDiagrams
            sage: from sage.misc.misc import SAGE_TMP

            sage: rep = os.path.join(SAGE_TMP, "cylinder_diagrams")
            sage: C = CylinderDiagrams(rep)
            sage: C.update(AbelianStratum(4))
            sage: import os
            sage: os.listdir(rep)
            ['4-hyp-2', '4-hyp-1', '4-hyp-3', '4-odd-2', '4-odd-3', '4-odd-1']
            sage: C.clean()
            sage: os.listdir(rep)
            []
        """
        assert not self.read_only
        for filename in os.listdir(self.path):
            os.remove(os.path.join(self.path,filename))

class IrregularComponentTwins(GenericRepertoryDatabase):
    r"""
    Twin data of generalized permutation of irregular components of strata of
    Abelian differentials.
    """
    default_name = "generalized_permutation_twins"
    default_path = FLAT_DB_HOME

    def __repr__(self):
        r"""
        String representation

        TESTS::

            sage: from sage.database.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: D.__repr__()  # random
            'Database of twins of irregular components at /path/to/database'
        """
        return "Database of twins of irregular components at %s"%self.path

    def filename(self, stratum):
        r"""
        Returns the name of the file for the given component.

        EXAMPLES::

            sage: from sage.database.flat_surfaces import IrregularComponentsTwins
            sage: D = IrregularComponentTwins()
            sage: D.filename(QuadraticStratum(12))
            'twins-12-irr'
            sage: D.filename(QuadraticStratum(2,2))
            Traceback (most recent call last):
            ...
            AssertionError: the stratum has no irregular component
        """
        from sage.dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        assert isinstance(stratum, QuadraticStratum)
        assert stratum.has_regular_and_irregular_components(), "the stratum has no irregular component"
        return ('twins-' +
                '_'.join(str(z) for z in stratum.zeros(poles=False)) +
                '_p' * stratum.nb_poles() +
                '-irr')

    def has_stratum(self, stratum):
        r"""
        Test whether the component is in the database.

        EXAMPLES::

            sage: from sage.database.flat_surfaces import IrregularComponentsTwins
            sage: D = IrregularComponentTwins()
            sage: D.has_stratum(QuadraticStratum(12))  # optional
            True
        """
        return os.path.isfile(os.path.join(self.path,self.filename(stratum)))

    def update(self, stratum=None):
        r"""
        Update the database with the component comp.
        If comp is None update the whole database.

        The database should be not in read only mode.
        """
        assert not self.read_only

        f = open(os.path.join(selfr.path, self.filename(stratum)))

        p = stratum.irregular_component().permutation_representative()
        res = []
        for q in p.rauzy_diagram(symmetric=True):
            if q.is_cylindric():
                res.append(q._twin)
        res.sort()

        for twin in res:
            f.write(str(twin) + '\n')

    def list_strata(self):
        r"""
        Returns the list of components for which the list of twins is stored.

        EXAMPLES::

            sage: from sage.databases.flat_surfaces import
            IrregularComponentTwins
            sage: G = IrregularComponentTwins()
            sage: G.list_strata()
            [Q_3(9, -1), Q_3(6, 3, -1), Q_4(12), Q_3(3^3, -1), Q_4(6^2), Q_4(9, 3), Q_4(6, 3^2), Q_4(3^4)]
        """
        from sage.dynamics.flat_surfaces.quadratic_strata import QuadraticStratum
        from sage.rings.all import Integer
        s = set()
        for f in os.listdir(self.path):
            if f.startswith('twins-'):
                g = f[6:]
                i = g.index('-')
                comp = g[:i].replace('p','-1')
                s.add(comp)

        return sorted(QuadraticStratum(map(Integer,g.split('_'))) for g in s)

    def get(self, stratum):
        r"""
        Get the list of twins for the stratum of quadratic differentials ``stratum``.

        EXAMPLES::

            sage: from sage.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: l = D.get(QuadraticStratum(6,3,-1))
            sage: l[0][0]
            [(1, 8), (1, 5), (1, 4), (0, 4), (0, 3), (1, 7), (1, 6)]
            sage: l[0][1][(1, 2),
            (1, 3), (1, 0), (1, 1), (0, 2), (0, 1), (0, 6), (0, 5), (0, 0)]
            sage: len(l)
            1634
        """
        assert self.has_stratum(stratum)
        f = open(os.path.join(self.path, self.filename(stratum)))

        res = []
        s = f.readline()
        while s:
            res.append(eval(s))
            s = f.readline()
        return res

    def count(self, stratum):
        r"""
        Returns the number of twins for that stratum.

        EXAMPLES::

            sage: from sage.databases.flat_surfaces import IrregularComponentTwins
            sage: D = IrregularComponentTwins()
            sage: Q = QuadraticStratum(12)
            sage: len(D.get(Q))
            6993
            sage: D.count(Q)
            6993
        """
        return line_count(os.path.join(self.path, self.filename(stratum)))

