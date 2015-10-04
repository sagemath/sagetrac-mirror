"""
Access the PALP database(s) of reflexive lattice polytopes

EXAMPLES::

    sage: from sage.geometry.polyhedron.palp_database import PALPreader
    sage: for lp in PALPreader(2):
    ....:     cone = Cone([(1,r[0],r[1]) for r in lp.vertices()])
    ....:     fan = Fan([cone])
    ....:     X = ToricVariety(fan)
    ....:     ideal = X.affine_algebraic_patch(cone).defining_ideal()
    ....:     print lp.n_vertices(), ideal.hilbert_series()
    3 (-t^2 - 7*t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    3 (-t^2 - t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    3 (t^2 + 6*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    3 (t^2 + 2*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    3 (t^2 + 4*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    4 (-t^2 - 5*t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    4 (-t^2 - 3*t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    4 (t^2 + 2*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    4 (t^2 + 6*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    4 (t^2 + 6*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    4 (t^2 + 2*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    4 (t^2 + 4*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    5 (-t^2 - 3*t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    5 (-t^2 - 5*t - 1)/(t^3 - 3*t^2 + 3*t - 1)
    5 (t^2 + 4*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
    6 (t^2 + 4*t + 1)/(-t^3 + 3*t^2 - 3*t + 1)
"""

from subprocess import Popen, PIPE

from sage.structure.sage_object import SageObject
from sage.matrix.all import matrix
from sage.rings.all import Integer, ZZ

from sage.geometry.polyhedron.ppl_lattice_polytope import LatticePolytope_PPL
from sage.geometry.polyhedron.constructor import Polyhedron



#########################################################################
class PALPreader(SageObject):
    """
    Read PALP database of polytopes.


    INPUT:

    - ``dim`` -- integer. The dimension of the poylhedra

    - ``data_basename`` -- string or ``None`` (default). The directory
      and database base filename (PALP usually uses ``'zzdb'``) name
      containing the PALP database to read. Defaults to the built-in
      database location.

    - ``output`` -- string. How to return the reflexive polyhedron
      data. Allowed values = ``'list'``, ``'Polyhedron'`` (default),
      ``'pointcollection'``, and ``'PPL'``. Case is ignored.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.palp_database import PALPreader
        sage: polygons = PALPreader(2)
        sage: [ (p.n_Vrepresentation(), len(p.integral_points())) for p in polygons ]
        [(3, 4), (3, 10), (3, 5), (3, 9), (3, 7), (4, 6), (4, 8), (4, 9),
         (4, 5), (4, 5), (4, 9), (4, 7), (5, 8), (5, 6), (5, 7), (6, 7)]

        sage: next(iter(PALPreader(2, output='list')))
        [[1, 0], [0, 1], [-1, -1]]
        sage: type(_)
        <type 'list'>

        sage: next(iter(PALPreader(2, output='Polyhedron')))
        A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        sage: type(_)
        <class 'sage.geometry.polyhedron.backend_ppl.Polyhedra_ZZ_ppl_with_category.element_class'>

        sage: next(iter(PALPreader(2, output='PPL')))
        A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
        sage: type(_)
        <class 'sage.geometry.polyhedron.ppl_lattice_polygon.LatticePolygon_PPL_class'>

        sage: next(iter(PALPreader(2, output='PointCollection')))
        [ 1,  0],
        [ 0,  1],
        [-1, -1]
        in Ambient free module of rank 2 over the principal ideal domain Integer Ring
        sage: type(_)
        <type 'sage.geometry.point_collection.PointCollection'>
    """

    def __init__(self, dim, data_basename=None, output='Polyhedron'):
        """
        The Python constructor

        See :class:`PALPreader` for documentation.

        TESTS::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
        """
        self._dim = dim
        if data_basename is not None:
            self._data_basename = data_basename
        else:
            import os
            from sage.env import SAGE_SHARE
            self._data_basename = os.path.join(SAGE_SHARE, 'reflexive_polytopes',
                                               'Full'+str(dim)+'d', 'zzdb')
            info = self._data_basename + '.info'
            if not os.path.exists(info):
                raise ValueError('Cannot find PALP database: '+info)
        from sage.geometry.polyhedron.parent import Polyhedra
        self._polyhedron_parent = Polyhedra(ZZ, dim)
        self._output = output.lower()

    def _palp_Popen(self):
        """
        Open PALP.

        OUTPUT:

        A PALP subprocess.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: polygons._palp_Popen()
            <subprocess.Popen object at 0x...>
        """
        return Popen(["class.x", "-b2a", "-di", self._data_basename], stdout=PIPE)

    def _read_vertices(self, stdout, rows, cols):
        """
        Read vertex data from the PALP output pipe.

        OUTPUT:

        A list of lists.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: palp = polygons._palp_Popen()
            sage: palp.stdout.readline()
            '2 3  \n'
            sage: polygons._read_vertices(palp.stdout, 2, 3)
            [[1, 0], [0, 1], [-1, -1]]
        """
        m = [ [] for col in range(0,cols) ]
        for row in range(0,rows):
            for col,x in enumerate(stdout.readline().split()):
                m[col].append(ZZ(x))
        return m

    def _read_vertices_transposed(self, stdout, rows, cols):
        """
        Read vertex data from the PALP output pipe.

        OUTPUT:

        A list of lists.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: palp = polygons._palp_Popen()
            sage: palp.stdout.readline()
            '2 3  \n'
            sage: polygons._read_vertices_transposed(palp.stdout, 2, 3)
            [[1, 0, -1], [0, 1, -1]]
        """
        m = []
        for row in range(0,rows):
            m.append( [ ZZ(x) for x in stdout.readline().split() ] )
        return m

    def _iterate_list(self, start, stop, step):
        """
        Iterate over the reflexive polytopes.

        INPUT:

        - ``start``, ``stop``, ``step`` -- integers specifying the
          range to iterate over.

        OUTPUT:

        A generator for vertex data as a list of lists.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: iter = polygons._iterate_list(0,4,2)
            sage: next(iter)
            [[1, 0], [0, 1], [-1, -1]]
        """
        if start is None:
            start = 0
        if step is None:
            step = 1
        palp = self._palp_Popen()
        try:
            palp_out = palp.stdout
            i = 0
            while True:
                l = palp_out.readline().strip()
                if l=='' or l.startswith('#'):
                    return  # EOF
                l=l.split()
                dim = ZZ(l[0]);  # dimension
                n = ZZ(l[1]);    # number of vertices
                if i>=start and (i-start) % step == 0:
                    if dim == self._dim:
                        vertices = self._read_vertices(palp_out, dim, n)
                    elif n == self._dim:  # PALP sometimes returns transposed data #@!#@
                        vertices = self._read_vertices_transposed(palp_out, dim, n)
                    else:
                        raise ValueError('PALP output dimension mismatch.')
                    yield vertices
                else:
                    for row in range(0,dim):
                        palp_out.readline()
                i += 1
                if stop is not None and i>=stop:
                    return
        finally:
            palp.poll()
            if palp.returncode is None:
                palp.terminate()
            palp.poll()
            if palp.returncode is None:
                palp.kill()


    def _iterate_Polyhedron(self, start, stop, step):
        """
        Iterate over the reflexive polytopes.

        INPUT:

        - ``start``, ``stop``, ``step`` -- integers specifying the
          range to iterate over.

        OUTPUT:

        A generator for lattice polyhedra.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: iter = polygons._iterate_Polyhedron(0,4,2)
            sage: next(iter)
            A 2-dimensional polyhedron in ZZ^2 defined as the convex hull of 3 vertices
        """
        parent = self._polyhedron_parent
        for vertices in self._iterate_list(start, stop, step):
            p = parent.element_class(parent, [vertices,[],[]], None)
            yield p

    def _iterate_PPL(self, start, stop, step):
        """
        Iterate over the reflexive polytopes.

        INPUT:

        - ``start``, ``stop``, ``step`` -- integers specifying the
          range to iterate over.

        OUTPUT:

        A generator for PPL-based lattice polyhedra.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: iter = polygons._iterate_PPL(0,4,2)
            sage: next(iter)
            A 2-dimensional lattice polytope in ZZ^2 with 3 vertices
        """
        for vertices in self._iterate_list(start, stop, step):
            yield LatticePolytope_PPL(*vertices)

    def _iterate_PointCollection(self, start, stop, step):
        """
        Iterate over the reflexive polytopes.

        INPUT:

        - ``start``, ``stop``, ``step`` -- integers specifying the
          range to iterate over.

        OUTPUT:

        A generator for PPL-based lattice polyhedra.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: iter = polygons._iterate_PointCollection(0,4,2)
            sage: next(iter)
            [ 1,  0],
            [ 0,  1],
            [-1, -1]
            in Ambient free module of rank 2 over the principal ideal domain Integer Ring
        """
        from sage.modules.free_module import FreeModule
        N = FreeModule(ZZ, self._dim)
        from sage.geometry.point_collection import PointCollection
        for vertices in self._iterate_list(start, stop, step):
            yield PointCollection(vertices, module=N)

    def _iterate(self, output=None):
        """
        Iterate over the reflexive polytopes.

        INPUT:

        - ``output`` -- as in the :class:`PALPreader` constructor.

        OUTPUT:

        A function generating lattice polytopes in the specified output format.

        EAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: func = polygons._iterate(output='list')
            sage: func
            <bound method PALPreader._iterate_list of <class 'sage.geometry.polyhedron.palp_database.PALPreader'>>
            sage: iter = func(0,1,1)
            sage: next(iter)
            [[1, 0], [0, 1], [-1, -1]]
        """
        if output is None:
            output = self._output
        if output == 'polyhedron':
            return self._iterate_Polyhedron
        elif output == 'ppl':
            return self._iterate_PPL
        elif output == 'pointcollection':
            return self._iterate_PointCollection
        elif output == 'list':
            return self._iterate_list
        else:
            raise TypeError('Unknown output format (='+str(self._output)+').')

    def __iter__(self):
        """
        Iterate over all polytopes.

        OUTPUT:

        An iterator for all polytopes.

        TESTS::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: polygons = PALPreader(2)
            sage: polygons.__iter__()
            <generator object _iterate_Polyhedron at 0x...>
        """
        iterator = self._iterate()
        return iterator(None, None, None)

    def __getitem__(self, item):
        """
        Return the polytopes(s) indexed by ``item``.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import PALPreader
            sage: palp = PALPreader(3)
            sage: list(palp[10:30:10])
            [A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices,
             A 3-dimensional polyhedron in ZZ^3 defined as the convex hull of 4 vertices]
        """
        iterator = self._iterate()
        if isinstance(item, slice):
            return iterator(item.start, item.stop, item.step)
        else:
            try:
                return next(iterator(item, item+1, 1))
            except StopIteration:
                raise IndexError('Index out of range.')



#########################################################################
class Reflexive4dHodge(PALPreader):
    """
    Read the PALP database for Hodge numbers of 4d polytopes.

    The database is very large and not installed by default. You can
    install it with the shell command ``sage -i polytopes_db_4d``.

    INPUT:

    - ``h11``, ``h21`` -- Integers. The Hodge numbers of the reflexive
      polytopes to list.

    Any additional keyword arguments are passed to
    :class:`PALPreader`.

    EXAMPLES::

        sage: from sage.geometry.polyhedron.palp_database import Reflexive4dHodge
        sage: ref = Reflexive4dHodge(1,101)             # optional - polytopes_db_4d
        sage: next(iter(ref)).Vrepresentation()         # optional - polytopes_db_4d
        (A vertex at (-1, -1, -1, -1), A vertex at (0, 0, 0, 1),
         A vertex at (0, 0, 1, 0), A vertex at (0, 1, 0, 0), A vertex at (1, 0, 0, 0))
    """
    def __init__(self, h11, h21, data_basename=None, **kwds):
        """
        The Python constructor.

        See :class:`Reflexive4dHodge` for documentation.

        TESTS::

        sage: from sage.geometry.polyhedron.palp_database import Reflexive4dHodge
        sage: Reflexive4dHodge(1,101)  # optional - polytopes_db_4d
        <class 'sage.geometry.polyhedron.palp_database.Reflexive4dHodge'>
        """
        dim = 4
        if data_basename is None:
            import os
            from sage.env import SAGE_SHARE
            data_basename = os.path.join(SAGE_SHARE, 'reflexive_polytopes',
                                         'Hodge4d', 'all')
            info = data_basename + '.vinfo'
            if not os.path.exists(info):
                raise ValueError('Cannot find PALP database: '+info+
                                 '. Did you install the polytopes_db_4d optional spkg?')
        PALPreader.__init__(self, dim, data_basename=data_basename, **kwds)
        self._h11 = h11
        self._h21 = h21

    def _palp_Popen(self):
        """
        Open PALP.

        OUTPUT:

        A PALP subprocess.

        EXAMPLES::

            sage: from sage.geometry.polyhedron.palp_database import Reflexive4dHodge
            sage: polygons = Reflexive4dHodge(1, 101)   # optional - polytopes_db_4d
            sage: polygons._palp_Popen()                # optional - polytopes_db_4d
            <subprocess.Popen object at 0x...>
        """
        return Popen(['class-4d.x', '-He',
                      'H'+str(self._h21)+':'+str(self._h11)+'L100000000',
                      '-di', self._data_basename], stdout=PIPE)



