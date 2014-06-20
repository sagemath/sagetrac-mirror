r"""
Lattice Polytopes in the GRDB, Library

Here is the documentation on all the functions mentioned
in :mod:`lattice_polytopes<sage.geometry.lattice_polytopes>`.
For a complete explanation of this library of functions
see the :mod:`lattice_polytopes<sage.geometry.lattice_polytopes>`
documentation.

These functions provide access to databases that form part of
the Graded Ring Database.  If you use these databases as part
of published research, please cite the appropriate papers as
listed at http://grdb.lboro.ac.uk/
"""
from sage.rings.all import Integer
from sage.misc.misc import SAGE_SHARE
from sage.geometry.lattice_polytope import LatticePolytope, all_cached_data
import os

GRDB_location = os.path.join(SAGE_SHARE,'grdb_polytopes')

def _GRDBBlockReader(n,databasename,blocksize,compute_properties=False):
    r"""
    Returns a list of lattice polytopes of the ``k``-th polytopes from
    for ``k`` in ``n`` from the database ``databasename`` which has blocks
    of size ``blocksize``.

    INPUT:

    - ``n`` -- a list of integers, each greater than or equal to one, corresponding
      to the Polytope IDs in the database
    - ``databasename`` -- string representing the name of the database
      to read in grdb_polytopes data
    - ``blocksize`` -- an integer, of the size of a block in ``databasename``
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A list of lattice polytopes.

    .. SEEALSO::

        :class:`sage.geometry.lattice_polytope.LatticePolytope`

    EXAMPLES::
        sage: import sage.geometry.lattice_polytopes_backend
        sage: M=sage.geometry.lattice_polytopes_backend._GRDBBlockReader([4654,31],"canonical3",5000); M # optional - grdb_polytopes
        [A lattice polytope: 3-dimensional, 12 vertices., A lattice polytope: 3-dimensional, 14 vertices.]

        sage: M=sage.geometry.lattice_polytopes_backend._GRDBBlockReader([6],"canonical2",16,True); M # optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 4 vertices.]

    """
    #prepare n to sort it
    n=zip(range(len(n)),n)
    n.sort(key=lambda s: s[1])


    polytopes = []
    block_and_line = [((k[1] - 1) // blocksize,((k[1] - 1) % blocksize) + 1) for k in n]
    current_block=block_and_line[0][0]
    current_line=0
    #Open the appropriate block file. The first line tells us the base the
    #data is encoded in. Extract that, then fetch the required line.
    filename = os.path.join(GRDB_location,databasename,"block" + str(current_block))
    try:
        fh = open(filename,"r")
    except(IOError):
        raise NotImplementedError("Cannot find database! Is the grdb_polytopes package installed?")
    base = fh.readline()
    #Open the appropriate block file. The first line tells us the base the
    #data is encoded in. Extract that, then fetch the required line.
    for block,line_number in block_and_line:
        if not block == current_block:
            current_block=block
            current_line=0
            filename = os.path.join(GRDB_location,databasename,"block" + str(block))
            try:
                fh = open(filename,"r")
            except(IOError):
                raise NotImplementedError("Cannot find database! Is the grdb_polytopes package installed?")
            base = fh.readline()
        while not current_line == line_number:
            line = fh.readline()
            current_line += 1

        #Convert the base and line into integers
        base_number = int(base)
        line = Integer(line)

        #Now unpack the integer into a sequence
        coeffs = line.digits(base=base_number)

        #The first entry in the sequence is the dimension of the polytope, the
        #second entry is a shift value we'll need
        dim = coeffs[0]
        shift = coeffs[1]

        #We need to subtract "shift" from the remaining values.
        coeffs = [ coeffs[i] - shift for i in range(2, len(coeffs))]

        #Finally collect the coefficients together into groups of the correct
        #dimension
        vertices = [ [ coeffs[dim * i + j] for j in range(dim) ] for
                    i in range(len(coeffs) // dim )]

        #We can create the polytopes without computing vertices: the database, guarantees full
        #dimension and no interior points
        polytopes.append(LatticePolytope(vertices,compute_vertices=False,copy_vertices=False))

    if compute_properties:
        all_cached_data(polytopes)

    #re-sort data
    polytopes=[x for (y,x) in sorted(zip([i[0] for i in n],polytopes))]

    return polytopes

def _GRDBlistpacking(n,databasename,blocksize,compute_properties,maximum,offset=0):
    r"""
    Acts as a wrapper around :func:`_GRDBBlockReader`, adding
    out of range error handling and converting ``n`` from
    an integer to a list and vice versa on output.

    INPUT:

    - ``n`` -- An integer or list of integers, each greater than or
      equal to one, corresponding to the Polytope IDs in the database.
    - ``databasename`` -- A string representing the name of the database
      to read in grdb_polytopes data.
    - ``blocksize`` -- an integer, of the size of a block in ``databasename``
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.
    - ``maximum`` -- An integer containing the total number of polytopes in this
      database.
    - ``offset`` -- An integer indicating a shift from which line to begin
      searching the database from.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :class:`sage.geometry.lattice_polytope.LatticePolytope`
        :func:`_GRDBBlockReader`

    EXAMPLES::
        sage: import sage.geometry.lattice_polytopes_backend
        sage: M=sage.geometry.lattice_polytopes_backend._GRDBlistpacking([4654,31],"canonical3",5000,False,674688); M # optional - grdb_polytopes
        [A lattice polytope: 3-dimensional, 12 vertices., A lattice polytope: 3-dimensional, 14 vertices.]

        sage: M=sage.geometry.lattice_polytopes_backend._GRDBlistpacking(6,"canonical2",16,True,16); M # optional - grdb_polytopes
        A lattice polytope: 2-dimensional, 4 vertices.


    """
    input_was_list = (type(n) == list)
    if not input_was_list:
        #pack into a list if not a list
        n=[n]

    for l,k in enumerate(n):
        if k < 1 or k > maximum:
            raise ValueError("there is an invalid ID, all IDs must be in range 1 to "+str(maximum))
        n[l]+=offset


    out=_GRDBBlockReader(n,databasename,blocksize,compute_properties)

    if not input_was_list:
        #unpack list if input wasn't a list
        out=out[0]

    return out

def CanonicalFano(dim,n,compute_properties=False):
    r"""
    Returns ``n``-th canonical Fano polytope of dimension ``dim`` from the database
    of 2 and 3 dimensional canonical Fano polytopes.

    A Fano polytope is a lattic polytope of maximal dimension
    containing the origin in its strict interior and whose vertices
    are primitive.

    A canonical Fano polytope is a Fano polytope with the only
    strict interior lattice point being the origin.

    INPUT:

    - ``dim`` -- An integer specifying the dimension of the polytope.
    - ``n`` -- An integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`


    EXAMPLES::

        sage: C = lattice_polytopes.CanonicalFano(3,5764); C # optional - grdb_polytopes
        A lattice polytope: 3-dimensional, 12 vertices.

        sage: C =lattice_polytopes.CanonicalFano(2,[10,14],True); C #optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 4 vertices., A lattice polytope: 2-dimensional, 3 vertices.]

    Any list of integers can be requested if all the entries are valid ID's::

        sage: C = lattice_polytopes.CanonicalFano(3,range(1,1000)) # optional - grdb_polytopes

    Watch out for possible errors::

        sage: C = lattice_polytopes.CanonicalFano(3,42342344) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 674688


    NOTES:

        #. Numeration starts with one, where for
           ``dim`` = 2 : `1 \leq k \leq 5` and for ``dim`` = 3
           : `1 \leq k \leq 674688` for each `k` in ``n``
           These IDs correspond with the IDs in the Graded
           Ring Database.


    REFERENCES:

    .. [Kas10]  Alexander M. Kasprzyk, "Canonical toric Fano threefolds",
        Canadian Journal of Mathematics, 62 (2010), no. 6, 1293-1309.
    .. [GRDBfc] http://grdb.lboro.ac.uk/forms/toricf3c
    """
    if dim==2:
        out=_GRDBlistpacking(n,"canonical2",16,compute_properties,16)
    elif dim==3:
        out=_GRDBlistpacking(n,"canonical3",5000,compute_properties,674688)
    else:
        raise NotImplementedError("only dimensions 2 and 3 are implemented.")
    return out

def TerminalFano(dim,n,compute_properties=False):
    r"""
    Returns ``n``-th terminal Fano polytope of
    dimension ``dim`` from the database of 2 and 3 dimensional
    terminal Fano polytopes.

    A Fano polytope is a lattic polytope of maximal dimension
    containing the origin in its strict interior and whose vertices
    are primitive.

    A terminal Fano polytope is a Fano polytope with the
    only boundary lattice points being vertices.

    INPUT:

    - ``dim`` -- An integer specifying the dimension of the polytope.
    - ``n`` -- An integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`

    EXAMPLES::

        sage: C = lattice_polytopes.TerminalFano(3,312); C # optional - grdb_polytopes
        A lattice polytope: 3-dimensional, 8 vertices.

        sage: C =lattice_polytopes.TerminalFano(2,[3,4],True); C #optional - grdb_polytopes
            [A lattice polytope: 2-dimensional, 4 vertices., A lattice polytope: 2-dimensional, 4 vertices.]

    Any list of integers can be requested if all the entries are valid ID's::

        sage: C = lattice_polytopes.TerminalFano(3,range(1,300)) # optional - grdb_polytopes

    Watch out for possible errors::

        sage: C = lattice_polytopes.TerminalFano(2,17) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 5


    NOTES:

        #. Numeration starts with one, where for
           ``dim`` = 2 : `1 \leq k \leq 16` and for ``dim`` = 3
           : `1 \leq k \leq 634` for each `k` in ``n``
           These IDs correspond with the IDs in the Graded
           Ring Database.
        #. For the terminal polygons we use the result that
           these are exactly the smooth Fano polygons.


    REFERENCES:

    .. [Kas06]  Alexander M. Kasprzyk, "Toric Fano threefolds with terminal
        singularities", Tohoku Mathematical Journal, 58 (2006),
        no. 1, 101-121.
    .. [GRDBtf] http://grdb.lboro.ac.uk/forms/toricf3t

    """

    if dim==2:
            out=SmoothFano(n,2,compute_properties)
    elif dim==3:
        out=_GRDBlistpacking(n,"terminal3",50,compute_properties,634)
    else:
        raise NotImplementedError("only dimensions 2 and 3 are implemented.")
    return out

def SmoothFano(n,dim=None,compute_properties=False):
    r"""
    Returns ``n``-th smooth Fano polytope from the database of 2- to
    8-dimensional smooth Fano polytopes if no dimension ``dim`` is specified;
    otherwise returns the ``n``-th smooth Fano polytope of dimension ``dim``

    A Fano polytope is a lattic polytope of maximal dimension
    containing the origin in its strict interior and whose vertices
    are primitive.

    A smooth Fano polytope is a Fano polytope such that for
    each facet the vertices of that facet form a basis for the lattice.


    INPUT:

    - ``n`` -- An integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``dim`` -- An integer specifying the dimension of the polytope, or
      ``None`` (by default) if the indexing the whole database of 2- to
      8-dimensional polytopes.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.



    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`

    EXAMPLES:

    ::

        sage: S = lattice_polytopes.SmoothFano(27958); S # optional - grdb_polytopes
        A lattice polytope: 7-dimensional, 13 vertices.


        sage: S = lattice_polytopes.SmoothFano([4,545,32],compute_properties=True); S # optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 4 vertices., A lattice polytope: 5-dimensional,
            10 vertices., A lattice polytope: 4-dimensional, 8 vertices.]

        sage: S = lattice_polytopes.SmoothFano(534,7); S # optional - grdb_polytopes
        A lattice polytope: 7-dimensional, 13 vertices.


    There are a few possible errors that can be thrown:

    ::


        sage: S = lattice_polytopes.SmoothFano(42342344) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 830783

        sage: S = lattice_polytopes.SmoothFano(21,2) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 5


    NOTES:
        #. Numeration starts with one, where for
           ``dim`` = 2 : `1 \leq k \leq 16` and for ``dim`` = 3
           : `1 \leq k \leq 634` for each `k` in ``n``
           These IDs correspond with the IDs in the Graded
           Ring Database.
        #. Numeration starts with one where, if ``dim`` is
           ``None``, `1 \leq k \leq 830783` for each `k` in ``n``
           For a specified dimension, for each `k` in ``n``:

            - ``dim`` = 2 then `1 \leq k \leq 5`
            - ``dim`` = 3 then `1 \leq k \leq 18`
            - ``dim`` = 4 then `1 \leq k \leq 124`
            - ``dim`` = 5 then `1 \leq k \leq 866`
            - ``dim`` = 6 then `1 \leq k \leq 7622`
            - ``dim`` = 7 then `1 \leq k \leq 72256`
            - ``dim`` = 8 then `1 \leq k \leq 749891`

           These IDs correspond with the IDs in the Graded
           Ring Database.
        #. The PALP backend cannot compute properties
           for smooth Fano polytopes of dimension 6 or greater,
           if ``compute_properties`` is ``True`` it is likely to throw an error.


    REFERENCES:

    .. [Obr]  Mikkel \Obro, "An algorithm for the classification of smooth Fano
        polytopes", arXiv:0704.0049v1.
    .. [GRDBts] http://grdb.lboro.ac.uk/forms/toricsmooth


    """
    if dim:
        if not 2<=dim<=8:
            raise NotImplementedError("only dimensions 2 to 8 can be specified!")
        dimension_props={2:(0,"smoothfano",250),3:(5,"smoothfano",250),4:(23,"smoothfano",250),
                         5:(147,"smoothfano",250),6:(1013,"smoothfano",250),7:(8635,"smoothfano7",722),
                         8:(80891,"smoothfano8",7498),9:(830782,)}
        shift=dimension_props[dim][0]
        maximum=dimension_props[dim+1][0]-dimension_props[dim][0]
        return _GRDBlistpacking(n,dimension_props[dim][1],dimension_props[dim][2],compute_properties,maximum,shift)
    else:
    #must split list to read from appropriate databases

        input_was_list = (type(n) == list)
        if not input_was_list:
            #pack into a list if not a list
            n=[n]

        out=[]
        list_2d_to6d=[]
        list_7d=[]
        list_8d=[]
        for k in n:
            if k < 1:
                raise ValueError("there is an invalid ID, enumeration begins from 1!")
            elif k < 8636:
                list_2d_to6d.append(k)
            elif k < 80892:
                list_7d.append(k-8635)
            elif k < 830784:
                list_8d.append(k-80891)
            elif k > 830783:
                raise ValueError("there is an invalid ID, all IDs must be in range 1 to 830783")

    if list_2d_to6d: out.extend(_GRDBBlockReader(list_2d_to6d,"smoothfano",250,compute_properties))
    if list_7d: out.extend(_GRDBBlockReader(list_7d,"smoothfano7",722,compute_properties))
    if list_8d: out.extend(_GRDBBlockReader(list_8d,"smoothfano8",7498,compute_properties))

    if not input_was_list:
        #unpack list if input wasn't a list
        out=out[0]

    return out

def SmallPolygon(n,compute_properties=False):
    r"""
    Returns the ``n``-th small polygon, from the database
    of all lattice polygons, up to translation and change
    of basis, that can lie in a `[0,7] \times [0,7]` box

    INPUT:

    - ``n`` -- an integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`

    EXAMPLES:

    ::

        sage: S = lattice_polytopes.SmallPolygon(574); S # optional - grdb_polytopes
        A lattice polytope: 2-dimensional, 5 vertices.

        sage: S = lattice_polytopes.SmallPolygon([23,133,44,65],True); S # optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 6 vertices., A lattice polytope: 2-dimensional, 4 vertices.,
            A lattice polytope: 2-dimensional, 5 vertices., A lattice polytope: 2-dimensional, 5 vertices.]

    There are a few possible errors that can be thrown:

    ::


        sage: S = lattice_polytopes.SmallPolygon(42342344) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 1249439


    NOTES:

        #. Numeration starts with one, with
           `1 \leq k \leq 1249439` for each `k` in ``n``



    REFERENCES:

    .. [BK] Gavin Brown and Alexander M. Kasprzyk, "Small polygons
        and toric codes", Journal of Symbolic Computation (2012),
        doi:10.1016/j.jsc.2012.07.001.
    .. [GRDBbp] http://grdb.lboro.ac.uk/forms/boxpoly


    """

    return _GRDBlistpacking(n,"smallpolygons",4861,compute_properties,1249439)

def lReflexive(n,compute_properties=False):
    r"""
    Returns the ``n``-th `l`-reflexive polygon from the
    database of `l`-reflexive polygons for `1\leq l \leq 200`.

    A lattice polygon is `l`-reflexive if and only if
    its polar, scaled by `l`, is reflexive.

    INPUT:

    - ``n`` -- an integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`
        :func:`lattice_polytopes.ReflexiveFano<ReflexiveFano>`

    EXAMPLES:

    ::

        sage: L = lattice_polytopes.lReflexive(574); L # optional - grdb_polytopes
        A lattice polytope: 2-dimensional, 3 vertices.

        sage: L = lattice_polytopes.lReflexive([3,654,574],True); L # optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 5 vertices., A lattice polytope: 2-dimensional, 4 vertices.,
            A lattice polytope: 2-dimensional, 3 vertices.]

    There are a few possible errors that can be thrown:

    ::

        sage: L = lattice_polytopes.lReflexive(42342344) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 41458


    NOTES:

        #. Numeration starts with one, with
           `1 \leq k \leq 41458` for each `k` in ``n``



    REFERENCES:

    .. [KN]   Alexander M. Kasprzyk and Benjamin Nill, "Reflexive polytopes
        of higher index and the number 12", Electronic Journal of Combinatorics,
        19 (2012), no. 3, P9.
    .. [GRDBlr] http://grdb.lboro.ac.uk/forms/toriclr2

    """

    return _GRDBlistpacking(n,"lreflexive2",500,compute_properties,41458)

def ReflexiveFano(dim,n,compute_properties=False):
    r"""
    Returns ``n``-th reflexive Fano polytope of dimension ``dim`` from the database
    of 2 and 3 dimensional reflexive Fano polytopes.

    A Fano polytope is a lattic polytope of maximal dimension
    containing the origin in its strict interior and whose vertices
    are primitive.

    A reflexive Fano lattice polytope is a Fano polytope such that its
    polar is also a Fano polytope.

    INPUT:

    - ``dim`` -- An integer specifying the dimension of the polytope.
    - ``n`` -- An integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes.

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`
        :func:`lReflexive`

    EXAMPLES::

        sage: C = lattice_polytopes.ReflexiveFano(3,432); C # optional - grdb_polytopes
        A lattice polytope: 3-dimensional, 4 vertices.

        sage: C =lattice_polytopes.ReflexiveFano(2,[3,4],True); C #optional - grdb_polytopes
            [A lattice polytope: 2-dimensional, 5 vertices., A lattice polytope: 2-dimensional, 5 vertices.]

    Any list of integers can be requested if all the entries are valid ID's::

        sage: C = lattice_polytopes.ReflexiveFano(3,range(1,300)) # optional - grdb_polytopes

    Watch out for possible errors::

        sage: C = lattice_polytopes.ReflexiveFano(2,17) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 16


    NOTES:

        #. Numeration starts with one, where for
           ``dim`` = 2 : `1 \leq k \leq 16` and for ``dim`` = 3
           : `1 \leq k \leq 4319` for each `k` in ``n``
           These IDs correspond with the IDs in the Graded
           Ring Database.
        #. For the reflexive Fano polygons
           we use the result that these are exactly the
           canonical Fano polygons.
        #. The Fano reflexive polytopes are the reflexive
           polytopes returned by :func:`sage.geometry.lattice_polytope.ReflexivePolytope`
           up to change of lattice basis and difference in
           Database IDs
    """
    if dim==2:
            out=CanonicalFano(2,n,compute_properties)
    elif dim==3:
        out=_GRDBlistpacking(n,"reflexive3",4319,compute_properties,4319)
    else:
        raise NotImplementedError("only dimensions 2 and 3 are implemented.")
    return out

def LDP(n,compute_properties=False):
    r"""
    Returns the ``n``-th LDP polygon from the
    database of LDP polygons.

    A lattice polygon `P` is LDP if and only if
    the vertices are primitive and the origin is contained
    in the strict interior of `P`.

    These correspond to log del Pezzo surfaces, in that
    the toric variety generated by the fan of `P` is log
    del Pezzo if and only if `P` is LDP.

    INPUT:

    - ``n`` -- an integer or list of integers all greater than or equal to one,
      corresponding to the Polytope IDs in the Graded Ring Database.
    - ``compute_properties`` -- A Boolean, that if true calculates a range
      of properties by use of the all* functions in
      :mod:`sage.geometry.lattice_polytope` for each polytope, and ensures this data
      is saved in the lattice polytopes that are returned.


    OUTPUT:
    A lattice polytope or a list of lattice polytopes

    .. SEEALSO::

        :func:`sage.geometry.lattice_polytope.LatticePolytope`

    EXAMPLES:

    ::

        sage: L = lattice_polytopes.LDP(574); L # optional - grdb_polytopes
        A lattice polytope: 2-dimensional, 4 vertices.

        sage: L = lattice_polytopes.LDP([3,654,574],True); L # optional - grdb_polytopes
        [A lattice polytope: 2-dimensional, 5 vertices., A lattice polytope: 2-dimensional,
            4 vertices., A lattice polytope: 2-dimensional, 4 vertices.]


    There are a few possible errors that can be thrown:

    ::

        sage: L = lattice_polytopes.LDP(42342344) # optional - grdb_polytopes
        Traceback (most recent call last):
        ...
        ValueError: there is an invalid ID, all IDs must be in range 1 to 15346


    NOTES:

        #. Numeration starts with one, with
           `1 \leq k \leq 15346` for each `k` in ``n``



    REFERENCES:

    .. [KKN]  Alexander M. Kasprzyk, Maximilian Kreuzer, Benjamin Nill, "On the
        combinatorial classification of toric log del Pezzo surfaces",
        MS Journal of Computation and Mathematics, 13 (2010), 33-46.
    .. [GRDBldp] http://grdb.lboro.ac.uk/forms/toricldp

    """
    return _GRDBlistpacking(n,"ldp",250,compute_properties,15346)
