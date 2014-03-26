r"""
Labelled matrices

Matrices and (two dimensional) arrays are typically indexed by the integers. The
class :class:`LabelledMatrix` allows two dimensional arrays to be indexed
by arbitrary sage objects. The purpose of this class is to do the following::

#.  display ``labelled matrices`` nicely, together their row and column
indices,

#.  allow a more intuitive way of accessing and setting the entries of such
matrices uses their natural indexing sets (rather than just using the integers).

#. allow the row and column spaces of a matrix to be naturally identified
with a :class:`CombinatorialFreeModule`.

These *labelled* arrays and matrices should be useful for decomposition
matrices, character tables, Cayley tables, homomorphisms of combinatorial free
modules, ....

More precisely, an instance of the class :class:`LabelledMatrix` is an
array with specified row and column labels and which have arbitrary (not
necessarily numerical) entries.  The indexing sets for the row and column labels
do not have to be the same and the entries of the array can, in principle, be
arbitrary sage objects, although typically they would be elements of some ring.
The entries in the underlying array can also be ``None`` in which case they are
treated as unknown quantities which, potentially, sage or the user is in the
process of calculating.

    sage: from sage.matrix.labelled_matrix import LabelledMatrix
    sage: LabelledMatrix(Partitions(5)).display()
    [5]             | ? ? ? ? ? ? ?
    [4, 1]          | ? ? ? ? ? ? ?
    [3, 2]          | ? ? ? ? ? ? ?
    [3, 1, 1]       | ? ? ? ? ? ? ?
    [2, 2, 1]       | ? ? ? ? ? ? ?
    [2, 1, 1, 1]    | ? ? ? ? ? ? ?
    [1, 1, 1, 1, 1] | ? ? ? ? ? ? ?
    sage: LabelledMatrix(Partitions(5), row_label='compact_high', default_value=0).display()
    5     | . . . . . . .
    4,1   | . . . . . . .
    3,2   | . . . . . . .
    3,1^2 | . . . . . . .
    2^2,1 | . . . . . . .
    2,1^3 | . . . . . . .
    1^5   | . . . . . . .
    sage: d=LabelledMatrix(Partitions(5), row_label='compact_high', default_value=0, unitriangular=True);d.display()
    5     | 1
    4,1   | . 1
    3,2   | . . 1
    3,1^2 | . . . 1
    2^2,1 | . . . . 1
    2,1^3 | . . . . . 1
    1^5   | . . . . . . 1
    sage: d[[3,2],[5]]=1; d.display()
    5     | 1
    4,1   | . 1
    3,2   | 1 . 1
    3,1^2 | . . . 1
    2^2,1 | . . . . 1
    2,1^3 | . . . . . 1
    1^5   | . . . . . . 1
    sage: d.options(row_order=lambda a,b: 1 if a>b else -1, col_label='compact_low',triangular='upper')
    sage: d.display()
    <BLANKLINE>
    -------------------------------------------------
        |   1^5 1^3,2 1,2^2 1^2,3   2,3   1,4     5
    ------|------------------------------------------
    1^5   |     1     .     .     .     .     .     .
    2,1^3 |           1     .     .     .     .     .
    2^2,1 |                 1     .     .     .     .
    3,1^2 |                       1     .     .     .
    3,2   |                             1     .     1
    4,1   |                                   1     .
    5     |                                         1

The indexing set for a :class:`LabelledMatrix` can be any iterable sage
object or a :class:`CombinatorialFreeModule`, in which case the rows and columns
of the array/matrix/table can be interpreted as elements of this module.

If an instance of :class:`LabelledMatrix` is (known to be) a matrix then
any matrix methods will automatically be applied to the underlying matrix. If
the result of the method is another matrix of the same size then the
corresponding :class:`LabelledMatrix` will be returned.

The indexing sets for the rows and columns can be quite general::

    sage: d=LabelledMatrix(rows=Partitions(5), columns=['dim S(mu)','dim D(mu)'],default_value=0)
    sage: for mu in d.rows(): d[mu,0]=StandardTableaux(mu).cardinality()
    sage: for mu in d.rows(): d[mu,1]=GradedSpechtModule(3,mu).simple_dimension() # TODO: not implemented
    sage: d[:,1]=[0,0,4,6,1,4,1]     # a work around until graded Specht modules are merged
    sage: d.options(row_label='compact_high',col_label='repr', name='Dimensions in characteristic 3')
    sage: d.display()
    Dimensions in characteristic 3
    ---------------------------
          | dim S(mu) dim D(mu)
    ------|--------------------
    5     |         1         .
    4,1   |         4         .
    3,2   |         5         4
    3,1^2 |         6         6
    2^2,1 |         5         1
    2,1^3 |         4         4
    1^5   |         1         1

    sage: d.transpose().display()
    Dimensions in characteristic 3
    -----------------------------------------------------
              |     5   4,1   3,2 3,1^2 2^2,1 2,1^3   1^5
    ----------|------------------------------------------
    dim S(mu) |     1     4     5     6     5     4     1
    dim D(mu) |     .     .     4     6     1     4     1



Labelled matrices and matrices
-----------------------------------

Not every :class:`LabelledMatrix` will be a matrix, however, if a
particular :class:`LabelledMatrix` is a matrix then matrix methods can be
applied to the underlying matrix -- although, tab-completion will not work for
most of these methods. ::


Indexing of row and column entries
----------------------------------


Row and column labels
---------------------


Reordering rows and columns
---------------------------


Indexing rows or columns by a :class:`CombinatorialFreeModule`
--------------------------------------------------------------


Options for a :class:`LabelledMatrix`
------------------------------------------

A :class:`LabelledMatrix` has many options which can be set when the
matrix is created or changed while working with the matrix.::


Extending the :class:`LabelledMatrix` class
------------------------------------------------

The :class:`LabelledMatrix` is not really intended to be used from the
command line. Instead, it is envisaged that this class will be used as the
base class by developers for constructing combinorially labelled matrices
which would typically have a richer feature set than an unlabelled matrix.
As an example, a rudimentary :class:`CharacterTable` has been implemented.::

    sage: G=SymmetricGroup(5); G.character_table()
    [ 1 -1  1  1 -1 -1  1]
    [ 4 -2  0  1  1  0 -1]
    [ 5 -1  1 -1 -1  1  0]
    [ 6  0 -2  0  0  0  1]
    [ 5  1  1 -1  1 -1  0]
    [ 4  2  0  1 -1  0 -1]
    [ 1  1  1  1  1  1  1]
    sage: from sage.matrix.labelled_matrix import CharacterTable
    sage: ct=CharacterTable(G); ct.display()
    Character table of Symmetric group of order 5! as a permutation group
    -------------------------
       | 1a 2a 2b 3a 6a 4a 5a
    ---|---------------------
    X0 |  1 -1  1  1 -1 -1  1
    X1 |  4 -2  0  1  1  0 -1
    X2 |  5 -1  1 -1 -1  1  0
    X3 |  6  0 -2  0  0  0  1
    X4 |  5  1  1 -1  1 -1  0
    X5 |  4  2  0  1 -1  0 -1
    X6 |  1  1  1  1  1  1  1
    sage: ct['X3']
    [6, 0, -2, 0, 0, 0, 1]
    sage: ct[3]
    [6, 0, -2, 0, 0, 0, 1]
    sage: ct[:,'3a']
    [1, 1, -1, 0, -1, 1, 1]

Apart from the code it inherits from :class:`LabelledMatrix` and a call to
:meth:`character_table`, the implementation of :class:`CharacterTable` is almost
trivial. The only extra lines of `code` that are needed are a few lines to work
out the labelling of the columns for the conjugacy classes. In order to make
this into a useful class some character theoretic methods would need to be added
to the class, such as methods to decompose characters into a linear combination
of irreducible characters.  It would also be good to make the row space of a
character table equal to a better implementation of the character ring of the
corresponding group as a :class:`CombinatorialFreeModule`.  Finally, as in Gap,
a better way of displaying non-integral character values is needed, as the
following example indicates::

    sage: from sage.matrix.labelled_matrix import CharacterTable
    sage: CharacterTable(CyclicPermutationGroup(6)).display()
    Character table of Cyclic group of order 6 as a permutation group
    ----------------------------------------------------------------------
       |         1a         6a         3a         2a         3b         6b
    ---|------------------------------------------------------------------
    X0 |          1          1          1          1          1          1
    X1 |          1         -1          1         -1          1         -1
    X2 |          1 -zeta3 - 1      zeta3          1 -zeta3 - 1      zeta3
    X3 |          1  zeta3 + 1      zeta3         -1 -zeta3 - 1     -zeta3
    X4 |          1      zeta3 -zeta3 - 1          1      zeta3 -zeta3 - 1
    X5 |          1     -zeta3 -zeta3 - 1         -1      zeta3  zeta3 + 1


TODO:

- add mechanism for saving and restoring combinatorial matrices - and possibly
  linking with a database


AUTHORS:

- Andrew Mathas (2012-12-17): initial version

"""

#*****************************************************************************
#  Copyright (C) 2012 Andrew Mathas andrew dot mathas at sydney dot edu dot au
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import types

from __builtin__ import classmethod
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.categories.modules_with_basis import ModulesWithBasis
from sage.categories.rings import Rings
from sage.categories.sets_cat import Sets
from sage.combinat.combinat import CombinatorialObject
from sage.graphs.graph import Graph
from sage.interfaces.html_display import HTMLDisplay, HTMLTable
from sage.matrix.constructor import matrix
from sage.misc.classcall_metaclass import ClasscallMetaclass
from sage.misc.flatten import flatten
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.rings.semirings.non_negative_integer_semiring import NN
from sage.sets.family import LazyFamily
from sage.structure.element import is_Element, is_Matrix
from sage.structure.parent import Parent
from sage.structure.sage_object import SageObject



class LabelledMatrix(SageObject):
    r"""
    A :class:`LabelledMatrix` provides a mechanism for storing data in an
    array, or matrix, which is indexed by sage objects, which can be tuples,
    lists, sets and so on but will typically be a :class:`CombinatorialObject`
    or a (basis for a) :class:`CombinatorialFreeModule`.

    This class is intended to be used as a base class by other applications.

    INPUTS:

    - rows        keys to index the rows.                     (required)
    - row_label  name of the method to print row labels     (default: repr)
    - row_order   comparison function for ordering rows       (default: None)
    - columns     keys to index the columns                   (default: rows)
    - col_label  name of the method to print column labels  (default: None)
    - col_order   comparison function for ordering columns    (default: None)
    - array       the underlying array/matrix/table           (default: [])
    - is_matrix   ``True`` if the array is a matrix           (default: False)
    - triangular  False, ``upper`` or ``lower``               (default: False)
    - default     default value                               (default: None)
    - base_ring   default base ring for the entries           (default: ZZ)
    - name        a descriptive identifier for the array      (default: None)
    - unknown_string a string for displaying unknown entries  (default: ?)

    If ``triangular``  is ``upper`` or ``lower`` then ``columns`` should be a
    subset of ``rows``, or ``rows`` a subset of ``columns``, although this is
    not enforced. It is necessary, however, that comparisons such as `row<col`
    make sense for ``row`` in ``rows`` and ``col`` in ``columns``.

    EXAMPLES::

        sage: from sage.matrix.labelled_matrix import LabelledMatrix
        sage: P4=CombinatorialFreeModule(ZZ,Partitions(4),prefix='P')
        sage: K4=Partitions(4, min_slope=-2)  # 3-regular partitions
        sage: d=LabelledMatrix(rows=P4,columns=K4, unitriangular=True, default_value=0)
        sage: d.display()
        <BLANKLINE>
        ------------------------------------------------------------------
                     |       [3, 1]       [2, 2]    [2, 1, 1] [1, 1, 1, 1]
        -------------|----------------------------------------------------
        [4]          |                                                    
        [3, 1]       |            1                                       
        [2, 2]       |            .            1                          
        [2, 1, 1]    |            .            .            1             
        [1, 1, 1, 1] |            .            .            .            1
        sage: d.options(row_label='compact_high'); d.display()
        <BLANKLINE>
        -----------------------------------------------------------
              |       [3, 1]       [2, 2]    [2, 1, 1] [1, 1, 1, 1]
        ------|----------------------------------------------------
        4     |                                                    
        3,1   |            1                                       
        2^2   |            .            1                          
        2,1^2 |            .            .            1             
        1^4   |            .            .            .            1
        sage: d[[1,1,1,1]]
        P[[1, 1, 1, 1]]
        sage: d[[2,2],[1,1,1,1]]=2
        sage: d[[2,2]]
        2*P[[1, 1, 1, 1]] + P[[2, 2]]
        sage: d[[2,2]] in P4
        True

    By default, a unitriangular array is assumed to be lower triangular.
    """

    # is_matrix is used by __getattr__ so we set it to it's default now
    # so as to avoid any possible problems
    _is_matrix=False

    def __init_indices(self, index, index_set):
        r"""
        The rows and columns of a :class:`LabelledMatrix` can be indexed
        by many different ``sage`` objects and the way in which the elements
        of these index sets are accessed can vary. To get around this, we
        define several new methods for the matrix ``self`` which allow the
        elements of the indexing sets to be accessed in a uniform way.

        INPUT:

        - ``index``      either 'row' or 'col', corresponding to the index set
                         for the rows and columns

        - ``index_set``  an indexing *set* which can be any iteratable object,
                         but will typically be either a :class:`Parent` or a
                         :class:`CombinatorialFreeModule`

        This method defines the following methods and attributes for accessing
        the relevant indices. If ``index==row`` then these methods are:

        - ``_rows``       the indexing set/parent/...

        - ``_row_length`` the number of rows

        - ``_row_rank``   returns the rank of an index in ``_rows``. That is,
                          if ``t`` is in ``_rows`` then
                          ``_rows[ _row_rank(t) ] ==t``, at least in those
                          cases where ``_rows[t]`` makes sense.

        - ``_row_unrank`` returns the element of ``_rows`` of the specified
                          rank.

        If ``index_set`` is a  :class:`CombinatorialFreeModule` then the
        attribute and method are also defined:

        - ``_row_space``  The associated row space: ``self._row_space`` is
                          ``index_set``.

        - ``_row_sum``    this is rows.sum_of-terms() which is the function
                          which takes as it argument a list of 2-tuples
                          (ind,coeff), where ``ind`` is a row index and
                          ``coeff`` is its coefficient.  This is just a
                          convenient shorthand.

        If ``index==col`` then the analogues methods ``_cols``,
        ``_col_length,`` _col_rank``, ``_col_unrank`` and ``_col_sum`` are
        created. The definitions of these methods depends upon what type of
        object ``index_set`` is.

        EXAMPLES::

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: LabelledMatrix(rows=Family([1,2,3])).display()
            1 | ? ? ?
            2 | ? ? ?
            3 | ? ? ?
            sage: LabelledMatrix(rows=CombinatorialFreeModule(ZZ,[1,2,3])).display()
            1 | ? ? ?
            2 | ? ? ?
            3 | ? ? ?
            sage: LabelledMatrix(rows=[1,2,3],unitriangular=True,default_value=0).display()
            1 | 1
            2 | . 1
            3 | . . 1
            sage: LabelledMatrix(rows=(2*i for i in range(3)),unitriangular=True,triangular='upper',default_value=0).display()
            0 | 1 . .
            2 |   1 .
            4 |     1
        """
        if index_set in FiniteEnumeratedSets:
            setattr(self,'_'+index+'s', index_set)
            setattr(self,'_'+index+'_length',index_set.cardinality())
            setattr(self,'_'+index+'_rank',index_set.rank)
            setattr(self,'_'+index+'_unrank',index_set.unrank)

        elif index_set in ModulesWithBasis:
            if index_set.dimension() in NN:
                setattr(self,'_'+index+'s', index_set._basis_keys)
                setattr(self,'_'+index+'_length',len(index_set._basis_keys))
                setattr(self,'_'+index+'_rank',index_set._basis_keys.rank)
                setattr(self,'_'+index+'_unrank',index_set._basis_keys.unrank)
                setattr(self,'_'+index+'_space', index_set)
                setattr(self,'_'+index+'_sum',index_set.sum_of_terms)
            else:
                raise ValuError('module must be finite dimensional\n')

        elif index_set in Sets or isinstance(index_set,(list,tuple)):
            setattr(self,'_'+index+'s',list(index_set))
            setattr(self,'_'+index+'_length',len(index_set))
            setattr(self,'_'+index+'_rank',index_set.index)
            setattr(self,'_'+index+'_unrank',lambda r:index_set[r])

        elif isinstance(index_set, types.GeneratorType):
            real_index=[ind for ind in index_set]
            setattr(self,'_'+index+'s',real_index)
            setattr(self,'_'+index+'_length',len(real_index))
            setattr(self,'_'+index+'_rank',real_index.index)
            setattr(self,'_'+index+'_unrank',lambda r:real_index[r])

        else:
            raise ValueError('unrecognised %s index\n'%index)


    def __init__(self, rows, **options):
        r"""
        The main data stored in a :class:`LabelledMatrix` is the following::
        - _array  the actual array of data, which can be arbitrary
        - _rows   the indexing set for the rows
        - _cols   the indexing set for the columns


        """
        # Set the default options. We do this in __init__ so that different
        # instances of the class do not share the same defaults in memory.
        self._default_options = dict(
                         base_ring=ZZ,
                         col_label=None,
                         col_labels_html=None,
                         col_label_math_mode=True,
                         col_order=None,
                         default_value=None,
                         is_matrix=False,
                         name='',
                         row_label='repr',
                         row_labels_html=None,
                         row_label_math_mode=True,
                         row_order=None,
                         triangular=False,
                         unitriangular=False,
                         unknown_string='?',
                         zero_string='.',
                         entries_math_mode=True
        )

        # we start by setting the row labels and associated helper functions
        self.__init_indices('row', rows)

        # if the column labels aren't specified they default to the row labels
        if 'columns' in options:
            self.__init_indices('col', options['columns'])
            options.__delitem__('columns')
            if not 'col_label' in options:
                options['col_label']='repr'
        else:
            self.__init_indices('col', self._rows)

        # To process options we first set all of the default options and
        # then overwrite them with the specified options
        self.options(**self._default_options)
        if options!={}: self.options(**options)

        # Finally if the array was not given as an argument then set it to the
        # default. If triangular is not None then we set the corresponding off
        # diagonal elements equal to zero.
        if hasattr(self,'_array'):
            is_matrix=is_Matrix(self._array)
            #print 'array already set!: %s' % is_matrix
            if self._is_matrix and not is_matrix:
                raise ValueError('is_matrix==True, however, the supplied array is not a matrix\n')
            self._is_matrix=is_matrix
            if self._is_matrix:
                self._base_ring=self._array.base_ring()

        else:
            # as the array has not been specified we fill it with the default value
            if not self._triangular:
                self._array=[[self._default_value]*len(self._cols)
                                   for row in xrange(len(self._rows))]
            else:
                if self._triangular=='upper':
                    self._array=[[0 if col>row else self._default_value
                                  for col in self._cols] for row in self._rows]
                elif self._triangular=='lower':
                    self._array=[[0 if col<row else self._default_value
                                  for col in self._cols] for row in self._rows]
                else:
                    raise ValueError('triangular must be one of None, upper or lower\n')

                if self._unitriangular:
                    if len(self._rows)>=len(self._cols):
                        for row in self._rows:
                            if row in self._cols:self[row,row]=1
                    else:
                        for col in self._cols:
                            if col in self._rows: self[col,col]=1

            # if the default value is a ring element then we have a matrix
            if is_Element(self._default_value) and self._default_value.parent() in Rings():
                if self._base_ring is None:
                    self._base_ring=self._default_value.parent()
                self._is_matrix=True
                self._array=matrix(self._base_ring, self._array)


    def options(self, *get_option, **set_option):
        r"""
        Sets and displays the options for the array. When the options are set
        we also do some post processing to ensure that they do what we want.
        """
        # print the current settings for the options in ``option``
        if get_option==() and set_option=={}:
            print 'Current options for this array:'
            for opt in sorted(self._default_options.keys()):
                if opt in ['col_label','row_label']:
                    print '  - %s = %s' % (opt, getattr(self, '_'+opt+'_str'))
                else:
                    print '  - %s = %s' % (opt, getattr(self, '_'+opt))

        # print the current settings for the options in ``option``
        for opt in get_option:
            if not opt in self._default_options:
                raise ValueError('%s is not a valid option\n'%opt)
            elif opt[:10] in ['col_label','row_label']:
                print 'The current value of  %s = %s' % (opt, getattr(self, '_'+opt+'_str'))
            else:
                print 'The current value of %s is %s' %(opt, getattr(self,'_'+opt))

        # set the specified options, but process them in alphabetical order as
        # some options need to be set before others and conveniently
        # alphabetical order is compatible with this.
        for opt in sorted(set_option.keys()):
            if opt =='array':
                self._array=set_option[opt]
                # the documentation for ``self._matrix`` explains why this is
                # NOT self._matrix
                if is_Matrix(self._array):
                    self._is_matrix=True
                    self._base_ring=self._array.base_ring()
                    if self._array.dimensions()!=(len(self._rows),len(self._cols)):
                        raise ValueError('The sizes of the array and the row and columns indices do not agree\n')
                elif not (len(self._array)==len(self._rows) and
                          all(len(row)==len(self._cols) for row in self._array)):
                    raise ValueError('The sizes of the array and the row and columns indices do not agree\n')

            elif opt[:10] in ['col_label','row_label']:
                # we have to make a function from the name of the method
                setattr(self, '_'+opt+'_str',set_option[opt])  # store the descriptor for later use
                if set_option[opt]=='repr':
                    setattr(self, '_'+opt,lambda t: '%s'%t)
                elif set_option[opt] is None:
                    if opt=='row_label':
                        raise ValueError('row_label cannot be set to None\n')
                    setattr(self, '_'+opt,None)
                else:
                    index=getattr(self,'_'+opt[:3]+'s')
                    try:
                        rep=getattr(index[0].__class__,'_repr_'+set_option[opt])
                        setattr(self, '_'+opt, rep)
                    except AttributeError:
                        raise ValueError('%s has no method %s\n'%(index, '_'+set_option[opt]))

            elif opt in ['row_order','col_order']:
                if not (set_option[opt] is None or isinstance(set_option[opt],types.FunctionType)):
                    raise ValueError('%s must be a comparison function\n'%opt)

                setattr(self, '_'+opt, set_option[opt])
                if opt=='row_order':
                    # by default, changing the row order also changes the column order
                    #  if the rows and column indices are the same
                    if self._col_order is None and self._rows==self._cols:
                        setattr(self, '_col_order', set_option[opt])

            elif opt=='is_matrix' and hasattr(self,'_matrix'):
                # if setting _is_matrix to True we check that this is OK
                self.is_matrix(set_option[opt])

            elif not opt in self._default_options:
                raise ValueError('%s is not a valid option\n'%opt)

            else:
                setattr(self, '_'+opt, set_option[opt])

                # finally, check the options which change other options
                if opt=='default_value' and is_Element(self._default_value) and self._base_ring is None:
                    self._base_ring=self._default_value.parent()
                elif opt=='unitriangular':
                    # By default, a unitriangular array is assumed to be lower triangular.
                    if self._unitriangular and not self._triangular:
                        self._triangular='lower'


    def _current_options(self, transpose=False):
        """
        Return a dictionary of the current **non-default** options. This is
        used when constructing related :class:`LabelledMatrix` from ``self``
        such as the inverse and transpose. If the optional argument
        ``transpose`` is set to ``True`` then the options for constructing the
        transposed labelled matrix are returned.

        EXAMPLES::

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: LabelledMatrix([1,2,3],unitriangular=True,default_value=0)._current_options()
            {'default_value': 0,
            'rows': [1, 2, 3],
            'unitriangular': True,
            'triangular': 'lower',
            'is_matrix': True}
            sage: LabelledMatrix([1,2,3,4],columns=[4,5],unitriangular=True)._current_options(transpose=True)
            {'rows': [4, 5],
            'unitriangular': True,
            'triangular': 'upper',
            'columns': [1, 2, 3, 4]}
        """
        options={}
        for opt in self._default_options:
            option=getattr(self, '_'+opt)
            if transpose and opt[:4] in ['row_','col_']:
                target_opt='%s_%s'%('row' if opt[:3]=='col' else 'col', opt[4:])
            else:
                target_opt=opt
            if option!=self._default_options[target_opt]:
                if opt in ['row_label','col_label']:
                    label=getattr(self,'_'+opt+'_str')
                    if not label in [None,'repr']: options[target_opt]=label
                else:
                    options[target_opt]=option

        # the row and column labels depend partly on how the matrix was constructed
        rows=self._row_space if hasattr(self,'_row_space') else self._rows
        cols=self._col_space if hasattr(self,'_col_space') else self._cols

        if transpose:
            options['rows']=cols
            if cols!=rows: options['columns']=rows

            # transpose upper and lower triangular flags
            if 'triangular' in options:
                if options['triangular']=='lower':
                    options['triangular']='upper'
                elif options['triangular']=='upper':
                    options['triangular']='lower'

        else:
            options['rows']=rows
            if rows!=cols: options['columns']=cols

        return options


    def __repr__(self):
        """
        Returns a string for printing the array
        """
        if self._name is not None: return self._name
        elif self._rows is self._cols: return 'Array with rows and columns indexed by %s' %(self._rows,)
        else: return 'Array with rows indexed by %s and columns by %s' %(self._rows, self._cols)


    def _matrix_(self, ring=None):
        r"""
        Returns ``self._array`` so that ``matrix(self)`` will return a matrix
        representation of the array. There is a default argument of ``ring``
        which causes the entries of ``self._array`` to be coerced into the
        :class:`~sage.categories.rings.Ring` ``ring``.

        EXAMPLES::

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: d=LabelledMatrix(Partitions(3), row_label='compact_high', default_value=0, unitriangular=True);
            sage: matrix(d)
            [1 0 0]
            [0 1 0]
            [0 0 1]
            sage: matrix(d, ring=PolynomialRing(ZZ,'q'))
            [1 0 0]
            [0 1 0]
            [0 0 1]
        """
        if ring is None:
            return self._array
        else:
            return [[ring(c) for c in row] for row in self._array]

    @classmethod
    def _new_labelled_matrix(cls, *args, **kw_args):
        r"""
        Return a new labelled matrix of the specified class.

        Many instances of :class:`LabelledMatrix`

        EXAMPLES::

        """
        return cls(*args, **kw_args)

    def _display_row(self, row, row_order, col_order, width, sep):
        """
        INPUTS:
        - ``self``: the matrix
        - ``row``:  an integer for indexing the row being printed
        - ``col_order``
        - ``width`` the width of each column
        - ``sep`` the string which separates each column (typically, ``sep`` is
          equal to ` ` or `&` (for latex)

        Return a string for pretty printing the ``row``th row of the array, with
        suitable defaults.
        """
        if self._array[row]==[]:
            return sep.join('%*s' % (width, self._unknown_string))

        if self._triangular is not None:
            # For each column, col, we will need to determine whether row<col or row>col.
            # We define a function row_less to do this. Note that the ordering
            # of the rows and columns is already determined by row_order and
            # col_order.
            row_unrank=self._row_unrank(row)  # the index of row in self._rows
            if row_unrank in self._cols:
                row_col_index=col_order.index(self._col_rank(row_unrank))# position in col_order
                row_less=lambda col: row_col_index<=col_order.index(col) # use order from col_order
            else:
                # since row is not in cols, col must be in rows so we use row_order
                row_row_index=row_order.index(row)
                row_less=lambda col: row_row_index<=row_order.index(self._row_rank(self._col_unrank(col)))

        row_str=''
        for col in col_order:
            if self._array[row][col] is None:
                row_str+='%s%*s'%(sep,width,self._unknown_string)
            elif self._array[row][col]==0:
                if ((self._triangular=='lower' and row_less(col))
                        or (self._triangular=='upper' and not row_less(col))):
                  row_str+='%s%*s'%(sep,width,'')
                else:
                  row_str+='%s%*s'%(sep, width,self._zero_string)
            else:
                row_str+='%s%*s'%(sep,width,self._array[row][col])
        return row_str


    def _display(self, display_format='text'):
        """
        Returns a string for displaying/pretty printing the array either to the
        screen or as latex.

        """
        # We will return a string ``display`` for the table. First some
        # defaults which set the separators depending upon the format of the
        # output, as determined by ``display_format``.
        if display_format=='latex':
            sep=hsep='&'   # separator between columns, hsep is for the header
            bar=hbar=''    # line between row labels and row
            start_row=''   # starts every row
            end_row=r'\\'  # ends every row
            if self._name=='':
                display=''
            else:
                display='\\textbf{ %s }\n\n'%self._name
            display+='\\begin{array}{r|*{%s}c}\n'% self._col_length # start of table
            end_display='\n\\end{array}'  # the end of the display!
        elif display_format=='html':
            sep='</td><td>'
            bar='</th><td>'
            hbar='</th><th>'
            start_row='<tr><th>'
            end_row='</td></tr>'
            if self._name is None:
                display='<table>\n'
            else:
                display='<h2>%s</h2>\n<table>\n'%self._name
            end_display='\n</table>'
        else:
            sep=' '
            bar=hbar='|'
            start_row=end_row=''
            if self._name is None:
                display=''
            else:
                display='%s\n'%self._name
            end_display=''

        # the row labels
        rows=['%s' % self._row_label(row) for row in self._rows]
        # the maximum width of a row label
        row_max=max(len(row) for row in rows)

        # determine the maximum width of a column entry
        if is_Matrix(self._array):
            col_width=max(1,max(len('%s'%(a is None and self._unknown_string or a))
                                 for a in set(self._array.list())))
        else:
            col_width=max(1,max(len('%s'%(a is None and self._unknown_string or a))
                                 for a in set(flatten(list(self._array)))))

        row_order=range(len(self._rows))
        if self._row_order is not None:
            row_order.sort(cmp = lambda s,t: int(self._row_order(self._row_unrank(s),self._row_unrank(t))))
        col_order=range(len(self._cols))
        if self._col_order is not None:
            col_order.sort(cmp = lambda s,t: int(self._col_order(self._col_unrank(s),self._col_unrank(t))))

        # the column labels
        if self._col_label is not None:
            # add the row of column labels to display
            columns=['%s' % self._col_label(col) for col in self._cols]
            col_width=max(col_width, max(len(col) for col in columns))
            if display_format=='text' and self._name is not None: # ------
                display+='%s\n'%('-'*(row_max+2+(1+col_width)*len(self._cols)))
            display+='%s%*s %s %s%s\n'%(start_row,row_max,'', hbar,
                    sep.join('%*s'%(col_width,self._col_label(self._col_unrank(col)))
                              for col in col_order), end_row)
            if display_format=='text':
                display+='%s-|%s\n'%('-'*row_max,'-'*(1+col_width)*len(self._cols))
            elif display_format=='latex':
                display+='\\\\\hline\n'

        return display+'\n'.join('%s%-*s %s%s%s' %(start_row, int(row_max),
                             self._row_label(self._row_unrank(row)), bar,
                             self._display_row(row,row_order,col_order,col_width,sep), end_row)
                                 for row in row_order)+end_display

    def _latex_(self):
        """
        Return a string for printing the tableau in latex.
        """
        return self._display('latex')

    def latex(self):
        return self._latex_()

    def display(self):
        print self._display()

    def _html_display(self, filename, url_include):
        r"""
        This is a dummy function which caches the html table being constructed
        by the matrix.
        """
        if not hasattr(self, '__html_display'):
            self.__html_display=HTMLDisplay(title=self._name, html_file=filename, url_include=url_include)
        return self.__html_display

    def _html_table(self):
        r"""
        Display the matrix in a web browser.

        EXAMPLES:

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: d=LabelledMatrix(Partitions(3), default_value=0, base_ring=ZZ['q'])
            sage: d.display()
            [3]       | . . .
            [2, 1]    | . . .
            [1, 1, 1] | . . .
            sage: d.html_display()   # not tested - opens in browser
        """
        # determine row and column orderings
        row_order=range(len(self._rows))
        if self._row_order is not None:
            row_order.sort(cmp = lambda s,t: int(self._row_order(self._row_unrank(s),self._row_unrank(t))))
        col_order=range(len(self._cols))
        if self._col_order is not None:
            col_order.sort(cmp = lambda s,t: int(self._col_order(self._col_unrank(s),self._col_unrank(t))))

        # Determine the row and column labels to be used in the html table
        if self._row_labels_html is None:
            row_label=self._rows
        else:
            row_label=[self._row_labels_html(row) for row in self._rows]

        if self._col_labels_html is None:
            if self._rows==self._cols:
                col_label=row_label
            elif self._col_label is not None:
                col_label=['%s' % self._col_label(col) for col in self._cols]
            else:
                col_label=self._cols
        else:
            col_label=[self._col_labels_html(col) for col in self._cols]

        # now construct the array to be displayed, prettifying the entries as necessary
        table=[ [] for row in row_order ]
        for row in xrange(len(row_order)):
            table+=[]
            for col in xrange(len(col_order)):
                if (self._triangular=='lower' and row<col) or (self._triangular=='upper' and row>col):
                    table[row].append('')
                elif self._zero_string is not None and self[row][col]==0:
                    table[row].append(self._zero_string)
                elif self._entries_math_mode:
                    table[row].append('\({}\)'.format(latex(self[row][col])))
                else:
                    table[row].append('{}'.format(self[row][col]))

        return HTMLTable(table=table, col_label=col_label, row_label=row_label,
                         col_label_math_mode=self._col_label_math_mode,
                         row_label_math_mode=self._row_label_math_mode,
                         entries_math_mode=False  # already taken care of
                        )


    def html_display(self, filename=None, url_include=None):
        html=self._html_display(filename, url_include)
        html.add_section(heading=self._name if self._name!='' else 'Table', 
                         section=self._html_table())
        html.html_display()


    def __getattr__(self,attr):
        """
        Any unrecognised method calls/attributes get passed on to the corresponding matrix methods.
        """
        if self._is_matrix:
            try:
                return getattr(self._array,attr)
            except AttributeError:
                pass
        raise AttributeError, '%s has no attribute %s\n'%(self,attr)


    def __getitem__(self,key):
        r"""
        We overload :meth:`__getitem__` so that the entries of the array
        ``self`` can be accessed using the row and column indices. We also
        implement numpy style indexing so the entry in row ``r`` and column ``c``
        can be accessed as ``self[r,c]``.

        """
        # replace ``self[row]`` calls with ``self[row,:]``
        if not isinstance(key,tuple):
            key=(key, slice(None,None,None))

        rows,cols=key
        # first ensure that rows and cols are integers or integer slices
        if isinstance(rows,slice):
            if rows==slice(None,None,None):
                rows=slice(0,self._row_length)
            elif rows.start in self._rows:
                rows=slice(self._row_rank(rows.start),self._row_rank(rows.stop),rows.step)
        elif rows in self._rows: 
            rows=self._row_rank(rows)

        if isinstance(cols,slice):
            if cols==slice(None,None,None):
                cols=slice(0,self._col_length)
            elif cols.start in self._cols:
                cols=slice(self._col_rank(cols.start),self._col_rank(cols.stop),cols.step)
        elif cols in self._cols: 
            cols=self._col_rank(cols)

        if isinstance(rows,(int,Integer)) and isinstance(cols,(int,Integer)):
            # return a single matrix entry
            return self._array[rows][cols]

        elif isinstance(rows,slice) and isinstance(cols,slice):
            # return the submatrix as a LabelledMatrix
            options=self.__current_options()
            if 'columns' in options: options.__delitem__('columns')
            return LabelledMatrix(row=self._rows[rows], columns=self._cols[cols],
                       array=[row[cols] for row in self._array[rows]], **options)

        elif isinstance(rows,slice):
            # return the corresponding column -- by the last case, cols cannot be a slice
            if hasattr(self, '_col_sum'):
                # return an element of the column space
                return self._col_sum([(self._row_unrank(row),self._array[row][cols])
                                     for row in range(self._row_length)[rows]], distinct=True)
            else:
                return [self._array[row][cols] for row in range(len(self._rows))[rows]]

        elif isinstance(cols,slice):
            # return the corresponding row -- again, rows cannot be a slice
            if hasattr(self, '_row_sum'):
                # return an element of the row space
                return self._row_sum([(self._col_unrank(col), self._array[rows][col])
                                     for col in range(self._col_length)[cols]], distinct=True)
            else:
                return [self._array[rows][col] for col in range(self._col_length)[cols]]

        else:
            raise ValueError('%s is not a valid index for %s\n'%(key,self))


    def __setitem__(self, key, value):
        r"""
        Sets entries, rows and columns equal to ``value``. The index ``key`` can
        be either be:

        - a pair ``(r,c)`` where ``r`` and ``c`` are either integers or row and
          column  indices. This sets the corresponding ``(r,c)`` entry equal to
          ``value``

        - a single index ``r``, or slice ``(r,:)``, which is either an integer
          or a row index. This sets the entire row equal to ``value``

        - the slice ``(:,c)`` where ``c`` is an integer or a column index. This
          sets the entire row equal to ``value``

        - the slice ``:``, in which case the underlying matrix is set to
          ``value``.

        The matrix ``value`` is always assumed to have the right `shape` to be
        inserted into the matrix. It can either be a single matrix entry (which
        can be any ``sage`` object, a list or tuple of the right length to be a
        row or a column, or an element of the associated row or column space (if
        it is defined).

        EXAMPLES::

        """
        # replace ``self[row]=value`` calls with ``self[row,:]=value``
        if not isinstance(key,tuple):
            if key in self._rows:
                key=(key, slice(None,None,None))
            else:
                raise ValueError('%s is not a valid index for %s\n'%(key,self))

        row=self._row_rank(key[0]) if key[0] in self._rows else key[0]
        col=self._col_rank(key[1]) if key[1] in self._cols else key[1]
        if row==slice(None,None,None) and col==slice(None,None,None):
            self._array=value
        elif row==slice(None,None,None):
            # insert a column - but allow ``value`` to belong to self._col_space
            if hasattr(self,'_col_space') and value in self._col_space:
                unrank=lambda row: row
            else:
                unrank=lambda row: self._row_rank(row)
                if len(value)!=self._row_length:
                    raise ValueError('%s does not have the right length for a column  of %s\n'%(value,self))
            if is_Matrix(self._array):
                for row in self._rows:
                    self._array[self._row_rank(row),col]=value[unrank(row)]
            else:
                for row in self._rows:
                    self._array[self._row_rank(row)][col]=value[unrank(row)]
        elif col==slice(None,None,None):
            # insert a row - but allow ``value`` to belong to self._row_space
            if hasattr(self,'_row_space') and value in self._row_space:
                unrank=lambda col: col
            else:
                unrank=lambda col: self._col_rank(col)
                if len(value)!=self._row_length:
                    raise ValueError('%s does not have the right length to be a row %s\n'%(value,self))
            if is_Matrix(self._array):
                for col in self._cols:
                    print 'col=%s, unrank(col)=%s' %(col, unrank(col))
                    print 'col=%s' % col
                    print 'col_unrank(col)=%s' %(self._col_unrank(col))
                    self._array[row,self._col_rank(col)]=value[unrank(col)]
            else:
                for col in self._cols:
                    self._array[row][self._col_unrank(col)]=value[unrank(col)]
        elif row in xrange(self._row_length) and col in xrange(self._col_length):
            if is_Matrix(self._array):  # unfortunately, self._is_matrix is not reliable
                self._array[row,col]=value
            else:
                self._array[row][col]=value
        else:
            raise ValueError('%s is not a valid index for %s\n'%(key,self))


    def rows(self):
        """
        Return the indexing set for the rows of the array ``self``.

        EXAMPLES:

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: LabelledMatrix(Partitions(5)).rows()
            Partitions of the integer 5
        """
        return self._rows

    def columns(self):
        """
        Return the indexing set for the columns of the array ``self``.
        EXAMPLES:

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: LabelledMatrix(Partitions(5)).columns()
            Partitions of the integer 5
        """
        return self._cols

    def is_matrix(self, value=None):
        r"""
        Returns ``True`` or ``False`` depending upon whether the underling
        array in ``self`` is a matrix. With no option

        EXAMPLES:

            sage: from sage.matrix.labelled_matrix import LabelledMatrix
            sage: S=LabelledMatrix(rows=Partitions(5), unitriangular=True, default_value=0)
            sage: S.is_matrix()
            True
            sage: S.is_matrix(False)
            False
            sage: S.is_matrix('check')
            True
            sage: tri=LabelledMatrix(rows=Partitions(5), unitriangular=True)
            sage: tri.is_matrix()
            False
            sage: tri.display()
            [5]             | 1
            [4, 1]          | ? 1
            [3, 2]          | ? ? 1
            [3, 1, 1]       | ? ? ? 1
            [2, 2, 1]       | ? ? ? ? 1
            [2, 1, 1, 1]    | ? ? ? ? ? 1
            [1, 1, 1, 1, 1] | ? ? ? ? ? ? 1
            sage: tri.is_matrix('check')
            False
            sage: tri.is_matrix(True)
            Traceback (most recent call last):
            ...
            ValueError: This array is NOT a matrix
        """
        if value is False:
            self._is_matrix=value

        elif isinstance(value,bool):
            self._is_matrix=value
            if self._is_matrix and not is_Matrix(self._array):
                # check that it is actually a matrix
                try:
                    self._array=matrix(self._array)
                except (TypeError,ValueError):
                    raise ValueError('This array is NOT a matrix\n')

        elif value=="check":
            self._is_matrix=is_Matrix(self._array)

        return self._is_matrix


    def transpose(self):
        r"""
        Return the :class:`LabelledMatrix` obtained from ``self`` by
        taking the transpose of the underlying array.
        """
        options=self._current_options(transpose=True)
        if is_Matrix(self._array):
            return self._new_labelled_matrix(array=self._array.transpose(), **options)
        else:
            options['array']=[[self._array[r][c] for r in xrange(len(self._rows))] for c in xrange(len(self._cols))]
            return self._new_labelled_matrix(**options)


    def inverse(self):
        r"""
        If ``self.is_matrix()==True`` the inverse of the matrix is computed
        and returned as a :class:`LabelledMatrix`.
        """
        if self._is_matrix:
            try:
                return self._new_labelled_matrix(matrix=self._array.inverse(), **self._current_options())
            except ZeroDivisionError:
                raise ZeroDivisionError('input matrix must be nonsingular\n')
        else:
            raise ValueError("array is not known to be a matrix. Try self.is_matrix('check')\n")


    def blocks(self):
        """
        Return a list of the *blocks*, or *linkage classes*, of the matrix.
        """
        if not is_Matrix(self._array) or self._row_length!=self._col_length:
            raise TypeError('blocks only make sense for square matrices')
        G = Graph({i: self._array[i].nonzero_positions() for i in xrange(self._array.nrows())})
        blocks=[]
        options=self._current_options()
        for comp in G.connected_components():
            options['rows']=[self._rows[i] for i in comp]
            options['array']=matrix([[self._array[r][c] for c in comp] for r in comp])
            blocks.append(self._new_labelled_matrix(**options))
        return blocks


#--------------------------------------------------
# Example application of LabelledMatrix class
#--------------------------------------------------
class CharacterTable(LabelledMatrix):
    """
    This is a toy example which wraps the character table of a finite
    group inside a :class:`LabelledMatrix`. The irreducible characters are
    labelled as `X0`, `X1`, ..., `Xn`, where `n` is the number of conjugacy
    classes and the conjugacy classes are labelled as '<order><letter>', where
    `order` is the order of any element of the conjugacy class and `letter` runs
    through the letters of the alphabet. The order of the irreducible caraters
    and the conjugacy classes is as determined by
    :meth:

    EXAMPLES:

        sage: from sage.matrix.labelled_matrix import CharacterTable
        sage: ct=CharacterTable( SymmetricGroup(4) ); ct.display()
        Character table of Symmetric group of order 4! as a permutation group
        -------------------
           | 1a 2a 2b 3a 4a
        ---|---------------
        X0 |  1 -1  1  1 -1
        X1 |  3 -1 -1  0  1
        X2 |  2  0  2 -1  0
        X3 |  3  1 -1  0 -1
        X4 |  1  1  1  1  1
        sage: ct['X0'] == ct[0]
        True
        sage: ct['X1']+2*ct['X3']  # unfortunately these are lists not characters!
        [3, -1, -1, 0, 1, 3, 1, -1, 0, -1, 3, 1, -1, 0, -1]

    Some work needs to be done to make this into a useful class. For example,
    some character theoretic methods should be added to the class, such as
    methods to decompose an arbitrary character into a linear combination of
    irreducible characters.  It would also be good to define the row space of a
    character table to be the character ring. Finally, as in Gap, it would be
    good to find a sensible way to display non-integral character values.
    """
    def __init__(self, G):
        r"""
        Returns the character table for the group ``G`` wrapped as a
        :class:`LabelledMatrix`.

        EXAMPLES::

            sage: from sage.matrix.labelled_matrix import CharacterTable
            sage: CharacterTable( DihedralGroup(4) ).display()
            Character table of Dihedral group of order 8 as a permutation group
            -------------------
               | 1a 2a 2b 4a 2c
            ---|---------------
            X0 |  1  1  1  1  1
            X1 |  1 -1 -1  1  1
            X2 |  1 -1  1 -1  1
            X3 |  1  1 -1 -1  1
            X4 |  2  0  0  0 -2
        """
        # The conjugacy classes are labelled as <order><letter>, where ``order``
        # is the order of a element of the conjugacy class and ``letter`` is a
        # letter of the alphabet
        alphabet='abcdefghijklnmopqrstuvwxyz'
        cc_orders=[c.order() for c in G.conjugacy_classes_representatives()]
        cc_labels=['%d%s' %(cc_orders[c], alphabet[cc_orders[:c].count(cc_orders[c])])
                                   for c in xrange(len(cc_orders))]
        irreds=['X%d' % irr for irr in xrange(len(cc_labels))]
        super(CharacterTable, self).__init__(rows=irreds,
                                             columns=cc_labels,
                                             array=G.character_table(),
                                             is_matrix=True,
                                             zero_string='0',
                                             name='Character table of %s'%G)

