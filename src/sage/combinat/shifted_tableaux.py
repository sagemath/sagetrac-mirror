from sage.structure.element import Element
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent import Parent
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.functions.other import factorial
from sage.misc.misc_c import prod
from sage.misc.prandom import randrange

class ShiftedTableau(Element):
    def __init__(self, parent, tableau):
        """Initialize the shifted tableau with a list of lists."""
        Element.__init__(self,parent)
        self.tableau = tableau
#        self._shape = [len(row) for row in tableau]

    def __getitem__(self, i):
        """
        Return ``i``-th row of ``self``.

        EXAMPLES::
        
            sage: T = ShiftedTableaux([4,2])
            sage: T([[1,2,4,5],[3,6]])[1]
            [3, 6]
        """
        return self.tableau[i]

    def __len__(self):
        """
        Return the length of ``self``.

        EXAMPLES::
        
            sage: T = ShiftedTableaux([4,2])
            sage: len(T([[1,2,4,5],[3,6]]))
            2
        """
        return len(self.tableau)

    def _repr_(self):
        """
        Return a string representation of ``self``.

        EXAMPLES::
        
            sage: T = ShiftedTableaux([4,2])
            sage: T([[1,2,4,5],[3,6]])
            [[1, 2, 4, 5], [3, 6]]
        """        
        return repr(self.tableau)


    def _latex_(self):
        """
        Return LaTex code for ``self``.

        EXAMPLES::
        
            sage: T = ShiftedTableaux([4,2])
            sage: latex(T([[1,2,4,5],[3,6]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}
            {\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{2}&\lr{4}&\lr{5}\\\cline{1-4}
            &\lr{3}&\lr{6}\\\cline{2-3}
            \end{array}$}
            }
        """         
        from sage.combinat.output import tex_from_array
        L = [[None]*i + row for (i,row) in enumerate(self)]
        return tex_from_array(L)

    def pp(self):
        """Print out a nice version of ``self``.
        
        EXAMPLE::
        
            sage: T = ShiftedTableaux([4,2])
            sage: T([[1,2,4,5],[3,6]]).pp() 
            1  2  4  5
               3  6
            
        """
        for row_num, row in enumerate(self.tableau):
            print "   "*row_num + "".join(repr(n).rjust(3) for n in row)

    def yang_baxter_moves(self):
        """Count the number of Yang-Baxter-type configurations in ``self``.
        These are configurations of the form::

                     i  i+1
                        i+2
                        
        along the superdiagonal.
        
        EXAMPLE::
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.yang_baxter_positions()
            1
            
        """
        return sum( self[r+1][0] == self[r][0] + 2
                    for r in range(len(self) - 1) )

    def yang_baxter_positions(self):
        """Return a vector of ``True`` and ``False`` values giving the locations of
        Yang-Baxter-type configurations in ``self``.
        
        EXAMPLE::

            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.yang_baxter_positions()
            [False, True, False]
        
        """

        return [ self[r+1][0] == self[r][0] + 2
                    for r in range(len(self) - 1) ]

    def position(self,value):
        """Return the coordinates of the cell of ``self`` labelled value.
        
        EXAMPLE::
            
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.position(6)
            (1,1)
            sage: t.position(12)
            (3,0)
            sage: t.position(22)
            Exception: Didn't find
        
        """
        for i in range(len(self)):
            if value in self[i]:
                return (i,self[i].index(value))
        raise Exception("Didn't find")
        
    def shape(self):
        """Return the shape of ``self``.
        
        EXAMPLE::
            
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: t = T([[1, 2, 3, 4, 9, 13], [5, 6, 8, 11], [7, 10, 14], [12]])
            sage: t.shape()
            [6, 4, 3, 1]
            
        """
        return self._shape


class ShiftedTableaux(UniqueRepresentation,Parent):
    @staticmethod
    def __classcall_private__(cls, shape):
        shape = tuple(shape)
        return super(ShiftedTableaux, cls).__classcall__(cls, shape)
    
    def __init__(self, shape):
        self._shape = shape
        Parent.__init__(self, category=FiniteEnumeratedSets())
    
    def __iter__(self):
        """ Iterative over self.
        
        EXAMPLE::
            
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T[0:4]
            [[[1, 2, 3, 4, 5, 6], [7, 8, 9, 10], [11, 12, 13], [14]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 10], [11, 12, 14], [13]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 11], [10, 12, 13], [14]],
             [[1, 2, 3, 4, 5, 6], [7, 8, 9, 11], [10, 12, 14], [13]]]
            
        """
        P = shape_to_poset(self._shape)
        L = P.linear_extensions()
        for l in L:
            t = linear_extension_to_tableau(l,self._shape)
            yield self.element_class(self,t)
    
    def cardinality(self):
        """ Return the number of shifted tableaux with the shape of ``shape``.
        This uses the shifted hook-length formula originally due to Thrall in 1952.
        
        EXAMPLE::
            
            sage: T = ShiftedTableaux([4,2])
            sage: T.cardinality()
            5
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.cardinality()
            1716
            
        """
        
        n = sum(self._shape)
        positions = []
        for i in range(len(self._shape)):
            for j in range(self._shape[i]):
                positions.append((i,j))
        hook_lengths = [len(hook(self._shape,pos))+1 for pos in positions]
        return factorial(n)/prod(hook_lengths)
    
    def random_element(self):
        """Return a random shifted tableau with the shape of ``self``. 
        This uses a variant of the Greene-Nijenhuis-Wilf hook-walk algorithm due to Sagan 1980.
        
        EXAMPLE::
            
            sage: T = ShiftedTableaux([4,2])
            sage: T.random_element()
            [[1, 2, 4, 5], [3, 6]]
            sage: T = ShiftedTableaux([6,4,3,1])
            sage: T.random_element()
            [[1, 2, 3, 5, 7, 10], [4, 6, 9, 12], [8, 11, 14], [13]]
            
        """
        tableau = [ [0] * row_length for row_length in self._shape ]
        N = sum(self._shape)
        current_shape = list(self._shape)
        for next_number in reversed(range(1, N+1)):
            cell = rand_cell(current_shape)
            hook_coords = hook(current_shape, cell)
            while hook_coords:
                cell = hook_coords[ randrange(len(hook_coords)) ]
                hook_coords = hook(current_shape, cell)
            # done with hook walk, so place next_number in tableau
            tableau[cell[0]][cell[1]] = next_number
            current_shape[cell[0]] -= 1

        return self.element_class(self,tableau)
        
        
    Element = ShiftedTableau

    
def hook(shape, pos):
    """ Determines the coordinates of the cells in the hook corresponding to the cell indexed by the tuple ``pos`` in the shifted shape ``shape``.
    
    EXAMPLES::
    
        sage: hook([6,4,3,1],(1,1))
        [(1, 2), (1, 3), (2, 0), (3, 0)]
    
    """
    # arm first:
    hook_coordinates = [(pos[0], col) for col in range(pos[1]+1, shape[pos[0]])]

    # now the leg:
    (r, c) = pos
    r, c = r+1, c-1
    while c >= 0 and shape_contains(shape, (r, c)):
        hook_coordinates.append( (r, c) )
        r, c = r+1, c-1
    if c < 0 and r < len(shape):
        hook_coordinates.extend( (r, j) for j in range(shape[r]) )
    return hook_coordinates    
    
def rand_cell(shape):
    """Return a uniformly random cell from the given shifted tableau shape.
    
    EXAMPLES:
    
        sage: [rand_cell([6,4,3,1]) for i in range(6)]
        [(2, 2), (1, 1), (1, 2), (1, 1), (2, 1), (3, 0)]
    
    """

    rnd = randrange(0, sum(shape))
    row = 0
    cells_so_far = 0
    while cells_so_far + shape[row] <= rnd:
        cells_so_far += shape[row]
        row += 1
    return (row, rnd - cells_so_far)
    
def shape_contains(shape, pos):
    """ Return ``True`` if the given shifted shape contains the given position.
    Here, ``pos`` must have both coordinates positive.
    
    EXAMPLES::
    
        sage: shape_contains([6,4,3,1],(2,2))
        True
        sage: shape_contains([6,4,3,1],(4,2))
        False
        
    """
    
    if pos[0] >= len(shape):
        return False
    return pos[1] < shape[pos[0]]
    
def shape_to_poset(shape):
    """ Converts the shifted shape ``shape`` into a SageMath poset with elements `1,2,...,n`.
    
    EXAMPLES::
        
        sage: shape_to_poset([4,2]).cover_relations()
        [[1, 2], [2, 3], [2, 5], [3, 4], [3, 6], [5, 6]]      
    
    """
    n = sum(shape)
    elts = range(1, n+1)  # elements of the poset
    cells = [] # elements of the poset as tuples
    for i in range(len(shape)):
        cells.extend([(i,j) for j in range(shape[i])])
    rels = []    # poset relations
    for p in range(n):
        for q in range(n):
            if cells[p][0] == cells[q][0]:  # same row
                if cells[p][1] < cells[q][1]: # p before q
                    rels.append((p+1,q+1))
            if sum(cells[p]) == sum(cells[q]): # same column
                if cells[p][0] < cells[q][0]: # p above q
                    rels.append((p+1,q+1))
    return Poset((elts, rels), cover_relations = False)
    
def linear_extension_to_tableau(linear_extension, shape):
    """ Converts a linear extension for the poset ``shape_to_poset(shape)`` into a shifted tableau of shape ``shape``.
    
    EXAMPLES::
        
        sage: P = shape_to_poset([4,2])
        sage: [linear_extension_to_tableau(L,[4,2]) for L in P.linear_extensions()]
        [[[1, 2, 3, 4], [5, 6]],
         [[1, 2, 3, 5], [4, 6]],
         [[1, 2, 3, 6], [4, 5]],
         [[1, 2, 4, 5], [3, 6]],
         [[1, 2, 4, 6], [3, 5]]]
    
    """
    
    L = list(Permutation(linear_extension).inverse()) # tableau in a list, rows not separate
    A = [sum(shape[:i]) for i in range(len(shape)+1)] # determine indices where new rows start
    return [L[A[i]:A[i+1]] for i in range(len(A)-1)]