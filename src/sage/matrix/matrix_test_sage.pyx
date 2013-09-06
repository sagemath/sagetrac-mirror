

from sage.matrix.matrix_generic_dense cimport Matrix_generic_dense
from sage.categories.discrete_valuation import DiscreteValuationRings
#from sage.matrix.matrix0 import add_multiple_of_row
from copy import copy

cdef class Matrix_test_sage(Matrix_generic_dense):# dense or sparse
    cpdef test_func(self):
        return self[0][0].valuation()
    
    cpdef row_echelon_test(self, pivot_choice = False,  l = None, compute_pivots = False):
        #M = self.copy()
        M = copy(self)
        entries = M
        m = self.nrows()
        n = self.ncols()
        nmin = min(m,n)
        rows = range(m)
        if compute_pivots: rows_pivot = [ ]
        j = 0
        if l is None:
            l = n
        while (j<l) and len(rows)>0:           
                #looking for a first real non-zero term in the i column
                pivot_found = False
                for k0 in rows:
                    entry = entries[k0][j]
                    if entry <> 0:
                        # first time we find a pivot. Its valuation is the smallest
                        # one for now.
                        if not(pivot_found):
                            minval = entry.valuation()
                            rpiv = k0
                            piv = entry
                        
                        pivot_found = True
                        if not(pivot_choice):
                            piv = entry
                            rpiv = k0
                            minval = entry.valuation()
                            break
                        else:
                            curval = entry.valuation()
                            if curval < minval:
                                piv = entry
                                rpiv = k0
                                minval = curval
                #performing the elimination
                if not(pivot_found):
                    print("Watch out, a problem appeared with false zeroes")
                else:
                    rows.remove(rpiv)
                    if compute_pivots: 
                        rows_pivot.append(rpiv)
                    if not(pivot_choice):
                        for k1 in range(rpiv+1,m):
                           if entries[k1][j] <> 0:
                                M.add_multiple_of_row(k1, rpiv, -entries[k1][j]/piv)
                    else:    
                        for k1 in rows:
                            #here might happen the loss of precision...
                            if entries[k1][j] <> 0:
                                M.add_multiple_of_row(k1, rpiv, -entries[k1][j]/piv)
                j += 1
                
        return M
    
    cdef alter_precision(self,i,j,col_ind):
        val_piv = self.get_unsafe(i,col_ind).valuation()
        for k in range(col_ind,self.ncols()):
            val_xk = self.get_unsafe(i,k).valuation()
            new_prec = self.get_unsafe(j,k).precision_absolute()-val_piv+val_xk
            self.set_unsafe(j,k,self.get_unsafe(j,k).add_bigoh(new_prec))
    
    cpdef row_echelon_test2(self , l = None, compute_pivots = False):
        #M = self.copy()
        M = copy(self)
        entries = M
        m = self.nrows()
        n = self.ncols()
        nmin = min(m,n)
        rows = range(m)
        if compute_pivots: rows_pivot = [ ]
        j = 0
        if l is None:
            l = n
        while (j<l) and len(rows)>0:           
                #looking for a first real non-zero term in the i column
                pivot_found = False
                for k0 in rows:
                    entry = entries[k0][j]
                    if entry <> 0:
                        if not(pivot_found):
                            minval = entry.valuation()
                        
                        pivot_found = True

                        curval = entry.val()
                        if curval < minval:
                            piv = entry
                            rpiv = k0
                            minval = curval
                #performing the elimination
                if not(pivot_found):
                    print("Watch out, a problem appeared with false zeroes")
                else:
                    rows.remove(rpiv)
                    if compute_pivots: 
                        rows_pivot.append(rpiv)
    
                    for k1 in rows:
                        #here might happen the loss of precision...
                        M.alter_precision(rpiv,k1,j)
                        M.set_unsafe(k1,j,0)
                                
                j += 1
                
        return M

 #   def __init__(self, parent):
 #       self._parent = parent
 #       self._base_ring = parent.base_ring()
 #       self._nrows = parent.nrows()
 #       self._ncols = parent.ncols()

 

