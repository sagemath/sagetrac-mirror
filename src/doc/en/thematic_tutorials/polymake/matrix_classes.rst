.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Dealing with matrices and vectors
---------------------------------

The following list gives an overview over the operators and basic
constructions for working with matrices in ``polymake``. Some of them
also apply to vectors. Operators and functions are described in
``C++``-syntax, but most of them can also be used in the interactive
shell on ``Perl`` side, compare
`release_docs/common <http://polymake.org/release_docs/2.12/common.html>`__.

Operations
~~~~~~~~~~

C++
^^^

::

   const GenericMatrix operator- (const GenericMatrix&); 
   const GenericMatrix operator+ (const GenericMatrix&, const GenericMatrix&); 
   const GenericMatrix operator- (const GenericMatrix&, const GenericMatrix&); 
   const GenericMatrix operator* (const Scalar&, const GenericMatrix&); 
   const GenericMatrix operator* (const GenericMatrix&, const Scalar&); 
   const GenericMatrix operator/ (const GenericMatrix&, const Scalar&); 
   const GenericMatrix operator* (const GenericMatrix&, const GenericMatrix&); 
   const GenericMatrix operator* (const GenericMatrix&, const GenericVector&); 
   const GenericMatrix operator* (const GenericVector&, const GenericMatrix&);

Do exactly what the mnemonics suggest. Every type combination is
allowed, including mixing of sparse and dense vectors and matrices.
ElementType of matrix and vector operands, and Scalar type can be
different, provided the arithmetic operator with corresponding arguments
exists.

The dimensions of the operands must match.

::

   _Matrix& GenericMatrix::negate(); 
   _Matrix& GenericMatrix::operator+= (const GenericMatrix&); 
   _Matrix& GenericMatrix::operator-= (const GenericMatrix&); 
   _Matrix& GenericMatrix::operator*= (const Scalar&); 
   _Matrix& GenericMatrix::operator/= (const Scalar&);

Corresponding assignment versions of the arithmetic operators.

::

   bool operator== (const GenericMatrix&, const GenericMatrix&); 
   bool operator!= (const GenericMatrix&, const GenericMatrix&); 
   bool operator< (const GenericMatrix&, const GenericMatrix&); 
   bool operator> (const GenericMatrix&, const GenericMatrix&); 
   bool operator<= (const GenericMatrix&, const GenericMatrix&); 
   bool operator>= (const GenericMatrix&, const GenericMatrix&);

Compare two matrices lexicographically (iterating over elements in row
order). The dimensions of the operands must match.

::

   GenericMatrix::operator bool() const; bool GenericMatrix::operator! () const;

Look for a non-zero element.

Perl
^^^^

The operators given above in ``C++``-syntax can be translated to
``Perl`` in the following way. Consider e.g. the operator

::

   const GenericMatrix operator+ (const GenericMatrix&, const GenericMatrix&);

In ``Perl`` you can apply this operator via


::

    polymake> $m1=new Matrix(3,3); $m2=new Matrix(3,3);
    ........> $m1+=$m2; # changing the variable $m1

or


::

    polymake> $m1=unit_matrix(4); $m2=3*unit_matrix(4);
    ........> $m=$m1+$m2; print $m; print dense($m);

All other operators are translated in a similar way.

Other constructions
~~~~~~~~~~~~~~~~~~~

C++
^^^

::

   GenericMatrix GenericMatrix::minor(const RowIndexSet& row_indices, const ColIndexSet& column_indices);

Select a matrix minor lying on the intersection of the given row and
column subsets. The const variant of this method creates an immutable
minor object.

The indices must lie in the valid ranges.

::

   GenericVector GenericMatrix::diagonal (int i=0); 
   GenericVector GenericMatrix::anti_diagonal (int i=0);

Select the diagonal or anti-diagonal of a matrix:

-  i=0 the main diagonal (anti-diagonal)

-  i>0 the i-th diagonal below the main

-  i<0 the (-i)-th diagonal above the main The const variants of these
   methods create immutable slice objects.

   GenericMatrix T (GenericMatrix&);

Transpose the matrix. The roles of Rows and Cols get swapped. (**Perl:**
``transpose()``)

::

   GenericMatrix operator/ (const GenericMatrix&, const GenericMatrix&); 
   GenericMatrix operator/ (const GenericVector&, const GenericMatrix&); 
   GenericMatrix operator/ (const GenericMatrix&, const GenericVector&); 
   GenericMatrix operator/ (const GenericVector&, const GenericVector&);

Create a block matrix, virtually appending the rows of the second matrix
after the last row of the first matrix. A vector argument is treated as
a matrix with one row. Column dimensions of the operands must be equal.
ElementType of the operands must be identical, otherwise you will get a
rather cryptic message from the compiler.

::

   GenericMatrix operator| (const GenericMatrix&, const GenericMatrix&); 
   GenericMatrix operator| (const GenericVector&, const GenericMatrix&); 
   GenericMatrix operator| (const GenericMatrix&, const GenericVector&);

Create a block matrix, virtually appending the columns of the second
matrix after the last column of the first matrix. A vector argument is
treated as a matrix with one column. Row dimensions of the operands must
be equal.

ElementType of the operands must be identical, otherwise you will get a
rather cryptic message from the compiler.

::

   Matrix& Matrix::operator/= (const GenericMatrix&); 
   Matrix& Matrix::operator/= (const GenericVector&); 
   Matrix& Matrix::operator|= (const GenericMatrix&); 
   Matrix& Matrix::operator|= (const GenericVector&);

Append rows or columns to the matrix. All constraints for the
non-assigning operators apply here too.

::

   const GenericMatrix diag(const GenericMatrix&, const GenericMatrix&); 
   const GenericMatrix diag(const GenericMatrix&, const GenericVector&); 
   const GenericMatrix diag(const GenericVector&, const GenericMatrix&);

Create a block-diagonal matrix. Vector arguments are treated as square
diagonal matrices.

::

   GenericMatrix vector2row(const GenericVector& v); 
   GenericMatrix vector2col(const GenericVector& v);

Disguise v as a matrix with 1 row (column).

::

   const GenericMatrix diag(const GenericVector& v);

Create a square diagonal matrix. The elements of v appear on the main
diagonal.

::

   const GenericMatrix repeat_row(const GenericVector& v, int n); 
   const GenericMatrix repeat_col(const GenericVector& v, int n);

Create a matrix with n rows (columns), each equal to v.

::

   const GenericMatrix same_element_matrix(const ElementType& x, int m, int n); 
   const GenericMatrix ones_matrix<ElementType>(int m, int n); 
   const GenericMatrix zero_matrix<ElementType>(int m, int n);

Create a dense matrix with m rows and n columns, with all elements equal
to x, 1, or 0, respectively. Note the obligatory explicit specification
of the template parameter ElementType, where it cannot be deduced from
the function arguments.

::

   const GenericMatrix same_element_sparse_matrix<ElementType>(const GenericIncidenceMatrix& M); 
   const GenericMatrix same_element_sparse_matrix(const GenericIncidenceMatrix& M, const ElementType& x);

Create a sparse matrix with the same structure as M, with elements equal
to x in true cells and 0 in false cells. If omitted, x is taken equal to
1; note the obligatory explicit template parameter specification in this
case.

::

   const GenericMatrix unit_matrix<ElementType>(int n);

Create a unit matrix n×n.

Perl
^^^^

Many of the ``C++``-functions also exist on ``Perl`` side (under the
same name), compare
`release_docs/common <http://polymake.org/release_docs/2.12/common.html>`__.
To append rows or columns in ``Perl``, you can also use the operators
``/`` and ``|``:


::

    polymake> $m1 = new Matrix(3,4); $v1 = new Vector(4);
    ........> $m1_ext = $m1/$v1;
    ........> $m2 = new Matrix(2,4);
    ........> $m = $m1_ext / $m2; print $m;
    ........> $zero_vec6 = ones_vector<Rational>(6);
    ........> print $zero_vec6|$m;


