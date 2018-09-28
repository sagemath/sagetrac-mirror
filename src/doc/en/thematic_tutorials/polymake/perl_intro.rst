.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


polymake for the Perl Newbie
============================

The language that the interactive version of ``polymake`` speaks is a
dialect of Perl that we refer to as ``polymake``/Perl. See
`www.perl.org <http://www.perl.org>`__ for comprehensive Perl
information. Note also that the ordinary Perl manual pages are
particularly useful, especially the perlintro man page which is also
available on `perldoc <http://perldoc.perl.org/perlintro.html>`__. This
short section here cannot be a replacement for a thorough introduction
to this language, but we want to focus on a few key points that are
relevant to ``polymake``.

Standard data structures
~~~~~~~~~~~~~~~~~~~~~~~~

The Perl programming language originally provides three different data
structures,
scalars(\ `), arrays(@), and hashes(%). The user always has to specify the type of a variable using the appropriate symbol ``\ ``,``\ @\ ``, or``\ %`.
If you forget to do so, you will receive the following error message:

::

    polytope > i=5;
   polymake:  ERROR: Unquoted string "i" may clash with future reserved word.
   </code>


   Here are some simple commands illustrating how to use the different data structures:
   ==Scalars==
   <code>


::

    polymake> $i=5;
    ........> $j=6;
    ........> $sum=$i+$j; print $sum;

Arrays
''''''


::

    polymake> @array=("a","b","c"); print scalar(@array);
    ........> push(@array,"d"); print "@array"; 
    ........> $first_entry=$array[0]; print $first_entry;
    ........> print join("\n",@array);
    ........> @array2=(3,1,4,2);
    ........> print sort(@array2);

Hashes
''''''


::

    polymake> %hash=();
    ........> $hash{"zero"}=0;
    ........> $hash{"four"}=4;
    ........> print keys %hash;
    ........> print join(", ",keys %hash);
    ........> print join(", ",values %hash);
    ........> %hash=("one",1,"two",2);
    ........> %hash=("one"=>1,"two"=>2);

``polymake``-Perl
~~~~~~~~~~~~~~~~~

In addition to the three standard data structures, the enriched version
of ``Perl`` used in ``polymake`` also provides special data structures
for dealing with more complicated structures. ``polymake``\ ’s object
hierarchy is completely reflected on the Perl side. Let us create a
small polytope as an example object.


::

    polymake> $p = new Polytope(POINTS=>[[1,0,1],[1,0,-1],[1,1,0],[1,-1,0]]);

Note that the ``Perl``-type of the variable ``$p`` is ``Scalar``, as the
variable is internally treated as a reference to a ``C++``-object. The
true nature of the object becomes visible if it is printed:


::

    polymake> print $p;
    Polymake::polytope::Polytope__Rational=ARRAY(0x2f2f1c0)
    





In this case it is a ``polymake`` object from the application
``polytope``, and it happens to be of type \`Polytope. Technically, $p
is a reference to an array (but it should be never treated as an array
unless you are deliberately trying to crash polymake). If you want less
technical information on the type of your object, use this:


::

    polymake> print $p->type->full_name;
    Polytope<Rational>
    





“Small objects”: Data structures inherited from C++
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use objects that are inherited from the ``C++``-side of
``polymake`` in the interactive shell. A complete list of so-called
“small objects” can be found in the `online
documentation </release_docs/latest/common.html>`__ under the heading
“Property types”. Here is a selection of three different structures that
facilitate everyday work with ``polymake``: ##### Arrays

The small object ``Array`` can be initialized in different ways and with
different template parameters:


::

    polymake> @array=("a","b","c");
    ........> $arr1=new Array<String>(\@array); print $arr1;
    ........> $arr2=new Array<Int>([3,2,5]); print $arr2;
    ........> $arr3=new Array<Int>(0,1,2,3); print $arr3;
    ........> $arr4=new Array<Int>(0..4); print $arr4;
    ........> $arr5=new Array<Int>($arr4); print $arr5;

You have random access:


::

    polymake> $arr5->[0] = 100;
    ........> print $arr5;
    100 1 2 3 4
    





It is also possible to convert the ``C++``-object ``Array`` into a
``Perl``-array by writing > @arr4=@{$arr4}; print
`arr2; </code> or simply<code> > @arr4=@`\ arr4; ##### Sets

On ``C++``-side sets are stored in a balanced binary search (AVL) tree.
For more information see the
`PTL-documentation <https///polymake.org/release_docs/master/PTL/classpm_1_1Set.html>`__.
In many cases, the small objects can be converted into ``Perl``-types in
the expected way: > $set=new Set(3,2,5); print $set; > print
`set->size; > @array_from_set=@`\ set; ##### Matrices

Here is a simple way to initialize a matrix: > $mat=new
Matrix([2,1,4,0,0],[3,1,5,2,1],`1,0,4,0,6 <2,1,4,0,0%5D,%5B3,1,5,2,1%5D,%5B1,0,4,0,6>`__);
> print $mat; You could also define it by passing a reference to an
(``Perl``-)array of ``Vectors``. The single entries are interpreted as
different rows: > $row1=new Vector([2,1,4,0,0]); > $row2=new
Vector([3,1,5,2,1]); >
`row3=new Vector<Rational>([1,0,4,0,6]); > @matrix_rows=(`\ row1,\ `row2,`\ row3);
> $matrix_from_array=new Matrix(@matrix_rows); You can change a single
entry of a matrix in the following way (if it is not already assigned to
an immutable property like ``VERTICES``!): > $mat->row(1)->[1]=7; print
$mat->row(1)->[1]; > print $mat; > $mat->(1,2)=8; print $mat; A unit
matrix of a certain dimension can be defined via the user-function
``unit_matrix<COORDINATE_TYPE>(.)``: > $unit_mat=4\ *unit_matrix(3); >
print
`unit_mat; </code> The reason for the "strange output" is the implementation as *sparse matrix*: <code> > print ref(`\ unit_mat);
However, some functions cannot deal with this special type of matrix. In
this case it is necessary to transform the sparse matrix into a dense
matrix first via: > `dense=new Matrix<Rational>(`\ unit_mat);print
$dense; or just > `dense2=dense(`\ unit_mat);print $dense2; You
can also work with matrices that have different types of coordinates
like ``Rational``, ``Float``, or ``Int``: > $m_rat=new
Matrix(3/5*\ unit_matrix(5)); print $m_rat; > `m2=`\ mat/$m_rat;
print $m2; > $m_int=new Matrix(unit_matrix(5)); print $m_int; >
`m3=`\ m_rat/$m_int; The error message polymake: ERROR: undefined
operator Matrix / Matrix at input line 1. indicates that you need to
convert the integer matrix to a rational matrix first: >
`m3=`\ m_rat/convert_to($m_int); print $m3; By “/” you can add
rows to a matrix, whereas “\|” adds columns. By the way, this also works
for ``Vector``. > `z_vec=zero_vector<Int>(`\ m_int->rows); >
`extended_matrix=(`\ z_vec|$m_int); print $extended_matrix; It is
also possible to nest template parameters in any way you like, e.g.


::

    polymake> $set=new Set<Int>(3,2,5);
    ........> $template_Ex=new Array<Set<Int>>(new Set<Int>(5,2,6), $set); print $template_Ex; print ref($template_Ex);

However, if you use a template combination, you have never used before,
it may take some time until you see the result. This is due to the fact
that ``polymake`` compiles your new combination *on the fly*. But this
is only a one-time effect, and next time you use this combination it
will work without delay.

“Big Objects”: Objects with properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A big object is an instance of a data type which represents a
mathematical concept with clear semantics. They may have template
parameters.


::

    polymake> $p=new Polytope<Rational>(POINTS=>cube(4)->VERTICES);
    ........> $lp=new LinearProgram<Rational>(LINEAR_OBJECTIVE=>[0,1,1,1,1]);

Big objects have properties which come with a type, which is either
built-in or a small object type or a big object type, and which can be
accessed using the :literal:`-``>` operator.


::

    polymake> # access the property named `LP`:
    ........> $p->LP=$lp;
    ........> # properties can have properties themselves.
    ........> print $p->LP->MAXIMAL_VALUE;

Scalar properties can be used in arithmetic expressions right away.


::

    polymake> $i = ($p->N_FACETS * $p->N_FACETS) * 15;




::

    polymake> print $i;
    2940
    





Check out the tutorial on `properties <properties>`__ to learn more
about the way properties are used and computed.

A small example script…
~~~~~~~~~~~~~~~~~~~~~~~

…to demonstrate the usage of ``polymake``/Perl. You can download the
matrix file {{:points.demo\| here}}.


::

    polymake> ### load matrix from file
    ........> open(INPUT, "< demo/Workshop2011/points.demo");
    ........> $matrix=new Matrix<Rational>(<INPUT>);
    ........> close(INPUT);
    ........> print $matrix;
    > 





::

    polymake> ### create a polytope from the matrix
    ........> $p=new Polytope<Rational>(POINTS=>$matrix);
    ........> print $p->FACETS;
    ........> print $p->DIM;
    ........> print $p->VERTEX_SIZES;
    > 





::

    polymake> ### print "simple" vertices
    ........> for(my $i=0;$i<scalar(@{$p->VERTEX_SIZES});$i++){
    ........>     if($p->VERTEX_SIZES->[$i]==$p->DIM){
    ........>     print $i.": ".$p->VERTICES->row($i)."\n";
    ........>     }
    ........> }
    > 





::

    polymake> ### put their indices in a set
    ........> $s=new Set<Int>();
    ........> for(my $i=0;$i<scalar(@{$p->VERTEX_SIZES});$i++){
    ........>     if($p->VERTEX_SIZES->[$i]==$p->DIM){
    ........>     $s+=$i;
    ........>     }
    ........> }
    > 





::

    polymake> ### iterate the set in two different ways
    ........> foreach(@{$s}){
    ........>     print $p->VERTICES->row($_)."\n";
    ........> }
    ........> foreach my $index(@{$s}){
    ........>     print $p->VERTICES->row($index)."\n";
    ........> }
    > 





::

    polymake> ### create a minor of the vertices matrix that only contains the simple ones
    ........> $special_points=$p->VERTICES->minor($s,All); print $special_points;

Writing scripts
~~~~~~~~~~~~~~~

Comprehensive information on how to use scripts within ``polymake`` can
be found `here <scripting/start>`__.
