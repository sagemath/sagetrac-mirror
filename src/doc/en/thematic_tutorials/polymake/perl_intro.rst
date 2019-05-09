.. -*- coding: utf-8 -*-
.. escape-backslashes
.. default-role:: math


Using Perl within polymake
==========================

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
structures, scalars(\ ``$``), arrays(\ ``@``), and hashes(\ ``%``). The
user always has to specify the type of a variable using the appropriate
symbol ``$``, ``@``, or ``%``. If you forget to do so, you will receive
the following error message:

::

   > i=5;
   polymake:  ERROR: Unquoted string "i" may clash with future reserved word.

Here are some simple commands illustrating how to use the different data
structures: ##### Scalars


::

    polymake> $i=5;
    ........> $j=6;
    ........> $sum=$i+$j; print $sum;
    11




Arrays
''''''


::

    polymake> @array=("a","b","c"); print scalar(@array);
    ........> push(@array,"d"); print "@array"; 
    ........> $first_entry=$array[0]; print $first_entry;
    ........> print join("\n",@array);
    ........> @array2=(3,1,4,2);
    ........> print sort(@array2);
    3a b c daa
    b
    c
    d1234




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
    zerofourzero, four0, 4




``polymake``-Perl
~~~~~~~~~~~~~~~~~

In addition to the three standard data structures, the enriched version
of ``Perl`` used in ``polymake`` also provides special data structures
for dealing with more complicated concepts. For an introduction to the
polymake object model see `here <.properties#objects>`__.

``polymake``\ ’s object hierarchy is completely reflected on the Perl
side. Let us create a small polytope as an example object.


::

    polymake> $p = new Polytope(POINTS=>[[1,0,1],[1,0,-1],[1,1,0],[1,-1,0]]);

Note that the ``Perl``-type of the variable ``$p`` is ``Scalar``, as the
variable is internally treated as a reference to a ``C++``-object. The
true nature of the object becomes visible if it is printed:


::

    polymake> print $p;
    Polymake::polytope::Polytope__Rational=ARRAY(0x55f7a179be78)




In this case it is a ``polymake`` object from the application
``polytope``, and it happens to be of type ``Polytope<Rational>``.
Technically, ``$p`` is a reference to an array (but it should be never
treated as an array unless you are deliberately trying to crash
``polymake``). If you want less technical information on the type of
your object, use this:


::

    polymake> print $p->type->full_name;
    Polytope<Rational>




“Small objects”: Data structures inherited from C++
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use objects that are inherited from the ``C++``-side of
``polymake`` in the interactive shell. A complete list of so-called
“small objects” can be found in the `online
documentation <https://polymake.org/release_docs/latest/common.html>`__
under the heading “Property types”. Here is a selection of three
different structures that facilitate everyday work with ``polymake``:
##### Arrays

The small object ``Array`` can be initialized in different ways and with
different template parameters:


::

    polymake> @array=("a","b","c");
    ........> $arr1=new Array<String>(\@array); print $arr1;
    ........> $arr2=new Array<Int>([3,2,5]); print $arr2;
    ........> $arr3=new Array<Int>(0,1,2,3); print $arr3;
    ........> $arr4=new Array<Int>(0..4); print $arr4;
    ........> $arr5=new Array<Int>($arr4); print $arr5;
    a b c3 2 50 1 2 30 1 2 3 40 1 2 3 4




You have random access:


::

    polymake> $arr5->[0] = 100;
    ........> print $arr5;
    100 1 2 3 4




It is also possible to convert the ``C++``-object ``Array`` into a
``Perl``-array by writing


::

    polymake> @arr4=@{$arr4}; print $arr2;
    3 2 5




or simply


::

    polymake> @arr4=@$arr4;

Sets
''''

On ``C++``-side sets are stored in a balanced binary search (AVL) tree.
For more information see the
`PTL-documentation <https://polymake.org/release_docs/master/PTL/classpm_1_1Set.html>`__.
In many cases, the small objects can be converted into ``Perl``-types in
the expected way:


::

    polymake> $set=new Set<Int>(3,2,5); print $set;
    ........> print $set->size;
    ........> @array_from_set=@$set;
    {2 3 5}3




Matrices
''''''''

Here is a simple way to initialize a matrix:


::

    polymake> $mat=new Matrix<Rational>([[2,1,4,0,0],[3,1,5,2,1],[1,0,4,0,6]]);
    ........> print $mat;
    2 1 4 0 0
    3 1 5 2 1
    1 0 4 0 6





You could also define it by passing a reference to an (``Perl``-)array
of ``Vectors``. The single entries are interpreted as different rows:


::

    polymake> $row1=new Vector<Rational>([2,1,4,0,0]);
    ........> $row2=new Vector<Rational>([3,1,5,2,1]);
    ........> $row3=new Vector<Rational>([1,0,4,0,6]);
    ........> @matrix_rows=($row1,$row2,$row3);
    ........> $matrix_from_array=new Matrix<Rational>(\@matrix_rows);

You can change a single entry of a matrix in the following way (if it is
not already assigned to an immutable property like ``VERTICES``!):


::

    polymake> $mat->row(1)->[1]=7;
    ........> print $mat->row(1)->[1], "\n";
    ........> print $mat, "\n";
    ........> $mat->elem(1,2)=8;
    ........> print $mat;
    7
    2 1 4 0 0
    3 7 5 2 1
    1 0 4 0 6
    
    2 1 4 0 0
    3 7 8 2 1
    1 0 4 0 6





A unit matrix of a certain dimension can be defined via the
user-function ``unit_matrix<COORDINATE_TYPE>(.)``:


::

    polymake> $unit_mat=4*unit_matrix<Rational>(3);
    ........> print $unit_mat;
    (3) (0 4)
    (3) (1 4)
    (3) (2 4)





The reason for the “strange output” is the implementation as *sparse
matrix*:


::

    polymake> print ref($unit_mat);
    Polymake::common::SparseMatrix_A_Rational_I_NonSymmetric_Z




However, some functions cannot deal with this special type of matrix. In
this case it is necessary to transform the sparse matrix into a dense
matrix first via:


::

    polymake> $dense=new Matrix<Rational>($unit_mat);print $dense;
    4 0 0
    0 4 0
    0 0 4





or just


::

    polymake> $dense2=dense($unit_mat);print $dense2;
    4 0 0
    0 4 0
    0 0 4





You can also work with matrices that have different types of coordinates
like ``Rational``, ``Float``, or ``Int``:


::

    polymake> $m_rat=new Matrix<Rational>(3/5*unit_matrix<Rational>(5)); print $m_rat, "\n"; 
    ........> $m2=$mat/$m_rat; print $m2, "\n";
    ........> $m_int=new Matrix<Int>(unit_matrix<Rational>(5)); print $m_int, "\n";
    3/5 0 0 0 0
    0 3/5 0 0 0
    0 0 3/5 0 0
    0 0 0 3/5 0
    0 0 0 0 3/5
    
    2 1 4 0 0
    3 7 8 2 1
    1 0 4 0 6
    3/5 0 0 0 0
    0 3/5 0 0 0
    0 0 3/5 0 0
    0 0 0 3/5 0
    0 0 0 0 3/5
    
    1 0 0 0 0
    0 1 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 1
    





Sometimes there is incompatible types:

::

   > $m3=$m_rat/$m_int;

C++/perl Interface module compilation failed; most likely due to a type
mismatch. Set the variable $Polymake::User::Verbose::cpp to a positive
value and repeat for more details.

The error message indicates that you need to convert the integer matrix
to a rational matrix first:


::

    polymake> $m3=$m_rat/(convert_to<Rational>($m_int)); print $m3;
    3/5 0 0 0 0
    0 3/5 0 0 0
    0 0 3/5 0 0
    0 0 0 3/5 0
    0 0 0 0 3/5
    1 0 0 0 0
    0 1 0 0 0
    0 0 1 0 0
    0 0 0 1 0
    0 0 0 0 1





By “/” you can add rows to a matrix, whereas “\|” adds columns. By the
way, this also works for ``Vector``.


::

    polymake> $z_vec=zero_vector<Int>($m_int->rows);
    ........> $extended_matrix=($z_vec|$m_int); print $extended_matrix;
    0 1 0 0 0 0
    0 0 1 0 0 0
    0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1





It is also possible to nest template parameters in any way you like,
e.g.


::

    polymake> $set=new Set<Int>(3,2,5);
    ........> $template_Ex=new Array<Set<Int>>((new Set<Int>(5,2,6)),$set); print $template_Ex; print ref($template_Ex);
    {2 5 6}
    {2 3 5}
    Polymake::common::Array__Set__Int




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
    4




Scalar properties can be used in arithmetic expressions right away.


::

    polymake> $i = ($p->N_FACETS * $p->N_FACETS) * 15;




::

    polymake> print $i;
    960




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
    ........> 
    ........> 
    ........> ### create a polytope from the matrix
    ........> $p=new Polytope<Rational>(POINTS=>$matrix);
    ........> print $p->FACETS;
    ........> print $p->DIM;
    ........> print $p->VERTEX_SIZES;
    ........> 
    ........> 
    ........> ### print "simple" vertices
    ........> for(my $i=0;$i<scalar(@{$p->VERTEX_SIZES});$i++){
    ........>     if($p->VERTEX_SIZES->[$i]==$p->DIM){
    ........>     print $i.": ".$p->VERTICES->row($i)."\n";
    ........>     }
    ........> }
    ........> 
    ........> 
    ........> ### put their indices in a set
    ........> $s=new Set<Int>();
    ........> for(my $i=0;$i<scalar(@{$p->VERTEX_SIZES});$i++){
    ........>     if($p->VERTEX_SIZES->[$i]==$p->DIM){
    ........>     $s+=$i;
    ........>     }
    ........> }
    ........> 
    ........> 
    ........> ### iterate the set in two different ways
    ........> foreach(@{$s}){
    ........>     print $p->VERTICES->row($_)."\n";
    ........> }
    ........> foreach my $index(@{$s}){
    ........>     print $p->VERTICES->row($index)."\n";
    ........> }
    ........> 
    ........> 
    ........> ### create a minor of the vertices matrix that only contains the simple ones
    ........> $special_points=$p->VERTICES->minor($s,All); print $special_points;
    -1




Writing scripts
~~~~~~~~~~~~~~~

Comprehensive information on how to use scripts within ``polymake`` can
be found `here <:user_guide:howto:scripting>`__.
