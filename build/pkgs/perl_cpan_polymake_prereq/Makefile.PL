## The missing Makefile.PL of polymake.
## We use it to declare the dependencies.
require 5.006;
use ExtUtils::MakeMaker;
WriteMakefile
(
  'NAME' => 'Polymake',
  'VERSION' => '0',
  'PREREQ_PM' => {
      'Types::Standard' => 0,
          'Type::Library' => 0,
          'XML::LibXML' => 0,
          'XML::Writer' => 0,
          'XML::LibXSLT' => 0,
          'File::Slurp' => 0,
          'JSON' => 0,
          'SVG' => 0,
          'MongoDB' => 0
 }
)
;
