## The missing Makefile.PL of polymake.
## We use it to declare the dependencies.
require 5.006;
use ExtUtils::MakeMaker;
WriteMakefile
(
  'NAME' => 'PolymakeStage1',
  'VERSION' => '0',
 'PREREQ_PM' => {
     'Moose' => 2.20
         # A MongoDB dependency causes Module::Runtime to be updated to >= 0.14,
         # which is incompatible with Moose from the perl 5.18 distribution on macOS;
         # this causes the Type-Tiny installation to fail.
         # http://blogs.perl.org/users/ovid/2014/02/moose-dying-with-invalid-version-format-version-required.html
 }
)
;
