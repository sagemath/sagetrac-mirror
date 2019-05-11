SAGE_SPKG_CONFIGURE(
    [perl_cpan_polymake_prereq], [
    AX_PROG_PERL_MODULES(XML::Writer XML::LibXML XML::LibXSLT File::Slurp JSON SVG MongoDB,
      [],
      [sage_spkg_install_perl_cpan_polymake_prereq=yes])
])
