SAGE_SPKG_CONFIGURE([brial], [
  SAGE_SPKG_DEPCHECK([boost m4ri], [
    # If we're using the system m4ri and boost, ensure that we can
    # compile and run an executable linked against both libbrial and
    # libbrial_groebner (both are used by SageMath).
    AC_LANG_PUSH(C++)
    SAVED_LIBS=$LIBS
    LIBS="$LIBS -lbrial -lbrial_groebner"
    AC_MSG_CHECKING([if we can link against brial libraries])
    AC_RUN_IFELSE([
      AC_LANG_PROGRAM([
        #include <polybori.h>
        #include <polybori/groebner/groebner_alg.h>
        USING_NAMESPACE_PBORI
        USING_NAMESPACE_PBORIGB

	class MyConstant : public BooleConstant{
          public: void negate() { this->m_value = !this->m_value; }
        };
      ],[
        BoolePolyRing r = BoolePolyRing(2, COrderEnums::dlex);
        ReductionStrategy rs = ReductionStrategy(r);
        rs.llReduceAll(); // uses groebner lib
        if (2 != r.nVariables()) { return 1; }
        if (r.constant(true) == r.constant(false)) { return 2; }
        MyConstant f = MyConstant();
        f.negate(); // ensures v1.1.0+ if m_value isn't const
        if (!f.isOne()) { return 3; }
        return 0;
      ])
    ],
    [
      AC_MSG_RESULT([yes])
      sage_spkg_install_brial=no
    ],
    [
      AC_MSG_RESULT([no])
      sage_spkg_install_brial=yes
    ])
    LIBS=$SAVED_LIBS
    AC_LANG_POP
  ],
  [ # If we're installing sage's boost or m4ri, then we have to
    # install its BRiAl, too.
    sage_spkg_install_brial=yes
  ])
])
