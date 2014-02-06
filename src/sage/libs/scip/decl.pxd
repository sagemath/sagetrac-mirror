cdef extern from "scip/scip.h":
    ctypedef long SCIP_RETCODE
    cdef SCIP_RETCODE SCIP_OKAY

    ctypedef struct SCIP_VAR:
        char *name
        double obj
        int index

    ctypedef struct SCIP_QUADVARTERM:
        SCIP_VAR*  var
        double  lincoef
        double sqrcoef
        int nadjbilin
        int adjbilinsize
        int* adjbilin

    ctypedef struct SCIP_BILINTERM:
        SCIP_VAR*  var1
        SCIP_VAR*  var2
        double   coef

    ctypedef struct SCIP_CONS:
        pass

    ctypedef struct SCIP_SOL:
        pass

    ctypedef struct SCIP_c "SCIP":
       pass

    ctypedef long SCIP_OBJSENSE
    cdef SCIP_OBJSENSE SCIP_OBJSENSE_MAXIMIZE
    cdef SCIP_OBJSENSE SCIP_OBJSENSE_MINIMIZE

    ctypedef long SCIP_VARTYPE
    cdef SCIP_VARTYPE SCIP_VARTYPE_CONTINUOUS
    cdef SCIP_VARTYPE SCIP_VARTYPE_INTEGER
    cdef SCIP_VARTYPE SCIP_VARTYPE_BINARY

    cdef SCIP_RETCODE SCIPcreate(SCIP_c** scip)
    cdef SCIP_RETCODE SCIPfree(SCIP_c** scip)

    cdef SCIP_RETCODE SCIPincludeDefaultPlugins(SCIP_c *scip)
    cdef SCIP_RETCODE SCIPsetObjsense(SCIP_c *scip, SCIP_OBJSENSE)

    cdef SCIP_RETCODE SCIPcreateProb(SCIP_c *scip, char *name, void *, void *, void *, void *, void *, void *, void *)

    cdef char *SCIPgetProbName(SCIP_c *)
    cdef SCIP_RETCODE SCIPsetProbName(SCIP_c *, char *name)

    cdef SCIP_OBJSENSE SCIPgetObjsense(SCIP_c *)
    cdef long SCIPgetNVars(SCIP_c *)
    cdef long SCIPgetNConss(SCIP_c *)
    cdef SCIP_CONS ** SCIPgetConss(SCIP_c *)

    cdef SCIP_RETCODE SCIPwriteOrigProblem(SCIP_c *scip,
                                           char *filename,
                                           char *extension,
                                           bint genericnames)


    cdef double SCIPinfinity(SCIP_c *scip)

    # Parameters

    cdef enum SCIP_PARAMTYPE:
        SCIP_PARAMTYPE_BOOL
        SCIP_PARAMTYPE_INT
        SCIP_PARAMTYPE_LONGINT
        SCIP_PARAMTYPE_REAL
        SCIP_PARAMTYPE_CHAR
        SCIP_PARAMTYPE_STRING

    cdef SCIP_RETCODE SCIPreadParams(SCIP_c *, char *)
    cdef SCIP_RETCODE SCIPgetBoolParam(SCIP_c *scip, char *name, bint *value)
    cdef SCIP_RETCODE SCIPgetIntParam(SCIP_c *scip, char *name, int *value)
    cdef SCIP_RETCODE SCIPgetLongintParam(SCIP_c *scip, char *name, long *value)
    cdef SCIP_RETCODE SCIPgetRealParam(SCIP_c *scip, char *name, double *value)
    cdef SCIP_RETCODE SCIPgetCharParam(SCIP_c *scip, char *name, char *value)
    cdef SCIP_RETCODE SCIPgetStringParam(SCIP_c *scip, char *name, char **value)
    cdef SCIP_RETCODE SCIPsetBoolParam(SCIP_c *scip, char *name, bint value)
    cdef SCIP_RETCODE SCIPsetIntParam(SCIP_c *scip, char *name, int value)
    cdef SCIP_RETCODE SCIPsetLongintParam(SCIP_c *scip, char *name, long value)
    cdef SCIP_RETCODE SCIPsetRealParam(SCIP_c *scip, char *name, double value)
    cdef SCIP_RETCODE SCIPsetCharParam(SCIP_c *scip, char *name, char value)
    cdef SCIP_RETCODE SCIPsetStringParam(SCIP_c *scip, char *name, char *value)
    cdef SCIP_RETCODE SCIPreadParams(SCIP_c *scip, char *filename)
    cdef SCIP_RETCODE SCIPwriteParams(SCIP_c *scip, char *filename, bint comments, bint onlychanged)
    cdef SCIP_RETCODE SCIPresetParams(SCIP_c *scip)
    cdef int SCIPgetNParams(SCIP_c *scip)

    # Variables

    cdef SCIP_RETCODE SCIPcreateVar(SCIP_c *scip, SCIP_VAR **var, char *name,
                                    double lb, double ub,  double obj, SCIP_VARTYPE vartype,
                                    int initial, int removable,
                                    void *, void *, void *, void *, void *)

    cdef SCIP_RETCODE SCIPaddVar(SCIP_c *scip, SCIP_VAR *var)
    cdef SCIP_RETCODE SCIPreleaseVar(SCIP_c *scip, SCIP_VAR **var)
    cdef SCIP_VARTYPE SCIPvarGetType(SCIP_VAR *var)

    cdef bint SCIPvarIsTransformed(SCIP_VAR *var)
    cdef double SCIPvarGetLbGlobal(SCIP_VAR *var)
    cdef double SCIPvarGetUbGlobal(SCIP_VAR *var)
    cdef double SCIPvarGetLbOriginal(SCIP_VAR *var)
    cdef double SCIPvarGetUbOriginal(SCIP_VAR *var)

    cdef SCIP_RETCODE SCIPchgVarLb(SCIP_c *scip, SCIP_VAR *var, double val)
    cdef SCIP_RETCODE SCIPchgVarUb(SCIP_c *scip, SCIP_VAR *var, double val)

    cdef SCIP_RETCODE SCIPtransformVar(SCIP_c *scip,
                                       SCIP_VAR * 	var,
                                       SCIP_VAR ** 	transvar
                                       )

    cdef SCIP_RETCODE SCIPchgVarType(SCIP_c *scip,
                                     SCIP_VAR *var,
                                     SCIP_VARTYPE vartype,
                                     bint *infeasible)

    cdef SCIP_RETCODE SCIPchgVarObj(SCIP_c *scip,
                                    SCIP_VAR *var,
                                    double objval)

    # Constraints


    cdef SCIP_RETCODE SCIPreleaseCons(SCIP_c *scip, SCIP_CONS **cons)
    cdef SCIP_RETCODE SCIPaddCons(SCIP_c *scip, SCIP_CONS *cons)

    cdef char *SCIPconsGetName(SCIP_CONS *cons)

    cdef SCIP_RETCODE SCIPcreateConsAnd(SCIP_c *,
                                        SCIP_CONS **,
                                        char *,
                                        SCIP_VAR *res,
                                        int,
                                        SCIP_VAR **,
                                        bint 	initial,
                                        bint 	separate,
                                        bint 	enforce,
                                        bint 	check,
                                        bint 	propagate,
                                        bint 	local,
                                        bint 	modifiable,
                                        bint 	dynamic,
                                        bint 	removable,
                                        bint    stickingatnode
                                        )
    cdef int 	SCIPgetNVarsAnd (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR ** SCIPgetVarsAnd (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR *  SCIPgetResultantAnd (SCIP_c *scip, SCIP_CONS *cons)


    cdef SCIP_RETCODE SCIPcreateConsOr(SCIP_c *,
                                       SCIP_CONS **,
                                       char *,
                                       SCIP_VAR *res,
                                       int,
                                       SCIP_VAR **,
                                       bint 	initial,
                                       bint 	separate,
                                       bint 	enforce,
                                       bint 	check,
                                       bint 	propagate,
                                       bint 	local,
                                       bint 	modifiable,
                                       bint 	dynamic,
                                       bint 	removable,
                                       bint    stickingatnode
                                       )
    cdef int 	SCIPgetNVarsOr (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR ** SCIPgetVarsOr (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR *  SCIPgetResultantOr (SCIP_c *scip, SCIP_CONS *cons)

    cdef SCIP_RETCODE SCIPcreateConsXor(SCIP_c *,
                                        SCIP_CONS **,
                                        char *,
                                        int rhs,
                                        int,
                                        SCIP_VAR **,
                                        bint 	initial,
                                        bint 	separate,
                                        bint 	enforce,
                                        bint 	check,
                                        bint 	propagate,
                                        bint 	local,
                                        bint 	modifiable,
                                        bint 	dynamic,
                                        bint 	removable,
                                        bint    stickingatnode
                                        )

    cdef int SCIPgetNVarsXor (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR **SCIPgetVarsXor (SCIP_c *scip, SCIP_CONS *cons)
    cdef bint SCIPgetRhsXor (SCIP_c *scip, SCIP_CONS *cons)

    cdef SCIP_RETCODE SCIPcreateConsLinear(SCIP_c *,
                                           SCIP_CONS **,
                                           char *,
                                           int,
                                           SCIP_VAR **,
                                           double *,
                                           double,
                                           double,
                                           bint 	initial,
                                           bint 	separate,
                                           bint 	enforce,
                                           bint 	check,
                                           bint 	propagate,
                                           bint 	local,
                                           bint 	modifiable,
                                           bint 	dynamic,
                                           bint 	removable,
                                           bint         stickingatnode
                                           )

    cdef SCIP_RETCODE SCIPaddCoefLinear(SCIP_c *, SCIP_CONS *, SCIP_VAR *, double)
    cdef double SCIPgetLhsLinear(SCIP_c *scip, SCIP_CONS *cons)
    cdef double SCIPgetRhsLinear(SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_RETCODE SCIPchgLhsLinear(SCIP_c *scip, SCIP_CONS *cons, double lhs)
    cdef SCIP_RETCODE SCIPchgRhsLinear(SCIP_c *scip, SCIP_CONS *cons, double rhs)
    cdef int 	SCIPgetNVarsLinear(SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR **SCIPgetVarsLinear(SCIP_c *scip, SCIP_CONS *cons)
    cdef double *SCIPgetValsLinear(SCIP_c *scip, SCIP_CONS *cons)


    cdef SCIP_RETCODE SCIPcreateConsQuadratic(SCIP_c *,
                                              SCIP_CONS ** 	cons,
                                              char * 	name,
                                              int 	nlinvars,
                                              SCIP_VAR ** 	linvars,
                                              double * 	lincoefs,
                                              int 	nquadterms,
                                              SCIP_VAR ** 	quadvars1,
                                              SCIP_VAR ** 	quadvars2,
                                              double * 	quadcoefs,
                                              double 	lhs,
                                              double 	rhs,
                                              bint 	initial,
                                              bint 	separate,
                                              bint 	enforce,
                                              bint 	check,
                                              bint 	propagate,
                                              bint 	local,
                                              bint 	modifiable,
                                              bint 	dynamic,
                                              bint 	removable,
                                              )

    # Gets the number of variables in the linear part of a quadratic constraint.
    cdef int SCIPgetNLinearVarsQuadratic(SCIP_c* scip, SCIP_CONS* cons)

    # Gets the variables in the linear part of a quadratic constraint.
    # Length is given by SCIPgetNLinearVarsQuadratic.

    cdef SCIP_VAR** SCIPgetLinearVarsQuadratic(SCIP_c* scip, SCIP_CONS* cons)

    # Gets the coefficients in the linear part of a quadratic constraint.
    # Length is given by SCIPgetNQuadVarsQuadratic.

    cdef double* SCIPgetCoefsLinearVarsQuadratic(SCIP_c* scip, SCIP_CONS* cons)

    # Gets the number of quadratic variable terms of a quadratic constraint.
    cdef  int SCIPgetNQuadVarTermsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Gets the quadratic variable terms of a quadratic constraint.
    #  Length is given by SCIPgetNQuadVarTermsQuadratic.
    cdef  SCIP_QUADVARTERM* SCIPgetQuadVarTermsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Ensures that quadratic variable terms are sorted. */
    cdef SCIP_RETCODE SCIPsortQuadVarTermsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Finds the position of a quadratic variable term for a given variable.
    # @note If the quadratic variable terms have not been sorted before, then a
    # search may reorder the current order of the terms.
    cdef  SCIP_RETCODE SCIPfindQuadVarTermQuadratic(SCIP_c* scip, SCIP_CONS* cons, SCIP_VAR* var, int* pos)

    # Gets the number of bilinear terms of a quadratic constraint.
    cdef int SCIPgetNBilinTermsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Gets the bilinear terms of a quadratic constraint.
    # Length is given by SCIPgetNBilinTermQuadratic.
    cdef SCIP_BILINTERM* SCIPgetBilinTermsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Gets the left hand side of a quadratic constraint.
    cdef double SCIPgetLhsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Gets the right hand side of a quadratic constraint.
    cdef double SCIPgetRhsQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Check the quadratic function of a quadratic constraint for its semi-definiteness, if not done yet.
    cdef  SCIP_RETCODE SCIPcheckCurvatureQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Indicates whether the quadratic function of a quadratic constraint is (known to be) convex.
    cdef  bint SCIPisConvexQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Indicates whether the quadratic function of a quadratic constraint is (known to be) concave.
    cdef bint SCIPisConcaveQuadratic(SCIP_c *scip, SCIP_CONS* cons)

    # Gets the violation of a constraint by a solution. */
    cdef  SCIP_RETCODE SCIPgetViolationQuadratic(SCIP_c* scip, SCIP_CONS* cons, SCIP_SOL* sol, double* violation)

    # Indicates whether the quadratic constraint is local w.r.t. the current
    # local bounds.  That is, checks whether each variable with a square term is
    # fixed and for each bilinear term at least one variable is fixed.
    cdef  bint SCIPisLinearLocalQuadratic(SCIP_c *scip, SCIP_CONS* cons)
   
    cdef SCIP_RETCODE SCIPcreateConsLogicor(SCIP_c * 	scip,
                                            SCIP_CONS ** 	cons,
                                            char * 	name,
                                            int 	nvars,
                                            SCIP_VAR ** 	vars,
                                            bint 	initial,
                                            bint 	separate,
                                            bint 	enforce,
                                            bint 	check,
                                            bint 	propagate,
                                            bint 	local,
                                            bint 	modifiable,
                                            bint 	dynamic,
                                            bint 	removable,
                                            bint 	stickingatnode
                                            )

    cdef int SCIPgetNVarsLogicor(SCIP_c* scip, SCIP_CONS* cons)
    cdef SCIP_VAR** SCIPgetVarsLogicor(SCIP_c* scip, SCIP_CONS* cons)

    cdef SCIP_RETCODE 	SCIPcreateConsSetpart(SCIP_c *scip,
                                              SCIP_CONS **cons,
                                              char *name,
                                              int nvars,
                                              SCIP_VAR **vars,
                                              bint initial,
                                              bint separate,
                                              bint enforce,
                                              bint check,
                                              bint propagate,
                                              bint local,
                                              bint modifiable,
                                              bint dynamic,
                                              bint removable,
                                              bint stickingatnode)

    cdef SCIP_RETCODE 	SCIPcreateConsSetpack(SCIP_c *scip,
                                              SCIP_CONS **cons,
                                              char *name,
                                              int nvars,
                                              SCIP_VAR **vars,
                                              bint initial,
                                              bint separate,
                                              bint enforce,
                                              bint check,
                                              bint propagate,
                                              bint local,
                                              bint modifiable,
                                              bint dynamic,
                                              bint removable,
                                              bint stickingatnode)


    cdef SCIP_RETCODE 	SCIPcreateConsSetcover(SCIP_c *scip,
                                              SCIP_CONS **cons,
                                              char *name,
                                              int nvars,
                                              SCIP_VAR **vars,
                                              bint initial,
                                              bint separate,
                                              bint enforce,
                                              bint check,
                                              bint propagate,
                                              bint local,
                                              bint modifiable,
                                              bint dynamic,
                                              bint removable,
                                              bint stickingatnode)

    cdef int SCIPgetNVarsSetppc (SCIP_c *scip, SCIP_CONS *cons)
    cdef SCIP_VAR **SCIPgetVarsSetppc (SCIP_c *scip, SCIP_CONS *cons)

    # Solutions

    cdef SCIP_RETCODE SCIPsolve(SCIP_c *)
    cdef SCIP_RETCODE SCIPtransformProb(SCIP_c *)
    cdef bint SCIPisTransformed(SCIP_c *)


    cdef SCIP_RETCODE SCIPfreeTransform(SCIP_c *scip)

    cdef int SCIPgetNSols(SCIP_c *)
    cdef SCIP_SOL * SCIPgetBestSol(SCIP_c *)
    cdef SCIP_SOL ** SCIPgetSols(SCIP_c *)
    cdef double SCIPgetSolVal(SCIP_c *, SCIP_SOL *, SCIP_VAR *)
    cdef double SCIPgetSolOrigObj(SCIP_c *, SCIP_SOL *)

    ctypedef struct SCIP_HEUR:
        pass

    cdef SCIP_RETCODE SCIPcreateSol(SCIP_c *scip, SCIP_SOL **sol, SCIP_HEUR *heur)

    cdef SCIP_RETCODE SCIPsetSolVal(SCIP_c * 	scip,
                                    SCIP_SOL * 	sol,
                                    SCIP_VAR * 	var,
                                    double 	val)

    cdef SCIP_RETCODE SCIPtrySol(SCIP_c * 	scip,
                                 SCIP_SOL * 	sol,
                                 bint 	printreason,
                                 bint 	checkbounds,
                                 bint 	checkintegrality,
                                 bint 	checklprows,
                                 bint * 	stored
                                 )

    cdef SCIP_RETCODE SCIPfreeSol(SCIP_c *scip, SCIP_SOL **sol)

    # Message Handler

    cdef struct SCIP_MESSAGEHDLR:
        pass

    cdef SCIP_RETCODE SCIPsetMessagehdlr(SCIP_c *, SCIP_MESSAGEHDLR *)
    cdef SCIP_MESSAGEHDLR *SCIPgetMessagehdlr(SCIP_c *)
    cdef void  SCIPsetMessagehdlrQuiet(SCIP_c *, bint)

    # Conflicts

    ctypedef void SCIP_BDCHGIDX
    ctypedef int SCIP_BOUNDTYPE

    cdef SCIP_RETCODE 	SCIPinitConflictAnalysis (SCIP_c *scip)
    cdef SCIP_RETCODE 	SCIPaddConflictLb (SCIP_c *scip, SCIP_VAR *var, SCIP_BDCHGIDX *bdchgidx)
    cdef SCIP_RETCODE 	SCIPaddConflictUb (SCIP_c *scip, SCIP_VAR *var, SCIP_BDCHGIDX *bdchgidx)
    cdef SCIP_RETCODE 	SCIPaddConflictBd (SCIP_c *scip, SCIP_VAR *var, SCIP_BOUNDTYPE boundtype, SCIP_BDCHGIDX *bdchgidx)
    cdef SCIP_RETCODE 	SCIPaddConflictBinvar (SCIP_c *scip, SCIP_VAR *var)
    cdef SCIP_RETCODE 	SCIPanalyzeConflict (SCIP_c *scip, int validdepth, bint *success)
    cdef SCIP_RETCODE 	SCIPanalyzeConflictCons (SCIP_c *scip, SCIP_CONS *cons, bint *success)

cdef extern from "scip/cons_setppc.h":
    cdef enum SCIP_SetppcType:
        SCIP_SETPPCTYPE_PARTITIONING = 0
        SCIP_SETPPCTYPE_PACKING = 1
        SCIP_SETPPCTYPE_COVERING = 2

cdef extern from "scip/struct_paramset.h":
    cdef struct SCIP_Param:
        char *name
        SCIP_PARAMTYPE paramtype

cdef extern from "scip/scip.h":
    cdef SCIP_Param **SCIPgetParams(SCIP_c *scip)
    cdef SCIP_SetppcType  SCIPgetTypeSetppc (SCIP_c *scip, SCIP_CONS *cons)


# We must include this otherwise all kinds of weird things happen,
# e.g. the lhs and rhs of linear constraints are all garbled when
# printed in __repr__ (malb)

cdef extern from "scip/scipdefplugins.h":
    pass
