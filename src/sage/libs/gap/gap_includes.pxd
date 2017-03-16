###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


cdef extern from "<gap/system.h>":
    ctypedef char Char
    ctypedef int Int
    ctypedef unsigned char UChar

cdef extern from "<gap/libgap.h>":
    void libgap_initialize(int argc, char** argv)
    ctypedef void(*libgap_gasman_callback_ptr)()
    void libgap_set_gasman_callback(libgap_gasman_callback_ptr callback)
    ctypedef void(*libgap_error_func_ptr)(char* msg)
    void libgap_set_error_handler(libgap_error_func_ptr error_handler)
    void libgap_call_error_handler()
    void libgap_finalize()
    void libgap_start_interaction(char* inputline)
    char* libgap_get_output()
    char* libgap_get_error()
    void libgap_finish_interaction()
    void libgap_mark_stack_bottom()
    void libgap_enter()
    void libgap_exit()

cdef extern from "<gap/code.h>":
    ctypedef unsigned int Stat
    ctypedef Stat* PtrBody

cdef extern from "<gap/gap.h>":
    ctypedef unsigned int UInt
    ctypedef void* ExecStatus
    void ViewObjHandler(void*)
    void InitializeGap(int*, char** argv)
    void set_system_variables(char**, char**)
    cdef UInt Last
    cdef UInt Last2
    cdef UInt Last3
    cdef ExecStatus STATUS_END
    cdef ExecStatus STATUS_RETURN_VAL
    cdef ExecStatus STATUS_RETURN_VOID
    cdef ExecStatus STATUS_TNM
    cdef ExecStatus STATUS_QUIT
    cdef ExecStatus STATUS_EOF
    cdef ExecStatus STATUS_ERROR
    cdef ExecStatus STATUS_QQUIT

cdef extern from "<gap/objects.h>":
    ctypedef void* Obj
    Obj SHALLOW_COPY_OBJ(Obj obj)
    bint IS_INTOBJ(Obj obj)
    Obj INTOBJ_INT(Int)
    Int INT_INTOBJ(Obj)
    UInt TNUM_OBJ(Obj obj)
    char* TNAM_OBJ(Obj obj)
    cdef int FIRST_REAL_TNUM
    cdef int FIRST_CONSTANT_TNUM
    cdef int T_INT             "(FIRST_CONSTANT_TNUM+ 0)"
    cdef int T_INTPOS
    cdef int T_INTNEG
    cdef int T_RAT
    cdef int T_CYC
    cdef int T_FFE
    cdef int T_PERM2
    cdef int T_PERM4
    cdef int T_BOOL            "(FIRST_CONSTANT_TNUM+ 8)"
    cdef int T_CHAR            "(FIRST_CONSTANT_TNUM+ 9)"
    cdef int T_FUNCTION
    cdef int T_FLAGS
    cdef int T_MACFLOAT
    cdef int T_RESERVED_BY_GAP
    cdef int LAST_CONSTANT_TNUM
    cdef int IMMUTABLE
    cdef int FIRST_IMM_MUT_TNUM
    cdef int FIRST_RECORD_TNUM
    cdef int T_PREC
    cdef int LAST_RECORD_TNUM
    cdef int FIRST_LIST_TNUM
    cdef int FIRST_PLIST_TNUM
    cdef int T_PLIST
    cdef int T_PLIST_NDENSE
    cdef int T_PLIST_DENSE
    cdef int T_PLIST_DENSE_NHOM
    cdef int T_PLIST_DENSE_NHOM_SSORT
    cdef int T_PLIST_DENSE_NHOM_NSORT
    cdef int T_PLIST_EMPTY
    cdef int T_PLIST_HOM
    cdef int T_PLIST_HOM_NSORT
    cdef int T_PLIST_HOM_SSORT
    cdef int T_PLIST_TAB
    cdef int T_PLIST_TAB_NSORT
    cdef int T_PLIST_TAB_SSORT
    cdef int T_PLIST_TAB_RECT
    cdef int T_PLIST_TAB_RECT_NSORT
    cdef int T_PLIST_TAB_RECT_SSORT
    cdef int T_PLIST_CYC
    cdef int T_PLIST_CYC_NSORT
    cdef int T_PLIST_CYC_SSORT
    cdef int T_PLIST_FFE
    cdef int LAST_PLIST_TNUM
    cdef int T_RANGE_NSORT
    cdef int T_RANGE_SSORT
    cdef int T_BLIST
    cdef int T_BLIST_NSORT
    cdef int T_BLIST_SSORT
    cdef int T_STRING            "(FIRST_LIST_TNUM+50)"
    cdef int T_STRING_NSORT
    cdef int T_STRING_SSORT
    cdef int LAST_LIST_TNUM
    cdef int LAST_IMM_MUT_TNUM
    cdef int FIRST_EXTERNAL_TNUM
    cdef int T_COMOBJ
    cdef int T_POSOBJ
    cdef int T_DATOBJ
    cdef int T_WPOBJ
    cdef int LAST_EXTERNAL_TNUM
    cdef int LAST_REAL_TNUM
    cdef int LAST_VIRTUAL_TNUM
    cdef int FIRST_COPYING_TNUM
    cdef int COPYING
    cdef int LAST_COPYING_TNUM
    cdef int FIRST_TESTING_TNUM
    cdef int TESTING
    cdef int LAST_TESTING_TNUM

cdef extern from "<gap/read.h>":
    void* ReadEvalCommand(Obj context, UInt *dualSemicolon)
    void* ReadEvalFile()
    void* ReadEvalResult "TLS(ReadEvalResult)"
    bint READ_ERROR()

cdef extern from "<gap/scanner.h>":
    void ClearError()
    UInt NrError  "TLS(NrError)"
    UInt Symbol  "TLS(Symbol)"
    void GetSymbol()
    void Match (UInt symbol, char* msg, UInt skipto)
    int S_ILLEGAL
    int S_IDENT
    int S_UNBIND
    int S_ISBOUND
    int S_TRYNEXT
    int S_INFO
    int S_ASSERT
    int S_SAVEWS
    int S_LOADWS
    int S_LBRACK
    int S_LBRACE
    int S_BLBRACK
    int S_BLBRACE
    int S_RBRACK
    int S_RBRACE
    int S_DOT
    int S_BDOT
    int S_LPAREN
    int S_RPAREN
    int S_COMMA
    int S_DOTDOT
    int S_COLON
    int S_PARTIALINT
    int S_INT
    int S_TRUE
    int S_FALSE
    int S_CHAR
    int S_STRING
    int S_PARTIALSTRING
    int S_REC
    int S_FUNCTION
    int S_LOCAL
    int S_END
    int S_MAPTO
    int S_MULT
    int S_DIV
    int S_MOD
    int S_POW
    int S_PLUS
    int S_MINUS
    int S_EQ
    int S_LT
    int S_GT
    int S_NE
    int S_LE
    int S_GE
    int S_IN
    int S_NOT
    int S_AND
    int S_OR
    int S_ASSIGN
    int S_IF
    int S_FOR
    int S_WHILE
    int S_REPEAT
    int S_THEN
    int S_ELIF
    int S_ELSE
    int S_FI
    int S_DO
    int S_OD
    int S_UNTIL
    int S_BREAK
    int S_RETURN
    int S_QUIT
    int S_QQUIT
    int S_CONTINUE
    int S_SEMICOLON
    int S_EOF

cdef extern from "<gap/gvars.h>":
    UInt GVarName(char* name)
    void AssGVar(UInt gvar, Obj val)
    Obj VAL_GVAR(UInt gvar)

cdef extern from "<gap/string.h>":
    char* CSTR_STRING(Obj list)
    int GET_LEN_STRING(Obj list)
    bint IS_STRING(Obj obj)
    bint IsStringConv(Obj obj)
    bint ConvString(Obj obj)
    void C_NEW_STRING(Obj new_gap_string, int length, char* c_string)

cdef extern from "<gap/gasman.h>":
    void InitGlobalBag(Obj* addr, char* cookie)
    Obj NewBag(UInt type, UInt size)
    void CHANGED_BAG(Obj bag)
    void MARK_BAG(Obj bag)
    bint IS_MARKED_ALIVE(Obj bag)
    bint IS_MARKED_DEAD(Obj bag)
    bint IS_MARKED_HALFDEAD(Obj bag)
    cdef UInt NrAllBags
    cdef UInt SizeAllBags
    cdef UInt NrLiveBags
    cdef UInt SizeLiveBags
    cdef UInt NrDeadBags
    cdef UInt SizeDeadBags
    cdef UInt NrHalfDeadBags
    UInt CollectBags(UInt size, UInt full)
    void CallbackForAllBags(void (*func)(Obj))
    char* TNAM_BAG(Obj obj)
    UInt TNUM_BAG(Obj)
    UInt SIZE_BAG(Obj)
    void CheckMasterPointers()
    Obj* MptrBags
    Obj* YoungBags
    Obj* OldBags
    Obj* AllocBags
    Obj* MarkedBags
    Obj* ChangedBags

# in gasman.c but not declared in gasman.h
cdef extern Obj* StopBags
cdef extern Obj* EndBags

cdef extern from "<gap/ariths.h>":
    Obj SUM (Obj, Obj)
    Obj DIFF(Obj, Obj)
    Obj PROD(Obj, Obj)
    Obj QUO(Obj, Obj)
    Obj POW(Obj, Obj)
    Obj MOD(Obj, Obj)
    Obj CALL_0ARGS(Obj f)              # 0 arguments
    Obj CALL_1ARGS(Obj f, Obj a1)      # 1 argument
    Obj CALL_2ARGS(Obj f, Obj a1, Obj a2)
    Obj CALL_3ARGS(Obj f, Obj a1, Obj a2, Obj a3)
    Obj CALL_4ARGS(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4)
    Obj CALL_5ARGS(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5)
    Obj CALL_6ARGS(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5, Obj a6)
    Obj CALL_XARGS(Obj f, Obj args)   # more than 6 arguments
    bint EQ(Obj opL, Obj opR)
    bint LT(Obj opL, Obj opR)

cdef extern from "<gap/calls.h>":
    bint IS_FUNC(Obj)

cdef extern from "<gap/plist.h>":
    Obj NEW_PLIST(int type, int len)
    bint IS_PLIST(Obj lst)
    int LEN_PLIST(Obj lst)
    Obj ELM_PLIST(Obj lst, int pos)

cdef extern from "<gap/lists.h>":
    void UNB_LIST(Obj list, int pos)
    bint IS_LIST(Obj lst)
    int LEN_LIST(Obj lst)
    Obj ELM_LIST(Obj lst, int pos)

cdef extern from "<gap/listfunc.h>":
    void AddList(Obj list, Obj obj)
    void AddPlist(Obj list, Obj obj)

cdef extern from "<gap/records.h>":
    char* NAME_RNAM(UInt rnam)
    UInt RNamIntg(int i)
    bint IS_REC(Obj obj)
    Obj ELM_REC(Obj rec, UInt rnam)
    UInt RNamName(Char* name)

cdef extern from "<gap/precord.h>":
    Obj NEW_PREC(int len)
    int LEN_PREC(Obj rec)
    int GET_RNAM_PREC(Obj rec, int i)
    Obj GET_ELM_PREC(Obj rec, int i)
    void AssPRec(Obj rec, UInt rnam, Obj val)
    void UnbPRec(Obj rec, UInt rnam)
    bint IsbPRec(Obj rec, UInt rnam)
    Obj ElmPRec(Obj rec, UInt rnam)

cdef extern from "<gap/cyclotom.h>":
    pass

cdef extern from "<gap/bool.h>":
    cdef Obj True
    cdef Obj False

cdef extern from "<gap/vars.h>":
     cdef int T_LVARS
     Obj BottomLVars "TLS(BottomLVars)"





