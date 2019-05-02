#ifndef SAGE_FLINT_NTL_WRAP_H
#define SAGE_FLINT_NTL_WRAP_H
/*
 * Similar to flint_wrap.h but specifically for wrapping the flint-NTL
 * interface.  It is separate from flint_wrap.h so that not every module
 * which uses flint has to pull in NTL headers as well.
 */

#include "flint_wrap.h"
#define ulong mp_limb_t
#define slong mp_limb_signed_t

#include <flint/NTL-interface.h>

#undef ulong
#undef slong
#endif
