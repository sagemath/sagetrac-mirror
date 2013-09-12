#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#ifdef USE_GMP
#include <gmp.h>

// newer gmp does not define the __* variant.
// flint/longlong.h still uses it
#ifndef __GMP_BITS_PER_MP_LIMB
# define __GMP_BITS_PER_MP_LIMB GMP_LIMB_BITS
#endif

#else
#include <mpir.h>
#endif
