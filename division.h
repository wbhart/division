#include "mpir.h"
#include "gmp-impl.h"
#include "longlong.h"

#ifndef DIVISION_H
#define DIVISION_H

mp_limb_t div_preinv1(mp_limb_t d, mp_limb_t d2);

mp_limb_t div_basecase(mp_ptr q, mp_ptr a, mp_size_t m, 
                                 mp_srcptr b, mp_size_t n, mp_limb_t dinv);

#define divrem21_preinv(q, a_hi, a_lo, dinv) \
   do { \
      mp_limb_t __q2, __q3, __q4; \
      umul_ppmm((q), __q2, (a_hi), (dinv)); \
      umul_ppmm(__q3, __q4, (a_lo), (dinv)); \
      add_ssaaaa((q), __q2, (q), __q2, 0, __q3); \
      add_ssaaaa((q), __q2, (q), __q2, (a_hi), (a_lo)); \
   } while (0)

#endif
