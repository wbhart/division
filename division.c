#include <stdio.h>
#include "mpir.h"
#include "gmp-impl.h"
#include "division.h"

mp_limb_t div_preinv1(mp_limb_t d1, mp_limb_t d2)
{
   mp_limb_t q, r[2], p[2], cy;
   
   if (d2 + 1 == 0 && d1 + 1 == 0)
      return 0;

   if (d1 + 1 == 0)
      q = ~d1, r[1] = ~d2;
   else
      udiv_qrnnd(q, r[1], ~d1, ~d2, d1 + 1);

   r[0] = 0;

   if (d2 + 1 == 0)
      add_ssaaaa(cy, r[1], 0, r[1], 0, q);   
   else
   {
      umul_ppmm(p[1], p[0], q, ~d2 - 1);
      cy = mpn_add_n(r, r, p, 2);
   }
 
   p[0] = d2 + 1, p[1] = d1 + (d2 + 1 == 0);
   if (cy || mpn_cmp(r, p, 2) >= 0)
      q++;
   
   return q;
}

/* 
   Divide { a, m } by { b, n }, returning the high limb of the quotient
   (which will either be 0 or 1), storing the remainder in-place in
   { a, n } and the rest of the quotient in { q, m - n }.
   We require the most significant bit of { a, m } to be 1.
   dinv must be computed from b[n - 1], b[n - 2] by div_preinv1.
   Thus, currently we require n >= 2 and m > n.
*/
mp_limb_t div_basecase(mp_ptr q, mp_ptr a, mp_size_t m, 
                                 mp_srcptr b, mp_size_t n, mp_limb_t dinv)
{
   mp_limb_t ret;
   mp_size_t i;

   /* ensure { a + i, n } < { b, n } */
   if (ret = (mpn_cmp(a + m - n, b, n) >= 0))
      mpn_sub_n(a + m - n, a + m - n, b, n);
   
   for (i = m - 1; i >= n; i--)
   {
      divrem21_preinv(q[i - n], a[i], a[i - 1], dinv);
      a[i] -= mpn_submul_1(a + i - n, b, n, q[i - n]);

      if (mpn_cmp(a + i - n, b, n) >= 0 || a[i] != 0)
      {
         q[i - n]++;
         a[i] -= mpn_sub_n(a + i - n, a + i - n, b, n);
      }
   }

   return ret;
}
