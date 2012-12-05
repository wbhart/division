#include <stdio.h>
#include <stdlib.h>
#include "mpir.h"
#include "gmp-impl.h"
#include "division.h"

#define NONE 0
#define MPIR 1
#define THIS 2

#define TIME THIS

#define TEST (TIME == 0)

int main(void)
{
   mp_limb_t r1, r2, dinv1, dinv2, a1[80], a2[80], b[40], q1[80], q2[80];
   mp_size_t limbs = 10;
   long i, j;

   for (i = 0; i < 1000; i++)
   {
      mpn_random(a1, 2*limbs);
      mpn_random(b, limbs);
      mpn_copyi(a2, a1, 2*limbs);
      b[limbs - 1] |= GMP_NUMB_HIGHBIT; /* normalise b */
      
#if (TIME == MPIR) || TEST
      invert_1(dinv1, b[limbs - 1], b[limbs - 2]);
#endif

#if (TIME == THIS) || TEST
      dinv2 = div_preinv1(b[limbs - 1], b[limbs - 2]);
#endif

#if TEST
      r1 = mpn_sb_div_qr(q1, a1, 2*limbs, b, limbs, dinv1);
      r2 = div_basecase(q2, a2, 2*limbs, b, limbs, dinv2);
#endif

#if TIME
      for (j = 0; j < 1000; j++)
      {
#if TIME == MPIR
          r1 = mpn_sb_div_qr(q1, a1, 2*limbs, b, limbs, dinv1);
          mpn_copyi(a1, a2, 2*limbs);
#endif

#if TIME == THIS
          r2 = div_basecase(q2, a2, 2*limbs, b, limbs, dinv2);
          mpn_copyi(a2, a1, 2*limbs);
#endif
      }
#endif

#if TEST
      if (r1 != r2)
      {
         printf("Error in most significant limb\n", j);
         printf("%lu vs %lu\n", r1, r2);
         abort();
      }
      
      for (j = 0; j < limbs; j++)
      {
         if (q1[limbs - j - 1] != q2[limbs - j - 1])
         {
            printf("Error in limb %ld of quotient\n", limbs - j - 1);
            printf("%lu vs %lu\n", q1[limbs - j - 1], q2[limbs - j - 1]);
            abort();
         }
      }

      for (j = 0; j < limbs; j++)
      {
         if (a1[limbs - j - 1] != a2[limbs - j - 1])
         {
            printf("Error in limb %ld of remainder\n", limbs - j - 1);
            printf("%lu vs %lu\n", a1[limbs - j - 1], a2[limbs - j - 1]);
            abort();
         }
      }
#endif

   }

   printf("PASS\n");

   return 0;
}

