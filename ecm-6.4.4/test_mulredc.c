#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gmp.h>

#include "mulredc.h"

void mp_print(const mp_limb_t *x, const int N) 
{
  int i;
  for (i = 0; i < N; ++i)
    {
      if (i>0) 
        printf (" + ");
      printf("%lu", x[i]);
      if (i>0) 
        printf ("*2^%d", i*GMP_NUMB_BITS); 
    }
  printf("\n");
}

static mp_limb_t
call_mulredc (const int N, mp_limb_t *z, const mp_limb_t *x, const mp_limb_t *y, 
              const mp_limb_t *m, const mp_limb_t invm)
{
  mp_limb_t cy;

  switch (N) 
    {
     case 1:
      cy = mulredc1(z, x[0], y[0], m[0], invm);
      break;
     case 2:
      cy = mulredc2(z, x, y, m, invm);
      break;
     case 3:
      cy = mulredc3(z, x, y, m, invm);
      break;
     case 4:
      cy = mulredc4(z, x, y, m, invm);
      break;
     case 5:
      cy = mulredc5(z, x, y, m, invm);
      break;
     case 6:
      cy = mulredc6(z, x, y, m, invm);
      break;
     case 7:
      cy = mulredc7(z, x, y, m, invm);
      break;
     case 8:
      cy = mulredc8(z, x, y, m, invm);
      break;
     case 9:
      cy = mulredc9(z, x, y, m, invm);
      break;
     case 10:
      cy = mulredc10(z, x, y, m, invm);
      break;
     case 11:
      cy = mulredc11(z, x, y, m, invm);
      break;
     case 12:
      cy = mulredc12(z, x, y, m, invm);
      break;
     case 13:
      cy = mulredc13(z, x, y, m, invm);
      break;
     case 14:
      cy = mulredc14(z, x, y, m, invm);
      break;
     case 15:
      cy = mulredc15(z, x, y, m, invm);
      break;
     case 16:
      cy = mulredc16(z, x, y, m, invm);
      break;
     case 17:
      cy = mulredc17(z, x, y, m, invm);
      break;
     case 18:
      cy = mulredc18(z, x, y, m, invm);
      break;
     case 19:
      cy = mulredc19(z, x, y, m, invm);
      break;
     case 20:
      cy = mulredc20(z, x, y, m, invm);
      break;
     default:
      cy = mulredc20(z, x, y, m, invm);
    }
  return cy;
}

#if defined(HAVE_NATIVE_MULREDC1_N)
static mp_limb_t
call_mulredc1 (const int N, mp_limb_t *z, const mp_limb_t x, const mp_limb_t *y, 
               const mp_limb_t *m, const mp_limb_t invm)
{
  mp_limb_t cy;

  switch (N) 
    {
     case 1:
      cy = mulredc1(z, x, y[0], m[0], invm);
      break;
     case 2:
      cy = mulredc1_2(z, x, y, m, invm);
      break;
     case 3:
      cy = mulredc1_3(z, x, y, m, invm);
      break;
     case 4:
      cy = mulredc1_4(z, x, y, m, invm);
      break;
     case 5:
      cy = mulredc1_5(z, x, y, m, invm);
      break;
     case 6:
      cy = mulredc1_6(z, x, y, m, invm);
      break;
     case 7:
      cy = mulredc1_7(z, x, y, m, invm);
      break;
     case 8:
      cy = mulredc1_8(z, x, y, m, invm);
      break;
     case 9:
      cy = mulredc1_9(z, x, y, m, invm);
      break;
     case 10:
      cy = mulredc1_10(z, x, y, m, invm);
      break;
     case 11:
      cy = mulredc1_11(z, x, y, m, invm);
      break;
     case 12:
      cy = mulredc1_12(z, x, y, m, invm);
      break;
     case 13:
      cy = mulredc1_13(z, x, y, m, invm);
      break;
     case 14:
      cy = mulredc1_14(z, x, y, m, invm);
      break;
     case 15:
      cy = mulredc1_15(z, x, y, m, invm);
      break;
     case 16:
      cy = mulredc1_16(z, x, y, m, invm);
      break;
     case 17:
      cy = mulredc1_17(z, x, y, m, invm);
      break;
     case 18:
      cy = mulredc1_18(z, x, y, m, invm);
      break;
     case 19:
      cy = mulredc1_19(z, x, y, m, invm);
      break;
     case 20:
      cy = mulredc1_20(z, x, y, m, invm);
      break;
     default:
      cy = mulredc1_20(z, x, y, m, invm);
    }
  return cy;
}
#endif

void test(mp_size_t N, int k)
{
  mp_limb_t *x, *y, *yp, *z, *m, invm, cy, cy2, *tmp, *tmp2, *tmp3;
  int i, j;
  
  x = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  y = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  z = (mp_limb_t *) malloc((N+1)*sizeof(mp_limb_t));
  m = (mp_limb_t *) malloc(N*sizeof(mp_limb_t));
  tmp = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));
  tmp2 = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));
  tmp3 = (mp_limb_t *) malloc((2*N+2)*sizeof(mp_limb_t));

  if (x == NULL || y == NULL || z == NULL || m == NULL || tmp == NULL ||
      tmp2 == NULL || tmp3 == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in test_mulredc\n");
      exit (1);
    }
 
  mpn_random2(m, N);
  m[0] |= 1UL;
  if (m[N-1] == 0) 
    m[N-1] = 1UL;

  invm = 1UL;
  for (i = 0; i < 10; ++i)
    invm = (2*invm-m[0]*invm*invm);
  invm = -invm;

  assert( (invm*m[0] +1UL) == 0UL);
  
  yp = y;
  for (i=0; i < k; ++i) {
    /* Try a few special cases */
    if (i == 0)
    {
      /* Try all 0, product should be 0 */
      for (j = 0; j < N; j++)
        x[j] = y[j] = 0;
    }
    else if (i == 1)
    {
      /* Try all 1 */
      for (j = 0; j < N; j++)
        x[j] = y[j] = 1;
    }
    else if (i == 2)
    {
      /* Try all 2^wordsize - 1 */
      for (j = 0; j < N; j++)
        x[j] = y[j] = ~(0UL);
    } 
    else 
    {
      /* In the other cases, try random data */
      if (i % 2 == 0)
        {
          /* Try squaring */
          mpn_random2(x, N);
          yp = x;
        }
      else
        {
          /* Try multiplication */
          mpn_random2(x, N);
          mpn_random2(y, N);
        }
    }
    
    /* Mixed mul and redc */
    cy = call_mulredc (N, z, x, yp, m, invm);
    
    if (cy)
      printf("!");
    z[N] = cy;
    /* Check with pure gmp : multiply by 2^(N*GMP_NUMB_BITS) and compare. */
    for (j=0; j < N; ++j) {
      tmp[j] = 0;
      tmp[j+N] = z[j]; 
    }
    tmp[2*N] = z[N];
    mpn_tdiv_qr(tmp2, tmp3, 0, tmp, 2*N+1, m, N);
    for (j=0; j < N; ++j)
      z[j] = tmp3[j]; 

    mpn_mul_n(tmp, x, yp, N);
    mpn_tdiv_qr(tmp2, tmp3, 0, tmp, 2*N, m, N);
    
    assert(mpn_cmp(z, tmp3, N) == 0);

#if defined(HAVE_NATIVE_MULREDC1_N)
    /* Test mulredc1_n() */
    z[N] = call_mulredc1 (N, z, x[0], yp, m, invm);
    tmp[0] = 0;
    for (j=0; j <= N; ++j) /* Multiply by 2^GMP_NUMB_BITS */
      tmp[j+1] = z[j];
    mpn_tdiv_qr(tmp2, tmp3, 0, tmp, N+2, m, N);
    for (j=0; j < N; ++j)
      z[j] = tmp3[j]; 

    tmp[N] = mpn_mul_1 (tmp, yp, N, x[0]);
    mpn_tdiv_qr(tmp2, tmp3, 0, tmp, N+1, m, N);
    
    assert(mpn_cmp(z, tmp3, N) == 0);
#endif
  }
  
  free(tmp); free(tmp2); free(tmp3);
  free(x); free(y); free(z); free(m);
}
  


int main(int argc, char** argv)
{
  int i, len;

  if (argc > 1) /* Test a specific length */
  {
    len = atoi (argv[1]);
    for (i = 0; i < 1; i++)
      test (len, 1000000);
    return 0;
  }

  for (;;) {
    for (i = 1; i <= 20; ++i) {
      test(i, 1000);
    }
#if 0
    test(1, 1000);
    test(2, 1000);
    test(3, 1000);
    test(4, 1000);
    test(5, 1000);
    test(6, 1000);
    test(7, 1000);
    test(8, 1000);
    test(9, 1000);
    test(10, 1000);
    test(11, 1000);
    test(12, 1000);
    test(13, 100);
    test(14, 100);
    test(15, 100);
    test(16, 100);
    test(17, 100);
    test(18, 100);
    test(44, 10);
    test(45, 10);
    test(46, 10);
    test(47, 10);
    test(48, 10);
    test(49, 10);
#endif
    printf("."); fflush(stdout);
  }
#if 0  
  x[0] = 12580274668139321508UL;
  x[1] = 9205793975152560417UL;
  x[2] = 7857372727033793057UL;

  y[0] = 13688385828267279103UL;
  y[1] = 10575011835742767258UL;
  y[2] = 8802048318027595690UL;

  
  m[0] = 2981542467342508025UL;
  m[1] = 5964669706257742025UL;
  m[2] = 18446744073678090270UL;

  invm = 9419286575570128311UL;

  carry = mulredc(z, x, y, m, 3, invm);

  printf("%lu + 2^64*(%lu + 2^64*%lu), carry=%lu\n", z[0], z[1], z[2], carry);
#endif
  return 0;
}


#if 0

W := 2^64;

x0:= 12580274668139321508;
x1:= 9205793975152560417;
x2:= 7857372727033793057;
x := x0 + W*(x1 + W*x2);

y0:= 13688385828267279103;
y1:= 10575011835742767258;
y2:= 8802048318027595690;
y := y0 + W*(y1 + W*y2);
  
m0:= 2981542467342508025;
m1:= 5964669706257742025;
m2:= 18446744073678090270;
m := m0 + W*(m1 + W*m2);
  
invm := 9419286575570128311;



#endif
