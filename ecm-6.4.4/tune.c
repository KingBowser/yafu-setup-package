/* Tune program for GMP-ECM.

Copyright 2003, 2005, 2006, 2007, 2008, 2009, 2010, 2012 Paul Zimmermann,
Alexander Kruppa, Dave Newman and Jason Papadopoulos.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, see
http://www.gnu.org/licenses/ or write to the Free Software Foundation, Inc.,
51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA. */

#include <stdio.h>
#include <stdlib.h>
#include "ecm-gmp.h"
#include "ecm-impl.h"

/* 250ms, we (probably) don't need any more precision */
#define GRANULARITY 250
#define MAX_LOG2_LEN 18 /* 2 * 131072 */
#define MAX_LEN (1U << max_log2_len)
#define MAX_LOG2_MPZSPV_NORMALISE_STRIDE (MIN (12, max_log2_len))
#define M_str "95209938255048826235189575712705128366296557149606415206280987204268594538412191641776798249266895999715600261737863698825644292938050707507901970225804581"

#define ELAPSED elltime (__st, cputime () )
#define TUNE_FUNC_START(x)                   \
double x (size_t n)                          \
  { unsigned int __i, __k = 1; long __st;

/* Keep doubling the number of iterations until the timing is 
   at least GRANULARITY */
#define TUNE_FUNC_LOOP(x)                    \
  do {                                       \
    do {                                     \
      __st = cputime ();                     \
      for (__i = 0; __i < __k; __i++) { x; } \
      __k *= 2;                              \
    } while (ELAPSED < GRANULARITY);         \
    __k /= 2;                                \
    __st = ELAPSED;                          \
  } while (0)

#define TUNE_FUNC_END(x)                     \
  if (tune_verbose)                          \
    fprintf (stderr, #x "(%2ld) = %f\n", (long)n, (double) __k / (double) __st); \
  return (double) __k / (double) __st; }


/* Throughout, each function pointer points to a function
 * 
 *   double f0 (size_t n);
 *
 * that runs for at least GRANULARITY ms and then returns the number of
 * iterations performed per ms.
 *
 * X_Y_THRESHOLD denotes the threshold at which to start using Y for X. */


mpz_t M; /* yes, global variables */
gmp_randstate_t gmp_randstate;
size_t mp_size;
mpzspm_t mpzspm;
mpzv_t x, y, z, t;
spm_t spm;
spv_t spv;
mpzspv_t mpzspv;
int tune_verbose;
int max_log2_len = MAX_LOG2_LEN;
int min_log2_len = 3;

size_t MPZMOD_THRESHOLD;
size_t REDC_THRESHOLD;
size_t NTT_GFP_TWIDDLE_DIF_BREAKOVER = MAX_LOG2_LEN;
size_t NTT_GFP_TWIDDLE_DIT_BREAKOVER = MAX_LOG2_LEN;
size_t MUL_NTT_THRESHOLD;
size_t PREREVERTDIVISION_NTT_THRESHOLD;
size_t POLYINVERT_NTT_THRESHOLD;
size_t POLYEVALT_NTT_THRESHOLD;
size_t MPZSPV_NORMALISE_STRIDE = 256;

void
mpz_quick_random (mpz_t x, mpz_t M, unsigned long b)
{
  mpz_urandomb (x, gmp_randstate, b);
  if (mpz_cmp (x, M) >= 0)
    mpz_sub (x, x, M);
}


double
tune_mpres_mul (mp_size_t limbs, int repr)
{
  mpmod_t modulus;
  mpres_t x, y, z;
  mpz_t N, p, q;
  unsigned int __k = 1, __i;
  long __st;

  mpz_init (N);
  mpz_init (p);
  mpz_init (q);
  
  /* No need to generate a probable prime, just ensure N is not
     divisible by 2 or 3 */
  do
    {
      mpz_urandomb (N, gmp_randstate, limbs * GMP_NUMB_BITS);
      while (mpz_gcd_ui (NULL, N, 6) != 1)
        mpz_add_ui (N, N, 1);
    }
  while ((mp_size_t) mpz_size (N) != limbs);
  
  if (repr == ECM_MOD_MPZ)
    mpmod_init_MPZ (modulus, N);
  else if (repr == ECM_MOD_MODMULN)
    mpmod_init_MODMULN (modulus, N);
  else if (repr == ECM_MOD_REDC)
    mpmod_init_REDC (modulus, N);

  mpz_urandomm (p, gmp_randstate, N);
  mpz_urandomm (q, gmp_randstate, N);
  
  mpres_init (x, modulus);
  mpres_init (y, modulus);
  mpres_init (z, modulus);

  mpres_set_z (x, p, modulus);
  mpres_set_z (y, q, modulus);

  TUNE_FUNC_LOOP (mpres_mul (z, x, y, modulus));

  mpres_clear (x, modulus);
  mpres_clear (y, modulus);
  mpres_clear (z, modulus);
  mpmod_clear (modulus);
  mpz_clear (N);
  mpz_clear (p);
  mpz_clear (q);

  return (double) __k / (double) __st;
}

double
tune_mpres_sqr (mp_size_t limbs, int repr)
{
  mpmod_t modulus;
  mpres_t x, z;
  mpz_t N, p;
  unsigned int __k = 1, __i;
  long __st;

  mpz_init (N);
  mpz_init (p);
  
  /* No need to generate a probable prime, just ensure N is not
     divisible by 2 or 3 */
  do
    {
      mpz_urandomb (N, gmp_randstate, limbs * GMP_NUMB_BITS);
      while (mpz_gcd_ui (NULL, N, 6) != 1)
        mpz_add_ui (N, N, 1);
    }
  while ((mp_size_t) mpz_size (N) != limbs);
  
  if (repr == ECM_MOD_MPZ)
    mpmod_init_MPZ (modulus, N);
  else if (repr == ECM_MOD_MODMULN)
    mpmod_init_MODMULN (modulus, N);
  else if (repr == ECM_MOD_REDC)
    mpmod_init_REDC (modulus, N);

  mpz_urandomm (p, gmp_randstate, N);
  
  mpres_init (x, modulus);
  mpres_init (z, modulus);

  mpres_set_z (x, p, modulus);

  TUNE_FUNC_LOOP (mpres_sqr (z, x, modulus));

  mpres_clear (x, modulus);
  mpres_clear (z, modulus);
  mpmod_clear (modulus);
  mpz_clear (N);
  mpz_clear (p);

  return (double) __k / (double) __st;
}

double
tune_mpres_mul_mpz (size_t n)
{
  return tune_mpres_mul (n, ECM_MOD_MPZ);
}

double
tune_mpres_mul_modmuln (size_t n)
{
  return tune_mpres_mul (n, ECM_MOD_MODMULN);
}

double
tune_mpres_mul_redc (size_t n)
{
  return tune_mpres_mul (n, ECM_MOD_REDC);
}

TUNE_FUNC_START (tune_spv_ntt_gfp_dif)
  NTT_GFP_TWIDDLE_DIF_BREAKOVER = n;
  TUNE_FUNC_LOOP (spv_ntt_gfp_dif (spv, max_log2_len, spm));
TUNE_FUNC_END (tune_spv_ntt_gfp_dif)


TUNE_FUNC_START (tune_spv_ntt_gfp_dit)
  NTT_GFP_TWIDDLE_DIT_BREAKOVER = n;
  TUNE_FUNC_LOOP (spv_ntt_gfp_dit (spv, max_log2_len, spm));
TUNE_FUNC_END (tune_spv_ntt_gfp_dit_recursive)


TUNE_FUNC_START (tune_ntt_mul)
  MUL_NTT_THRESHOLD = 0;

  TUNE_FUNC_LOOP (ntt_mul (z, x, y, 1 << n, NULL, 1, mpzspm));
TUNE_FUNC_END (tune_ntt_mul)


TUNE_FUNC_START (tune_list_mul)
  TUNE_FUNC_LOOP (list_mul (z, x, 1 << n, 1, y, 1 << n, 1, t));
TUNE_FUNC_END (tune_list_mul)


TUNE_FUNC_START (tune_ntt_PrerevertDivision)
  PREREVERTDIVISION_NTT_THRESHOLD = 0;

  TUNE_FUNC_LOOP (ntt_PrerevertDivision (z, x, y, mpzspv, mpzspv,
    1 << n, t, mpzspm));
TUNE_FUNC_END (tune_ntt_PrerevertDivision)


TUNE_FUNC_START (tune_PrerevertDivision)
  TUNE_FUNC_LOOP (PrerevertDivision (z, x, y, 1 << n, t, mpzspm->modulus));
TUNE_FUNC_END (tune_PrerevertDivision)


TUNE_FUNC_START (tune_ntt_PolyInvert)
  POLYINVERT_NTT_THRESHOLD = 1 << n;
  
  TUNE_FUNC_LOOP (ntt_PolyInvert (z, x, 1 << n, t, mpzspm));
TUNE_FUNC_END (tune_ntt_PolyInvert)


TUNE_FUNC_START (tune_PolyInvert)
  
  TUNE_FUNC_LOOP (PolyInvert (z, x, 1 << n, t, mpzspm->modulus));
TUNE_FUNC_END (tune_PolyInvert)
  

TUNE_FUNC_START (tune_ntt_polyevalT)
  unsigned int i;
  mpzv_t *Tree = (mpzv_t *) malloc ((n + 1) * sizeof (mpzv_t));
  if (Tree == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in tune_ntt_polyevalT\n");
      exit (1);
    }
  
  for (i = 0; i <= n; i++)
    Tree[i] = x;

  POLYEVALT_NTT_THRESHOLD = 1 << n;

  TUNE_FUNC_LOOP (ntt_polyevalT (z, 1 << n, Tree, t, mpzspv, mpzspm, NULL));

  free (Tree);
TUNE_FUNC_END (tune_ntt_polyevalT) 


TUNE_FUNC_START (tune_polyevalT)
  unsigned int i;
  mpzv_t *Tree = (mpzv_t *) malloc ((n + 1) * sizeof (mpzv_t));
  if (Tree == NULL)
    {
      fprintf (stderr, "Cannot allocate memory in tune_polyevalT\n");
      exit (1);
    }

  for (i = 0; i <= n; i++)
    Tree[i] = x;

  TUNE_FUNC_LOOP (polyeval_tellegen (z, 1 << n, Tree, t, 3 * (1 << n),
	  x, mpzspm->modulus, NULL));

  free (Tree);
TUNE_FUNC_END (tune_polyevalT)


TUNE_FUNC_START (tune_mpzspv_normalise)
  MPZSPV_NORMALISE_STRIDE = 1 << n;
  
  TUNE_FUNC_LOOP (mpzspv_normalise (mpzspv, 0,
    1 << MAX_LOG2_MPZSPV_NORMALISE_STRIDE, mpzspm));
TUNE_FUNC_END (tune_mpzspv_normalise)


TUNE_FUNC_START (tune_ecm_mul_lo_n)
  mp_limb_t rp[2 * MPN_MUL_LO_THRESHOLD];
  mp_limb_t xp[MPN_MUL_LO_THRESHOLD];
  mp_limb_t yp[MPN_MUL_LO_THRESHOLD];

  if (n > 1 && n < (mp_size + 1) / 2)
    return 0.0;
  
  mpn_random (xp, mp_size);
  mpn_random (yp, mp_size);
  
  mpn_mul_lo_threshold[mp_size] = n;

  TUNE_FUNC_LOOP (ecm_mul_lo_n (rp, xp, yp, mp_size));
TUNE_FUNC_END (tune_ecm_mul_lo_n)

/* Return the lowest n with min_n <= n < max_n such that
 * f1(t) >= f0(t) for all t in [n, n + k), or return max_n if no such
 * n exists. This function will typically return high values if there
 * is no 'clean' threshold between f0(n) and f1(n). */
size_t
crossover2 (double (*f0)(size_t), double (*f1)(size_t),
    size_t min_n, size_t max_n, size_t k)
{
  size_t n = min_n;
  size_t t;
  
  while (n < max_n)
    {
      for (t = MIN (max_n, n + k); t > n; t--)
        {
	  if ((f0)(t - 1) > (f1)(t - 1))
            break;
	}

      if (t == n)
        return n;

      n = t;
    };

  return max_n;
}


/* Assume f0 and f1 are monotone decreasing. Return the first n in the range
 * [min_n, max_n) for which f1(n) >= f0(n), or return max_n if no such n
 * exists. We use a bisection algorithm so the function is fast but
 * may give slightly varied results. */
size_t
crossover (double (*f0)(size_t), double (*f1)(size_t),
    size_t min_n, size_t max_n)
{
  size_t mid_n;
  
#ifdef TUNE_SLOW
  return crossover2 (f0, f1, min_n, max_n, 1);
#endif
    
  if (min_n == max_n)
    return min_n;

  mid_n = (max_n + min_n) / 2;
  return ((f0)(mid_n) > (f1)(mid_n))
    ? crossover (f0, f1, mid_n + 1, max_n)
    : crossover (f0, f1, min_n, mid_n);
}


/* Return the n in the range [min_n, max_n) that maximises f(n).
 * We make no assumptions about the shape of f(n) and so evaluate
 * f at every point. */
size_t
maximise (double (*f)(size_t), size_t min_n, size_t max_n)
{
  size_t n, best_n = 0;
  double f_n, f_best_n = -1.0;

  for (n = min_n; n < max_n; n++)
    {
      f_n = f (n);
      if (f_n > f_best_n)
        {
	  f_best_n = f_n;
	  best_n = n;
	}
    }

  return best_n;
}

/* Debugging. Print the value of f0(n) and f1(n) and which is fastest. */
void
print_timings (double (*f0)(size_t), double (*f1)(size_t),
  size_t min_n, size_t max_n)
{
  size_t n;
  double f0_n, f1_n;
  
  for (n = min_n; n < max_n; n++)
    {
      f0_n = (f0)(n);
      f1_n = (f1)(n);
      printf ("n=%2ld: %8.2f %8.2f (f%d)\n",
          (long) n, f0_n, f1_n, (f0_n <= f1_n) ? 1 : 0);
    }
}

int
main (int argc, char **argv)
{
  spv_size_t i;
  unsigned long b;

  while (argc > 1)
    {
      if (strcmp (argv[1], "-v") == 0)
        {
          tune_verbose = 1;
          argc --;
          argv ++;
        }
      else if (argc > 2 && strcmp (argv[1], "-max_log2_len") == 0)
        {
          max_log2_len = atoi (argv[2]);
	  if (max_log2_len < min_log2_len)
	    max_log2_len = min_log2_len;
          argc -= 2;
          argv += 2;
        }
      else
        {
          fprintf (stderr, "Usage: tune [-v] [-max_log2_len nnn]\n");
          exit (1);
        }
    }
  
  gmp_randinit_default (gmp_randstate);
  mpz_init_set_str (M, M_str, 10);
  b = (unsigned long) mpz_sizeinbase (M, 2);

  x = init_list (MAX_LEN);
  y = init_list (MAX_LEN);
  z = init_list (MAX_LEN);
  t = init_list (list_mul_mem (MAX_LEN / 2) + 3 * MAX_LEN / 2);
  
  mpzspm = mpzspm_init (MAX_LEN, M);
  if (mpzspm == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in mpzspm_init\n");
      exit (1);
    }
  mpzspv = mpzspv_init (MAX_LEN, mpzspm);
  if (mpzspv == NULL)
    {
      fprintf (stderr, "Error, cannot allocate memory in mpzspv_init\n");
      exit (1);
    }
  mpzspv_random (mpzspv, 0, MAX_LEN, mpzspm);
  
  for (i = 0; i < MAX_LEN; i++)
    mpz_quick_random (x[i], M, b);
  for (i = 0; i < MAX_LEN; i++)
    mpz_quick_random (y[i], M, b);
  for (i = 0; i < MAX_LEN; i++)
    mpz_quick_random (z[i], M, b);    
  
  spm = mpzspm->spm[0];
  spv = mpzspv[0];
  
  MPZMOD_THRESHOLD = crossover2 (tune_mpres_mul_modmuln, tune_mpres_mul_mpz,
      1, 512, 10);
  
  printf ("#define MPZMOD_THRESHOLD %lu\n", (unsigned long) MPZMOD_THRESHOLD);
  
  REDC_THRESHOLD = crossover2 (tune_mpres_mul_mpz, tune_mpres_mul_redc,
      MPZMOD_THRESHOLD, 512, 10);
  
  printf ("#define REDC_THRESHOLD %lu\n", (unsigned long) REDC_THRESHOLD);

  mpn_mul_lo_threshold[0] = 0;
  mpn_mul_lo_threshold[1] = 0;

  printf ("#define MPN_MUL_LO_THRESHOLD_TABLE {0, 0");

  for (mp_size = 2; mp_size < MPN_MUL_LO_THRESHOLD; mp_size++)
    {
      mpn_mul_lo_threshold[mp_size] = maximise (tune_ecm_mul_lo_n, 0, mp_size);
      printf (", %lu", (unsigned long) mpn_mul_lo_threshold[mp_size]);
      fflush (stdout);
    }

  printf ("}\n");
	  
  NTT_GFP_TWIDDLE_DIF_BREAKOVER = maximise
      (tune_spv_ntt_gfp_dif, min_log2_len, max_log2_len);

  printf ("#define NTT_GFP_TWIDDLE_DIF_BREAKOVER %lu\n",
      (unsigned long) NTT_GFP_TWIDDLE_DIF_BREAKOVER);
   
  NTT_GFP_TWIDDLE_DIT_BREAKOVER = maximise
      (tune_spv_ntt_gfp_dit, min_log2_len, max_log2_len);

  printf ("#define NTT_GFP_TWIDDLE_DIT_BREAKOVER %lu\n",
      (unsigned long) NTT_GFP_TWIDDLE_DIT_BREAKOVER);
  
  MUL_NTT_THRESHOLD = 1 << crossover2 (tune_list_mul, tune_ntt_mul, 1,
      max_log2_len, 2);

  printf ("#define MUL_NTT_THRESHOLD %lu\n", (unsigned long) MUL_NTT_THRESHOLD);

  PREREVERTDIVISION_NTT_THRESHOLD = 1 << crossover2 (tune_PrerevertDivision,
      tune_ntt_PrerevertDivision, 1, max_log2_len, 2);

  printf ("#define PREREVERTDIVISION_NTT_THRESHOLD %lu\n",
      (unsigned long) PREREVERTDIVISION_NTT_THRESHOLD);

  POLYINVERT_NTT_THRESHOLD = 1 << crossover (tune_PolyInvert,
      tune_ntt_PolyInvert, 5, max_log2_len);

  printf ("#define POLYINVERT_NTT_THRESHOLD %lu\n", 
      (unsigned long) POLYINVERT_NTT_THRESHOLD);
  
  POLYEVALT_NTT_THRESHOLD = 1 << crossover (tune_polyevalT,
      tune_ntt_polyevalT, 5, max_log2_len);

  printf ("#define POLYEVALT_NTT_THRESHOLD %lu\n", 
      (unsigned long) POLYEVALT_NTT_THRESHOLD);
  
  MPZSPV_NORMALISE_STRIDE = 1 << maximise (tune_mpzspv_normalise,
      1, MAX_LOG2_MPZSPV_NORMALISE_STRIDE);
	  
  printf ("#define MPZSPV_NORMALISE_STRIDE %lu\n", 
      (unsigned long) MPZSPV_NORMALISE_STRIDE);

  mpzspv_clear (mpzspv, mpzspm);
  mpzspm_clear (mpzspm);
  
  clear_list (x, MAX_LEN);
  clear_list (y, MAX_LEN);
  clear_list (z, MAX_LEN);
  clear_list (t, list_mul_mem (MAX_LEN / 2) + 3 * MAX_LEN / 2);

  mpz_clear (M);
  
  gmp_randclear (gmp_randstate);
  
  return 0;
}
