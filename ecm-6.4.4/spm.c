/* spm.c - "small prime modulus" functions to precompute an inverse and a
   primitive root for a small prime

Copyright 2005, 2006, 2008, 2009, 2010, 2012 Dave Newman, Jason Papadopoulos,
Paul Zimmermann, Alexander Kruppa.

The SP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The SP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the SP Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,
MA 02110-1301, USA. */

#include <stdlib.h>
#include "sp.h"

/* Returns the exponent of $q$ in the factorisation of $n$ */
static int
exponent (const sp_t q, sp_t n)
{
  int i;
  for (i = 0; n % q == (sp_t) 0; i++, n /= q);
  return i;
}

/* Returns i so that ord(a) = q^i. This assumes that ord(a) is indeed
   a low power of q. */
static int
ordpow (const sp_t q, sp_t a, const sp_t sp, const sp_t mul_c)
{
  int i = 0;
  for (i = 0; a != (sp_t) 1; i++, a = sp_pow (a, q, sp, mul_c));
  return i;
}

/* initialize roots of unity and twiddle factors for one NTT.
   If successful, returns 1.
   If unsuccessful, returns 0 (and frees allocated memory) */
static int
nttdata_init (const sp_t sp, const sp_t mul_c, 
		const sp_t prim_root, const spv_size_t log2_len,
		sp_nttdata_t data, spv_size_t breakover)
{
  spv_t r, t;
  spv_size_t i, j, k;

  r = data->ntt_roots = 
	  (spv_t) sp_aligned_malloc ((log2_len + 1) * sizeof(sp_t));
  if (r == NULL)
    return 0;

  i = log2_len;
  r[i] = prim_root;
  for (i--; (int)i >= 0; i--)
    r[i] = sp_sqr (r[i+1], sp, mul_c);

  k = MIN(log2_len, breakover);
  t = data->twiddle = (spv_t) sp_aligned_malloc (sizeof(sp_t) << k);
  if (t == NULL)
    {
      sp_aligned_free (r);
      return 0;
    }
  data->twiddle_size = 1 << k;

  for (i = k; i; i--) 
    {
      sp_t w = r[i];
      for (j = t[0] = 1; j < ((spv_size_t) 1 << (i-1)); j++) 
      	t[j] = sp_mul (t[j-1], w, sp, mul_c);

      t += j;
    }
  return 1;
}

static void
nttdata_clear(sp_nttdata_t data)
{
  sp_aligned_free(data->ntt_roots);
  sp_aligned_free(data->twiddle);
}

/* Compute some constants, including a primitive n'th root of unity. 
   Returns NULL in case of error.
   k is the number of limbs of the number to factor
*/
spm_t
spm_init (spv_size_t n, sp_t sp, mp_size_t k)
{
  sp_t a, b, bd, sc;
  spv_size_t q, nc, ntt_power;
  spm_t spm = (spm_t) malloc (sizeof (__spm_struct));
  if (spm == NULL)
    return NULL;

  ASSERT (sp % (sp_t) n == (sp_t) 1);

  spm->sp = sp;
  sp_reciprocal (spm->mul_c, sp);

  /* compute spm->invm = -1/p mod B where B = 2^GMP_NUMB_BITS */
  a = sp_pow (2, GMP_NUMB_BITS, sp, spm->mul_c); /* a = B mod p */
  a = sp_inv (a, sp, spm->mul_c);                /* a = 1/B mod p */
  /* a = 1/B mod p thus B*a - 1 = invm*p */
  a --;
  b = GMP_NUMB_MASK;
#if SP_NUMB_BITS == W_TYPE_SIZE - 2
  a = (a << 2) + (b >> (GMP_NUMB_BITS - 2));
  b = (b << 2) & GMP_NUMB_MASK;
  udiv_qrnnd (bd, sc, a, b, sp << 2);
#else
  a = (a << 1) + (b >> (GMP_NUMB_BITS - 1));
  b = (b << 1) & GMP_NUMB_MASK;
  udiv_qrnnd (bd, sc, a, b, sp << 1);
#endif
  spm->invm = bd;

  /* compute spm->Bpow = B^(k+1) mod p */
  spm->Bpow = sp_pow (2, GMP_NUMB_BITS * (k + 1), sp, spm->mul_c);

  /* find an $n$-th primitive root $a$ of unity $(mod sp)$. */

  /* Construct a $b$ whose order $(mod sp)$ is equal to $n$.
     We try different $a$ values and test if the exponent of $q$ in $ord(a)$
     is at least as large as in $n$. If it isn't, we move to another $a$.
     If it is, we optionally exponentiate to make the exponents equal and
     test for the remaining $q$'s.
     We assume that the largest prime dividing $n$ is very small, 
     so no optimizations in factoring n are made. */
  a = 2;
  b = a;
  nc = n; /* nc is remaining cofactor of n */
  q = 2;
  sc = sp - 1;
#ifdef PARI
  printf ("/* spm_init */ n = %lu; sp = %lu; /* PARI */\n", n, sp);
  printf ("exponent(a,b) = {local(i); while(b%%a == 0,i++;b/=a); "
	  "return(i)} /* PARI */\n");
#endif
  for ( ; nc != (spv_size_t) 1; q++)
    {
      if (nc % q == (spv_size_t) 0)
        {
	  const int k = exponent (q, n); /* q^k || n */
	  sp_t d;
	  int l;
#ifdef PARI
	  printf ("exponent(%lu, n) == %d /* PARI */\n", q, k);
#endif
          /* Remove all factors of $q$ from $sp-1$ */
          for (d = sp - 1; d % q == (spv_size_t) 0; d /= q);
	  bd = sp_pow (b, d, sp, spm->mul_c);
          /* Now ord(bd) = q^l, q^l || ord(a) */
	  l = ordpow (q, bd, sp, spm->mul_c);
#ifdef PARI
	  printf ("exponent(%lu, znorder(Mod(%lu, sp))) == %d /* PARI */\n", 
		  q, b, l);
#endif
          if (l < k)
            {
              /* No good, q appears in ord(a) in a lower power than in n. 
		 Try next $a$ */
              a++;
              b = a;
              nc = n;
              q = 1; /* Loop increment following "continue" will make q=2 */
	      sc = sp - 1;
              continue;
            }
          else
            {
	      /* Reduce the exponent of $q$ in $ord(b)$ until is it 
		 equal to that in $n$ */
	      for ( ; l > k; l--)
		{
#ifdef PARI
		  printf ("Exponentiating %lu by %lu\n", b, q);
#endif
		  b = sp_pow (b, q, sp, spm->mul_c);
		}
#ifdef PARI
	      printf ("New b = %lu\n", b);
#endif
            }
	  do {nc /= q;} while (nc % q == 0); /* Divide out all q from nc */
	  while (sc % q == (sp_t) 0) /* Divide out all q from sc */
	    sc /= q;
        }
    }
  
  b = sp_pow (b, sc, sp, spm->mul_c);
#ifdef PARI
  printf ("znorder(Mod(%lu, sp)) == n /* PARI */\n", b, sp, n);
#endif

  /* turn this into a primitive n'th root of unity mod p */
  spm->prim_root = b;
  spm->inv_prim_root = sp_inv (b, sp, spm->mul_c);

  /* initialize auxiliary data for all supported power-of-2 NTT sizes */
  ntt_power = 0;
  while (1)
    {
      if (n & (1 << ntt_power))
        break;
      ntt_power++;
    }

  if (!nttdata_init (sp, spm->mul_c, 
		sp_pow (spm->prim_root, 
			n >> ntt_power, sp, spm->mul_c),
		ntt_power, spm->nttdata, 
		NTT_GFP_TWIDDLE_DIF_BREAKOVER))
    goto free_spm;
  if (!nttdata_init (sp, spm->mul_c, 
		sp_pow (spm->inv_prim_root, 
			n >> ntt_power, sp, spm->mul_c),
		ntt_power, spm->inttdata, 
		NTT_GFP_TWIDDLE_DIT_BREAKOVER))
    goto free_nttdata;
  spm->scratch = (spv_t) sp_aligned_malloc (
		  	MAX_NTT_BLOCK_SIZE * sizeof(sp_t));
  if (spm->scratch == NULL)
    goto free_inttdata;
    
  return spm;

  free_inttdata:
  nttdata_clear (spm->inttdata);
  
  free_nttdata:
  nttdata_clear (spm->nttdata);

  free_spm:
  free (spm);

  return NULL;
}

void
spm_clear (spm_t spm)
{
  nttdata_clear (spm->nttdata);
  nttdata_clear (spm->inttdata);
  sp_aligned_free (spm->scratch);
  free (spm);
}



