/*----------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Ben Buhrow. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

Some parts of the code (and also this header), included in this 
distribution have been reused from other sources. In particular I 
have benefitted greatly from the work of Jason Papadopoulos's msieve @ 
www.boo.net/~jasonp, Scott Contini's mpqs implementation, and Tom St. 
Denis Tom's Fast Math library.  Many thanks to their kind donation of 
code to the public domain.
       				   --bbuhrow@gmail.com 11/24/09
----------------------------------------------------------------------*/

#include "yafu.h"
#include "qs.h"
#include "factor.h"
#include "util.h"
#include "common.h"
#include "tdiv_macros_common.h"
#include "tdiv_macros_64k.h"

//#define SIQSDEBUG 1

/*
We are given an array of bytes that has been sieved.  The basic trial 
division strategy is as follows:

1) Scan through the array and 'mark' locations that meet criteria 
indicating they may factor completely over the factor base.  

2) 'Filter' the marked locations by trial dividing by small primes
that we did not sieve.  These primes are all less than 256.  If after
removing small primes the location does not meet another set of criteria,
remove it from the 'marked' list (do not subject it to further trial
division).

3) Divide out primes from the factor base between 256 and 2^13 or 2^14, 
depending on the version (2^13 for 32k version, 2^14 for 64k).  

4) Resieve primes between 2^{13|14} and 2^16, max.  

5) Primes larger than 2^16 will have been bucket sieved.  Remove these
by scanning the buckets for sieve hits equal to the current block location.

6) If applicable/appropriate, factor a remaining composite with squfof

this file contains code implementing 3)


*/

void tdiv_medprimes_64k(uint8 parity, uint32 poly_id, uint32 bnum, 
						 static_conf_t *sconf, dynamic_conf_t *dconf)
{
	//we have flagged this sieve offset as likely to produce a relation
	//nothing left to do now but check and see.
	int i;
	uint32 bound, tmp, prime, root1, root2, report_num;
	int smooth_num;
	uint32 *fb_offsets;
	sieve_fb *fb;
	sieve_fb_compressed *fbc;
	fb_element_siqs *fullfb_ptr, *fullfb = sconf->factor_base->list;
	uint32 block_loc;

#ifdef USE_8X_MOD_ASM
	uint16 *bl_sizes;
	uint16 *bl_locs;

	bl_sizes = (uint16 *)xmalloc_align(8 * sizeof(uint16));
	bl_locs = (uint16 *)xmalloc_align(8 * sizeof(uint16));

#endif

	fullfb_ptr = fullfb;
	if (parity)
	{
		fb = dconf->fb_sieve_n;
		fbc = dconf->comp_sieve_n;
	}
	else
	{
		fb = dconf->fb_sieve_p;
		fbc = dconf->comp_sieve_p;
	}

#ifdef QS_TIMING
	gettimeofday(&qs_timing_start, NULL);
#endif

	for (report_num = 0; report_num < dconf->num_reports; report_num++)
	{
#ifdef USE_YAFU_TDIV
		z32 *tmp32 = &dconf->Qvals32[report_num];
#endif

		if (!dconf->valid_Qs[report_num])
			continue;

		// for each report, we trial divide, and then either trial divide or
		// resieve.  for the first trial division step, we have either
		// unrolled C routines or SIMD assembly routines to choose from.  The
		// second trial division step only has a C routine - the more optimized
		// path is resieving.
		//
		// the basic idea of trial division is to test if the block location
		// in question lies on the arithmetic progression of a prime.  By the time
		// we get to this routine, the arithmetic progression has been reset for
		// the next sieve block, so we have to do a few manipulations to revert
		// to the "real" progression (add the blocksize back to the roots).  
		// there are various methods for doing the test depending on the size of
		// the primes (and therefore how many times it could have possibly landed
		// in the sieve block).  The most straightforward, and the method the SIMD
		// assembly uses, is to see if the roots (adjusted for the "real" progression)
		// minus the block location in question, divided by the prime, is zero.  if
		// so, this block location is on the arithmetic progression of that prime,
		// and we can proceed to trial divide the prime into Q(x) for this sieve hit.
		// 
		// the basic idea of resieving is to start from the end of the "real"
		// arithmetic progression, repeatedly subtract each prime, and test after
		// each subtraction if we've hit the sieve location in question.  if so,
		// we know this location is on the prime's arithmetic progression and can 
		// proceed to trial divide.  For "large enough" primes, this is very efficient
		// because we only need to do a few subtractions and tests instead of a 
		// division.  Since we are not really doing division tests, and instead are
		// doing multiplication by inverses, and futhermore since we might be doing those
		// multiplications 8 at a time using SIMD, resieving is only a win for the 
		// very largest primes less than 16 bits in size.  
		//
		// for OS/architecture/compilers where resieving isn't implemented, there are
		// further trial division steps instead.  These are more efficient than
		// the "check for exact division of a difference" method described above,
		// but are only implemented in portable C.  See code below for more detail.

		// pull the details of this report to get started.
		fb_offsets = &dconf->fb_offsets[report_num][0];
		smooth_num = dconf->smooth_num[report_num];
		block_loc = dconf->reports[report_num];

		//do the primes less than the blocksize.  primes bigger than the blocksize can be handled
		//even more efficiently.
		//a couple of observations from jasonp:
		//if a prime divides Q(x), then this index (j) and either
		//root1 or root2 are on the same arithmetic progression.  this we can
		//test with a single precision mod operation
		//do the first few until the rest can be done in batches of 4 that are aligned to 16 byte
		//boundaries.  this is necessary to use the SSE2 batch mod code, if it ever
		//becomes faster...

		i=sconf->sieve_small_fb_start;

		// the bound before resieving takes over
		bound = sconf->factor_base->fb_14bit_B;

		// although if we are using 8x ASM division,
		// use a lower bound for a while first
#ifdef USE_8X_MOD_ASM
		bound = sconf->factor_base->fb_10bit_B;
#endif
		
		// single-up test until i is a multiple of 8
		while ((uint32)i < bound && ((i & 7) != 0))
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			//tmp = distance from this sieve block offset to the end of the block
			tmp = 65536 - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - 65536;

			//tmp = advance the offset to where it should be after the interval, and
			//check to see if that's where either of the roots are now.  if so, then
			//this offset is on the progression of the sieve for this prime
			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#ifdef USE_8X_MOD_ASM

		// the 64k blocksize will overflow the 16 bit word
		// so we use blocksize-1 and then compensate in 
		// the assembly by adding back 1 when needed
		bl_sizes[0] = 65535;
		bl_sizes[1] = 65535;
		bl_sizes[2] = 65535;
		bl_sizes[3] = 65535;
		bl_sizes[4] = 65535;
		bl_sizes[5] = 65535;
		bl_sizes[6] = 65535;
		bl_sizes[7] = 65535;

		bl_locs[0] = block_loc;
		bl_locs[1] = block_loc;
		bl_locs[2] = block_loc;
		bl_locs[3] = block_loc;
		bl_locs[4] = block_loc;
		bl_locs[5] = block_loc;
		bl_locs[6] = block_loc;
		bl_locs[7] = block_loc;

		MOD_INIT_8X;

		while ((uint32)i < bound)
		{
			uint32 tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(8);
#else
			MOD_CMP_8X("8");
#endif

			//if ((((((uint32)fbc->root1[i] + 65536 - block_loc) % fbc->prime[i]) == 0) ||
			//	((((uint32)fbc->root2[i] + 65536 - block_loc) % fbc->prime[i]) == 0)) &&
			//	((tmp3 & 0x2) == 0))
			//	printf("index = %d, prime = %u, inv = %u, corr = %u, root1 = %u, root2 = %u, "
			//	"loc = %u, result = %u\n",
			//		i, fbc->prime[i], fullfb_ptr->small_inv[i], 
			//		fullfb_ptr->correction[i], fbc->root1[i], fbc->root2[i], block_loc, tmp3);

			if (tmp3 == 0)
			{
				i += 8;
				continue;
			}
			
			if (tmp3 & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (tmp3 & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (tmp3 & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (tmp3 & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (tmp3 & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (tmp3 & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (tmp3 & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (tmp3 & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;			

		}

		bound = sconf->factor_base->fb_12bit_B;
		while ((uint32)i < bound)
		{
			uint32 tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(10);
#else
			MOD_CMP_8X("10");
#endif

			if (tmp3 == 0)
			{
				i += 8;
				continue;
			}
			
			if (tmp3 & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (tmp3 & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (tmp3 & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (tmp3 & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (tmp3 & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (tmp3 & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (tmp3 & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (tmp3 & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;			

		}

		bound = sconf->factor_base->fb_14bit_B;
			
		while ((uint32)i < bound)
		{
			uint32 tmp3 = 0;

#ifdef _MSC_VER
			MOD_CMP_8X(12);
#else
			MOD_CMP_8X("12");
#endif

			if (tmp3 == 0)
			{
				i += 8;
				continue;
			}
			
			if (tmp3 & 0x2)
			{
				DIVIDE_RESIEVED_PRIME(0);
			}

			if (tmp3 & 0x8)
			{
				DIVIDE_RESIEVED_PRIME(1);
			}

			if (tmp3 & 0x20)
			{
				DIVIDE_RESIEVED_PRIME(2);
			}

			if (tmp3 & 0x80)
			{
				DIVIDE_RESIEVED_PRIME(3);
			}

			if (tmp3 & 0x200)
			{
				DIVIDE_RESIEVED_PRIME(4);
			}

			if (tmp3 & 0x800)
			{
				DIVIDE_RESIEVED_PRIME(5);
			}

			if (tmp3 & 0x2000)
			{
				DIVIDE_RESIEVED_PRIME(6);
			}

			if (tmp3 & 0x8000)
			{
				DIVIDE_RESIEVED_PRIME(7);
			}

			i += 8;			

		}

		TDIV_MED_CLEAN;

#else
		//now do things in batches of 4 which are aligned on 16 byte boundaries.
		while ((uint32)i < bound)
		{
			uint64 q64;
			uint32 tmp1 = 65536 - block_loc;
			uint32 tmp2 = 65536 - block_loc;
			uint32 tmp3 = 65536 - block_loc;
			uint32 tmp4 = 65536 - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - 65536;
			tmp2 = tmp2 - 65536;
			tmp3 = tmp3 - 65536;
			tmp4 = tmp4 - 65536;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}

		//now cleanup any that don't fit in the last batch
		while ((uint32)i < bound)
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			tmp = 65536 - block_loc;
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - 65536;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

		// single-up test until i is a multiple of 8
		while ((uint32)i < bound && ((i & 7) != 0))
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			//tmp = distance from this sieve block offset to the end of the block
			tmp = 65536 - block_loc;
	
			//tmp = tmp/prime + 1 = number of steps to get past the end of the sieve
			//block, which is the state of the sieve now.
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - 65536;

			//tmp = advance the offset to where it should be after the interval, and
			//check to see if that's where either of the roots are now.  if so, then
			//this offset is on the progression of the sieve for this prime
			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

		//now do things in batches of 4 which are aligned on 16 byte boundaries.
		while ((uint32)i < bound)
		{
			uint64 q64;
			uint32 tmp1 = 65536 - block_loc;
			uint32 tmp2 = 65536 - block_loc;
			uint32 tmp3 = 65536 - block_loc;
			uint32 tmp4 = 65536 - block_loc;

			tmp1 = tmp1 + fullfb_ptr->correction[i];
			q64 = (uint64)tmp1 * (uint64)fullfb_ptr->small_inv[i];
			tmp1 = q64 >> FOGSHIFT; 
			tmp1 = tmp1 + 1;
			tmp1 = block_loc + tmp1 * fullfb_ptr->prime[i];
			
			tmp2 = tmp2 + fullfb_ptr->correction[i+1];
			q64 = (uint64)tmp2 * (uint64)fullfb_ptr->small_inv[i+1];
			tmp2 = q64 >> FOGSHIFT; 
			tmp2 = tmp2 + 1;
			tmp2 = block_loc + tmp2 * fullfb_ptr->prime[i+1];

			tmp3 = tmp3 + fullfb_ptr->correction[i+2];
			q64 = (uint64)tmp3 * (uint64)fullfb_ptr->small_inv[i+2];
			tmp3 = q64 >> FOGSHIFT; 
			tmp3 = tmp3 + 1;
			tmp3 = block_loc + tmp3 * fullfb_ptr->prime[i+2];

			tmp4 = tmp4 + fullfb_ptr->correction[i+3];
			q64 = (uint64)tmp4 * (uint64)fullfb_ptr->small_inv[i+3];
			tmp4 = q64 >>  FOGSHIFT; 
			tmp4 = tmp4 + 1;
			tmp4 = block_loc + tmp4 * fullfb_ptr->prime[i+3];

			tmp1 = tmp1 - 65536;
			tmp2 = tmp2 - 65536;
			tmp3 = tmp3 - 65536;
			tmp4 = tmp4 - 65536;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp1 == root1 || tmp1 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp2 == root1 || tmp2 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp3 == root1 || tmp3 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			if (tmp4 == root1 || tmp4 == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;

		}

		//now cleanup any that don't fit in the last batch
		while ((uint32)i < bound)
		{
			prime = fbc->prime[i];
			root1 = fbc->root1[i];
			root2 = fbc->root2[i];

			tmp = 65536 - block_loc;
			tmp = 1+(uint32)(((uint64)(tmp + fullfb_ptr->correction[i])
					* (uint64)fullfb_ptr->small_inv[i]) >> FOGSHIFT_2); 
			tmp = block_loc + tmp*prime;
			tmp = tmp - 65536;

			if (tmp == root1 || tmp == root2)
			{
				//it will divide Q(x).  do so as many times as we can.
				DIVIDE_ONE_PRIME;
			}
			i++;
		}

#endif

		// either after 8x SSE2 ASM, or standard trial division, record
		// how many factors we've found so far
		dconf->smooth_num[report_num] = smooth_num;	

	}

#ifdef QS_TIMING
	gettimeofday (&qs_timing_stop, NULL);
	qs_timing_diff = my_difftime (&qs_timing_start, &qs_timing_stop);

	TF_STG2 += ((double)qs_timing_diff->secs + (double)qs_timing_diff->usecs / 1000000);
	free(qs_timing_diff);
#endif

#ifdef USE_8X_MOD_ASM
	align_free(bl_sizes);
	align_free(bl_locs);
#endif

	return;
}


