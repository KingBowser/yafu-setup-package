#ifndef __ASM_REDC_H__
#define __ASM_REDC_H__

#include "config.h"
#include <gmp.h>

/* Signals that we have assembly code for 1xN mul/redc */
#define HAVE_NATIVE_MULREDC1_N
/* Signals that we have assembly code for variable size redc */
#define HAVE_ASM_REDC3

/* Call the mulredc*() function with MS Windows parameter passing if
  WINDOWS64_ABI is defined. This is useful for testing the functions
  with Microsoft ABI under Linux */
#ifdef WINDOWS64_ABI
#define MULREDC_ABI __attribute__((ms_abi))
#else
#define MULREDC_ABI
#endif

extern void ecm_redc3(mp_limb_t *, const mp_limb_t *, mp_size_t, mp_limb_t) MULREDC_ABI;


/* WARNING: the size-1 version doesn't take pointers in input */
extern mp_limb_t mulredc1(mp_limb_t *, mp_limb_t, mp_limb_t, mp_limb_t, mp_limb_t) MULREDC_ABI;

extern mp_limb_t mulredc2(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc3(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc4(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc5(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc6(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc7(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc8(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc9(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc10(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc11(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc12(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc13(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc14(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc15(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc16(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc17(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc18(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc19(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc20(mp_limb_t *, const mp_limb_t *, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;

extern mp_limb_t mulredc1_2(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_3(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_4(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_5(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_6(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_7(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_8(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_9(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_10(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_11(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_12(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_13(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_14(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_15(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_16(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_17(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_18(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_19(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;
extern mp_limb_t mulredc1_20(mp_limb_t *, const mp_limb_t, const mp_limb_t *, const mp_limb_t *, mp_limb_t) MULREDC_ABI;

#endif
