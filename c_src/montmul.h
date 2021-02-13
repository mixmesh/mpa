#ifndef __MONTMUL_H__
#define __MONTMUL_H__

#include "erl_nif_big.h"

#define GLOBAL
#define LOCAL
#define PRIVATE
#define CONST   const
#define INLINE  inline

#define UINT_T  ErlNifBigDigit
#define UINTH_T ErlNifBigHalfDigit
#ifdef BIG_HAVE_DOUBLE_DIGIT
#define UINTD_T ErlNifBigDoubleDigit
#endif
#define INT_T  ErlNifBigSignedDigit

typedef enum {
    REDC_SOS,
    REDC_SPS,
    REDC_CIOS,
    REDC_FIPS,
    REDC_FIOS,
    REDC_CIHS,
} redc_type_t;

extern void big_print(UINT_T* x, int xl);
extern int big_mont_norm(UINT_T* r, int rl, UINT_T* n, int nl);

extern int big_mont_mul(redc_type_t redc_type, UINT_T* a, UINT_T* b,
			UINT_T* n, UINT_T* np, UINT_T* r, int s);
extern int big_mont_sqr(redc_type_t redc_type,
			UINT_T* a, UINT_T* n,  UINT_T* np, UINT_T* r, int s);
extern int big_mont_pow(redc_type_t redc_type,
			UINT_T* a,
			UINT_T* e, int el,
			UINT_T* p, UINT_T* n, UINT_T* np, UINT_T* r, int s);
// test
extern int big_mod2_sqr(UINT_T* a, UINT_T* r, int s);

#endif
