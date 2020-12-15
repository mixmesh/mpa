#ifndef __MONTMUL_H__
#define __MONTMUL_H__

#include "erl_nif_big.h"

typedef enum {
    REDC_SOS,
    REDC_SPS,
    REDC_CIOS,
    REDC_FIPS,
    REDC_FIOS,
    REDC_CIHS,
} redc_type_t;

extern int big_mont_mul(redc_type_t redc_type,
			ErlNifBigDigit* a, ErlNifBigDigit* b,
			ErlNifBigDigit* n, ErlNifBigDigit* np,
			ErlNifBigDigit* r, int s);
extern int big_mont_sqr(redc_type_t redc_type,
			ErlNifBigDigit* a,
			ErlNifBigDigit* n, ErlNifBigDigit* np,
			ErlNifBigDigit* r, int s);
extern int big_mont_pow(redc_type_t redc_type,
			ErlNifBigDigit* a, 
			ErlNifBigDigit* e, int el,
			ErlNifBigDigit* p,  ErlNifBigDigit* n, 
			ErlNifBigDigit* np, ErlNifBigDigit* r, int s);
#endif
