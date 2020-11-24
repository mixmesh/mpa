#ifndef __MONTMUL_H__
#define __MONTMUL_H__

#include "erl_nif_big.h"


typedef enum {
    REDC_DEFAULT,
    REDC_SOS,
    REDC_FIPS
} redc_type_t;

extern int big_mont_redc(redc_type_t redc_type,
			 ErlNifBigDigit* t, int tl, ErlNifBigDigit* n, int nl,
			 ErlNifBigDigit* np, int npl,
			 ErlNifBigDigit* r, int szr);
extern int big_mont_mul(redc_type_t redc_type,
			ErlNifBigDigit* a, int al,
			ErlNifBigDigit* b, int bl,
			ErlNifBigDigit* n, int nl,
			ErlNifBigDigit* np, int npl,
			ErlNifBigDigit* r, int rzr);
extern int big_mont_sqr(redc_type_t redc_type,
			ErlNifBigDigit* a, int al,
			ErlNifBigDigit* n, int nl,
			ErlNifBigDigit* np, int npl,
			ErlNifBigDigit* r, int szr);
extern int big_mont_pow(redc_type_t redc_type,
			ErlNifBigDigit* a, int al,
			ErlNifBigDigit* e, int el,
			ErlNifBigDigit* p, int pl,
			ErlNifBigDigit* n, int nl,
			ErlNifBigDigit* np, int npl,
			ErlNifBigDigit* r, int szr);
#endif
