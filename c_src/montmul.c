//
// unsigned montgomery multiplication
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define NDEBUG
#include <assert.h>

#include "montmul.h"

#define UINT_T  ErlNifBigDigit
#define UINTH_T ErlNifBigHalfDigit
#ifdef BIG_HAVE_DOUBLE_DIGIT
#define UINTD_T ErlNifBigDoubleDigit
#endif

#include "digit.i"
#include "big.i"

#include "redc.i"
#include "sos.i"
#include "sps.i"
#include "cios.i"
#include "fios.i"
#include "fips.i"
#include "cihs.i"

int big_mont_redc(redc_type_t redc_type,
		  UINT_T* t, int tl, UINT_T* n, int nl,
		  UINT_T* np, int npl, UINT_T* r, int szr)
{
    switch(redc_type) {
    case REDC_DEFAULT:
	return big_mont_redc_default(t, tl, n, nl, np, npl, r, szr);
    case REDC_SOS:
	return big_mont_redc_sos(t, tl, n, nl, np, npl, r, szr);
    case REDC_SPS:
	return big_mont_redc_sps(t, tl, n, nl, np, npl, r, szr);
    default:
	return -1;
    }
}

// al,bl < nl < k   R[al+bl]
int big_mont_mul(redc_type_t redc_type,
		 UINT_T* a, int al, UINT_T* b, int bl,
		 UINT_T* n, int nl,
		 UINT_T* np, int npl,
		 UINT_T* r, int szr)
{
    UINT_T T[2*nl];
    big_zero(T, 2*nl);
    big_mul(a, al, b, bl, T, BIGNUM_SIZE(T));
    return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
}

// al < nl < k
int big_mont_sqr(redc_type_t redc_type,
		 UINT_T* a, int al,
		 UINT_T* n, int nl,
		 UINT_T* np, int npl,
		 UINT_T* r, int szr)
{
    UINT_T T[2*nl];
    big_zero(T, 2*nl);
    big_sqr(a, al, T, BIGNUM_SIZE(T));
    return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
}

#define BIGPRINT1(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)

#if 0
#define BIGPRINT(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)
#endif
#define BIGPRINT(fmt, x, xl, args...)

// a^e (mod R)  (R=B^k > N)  (B = UINT_T base) r[2*al]
int big_mont_pow(redc_type_t redc_type,
		 UINT_T* a, int al,
		 UINT_T* e, int el,
		 UINT_T* p, int pl,
		 UINT_T* n, int nl,
		 UINT_T* np, int npl,
		 UINT_T* R, int szR)
{
    UINT_T P[2][2*nl+1];
    UINT_T A[2][2*nl+1];
    int u, v;
    int s, t;
    int pos, nbits;
    int rl;

    assert(al <= nl);
    assert(npl <= nl);
    assert(pl <= nl);
    
    u = 0; v = u^1;
    big_copy(P[u], p, pl);   // mont(1) !

    BIGPRINT("P%s=", P[u], pl, "");

    s = 0; t = s^1;
    big_copy(A[s], a, al);  // check al!

    BIGPRINT("A%s=", A[s], al, "");
    BIGPRINT("E%s=", e, el, "");

    nbits = big_bits(e, el)-1;
    // printf("E #bits = %d\r\n", nbits);
    for (pos = 0; pos < nbits; pos++) {
	int bit = big_bit_test(e, el, pos);
	// printf("E[%d] = %d\r\n", pos, bit);
	if (bit) {
	    // P' = A*P (mod R)
	    pl = big_mont_mul(redc_type,
			      A[s],al, P[u],pl, n,nl, np,npl,
			      P[v],BIGNUM_SIZE(P[v]));
	    u = v; v = u^1;
	    BIGPRINT("P%d = ", P[u], pl, pos);
	}
	// A' = A^2 (mod R)
	al = big_mont_sqr(redc_type, A[s], al,
			  n, nl, np, npl, A[t], BIGNUM_SIZE(A[t]));
	s = t; t = s^1;
	BIGPRINT("A%d = ", A[s], al, pos);
    }
    // r = A*P (mod R)
    //printf("al=%d, pl=%d, nl=%d, npl=%d Rsize=%d\r\n", al, pl, nl, npl, 2*nl);
    rl = big_mont_mul(redc_type, A[s], al, P[u], pl, n, nl, np, npl, R, szR);
    BIGPRINT("R%s = ", R, rl, "");
    return rl;
}
