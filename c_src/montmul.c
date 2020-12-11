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

void big_print(UINT_T* x, int xl);

#define BIGPRINT1(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)

#if 0
#define BIGPRINT(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)
#endif
#define BIGPRINT(fmt, x, xl, args...)


#include "digit.i"
#include "big.i"

#include "redc.i"
#include "sos.i"
#include "sps.i"
#include "cios.i"
#include "fios.i"
#include "fips.i"
#include "cihs.i"


void big_print(UINT_T* x, int xl)
{
    int i;
    printf("[%d]{%lu",xl,(unsigned long)x[xl-1]);
    for (i = xl-2; i >= 0; i--) printf(",%lu", (unsigned long)x[i]);
    printf("}");
}



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
    switch(redc_type) {
    case REDC_DEFAULT:
    case REDC_SOS:
    case REDC_SPS: {
	UINT_T T[2*nl+1];
	big_zero(T, BIGNUM_SIZE(T));
	big_mul(a, al, b, bl, T, BIGNUM_SIZE(T));
	return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
    }
    case REDC_CIOS: {
	UINT_T A[nl];
	UINT_T B[nl];
	UINT_T* Ap = a;
	UINT_T* Bp = b;
	int rl;
	if ((al > nl)||(bl > nl))
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	if (bl < nl) {
	    big_copyz(B, nl, b, bl);
	    Bp = B;
	}
	big_zero(r, szr);
	BIGPRINT1("%s: A=", Ap, nl, "cios");
	BIGPRINT1("%s: B=", Bp, nl, "cios");
	rl = big_mont_mul_cios(Ap, Bp, np, n, r, nl);
	BIGPRINT1("%s: A*B=", r, rl, "cios");
	return rl;
    }
    case REDC_FIPS: {
	UINT_T A[nl];
	UINT_T B[nl];
	UINT_T* Ap = a;
	UINT_T* Bp = b;
	if ((al > nl)||(bl > nl))
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	if (bl < nl) {
	    big_copyz(B, nl, b, bl);
	    Bp = B;
	}
	big_zero(r, szr);
	return big_mont_mul_fips(Ap, Bp, np, n, r, nl);
    }
    case REDC_FIOS:{
	UINT_T A[nl];
	UINT_T B[nl];
	UINT_T* Ap = a;
	UINT_T* Bp = b;
	if ((al > nl)||(bl > nl))
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	if (bl < nl) {
	    big_copyz(B, nl, b, bl);
	    Bp = B;
	}
	big_zero(r, szr);
	return big_mont_mul_fios(Ap, Bp, np, n, r, nl);
    }
    case REDC_CIHS:{
	UINT_T A[nl];
	UINT_T B[nl];
	UINT_T* Ap = a;
	UINT_T* Bp = b;
	if ((al > nl)||(bl > nl))
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	if (bl < nl) {
	    big_copyz(B, nl, b, bl);
	    Bp = B;
	}
	big_zero(r, szr);
	return big_mont_mul_cihs(Ap, Bp, np, n, r, nl);
    }
    default:
	return -1;
    }
}

// al < nl < k
int big_mont_sqr(redc_type_t redc_type,
		 UINT_T* a, int al,
		 UINT_T* n, int nl,
		 UINT_T* np, int npl,
		 UINT_T* r, int szr)
{
    switch(redc_type) {
    case REDC_DEFAULT:
    case REDC_SOS:
    case REDC_SPS: {	
	UINT_T T[2*nl+1];
	big_zero(T, BIGNUM_SIZE(T));
	big_sqr(a, al, T, BIGNUM_SIZE(T));
	return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
    }
    case REDC_CIOS: {
	UINT_T A[nl];
	UINT_T* Ap = a;
	if (al > nl)
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	big_zero(r, szr);
	return big_mont_mul_cios(Ap, Ap, np, n, r, nl);
    }
    case REDC_FIPS: {
	UINT_T A[nl];
	UINT_T* Ap = a;
	if (al > nl)
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	big_zero(r, szr);
	return big_mont_mul_fips(Ap, Ap, np, n, r, nl);
    }
    case REDC_FIOS:{
	UINT_T A[nl];
	UINT_T* Ap = a;
	if (al > nl)
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	big_zero(r, szr);
	return big_mont_mul_fios(Ap, Ap, np, n, r, nl);
    }
    case REDC_CIHS:{
	UINT_T A[nl];
	UINT_T* Ap = a;
	if (al > nl)
	    return -1;
	if (al < nl) {
	    big_copyz(A, nl, a, al);
	    Ap = A;
	}
	big_zero(r, szr);	
	return big_mont_mul_cihs(Ap, Ap, np, n, r, nl);
    }
    default:
	return -1;
    }
}


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
    big_copy(P[u], BIGNUM_SIZE(P[0]), p, pl);   // mont(1) !

    BIGPRINT("P%s=", P[u], pl, "");

    s = 0; t = s^1;
    big_copy(A[s], BIGNUM_SIZE(A[0]), a, al);  // check al!

    BIGPRINT("A%s=", A[s], al, "");
    BIGPRINT("E%s=", e, el, "");

    nbits = big_bits(e, el)-1;
    // printf("E #bits = %d\r\n", nbits);
    for (pos = 0; pos < nbits; pos++) {
	int bit = big_bit_test(e, el, pos);
	// printf("E[%d] = %d\r\n", pos, bit);
	if (bit) {
	    // FIXME: make sure al=pl==nl extend A and P
	    //   by zero MSB if needed
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
