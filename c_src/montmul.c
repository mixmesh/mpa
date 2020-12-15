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

int big_mont_norm(UINT_T* r, int rl, UINT_T* n, int nl)
{
    while((rl>1) && (r[rl-1]==0))
	rl--;
    return big_norm0(r, rl, n, nl);
}

// al,bl < nl < k   R[al+bl]
int big_mont_mul(redc_type_t redc_type, UINT_T* a, UINT_T* b,
		 UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    int rl;

    BIGPRINT1("%sa=", a, s, "");
    BIGPRINT1("%sb=", b, s, "");
    
    switch(redc_type) {
    case REDC_SOS:
	rl = big_mont_mul_sos(a, b, np, n, r, s);
	break;
	return big_mont_norm(r, rl, n, s);
    case REDC_SPS:
	rl = big_mont_mul_sps(a, b, np, n, r, s);
	break;
    case REDC_CIOS:
	big_zero(r, s+1);
	rl = big_mont_mul_cios(a, b, np, n, r, s);
	break;
    case REDC_FIPS:
	big_zero(r, s+1);
	rl = big_mont_mul_fips(a, b, np, n, r, s);
	break;
    case REDC_FIOS:
	big_zero(r, s+2);
	rl = big_mont_mul_fios(a, b, np, n, r, s);
	break;
    case REDC_CIHS:
	big_zero(r, 2*s);
	rl = big_mont_mul_cihs(a, b, np, n, r, s);
	break;
    default:
	return -1;
    }
    BIGPRINT1("%sa*b=", r, rl, "");    
    return big_mont_norm(r, rl, n, s);
}

// al < nl < k
int big_mont_sqr(redc_type_t redc_type,
		 UINT_T* a, UINT_T* n,  UINT_T* np, UINT_T* r, int s)
{
    int rl;
	
    BIGPRINT1("%sa=", a, s, "");
    switch(redc_type) {
    case REDC_SOS:
	rl = big_mont_mul_sos(a, a, np, n, r, s);
	break;
    case REDC_SPS:
	rl = big_mont_mul_sps(a, a, np, n, r, s);
	break;
    case REDC_CIOS:
	big_zero(r, s);
	rl = big_mont_mul_cios(a, a, np, n, r, s);
	break;
    case REDC_FIPS:
	big_zero(r, s);
	rl = big_mont_mul_fips(a, a, np, n, r, s);
	break;
    case REDC_FIOS:
	big_zero(r, s);
	rl = big_mont_mul_fios(a, a, np, n, r, s);
	break;
    case REDC_CIHS:
	big_zero(r, s);	
	rl = big_mont_mul_cihs(a, a, np, n, r, s);
	break;
    default:
	return -1;
    }
    BIGPRINT1("%sa^2=", r, rl, "");
    rl = big_mont_norm(r, rl, n, s);
    BIGPRINT1("%sa^2'=", r, rl, "");
    return rl;
}


// a^e (mod R)  (R=B^k > N)  (B = UINT_T base) r[2*al]
int big_mont_pow(redc_type_t redc_type,
		 UINT_T* a,
		 UINT_T* e, int el,
		 UINT_T* p, UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    UINT_T P[2][2*s+1];
    UINT_T A[2][2*s+1];
    int u, v;
    int c, d;
    int pos, nbits;
    int rl;

    u = 0; v = u^1;
    big_copy(P[u], p, s);   // = mont(1) !

    c = 0; d = c^1;
    big_copy(A[c], a, s);  // check al!

    nbits = big_bits(e, el)-1;

    for (pos = 0; pos < nbits; pos++) {
	int bit = big_test(e, el, pos);
	if (bit) {
	    big_mont_mul(redc_type,A[c],P[u],n,np,P[v],s);
	    u = v; v = u^1;
	}
	// A' = A^2 (mod R)
	big_mont_sqr(redc_type,A[c],n,np,A[d],s);
	c = d; d = c^1;
    }
    // r = A*P (mod R)
    rl = big_mont_mul(redc_type, A[c], P[u], n, np, r, s);
    return rl;
}
