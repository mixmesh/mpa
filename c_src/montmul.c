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

#define BIGPRINT0(fmt, x, xl, args...)

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
// for size to 4k use fips3 for 64-bit and fips4 32-bit (upto 2^18 64bitdigits)
#include "fips3.i"
#include "cihs.i"

void big_print(UINT_T* x, int xl)
{
    int i;
    printf("[%d]{\r\n", xl);
    i = 0;
    while(i+4 < xl) {
	printf("  ");
	printf("%016lx ", (unsigned long)x[i]);
	printf("%016lx ", (unsigned long)x[i+1]);
	printf("%016lx ", (unsigned long)x[i+2]);
	printf("%016lx ", (unsigned long)x[i+3]);
	i += 4;
	printf("\r\n");
    }
    if (i < xl) {
	printf("  ");
	while(i < xl) {
	    printf("%016lx ", (unsigned long)x[i]);
	    i++;
	}
    }
    printf("}\r\n");
}

int big_mont_norm(UINT_T* r, int rl, UINT_T* n, int nl)
{
    rl = big_trim(r, rl);
    // return big_norm1(r, n, nl);
    return big_norm0(r, rl, n, nl);
}

// al,bl < nl < k   R[al+bl]
int big_mont_mul(redc_type_t redc_type, UINT_T* a, UINT_T* b,
		 UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    int rl;

    BIGPRINT0("%sa=", a, s, "");
    BIGPRINT0("%sb=", b, s, "");
    
    switch(redc_type) {
    case REDC_SOS:
	rl = big_mont_mul_sos(a, b, n, np, r, s);
	break;
    case REDC_SPS:
	rl = big_mont_mul_sps(a, b, n, np, r, s);
	break;
    case REDC_CIOS:
	big_zero(r, s+1);
	rl = big_mont_mul_cios(a, b, n, np, r, s);
	break;
    case REDC_FIPS:
	big_zero(r, s+1);
	rl = big_mont_mul_fips(a, b, n, np, r, s);
	break;
    case REDC_FIOS:
	big_zero(r, s+2);
	rl = big_mont_mul_fios(a, b, n, np, r, s);
	break;
    case REDC_CIHS:
	big_zero(r, s+2);
	rl = big_mont_mul_cihs(a, b, n, np, r, s);
	break;
    default:
	return -1;
    }
    BIGPRINT0("%sa*b=", r, rl, "");
    rl = big_mont_norm(r, rl, n, s);
    BIGPRINT0("%sa*b'=", r, rl, "");
    return rl;
}

// al < nl < k
int big_mont_sqr(redc_type_t redc_type,
		 UINT_T* a, UINT_T* n,  UINT_T* np, UINT_T* r, int s)
{
    int rl;
	
    BIGPRINT0("%sa=", a, s, "");
    switch(redc_type) {
    case REDC_SOS:
	rl = big_mont_sqr_sos(a, n, np, r, s);
	break;
    case REDC_SPS:
	rl = big_mont_sqr_sps(a, n, np, r, s);
	break;
    case REDC_CIOS:
	big_zero(r, s+1);
	rl = big_mont_sqr_cios(a, n, np, r, s);
	break;
    case REDC_FIPS:
	big_zero(r, s+1);
	rl = big_mont_sqr_fips(a, n, np, r, s);
	break;
    case REDC_FIOS:
	big_zero(r, s+2);
	rl = big_mont_sqr_fios(a, n, np, r, s);
	break;
    case REDC_CIHS:
	big_zero(r, s+2);
	rl = big_mont_sqr_cihs(a, n, np, r, s);
	break;
    default:
	return -1;
    }
    BIGPRINT0("%sa^2=", r, rl, "");
    rl = big_mont_norm(r, rl, n, s);
    BIGPRINT0("%sa^2'=", r, rl, "");
    return rl;
}

// test function
// a[s] r[2s]
int big_mod2_sqr(UINT_T* a, UINT_T* r, int s)
{
    int rl;
    big_zero(r, 2*s);
    rl = big_n_sqr(a, r, s);
    return big_trim(r, rl);
}


// a^x (mod R)  (R=B^k > N)  (B = UINT_T base) r[s+2]
int big_mont_pow(redc_type_t redc_type,
		 UINT_T* a,
		 UINT_T* x, int xl,
		 UINT_T* p, UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    UINT_T P[2][s+2];
    UINT_T A[2][s+2];
    int u, v;
    int pos, xbits;
    int j;
    int rl;

    BIGPRINT1("%sN=", n, s, "");
    BIGPRINT1("%sN'=", np, s, "");

    u = 0;
    big_copy(P[u], p, s);   // = mont(1) !
    BIGPRINT1("%s1=P[u]=", P[u], s, "");
    
    v = 0;
    big_copy(A[v], a, s);  // check al!
    BIGPRINT1("%sA=A[u]=", A[v], s, "");

    xbits = big_bits(x, xl);

    for (pos=1, j=0; pos < xbits; pos += D_SIZE, j++) {
	UINT_T xj = x[j];
	int size = (pos+D_SIZE < xbits) ? D_SIZE : xbits-pos;
	int bit;
	for (bit = 0; bit < size; bit++) {
	    if (xj & 1) {
		big_mont_mul(redc_type,A[v],P[u],n,np,P[!u],s);
		BIGPRINT1("%sP[!u]=", P[!u], s, "");
		u = !u;
	    }
	    // A' = A^2 (mod R)
	    big_mont_sqr(redc_type,A[v],n,np,A[!v],s);
	    BIGPRINT1("%sA[!v]=", A[!v], s, "");
	    v = !v;
	    
	    xj >>= 1;
	}
    }
    // r = A*P (mod R)
    rl = big_mont_mul(redc_type,A[v],P[u],n,np,r,s);
    BIGPRINT1("%sR=", r, s, "");
    return rl;
}
