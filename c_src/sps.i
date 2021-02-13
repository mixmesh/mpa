// Separated Product Scanning
#ifndef __SPS_I__
#define __SPS_I__

#include "big3r.i"

// P[2s]=0, n[s], np[s], r[s+1]
static int big_mont_redc_sps(UINT_T* P,UINT_T* n,UINT_T* np,UINT_T* r,int s)
{
    DECL3(u);
    UINT_T p1,p0;
    int i, j;

    ZERO3(u);
    for (i = 0; i < s; i++) {
	for (j = 0; j < i; j++) {
	    MUL(r[j],n[i-j],&p1,&p0);
	    ADD32p(u,p1,p0);
	}
	ADD31p(u,P[i]);
	MUL0(ELEM3(u,0), np[0], &r[i]);
	MUL(r[i],n[0],&p1,&p0);
	ADD32p(u,p1,p0);
	SHR3(u);
    }
    for (i = s; i < 2*s-1; i++) {
	for (j = i-s+1; j < s; j++) {
	    MUL(r[j],n[i-j],&p1,&p0);
	    ADD32p(u,p1,p0);
	}
	ADD31p(u,P[i]);
	r[i-s] = ELEM3(u,0);
	SHR3(u);
    }
    ADD31p(u,P[2*s-1]);
    r[s-1] = ELEM3(u,0);
    r[s] = ELEM3(u,1);
    return s+1;
}
//
// a[s], b[s] r[s+1]
// product a[s]*b[s] result is placed in r[s+1] 
static int big_mont_mul_sps(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    UINT_T P[2*s];
    big_zero(P, BIGNUM_SIZE(P));
    big_n_mul(a, b, P, s);
    return big_mont_redc_sps(P, n, np, r, s);
}


static int big_mont_sqr_sps(UINT_T* a, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    UINT_T P[2*s];
    big_zero(P, BIGNUM_SIZE(P));
    big_n_sqr(a, P, s);
    // big_n_mul(a, a, P, s);
    BIGPRINT0("%sa^2,P=", P, 2*s, "");    
    return big_mont_redc_sps(P, n, np, r, s);
}

#endif
