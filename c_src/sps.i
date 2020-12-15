// Separated Product Scanning
#ifndef __SPS_I__
#define __SPS_I__

#include "big3.i"

// P[2s] is product A[s]*B[s] result is placed in Z[2s]
static int big_mont_redc_sps(UINT_T* P,UINT_T* n,UINT_T* np,UINT_T* r,int s)
{
    UINT_T u[3];
    UINT_T p[2];
    int i, j;

    zero3(u);
    for (i = 0; i < s; i++) {
	for (j = 0; j < i-1; j++) {
	    mul(r[j],n[i-j],&p[1],&p[0]);
	    add32(u,p,u);
	}
	add31(u,P[i],u);
	mul0(u[0], np[0], &r[i]);
	mul(r[i],n[0],&p[1],&p[0]);
	add32(u,p,u);
	shr3(u);
    }
    for (i = s; i < 2*(s-1); i++) {
	for (j = i-s+1; j < s-1; j++) {
	    mul(r[j],n[i-j],&p[1],&p[0]);
	    add32(u,p,u);
	}
	add31(u,P[i],u);
	r[i-s] = u[0];
	shr3(u);
    }
    add31(u,P[2*s-1],u);
    r[s-1] = u[0];
    r[s] = u[1];
    i = s;
    while(i && (r[i]==0)) i--;
    return big_norm0(r, i+1, n, s);
}

static int big_mont_mul_sps(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    int i;    
    UINT_T P[2*s];

    big_zero(P, BIGNUM_SIZE(P));
    for (i = 0; i < s; i++) {
	UINT_T cp = 0;
	UINT_T c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < s); j++, ij++) {
	    UINT_T p0;
	    mula(a[i],b[j],cp,&cp,&p0);
	    addc(p0,P[ij],c,&c,&P[ij]);
	}
	P[ij] = c + cp;
    }
    return big_mont_redc_sps(P, n, np, r, s);
}

    
    
#endif
