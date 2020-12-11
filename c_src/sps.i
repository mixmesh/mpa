// Separated Product Scanning
#ifndef __SPS_I__
#define __SPS_I__

#include "big3.i"

int big_mont_redc_sps(UINT_T* P, int pl,
		      UINT_T* n, int s,
		      UINT_T* np, int npl,
		      UINT_T* Z, int szZ)
{
    UINT_T u[3];
    UINT_T p[2];
    int i, j;

    zero3(u);
    for (i = 0; i < s; i++) {
	for (j = 0; j < i-1; j++) {
	    mul(Z[j],n[i-j],&p[1],&p[0]);
	    add32(u,p,u);
	}
	add31(u,P[i],u);
	mul0(u[0], np[0], &Z[i]);
	mul(Z[i],n[0],&p[1],&p[0]);
	add32(u,p,u);
	shr3(u);
    }
    for (i = s; i < 2*(s-1); i++) {
	for (j = i-s+1; j < s-1; j++) {
	    mul(Z[j],n[i-j],&p[1],&p[0]);
	    add32(u,p,u);
	}
	add31(u,P[i],u);
	Z[i-s] = u[0];
	shr3(u);
    }
    add31(u,P[2*s-1],u);
    Z[s-1] = u[0];
    Z[s] = u[1];
    i = s;
    while(i && (Z[i]==0)) i--;
    return big_norm0(Z, i+1, n, s);
}

#endif
