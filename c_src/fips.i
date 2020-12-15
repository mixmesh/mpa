// Finely Integrated Product Scanning (FIPS)
#ifndef __FIPS_I__
#define __FIPS_I__

#include "big3.i"

// a[s], b[s], m[s+2]! (s < W=sizeof(UINT_T)*8)
static int big_mont_mul_fips(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* m,int s)
{
    UINT_T t[3]; // FIXME: 2+logw(S) S=1024,W=32 => 2+2 = 4!
    UINT_T C;
    UINT_T S;
    int i;

    for(i=0; i < s; i++) {
	int j;
	for (j=0; j < i; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    big_addc(t, 1, C);
	    mula(m[j],n[i-j],S,&C,&S);
	}
	mula(a[i],b[0],t[0],&C,&S);
	big_addc(t, 1, C);
	mul0(S,np[0],&m[i]);
	mula(m[i],n[0], S, &C, &S);
	big_addc(t, 1, C);
	shr3(t);
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    big_addc(t, 1, C);
	    mula(m[j],n[i-j],S,&C,&t[0]);
	    big_addc(t,1,C);
	}
	m[i-s] = t[0];
	shr3(t);
    }
    return s;
}

#endif
