// Finely Integrated Product Scanning (FIPS)
#ifndef __FIPS_I__
#define __FIPS_I__

#include "big2r.i"

// acculmulator size 3 should work to 256 bits (32 bit)
// log64(64*63*128)  ~ 3
// log64(64*63*4096)  ~ 4
// a[s], b[s], m[s+2]!
static int big_mont_mul_fips(PRIVATE UINT_T* a, PRIVATE UINT_T* b,
			     CONST UINT_T* np, CONST UINT_T* n,
			     PRIVATE UINT_T* r,int s)
{
    DECL2(t);
    UINT_T ti;
    UINT_T C;
    UINT_T S;
    int i;

    ZERO2(t);
    ti = 0;
    for(i=0; i < s; i++) {
	int j;
	for (j=0; j < i; j++) {
	    MULA(a[j],b[i-j],ti,&C,&S);
	    ADD21p(t, C);
	    MULA(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    ADD21p(t, C);
	}
	MULA(a[i],b[0],ti,&C,&S);
	ADD21p(t, C);
	MUL0(S,np[0],&r[i]);
	MULA(r[i],n[0], S, &C, &S);
	ADD21p(t, C);
	ti = elem2(t,0);
	SHR2(t);
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    MULA(a[j],b[i-j],ti,&C,&S);
	    ADD21p(t, C);
	    MULA(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    ADD21p(t, C);
	}
	r[i-s] = ti;
	ti = ELEM2(t,0);
	SHR2(t);
    }
    r[s] = ti;
    return s+1;
}

static int big_mont_sqr_fips(PRIVATE UINT_T* a,
			     CONST UINT_T* np, CONST UINT_T* n,
			     PRIVATE UINT_T* r, int s)
{
    return big_mont_mul_fips(a, a, np, n, r, s);
}

#endif
