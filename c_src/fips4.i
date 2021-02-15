// Finely Integrated Product Scanning (FIPS)
#ifndef __FIPS_I__
#define __FIPS_I__

#include "big3r.i"

// acculmulator size 4 should work to 4096 bits
// a[s], b[s], r[s+2]! 
static int big_mont_mul_fips(PRIVATE UINT_T* a, PRIVATE UINT_T* b,
			     CONST UINT_T* n, CONST UINT_T* np, 
			     PRIVATE UINT_T* r,int s)
{
    DECL3(t);
    UINT_T ti;
    UINT_T C;
    UINT_T S;
    int i;

    ZERO3(t);
    ti = 0;
    for(i=0; i < s; i++) {
	int j;
	for (j=0; j < i; j++) {
	    MULA(a[j],b[i-j],ti,&C,&S);
	    ADD31P(t, C);
	    MULA(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    ADD31P(t, C);
	}
	MULA(a[i],b[0],ti,&C,&S);
	ADD31P(t, C);
	MUL0(S,np[0],&r[i]);
	MULA(r[i],n[0], S, &C, &S);
	ADD31P(t, C);
	ti = ELEM3(t,0);
	SHR3(t);
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    MULA(a[j],b[i-j],ti,&C,&S);
	    ADD31P(t, C);
	    MULA(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    ADD31P(t, C);
	}
	r[i-s] = ti;
	ti = ELEM3(t,0);	
	SHR3(t);
    }
    r[s] = ti;
    return s+1;
}

static int big_mont_sqr_fips(PRIVATE UINT_T* a,
			     CONST UINT_T* n, CONST UINT_T* np,
			     PRIVATE UINT_T* r, int s)
{
    return big_mont_mul_fips(a, a, n, np, r, s);
}

#endif
