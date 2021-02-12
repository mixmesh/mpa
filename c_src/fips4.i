// Finely Integrated Product Scanning (FIPS)
#ifndef __FIPS_I__
#define __FIPS_I__

#include "big3r.i"

// acculmulator size 4 should work to 4096 bits
// a[s], b[s], r[s+2]! 
static int big_mont_mul_fips(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* r,int s)
{
    decl3(t);
    UINT_T ti;
    UINT_T C;
    UINT_T S;
    int i;

    zero3(t);
    ti = 0;
    for(i=0; i < s; i++) {
	int j;
	for (j=0; j < i; j++) {
	    mula(a[j],b[i-j],ti,&C,&S);
	    add31p(t, C);
	    mula(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    add31p(t, C);
	}
	mula(a[i],b[0],ti,&C,&S);
	add31p(t, C);
	mul0(S,np[0],&r[i]);
	mula(r[i],n[0], S, &C, &S);
	add31p(t, C);
	ti = elem3(t,0);
	shr3(t);
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    mula(a[j],b[i-j],ti,&C,&S);
	    add31p(t, C);
	    mula(r[j],n[i-j],S,&C,&S);
	    ti = S;
	    add31p(t, C);
	}
	r[i-s] = ti;
	ti = elem3(t,0);	
	shr3(t);
    }
    r[s] = ti;
    return s+1;
}

static int big_mont_sqr_fips(UINT_T* a, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    return big_mont_mul_fips(a, a, np, n, r, s);
}

#endif
