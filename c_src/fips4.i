// Finely Integrated Product Scanning (FIPS)
#ifndef __FIPS_I__
#define __FIPS_I__

#include "big2.i"
#include "big3.i"
#include "big4.i"

// acculmulator size 4 should work to 4096 bits
// a[s], b[s], r[s+2]! 
static int big_mont_mul_fips(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* r,int s)
{
    UINT_T t[4];
    UINT_T C;
    UINT_T S;
    int i;

    zero4(t);
    for(i=0; i < s; i++) {
	int j;
	for (j=0; j < i; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    add31(&t[1], C, &t[1]);
	    mula(r[j],n[i-j],S,&C,&S);
	    t[0] = S;
	    add31(&t[1], C, &t[1]);	    
	}
	mula(a[i],b[0],t[0],&C,&S);
	add31(&t[1], C, &t[1]);
	mul0(S,np[0],&r[i]);
	mula(r[i],n[0], S, &C, &S);
	add31(&t[1], C, &t[1]);	
	shr4(t);
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    add31(&t[1], C, &t[1]);
	    mula(r[j],n[i-j],S,&C,&S);
	    t[0] = S;
	    add31(&t[1], C, &t[1]);
	}
	r[i-s] = t[0];
	shr4(t);
    }
    r[s] = t[0];
    return s+1;
}

static int big_mont_sqr_fips(UINT_T* a, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    return big_mont_mul_fips(a, a, np, n, r, s);
}

#endif
