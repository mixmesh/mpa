// Coarsely Integrated Hybrid Scanning (CIHS) Method
// From 1996/01/j37acmon.pdf (buggy)
#ifndef __CIHS_I__
#define __CIHS_I__

// a[s], b[s], r[s+2]
static int big_mont_mul_cihs(PRIVATE UINT_T* a, PRIVATE UINT_T* b,
			     CONST UINT_T* np, CONST UINT_T* n,
			     PRIVATE UINT_T* r, int s)
{
    int i;

    for (i=0; i<s; i++) {
	UINT_T C = 0;
	int j;
	for(j=0; j < s-i; j++) {
	    MULAB(a[j],b[i],r[i+j],C,&C,&r[i+j]);
	}
	ADD(r[s],C,&C,&r[s]);
	ADD(r[s+1],C,&C,&r[s+1]); // &C is not used
    }
    
    for (i=0; i < s; i++) {
	UINT_T C;
	UINT_T m;
	int j;
	MUL0(r[0],np[0],&m);
	MUL1A(m, n[0], r[0], &C);
	for(j=1; j < s; j++) {
	    MULAB(m,n[j],r[j],C,&C,&r[j-1]);
	}
	ADD(r[s],C,&C,&r[s-1]);
	ADD(r[s+1],C,&C,&r[s]);  // &C is not used
	r[s+1] = 0;
	
	for (j=i+1; j < s; j++) {
	    MULA(b[j], a[s-j+i], r[s-1], &C, &r[s-1]);
	    ADD(r[s], C, &C, &r[s]);
	    ADD(r[s+1], C, &C, &r[s+1]);  // &C is not used
	}
    }
    return s+1;
}

static int big_mont_sqr_cihs(PRIVATE UINT_T* a,
			     CONST UINT_T* np, CONST UINT_T* n,
			     PRIVATE UINT_T* r, int s)
{
    return big_mont_mul_cihs(a, a, np, n, r, s);
}

#endif
