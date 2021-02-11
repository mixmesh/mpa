// Coarsely Integrated Hybrid Scanning (CIHS) Method
// From 1996/01/j37acmon.pdf (buggy)
#ifndef __CIHS_I__
#define __CIHS_I__

// a[s], b[s], r[s+2]
static int big_mont_mul_cihs(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    int i;

    for (i=0; i<s; i++) {
	UINT_T C = 0;
	int j;
	for(j=0; j < s-i; j++) {
	    mulab(a[j],b[i],r[i+j],C,&C,&r[i+j]);
	}
	add(r[s],C,&C,&r[s]);
	add(r[s+1],C,&C,&r[s+1]); // &C is not used
    }
    
    for (i=0; i < s; i++) {
	UINT_T C;
	UINT_T m;
	int j;
	mul0(r[0],np[0],&m);
	mul1a(m, n[0], r[0], &C);
	for(j=1; j < s; j++) {
	    mulab(m,n[j],r[j],C,&C,&r[j-1]);
	}
	add(r[s],C,&C,&r[s-1]);
	add(r[s+1],C,&C,&r[s]);  // &C is not used
	r[s+1] = 0;
	
	for (j=i+1; j < s; j++) {
	    mula(b[j], a[s-j+i], r[s-1], &C, &r[s-1]);
	    add(r[s], C, &C, &r[s]);
	    add(r[s+1], C, &C, &r[s+1]);  // &C is not used
	}
    }
    return s+1;
}

static int big_mont_sqr_cihs(UINT_T* a, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    return big_mont_mul_cihs(a, a, np, n, r, s);
}

#endif
