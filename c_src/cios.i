// Coarsely Integrated Operand Scanning (CIOS) Method
#ifndef __CIOS_I__
#define __CIOS_I__

// a[s], b[s], r[s+2]
static int big_mont_mul_cios(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    int i;
    for (i=0; i < s; i++) {
	UINT_T C = 0;
	UINT_T m;
	int j;
	for(j=0; j < s; j++) {
	    MULAB(a[j],b[i],r[j],C,&C,&r[j]);
	}
	ADD(r[s],C,&C,&r[s]);
	r[s+1] = C;
	C = 0;
	MUL0(r[0],np[0],&m);
	MUL1A(m, n[0], r[0], &C);
	for(j=1; j<s; j++) {
	    MULAB(m,n[j],r[j],C,&C,&r[j-1]);
	}
	ADD(r[s], C, &C, &r[s-1]);
	r[s] = r[s+1] + C;
    }
    return s+1;
}

static int big_mont_sqr_cios(UINT_T* a, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    return big_mont_mul_cios(a, a, np, n, r, s);
}


#endif
