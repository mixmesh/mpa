// Coarsely Integrated Operand Scanning (CIOS) Method
#ifndef __CIOS_I__
#define __CIOS_I__

// a[s], b[s], r[s+2]
STATIC INLINE int big_mont_mul_cios(PRIVATE UINT_T* a, PRIVATE UINT_T* b,
				    CONST UINT_T* n, CONST UINT_T* np,
				    PRIVATE UINT_T* r, int s)
{
    int i;
UNROLL    
    for (i=0; i < s; i++) {
	UINT_T C = 0;
	UINT_T m;
	int j;
	UNROLL
	for(j=0; j < s; j++) {
	    MULAB(a[j],b[i],r[j],C,&C,&r[j]);
	}
	ADD(r[s],C,&C,&r[s]);
	r[s+1] = C;
	C = 0;
	MUL0(r[0],np[0],&m);
	MUL1A(m, n[0], r[0], &C);
	UNROLL
	for(j=1; j<s; j++) {
	    MULAB(m,n[j],r[j],C,&C,&r[j-1]);
	}
	ADD(r[s], C, &C, &r[s-1]);
	r[s] = r[s+1] + C;
    }
    return s+1;
}

STATIC INLINE int big_mont_sqr_cios(PRIVATE UINT_T* a,
				    CONST UINT_T* n, CONST UINT_T* np, 
				    PRIVATE UINT_T* r, int s)
{
    return big_mont_mul_cios(a, a, n, np, r, s);
}


#endif
