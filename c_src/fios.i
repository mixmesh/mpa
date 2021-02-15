// Finly Integrated Operand Scanning (FIOS) Method
#ifndef __FIOS_I__
#define __FIOS_I__

// a[s], b[s], r[s+2]!
STATIC INLINE int big_mont_mul_fios(PRIVATE UINT_T* a, PRIVATE UINT_T* b,
				    CONST UINT_T* n, CONST UINT_T* np,
				    PRIVATE UINT_T* r,int s)
{
    int i;
    for(i=0; i < s; i++) {
	UINT_T m;
	UINT_T C;
	UINT_T S;
	int j;

	MULA(a[0], b[i], r[0], &C, &S);
	big_addc(r, 1, s+2, C);
	MUL0(S, np[0], &m);
	MULA(m, n[0], S, &C, &S);

	for (j = 1; j < s; j++) {
	    MULAB(a[j],b[i],r[j],C,&C,&S);
	    big_addc(r, j+1, s+2, C);
	    MULA(m, n[j], S, &C, &r[j-1]);
	}
	ADD(r[s], C, &C, &r[s-1]);
	r[s] = r[s+1] + C;
	r[s+1] = 0;
    }
    return s+1;
}

STATIC INLINE int big_mont_sqr_fios(PRIVATE UINT_T* a,
				    CONST UINT_T* n, CONST UINT_T* np,
				    PRIVATE UINT_T* r, int s)
{
    return big_mont_mul_fios(a, a, n, np, r, s);
}

#endif
