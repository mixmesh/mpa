// Finly Integrated Operand Scanning (FIOS) Method
#ifndef __FIOS_I__
#define __FIOS_I__

// include "unroll.i"


// a[s], b[s], t[s+3]!
static int big_mont_mul_fios(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* t,int s)
{
    int i;
    for(i=0; i < s; i++) {
	UINT_T m;
	UINT_T C;
	UINT_T S;
	int j;

	mula(a[0], b[i], t[0], &C, &S);
	big_addc(t, 1, C);
	mul0(S, np[0], &m);
	mula(m, n[0], S, &C, &S);

	for (j = 1; j < s; j++) {
	    mulab(a[j],b[i],t[j],C,&C,&S);
	    big_addc(t, j+1, C);
	    mula(m, n[j], S, &C, &t[j-1]);
	}
	add(t[s], C, &C, &t[s-1]);
	t[s] = t[s+1] + C;
	t[s+1] = 0;
    }
    return s+1;
}

#endif
