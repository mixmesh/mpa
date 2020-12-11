// Coarsely Integrated Operand Scanning (CIOS) Method
#ifndef __CIOS_I__
#define __CIOS_I__

// a[s], b[s], t[s+2]!
static int big_mont_mul_cios(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* t,int s)
{
    int i;
    for (i=0; i < s; i++) {
	UINT_T C = 0;
	UINT_T m;
	int j;
	for(j=0; j < s; j++) {
	    mulab(a[j],b[i],t[j],C,&C,&t[j]);
	}
	add(t[s],C,&C,&t[s]);
	t[s+1] = C;
	C = 0;
	mul0(t[0],np[0],&m);
	for(j=0; j<s; j++) {
	    mulab(m,n[j],t[j],C,&C,&t[j]);
	}
	add(t[s], C, &C, &t[s]);
	t[s+1] += C;
	for(j=0; j <= s; j++)
	    t[j] = t[j+1];
	mul0(t[0],np[0], &m);
	mul1a(m, n[0], t[0], &C);
	for(j=1; j < s; j++) {
	    mulab(m, n[j], t[j], C, &C, &t[j-1]);
	}
	add(t[s], C, &C, &t[s-1]);
	t[s] = t[s+1] + C;
    }
    return s+1; // length?
}

#endif
