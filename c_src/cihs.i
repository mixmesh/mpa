// Coarsely Integrated Hybrid Scanning (CIHS) Method
#ifndef __CIHS_I__
#define __CIHS_I__

static int big_mont_mul_cihs(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* t, int s)
{
    UINT_T C;
    UINT_T S;
    int i;

    for(i=0; i<s;i++) {
	int j;
	C = 0;
	for(j=0; j<s; j++) {
	    mulab(a[j],b[i],t[i+j],C,&C,&S);
	    t[i+j] = S;
	}
	add(t[s],C,&C,&S);
	t[s] = S;
	t[s+1] = C;
    }
    
    for(i=0; i < s; i++) {
	int j;
	UINT_T m;
	mul0(t[0],np[0],&m);
	mula(m, n[0], t[0], &C, &S);
	for(j=1; j<s; j++) {
	    mulab(m,n[j],t[j],C,&C,&S);
	    t[j-1] = S;
	}
	add(t[s],C,&C,&S);
	t[s-1] = S;
	t[s] = t[s+1] + C;
	t[s+1] = 0;
	for (j=i+1; j < s; j++) {
	    mula(b[j], a[s-j+i], t[s-1], &C, &S);
	    t[s-1] = S;
	    add(t[s], C, &C, &S);
	    t[s] = S;
	    t[s+1] = C;
	}
    }
    return s+1;
}

#endif
