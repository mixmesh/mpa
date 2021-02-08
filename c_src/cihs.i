// Coarsely Integrated Hybrid Scanning (CIHS) Method
#ifndef __CIHS_I__
#define __CIHS_I__

// a[s], b[s], r[s+2]
static int big_mont_mul_cihs(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* r, int s)
{
    UINT_T C;
    UINT_T S;
    int i;

    for(i=0; i<s;i++) {
	int j;
	C = 0;
	for(j=0; j<s-i; j++) {
	    mulab(a[j],b[i],r[i+j],C,&C,&S);
	    r[i+j] = S;
	}
	add(r[s],C,&C,&S);
	r[s] = S;
	r[s+1] = C;
    }
    
    for(i=0; i < s; i++) {
	int j;
	UINT_T m;
	mul0(r[0],np[0],&m);
	mula(m, n[0], r[0], &C, &S);
	for(j=1; j<s; j++) {
	    mulab(m,n[j],r[j],C,&C,&S);
	    r[j-1] = S;
	}
	add(r[s],C,&C,&S);
	r[s-1] = S;
	r[s] = r[s+1] + C;
	r[s+1] = 0;
	for (j=i+1; j < s; j++) {
	    mula(b[j], a[s-j+i], r[s-1], &C, &S);
	    r[s-1] = S;
	    add(r[s], C, &C, &S);
	    r[s] = S;
	    r[s+1] = C;
	}
    }
//    printf("r[s-1] = %lu\r\n", r[s-1]);
//    printf("t[s] = %lu\r\n", r[s]);
//    printf("t[s+1] = %lu\r\n", r[s+1]);
//    printf("S = %lu\r\n", S);
//    printf("C = %lu\r\n", C);    
    return s+1;
}

#endif
