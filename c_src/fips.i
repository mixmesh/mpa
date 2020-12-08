// Finely Integrated Operand Scanning (COIS) Method
#ifndef __FIPS_I__
#define __FIPS_I__

// a[s], b[s], t[s+2]!
static int big_mont_mul_fips(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			     UINT_T* m,int s)
{
    UINT_T t[3];
    UINT_T C;
    UINT_T S;
    int i;
    for(i=0; i < s; i++) {
	int j;    
	for (j=0; j < i; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    big_addc(t, 1, C);
	    mula(m[j],n[i-j],S,&C,&S);
	}
	mula(a[i],b[0],t[0],&C,&S);
	big_addc(t, 1, C);
	mul0(S,np[0],&m[i]);
	mula(m[i],n[0], S, &C, &S);
	big_addc(t, 1, C);
	t[0] = t[1];
	t[1] = t[2];
	t[2] = 0;
    }

    for (i=s; i < 2*s; i++) {
	int j;
	for (j=i-s+1; j < s; j++) {
	    mula(a[j],b[i-j],t[0],&C,&S);
	    big_addc(t, 1, C);
	    mula(m[j],n[i-j],S,&C,&S);
	    t[0] = S;
	    big_addc(t,1,C);
	}
	m[i-s] = t[0];
	t[0] = t[1];
	t[1] = t[2];
	t[2] = 0;	
    }
    return s+1;
}

#endif
