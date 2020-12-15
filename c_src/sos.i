// Separated Operand Scanning
#ifndef __SOS_I
#define __SOS_I

// note that P is destructivly updated
static int big_mont_redc_sos(UINT_T* P, UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    int i, j;
    UINT_T t = 0;
    
    for (i = 0; i < s; i++) {
	UINT_T C = 0;
	UINT_T m;

	mul0(P[i],np[0],&m);

	for (j = 0; j < s; j++) {
	    mulab(n[j],m,P[i+j],C,&C,&P[i+j]);
	}
	addc(P[i+s],t,C,&C,&P[i+s]);
	t = C;
    }
    for (j = 0; j < s; j++)
	r[j] = P[j+s];
    r[s] = t;
    return s+1;
}

static int big_mont_mul_sos(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    int i;    
    UINT_T P[2*s];

    big_zero(P, BIGNUM_SIZE(P));
    for (i = 0; i < s; i++) {
	UINT_T cp = 0;
	UINT_T c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < s); j++, ij++) {
	    UINT_T p0;
	    mula(a[i],b[j],cp,&cp,&p0);
	    addc(p0,P[ij],c,&c,&P[ij]);
	}
	P[ij] = c + cp;
    }    
    return big_mont_redc_sos(P, n, np, r, s);
}

#endif
