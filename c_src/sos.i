// Separated Operand Scanning
#ifndef __SOS_I
#define __SOS_I

// note that P is destructivly updated (but kept as local array)
static int big_mont_redc_sos(UINT_T* P, UINT_T* n, UINT_T* np, UINT_T* r, int s)
{
    int i, j;
    UINT_T t = 0;
    
    for (i = 0; i < s; i++) {
	UINT_T C = 0;
	UINT_T m;
	int ij = i;

	MUL0(P[i],np[0],&m);

	for (j = 0; j < s; j++) {
	    MULAB(n[j],m,P[ij],C,&C,&P[ij]);
	    ij++;
	}
	ADDC(P[ij],t,C,&C,&P[ij]);
	t = C;
    }
    for (j = 0; j < s; j++)
	r[j] = P[j+s];
    r[s] = t;
    return s+1;
}

// a[s], b[s], r[s+1]
static int big_mont_mul_sos(UINT_T* a, UINT_T* b, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    UINT_T P[2*s];
    big_zero(P, BIGNUM_SIZE(P));
    big_n_mul(a, b, P, s);
    return big_mont_redc_sos(P, n, np, r, s);
}

static int big_mont_sqr_sos(UINT_T* a, UINT_T* np, UINT_T* n,
			    UINT_T* r,int s)
{
    UINT_T P[2*s];
    
    big_zero(P, BIGNUM_SIZE(P));
    big_n_sqr(a, P, s);
    // big_n_mul(a,a,P,s);
    BIGPRINT0("%sa^2,P=", P, 2*s, "");
    return big_mont_redc_sos(P, n, np, r, s);    
}


#endif
