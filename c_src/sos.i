// Separated Operand Scanning
#ifndef __SOS_I
#define __SOS_I

// note that P is destructivly updated
int big_mont_redc_sos(UINT_T* P, int pl,
		      UINT_T* n, int s,
		      UINT_T* np, int npl,
		      UINT_T* Z, int szZ)
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
	Z[j] = P[j+s];
    Z[s] = t;
    i = s;
    while(i && (Z[i]==0)) i--;
    return big_norm0(Z, i+1, n, s);
}

#endif
