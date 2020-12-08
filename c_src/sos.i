#ifndef __SOS_I
#define __SOS_I

// note that P is destructivly updated
int big_mont_redc_sos(UINT_T* P, int pl,
		      UINT_T* n, int s,
		      UINT_T* np, int npl,
		      UINT_T* Z, int szZ)
{
    int i, j;
    int Zl;
    UINT_T t = 0;
    
    for (i = 0; i < s; i++) {
	UINT_T u = 0;
	UINT_T v;
	UINT_T q;

	mul0(P[i],np[0],&q);

	for (j = 0; j < s; j++) {
	    mulab(n[j],q,P[i+j],u,&u,&v);
	    P[i+j] = v;
	}
	addc(P[i+s],t,u,&u,&P[i+s]);
	t = u;
    }
    for (j = 0; j < s; j++)
	Z[j] = P[j+s];
    Z[s] = t;
    i = s;
    while(i && (Z[i]==0)) i--;
    Zl = i+1;
    if (big_gt(Z, Zl, n, s)) // R>N?
	return big_sub(Z, Zl, n, s, Z, Zl);   // R = R - N
    return Zl;
}

#endif
