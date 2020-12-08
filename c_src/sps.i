#ifndef __SPS_I__
#define __SPS_I__


int big_mont_redc_sps(UINT_T* P, int pl,
		      UINT_T* n, int s,
		      UINT_T* np, int npl,
		      UINT_T* Z, int szZ)
{
    UINT_T t=0, u=0, v=0;
    UINT_T p0,p1;
    int i, j, Zl;
    
    for (i = 0; i < s; i++) {
	for (j = 0; j < i-1; j++) {
	    // (t,u,v) = (t,u,v) + Z[j]*n[i-j]
	    mul(Z[j],n[i-j],&p0,&p1);
	    add32(t,u,v,p1,p0,&t,&u,&v);
	}
	// (t,y,v) = (t,u,v) + P[i]	
	add31(t,u,v,P[i],&t,&u,&v);
	// Z[i] = v*m0' mod 2^w
	mul0(v, np[0], &Z[i]);
	// (t,u,v) = (t,u,v)+Z[i]*n[0]
	mul(Z[i],n[0],&p0,&p1);
	add32(t,u,v,p1,p0,&t,&u,&v);
	v=u; u=t; t=0;
    }
    for (i = s; i < 2*(s-1); i++) {
	for (j = i-s+1; j < s-1; j++) {
	    // (t,u,v) = (t,u,v) + Z[j]*n[i-j]
	    mul(Z[j],n[i-j],&p0,&p1);
	    add32(t,u,v,p1,p0,&t,&u,&v);
	}
	// (t,y,v) = (t,u,v) + P[i]
	add31(t,u,v,P[i],&t,&u,&v);
	Z[i-s] = v;
	v=u; u=t; t=0;
    }
    // (t,u,v) = (t,u,v) + P[2s-1]
    add31(t,u,v,P[2*s-1],&t,&u,&v);
    Z[s-1] = v;
    Z[s] = u;
    i = s;
    while(i && (Z[i]==0)) i--;
    Zl = i+1;    
    if (big_gt(Z, Zl, n, s)) // R>N?
	return big_sub(Z, Zl, n, s, Z, Zl);   // R = R - N
    return Zl;    
}

#endif
