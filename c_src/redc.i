#ifndef __REDC_I__
#define __REDC_I__

// r must be of size 2nl
int big_mont_redc_default(UINT_T* T, int tl,
			  UINT_T* n, int nl,
			  UINT_T* np, int npl,
			  UINT_T* r, int szr)
{
    UINT_T M[nl];
    int ml, rl, tk;

    tk = (tl < nl) ? tl : nl;
    big_zero(M, nl);
    ml = big_mul_k(T, tk, np, npl, M, nl); // M = T*N (mod B^k)
    big_zero(r, ml);
    rl = big_mul(M, ml, n, nl, r, szr);    // r = M*n
    rl = big_add(T, tl, r, rl, r, szr);    // r = t+(M*N)
    rl = big_shr(r, rl, nl);               // r = r >> k
    if (big_gt(r, rl, n, nl)) // r>N?
	return big_sub(r, rl, n, nl, r, szr);   // r = r - N
    return
	rl;
}

#endif
