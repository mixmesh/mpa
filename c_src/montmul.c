//
// unsigned montgomery multiplication
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define NDEBUG
#include <assert.h>

#include "montmul.h"

#define UINT_T  ErlNifBigDigit
#define UINTH_T ErlNifBigHalfDigit
#ifdef BIG_HAVE_DOUBLE_DIGIT
#define UINTD_T ErlNifBigDoubleDigit
#endif

#include "unroll.i"

#define D_EXP (__SIZEOF_POINTER__*8)
#define D_MASK ((ErlNifBigDigit)(-1))      /* D_BASE-1 */


// print hi->lo
void big_print(ErlNifBigDigit* x, int xl)
{
    int i;
    printf("[%d]{%lu",xl,(unsigned long)x[xl-1]);
    for (i = xl-2; i >= 0; i--) printf(",%lu", (unsigned long)x[i]);
    printf("}");
}

static void big_copy(ErlNifBigDigit* dst, ErlNifBigDigit* src, int n)
{
    memcpy(dst, src, n*sizeof(ErlNifBigDigit));
}

static void big_zero(ErlNifBigDigit* dst, int n)
{
    memset(dst, 0, n*sizeof(ErlNifBigDigit));
}

static int big_bits(ErlNifBigDigit* x, int xl)
{
    int n = 8*sizeof(ErlNifBigDigit)*(xl-1);
    ErlNifBigDigit h = x[xl-1];
    if ((n == 0) && (h == 0))
	return 1;
    while(h) {
	n++;
	h >>= 1;
    }
    return n;
}

static int big_bit_test(ErlNifBigDigit* x, int xl, unsigned pos)
{
    int d = pos / D_EXP; // digit
    pos %= D_EXP;      // bit
    if (d >= xl) return 0; // definied as zero
    return (x[d] & (1 << pos)) != 0;
}


#if 0
static int big_small(ErlNifBigDigit* x, ErlNifBigDigit d)
{
    x[0] = d;
    return 1;
}
#endif

static int big_comp(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl)
{
    if (xl < yl)
	return -1;
    else if (xl > yl)
	return 1;
    else {
	if (x == y)
	    return 0;
	x += (xl-1);
	y += (yl-1);
	while((xl > 0) && (*x == *y)) {
	    x--;
	    y--;
	    xl--;
	}
	if (xl == 0)
	    return 0;
	return (*x < *y) ? -1 : 1;
    }
}

static int big_gt(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl)
{
    return big_comp(x, xl, y, yl) > 0;
}

static int big_add(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		   ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit c = 0;
    int i = 0;

    while((i < xl) && (i < yl)) {
	addc(x[i],y[i],c,&c,&r[i]);
	i++;
    }
    while(i < xl) {
	addc(x[i],0,c,&c,&r[i]);
	i++;
    }
    while(i < yl) {
	addc(0,y[i],c,&c,&r[i]);	
	i++;
    }
    if (c)
	r[i++] = c;
    return i;
}

// x >= y
static int big_sub(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		   ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit b = 0;
    int i = 0;

    while(i < yl) {
	subb(x[i],y[i],b,&b,&r[i]);
	i++;
    }
    while(i < xl) {
	subb(x[i],0,b,&b,&r[i]);
	i++;
    }
    do {
	i--;
    } while((i>0) && (r[i] == 0));
    return i+1;
}

// inline shift x:xl k digits to the right
static int big_shr(ErlNifBigDigit* x, int xl, int k)
{
    int n = xl - k;
    int i;
    for (i = 0; i < n; i++)
	x[i] = x[i+k];
    return n;
}

// multiply x and y module (B^k) (B is 2^w) (w = machine word size)
static int inline big_mul_k(ErlNifBigDigit* x, int xl,
			    ErlNifBigDigit* y, int yl,
			    ErlNifBigDigit* r, int k)
{
    int i;

    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl) && (ij < k); j++, ij++) {
	    ErlNifBigDigit p;
	    mula(x[i],y[j],cp,&cp,&p);
	    addc(p,r[ij],c,&c,&r[ij]);
	}
	if (ij < k)
	    r[ij] = c + cp;
    }
    i = ((xl+yl-1) < k) ? (xl+yl-1) : k-1;
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

static int inline big_mul(ErlNifBigDigit* x, int xl,
			  ErlNifBigDigit* y, int yl,
			  ErlNifBigDigit* r, int szr)
{
    int i;

    assert(xl+yl <= szr);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl); j++, ij++) {
	    ErlNifBigDigit p;
	    mula(x[i],y[j],cp,&cp,&p);
	    addc(p,r[ij],c,&c,&r[ij]);
	}
	r[ij] = c + cp;
    }
    i = xl+yl-1;
    assert(i < szr);    
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

static int inline big_sqr(ErlNifBigDigit* x, int xl,
			  ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit d;
    int i, n;
    int ri, si;

    assert(xl+xl <= szr);

    ri = si = i = 0;
    n = xl;
    while(n--) {
	ErlNifBigDigit y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0;
	ErlNifBigDigit b0, b1;
	ErlNifBigDigit z0, z1, z2;
	int m = n;
	int ij;

	si = ri;
	
	d = x[i++];
	sqr(d, &z1, &b0);
	// mul(d, d, &z1, &b0);
	addc(r[si],b0,y_3,&y_3,&r[si]);
	si++;

	ij = i;
	while(m--) {
	    mul(d, x[ij], &b1, &b0);
	    ij++;
	    addc(b0, b0, y_0, &y_0, &z0);
	    addc(z0, z1, y_2, &y_2, &z2);
	    addc(r[si],z2,y_3,&y_3,&r[si]);
	    si++;
	    addc(b1, b1, y_1, &y_1, &z1);
	}
	z0 = y_0;
	addc(z0, z1, y_2, &y_2, &z2);
	assert(si < szr);
	addc(r[si], z2, y_3, &y_3, &r[si]);
	if (n != 0) {
	    si++;
	    r[si] = (y_1+y_2+y_3);
	    ri += 2;
	}
	else {
	    assert((y_1+y_2+y_3) == 0);
	}
    }
    assert(si < szr);    
    if (r[si] == 0)
	return si;
    else
	return si + 1;    
}

// r must be of size 2nl
int big_mont_redc_default(ErlNifBigDigit* T, int tl,
			  ErlNifBigDigit* n, int nl,
			  ErlNifBigDigit* np, int npl,
			  ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit M[nl];
    int ml, rl, tk;

    assert(tl <= 2*nl);
    assert(npl <= nl);

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

// note that P is destructivly updated
int big_mont_redc_sos(ErlNifBigDigit* P, int pl,
		      ErlNifBigDigit* n, int s,
		      ErlNifBigDigit* np, int npl,
		      ErlNifBigDigit* Z, int szZ)
{
    int i, j;
    int Zl;
    ErlNifBigDigit t = 0;
    
    for (i = 0; i < s; i++) {
	ErlNifBigDigit u = 0;
	ErlNifBigDigit v;
	ErlNifBigDigit q;

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

static void inline add31(ErlNifBigDigit a2,ErlNifBigDigit a1,ErlNifBigDigit a0,
			 ErlNifBigDigit b0,
			 ErlNifBigDigit* d2,ErlNifBigDigit* d1, ErlNifBigDigit* d0)
{
    ErlNifBigDigit c;
    add(a0, b0, &c, d0);
    add(a1, c,  &c, d1);
    add(a2, c,  &c, d2);
}

static void inline add32(ErlNifBigDigit a2,ErlNifBigDigit a1,ErlNifBigDigit a0,
			 ErlNifBigDigit b1,ErlNifBigDigit b0,
			 ErlNifBigDigit* d2,ErlNifBigDigit* d1,ErlNifBigDigit* d0)
{
    ErlNifBigDigit c;
    add(a0, b0, &c,d0);
    addc(a1,b1,c,&c,d1);
    add(a2,c,&c,d2);
}

int big_mont_redc_sps(ErlNifBigDigit* P, int pl,
		      ErlNifBigDigit* n, int s,
		      ErlNifBigDigit* np, int npl,
		      ErlNifBigDigit* Z, int szZ)
{
    ErlNifBigDigit t=0, u=0, v=0;
    ErlNifBigDigit p0,p1;
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

int big_mont_redc(redc_type_t redc_type,
		  ErlNifBigDigit* t, int tl, ErlNifBigDigit* n, int nl,
		  ErlNifBigDigit* np, int npl, ErlNifBigDigit* r, int szr)
{
    switch(redc_type) {
    case REDC_DEFAULT:
	return big_mont_redc_default(t, tl, n, nl, np, npl, r, szr);
    case REDC_SOS:
	return big_mont_redc_sos(t, tl, n, nl, np, npl, r, szr);
    case REDC_SPS:
	return big_mont_redc_sps(t, tl, n, nl, np, npl, r, szr);
    default:
	return -1;
    }
}

// al,bl < nl < k   R[al+bl]
int big_mont_mul(redc_type_t redc_type,
		 ErlNifBigDigit* a, int al, ErlNifBigDigit* b, int bl,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl,
		 ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit T[2*nl];
    // int tl;
    assert(al <= nl);
    assert(bl <= nl);
    assert(npl <= nl);
    // fixme: the product must be 2*nl number of digits!)
    big_zero(T, 2*nl);
    big_mul(a, al, b, bl, T, BIGNUM_SIZE(T));
    return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
}

// al < nl < k
int big_mont_sqr(redc_type_t redc_type,
		 ErlNifBigDigit* a, int al,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl,
		 ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit T[2*nl];
    // int tl;
    assert(al <= nl);
    assert(npl <= nl);
    big_zero(T, 2*nl);
    big_sqr(a, al, T, BIGNUM_SIZE(T));
    return big_mont_redc(redc_type, T, 2*nl, n, nl, np, npl, r, szr);
}

#define BIGPRINT1(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)

#if 0
#define BIGPRINT(fmt, x, xl, args...) do { 			\
    printf(fmt, args); big_print((x),(xl)); printf("\r\n");	\
    } while(0)
#endif
#define BIGPRINT(fmt, x, xl, args...)

// a^e (mod R)  (R=B^k > N)  (B = ErlNifBigDigit base) r[2*al]
int big_mont_pow(redc_type_t redc_type,
		 ErlNifBigDigit* a, int al,
		 ErlNifBigDigit* e, int el,
		 ErlNifBigDigit* p, int pl,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl,
		 ErlNifBigDigit* R, int szR)
{
    ErlNifBigDigit P[2][2*nl+1];
    ErlNifBigDigit A[2][2*nl+1];
    int u, v;
    int s, t;
    int pos, nbits;
    int rl;

    assert(al <= nl);
    assert(npl <= nl);
    assert(pl <= nl);
    
    u = 0; v = u^1;
    big_copy(P[u], p, pl);   // mont(1) !

    BIGPRINT("P%s=", P[u], pl, "");

    s = 0; t = s^1;
    big_copy(A[s], a, al);  // check al!

    BIGPRINT("A%s=", A[s], al, "");
    BIGPRINT("E%s=", e, el, "");

    nbits = big_bits(e, el)-1;
    // printf("E #bits = %d\r\n", nbits);
    for (pos = 0; pos < nbits; pos++) {
	int bit = big_bit_test(e, el, pos);
	// printf("E[%d] = %d\r\n", pos, bit);
	if (bit) {
	    // P' = A*P (mod R)
	    pl = big_mont_mul(redc_type,
			      A[s],al, P[u],pl, n,nl, np,npl,
			      P[v],BIGNUM_SIZE(P[v]));
	    u = v; v = u^1;
	    BIGPRINT("P%d = ", P[u], pl, pos);
	}
	// A' = A^2 (mod R)
	al = big_mont_sqr(redc_type, A[s], al,
			  n, nl, np, npl, A[t], BIGNUM_SIZE(A[t]));
	s = t; t = s^1;
	BIGPRINT("A%d = ", A[s], al, pos);
    }
    // r = A*P (mod R)
    //printf("al=%d, pl=%d, nl=%d, npl=%d Rsize=%d\r\n", al, pl, nl, npl, 2*nl);
    rl = big_mont_mul(redc_type, A[s], al, P[u], pl, n, nl, np, npl, R, szR);
    BIGPRINT("R%s = ", R, rl, "");
    return rl;
}
