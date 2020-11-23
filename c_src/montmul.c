//
// unsigned montgomery multiplication
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#define NDEBUG
#include <assert.h>

#if (__SIZEOF_POINTER__ == 4) && defined(__SIZEOF_LONG_LONG__) && (__SIZEOF_LONG_LONG__ == 8)
/* Assume 32-bit machine with long long support */
typedef uint64_t   ErlNifBigDoubleDigit;
typedef uint32_t   ErlNifBigDigit;
typedef uint16_t   ErlNifHalfBigDigit;
#define BIG_HAVE_DOUBLE_DIGIT 1

#elif (__SIZEOF_POINTER__ == 4)
/* Assume 32-bit machine with no long support */
#undef  BIG_HAVE_DOUBLE_DIGIT
typedef uint32_t   ErlNifBigDigit;
typedef uint16_t  ErlNifHalfBigDigit;

#elif (__SIZEOF_POINTER__ == 8)
typedef uint64_t ErlNifBigDigit;
typedef uint32_t ErlNifHalfBigDigit;
/* Assume 64-bit machine, does it exist 128 bit long long long ? */
#ifdef __SIZEOF_FLOAT128__
typedef __uint128_t  ErlNifBigDoubleDigit;
#define BIG_HAVE_DOUBLE_DIGIT 1
#else
#undef  BIG_HAVE_DOUBLE_DIGIT
#endif
#else
#error "cannot determine machine size"
#endif

#define BIGNUM_SIZE(arr)  (sizeof((arr))/sizeof(ErlNifBigDigit))

#define MIN(a,b) (((a)<(b)) ? (a) : (b))
#define MAX(a,b) (((a)>(b)) ? (a) : (b))

/* add a and b with carry in + out */
#define DSUMc(a,b,c,s) do {						\
	ErlNifBigDigit ___cr = (c);					\
	ErlNifBigDigit ___xr = (a)+(___cr);				\
	ErlNifBigDigit ___yr = (b);					\
	___cr = (___xr < ___cr);					\
	___xr = ___yr + ___xr;						\
	___cr += (___xr < ___yr);					\
	s = ___xr;							\
	c = ___cr;							\
    }  while(0)

/* add a and b with carry out */
#define DSUM(a,b,c,s) do {						\
	ErlNifBigDigit ___xr = (a);					\
	ErlNifBigDigit ___yr = (b);					\
	___xr = ___yr + ___xr;						\
	s = ___xr;							\
	c = (___xr < ___yr);						\
    }  while(0)

#define DSUBb(a,b,r,d) do {						\
	ErlNifBigDigit ___cr = (r);					\
	ErlNifBigDigit ___xr = (a);					\
	ErlNifBigDigit ___yr = (b)+___cr;				\
	___cr = (___yr < ___cr);					\
	___yr = ___xr - ___yr;						\
	___cr += (___yr > ___xr);					\
	d = ___yr;							\
	r = ___cr;							\
    } while(0)

#define DSUB(a,b,r,d) do {				\
	ErlNifBigDigit ___xr = (a);			\
	ErlNifBigDigit ___yr = (b);			\
	___yr = ___xr - ___yr;				\
	r = (___yr > ___xr);				\
	d = ___yr;					\
    } while(0)

#define DCONST(n) ((ErlNifBigDigit)(n))

#ifdef BIG_HAVE_DOUBLE_DIGIT

#define DLOW(x)        ((ErlNifBigDigit)(x))
#define DHIGH(x)       ((ErlNifBigDigit)(((ErlNifBigDoubleDigit)(x)) >> D_EXP))

#define DLOW2HIGH(x)   (((ErlNifDoubleDigit)(x)) << D_EXP)
#define DDIGIT(a1,a0)  (DLOW2HIGH(a1) + (a0))

#define DMULc(a,b,c,p) do {			       \
        ErlNifBigDoubleDigit _t = ((ErlNifBigDoubleDigit)(a))*(b) + (c); \
	p = DLOW(_t);						\
	c = DHIGH(_t);						\
    } while(0)

#define DMUL(a,b,c1,c0) do {						\
	ErlNifBigDoubleDigit _t = ((ErlNifBigDoubleDigit)(a))*(b);	\
	c0 = DLOW(_t);					\
	c1 = DHIGH(_t);					\
    } while(0)

#else

#define H_EXP (D_EXP >> 1)
#define LO_MASK ((ErlNifBigDigit)((DCONST(1) << H_EXP)-1))
#define HI_MASK ((ErlNifBigDigit)(LO_MASK << H_EXP))

/* Calculate a*b + d1 and store double prec result in d1, d0 */
#define DMULc(a,b,d1,d0) do {					\
	ErlNifHalfBigDigit __a0 = (a);				\
	ErlNifHalfBigDigit __a1 = ((a) >> H_EXP);			\
	ErlNifHalfBigDigit __b0 = (b);				\
	ErlNifHalfBigDigit __b1 = ((b) >> H_EXP);			\
	ErlNifBigDigit __a0b0 = (ErlNifBigDigit)__a0*__b0;		\
	ErlNifBigDigit __a0b1 = (ErlNifBigDigit)__a0*__b1;		\
	ErlNifBigDigit __a1b0 = (ErlNifBigDigit)__a1*__b0;		\
	ErlNifBigDigit __a1b1 = (ErlNifBigDigit)__a1*__b1;		\
	ErlNifBigDigit __p0,__p1,__p2,__c0;				\
	DSUM(__a0b0,d1,__c0,__p0);				\
	DSUM((__c0<<H_EXP),(__p0>>H_EXP),__p2,__p1);		\
	DSUM(__p1,__a0b1,__c0,__p1);				\
	__p2 += __c0;						\
	DSUM(__p1,__a1b0,__c0,__p1);				\
	__p2 += __c0;						\
	DSUM(__p1,__a1b1<<H_EXP,__c0,__p1);			\
	__p2 += __c0;						\
	DSUM(__a1b1, (__p2<<H_EXP),__c0,__p2);			\
	d1 = (__p2 & HI_MASK) | (__p1 >> H_EXP);		\
	d0 = (__p1 << H_EXP) | (__p0 & LO_MASK);		\
    } while(0)

#define DMUL(a,b,d1,d0) do {				\
	ErlNifBigDigit _ds = 0;				\
	DMULc(a,b,_ds,d0);				\
	d1 = _ds;					\
    } while(0)

#endif

#define D_EXP (__SIZEOF_POINTER__*8)
#define D_MASK ((ErlNifBigDigit)(-1))      /* D_BASE-1 */

#if 0
static inline ErlNifBigDigit digit_mul_c(ErlNifBigDigit a, ErlNifBigDigit b,
					 ErlNifBigDigit* carry)
{
    ErlNifBigDigit c = *carry;
    ErlNifBigDigit p;
    DMULc(a,b,c,p);
    *carry = c;
    return p;
}

static inline ErlNifBigDigit digit_mul(ErlNifBigDigit a, ErlNifBigDigit b,
				       ErlNifBigDigit* r0)
{
    ErlNifBigDigit c1,c0;
    DMUL(a,b,c1,c0);
    *r0 = c0;
    return c1;
}

static inline ErlNifBigDigit digit_add_c(ErlNifBigDigit a, ErlNifBigDigit b,
					 ErlNifBigDigit* ci)
{
    ErlNifBigDigit c = *ci;
    ErlNifBigDigit s;
    DSUMc(a,b,c,s);
    *ci = c;
    return s;
}

static inline ErlNifBigDigit digit_sub_b(ErlNifBigDigit a, ErlNifBigDigit b,
					 ErlNifBigDigit* bi)
{
    ErlNifBigDigit r = *bi;
    ErlNifBigDigit d;
    DSUBb(a,b,r,d);
    *bi = r;
    return d;
}
#endif

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
	assert(i < szr);
	// r[i] = digit_add_c(x[i], y[i], &c);
	DSUMc(x[i],y[i],c,r[i]);
	i++;
    }
    while(i < xl) {
	assert(i < szr);
	// r[i] = digit_add_c(x[i], 0, &c);
	DSUMc(x[i],0,c,r[i]);
	i++;
    }
    while(i < yl) {
	assert(i < szr);	
	// r[i] = digit_add_c(0, y[i], &c);
	DSUMc(0,y[i],c,r[i]);	
	i++;
    }
    if (c) {
	assert(i < szr);
	r[i++] = c;
    }
    return i;
}

// x >= y
static int big_sub(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		   ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit b = 0;
    int i = 0;

    while(i < yl) {
	assert(i < szr);
	// r[i] = digit_sub_b(x[i], y[i], &b);
	DSUBb(x[i],y[i],b,r[i]);
	i++;
    }
    while(i < xl) {
	assert(i < szr);
	// r[i] = digit_sub_b(x[i], 0, &b);
	DSUBb(x[i],0,b,r[i]);
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

    big_zero(r, k);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl) && (ij < k); j++, ij++) {
	    ErlNifBigDigit p;
	    // ErlNifBigDigit p = digit_mul_c(x[i],y[j],&cp);
	    DMULc(x[i],y[j],cp,p);
	    // r[ij] = digit_add_c(p,r[ij],&c);
	    DSUMc(p,r[ij],c,r[ij]);
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
    big_zero(r, xl);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl); j++, ij++) {
	    ErlNifBigDigit p;
	    // ErlNifBigDigit p = digit_mul_c(x[i],y[j],&cp);
	    DMULc(x[i],y[j],cp,p);
	    assert(ij < szr);
	    // r[ij] = digit_add_c(p,r[ij],&c);
	    DSUMc(p,r[ij],c,r[ij]);
	}
	assert(ij < szr);
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
    big_zero(r, xl+1);

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
	// z1 = digit_mul(d, d, &b0);
	DMUL(d, d, z1, b0);
	assert(si < szr);
	// r[si] = digit_add_c(r[si], b0, &y_3);
	DSUMc(r[si],b0,y_3,r[si]);
	si++;

	ij = i;
	while(m--) {
	    assert(ij < xl);
	    // b1 = digit_mul(d, x[ij], &b0);
	    DMUL(d, x[ij], b1, b0);
	    ij++;
	    // z0 = digit_add_c(b0,   b0, &y_0);
	    DSUMc(b0, b0, y_0, z0);
	    // z2 = digit_add_c(z0,   z1, &y_2);
	    DSUMc(z0, z1, y_2, z2);
	    assert(si < szr);
	    // r[si] = digit_add_c(r[si],z2, &y_3);
	    DSUMc(r[si],z2,y_3,r[si]);
	    si++;
	    // z1 = digit_add_c(b1, b1, &y_1);
	    DSUMc(b1, b1, y_1, z1);
	}
	z0 = y_0;
	// z2 = digit_add_c(z0, z1, &y_2);
	DSUMc(z0, z1, y_2, z2);
	assert(si < szr);
	// r[si] = digit_add_c(r[si], z2, &y_3);
	DSUMc(r[si], z2, y_3, r[si]);
	if (n != 0) {
	    si++;
	    assert(si < szr);
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
int big_mont_redc(ErlNifBigDigit* t, int tl, ErlNifBigDigit* n, int nl,
		  ErlNifBigDigit* np, int npl, ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit M[nl];
    int ml, rl, tk;

    assert(tl <= 2*nl);
    assert(npl <= nl);

    tk = (tl < nl) ? tl : nl;
    ml = big_mul_k(t, tk, np, npl, M, nl); // M = T*N (mod B^k)
    rl = big_mul(M, ml, n, nl, r, szr);    // R = M*N
    rl = big_add(t, tl, r, rl, r, szr);    // R = T+M*N
    rl = big_shr(r, rl, nl);               // R = R >> k
    if (big_gt(r, rl, n, nl)) // R>N?
	return big_sub(r, rl, n, nl, r, szr);   // R = R - N
    return
	rl;
}

#if NOT_YET
int big_mont_redc(ErlNifBigDigit* t, int tl, ErlNifBigDigit* n, int nl,
		  ErlNifBigDigit* np, int npl, ErlNifBigDigit* r, int szr)
{
    int i;
    ErlNifBigDigit c = 0;
    
    for (i = 0; i < nl; i++) {
	ErlNifBigDigit u = 0;
	ErlNifBigDigit q = t[i]*np0;
	int j;
	for (j = 0; j < nl; j++) {
	    (u,v) = n[j]*q + t[i+j] + u;
	    t[i+j] = v;
	}
	(u,v) = t[i+nl] + u + c;
	t[i+nl] = v;
	c = u;
    }
    for (j = 0; j < nl; j++)
	r[j] = t[j+nl];
    r[nl] = t;
    rl = nl;
    if (big_gt(r, rl, n, nl)) // R>N?
	return big_sub(r, rl, n, nl, r, szr);   // R = R - N
    return nl;
}
#endif


// al,bl < nl < k   R[al+bl]
int big_mont_mul(ErlNifBigDigit* a, int al, ErlNifBigDigit* b, int bl,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl,
		 ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit T[2*nl];
    int tl;
    assert(al <= nl);
    assert(bl <= nl);
    assert(npl <= nl);
    tl = big_mul(a, al, b, bl, T, BIGNUM_SIZE(T));
    return big_mont_redc(T, tl, n, nl, np, npl, r, szr);
}

// al < nl < k
int big_mont_sqr(ErlNifBigDigit* a, int al,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl,
		 ErlNifBigDigit* r, int szr)
{
    ErlNifBigDigit T[2*nl];
    int tl;
    assert(al <= nl);
    assert(npl <= nl);
    tl = big_sqr(a, al, T, BIGNUM_SIZE(T));
    return big_mont_redc(T, tl, n, nl, np, npl, r, szr);
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
int big_mont_pow(ErlNifBigDigit* a, int al,
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
	    pl = big_mont_mul(A[s],al, P[u],pl, n,nl, np,npl,
			      P[v],BIGNUM_SIZE(P[v]));
	    u = v; v = u^1;
	    BIGPRINT("P%d = ", P[u], pl, pos);
	}
	// A' = A^2 (mod R)
	al = big_mont_sqr(A[s], al, n, nl, np, npl, A[t], BIGNUM_SIZE(A[t]));
	s = t; t = s^1;
	BIGPRINT("A%d = ", A[s], al, pos);
    }
    // r = A*P (mod R)
    //printf("al=%d, pl=%d, nl=%d, npl=%d Rsize=%d\r\n", al, pl, nl, npl, 2*nl);
    rl = big_mont_mul(A[s], al, P[u], pl, n, nl, np, npl, R, szR);
    BIGPRINT("R%s = ", R, rl, "");
    return rl;
}
