//
// unsigned montgomery multiplication
//
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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
/* Assume 64-bit machine, does it exist 128 bit long long long ? */
#undef  BIG_HAVE_DOUBLE_DIGIT
typedef uint64_t ErlNifBigDigit;
typedef uint32_t ErlNifHalfBigDigit;
#else
#error "cannot determine machine size"
#endif

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
#define DMUL(a,b,c1,c0) do { \
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


static inline ErlNifBigDigit digit_mul_c(ErlNifBigDigit a, ErlNifBigDigit b,
					 ErlNifBigDigit* carry)
{
    ErlNifBigDigit c = *carry;
    ErlNifBigDigit p;
    DMULc(a,b,c,p);
    *carry = c;
    return p;
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

// print hi->lo
void big_print(ErlNifBigDigit* x, int xl)
{
    int i;
    printf("{%lu",(unsigned long)x[xl-1]);
    for (i = xl-2; i >= 0; i--) printf(",%lu", (unsigned long)x[i]);
    printf("}\n");
}

static int big_bit_test(ErlNifBigDigit* x, int xl, unsigned pos)
{
    int d = pos / D_EXP; // digit
    pos %= D_EXP;      // bit
    if (d >= xl) return 0; // definied as zero
    return (x[d] & (1 << pos)) != 0;
}

static int big_small(ErlNifBigDigit* x, ErlNifBigDigit d)
{
    x[0] = d;
    return 1;
}

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
		   ErlNifBigDigit* r)
{
    ErlNifBigDigit c = 0;
    int i = 0;

    while((i < xl) && (i < yl)) {
	r[i] = digit_add_c(x[i], y[i], &c);
	i++;
    }
    while(i < xl) {
	r[i] = digit_add_c(x[i], 0, &c);
	i++;
    }
    while(i < yl) {
	r[i] = digit_add_c(0, y[i], &c);
	i++;
    }
    if (c)
	r[i++] = c;
    return i;
}

// x > y
static int big_sub(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		   ErlNifBigDigit* r)
{
    ErlNifBigDigit b = 0;
    int i = 0;

    while(i < yl) {
	r[i] = digit_sub_b(x[i], y[i], &b);
	i++;
    }
    while(i < xl) {
	r[i] = digit_sub_b(x[i], 0, &b);
	i++;
    }
    do {
	i--;
    } while(i && (r[i] == 0));
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

// regular multiply, multiply upto digit k
// multiply x[xl] and y[yl] and return in r[xl+ýl] (return actual size)
static int big_mul_k(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		     ErlNifBigDigit* r, int k)
{
    int i;

    big_zero(r, k);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j;
	for (j = 0; (j < yl) && ((i+j) < k); j++) {
	    ErlNifBigDigit p = digit_mul_c(x[i],y[j],&cp);
	    r[i+j] = digit_add_c(p,r[i+j],&c);
	}
	if ((i+j) < k)
	    r[i+j] = c + cp;
    }
    i = ((xl+yl-1) < k) ? (xl+yl-1) : k-1;
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

static int big_mul(ErlNifBigDigit* x, int xl, ErlNifBigDigit* y, int yl,
		   ErlNifBigDigit* r)
{
    int i;

    big_zero(r, xl+yl);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j;
	for (j = 0; (j < yl); j++) {
	    ErlNifBigDigit p = digit_mul_c(x[i],y[j],&cp);
	    r[i+j] = digit_add_c(p,r[i+j],&c);
	}
	r[i+j] = c + cp;
    }
    i = xl+yl-1;
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

static int big_sqr(ErlNifBigDigit* x, int xl, ErlNifBigDigit* r)
{
    int i;

    big_zero(r, xl+xl);
    
    for (i = 0; i < xl; i++) {
	ErlNifBigDigit cp = 0;
	ErlNifBigDigit c = 0;
	int j;
	for (j = 0; (j < xl); j++) {
	    ErlNifBigDigit p = digit_mul_c(x[i],x[j],&cp);
	    r[i+j] = digit_add_c(p,r[i+j],&c);
	}
	r[i+j] = c + cp;
    }
    i = xl+xl-1;
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

int big_mont_redc(ErlNifBigDigit* t, int tl, ErlNifBigDigit* n, int nl,
		  ErlNifBigDigit* np, int npl, ErlNifBigDigit* r)
{
    ErlNifBigDigit M[nl];
    int ml, rl, tk;

    assert(tl <= 2*nl);
    assert(npl <= nl);

    tk = (tl < nl) ? tl : nl;
    ml = big_mul_k(t, tk, np, npl, M, nl);  // M = T*N (mod B^k)
    rl = big_mul(M, ml, n, nl, r);         // R = M*N
    rl = big_add(t, tl, r, rl, r);         // R = T+M*N
    rl = big_shr(r, rl, nl);                // R = R >> k
    if (big_gt(r, rl, n, nl))
	return big_sub(r, rl, n, nl, r);
    return
	rl;
}

// al,bl < nl < k   R[al+bl]
int big_mont_mul(ErlNifBigDigit* a, int al, ErlNifBigDigit* b, int bl,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl, ErlNifBigDigit* r)
{
    ErlNifBigDigit T[2*nl];
    int tl;
    assert(al <= nl);
    assert(bl <= nl);
    assert(npl <= nl);
    tl = big_mul(a, al, b, bl, T);
    return big_mont_redc(T, tl, n, nl, np, npl, r);
}

// al < nl < k
int big_mont_sqr(ErlNifBigDigit* a, int al,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl, ErlNifBigDigit* r)
{
    ErlNifBigDigit T[2*nl];
    int tl;
    assert(al <= nl);
    assert(npl <= nl);
    tl = big_sqr(a, al, T);
    return big_mont_redc(T, tl, n, nl, np, npl, r);
}

// a^e (mod R)  (R=B^k > N)  (B = ErlNifBigDigit base) r[2*al]
int big_mont_pow(ErlNifBigDigit* a, int al, ErlNifBigDigit* e, int el,
		 ErlNifBigDigit* n, int nl,
		 ErlNifBigDigit* np, int npl, ErlNifBigDigit* r)
{
    ErlNifBigDigit P[2][2*nl];
    ErlNifBigDigit A[2][2*nl];
    int u, v;
    int s, t;
    int pl;
    int pos, nbits;

    assert(al <= nl);
    assert(npl <= nl);
    
    u = 0; v = u^1;
    pl = big_small(P[u], 1);  // P = 1

    s = 0; t = s^1;
    big_copy(A[s], a, al);  // check al!

    nbits = big_bits(e, el);
    for (pos = 0; pos < nbits; pos++) {
	if (big_bit_test(e, el, pos)) {
	    // P' = A*P (mod R)
	    pl = big_mont_mul(A[s], al, P[u], pl, n, nl, np, npl, P[v]);
	    u = v; v = u^1;
	}
	// A' = A^2 (mod R)
	al = big_mont_sqr(A[s], al, n, nl, np, npl, A[t]);
	s = t; t = s^1;
    }
    // r = A*P (mod R)
    return big_mont_mul(A[s], al, P[u], pl, n, nl, np, npl, r);
}

#ifdef TEST
int main(int argc, char** argv)
{
    // test
    ErlNifBigDigit a[4];
    ErlNifBigDigit b[4];
    ErlNifBigDigit r[8];
    int i;
    int n;
    
    for (i = 0; i < 4; i++)
	a[i] = b[i] = 0x9;
    big_zero(r, 8);
    
    n = big_mul(a, 4, b, 4, r);

    printf("a:4 = "); big_print(a, 4);
    printf("b:4 = "); big_print(b, 4);
    printf("r:%d = ", n); big_print(r, n);
    exit(0);
}
#endif
