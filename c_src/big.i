#ifndef __BIG_I__
#define __BIG_I__

#include "big3.i"

#define D_EXP (__SIZEOF_POINTER__*8)
#define D_MASK ((ErlNifBigDigit)(-1))      /* D_BASE-1 */
#define D_SIZE (sizeof(UINT_T)*8)

STATIC void big_zero(UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = 0;
}

STATIC void big_copy(UINT_T* dst, UINT_T* src, int s)
{
    int i;
    for (i = 0; i < s; i++)
	dst[i] = src[i];
}

STATIC int big_bits(UINT_T* x, int xl)
{
    int n = 8*sizeof(UINT_T)*(xl-1);
    UINT_T h = x[xl-1];
    if ((n == 0) && (h == 0))
	return 1;
    while(h) {
	n++;
	h >>= 1;
    }
    return n;
}

STATIC int big_test(UINT_T* x, int xl, unsigned pos)
{
    int d = pos / D_EXP; // digit
    pos %= D_EXP;      // bit
    if (d >= xl) return 0; // definied as zero
    return (x[d] & (1 << pos)) != 0;
}

STATIC int big_comp(UINT_T* x, int xl, UINT_T* y, int yl)
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

STATIC int big_gt(UINT_T* x, int xl, UINT_T* y, int yl)
{
    return big_comp(x, xl, y, yl) > 0;
}

STATIC int big_addc(UINT_T* t, int i, int n, UINT_T c)
{
    while(c && (i < n)) {
	ADDC(t[i],0,c,&c,&t[i]);
	i++;
    }
    return i;
}

STATIC int big_add(UINT_T* x, int xl, UINT_T* y, int yl,
		   UINT_T* r, int szr)
{
    UINT_T c = 0;
    int i = 0;

    while((i < xl) && (i < yl)) {
	ADDC(x[i],y[i],c,&c,&r[i]);
	i++;
    }
    while(i < xl) {
	ADDC(x[i],0,c,&c,&r[i]);
	i++;
    }
    while(i < yl) {
	ADDC(0,y[i],c,&c,&r[i]);	
	i++;
    }
    if (c)
	r[i++] = c;
    return i;
}

// x >= y
STATIC int big_sub(UINT_T* x, int xl, UINT_T* y, int yl,
		   UINT_T* r, int szr)
{
    UINT_T b = 0;
    int i = 0;

    while(i < yl) {
	SUBB(x[i],y[i],b,&b,&r[i]);
	i++;
    }
    while(i < xl) {
	SUBB(x[i],0,b,&b,&r[i]);
	i++;
    }
    do {
	i--;
    } while((i>0) && (r[i] == 0));
    return i+1;
}

// subtract and return borrow
STATIC UINT_T big_subb(UINT_T* x, UINT_T* y, UINT_T* r, int s)
{
    UINT_T B = 0;
    int i;
    for (i = 0; i < s; i++)
	SUBB(x[i],y[i],B,&B,&r[i]);
    SUBB(x[s],0,B,&B,&r[s]);
    return B;
}

// trim MSB=0
STATIC INLINE int big_trim(UINT_T* x, int xl)
{
    while((xl>1) && (x[xl-1]==0))
	xl--;
    return xl;
}

//
// if (X > N) X = X - N
//  
STATIC INLINE int big_norm0(UINT_T* x, int xl, UINT_T* n, int nl)
{
    if (big_comp(x, xl, n, nl) > 0)
	return big_sub(x, xl, n, nl, x, xl);
    else
	return xl;
}

// version of above but wihtout comparsion x must be of size xl+1!
#define NSHIFT (bit_sizeof(UINT_T)-1)
STATIC INLINE int big_norm1(UINT_T* z, UINT_T* n, int s)
{
    UINT_T mask = ((-(INT_T)z[s]) >> NSHIFT);
    int i;
    UINT_T B;

    SUB(z[0], (n[0] & mask), &B, &z[0]);
    for (i = 1; i < s; i++) {
	SUBB(z[i], (n[i] & mask), B, &B, &z[i]);
    }
    return s;
}

// inline shift x:xl k digits to the right
STATIC INLINE int big_shr(UINT_T* x, int xl, int k)
{
    int n = xl - k;
    int i;
    for (i = 0; i < n; i++)
	x[i] = x[i+k];
    return n;
}

// multiply x and y module (B^k) (B is 2^w) (w = 8*sizeof(UINT_T))
STATIC INLINE int big_mul_k(UINT_T* x, int xl,
			    UINT_T* y, int yl,
			    UINT_T* r, int k)
{
    int i;
    for (i = 0; i < xl; i++) {
	UINT_T cp = 0;
	UINT_T c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl) && (ij < k); j++, ij++) {
	    UINT_T p;
	    MULA(x[i],y[j],cp,&cp,&p);
	    ADDC(p,r[ij],c,&c,&r[ij]);
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

// multiply x[n] and y[n] produce result r[2n]
// the result of x[n]*y[n] is added to r[2n] so it must be
// set to zero before calling big_n_mul, unless intensional
//
STATIC INLINE void big_n_mul(UINT_T* x, UINT_T* y, UINT_T* r, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	UINT_T c = 0;
	int j, ij=i;
	for (j = 0; j < n; j++) {
	    MULAB(x[i],y[j],r[ij],c,&c,&r[ij]);
	    ij++;
	}
	r[ij] = c;  // top index is 2n-1
    }
}

STATIC INLINE int big_mul(UINT_T* x, int xl,
			  UINT_T* y, int yl,
			  UINT_T* r, int szr)
{
    int i;
    for (i = 0; i < xl; i++) {
	UINT_T cp = 0;
	UINT_T c = 0;
	int j, ij;
	for (j = 0, ij=i; (j < yl); j++, ij++) {
	    UINT_T p;
	    MULA(x[i],y[j],cp,&cp,&p);
	    ADDC(p,r[ij],c,&c,&r[ij]);
	}
	r[ij] = c + cp;
    }
    i = xl+yl-1;
    if (r[i] == 0)
	return i;
    else
	return i+1;
}

// x[s]  r[2n]
STATIC INLINE int big_n_sqr(UINT_T* x, UINT_T* r, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	DECL3(c);
	UINT_T Co;
	int j;
	ZERO3(c);
	SQRA(x[i], r[i+i], &Co, &r[i+i]);
	ELEM3(c,0) = Co;
	for (j = i+1; j < n; j++) {
	    UINT_T b1,b0;
	    MUL(x[i],x[j],&b1,&b0);   // (b1,b0) = xi*xj
	    ADD32p(c, b1,b0);         // (c2,c1,c0) += (b1,b0)
	    ADD32p(c, b1,b0);         // (c2,c1,c0) += (b1,b0)
	    ADD31p(c, r[i+j]);        // (c2,c1,c0) += r[i+j]
	    r[i+j] = ELEM3(c,0);
	    SHR3(c);
	}
	ADD(r[i+n], ELEM3(c,0), &Co, &r[i+n]);
	ADDC(r[i+n+1], ELEM3(c,1), Co, &Co, &r[i+n+1]);
    }
    return n+n;
}

STATIC INLINE int big_sqr(UINT_T* x, int xl,
			  UINT_T* r, int szr)
{
    UINT_T d;
    int i, n;
    int ri, si;

    ri = si = i = 0;
    n = xl;
    while(n--) {
	UINT_T y_0 = 0, y_1 = 0, y_2 = 0, y_3 = 0;
	UINT_T b0, b1;
	UINT_T z0, z1, z2;
	int m = n;
	int ij;

	si = ri;
	
	d = x[i++];
	SQR(d, &z1, &b0);
	ADDC(r[si],b0,y_3,&y_3,&r[si]);
	si++;

	ij = i;
	while(m--) {
	    MUL(d, x[ij], &b1, &b0);
	    ij++;
	    ADDC(b0, b0, y_0, &y_0, &z0);
	    ADDC(z0, z1, y_2, &y_2, &z2);
	    ADDC(r[si],z2,y_3,&y_3,&r[si]);
	    si++;
	    ADDC(b1, b1, y_1, &y_1, &z1);
	}
	z0 = y_0;
	ADDC(z0, z1, y_2, &y_2, &z2);
	ADDC(r[si], z2, y_3, &y_3, &r[si]);
	if (n != 0) {
	    si++;
	    r[si] = (y_1+y_2+y_3);
	    ri += 2;
	}
    }
    if (r[si] == 0)
	return si;
    else
	return si + 1;    
}

#endif
