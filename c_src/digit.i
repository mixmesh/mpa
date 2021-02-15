#ifndef __DIGIT_I__
#define __DIGIT_I__

#ifndef UINT_T
#error "UINT_T not defined"
#endif
#ifndef UINTH_T
#error "UINTH_T not defined"
#endif
#ifndef UINTD_T
#warning "UINTD_T not defined"
#endif

#ifdef __OPENCL_VERSION__
#define CARRY(x,y) ((x) < (y))
#else
#define CARRY(x,y) ((x) < (y))
#endif

#define bit_sizeof(T) (sizeof(T)*8)

// ADD0(X,Y,S)   S = X + Y

#define ADD0(X, Y, S) do {	\
	*(S) = (X) + (Y);	\
    } while(0)

// ADD(X,Y,C,S)   (C,S) = X + Y

#define ADD(X, Y, C, S) do {			\
	UINT_T Xl = (X);			\
	Xl += (Y);				\
	*(S) = Xl;				\
	*(C) = CARRY(Xl,(Y));			\
    } while(0)

// ADDC(X,Y,Ci,C,S)  (C,S) = X + Y + Ci

#define ADDC(X, Y, Ci, C, S) do {		\
	UINT_T Xl = (X);			\
	UINT_T CIl = (Ci);			\
	Xl += CIl;				\
	CIl = CARRY(Xl,CIl);			\
	Xl += (Y);				\
	*(S) = Xl;				\
	*(C) = CIl + CARRY(Xl,(Y));		\
    } while(0)

// SUB0(X,Y,D)  D = X - Y

#define SUB0(X, Y, D) do {    \
	*(D) = (X) - (Y);     \
    } while(0)


// SUB(X,Y,B,D) (B,D) = X - Y

#define SUB(X, Y, B, D) do {	       \
	UINT_T Yl = (Y);	       \
	Yl = (X) - Yl;		       \
	*(D) = Yl;		       \
	*(B) = CARRY((X),Yl);	       \
    } while(0)

// SUB(X,Y,Bi,B,D) (B,D) = X - Y - Bi

#define SUBB(X, Y, Bi, B, D) do {		\
	UINT_T Yl = (Y);			\
	UINT_T BIl = (Bi);			\
	Yl += (Bi);				\
	BIl = CARRY(Yl,BIl);			\
	Yl = (X) - Yl;				\
	*(B) = BIl + CARRY((X),Yl);		\
	*(D) = Yl;				\
    } while(0)

// MUL0(X,Y,P)  P = X*Y

#define MUL0(X,Y,p0) do {			\
	*(p0) = (X) * (Y);			\
    } while(0)

#ifdef UINTD_T

#define HSHIFT bit_sizeof(UINT_T)

// MUL(X,Y,P1,P0)  (P1,P0) = X*Y

#define MUL(X,Y,P1,P0) do {			\
	UINTD_T TMPl = ((UINTD_T)(X))*(Y);	\
	*(P1) = TMPl >> HSHIFT;			\
	*(P0) = TMPl;				\
    } while(0)

// MUL1(X,Y,P1)  (P1,_) = X*Y

#define MUL1(X, Y, P1) do {			\
	UINTD_T Tl = ((UINTD_T)(X))*(Y);	\
	*(P1) = Tl >> HSHIFT;			\
    } while(0)

// MULA(X,Y,A,P1,P0)  (P1,P0) = X*Y + A

#define MULA(X, Y, A, P1, P0) do {			\
	UINTD_T Tl = ((UINTD_T)(X))*(Y) + (A);		\
	*(P1) = Tl >> HSHIFT;				\
	*(P0) = Tl;					\
    } while(0)

// MULA(X,Y,A,P1)  (P1,_) = X*Y + A

#define MUL1A(X, Y, A, P1) do {				\
	UINTD_T Tl = ((UINTD_T)(X))*(Y) + (A);		\
	*(P1) = Tl >> HSHIFT;				\
    } while(0)

// SQRA(X,A,P1,P0)  (P1,P0) = X*X + A

#define SQRA(X, A, P1, P0) do {				\
	UINTD_T Tl = ((UINTD_T)(X))*(X) + (A);		\
	*(P1) = Tl >> HSHIFT;				\
	*(P0) = Tl;					\
    } while(0)

// SQR(X,P1,P0)  (P1,P0) = X*X

#define SQR(X, P1, P0) do {			\
	UINTD_T Tl = ((UINTD_T)(X))*(X);	\
	*(P1) = Tl >> HSHIFT;			\
	*(P0) = Tl;				\
    } while(0)

// MULAB(X,Y,A,B,P1,P0)  (P1,P0) = X*X + A + B

#define MULAB(X, Y, A, B, P1, P0) do {				\
	UINTD_T Tl = ((UINTD_T)(X))*(Y) + (A) + (B);		\
	*(P1) = Tl >> HSHIFT;					\
	*(P0) = Tl;						\
    } while(0)

#elif defined(UINTH_T)

#define HSHIFT bit_sizeof(UINTH_T)
#define LMASK  ((UINT_T)((((UINT_T)1) << HSHIFT)-1))
#define HMASK  (LMASK << HSHIFT)

STATIC INLINE void mula(UINT_T src1, UINT_T src2, UINT_T a,
			UINT_T* prod1, UINT_T* prod0)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINTH_T b0 = src2;
    UINTH_T b1 = src2 >> HSHIFT;
    UINT_T a0b0 = ((UINT_T)a0)*b0;
    UINT_T a0b1 = ((UINT_T)a0)*b1;
    UINT_T a1b0 = ((UINT_T)a1)*b0;
    UINT_T a1b1 = ((UINT_T)a1)*b1;
    UINT_T p0,p1,p2,c0;

    ADD(a0b0,a,&c0,&p0);
    ADD(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    ADD(p1,a0b1,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b0,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    ADD(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

STATIC INLINE void sqra(UINT_T src1, UINT_T a,
			UINT_T* prod1, UINT_T* prod0)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINT_T a0a0 = ((UINT_T)a0)*a0;
    UINT_T a0a1 = ((UINT_T)a0)*a1;
    UINT_T a1a1 = ((UINT_T)a1)*a1;
    UINT_T p0,p1,p2,c0;

    ADD(a0a0,a,&c0,&p0);
    ADD(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    ADD(p1,a0a1,&c0,&p1);
    p2 += c0;
    ADD(p1,a0a1,&c0,&p1);
    p2 += c0;
    ADD(p1,a1a1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    ADD(a1a1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

// multiply add and return MSB
// (prod1,_) = src1*src2 + a
STATIC INLINE void mul1a(UINT_T src1, UINT_T src2, UINT_T a, UINT_T* prod1)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINTH_T b0 = src2;
    UINTH_T b1 = src2 >> HSHIFT;
    UINT_T a0b0 = ((UINT_T)a0)*b0;
    UINT_T a0b1 = ((UINT_T)a0)*b1;
    UINT_T a1b0 = ((UINT_T)a1)*b0;
    UINT_T a1b1 = ((UINT_T)a1)*b1;
    UINT_T p0,p1,p2,c0;

    ADD(a0b0,a,&c0,&p0);
    ADD(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    ADD(p1,a0b1,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b0,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 += c0;
    ADD(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
}

// multiply and return MSB
// (prod1,_) = src1*src2
STATIC INLINE void mul1(UINT_T src1, UINT_T src2, UINT_T* prod1)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINTH_T b0 = src2;
    UINTH_T b1 = src2 >> HSHIFT;
    UINT_T a0b0 = ((UINT_T)a0)*b0;
    UINT_T a0b1 = ((UINT_T)a0)*b1;
    UINT_T a1b0 = ((UINT_T)a1)*b0;
    UINT_T a1b1 = ((UINT_T)a1)*b1;
    UINT_T p0,p1,p2,c0;

    ADD(a0b0,a,&c0,&p0);
    ADD(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    ADD(p1,a0b1,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b0,&c0,&p1);
    p2 += c0;
    ADD(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    ADD(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
}

STATIC INLINE void sqr(UINT_T src1,UINT_T* prod1, UINT_T* prod0)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINT_T a0a0 = ((UINT_T)a0)*a0;
    UINT_T a0a1 = ((UINT_T)a0)*a1;
    UINT_T a1a1 = ((UINT_T)a1)*a1;
    UINT_T p0=a0a0,p1=p0>>HSHIFT,p2=0,c0;

    ADD(p1,a0a1,&c0,&p1);
    p2 += c0;
    ADD(p1,a0a1,&c0,&p1);
    p2 += c0;    
    ADD(p1,a1a1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    ADD(a1a1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

STATIC INLINE void mulab(UINT_T src1, UINT_T src2, UINT_T a, UINT_T b,
			 UINT_T* prod1, UINT_T* prod0)
{
    UINT_T c0;
    mula(src1, src2, a, prod1, prod0);
    ADD(*prod0, b, &c0, prod0);
    ADD(*prod1, c0, &c0, prod1);  // result is ignored
    assert(c0 == 0);
}

STATIC INLINE void mul(UINT_T src1, UINT_T src2,
		       UINT_T* prod1, UINT_T* prod0)
{
    mula(src1, src2, 0, prod1, prod0);
}

#endif
#endif
