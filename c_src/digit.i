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

// MULA(X,Y,A,P1,P0)  (P1,P0) = X*Y + A

STATIC INLINE void MULA(UINT_T X, UINT_T Y, UINT_T A,
			UINT_T* P1, UINT_T* P0)
{
    UINTH_T A0l = (X);
    UINTH_T A1l = (X) >> HSHIFT;
    UINTH_T B0l = (Y);
    UINTH_T B1l = (Y) >> HSHIFT;
    UINT_T A0B0l = ((UINT_T)A0l)*B0l;
    UINT_T A0B1l = ((UINT_T)A0l)*B1l;
    UINT_T A1B0l = ((UINT_T)A1l)*B0l;
    UINT_T A1B1l = ((UINT_T)A1l)*B1l;
    UINT_T P0l,P1l,P2l,C0l;

    ADD(A0B0l,A,&C0l,&P0l);
    ADD(C0l<<HSHIFT,P0l>>HSHIFT,&P2l,&P1l);
    ADD(P1l,A0B1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B0l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B1l<<HSHIFT,&C0l,&P1l);
    P2l +=C0l;
    ADD(A1B1l,P2l<<HSHIFT,&C0l,&P2l);
    *P1 = (P2l & HMASK) | (P1l >> HSHIFT);
    *P0 = (P1l << HSHIFT) | (P0l & LMASK);
}

STATIC INLINE void SQRA(UINT_T X, UINT_T A,
			UINT_T* P1, UINT_T* P0)
{
    UINTH_T A0l = (X);
    UINTH_T A1l = (X) >> HSHIFT;
    UINT_T A0A0l = ((UINT_T)A0l)*A0l;
    UINT_T A0A1l = ((UINT_T)A0l)*A1l;
    UINT_T A1A1l = ((UINT_T)A1l)*A1l;
    UINT_T P0l,P1l,P2l,C0l;

    ADD(A0A0l,A,&C0l,&P0l);
    ADD(C0l<<HSHIFT,P0l>>HSHIFT,&P2l,&P1l);
    ADD(P1l,A0A1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A0A1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1A1l<<HSHIFT,&C0l,&P1l);
    P2l +=C0l;
    ADD(A1A1l,P2l<<HSHIFT,&C0l,&P2l);
    *P1 = (P2l & HMASK) | (P1l >> HSHIFT);
    *P0 = (P1l << HSHIFT) | (P0l & LMASK);
}

// multiply add and return MSB
// (P1,_) = X*Y + a
STATIC INLINE void MUL1A(UINT_T X, UINT_T Y, UINT_T A, UINT_T* P1)
{
    UINTH_T A0l = (X);
    UINTH_T A1l = (X) >> HSHIFT;
    UINTH_T B0l = (Y);
    UINTH_T B1l = (Y) >> HSHIFT;
    UINT_T A0B0l = ((UINT_T)A0l)*B0l;
    UINT_T A0B1l = ((UINT_T)A0l)*B1l;
    UINT_T A1B0l = ((UINT_T)A1l)*B0l;
    UINT_T A1B1l = ((UINT_T)A1l)*B1l;
    UINT_T P0l,P1l,P2l,C0l;

    ADD(A0B0l,A,&C0l,&P0l);
    ADD(C0l<<HSHIFT,P0l>>HSHIFT,&P2l,&P1l);
    ADD(P1l,A0B1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B0l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B1l<<HSHIFT,&C0l,&P1l);
    P2l += C0l;
    ADD(A1B1l,P2l<<HSHIFT,&C0l,&P2l);
    *P1 = (P2l & HMASK) | (P1l >> HSHIFT);
}

// multiply and return MSB
// (P1,_) = X*Y
STATIC INLINE void MUL1(UINT_T X, UINT_T Y, UINT_T* P1)
{
    UINTH_T A0l = (X);
    UINTH_T A1l = (X) >> HSHIFT;
    UINTH_T B0l = (Y);
    UINTH_T B1l = (Y) >> HSHIFT;
    UINT_T A0B0l = ((UINT_T)A0l)*B0l;
    UINT_T A0B1l = ((UINT_T)A0l)*B1l;
    UINT_T A1B0l = ((UINT_T)A1l)*B0l;
    UINT_T A1B1l = ((UINT_T)A1l)*B1l;
    UINT_T P0l=A0B0l,P1l=P0l>>HSHIFT,P2l=0,C0l;
    
    ADD(P1l,A0B1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B0l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A1B1l<<HSHIFT,&C0l,&P1l);
    P2l += C0l;
    ADD(A1B1l,P2l<<HSHIFT,&C0l,&P2l);
    *P1 = (P2l & HMASK) | (P1l >> HSHIFT);
}

STATIC INLINE void SQR(UINT_T X,UINT_T* P1, UINT_T* P0)
{
    UINTH_T A0l = (X);
    UINTH_T A1l = (X) >> HSHIFT;
    UINT_T A0A0l = ((UINT_T)A0l)*A0l;
    UINT_T A0A1l = ((UINT_T)A0l)*A1l;
    UINT_T A1A1l = ((UINT_T)A1l)*A1l;
    UINT_T P0l=A0A0l,P1l=P0l>>HSHIFT,P2l=0,C0l;

    ADD(P1l,A0A1l,&C0l,&P1l);
    P2l += C0l;
    ADD(P1l,A0A1l,&C0l,&P1l);
    P2l += C0l;    
    ADD(P1l,A1A1l<<HSHIFT,&C0l,&P1l);
    P2l +=C0l;
    ADD(A1A1l,P2l<<HSHIFT,&C0l,&P2l);
    *P1 = (P2l & HMASK) | (P1l >> HSHIFT);
    *P0 = (P1l << HSHIFT) | (P0l & LMASK);
}

STATIC INLINE void MULAB(UINT_T X, UINT_T Y, UINT_T A, UINT_T B,
			 UINT_T* P1, UINT_T* P0)
{
    UINT_T C0l;
    MULA(X, Y, A, P1, P0);
    ADD(*P0, B, &C0l, P0);
    ADD(*P1, C0l, &C0l, P1);  // result is ignored
    assert(C0l == 0);
}

STATIC INLINE void MUL(UINT_T X, UINT_T Y,
		       UINT_T* P1, UINT_T* P0)
{
    MULA(X, Y, 0, P1, P0);
}

#endif
#endif
