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

#define bit_sizeof(T) (sizeof(T)*8)

#define ADD0(src1, src2, sum) do {		\
	*(sum) = (src1) + (src2);		\
    } while(0)

#ifdef USE_INLINE
static INLINE void add0(UINT_T src1, UINT_T src2, UINT_T* sum)
{
    *sum = src1 + src2;
}
#endif

#define ADD(src1, src2, co, sum) do {		\
    UINT_T SRC1l = (src1);			\
    SRC1l += (src2);				\
    *(sum) = SRC1l;				\
    *(co) = (SRC1l < (src2));			\
    } while(0)

#ifdef USE_INLINE
static INLINE void add(UINT_T src1, UINT_T src2, UINT_T* co, UINT_T* sum)
{
    src1 += src2;
    *sum = src1;
    *co = (src1 < src2);
}
#endif

#define ADDC(src1, src2, ci, co, sum) do {	\
    UINT_T SRC1l = (src1);			\
    UINT_T CIl = (ci);				\
    SRC1l += (ci);				\
    CIl = (SRC1l < (CIl));			\
    SRC1l += (src2);				\
    *(sum) = SRC1l;				\
    *(co) = CIl + (SRC1l < (src2));		\
    } while(0)


#ifdef USE_INLINE
static INLINE void addc(UINT_T src1, UINT_T src2, UINT_T ci, UINT_T* co, UINT_T* sum)
{
    src1 += ci;
    ci = (src1 < ci);
    src1 += src2;
    *sum = src1;
    *co = ci + (src1 < src2);
}
#endif

#define SUB0(src1, src2, diff) do { \
    *(diff) = (src1) - (src2);	    \
    } while(0)

#ifdef USE_INLINE
static INLINE void sub0(UINT_T src1, UINT_T src2, UINT_T* diff)
{
    *diff = src1 - src2;
}
#endif

#define SUB(src1, src2, bo, diff) do { \
    UINT_T SRC2l = (src2);	       \
    SRC2l = (src1) - SRC2l;	       \
    *(diff) = SRC2l;		       \
    *(bo) = (SRC2l > (src1));	       \
    } while(0)

#ifdef USE_INLINE
static INLINE void sub(UINT_T src1, UINT_T src2, UINT_T* bo, UINT_T* diff)
{
    src2 = src1 - src2;
    *diff = src2;
    *bo = (src2 > src1);
}
#endif

#define SUBB(src1, src2, bi, bo, diff) do {	\
    UINT_T SRC2l = (src2);			\
    UINT_T BIl = (bi);				\
    SRC2l += (bi);				\
    BIl = (SRC2l < BIl);			\
    SRC2l = (src1) - SRC2l;			\
    *(bo) = BIl + (SRC2l > (src1));		\
    *(diff) = SRC2l;				\
    } while(0)

#ifdef USE_INLINE
static INLINE void subb(UINT_T src1, UINT_T src2, UINT_T bi, UINT_T* bo, UINT_T* diff)
{
    src2 += bi;
    bi = (src2 < bi);
    src2 = src1 - src2;
    *bo = bi + (src2 > src1);
    *diff = src2;
}
#endif

// multiply and return LSB
#ifdef USE_INLINE
static INLINE void mul0(UINT_T src1, UINT_T src2, UINT_T* p0)
{
    *p0 = src1 * src2;
}
#endif

#define MUL0(src1,src2,p0) do { \
    *(p0) = (src1) * (src2);	\
    } while(0)


#ifdef UINTD_T

#define HSHIFT bit_sizeof(UINT_T)

#define MUL(src1,src2,p1,p0) do {		\
    UINTD_T TMPl = ((UINTD_T)(src1))*(src2);	\
    *(p1) = TMPl >> HSHIFT;			\
    *(p0) = TMPl;				\
    } while(0)
	  

#ifdef USE_INLINE
static INLINE void mul(UINT_T src1, UINT_T src2, UINT_T* p1, UINT_T* p0)
{
    UINTD_T t = ((UINTD_T)src1)*src2;
    *p1 = t >> HSHIFT;
    *p0 = t;
}
#endif

#define MUL1(src1, src2, p1) do {		\
    UINTD_T Tl = ((UINTD_T)(src1))*(src2);	\
    *(p1) = Tl >> HSHIFT;			\
    } while(0)

// multiply and return MSB
#ifdef USE_INLINE
static INLINE void mul1(UINT_T src1, UINT_T src2, UINT_T* p1)
{
    UINTD_T t = ((UINTD_T)src1)*src2;
    *p1 = t >> HSHIFT;
}
#endif

#define MULA(src1, src2, a, prod1, prod0) do {	\
    UINTD_T Tl = ((UINTD_T)(src1))*(src2) + (a);	\
    *(prod1) = Tl >> HSHIFT;				\
    *(prod0) = Tl;					\
    } while(0)

#ifdef USE_INLINE
static INLINE void mula(UINT_T src1, UINT_T src2, UINT_T a,
			UINT_T* prod1, UINT_T* prod0)
{
    UINTD_T t = ((UINTD_T)src1)*src2 + a;
    *prod1 = t >> HSHIFT;
    *prod0 = t;
}
#endif

#define MUL1A(src1, src2, a, p1) do {		\
    UINTD_T Tl = ((UINTD_T)(src1))*(src2) + (a);	\
    *(p1) = Tl >> HSHIFT;				\
    } while(0)

// multiply and return MSB
#ifdef USE_INLINE
static INLINE void mul1a(UINT_T src1, UINT_T src2, UINT_T a, UINT_T* p1)
{
    UINTD_T t = ((UINTD_T)src1)*src2 + a;
    *p1 = t >> HSHIFT;
}
#endif

#define SQRA(src1, a, prod1, prod0) do {	\
    UINTD_T Tl = ((UINTD_T)(src1))*(src1) + (a);	\
    *(prod1) = Tl >> HSHIFT;				\
    *(prod0) = Tl;					\
    } while(0)

#ifdef USE_INLINE
static INLINE void sqra(UINT_T src1, UINT_T a, UINT_T* prod1, UINT_T* prod0)
{
    UINTD_T t = ((UINTD_T)src1)*src1 + a;
    *prod1 = t >> HSHIFT;
    *prod0 = t;
}
#endif

#define SQR(src1, prod1, prod0) do {		\
    UINTD_T Tl = ((UINTD_T)(src1))*(src1);	\
    *(prod1) = Tl >> HSHIFT;			\
    *(prod0) = Tl;				\
    } while(0)

#ifdef USE_INLINE
static INLINE void sqr(UINT_T src1, UINT_T* prod1, UINT_T* prod0)
{
    UINTD_T t = ((UINTD_T)src1)*src1;
    *prod1 = t >> HSHIFT;
    *prod0 = t;
}
#endif

#define MULAB(src1, src2, a, b, prod1, prod0) do {	\
    UINTD_T Tl = ((UINTD_T)(src1))*(src2) + (a) + (b);	\
    *(prod1) = Tl >> HSHIFT;				\
    *(prod0) = Tl;					\
    } while(0)

#ifdef USE_INLINE
static INLINE void mulab(UINT_T src1, UINT_T src2, UINT_T a, UINT_T b,
			 UINT_T* prod1, UINT_T* prod0)
{
    UINTD_T t = ((UINTD_T)src1)*src2 + a + b;
    *prod1 = t >> HSHIFT;
    *prod0 = t;
}
#endif

#elif defined(UINTH_T)

#define HSHIFT bit_sizeof(UINTH_T)
#define LMASK  ((UINT_T)((((UINT_T)1) << HSHIFT)-1))
#define HMASK  (LMASK << HSHIFT)

static INLINE void mula(UINT_T src1, UINT_T src2, UINT_T a,
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

    add(a0b0,a,&c0,&p0);
    add(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    add(p1,a0b1,&c0,&p1);
    p2 += c0;
    add(p1,a1b0,&c0,&p1);
    p2 += c0;
    add(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    add(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

static INLINE void sqra(UINT_T src1, UINT_T a,
			UINT_T* prod1, UINT_T* prod0)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINT_T a0a0 = ((UINT_T)a0)*a0;
    UINT_T a0a1 = ((UINT_T)a0)*a1;
    UINT_T a1a1 = ((UINT_T)a1)*a1;
    UINT_T p0,p1,p2,c0;

    add(a0a0,a,&c0,&p0);
    add(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    add(p1,a0a1,&c0,&p1);
    p2 += c0;
    add(p1,a0a1,&c0,&p1);
    p2 += c0;
    add(p1,a1a1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    add(a1a1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

// multiply add and return MSB
// (prod1,_) = src1*src2 + a
static INLINE void mul1a(UINT_T src1, UINT_T src2, UINT_T a, UINT_T* prod1)
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

    add(a0b0,a,&c0,&p0);
    add(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    add(p1,a0b1,&c0,&p1);
    p2 += c0;
    add(p1,a1b0,&c0,&p1);
    p2 += c0;
    add(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 += c0;
    add(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
}

// multiply and return MSB
// (prod1,_) = src1*src2
static INLINE void mul1(UINT_T src1, UINT_T src2, UINT_T* prod1)
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

    add(a0b0,a,&c0,&p0);
    add(c0<<HSHIFT,p0>>HSHIFT,&p2,&p1);
    add(p1,a0b1,&c0,&p1);
    p2 += c0;
    add(p1,a1b0,&c0,&p1);
    p2 += c0;
    add(p1,a1b1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    add(a1b1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
}

static INLINE void sqr(UINT_T src1,UINT_T* prod1, UINT_T* prod0)
{
    UINTH_T a0 = src1;
    UINTH_T a1 = src1 >> HSHIFT;
    UINT_T a0a0 = ((UINT_T)a0)*a0;
    UINT_T a0a1 = ((UINT_T)a0)*a1;
    UINT_T a1a1 = ((UINT_T)a1)*a1;
    UINT_T p0=a0a0,p1=p0>>HSHIFT,p2=0,c0;

    add(p1,a0a1,&c0,&p1);
    p2 += c0;
    add(p1,a0a1,&c0,&p1);
    p2 += c0;    
    add(p1,a1a1<<HSHIFT,&c0,&p1);
    p2 +=c0;
    add(a1a1,p2<<HSHIFT,&c0,&p2);
    *prod1 = (p2 & HMASK) | (p1 >> HSHIFT);
    *prod0 = (p1 << HSHIFT) | (p0 & LMASK);
}

static INLINE void mulab(UINT_T src1, UINT_T src2, UINT_T a, UINT_T b,
			 UINT_T* prod1, UINT_T* prod0)
{
    UINT_T c0;
    mula(src1, src2, a, prod1, prod0);
    add(*prod0, b, &c0, prod0);
    add(*prod1, c0, &c0, prod1);  // result is ignored
    assert(c0 == 0);
}

static INLINE void mul(UINT_T src1, UINT_T src2,
		       UINT_T* prod1, UINT_T* prod0)
{
    mula(src1, src2, 0, prod1, prod0);
}

#endif

#endif
