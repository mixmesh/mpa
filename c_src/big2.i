#ifndef __BIG2_I__
#define __BIG2_I__

#define DECL2(x) UINT_T x[2]
#define ELEM2(x,i) x[(i)]

#ifdef USE_INLINE
static INLINE void zero2(UINT_T* a)
{
    a[0] = 0;
    a[1] = 0;
}
#endif

#define ZERO2(a) do { \
	(a)[0] = 0; (a)[1] = 0;			\
    } while(0)

#ifdef USE_INLINE
static INLINE void add21(UINT_T* a, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add0(a[1], c, &d[1]);
}
#endif

#define ADD21(a, b0, d) do { \
    UINT_T c;					\
    ADD(a[0], b0, &c, &d[0]);			\
    ADD0(a[1], c, &d[1]);			\
    } while(0)

#ifdef USE_INLINE
static INLINE void add21p(UINT_T* a, UINT_T b0)
{
    UINT_T c;
    add(a[0], b0, &c, &a[0]);
    add0(a[1], c, &a[1]);
}
#endif

#define ADD21p(a, b0) do { \
    UINT_T c;					\
    ADD((a)[0], b0, &c, &(a)[0]);		\
    ADD0((a)[1], c, &(a)[1]);			\
    } while(0)

#ifdef USE_INLINE
static INLINE void add22(UINT_T* a, UINT_T b1, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0],b0,&c,&d[0]);
    addc(a[1],b1,c,&c,&d[1]);
}
#endif

#define ADD22(a, b1, b0, d) do {		\
	UINT_T c;				\
	ADD((a)[0],(b0),&c,&(d)[0]);		\
	ADDC((a)[1],(b1),c,&c,&(d)[1]);		\
    } while(0)

#ifdef USE_INLINE
static INLINE void add22p(UINT_T* a, UINT_T b1, UINT_T b0)
{
    UINT_T c;
    add(a[0],b0,&c,&a[0]);
    addc(a[1],b1,c,&c,&a[1]);
}
#endif

#define ADD22p(a, b1, b0) do {			\
	UINT_T c;				\
	ADD((a)[0],(b0),&c,&(a)[0]);		\
	ADDC((a)[1],(b1),c,&c,&(a)[1]);		\
    } while(0)

// shift right one digit
#ifdef USE_INLINE
static INLINE void shr2(UINT_T* a)
{
    a[0] = a[1];
    a[1] = 0;
}
#endif

#define SHR2(a) do {				\
    (a)[0] = (a)[1];				\
    (a)[1] = 0;					\
    } while(0)

#endif
