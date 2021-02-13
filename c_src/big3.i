#ifndef __BIG3_I__
#define __BIG3_I__

#define DECL3(x) UINT_T x[3]
#define ELEM3(x,i) x[(i)]

#ifdef USE_INLINE
static INLINE void zero3(UINT_T* a)
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
}
#endif

#define ZERO3(a) do { \
	(a)[0] = 0; (a)[1] = 0; (a)[2] = 0;	\
    } while(0)

#ifdef USE_INLINE
static INLINE void add31(UINT_T* a, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add(a[1], c,  &c, &d[1]);
    add0(a[2], c,  &d[2]);
}
#endif

#define ADD31(a, b0, d) do {				\
    UINT_T Cl;						\
    ADD(a[0], (b0), &Cl, &d[0]);			\
    ADD(a[1], Cl,  &Cl, &d[1]);				\
    ADD0(a[2], Cl, &d[2]);				\
	 } while(0)


#ifdef USE_INLINE
static INLINE void add31p(UINT_T* a, UINT_T b0)
{
    UINT_T c;
    add(a[0], b0, &c, &a[0]);
    add(a[1], c,  &c, &a[1]);
    add0(a[2], c,  &a[2]);
}
#endif

#define ADD31p(a, b0) do {			\
	UINT_T c;				\
	ADD(a[0], (b0), &c, &a[0]);		\
	ADD(a[1], c,  &c, &a[1]);		\
	ADD0(a[2], c, &a[2]);			\
    } while(0)

#ifdef USE_INLINE
static INLINE void add32(UINT_T* a, UINT_T b1, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0],b0,&c,&d[0]);
    addc(a[1],b1,c,&c,&d[1]);
    add0(a[2],c,&d[2]);
}
#endif

#define ADD32(a, b1, b0, d) do {		\
	UINT_T Cl;				\
	ADD(a[0],(b0),&Cl,&d[0]);		\
	ADDC(a[1],(b1),Cl,&Cl,&d[1]);		\
	ADD0(a[2],Cl,&d[2]);			\
    } while(0)

// d = a  destructive update
#ifdef USE_INLINE
static INLINE void add32p(UINT_T* a, UINT_T b1, UINT_T b0)
{
    UINT_T c;
    add(a[0],b0,&c,&a[0]);
    addc(a[1],b1,c,&c,&a[1]);
    add0(a[2],c,&a[2]);
}
#endif

#define ADD32p(a, b1, b0) do {			\
    UINT_T c;					\
    ADD((a)[0],(b0),&c,&(a)[0]);		\
    ADDC((a)[1],(b1),c,&c,&(a)[1]);		\
    ADD0((a)[2],c,&(a)[2]);			\
    } while(0)


// shift right one digit
#ifdef USE_INLINE
static INLINE void shr3(UINT_T* a)
{
    a[0] = a[1];
    a[1] = a[2];
    a[2] = 0;
}
#endif

#define SHR3(a) do {				\
    (a)[0] = (a)[1];				\
    (a)[1] = (a)[2];				\
    (a)[2] = 0;					\
    } while(0)

#endif
