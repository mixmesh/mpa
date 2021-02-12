#ifndef __BIG2_I__
#define __BIG2_I__

#ifndef INLINE
#define INLINE inline
#endif

#define decl2(x) UINT_T x[2]
#define elem2(x,i) x[(i)]

static INLINE void zero2(UINT_T* a)
{
    a[0] = 0;
    a[1] = 0;
}

static INLINE void add21(UINT_T* a, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add0(a[1], c, &d[1]);
}

static INLINE void add21p(UINT_T* a, UINT_T b0)
{
    UINT_T c;
    add(a[0], b0, &c, &a[0]);
    add0(a[1], c, &a[1]);
}

static INLINE void add22(UINT_T* a, UINT_T b1, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0],b0,&c,&d[0]);
    addc(a[1],b1,c,&c,&d[1]);
}

static INLINE void add22p(UINT_T* a, UINT_T b1, UINT_T b0)
{
    UINT_T c;
    add(a[0],b0,&c,&a[0]);
    addc(a[1],b1,c,&c,&a[1]);
}

// shift right one digit
static INLINE void shr2(UINT_T* a)
{
    a[0] = a[1];
    a[1] = 0;
}

#endif
