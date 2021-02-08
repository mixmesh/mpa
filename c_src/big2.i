#ifndef __BIG2_I__
#define __BIG2_I__

#ifndef INLINE
#define INLINE inline
#endif

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

static INLINE void add22(UINT_T* a, UINT_T b[2], UINT_T* d)
{
    UINT_T c;
    add(a[0],b[0],&c,&d[0]);
    addc(a[1],b[1],c,&c,&d[1]);
}

// shift right one digit
static INLINE void shr2(UINT_T* a)
{
    a[0] = a[1];
    a[1] = 0;
}

#endif
