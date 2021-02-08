#ifndef __BIG4_I__
#define __BIG4_I__

#ifndef INLINE
#define INLINE inline
#endif

static INLINE void zero4(UINT_T* a)
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
}

static INLINE void add41(UINT_T* a, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add(a[1], c,  &c, &d[1]);
    add(a[2], c,  &c, &d[2]);
    add0(a[3], c, &d[3]);
}

static INLINE void add42(UINT_T* a, UINT_T b[2], UINT_T* d)
{
    UINT_T c;
    add(a[0],b[0],&c,&d[0]);
    addc(a[1],b[1],c,&c,&d[1]);
    add(a[2],c,&c,&d[2]);
    add0(a[3],c,&d[3]);
}

// shift right one digit
static INLINE void shr4(UINT_T* a)
{
    a[0] = a[1];
    a[1] = a[2];
    a[2] = a[3];
    a[3] = 0;
}

#endif
