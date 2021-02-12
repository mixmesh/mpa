#ifndef __BIG3R_I__
#define __BIG3R_I__

#ifndef INLINE
#define INLINE inline
#endif

#define decl3(x) UINT_T x[3]
#define elem3(x,i) x[(i)]

static INLINE void zero3(UINT_T* a)
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
}

static INLINE void add31(UINT_T* a, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add(a[1], c,  &c, &d[1]);
    add0(a[2], c,  &d[2]);
}

static INLINE void add31p(UINT_T* a, UINT_T b0)
{
    UINT_T c;
    add(a[0], b0, &c, &a[0]);
    add(a[1], c,  &c, &a[1]);
    add0(a[2], c,  &a[2]);
}

static INLINE void add32(UINT_T* a, UINT_T b1, UINT_T b0, UINT_T* d)
{
    UINT_T c;
    add(a[0],b0,&c,&d[0]);
    addc(a[1],b1,c,&c,&d[1]);
    add0(a[2],c,&d[2]);
}

// d = a  destructive update
static INLINE void add32p(UINT_T* a, UINT_T b1, UINT_T b0)
{
    UINT_T c;
    add(a[0],b0,&c,&a[0]);
    addc(a[1],b1,c,&c,&a[1]);
    add0(a[2],c,&a[2]);
}

// shift right one digit
static INLINE void shr3(UINT_T* a)
{
    a[0] = a[1];
    a[1] = a[2];
    a[2] = 0;
}

#endif
