#ifndef __BIG3_I__
#define __BIG3_I__

#ifndef INLINE
#define INLINE inline
#endif

// shift right one digit
static INLINE void zero3(UINT_T a[3])
{
    a[0] = 0;
    a[1] = 0;
    a[2] = 0;
}

static INLINE void add31(UINT_T a[3], UINT_T b0, UINT_T d[3])
{
    UINT_T c;
    add(a[0], b0, &c, &d[0]);
    add(a[1], c,  &c, &d[1]);
    add(a[2], c,  &c, &d[2]);
}

static INLINE void add32(UINT_T a[3], UINT_T b[2], UINT_T d[3])
{
    UINT_T c;
    add(a[0],b[0],&c,&d[0]);
    addc(a[1],b[1],c,&c,&d[1]);
    add(a[2],c,&c,&d[2]);
}

// shift right one digit
static INLINE void shr3(UINT_T a[3])
{
    a[0] = a[1];
    a[1] = a[2];
    a[2] = 0;
}


#endif
