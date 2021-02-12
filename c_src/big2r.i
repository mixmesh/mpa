#ifndef __BIG2R_I__
#define __BIG2R_I__

#ifndef INLINE
#define INLINE inline
#endif

#define decl2(x) UINT_T x##_0,x##_1
#define elem2(x,i) x##_##i
#define zero2(x) x##_0 = x##_1 = 0

// d/2 = r/2 + b
#define add21(r, b0, d) do { \
    UINT_T c;		      \
    add(r##_0,  b0, &c, &d[0]); \
    add0(r##_1,  c,  &d[1]); \
    } while(0)

// r/2 += b
#define add21p(r, b0) do { \
    UINT_T c;		      \
    add(r##_0,  b0, &c, &r##_0); \
    add0(r##_1,  c,  &r##_1); \
    } while(0)

// d/2 = r/2 + b/2
#define add22(r, b1, b0, d) do {		\
    UINT_T c;					\
    add(r##_0,b0,&c,&d[0]);			\
    add0(r##_1,c,&d[1]);			\
    } while(0)

// r/2 += b/2
#define add22p(r, b1, b0) do {			\
    UINT_T c;					\
    add(r##_0,b0,&c,&r##_0);			\
    add0(r##_1,c,&r##_1);			\
    } while(0)

// shift right one digit
#define shr2(r) r##_0 = r##_1, r##_1 = 0

#endif
