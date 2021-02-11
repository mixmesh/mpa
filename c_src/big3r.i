#ifndef __BIG3R_I__
#define __BIG3R_I__

#ifndef INLINE
#define INLINE inline
#endif

#define decl3r(r) UINT_T n##_0,n##_1,n##_2

#define zero3r(r) r##_0 = r##_1 = r##_2 = 0

// d/3 = r/3 + b
#define add31r(r, b0, d) do { \
    UINT_T c;		      \
    add(r##_0,  b0, &c, &d[0]); \
    add(r##_1,  c,  &c, &d[1]); \
    add0(r##_2, c,  &d[2]);	\
    } while(0)

// d/3 = r/3 + b/2
#define add32r(r, b1, b0, d) do {		\
    UINT_T c;					\
    add(r##_0,b0,&c,&d[0]);			\
    addc(r##_1,b1,c,&c,&d[1]);		\
    add0(r##_2,c,&d[2]);			\
    } while(0)

// r/3 += b/2
#define add32rp(r, b1, b0) do {			\
    UINT_T c;					\
    add(r##_0,b0,&c,&r##_0);			\
    addc(r##_1,b1,c,&c,&r##_1);		\
    add0(r##_2,c,&r##_2);			\
    } while(0)

// shift right one digit
#define shr3r(r) r##_0 = r##_1, r##_1 = r##_2, r##_2 = 0

#endif
