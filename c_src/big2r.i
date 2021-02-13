#ifndef __BIG2R_I__
#define __BIG2R_I__

#define DECL2(x) UINT_T x##_0,x##_1
#define ELEM2(x,i) x##_##i
#define ZERO2(x) x##_0 = x##_1 = 0

// d/2 = r/2 + b
#define ADD21(r, b0, d) do { \
    UINT_T Cl;		      \
    ADD(r##_0,  b0, &Cl, &d[0]); \
    ADD0(r##_1,  Cl,  &d[1]); \
    } while(0)

// r/2 += b
#define ADD21P(r, b0) do { \
    UINT_T Cl;		      \
    ADD(r##_0,  b0, &Cl, &r##_0); \
    ADD0(r##_1,  Cl,  &r##_1);					\
    } while(0)

// d/2 = r/2 + b/2
#define ADD22(r, b1, b0, d) do {		\
    UINT_T Cl;					\
    ADD(r##_0,b0,&Cl,&d[0]);			\
    ADD0(r##_1,Cl,&d[1]);			\
    } while(0)

// r/2 += b/2
#define ADD22P(r, b1, b0) do {			\
    UINT_T Cl;					\
    ADD(r##_0,b0,&Cl,&r##_0);			\
    ADD0(r##_1,Cl,&r##_1);			\
    } while(0)

// shift right one digit
#define SHR2(r) r##_0 = r##_1, r##_1 = 0

#endif
