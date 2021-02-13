#ifndef __BIG3_I__
#define __BIG3_I__

#define DECL3(x) UINT_T x##_0,x##_1,x##_2
#define ELEM3(x,i) x##_##i
#define ZERO3(x) x##_0 = x##_1 = x##_2 = 0

// d/3 = r/3 + b
#define ADD31(r, b0, d) do { \
    UINT_T c;		      \
    add(r##_0,  b0, &c, &d[0]); \
    add(r##_1,  c,  &c, &d[1]); \
    add0(r##_2, c,  &d[2]);	\
    } while(0)

// d/3 = r/3 + b/2
#define ADD32(r, b1, b0, d) do {		\
    UINT_T c;					\
    add(r##_0,b0,&c,&d[0]);			\
    addc(r##_1,b1,c,&c,&d[1]);		\
    add0(r##_2,c,&d[2]);			\
    } while(0)

// r/3 += b/2
#define ADD32P(r, b1, b0) do {			\
    UINT_T c;					\
    add(r##_0,b0,&c,&r##_0);			\
    addc(r##_1,b1,c,&c,&r##_1);			\
    add0(r##_2,c,&r##_2);			\
    } while(0)

// shift right one digit
#define SHR3(r) r##_0 = r##_1, r##_1 = r##_2, r##_2 = 0

#endif
