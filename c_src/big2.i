#ifndef __BIG2_I__
#define __BIG2_I__

#define DECL2(x) UINT_T x[2]
#define ELEM2(x,i) x[(i)]

// (A1,A0) = (0,0,0)
#define ZERO2(a) do { \
	(a)[0] = 0; (a)[1] = 0;			\
    } while(0)

// (D1,D0) = (A1,A0) + (B0)
#define ADD21(a, b0, d) do { \
    UINT_T c;					\
    ADD(a[0], b0, &c, &d[0]);			\
    ADD0(a[1], c, &d[1]);			\
    } while(0)

// (A1,A0) += (B0)
#define ADD21p(a, b0) do { \
    UINT_T c;					\
    ADD((a)[0], b0, &c, &(a)[0]);		\
    ADD0((a)[1], c, &(a)[1]);			\
    } while(0)

// (D1,D0) = (A1,A0) + (B1,B0)
#define ADD22(a, b1, b0, d) do {		\
	UINT_T c;				\
	ADD((a)[0],(b0),&c,&(d)[0]);		\
	ADDC((a)[1],(b1),c,&c,&(d)[1]);		\
    } while(0)

// (A1,A0) += (B1,B0)
#define ADD22p(a, b1, b0) do {			\
	UINT_T c;				\
	ADD((a)[0],(b0),&c,&(a)[0]);		\
	ADDC((a)[1],(b1),c,&c,&(a)[1]);		\
    } while(0)

// (A1,A0) = (0,A1)
#define SHR2(a) do {				\
    (a)[0] = (a)[1];				\
    (a)[1] = 0;					\
    } while(0)

#endif
