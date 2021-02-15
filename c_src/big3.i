#ifndef __BIG3_I__
#define __BIG3_I__

#define DECL3(x) UINT_T x[3]
#define ELEM3(x,i) x[(i)]

// (A2,A1,A0) = (0,0,0)
#define ZERO3(a) do {				\
	(a)[0] = 0; (a)[1] = 0; (a)[2] = 0;	\
    } while(0)

// (D2,D1,D0) = (A2,A1,A0) + (B0)
#define ADD31(a, b0, d) do {				\
	UINT_T Cl;					\
	ADD(a[0], (b0), &Cl, &d[0]);			\
	ADD(a[1], Cl,  &Cl, &d[1]);			\
	ADD0(a[2], Cl, &d[2]);				\
    } while(0)

// (A2,A1,A0) += (B0)
#define ADD31p(a, b0) do {			\
	UINT_T Cl;				\
	ADD(a[0], (b0), &Cl, &a[0]);		\
	ADD(a[1], Cl,  &Cl, &a[1]);		\
	ADD0(a[2], Cl, &a[2]);			\
    } while(0)

// (D2,D1,D0) = (A2,A1,A0) + (B1,B0)
#define ADD32(a, b1, b0, d) do {		\
	UINT_T Cl;				\
	ADD(a[0],(b0),&Cl,&d[0]);		\
	ADDC(a[1],(b1),Cl,&Cl,&d[1]);		\
	ADD0(a[2],Cl,&d[2]);			\
    } while(0)

// (A2,A1,A0) += (B1,B0)
#define ADD32p(a, b1, b0) do {			\
	UINT_T Cl;				\
	ADD((a)[0],(b0),&Cl,&(a)[0]);		\
	ADDC((a)[1],(b1),Cl,&Cl,&(a)[1]);	\
	ADD0((a)[2],Cl,&(a)[2]);		\
    } while(0)

// (A2,A1,A0) = (0,A2,A1)
#define SHR3(a) do {				\
    (a)[0] = (a)[1];				\
    (a)[1] = (a)[2];				\
    (a)[2] = 0;					\
    } while(0)

#endif
