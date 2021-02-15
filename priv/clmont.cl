// -*- c -*-
// montgomery operations

// Clang definitions
// #ifdef __clang__
// #define __global __attribute__((address_space(1)))
// int get_global_id(int index);
// #endif

// unrolling
//#define UNROLL __attribute__((opencl_unroll_hint))
//#define UNROLLN(n) __attribute__((opencl_unroll_hint(n)))

#define UNROLL
#define UNROLLN(n)

#define INT_T   int
#define UINT_T  unsigned int
#define UINTD_T unsigned long
#define UINTH_T unsigned short

#define STATIC
#define INLINE inline
#define GLOBAL __global
#define CONST  __constant
#define LOCAL  __local
#define PRIVATE __private

#define bit_sizeof(T) (sizeof(T)*8)
#define D_SIZE (sizeof(UINT_T)*8)

// GENERATED: K,S,N,Np,1
// generate:  mpz:format_mont(mpz:mont_w(Type,32)).


#include "digit.i"

#include "cios.i"
CONST int mont_K = 1056;
#define mont_S 33
CONST UINT_T mont_N[33] = {0xA7D5BFEB,0x09AFCC70,0xAA7AAC04,0x4A922B39,0x12026166,0xE2D2F727,0xA6B7A79F,0x847FEF21,0xCC08533B,0xF7ECBEFB,0xE1138146,0xBE94E560,0xD6CA5ADB,0x2FD38236,0x20084F6D,0x60663161,0x8155C999,0xCF71EAFA,0x5790E9CB,0x2EFD631C,0xCEC4921B,0x9946F63D,0x7A4F1B62,0x19085E86,0x7A881BA5,0x5C6B3C8C,0x782E32AF,0xD00AF4C2,0xD9BF852C,0x7E3D4768,0xFFD2F4DF,0xA10AEC7B,0x00000006};
CONST UINT_T mont_Np[33] = {0x69578F3D,0xBC1EBB24,0x847D9C4E,0x5F6730B1,0x1F7AD209,0x4910512F,0x5E6233E9,0xC2BCDC19,0xC4A019A2,0x786271B4,0xD9F90959,0xD12E22CC,0x469B740E,0x7C5112AE,0x1E14E970,0x2FCA4DC5,0x39D33F80,0xD707B247,0xC836F87C,0x3D624C4F,0x600C1F1F,0x2B3F141A,0xC14F9DF5,0xD32D8B3B,0x7261FE5A,0x43DC4F4C,0x3362EF4E,0x1339EF57,0x59C7B1DD,0x868A70DC,0x1F243EAE,0x4E1A74E9,0x9C8B6ADE};
CONST UINT_T mont_1[33] = {0x09F55EBD,0xAA07C6B7,0x08394792,0x718C7E6C,0xB741ABB4,0x5F1C0B06,0xB05010D1,0xC47D70D6,0x62F378C8,0xE4A2BF5E,0xDE85C33C,0x3F83BEF4,0xB77C760A,0x6A81ADCF,0x696E1631,0x8C60D7E9,0x87167AE5,0x964EFF1F,0xDC9E95DA,0x5BE3C0DC,0x57A636F8,0xED1D6020,0x6B8C446E,0x93857C3A,0xEB7FD448,0x9C5A30C3,0x9AE73AF1,0x216B7F49,0x11F5A373,0xD4361B05,0x628F34FB,0x5B02131E,0x00000004};

#define MMUL(a,b,n,np,r,s) big_mont_mul_cios((a),(b),(n),(np),(r),(s))
#define MSQR(a,n,np,r,s)   big_mont_sqr_cios((a),(n),(np),(r),(s))

// subtract n if needed
STATIC INLINE int norm(PRIVATE UINT_T* z, CONST UINT_T* n, int s)
{
    if (z[s]) {
	int i;
	UINT_T B;
	SUB(z[0], n[0], &B, &z[0]);
	for (i = 1; i < s; i++) {
	    SUBB(z[i], n[i], B, &B, &z[i]);
	}
    }
    return s;
}

// multiply r = x*y x[n], y[n] r[n]! result is same size as operands
STATIC INLINE void big_mul0(PRIVATE UINT_T* x, PRIVATE UINT_T* y,
			    PRIVATE UINT_T* r, int n)
{
    int i;
    for (i = 0; i < n; i++) {
	UINT_T c = 0;
	int j, ij=i;
	for (j = 0; i+j < n; j++) {
	    MULAB(x[i],y[j],r[ij],c,&c,&r[ij]);
	    ij++;
	}
	r[ij] = c;  // top index is n
    }
}


static void constant_to_private(CONST UINT_T* src, PRIVATE UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = src[i];
}

static void constant_to_global(CONST UINT_T* src, GLOBAL UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = src[i];
}


static void global_to_private(GLOBAL UINT_T* src, PRIVATE UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = src[i];
}

static void private_to_global(PRIVATE UINT_T* src, GLOBAL UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = src[i];
}

static void zero_to_private(PRIVATE UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = 0;
}

static void zero_to_global(GLOBAL UINT_T* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = 0;
}

__kernel void mul0(GLOBAL UINT_T* a,   // in a[i]
		  GLOBAL UINT_T* b,   // in b[i]
		  GLOBAL UINT_T* r,   // out r[i]
		  const unsigned int n)
{
    UINT_T A[mont_S];
    UINT_T B[mont_S];
    UINT_T R[mont_S];
    const int i = get_global_id(0);

    if (i >= n) return;

    global_to_private(a, A, mont_S);
    global_to_private(b, B, mont_S);
    zero_to_private(R, mont_S);
    big_mul0(A, B, R, mont_S);
    private_to_global(R, r, mont_S);
}

void print_global(GLOBAL UINT_T* a, int n)
{
    printf("{");
    if (n > 0) {
	int i;
	printf(" 0x%08x", a[0]);
	for (i = 1; i < n; i++) {
	    printf(",0x%08x", a[i]);
	}
    }
    printf("}\n");
}

void print_private(PRIVATE UINT_T* a, int n)
{
    printf("{");
    if (n > 0) {
	int i;
	printf(" 0x%08x", a[0]);
	for (i = 1; i < n; i++) {
	    printf(",0x%08x", a[i]);
	}
    }
    printf("}\n");
}

// a = (a^x) mod N
__kernel void montpow(GLOBAL UINT_T* a,   // in/out a[i]
		      GLOBAL UINT_T* x,   // in
		      const unsigned int xbits,  // in (number of bits in e)
		      GLOBAL UINT_T* r,
		      const unsigned int n) // in number of a's and r's

{
    UINT_T P[2][mont_S+2];
    UINT_T A[2][mont_S+2];
    int u, v;
    int pos,j;
    const int i = get_global_id(0);

    if (i >= n) return;
    
    u = 0;
    constant_to_private(mont_1, P[u], mont_S);   // = mont(1) !

    v = 0;
    global_to_private(a+mont_S*i, A[v], mont_S);

    for (pos=1, j=0; pos < xbits; pos += D_SIZE, j++) {
	UINT_T xj = x[j];
	int size = (pos+D_SIZE < xbits) ? D_SIZE : xbits-pos;
	int bit;
	for (bit = 0; bit < size; bit++) {
	    if (xj & 1) {
		zero_to_private(P[!u], mont_S+2);
		MMUL(A[v],P[u],mont_N,mont_Np,P[!u],mont_S);
		norm(P[!u], mont_N, mont_S);
		printf("P[!u] = "); print_private(P[!u], mont_S);
		u = !u;
	    }
	    // A' = A^2 (mod R)
	    zero_to_private(A[!v], mont_S+2);
	    MSQR(A[v],mont_N,mont_Np,A[!v],mont_S);
	    norm(A[!v], mont_N, mont_S);
	    printf("A[!v] = "); print_private(A[!v], mont_S);
	    v = !v;
	    
	    xj >>= 1;
	}
    }
    zero_to_private(P[!u], mont_S+2);
    MMUL(A[v],P[u],mont_N,mont_Np,P[!u],mont_S);
    norm(P[!u], mont_N, mont_S);
    private_to_global(P[!u], r+mont_S*i, mont_S);
    printf("r = "); print_global(r, mont_S);
}
