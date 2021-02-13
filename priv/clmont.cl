// montgomery operations

#define INT_T   int
#define UINT_T  uint
#define UINTD_T ulong
#define UINTH_T ushort

#define INLINE inline
#define GLOBAL __global
#define CONST  __constant
#define LOCAL  __local
#define PRIVATE __private

#define bit_sizeof(T) (sizeof(T)*8)
#define D_SIZE (sizeof(UINT_T)*8)

// GENERATED: K,S,N,Np,1
// generate:  mpz:format_mont(mpz:mont_w_fips(32)).

CONST int mont_K = 1056;
CONST int mont_S = 33;
CONST UINT_T mont_N[33] = {0xA7D5BFEB,0x09AFCC70,0xAA7AAC04,0x4A922B39,0x12026166,0xE2D2F727,0xA6B7A79F,0x847FEF21,0xCC08533B,0xF7ECBEFB,0xE1138146,0xBE94E560,0xD6CA5ADB,0x2FD38236,0x20084F6D,0x60663161,0x8155C999,0xCF71EAFA,0x5790E9CB,0x2EFD631C,0xCEC4921B,0x9946F63D,0x7A4F1B62,0x19085E86,0x7A881BA5,0x5C6B3C8C,0x782E32AF,0xD00AF4C2,0xD9BF852C,0x7E3D4768,0xFFD2F4DF,0xA10AEC7B,0x00000006};
CONST UINT_T mont_Np[33] = {0x69578F3D,0xBC1EBB24,0x847D9C4E,0x5F6730B1,0x1F7AD209,0x4910512F,0x5E6233E9,0xC2BCDC19,0xC4A019A2,0x786271B4,0xD9F90959,0xD12E22CC,0x469B740E,0x7C5112AE,0x1E14E970,0x2FCA4DC5,0x39D33F80,0xD707B247,0xC836F87C,0x3D624C4F,0x600C1F1F,0x2B3F141A,0xC14F9DF5,0xD32D8B3B,0x7261FE5A,0x43DC4F4C,0x3362EF4E,0x1339EF57,0x59C7B1DD,0x868A70DC,0x1F243EAE,0x4E1A74E9,0x9C8B6ADE};
CONST UINT_T mont_1[33] = {0x09F55EBD,0xAA07C6B7,0x08394792,0x718C7E6C,0xB741ABB4,0x5F1C0B06,0xB05010D1,0xC47D70D6,0x62F378C8,0xE4A2BF5E,0xDE85C33C,0x3F83BEF4,0xB77C760A,0x6A81ADCF,0x696E1631,0x8C60D7E9,0x87167AE5,0x964EFF1F,0xDC9E95DA,0x5BE3C0DC,0x57A636F8,0xED1D6020,0x6B8C446E,0x93857C3A,0xEB7FD448,0x9C5A30C3,0x9AE73AF1,0x216B7F49,0x11F5A373,0xD4361B05,0x628F34FB,0x5B02131E,0x00000004};


#include "digit.i"
#include "fips3.i"

// version of above but wihtout comparsion x must be of size xl+1!
#define NSHIFT bit_sizeof(UINT_T)  // need -1 ?
static INLINE int big_norm1(GLOBAL UINT_T* z, CONST UINT_T* n, int s)
{
    UINT_T mask = ((-(INT_T)z[s]) >> NSHIFT);
    int i;
    UINT_T B;

    SUB(z[0], (n[0] & mask), &B, &z[0]);
    for (i = 1; i < s; i++) {
	SUBB(z[i], (n[i] & mask), B, &B, &z[i]);
    }
    return s;
}

static void constant_to_private(CONST UINT_T* src, PRIVATE UINT_T* dst, int n)
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

// a = (a^x) mod N
__kernel void pow(GLOBAL UINT_T* a,   // in/out a[i]
		  GLOBAL UINT_T* x,   // in
		  const unsigned int xbits,  // in (number of bits in e)
		  const unsigned int n) // in number of a's and r's

{
    UINT_T P[2][mont_S+2];
    UINT_T A[2][mont_S+2];
    int i;

    i = get_global_id(0);
    if (i < n) {
	int u, v;
	int c, d;
	int pos,j;

	u = 0; v = u^1;
	constant_to_private(mont_1, P[u], mont_S);   // = mont(1) !

	c = 0; d = c^1;
	global_to_private(&a[mont_S*i], A[c], mont_S);

	for (pos=1, j=0; pos < xbits; pos += D_SIZE, j++) {
	    UINT_T xj = x[j];
	    int size = (pos+D_SIZE < xbits) ? D_SIZE : xbits-pos;
	    int bit;
	    for (bit = 0; bit < size; bit++) {
		if (xj & 1) {
		    zero_to_private(P[v], mont_S+2);
		    big_mont_mul_fips(A[c],P[u],mont_N,mont_Np,P[v],mont_S);
		    u = v; v = u^1;
		}
		// A' = A^2 (mod R)
		zero_to_private(A[d], mont_S+2);		
		big_mont_sqr_fips(A[c],mont_N,mont_Np,A[d],mont_S);
		c = d; d = c^1;
		xj >>= 1;
	    }
	}
	zero_to_private(P[v], mont_S+2);
	big_mont_mul_fips(A[c],P[u],mont_N,mont_Np,P[v],mont_S);
	private_to_global(P[v], &a[mont_S*i], mont_S);
    }
}
