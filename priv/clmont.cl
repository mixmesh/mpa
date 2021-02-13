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

CONST int mont_K = 64;
CONST int mont_S = 2;
CONST UINT_T mont_N[2] = {0x0000000D,0x80000000};
CONST UINT_T mont_Np[2] = {0x3B13B13B,0x313B13B1};
CONST UINT_T mont_1[2] = {0xFFFFFFF3,0x7FFFFFFF};

#include "digit.i"
#include "fips3.i"

// version of above but wihtout comparsion x must be of size xl+1!
#define NSHIFT (bit_sizeof(UINT_T)-1)  // need -1 ?
static INLINE int big_norm1(PRIVATE UINT_T* z, CONST UINT_T* n, int s)
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
		    big_norm1(P[v], mont_N, mont_S);
		    u = v; v = u^1;
		}
		// A' = A^2 (mod R)
		zero_to_private(A[d], mont_S+2);		
		big_mont_sqr_fips(A[c],mont_N,mont_Np,A[d],mont_S);
		big_norm1(A[d], mont_N, mont_S);
		c = d; d = c^1;
		xj >>= 1;
	    }
	}
	zero_to_private(P[v], mont_S+2);
	big_mont_mul_fips(A[c],P[u],mont_N,mont_Np,P[v],mont_S);
	big_norm1(P[v], mont_N, mont_S);
	private_to_global(P[v], &a[mont_S*i], mont_S);
    }
}
