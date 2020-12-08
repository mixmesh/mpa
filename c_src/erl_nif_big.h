#ifndef __ERL_NIF_BIG_H__
#define __ERL_NIF_BIG_H__

#include <stdint.h>
#include "erl_nif.h"

typedef ERL_NIF_TERM ErlNifBigDigit;

#if (__SIZEOF_POINTER__ == 4) && defined(__SIZEOF_LONG_LONG__) && (__SIZEOF_LONG_LONG__ == 8)
/* Assume 32-bit machine with long long support */
typedef uint64_t   ErlNifBigDoubleDigit;
typedef uint16_t   ErlNifBigHalfDigit;
#define BIG_HAVE_DOUBLE_DIGIT 1

#elif (__SIZEOF_POINTER__ == 4)
/* Assume 32-bit machine with no long support */
#undef  BIG_HAVE_DOUBLE_DIGIT
// typedef uint32_t   ErlNifBigDigit;
typedef uint16_t  ErlNifBigHalfDigit;

#elif (__SIZEOF_POINTER__ == 8)
// typedef uint64_t ErlNifBigDigit;
typedef uint32_t ErlNifBigHalfDigit;
/* Assume 64-bit machine, does it exist 128 bit long long long ? */
#ifdef __SIZEOF_INT128__
typedef __uint128_t  ErlNifBigDoubleDigit;
#define BIG_HAVE_DOUBLE_DIGIT 1
#else
#undef  BIG_HAVE_DOUBLE_DIGIT
#endif
#else
#error "cannot determine machine size"
#endif

#define NUM_TMP_DIGITS 4
#define DIGIT_BITS (sizeof(ErlNifBigDigit)*8)

typedef struct 
{
    unsigned size;          // number of digits 
    unsigned sign;          // 1= negative, 0=none-negative
    ErlNifBigDigit* digits;  // least significant digit first D0 D1 .. Dsize-1
    ErlNifBigDigit  ds[NUM_TMP_DIGITS];
} ErlNifBignum;

#define BIGNUM_SIZE(arr)  (sizeof((arr))/sizeof(ErlNifBigDigit))

extern int enif_is_big(ErlNifEnv* env, ERL_NIF_TERM big_term);
extern int enif_inspect_big(ErlNifEnv* env,ERL_NIF_TERM big_term,ErlNifBignum* big);
extern ERL_NIF_TERM enif_make_number(ErlNifEnv* env, ErlNifBignum* big);

extern int enif_get_number(ErlNifEnv* env, ERL_NIF_TERM t, ErlNifBignum* big);
extern int enif_copy_number(ErlNifEnv* env, ErlNifBignum* big, size_t min_size);
extern void enif_release_number(ErlNifEnv* env, ErlNifBignum* big);
extern int enif_get_copy_number(ErlNifEnv* env, ERL_NIF_TERM t, 
				ErlNifBignum* big,size_t min_size);
#endif
