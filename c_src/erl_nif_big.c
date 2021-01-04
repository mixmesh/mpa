
#include "erl_nif_big.h"

#define TAG_PRIMARY_HEADER	0x0

#define _TAG_PRIMARY_SIZE	2
#define _TAG_PRIMARY_MASK	0x3
#define TAG_PRIMARY_BOXED	0x2

#define POS_BIG_SUBTAG          (0x2 << _TAG_PRIMARY_SIZE)
#define NEG_BIG_SUBTAG          (0x3 << _TAG_PRIMARY_SIZE)
#define _BIG_SIGN_BIT		(0x1 << _TAG_PRIMARY_SIZE)

#define _TAG_HEADER_POS_BIG	(TAG_PRIMARY_HEADER|POS_BIG_SUBTAG)
#define _TAG_HEADER_NEG_BIG	(TAG_PRIMARY_HEADER|NEG_BIG_SUBTAG)

#define _TAG_HEADER_MASK	0x3F
#define _HEADER_SUBTAG_MASK	0x3C	/* 4 bits for subtag */
#define _HEADER_ARITY_OFFS	6

#define EXPAND_POINTER(AnEterm) ((ERL_NIF_TERM) (AnEterm))

#define _is_bignum_header(x)	(((x) & (_TAG_HEADER_MASK-_BIG_SIGN_BIT)) == _TAG_HEADER_POS_BIG)
#define boxed_val(x) ((ERL_NIF_TERM*) EXPAND_POINTER(((x) - TAG_PRIMARY_BOXED)))
#define ptr_val(x)   ((ERL_NIF_TERM*)((x) & ~((unsigned long) 0x3)))

#define is_boxed(x)	(((x) & _TAG_PRIMARY_MASK) == TAG_PRIMARY_BOXED)
#define is_big(x)	(is_boxed((x)) && _is_bignum_header(*boxed_val((x))))

#define _header_arity(x)	((x) >> _HEADER_ARITY_OFFS)


// copy sl, or at most dl bytes, from src to dst
// return -1 if not all byte where copied
int enif_digits_copy(ErlNifBigDigit* dst, int dl, ErlNifBigDigit* src, int sl)
{
    int i;
    int n = (sl > dl) ? dl : sl;
    for (i = 0; i < n; i++)
	dst[i] = src[i];
    return (sl > dl) ? -1 : sl;
}

int enif_digits_small(ErlNifBigDigit* x, ErlNifBigDigit d)
{
    x[0] = d;
    return 1;
}

void enif_digits_zero(ErlNifBigDigit* dst, int n)
{
    int i;
    for (i = 0; i < n; i++)
	dst[i] = 0;
}

// #define USE_MEMCPY
// #define USE_MEMSET

// copy and zero pad MSB or truncate if needed
int enif_digits_copyz(ErlNifBigDigit* dst, int dl, ErlNifBigDigit* src, int sl)
{
    int i;
    int n = (sl > dl) ? dl : sl;

#ifdef USE_MEMCPY
    memcpy(dst, src, n*sizeof(ErlNifBigDigit));
#else
    for (i = 0; i < n; i++) dst[i] = src[i];
#endif
    if (sl > dl)
	return -1;  // truncated
#ifdef USE_MEMSET
    memset(dst+i, 0, (dl-sl)*sizeof(ErlNifBigDigit));
#else
    for (i = n; i < dl; i++) dst[i] = 0;
#endif
    return dl;
}

int enif_is_big(ErlNifEnv* env, ERL_NIF_TERM big_term)
{
    return is_big(big_term);
}

int enif_inspect_big(ErlNifEnv* env,ERL_NIF_TERM big_term,ErlNifBignum* big)
{
    unsigned long* ptr;

    if (!is_boxed(big_term))
	return 0;
    ptr = boxed_val(big_term);
    if (!_is_bignum_header(*ptr))
	return 0;
    big->sign   = *ptr & _BIG_SIGN_BIT;
    big->size   = _header_arity(*ptr);
    big->digits = ptr + 1;
    return 1;
}

// Create a new number (big or small)
// We fake this by making a tuple, and patch it
ERL_NIF_TERM enif_make_number(ErlNifEnv* env, ErlNifBignum* big)
{
    ERL_NIF_TERM t;
    ERL_NIF_TERM* ptr;
    unsigned size = big->size;
    
    // trim off zeros from high digits
    while((size > 1) && !big->digits[size-1])
	size--;
    if (size == 1) {
	// use make_uint64 / int64 to cover "all" cases
	if (!big->sign)
	    return enif_make_uint64(env, big->digits[0]);
	else if (!(big->digits[0] >> (DIGIT_BITS - 1))) {
	    ErlNifSInt64 d = (ErlNifSInt64) big->digits[0];
	    return enif_make_int64(env, -d);
	}
    }
    // this works since the elements are not checked!!!
    t = enif_make_tuple_from_array(env,(ERL_NIF_TERM*)big->digits,size);
    ptr = ptr_val(t);
    *ptr = (*ptr & ~_TAG_HEADER_MASK) | 
	(big->sign ? _TAG_HEADER_NEG_BIG : _TAG_HEADER_POS_BIG);
    return t;
}

// Load numer, big or small as ErlNigBignum
int enif_get_number(ErlNifEnv* env, ERL_NIF_TERM t, ErlNifBignum* big)
{
    if (enif_inspect_big(env, t, big))
	return 1;
    else {
	ErlNifSInt64 digit;
	if (enif_get_int64(env, t, &digit)) {
	    big->size = 1;
	    if (digit < 0) {
		big->sign = 1;
		big->ds[0] = -digit;
	    }
	    else {
		big->sign = 0;
		big->ds[0] = digit;
	    }
	    big->digits = big->ds;
	    return 1;
	}
    }
    return 0;
}

// get number with fixed number of digits, set MSB to zero if needed
// if big->size == s then the digits are not touched or patched
// if big->size < s then ds digits are used
int enif_get_number_ds(ErlNifEnv* env, ERL_NIF_TERM t,
		       ErlNifBignum* big, ErlNifBigDigit* ds, int s)
{
    if (!enif_get_number(env, t, big))
	return 0;
    if (big->size > s)
	return 0;
    if (big->size < s) {
	enif_digits_copyz(ds, s, big->digits, big->size);
	big->digits = ds;
	big->size = s;
    }
    return 1;
}

// Copy a bignum to a "safe" location, allow modifications of big num digits
int enif_copy_number(ErlNifEnv* env, ErlNifBignum* big, size_t min_size)
{
    size_t n;
    ErlNifBigDigit* digits;
    int i;

    n = (min_size > big->size) ? min_size : big->size;
    if (n <= NUM_TMP_DIGITS)
	digits = &big->ds[0];
    else
	digits = (ErlNifBigDigit*) enif_alloc(sizeof(ErlNifBigDigit)*n);
    if (!digits)
	return 0;
    i = 0;
    while(i < big->size) {
	digits[i] = big->digits[i];
	i++;
    }
    while(i < n) {
	digits[i] = 0;
	i++;
    }
    big->size   = n;
    big->digits = digits;
    return 1;
}

//  Release temporary memory associated with ErlNifBignum
// FIXME: debug: make sure (flag) digits are allocated (using copy_number)
void enif_release_number(ErlNifEnv* env, ErlNifBignum* big)
{
    if (big->digits && (big->digits != big->ds))
	enif_free(big->digits);
}

// Get and copy
int enif_get_copy_number(ErlNifEnv* env, ERL_NIF_TERM t, ErlNifBignum* big,
			 size_t min_size)
{
    if (!enif_get_number(env, t, big))
	return 0;
    return enif_copy_number(env, big, min_size);
}
