#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdbool.h>
#include <time.h>
#include <gmp.h>
#include <erl_nif.h>
#include "dloglib.h"

#include "erl_nif_big.h"
#include "montmul.h"

#ifdef DEBUG
FILE *_log_file;
char _temp_str[4096];
#define GET_STR(_value) mpz_get_str(_temp_str, 10, _value)
#define LOG(...) \
do { \
  fprintf(_log_file, ##__VA_ARGS__); \
  fflush(_log_file); \
} while (0)
#define OPEN_LOG \
do { \
  if ((_log_file = fopen("/tmp/gmp.log", "a")) == NULL) { \
    return enif_make_tuple2(env, \
                            enif_make_atom(env, "error"), \
                            enif_make_string(env, strerror(errno), \
                                             ERL_NIF_LATIN1)); \
  } \
} while (0)
#define CLOSE_LOG fclose(_log_file)
#else
#define GET_STR(value) (void)0
#define LOG(...) (void)0
#define OPEN_LOG (void)0
#define CLOSE_LOG (void)0
#endif

bool export_binary(ErlNifBinary *bin, mpz_t value);
void import_binary(mpz_t value, ErlNifBinary *bin);

gmp_randstate_t rand_state;

/*
 * dlog
 */

static ERL_NIF_TERM _dlog(ErlNifEnv* env, int argc,
                          const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** dlog\n");

  ErlNifBinary h_bin;
  if (!enif_inspect_binary(env, argv[0], &h_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary g_bin;
  if (!enif_inspect_binary(env, argv[1], &g_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary p_bin;
  if (!enif_inspect_binary(env, argv[2], &p_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t h;
  mpz_init(h);
  import_binary(h, &h_bin);
  LOG("h = %s\n", GET_STR(h));

  mpz_t g;
  mpz_init(g);
  import_binary(g, &g_bin);
  LOG("g = %s\n", GET_STR(g));

  mpz_t p;
  mpz_init(p);
  import_binary(p, &p_bin);
  LOG("p = %s\n", GET_STR(p));

  mpz_t rop;
  mpz_init(rop);
  dlog(rop, h, g, p);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(h);
  mpz_clear(g);
  mpz_clear(p);
  
  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&h_bin);
    enif_release_binary(&g_bin);
    enif_release_binary(&p_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&h_bin);
    enif_release_binary(&g_bin);
    enif_release_binary(&p_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * generate_safe_prime
 */

// http://cs.uccs.edu/~gsc/pub/master/cmccullo/src/sources/Common/Common.cpp

static ERL_NIF_TERM _generate_safe_prime(ErlNifEnv* env, int argc,
                                         const ERL_NIF_TERM argv[]) {

  OPEN_LOG;
  LOG("**** generate_safe_prime\n");

  long len;
  if (!enif_get_long(env, argv[0], &len)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  LOG("len = %ld\n", len);

  mpz_t p;
  mpz_init(p);
  mpz_t q;
  mpz_init(q);
  mpz_t size; // 2 ^ len
  mpz_init(size);
  mpz_set_ui(size, 1 << len);
  mpz_t rop;
  mpz_init(rop);

  do {
    LOG(".");
    mpz_urandomb(p, rand_state, len);
    LOG("RANDOM: %s\n", GET_STR(p));
    mpz_add(p, p, size);
    LOG("RANDOM1: %s\n", GET_STR(p));
    mpz_nextprime(q, p);
    LOG("RANDOM2: %s\n", GET_STR(p));
    mpz_mul_ui(p, q, 2);
    LOG("RANDOM3: %s\n", GET_STR(p));
    mpz_add_ui(p, p, 1);

    if (mpz_probab_prime_p(p, 10) > 0) {
       LOG("rop = %s is prime\n", GET_STR(rop));
       mpz_set(rop, p);
       break;
    } else {
      LOG("rop = %s is not prime\n", GET_STR(rop));
    }
    mpz_sub_ui(p, q, 1);
    mpz_tdiv_q_ui(p, p, 2);
    if (mpz_probab_prime_p(p, 10) > 0) {
       LOG("rop = %s is prime\n", GET_STR(rop));
       mpz_set(rop, q);
       break;
    } else {
      LOG("rop = %s is not prime\n", GET_STR(rop));
    }
   } while (true);

  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(p);
  mpz_clear(q);
  mpz_clear(size);  

  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * mpz_gcd
 */

static ERL_NIF_TERM _mpz_gcd(ErlNifEnv* env, int argc,
                             const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_gcd\n");

  ErlNifBinary op1_bin;
  if (!enif_inspect_binary(env, argv[0], &op1_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary op2_bin;
  if (!enif_inspect_binary(env, argv[1], &op2_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t op1;
  mpz_init(op1);
  import_binary(op1, &op1_bin);
  LOG("op1 = %s\n", GET_STR(op1));

  mpz_t op2;
  mpz_init(op2);
  import_binary(op2, &op2_bin);
  LOG("op2 = %s\n", GET_STR(op2));

  mpz_t rop;
  mpz_init(rop);
  mpz_gcd(rop, op1, op2);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(op1);
  mpz_clear(op2);
  
  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);    
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * mpz_invert
 */

static ERL_NIF_TERM _mpz_invert(ErlNifEnv* env, int argc,
                                const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_invert\n");

  ErlNifBinary op1_bin;
  if (!enif_inspect_binary(env, argv[0], &op1_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary op2_bin;
  if (!enif_inspect_binary(env, argv[1], &op2_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t op1;
  mpz_init(op1);
  import_binary(op1, &op1_bin);
  LOG("op1 = %s\n", GET_STR(op1));

  mpz_t op2;
  mpz_init(op2);
  import_binary(op2, &op2_bin);
  LOG("op2 = %s\n", GET_STR(op2));

  mpz_t rop;
  mpz_init(rop);
  mpz_invert(rop, op1, op2);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(op1);
  mpz_clear(op2);

  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * mpz_lcm
 */

static ERL_NIF_TERM _mpz_lcm(ErlNifEnv* env, int argc,
                             const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_lcm\n");

  ErlNifBinary op1_bin;
  if (!enif_inspect_binary(env, argv[0], &op1_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary op2_bin;
  if (!enif_inspect_binary(env, argv[1], &op2_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t op1;
  mpz_init(op1);
  import_binary(op1, &op1_bin);
  LOG("op1 = %s\n", GET_STR(op1));

  mpz_t op2;
  mpz_init(op2);
  import_binary(op2, &op2_bin);
  LOG("op2 = %s\n", GET_STR(op2));

  mpz_t rop;
  mpz_init(rop);
  mpz_lcm(rop, op1, op2);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(op1);
  mpz_clear(op2);  

  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&op1_bin);
    enif_release_binary(&op2_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * mpz_powm
 */

static ERL_NIF_TERM _mpz_powm(ErlNifEnv* env, int argc,
                              const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_powm\n");

  ErlNifBinary base_bin;
  if (!enif_inspect_binary(env, argv[0], &base_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary exp_bin;
  if (!enif_inspect_binary(env, argv[1], &exp_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  ErlNifBinary mod_bin;
  if (!enif_inspect_binary(env, argv[2], &mod_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t base;
  mpz_init(base);
  import_binary(base, &base_bin);
  LOG("base = %s\n", GET_STR(base));

  mpz_t exp;
  mpz_init(exp);
  import_binary(exp, &exp_bin);
  LOG("exp = %s\n", GET_STR(exp));

  mpz_t mod;
  mpz_init(mod);
  import_binary(mod, &mod_bin);
  LOG("mod = %s\n", GET_STR(mod));

  mpz_t rop;
  mpz_init(rop);
  mpz_powm(rop, base, exp, mod);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(base);
  mpz_clear(exp);
  mpz_clear(mod);

  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&base_bin);
    enif_release_binary(&exp_bin);
    enif_release_binary(&mod_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&base_bin);
    enif_release_binary(&exp_bin);
    enif_release_binary(&mod_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}


/*
 * mpz_pow_ui
 */

static ERL_NIF_TERM _mpz_pow_ui(ErlNifEnv* env, int argc,
                                const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_pow_ui\n");

  ErlNifBinary base_bin;
  if (!enif_inspect_binary(env, argv[0], &base_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t base;
  mpz_init(base);
  import_binary(base, &base_bin);
  LOG("base = %s\n", GET_STR(base));

  long exp;
  if (!enif_get_long(env, argv[1], &exp)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }
  LOG("exp = %ld\n", exp);

  mpz_t rop;
  mpz_init(rop);
  mpz_pow_ui(rop, base, exp);
  LOG("rop = %s\n", GET_STR(rop));

  mpz_clear(base);

  ErlNifBinary rop_bin;
  if (!export_binary(&rop_bin, rop)) {
    enif_release_binary(&base_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_badarg(env);
  } else {
    enif_release_binary(&base_bin);
    mpz_clear(rop);
    CLOSE_LOG;
    return enif_make_binary(env, &rop_bin);
  }
}

/*
 * mpz_probab_prime_p
 */

static ERL_NIF_TERM _mpz_probab_prime_p(ErlNifEnv* env, int argc,
                                        const ERL_NIF_TERM argv[]) {
  OPEN_LOG;
  LOG("**** mpz_probab_prime_p\n");

  ErlNifBinary n_bin;
  if (!enif_inspect_binary(env, argv[0], &n_bin)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  int reps;
  if (!enif_get_int(env, argv[1], &reps)) {
    CLOSE_LOG;
    return enif_make_badarg(env);
  }

  mpz_t n;
  mpz_init(n);
  import_binary(n, &n_bin);
  LOG("n = %s\n", GET_STR(n));

  LOG("mpz_probab_prime_p(_, %d)\n", reps);
  int result = mpz_probab_prime_p(n, reps);
  LOG("result = %d\n", result);

  enif_release_binary(&n_bin);
  mpz_clear(n);
  CLOSE_LOG;
  return enif_make_int(env, result);
}

void import_big(mpz_t n, ErlNifBignum *big)
{
    mpz_import(n, big->size, -1, sizeof(ErlNifBigDigit), 0, 0, big->digits);
}

ERL_NIF_TERM make_big(ErlNifEnv* env, mpz_t zr)
{
    ErlNifBignum r;
    size_t size = (mpz_sizeinbase(zr, 2) + (8*sizeof(ErlNifBigDigit)) -1) /
	(8*sizeof(ErlNifBigDigit));
    size_t ndigits;
    ErlNifBigDigit ds[size]; // GCC stack allocation
    
    mpz_export((void *)ds, &ndigits, -1, sizeof(ErlNifBigDigit), 0, 0, zr);
    r.size = ndigits;
    r.sign = 0;
    r.digits = ds;
    return enif_make_number(env, &r);
}

/*
 * big_powm using erl_nif_big for marshalling (remove some overhead)
 */

static ERL_NIF_TERM _big_powm(ErlNifEnv* env, int argc,
			      const ERL_NIF_TERM argv[])
{
    ErlNifBignum a;
    ErlNifBignum n;
    ErlNifBignum mod;
    mpz_t za;
    mpz_t zn;
    mpz_t zmod;
    mpz_t zr;
    ERL_NIF_TERM r;
    
    if (!enif_get_number(env, argv[0], &a))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[1], &n))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[2], &mod))
	return enif_make_badarg(env);

    mpz_init(za);
    import_big(za, &a);

    mpz_init(zn);
    import_big(zn, &n);

    mpz_init(zmod);
    import_big(zmod, &mod);

    mpz_init(zr);
    mpz_powm(zr, za, zn, zmod);

    r = make_big(env, zr);
    mpz_clear(za);
    mpz_clear(zn);
    mpz_clear(zmod);
    mpz_clear(zr);
    return r;
}

static ERL_NIF_TERM _big_size(ErlNifEnv* env, int argc,
			      const ERL_NIF_TERM argv[])
{
    ErlNifBignum t;
    if (!enif_get_number(env, argv[0], &t))
	return enif_make_badarg(env);
    return enif_make_int(env,t.size);
}

static ERL_NIF_TERM _big_bits(ErlNifEnv* env, int argc,
			      const ERL_NIF_TERM argv[])
{
    ErlNifBignum t;
    ErlNifBigDigit h;
    int n;
    if (!enif_get_number(env, argv[0], &t))
	return enif_make_badarg(env);
    n = sizeof(ErlNifBigDigit)*8*(t.size - 1);
    h = t.digits[t.size-1];
    if ((n == 0) && (h == 0))
	n = 1;
    else {
	while(h) {
	    n++;
	    h >>= 1;
	}
    }
    return enif_make_int(env,n);
}

static int get_redc_type(ErlNifEnv* env, ERL_NIF_TERM arg, redc_type_t* type)
{
    if (arg == enif_make_atom(env, "default"))
	*type = REDC_DEFAULT;
    else if (arg == enif_make_atom(env, "sos"))
	*type = REDC_SOS;
    else if (arg == enif_make_atom(env, "fips"))
	*type = REDC_FIPS;
    else
	return 0;
    return 1;
}


// args mont_redc(T, N, Np, K)
static ERL_NIF_TERM _big_mont_redc(ErlNifEnv* env, int argc,
				   const ERL_NIF_TERM argv[])
{
    ErlNifBignum t;
    ErlNifBignum n;
    ErlNifBignum np;
    redc_type_t redc_type;

    if (!get_redc_type(env, argv[0], &redc_type))
	return enif_make_badarg(env);    
    if (!enif_get_number(env, argv[1], &t))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[2], &n))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[3], &np))
	return enif_make_badarg(env);
    {
	ErlNifBignum r;
	ErlNifBigDigit R[2*n.size+1];
	int rl;
	rl = big_mont_redc(redc_type,
			   t.digits, t.size,
			   n.digits, n.size,
			   np.digits, np.size,
			   R, BIGNUM_SIZE(R));
	r.size = rl;
	r.sign = 0;
	r.digits = R;
	return enif_make_number(env, &r);
    }
}

// args mont_redc(A, B, N, Np, K)
static ERL_NIF_TERM _big_mont_mul(ErlNifEnv* env, int argc,
				   const ERL_NIF_TERM argv[])
{
    ErlNifBignum a;
    ErlNifBignum b;
    ErlNifBignum n;
    ErlNifBignum np;
    redc_type_t redc_type;

    if (!get_redc_type(env, argv[0], &redc_type))
	return enif_make_badarg(env);    
    if (!enif_get_number(env, argv[1], &a))
	return enif_make_badarg(env);    
    if (!enif_get_number(env, argv[2], &b))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[3], &n))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[4], &np))
	return enif_make_badarg(env);
    
    {
	ErlNifBignum r;
	ErlNifBigDigit R[2*n.size+1];
	int rl;
	rl = big_mont_mul(redc_type,
			  a.digits, a.size,
			  b.digits, b.size,
			  n.digits, n.size,
			  np.digits, np.size,
			  R, BIGNUM_SIZE(R));
	r.size = rl;
	r.sign = 0;
	r.digits = R;
	return enif_make_number(env, &r);
    }
}

// args mont_redc(A, N, Np)
static ERL_NIF_TERM _big_mont_sqr(ErlNifEnv* env, int argc,
				  const ERL_NIF_TERM argv[])
{
    ErlNifBignum a;
    ErlNifBignum n;
    ErlNifBignum np;
    redc_type_t redc_type;

    if (!get_redc_type(env, argv[0], &redc_type))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[1], &a))
	return enif_make_badarg(env);    
    if (!enif_get_number(env, argv[2], &n))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[3], &np))
	return enif_make_badarg(env);
    
    {
	ErlNifBignum r;
	ErlNifBigDigit R[2*n.size+1];
	int rl;
	
	rl = big_mont_sqr(redc_type,
			  a.digits, a.size,
			  n.digits, n.size,
			  np.digits, np.size,
			  R, BIGNUM_SIZE(R));
	r.size = rl;
	r.sign = 0;
	r.digits = R;
	return enif_make_number(env, &r);
    }    
}


// args mont_pow(A, E, P, N, Np)
static ERL_NIF_TERM _big_mont_pow(ErlNifEnv* env, int argc,
				  const ERL_NIF_TERM argv[])
{
    ErlNifBignum a;
    ErlNifBignum e;
    ErlNifBignum p;
    ErlNifBignum n;
    ErlNifBignum np;
    redc_type_t redc_type;

    if (!get_redc_type(env, argv[0], &redc_type))
	return enif_make_badarg(env);        
    if (!enif_get_number(env, argv[1], &a))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[2], &e))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[3], &p))
	return enif_make_badarg(env);
    if (!enif_get_number(env, argv[4], &n))
	return enif_make_badarg(env);    
    if (!enif_get_number(env, argv[5], &np))
	return enif_make_badarg(env);

    {
	ErlNifBignum r;
	ErlNifBigDigit R[2*n.size+1];
	int rl;
	rl = big_mont_pow(redc_type,
			  a.digits, a.size,
			  e.digits, e.size,
			  p.digits, p.size,
			  n.digits, n.size,
			  np.digits, np.size,
			  R, BIGNUM_SIZE(R));
	r.size = rl;
	r.sign = 0;
	r.digits = R;
	return enif_make_number(env, &r);
    }
}


static ErlNifFunc nif_funcs[] = {
  {"dlog", 3, _dlog, ERL_NIF_DIRTY_JOB_CPU_BOUND},
  {"generate_safe_prime", 1, _generate_safe_prime, ERL_NIF_DIRTY_JOB_CPU_BOUND},
  {"mpz_gcd", 2, _mpz_gcd},
  {"mpz_invert", 2, _mpz_invert},
  {"mpz_lcm", 2, _mpz_lcm},
  {"mpz_powm", 3, _mpz_powm},
  {"mpz_pow_ui", 2, _mpz_pow_ui},
  {"mpz_probab_prime_p", 2, _mpz_probab_prime_p},
  {"big_size", 1, _big_size},
  {"big_bits", 1, _big_bits},
  {"big_powm", 3, _big_powm},
  {"big_mont_redc", 4, _big_mont_redc},
  {"big_mont_mul", 5, _big_mont_mul},
  {"big_mont_sqr", 4, _big_mont_sqr},
  {"big_mont_pow", 6, _big_mont_pow},  
};

static int load(ErlNifEnv* env, void** priv_data, ERL_NIF_TERM load_info) {
  time_t timeseed;
  time(&timeseed);
  gmp_randinit_default(rand_state);
  gmp_randseed_ui(rand_state, timeseed);
  return 0;
}

ERL_NIF_INIT(gmp_nif, nif_funcs, &load, NULL, NULL, NULL)

// Export and export of binaries

bool export_binary(ErlNifBinary *bin, mpz_t n) {
  size_t size = (mpz_sizeinbase(n, 2) + 7) / 8;
  LOG("size = %zu\n", size);
  if (!enif_alloc_binary(size, bin)) {
    return false;
  }
  size_t countp;
  mpz_export((void *)bin->data, &countp, 1, 1, 0, 0, n);
  // empty loop is probably removed, but you never know...
#ifdef DEBUG  
  for (int i = 0; i < countp; i++) {
    LOG("[%d] -> %d\n", i, bin->data[i]);
  }
#endif
  return true;
}

void import_binary(mpz_t n, ErlNifBinary *bin) {
  mpz_import(n, bin->size, 1, 1, 0, 0, bin->data);
}
