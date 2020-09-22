#include <stdio.h>
#include <time.h>
#include <gmp.h>

// https://www.luke.maurits.id.au/files/misc/honours_thesis.pdf
// Compile: gcc -std=c99 -Wall -Werror dlog.c -lgmp -o dlog

/* Function prototypes */

void pollardrho(mpz_t * result, mpz_t h, mpz_t g, mpz_t p);
void iterate(mpz_t * x, mpz_t * a, mpz_t * b);

/* Global variables */

mpz_t h, g, p, pminus1;
gmp_randstate_t state;

/****************
 * Main program *
 ****************/

int main(int argc, char *argv[]) {

  /* Initialise and set global variables */

  // h = 12, g = 3, p = 23 results in 4
/*
  mpz_init_set_str(h, argv[1], 10);
  mpz_init_set_str(g, argv[2], 10);
  mpz_init_set_str(p, argv[3], 10);
*/

  // NOTE: Solving the discrete logarithm for the follwoing H, G and P takes
  // forever. Solving it with elgamal:brute_force_dlog/2 just takes a few
  // seconds. Something is indeed broken!
  mpz_init_set_str(h, "174757796725793774999567027597295197505665569273004209469970148397702384019044105191366253937185396858958956518015152796877013077655601176287549192437699162266404350225879496567141277628891297637669709009111308095675577974172825102201098090999441627405337875479851790438923887761011101214465827045413568965013", 10);
  mpz_init_set_str(g, "7", 10);
  mpz_init_set_str(p, "179769313486231590772930519078902473361797697894230657273430081157732675805500963132708477322407536021120113879871393357658789768814416622492847430639474124377767893424865485276302219601246094119453082952085005768838150682342462881473913110540827237163350510684586298239947245938479716304835356329624223043879", 10);

  mpz_init(pminus1);
  mpz_sub_ui(pminus1, p, 1);

  /* Prepare random number generator and seed by time */

  time_t timeseed;
  time(&timeseed);
  gmp_randinit_default(state);
  gmp_randseed_ui(state, timeseed);

  /* Compute, store and display logarithm */

  mpz_t result;
  mpz_init(result);
  pollardrho(&result, h, g, p);
  gmp_printf("Logarithm of %Zd to base %Zd in Z*(%Zd) is %Zd.\n",
             h, g, p, result);
}

/*********************************
 * Pollard’s rho-method function *
 *********************************/

void pollardrho(mpz_t * result, mpz_t h, mpz_t g, mpz_t p) {

  /* Randomly set initial sequence terms */

  mpz_t a1, b1, x1, a2, b2, x2, temp;
  mpz_init(a1);
  mpz_init(b1), mpz_init(x1);
  mpz_init(temp);
  mpz_urandomm(a1, state, p);
  mpz_urandomm(b1, state, p);
  mpz_powm(x1, g, a1, p);
  mpz_powm(temp, h, b1, p);
  mpz_mul(x1, x1, temp);
  mpz_mod(x1, x1, p);
  mpz_clear(temp);
  mpz_init(a2);
  mpz_init(b2);
  mpz_init(x2);
  mpz_set(a2, a1);
  mpz_set(b2, b1);
  mpz_set(x2, x1);

  /* Use Floyd’s algorithm to find sequence colliison */

  do {
    iterate(&x1, &a1, &b1);
    iterate(&x2, &a2, &b2);
    iterate(&x2, &a2, &b2);
  } while (mpz_cmp(x1, x2) != 0);

  /* Construct linear congruence */

  mpz_t a, b;    /* a(log) = b (mod p-1) */
  mpz_init(a);
  mpz_init(b);
  mpz_sub(b, a2, a1);
  mpz_mod(b, b, pminus1);
  mpz_sub(a, b1, b2);
  mpz_mod(a, a, pminus1);
  mpz_clear(a1);
  mpz_clear(a2);
  mpz_clear(b1);
  mpz_clear(b2);

  /* Ensure congruence is non-trivial and retry if not */

  if (mpz_cmp_ui(a, 0) == 0) {
    pollardrho(result, h, g, p);
  }

  /* Solve linear congruence if possible and retry otherwise */

  mpz_t gcd, soln;
  mpz_init(gcd);
  mpz_init(soln);
  mpz_gcdext(gcd, soln, NULL, a, pminus1);
  if (mpz_divisible_p(b, gcd)) {
    /* Congruence has a solution, solve it */
    mpz_mul(soln, b, soln);
    mpz_div(soln, soln, gcd);
    mpz_mod(soln, soln, pminus1);
  } else {
    /* Congruence has no solution, try again */
    pollardrho(result, h, g, p);
  }

  /* Perform trial exponentiation on candidates */

  mpz_t step, trial;
  mpz_init(step);
  mpz_init(trial);
  mpz_div(step, pminus1, gcd);
  mpz_mod(soln, soln, step);

  while (1) {
    mpz_powm(trial, g, soln, p);
    if (mpz_cmp(trial, h) == 0)
      break;
    else
      mpz_add(soln, soln, step);
  }

  /* Store result */

  mpz_set(*result, soln);
}

/************************************
 * Rho sequence updating subroutine *
 ************************************/

void iterate(mpz_t * x, mpz_t * a, mpz_t * b)
{
  if (mpz_congruent_ui_p(*x, 0, 3)) {
    mpz_mul(*x, *x, *x);
    mpz_mod(*x, *x, p);
    mpz_mul_ui(*a, *a, 2);
    mpz_mod(*a, *a, pminus1);
    mpz_mul_ui(*b, *b, 2);
    mpz_mod(*b, *b, pminus1);
  } else if (mpz_congruent_ui_p(*x, 1, 3)) {
    mpz_mul(*x, g, *x);
    mpz_mod(*x, *x, p);
    mpz_add_ui(*a, *a, 1);
    mpz_mod(*a, *a, pminus1);
  } else if (mpz_congruent_ui_p(*x, 2, 3)) {
    mpz_mul(*x, h, *x);
    mpz_mod(*x, *x, p);
    mpz_add_ui(*b, *b, 1);
    mpz_mod(*b, *b, pminus1);
  }
}
