# Multiple precision arithmetic

An Erlang NIF library to provide access to https://gmplib.org/ with additional functionality to solve the discrete logarithm.

## Files

<dl>
  <dt>./src/mpz.erl</dt>
  <dd>The API module</dd>
  <dt>./src/gmp_nif.erl</dt>
  <dd>The Erlang part of the NIF adaption, i.e. used by mpz</dd>
  <dt>./src/gmp_nif.c</dt>
  <dd>The C part of the NIF adaption used by the gmp_nif module</dd>
  <dt>./src/dloglib.c</dt>
  <dd>Additional C code to solve the discrete logarithm using the Pollard's rho-method</dd>
</dl>

## Additional information

The mpz module makes a number of functions available from GMP:

* generate_safe_prime
* mpz_gcd
* mpz_invert
* mpz_lcm
* mpz_powm
* mpz_pow_ui
* mpz_probab_prime_p
 
Additionally a dlog function is made available to perform fast
calculation of the discrete logarithm using Pollard’s rho-method as
described in https://www.luke.maurits.id.au/files/misc/honours_thesis.pdf

Note: To debug the dlog stuff use the standalone lib/sstdlib/c_src/dlog.c
