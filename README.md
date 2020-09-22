# Multiple precision arithmetic library 

```
./lib/sstdlib/src/mpz.erl
./lib/sstdlib/src/gmp_nif.erl
./lib/sstdlib/c_src/gmp_nif.c
./lib/sstdlib/c_src/dloglib.c
```

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
