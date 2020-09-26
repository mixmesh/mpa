%% -*- erlang -*-
{application, mpa,
 [{description,"Multiple precision arithmetic"},
  {vsn, "1.0"},
  {modules, [gmp_nif, mpz]},
  {applications, [kernel, stdlib]}]}.
