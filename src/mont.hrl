-ifndef(__MONT_HRL__).
-define(__MONT_HRL__, true).

-type mont_meth() :: sos | sps | cios | fips | fios | cihs.
-type unsigned() :: non_neg_integer().

-define(defmeth, sos).
-define(allmeth, [sos,sps,cios,fips,fios,cihs]).
-define(imeth, [cios,fips,fios,cihs]).  %% integrated
-define(smeth, [sos,sps]).      %% separate

-define(is_meth(M), 
	((((M)=:=sos) orelse ((M)=:=sps)
	  orelse ((M)=:=cios) orelse ((M)=:=fips) orelse ((M)=:=fios)
	  orelse ((M)=:=cihs)))).

-record(mont,
	{
	 meth = ?defmeth :: mont_meth(),
	 k,   %% size of n in number of bignum digits
	 m1,  %% mont(1)  1 in montgomery space
	 n,   %% modulus
	 np,  %% -N^-1 (mod 2^k)
	 ri   %% R^-1 (mod n)  (r = 1 << k)
	}).
-endif.
