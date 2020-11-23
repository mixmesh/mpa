-module(mpz).
-export([dlog/3,
         generate_safe_prime/1,
         invert/2,
         gcd/2,
         lcm/2,
         powm/3,
         pow_ui/2,
         probab_prime_p/2]).
-export([big_powm/3]).
-export([big_size/1, big_bits/1]).
-export([big_mont_redc/3, big_mont_redc/2]).
-export([big_mont_mul/4, big_mont_mul/3]).
-export([big_mont_sqr/3, big_mont_sqr/2]).
-export([big_mont_pow/5, big_mont_pow/3]).
-export([mont/1, redc/2]).
-export([to_mont/2, from_mont/2]).
-export([to_mont/3, from_mont/3]).
%% test
-export([test1/0, test2/0, test3/0, test4/0]).
-export([test11/0, test12/0, test13/0, test14/0]).
-export([test21/0, test22/0, test23/0, test24/0, test25/0]).
-export([mont_pow/3, mont_mul/3, mont_sqr/2]).
-export([test_sqr_random/2, test_mul_random/2, test_pow_random/2]).


%% Exported: dlog

dlog(H, G, P) ->
  binary:decode_unsigned(
    gmp_nif:dlog(binary:encode_unsigned(H),
                 binary:encode_unsigned(G),
                 binary:encode_unsigned(P))).

%% Exported: generate_safe_prime

generate_safe_prime(Len) ->
  binary:decode_unsigned(gmp_nif:generate_safe_prime(Len)).

%% Exported: invert

invert(Op1, Op2) ->
  binary:decode_unsigned(
    gmp_nif:mpz_invert(binary:encode_unsigned(Op1),
                       binary:encode_unsigned(Op2))).

%% Exported: gcd

gcd(Op1, Op2) ->
  binary:decode_unsigned(
    gmp_nif:mpz_gcd(binary:encode_unsigned(Op1),
                    binary:encode_unsigned(Op2))).

%% Exported: lcm

lcm(Op1, Op2) ->
  binary:decode_unsigned(
    gmp_nif:mpz_lcm(binary:encode_unsigned(Op1),
                    binary:encode_unsigned(Op2))).

%% Exported: powm

powm(Base, Exp, Mod) ->
  binary:decode_unsigned(
    gmp_nif:mpz_powm(binary:encode_unsigned(Base),
                     binary:encode_unsigned(Exp),
                     binary:encode_unsigned(Mod))).

big_powm(Base, Exp, Mod) ->
    gmp_nif:big_powm(Base, Exp, Mod).

big_size(X) -> gmp_nif:big_size(X).
big_bits(X) -> gmp_nif:big_bits(X).

big_mont_redc(T, {_K,_M1,N,Np,_Ri}) ->
    big_mont_redc(T, N, Np).
big_mont_redc(T, N, Np) ->
    gmp_nif:big_mont_redc(T, N, Np).

big_mont_mul(A, B, {_K,_M1,N,Np,_Ri}) ->
    big_mont_mul(A, B, N, Np).
big_mont_mul(A, B, N, Np) ->
    gmp_nif:big_mont_mul(A, B, N, Np).

big_mont_sqr(A, {_K,_M1,N,Np,_Ri}) ->
    big_mont_sqr(A, N, Np).
big_mont_sqr(A, N, Np) ->
    gmp_nif:big_mont_sqr(A, N, Np).

big_mont_pow(A, E, {_K,M1,N,Np,_Ri}) ->
    big_mont_pow(A, E, M1, N, Np).
big_mont_pow(A, E, P1, N, Np) ->
    gmp_nif:big_mont_pow(A, E, P1, N, Np).

%% Exported: pow_ui

pow_ui(Base, Exp) ->
  binary:decode_unsigned(
    gmp_nif:mpz_pow_ui(binary:encode_unsigned(Base), Exp)).

%% Exported: probab_prime_p

probab_prime_p(N, Reps) ->
  gmp_nif:mpz_probab_prime_p(binary:encode_unsigned(N), Reps).

%% calculate mongomery paramters Ri,Np
mont(N) when is_integer(N), N>0, N band 1 =:= 1 ->
    W = 8*erlang:system_info(wordsize),
    K = big_size(N)*W,
    R = (1 bsl K),
    {1,{S,T}} = egcd(R, N),
    Ri = mod(S, N),
    Np = mod(-T, R),
    M1 = (1 bsl K) rem N,
    {K,M1,N,Np,Ri}.

to_mont(X, {K,_M1,N,_Np,_Ri}) -> to_mont(X, K, N).
to_mont(X, K, N) -> (X bsl K) rem N.

from_mont(Y, {_K,_M1,N,_Np,Ri}) -> from_mont(Y,Ri,N).
from_mont(Y,Ri,N) -> (Y*Ri) rem N.

%%  Montgomery reduction
redc(T, {K,_M1,N,Np,_Ri}) ->
    R1 = ((1 bsl K)-1),
    M = ((T band R1)*Np) band R1,
    V = (T + M*N) bsr K,
    if V > N -> V - N;
       true -> V
    end.

mod(A,N) when A>0, N>0, A < N -> A;
mod(A,N) ->
    A1 = A rem N,
    if A1 < 0 -> 
	    A1 + N;
       true ->
	    A1
    end.

egcd(R, Q) when is_integer(R), is_integer(Q) ->
    R1 = abs(R),
    Q1 = abs(Q),
    if Q1 < R1 -> egcd_(Q1,R1,1,0,0,1);
       true -> egcd_(R1,Q1,0,1,1,0)
    end.

egcd_(0,Q,_,_,Q1,Q2) -> {Q, {Q2,Q1}};
egcd_(R,Q,R1,R2,Q1,Q2) ->
    D = Q div R,
    egcd_(Q rem R, R, Q1-D*R1, Q2-D*R2, R1, R2).

%% montgomery test

test1() ->
    P = 101,
    M = mont(P),
    A = 79,
    Am = to_mont(A, M),
    B = 33,
    Bm = to_mont(B, M),
    Rm = mpz:redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    C = 82,
    ok.

test2() ->
    P = 101,
    M = mont(P),
    A = 79,
    Am = to_mont(A, M),
    B = 33,
    Bm = to_mont(B, M),
    Rm = big_mont_redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    C = 82,
    ok.    

test3() ->
    P = 101,
    M = mont(P),
    A = 79,
    Am = to_mont(A, M),
    B = 33,
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    C = 82,
    ok.    

test4() ->
    P = 101,
    M = mont(P),
    A = 79,
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    C = (A*A) rem P,
    C = 80,
    ok.    

test11() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = mpz:redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    C = 3060820620989551345058379044987056313,
    ok.

test12() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = big_mont_redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    ok.    

test13() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem P,
    ok.    

test14() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    C = (A*A) rem P,
    ok.    

-define(P, ((1 bsl 1024) - 1093337)).
-define(G, 7).
-define(X, ((1 bsl 512) div 5)).

test21() ->
    M = mont(?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    Am = to_mont(A, M),
    B = 124060833865307543493568060615888984483974452188008593376614358923889775329311551928017371949352869174793233582962966084426368046363097562383329228081040634647665584981130599395343295428306766898340404601804058174196606881036706450931509901449845820591347944486874994264576959961749516461160241513031362543684,
    Bm = to_mont(B, M),
    Rm = mpz:redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem ?P,
    ok.

test22() ->
    M = mont(?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    Am = to_mont(A, M),
    B = 124060833865307543493568060615888984483974452188008593376614358923889775329311551928017371949352869174793233582962966084426368046363097562383329228081040634647665584981130599395343295428306766898340404601804058174196606881036706450931509901449845820591347944486874994264576959961749516461160241513031362543684,
    Bm = to_mont(B, M),
    Rm = big_mont_redc(Am*Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem ?P,
    ok.    

test23() ->
    M = mont(?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    Am = to_mont(A, M),
    B = 124060833865307543493568060615888984483974452188008593376614358923889775329311551928017371949352869174793233582962966084426368046363097562383329228081040634647665584981130599395343295428306766898340404601804058174196606881036706450931509901449845820591347944486874994264576959961749516461160241513031362543684,
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem ?P,
    ok.    


test24() ->
    M = mont(?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    C = (A*A) rem ?P,
    ok.

test25() ->
    M = mont(?P),
    Gm = to_mont(?G, M),
    test25_(M, Gm, 2).

test25_(_M, _Gm, X) when X > 1000 -> 
    ok;
test25_(M, Gm, X) ->
    %% io:format("X = ~w\n", [X]),
    Rm = big_mont_pow(Gm, X, M),
    C = from_mont(Rm, M),
    C = imath:pow(?G, X, ?P),
    test25_(M, Gm, X+1).


test_sqr_random(N, Range) -> 
    rand:seed(exsss, erlang:system_time()),
    M = mont(?P),
    test_sqr_random_(N, Range, M).

test_sqr_random_(0, _Range, _M) -> 
    ok;
test_sqr_random_(I, Range, M) ->
    A = random(Range),
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    C = (A*A) rem ?P,
    test_sqr_random_(I-1, Range, M).


test_mul_random(N, Range) -> 
    rand:seed(exsss, erlang:system_time()),
    M = mont(?P),
    test_mul_random_(N, Range, M).

test_mul_random_(0, _Range, _M) -> 
    ok;
test_mul_random_(I, Range, M) ->
    A = random(Range),
    Am = to_mont(A, M),
    B = random(Range),
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    C = (A*B) rem ?P,
    test_mul_random_(I-1, Range, M).

test_pow_random(N, Range) -> 
    rand:seed(exsss, erlang:system_time()),
    M = mont(?P),
    test_pow_random_(N, Range, M).

test_pow_random_(0, _Range, _M) -> 
    ok;
test_pow_random_(I, Range, M) ->
    A = random(Range),
    Am = to_mont(A, M),
    X = random(element(1,M)),
    Rm = big_mont_pow(Am, X, M),
    C = from_mont(Rm, M),
    C = imath:pow(A, X, ?P),
    test_pow_random_(I-1, Range, M).

%%
mont_mul(A, B, N) ->
    M = mont(N),
    Am = to_mont(A, M),
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    from_mont(Rm, M).

mont_sqr(A, N) ->
    M = mont(N),
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    from_mont(Rm, M).

mont_pow(A, B, N) ->
    M = mont(N),
    Am = to_mont(A, M),
    Rm = big_mont_pow(Am, B, M),
    from_mont(Rm, M).

random({Min,Max}) ->
    Min + rand:uniform((Max-Min)+1) - 1;
random(Max) ->
    rand:uniform(Max).
