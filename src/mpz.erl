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
-export([big_mont_pow/4, big_mont_pow/3]).
-export([mont/1, redc/2]).
-export([to_mont/2, from_mont/2]).
-export([to_mont/3, from_mont/3]).
%% test
-export([test1/0, test2/0, test3/0]).
-export([test11/0, test12/0, test13/0]).
-export([mont_pow/3]).

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

big_mont_redc(T, {_K,N,Np,_Ri}) ->
    big_mont_redc(T, N, Np).
big_mont_redc(T, N, Np) ->
    gmp_nif:big_mont_redc(T, N, Np).

big_mont_mul(A, B, {_K,N,Np,_Ri}) ->
    big_mont_mul(A, B, N, Np).
big_mont_mul(A, B, N, Np) ->
    gmp_nif:big_mont_mul(A, B, N, Np).

big_mont_sqr(A, {_K,N,Np,_Ri}) ->
    big_mont_sqr(A, N, Np).
big_mont_sqr(A, N, Np) ->
    gmp_nif:big_mont_sqr(A, N, Np).

big_mont_pow(A, E, {_K,N,Np,_Ri}) ->
    big_mont_pow(A, E, N, Np).
big_mont_pow(A, E, N, Np) ->
    gmp_nif:big_mont_pow(A, E, N, Np).

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
    {K,N,Np,Ri}.

to_mont(X, {K,N,_Np,_Ri}) -> to_mont(X, K, N).
to_mont(X, K, N) -> (X bsl K) rem N.

from_mont(Y, {_K,N,_Np,Ri}) -> from_mont(Y,Ri,N).
from_mont(Y,Ri,N) -> (Y*Ri) rem N.

%%  Montgomery reduction
redc(T, {K,N,Np,_Ri}) ->
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
    C = mpz:from_mont(Rm, M),
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
    C = mpz:from_mont(Rm, M),
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
    C = mpz:from_mont(Rm, M),
    C = (A*B) rem P,
    C = 82,
    ok.    


test11() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = mpz:redc(Am*Bm, M),
    C = mpz:from_mont(Rm, M),
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
    C = mpz:from_mont(Rm, M),
    C = (A*B) rem P,
    C = 3060820620989551345058379044987056313,
    ok.    

test13() ->
    P = 1262773213764120865151395821008507246189,
    M = mont(P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = mpz:from_mont(Rm, M),
    C = (A*B) rem P,
    C = 3060820620989551345058379044987056313,
    ok.    

mont_pow(A, B, N) ->
    M = mont(N),
    Am = to_mont(A, M),
    Rm = big_mont_pow(Am, B, M),
    mpz:from_mont(Rm, M).

    
    
