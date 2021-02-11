%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2020, Tony Rogvall
%%% @doc
%%%    Power wrench
%%% @end
%%% Created : 19 Sep 2020 by Tony Rogvall <tony@rogvall.se>

-module(pow_bench).

-compile(export_all).

-define(P, ((1 bsl 1024) - 1093337)).
-define(G, 7).
-define(X, ((1 bsl 512) div 5)).

pow_1() ->
    benchmark:run(fun (X) -> pow(?G, X, ?P) end, [?X]).

pow_2() ->
    benchmark:run(fun (X) -> binary:decode_unsigned(crypto:mod_pow(?G, X, ?P)) end, [?X]).

pow_3() ->
    benchmark:run(fun (X) -> mpz:powm(?G, X, ?P) end, [?X]).

pow_4() ->
    benchmark:run(fun (X) -> mpz:big_powm(?G, X, ?P) end, [?X]).

pow_5() -> pow_5(sos).
    
pow_5(Type) ->
    Mp = mpz:mont(Type,?P),
    Gh = mpz:to_mont(?G, Mp),
    benchmark:run(fun (X) -> mpz:big_mont_pow(Gh, X, Mp) end, [?X]).

%% special!
pow_5_bench() ->
    pow_5_bench(sos).

pow_5_bench(Type) ->
    pow_5_bench(Type,10000).
pow_5_bench(Type,Laps) ->
    Mp = mpz:mont(Type,?P),    
    Gh = mpz:to_mont(?G, Mp),
    T0 = erlang:monotonic_time(),    
    pow_5_loop(Laps, Gh, ?X, Mp),
    T1 = erlang:monotonic_time(),
    Time = erlang:convert_time_unit(T1-T0,native,microsecond),
    [{time,Time},{time_avg,Laps/(Time/1000000)}].

pow_5_loop(0, _Gh, _X, _Mp) ->
    ok;
pow_5_loop(I, Gh, X, Mp) ->
    mpz:big_mont_pow(Gh, ?X, Mp),
    pow_5_loop(I-1, Gh, X, Mp).



-define(is_non_neg_integer(X), (is_integer((X)) andalso ((X)>=0))).

%%
%% calculate A^B mod M
%%

-spec pow(A::non_neg_integer(), B::non_neg_integer(), M::non_neg_integer()) ->
		 non_neg_integer().

pow(0, B, M) when ?is_non_neg_integer(B),?is_non_neg_integer(M) -> 0;
pow(A, B, 0) when ?is_non_neg_integer(B),?is_non_neg_integer(A) -> 0;
pow(A, 0, M) when is_integer(A),?is_non_neg_integer(M) -> 1;
pow(1, B, M) when ?is_non_neg_integer(B),?is_non_neg_integer(M) -> 1;
pow(-1,B, M) when ?is_non_neg_integer(B),?is_non_neg_integer(M) ->
    (1 - 2*(B band 1));
pow(A, B, M) when is_integer(A),?is_non_neg_integer(B),?is_non_neg_integer(M) ->
    pow_(mod(A, M), B, M, 1).

pow_(A, 1, M, P) ->
    mod(A*P,M);
pow_(A, B, M, P)  ->
    B1 = B bsr 1,
    A1 = mod(A*A, M),
    if B - B1 =:= B1 ->
	    pow_(A1, B1, M, P);
       true ->
	    pow_(A1, B1, M, mod(A*P, M))
    end.

mod(X, N) when X < N -> X;
mod(X, N) -> X rem N.
