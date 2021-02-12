%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2020, Tony Rogvall
%%% @doc
%%%    Power wrench
%%% @end
%%% Created : 19 Sep 2020 by Tony Rogvall <tony@rogvall.se>

-module(pow_bench).

-compile(export_all).

-define(P, 1191703890297837857254846218124820162520314254482239260141586246493315566589245659462156276340012962327654624865776671922725912417154643528357403702766406672783187741039499777500937664819366321506835371609274218842538110523885904400885445461904752292635899168049169243216400297218378136654191604761801220538347).
%% -define(P, ((1 bsl 1024) - 1093337)).
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
pow_mont() ->
    pow_mont(sos),
    pow_mont(sps),
    pow_mont(cios),
    pow_mont(fips),
    pow_mont(fios),
    pow_mont(cihs),
    ok.

pow_mont(Type) ->
    pow_mont(Type,1000).

pow_mont(Type,Laps) ->
    Mp = mpz:mont(Type,?P),
    Gh = mpz:to_mont(?G, Mp),
    T0 = erlang:monotonic_time(),
    pow_mont_loop(Laps, Gh, ?X, Mp),
    T1 = erlang:monotonic_time(),
    Time = erlang:convert_time_unit(T1-T0,native,microsecond),
    TimeS = Time/1000000,
    PPS = Laps/TimeS,
    io:format("~s: time=~f,  #pow/s = ~f\n", [Type,TimeS,PPS]),
    [{time,Time},{time_avg,PPS}].

pow_mont_loop(0, _Gh, _X, _Mp) ->
    ok;
pow_mont_loop(I, Gh, X, Mp) ->
    mpz:big_mont_pow(Gh, ?X, Mp),
    pow_mont_loop(I-1, Gh, X, Mp).



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
