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
-export([big_mod2_sqr/2]).
-export([big_mont_mul/5, big_mont_mul/3]).
-export([big_mont_sqr/4, big_mont_sqr/2]).
-export([big_mont_pow/6, big_mont_pow/3]).
-export([mont/1, mont/2, mont_w/2, mont_w/3, redc/2]).
-export([mont_w_fips/1, mont_w_fips/2]).

-export([to_mont/2, from_mont/2]).
-export([to_mont/3, from_mont/3]).
-export([format_mont/1,format_mont/2]).
-export([format_cnum/3]).
%% test
-export([all/0]).
-export([test1/0, test3/0, test4/0]).
-export([test11/0, test13/0, test14/0]).
-export([test21/0, test23/0, test24/0, test25/0]).
-export([test26/0]).
-export([test1/1, test3/1, test4/1]).
-export([test11/1, test13/1, test14/1]).
-export([test21/1, test23/1, test24/1, test25/1]).
-export([test26/1]).

-export([mont_pow/4, mont_mul/4, mont_sqr/3]).
-export([test_sqr_random/2, test_sqr_random/3]).
-export([test_mul_random/2, test_mul_random/3]).
-export([test_pow_random/2, test_pow_random/3]).
-export([fips_accumulator_size/1, fips_accumulator_size/2]).

-include("mont.hrl").

-define(is_odd(P), (((P) band 1) =:= 1)).

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

big_mod2_sqr(A, K) -> gmp_nif:big_mod2_sqr(A, K).
    
big_mont_mul(Am,Bm,#mont{meth=Meth,n=N,np=Np}) ->
    big_mont_mul(Meth,Am,Bm,N,Np).
big_mont_mul(Meth,Am,Bm,N,Np) ->
    gmp_nif:big_mont_mul(Meth,Am,Bm,N,Np).

big_mont_sqr(Am, #mont{meth=Meth,n=N,np=Np}) ->
    big_mont_sqr(Meth,Am,N,Np).
big_mont_sqr(Meth,Am,N,Np) ->
    gmp_nif:big_mont_sqr(Meth,Am,N,Np).

big_mont_pow(Am, X, #mont{meth=Meth,m1=M1,n=N,np=Np}) ->
    big_mont_pow(Meth,Am,X,M1,N,Np).
big_mont_pow(Meth,Am,X,M1,N,Np) ->
    gmp_nif:big_mont_pow(Meth,Am,X,M1,N,Np).

%% Exported: pow_ui

pow_ui(Base, Exp) ->
  binary:decode_unsigned(
    gmp_nif:mpz_pow_ui(binary:encode_unsigned(Base), Exp)).

%% Exported: probab_prime_p

probab_prime_p(N, Reps) ->
  gmp_nif:mpz_probab_prime_p(binary:encode_unsigned(N), Reps).

%% calculate mongomery paramters Ri,Np
-define(P_1024, 1191703890297837857254846218124820162520314254482239260141586246493315566589245659462156276340012962327654624865776671922725912417154643528357403702766406672783187741039499777500937664819366321506835371609274218842538110523885904400885445461904752292635899168049169243216400297218378136654191604761801220538347).

mont(N) ->
    mont(?defmeth,N).

-spec mont(Meth::mont_meth(), N::unsigned()) -> #mont{}.
mont(Meth,N) when is_integer(N), N>0, ?is_odd(N) ->
    W = 8*erlang:system_info(wordsize),
    mont_w(Meth, N, W).

mont_w_fips(W) ->
    mont_w(fips,?P_1024,W).

mont_w_fips(N, W) ->
    mont_w(fips,N,W).

mont_w(Meth,W) ->
    mont_w(Meth,?P_1024,W).
mont_w(Meth,N,W) when is_integer(N), N>0, ?is_odd(N) ->
    Ns = big_bits(N),
    K = ((Ns + W - 1) div W)*W,
    R = (1 bsl K),
    {1,{S,T}} = egcd(R, N),
    Ri = mod(S, N),
    Np = mod(-T, R),
    M1 = (1 bsl K) rem N,
    #mont{meth=Meth,k=K,m1=M1,n=N,np=Np,ri=Ri}.

%% given S number of words
fips_accumulator_size(S) ->
    W = 8*erlang:system_info(wordsize),
    fips_accumulator_size(W, S).

fips_accumulator_size(W, S) ->
    trunc(math:log(W*(W-1)*S) / math:log(W)).

to_mont(X, #mont{k=K,n=N}) -> to_mont(X, K, N).
to_mont(X, K, N) -> (X bsl K) rem N.

from_mont(Y, #mont{n=N,ri=Ri}) -> from_mont(Y,Ri,N).
from_mont(Y,Ri,N) -> (Y*Ri) rem N.

%% Output current mont contants i C code format
%% to be used by C/openCL etc
%% W must be equal to sizeof(UINT_T) on target
%% W is typically 32 for openCL (int), and 32 or 64 for C (int or long)
format_mont(M) -> format_mont(M,32).
format_mont(M,W) ->
    S = (M#mont.k+W-1) div W,
    io:format("#include \"~s.i\"\n", [M#mont.meth]),
    io:format("CONST int mont_K = ~w;\n", [M#mont.k]),
    io:format("CONST int mont_S = ~w;\n", [S]),
    io:format("CONST UINT_T mont_N[~w] = ~s;\n", 
	      [S,format_cnum(W,M#mont.k,M#mont.n)]),
    io:format("CONST UINT_T mont_Np[~w] = ~s;\n", 
	      [S,format_cnum(W,M#mont.k,M#mont.np)]),
    io:format("CONST UINT_T mont_1[~w] = ~s;\n", 
	      [S,format_cnum(W,M#mont.k,M#mont.m1)]),
    ok.

format_cnum(W,K,N) ->
    Ds = format_cnum_(W,(1 bsl W)-1,K,N,[]),
    DsL = ["0x"++tl(integer_to_list(D+(1 bsl W), 16)) || D <- Ds],
    ["{",string:join(DsL,","),"}"].

format_cnum_(_W, _Wm, K, _N, Acc) when K =< 0 ->
    lists:reverse(Acc);
format_cnum_(W, Wm, K, N, Acc) ->
    format_cnum_(W, Wm, K-W, N bsr W, [(N band Wm) | Acc]).

%%  Montgomery reduction
redc(T, #mont{k=K,n=N,np=Np}) ->
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

%% invert(X,P) when is_integer(X), X > 0, 
%% 		 is_integer(P), P > 0, X < P ->
%%    case egcd(X, P) of
%%	{1,{A,_B}} when A < 0 -> P + A;
%%	{1,{A,_B}} -> A;
%%	_ ->
%%	    erlang:error(not_defined)
%%    end.

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


ipowm(A, B, M) ->
    ipowm_(A rem M, B, M, 1).

ipowm_(A, 1, M, P) ->
     (A*P) rem M;
ipowm_(A, B, M, P)  ->
    B1 = B bsr 1,
    A1 = sqr(A) rem M,
    if B - B1 =:= B1 ->
 	    ipowm_(A1, B1, M, P);
       true ->
 	    ipowm_(A1, B1, M, (A*P) rem M)
    end.


sqr(A) when is_integer(A) ->
    A*A.

%% montgomery test

-define(assert(Value, Expr),
	case (Expr) of
	    Value -> ok;
	    _OtherValue ->
		io:format("FAIL: ~s != ~w\n", [??Expr, Value]),
		throw(fail)
	end).

all() ->
    Tests =
	[fun test1/0, fun test3/0, fun test4/0,
	 fun test11/0, fun test13/0, fun test14/0,
	 fun test21/0, fun test23/0, fun test24/0,
	 fun test25/0, fun test26/0 ],
    lists:foreach(
      fun(Fun) ->
	      try Fun() of
		  _ -> ok
	      catch
		  throw:fail ->
		      error
	      end
      end, Tests).


test1() ->
    test1(?defmeth).

test1(Meth) when ?is_meth(Meth) ->
    io:format("test1: ~w ", [Meth]),
    P = 101,
    M = mont(Meth,P),
    A = 79,
    Am = to_mont(A, M),
    B = 33,
    Bm = to_mont(B, M),
    Rm = redc(Am*Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem P),
    %% ?assert(C, 82),
    io:format("OK\n"),
    ok.

test3() ->
    [test3(M) || M <- ?allmeth].
test3(Meth) when ?is_meth(Meth) ->
    io:format("test3: ~w ", [Meth]),
    P = 101,
    M = mont(Meth,P),
    A = 79,
    Am = to_mont(A, M),
    B = 33,
    Bm = to_mont(B, M),
    %% io:format("P=~w, A=~w,Am=~w, B=~w,Bm=~w\n", [P,A,Am,B,Bm]),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem P),
    %% ?assert(C, 82),
    io:format("OK\n"),
    ok.    

test4() ->
    [test4(M) || M <- ?allmeth].
test4(Meth) when ?is_meth(Meth) ->
    io:format("test4: ~w ", [Meth]),
    P = 101,
    M = mont(Meth,P),
    A = 79,
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*A) rem P),
    %% ?assert(C, 80),
    io:format("OK\n"),
    ok.    

test11() ->
    test11(?defmeth).
test11(Meth) when ?is_meth(Meth) ->
    io:format("test11: ~w ", [Meth]),
    P = 1262773213764120865151395821008507246189,
    M = mont(Meth,P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = redc(Am*Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem P),
    %% ?assert(C, 3060820620989551345058379044987056313),
    io:format("OK\n"),
    ok.

test13() ->
    [test13(M) || M <- ?allmeth].
test13(Meth) when ?is_meth(Meth) ->
    io:format("test13: ~w ", [Meth]),
    P = 1262773213764120865151395821008507246189,
    M = mont(Meth,P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    B = 1491922486076647757424410593223,
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem P),
    io:format("OK\n"),
    ok.    

test14() ->
    [test14(M) || M <- ?allmeth].
test14(Meth) when ?is_meth(Meth) ->
    io:format("test14: ~w ", [Meth]),
    P = 1262773213764120865151395821008507246189,
    M = mont(Meth,P),
    A = 2209866513432185383910552416615,
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*A) rem P),
    io:format("OK\n"),
    ok.    

-define(P, ((1 bsl 1024) - 1093337)).
-define(G, 7).
-define(X, ((1 bsl 512) div 5)).

test21() ->
    test21(?defmeth).
test21(Meth) when ?is_meth(Meth) ->
    io:format("test21: ~w ", [Meth]),
    M = mont(Meth,?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    Am = to_mont(A, M),
    B = 124060833865307543493568060615888984483974452188008593376614358923889775329311551928017371949352869174793233582962966084426368046363097562383329228081040634647665584981130599395343295428306766898340404601804058174196606881036706450931509901449845820591347944486874994264576959961749516461160241513031362543684,
    Bm = to_mont(B, M),
    Rm = redc(Am*Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem ?P),
    io:format("OK\n"),
    ok.

test23() ->
    [test23(M) || M <- ?allmeth].
test23(Meth) when ?is_meth(Meth) ->
    io:format("test23: ~w ", [Meth]),
    M = mont(Meth,?P),
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    %% io:format("log2(A) = ~w\n", [imath:ilog2(A)]),
    Am = to_mont(A, M),
    B = 124060833865307543493568060615888984483974452188008593376614358923889775329311551928017371949352869174793233582962966084426368046363097562383329228081040634647665584981130599395343295428306766898340404601804058174196606881036706450931509901449845820591347944486874994264576959961749516461160241513031362543684,
    %% io:format("log2(B) = ~w\n", [imath:ilog2(B)]),
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem ?P),
    io:format("OK\n"),
    ok.    

test24() ->
    [test24(M) || M <- ?allmeth].
test24(Meth) when ?is_meth(Meth) ->
    io:format("test24: ~w ", [Meth]),
    Mp = mont(Meth,?P),
    %% A is 1024 bits
    A = 129015633621243001913449155039089342460245014035423996408293170177875331730667287169606301577788634719390344399495012523585332459829150811854319143390944179836371880331297042928788452376675529244126810320465199234812830335213501532015809681469166815018523250762130729029190848536853958461583303676087231435574,
    %% io:format("log2(A) = ~w\n", [imath:ilog2(A)]),
    Am = to_mont(A, Mp),
    Rm = big_mont_sqr(Am, Mp),
    C = from_mont(Rm, Mp),
    ?assert(C, (A*A) rem Mp#mont.n),
    io:format("OK\n"),
    ok.

test25() ->
    [test25(Meth) || Meth <- ?allmeth].

test25(Meth) when ?is_meth(Meth) ->
    io:format("test25: ~w ", [Meth]),
    Mp = mont(Meth,?P),
    Gm = to_mont(?G, Mp),
    test25_(Mp, Gm, 999).

test25_(_Mp, _Gm, X) when X > 1000 -> 
    io:format("OK\n"),
    ok;
test25_(Mp, Gm, X) ->
    %% io:format("X = ~w\n", [X]),
    Rm = big_mont_pow(Gm, X, Mp),
    C = from_mont(Rm, Mp),
    ?assert(C, ipowm(?G, X, Mp#mont.n)),
    test25_(Mp, Gm, X+1).

%% test vectors from:
%% https://www.researchgate.net/publication/4107322_Montgomery_modular_multiplication_architecture_for_public_key_cryptosystems

test26() ->
    [test26(Meth) || Meth <- ?allmeth].

test26(Meth) ->
    io:format("test26: ~w ", [Meth]),
    test26_1(Meth),
    test26_2(Meth),
    test26_3(Meth),
    test26_4(Meth),
    test26_5(Meth),
    io:format("OK\n").
    
test26_1(Meth) ->
    P  = 16#eee74404d129949520704c5bf5814703,
    M = mont(Meth, P),
    Am  = 16#223520375bd184b2bac64c9d1a6c55fa,
    Bm  = 16#d857ed0720d590d61f05c150e1e40917,
    Cm = big_mont_mul(Am, Bm, M),
    Cm = 16#e3bd635debc8021ea0208d75df078ea6,
    ok.

test26_2(Meth) ->
    P  = 16#f0fe9f3c608d779379bb3676fdb85071,
    M = mont(Meth, P),
    Am  = 16#adfdb6089758064a69aad900ad18274b,
    Bm  = 16#828994fe60ddcab48671399fd349b0ff,
    Cm = big_mont_mul(Am, Bm, M),
    ?assert(Cm, 16#d848da961b4c3092def8bdaca7a73bee),
    ok.

test26_3(Meth) ->
    P  = 16#dd54fd6aa41a0dbcb7550b284862b7a5,
    M = mont(Meth, P),
    Am  = 16#c9d5398bf1cec1342f3c7cca33ea04a6,
    Bm  = 16#2f5bf07df1f473c46318eea49d7f1de3,
    Cm = big_mont_mul(Am, Bm, M),
    ?assert(Cm, 16#7dccdc7deb3d1a1b3afb2b7c0ca5c53a),
    ok.    

test26_4(Meth) ->
    P  = 16#9848765c71a41c12921ca63ca42a203d,
    M = mont(Meth, P),
    Am  = 16#c52c1e570e620174fcbb51063ecbfcb1,
    Bm  = 16#4b1a9caa8e183aabd89e8f3ab7561273,
    Cm = big_mont_mul(Am, Bm, M),
    ?assert(Cm,16#55352eb5475ce94fefcce0a8b687044f),
    ok.

test26_5(Meth) ->
    P  = 16#c71a6ffca861df3a175c6eb226581289,
    M = mont(Meth, P),
    Am  = 16#b253687d5f73bacf1fbd542d5d604272,
    Bm  = 16#af9764ea40aa14153e8af13177851ab1,
    Cm = big_mont_mul(Am, Bm, M),
    ?assert(Cm, 16#8b89addf457e084b752d716c73d29821),
    ok.


test_mul_random(N, Range) -> 
    [test_mul_random(M, N, Range) || M <- ?allmeth].
test_mul_random(Meth, N, Range) -> 
    io:format("test_mul_random: ~s ", [Meth]),
    rand:seed(exsss, erlang:system_time()),
    M = mont(Meth, ?P),
    try test_mul_random_(N, Range, M) of
	ok -> io:format("OK\n")
    catch
	throw:fail ->
	    io:format("FAIL\n")
    end.


test_mul_random_(0, _Range, _M) -> 
    ok;
test_mul_random_(I, Range, M) ->
    A = random(Range, M#mont.n),
    Am = to_mont(A, M),
    B = random(Range, M#mont.n),
    Bm = to_mont(B, M),
    Rm = big_mont_mul(Am, Bm, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*B) rem M#mont.n),
    test_mul_random_(I-1, Range, M).

test_sqr_random(N, Range) -> 
    [test_sqr_random(M, N, Range) || M <- ?allmeth].
test_sqr_random(Meth, N, Range) -> 
    io:format("test_sqr_random: ~s ", [Meth]),
    rand:seed(exsss, erlang:system_time()),
    M = mont(Meth, ?P),
    try test_sqr_random_(N, Range, M) of
	ok -> io:format("OK\n")
    catch
	throw:fail ->
	    io:format("FAIL\n")
    end.

test_sqr_random_(0, _Range, _M) -> 
    ok;
test_sqr_random_(I, Range, M) ->
    A = random(Range, M#mont.n),
    Am = to_mont(A, M),
    Rm = big_mont_sqr(Am, M),
    C = from_mont(Rm, M),
    ?assert(C, (A*A) rem M#mont.n),
    test_sqr_random_(I-1, Range, M).


test_pow_random(N, Range) -> 
    [test_pow_random(M, N, Range) || M <- ?allmeth].
test_pow_random(Meth, N, Range) ->
    io:format("test_pow_random: ~s", [Meth]),
    rand:seed(exsss, erlang:system_time()),
    M = mont(Meth, ?P),
    try test_pow_random_(N, Range, M) of
	ok -> io:format(" OK\n")
    catch
	throw:fail ->
	    io:format(" FAIL\n")
    end.

test_pow_random_(0, _Range, _M) -> 
    ok;
test_pow_random_(I, Range, M) ->
    A = random(Range, M#mont.n),
    Am = to_mont(A, M),
    X = random(M#mont.k, M#mont.n),
    Rm = big_mont_pow(Am, X, M),
    C = from_mont(Rm, M),
    ?assert(C, ipowm(A, X, M#mont.n)),
    test_pow_random_(I-1, Range, M).

%% A*B mod P (using method Meth)
-spec mont_mul(Meth::mont_meth(),
	       A::unsigned(), B::unsigned(), P::unsigned()) ->
	  unsigned().

mont_mul(Meth, A, B, P) when ?is_meth(Meth), ?is_odd(P) ->
    M = mont(Meth, P),
    Am = to_mont(A, M),
    Bm = to_mont(B, M),
    Cm = big_mont_mul(Am, Bm, M),
    from_mont(Cm, M).

-spec mont_sqr(Meth::mont_meth(),A::unsigned(),
	       P::unsigned()) -> unsigned().
%% A^2 mod P (using method Meth)
mont_sqr(Meth, A, P) when ?is_meth(Meth), ?is_odd(P) ->
    M = mont(Meth, P),
    Am = to_mont(A, M),
    Cm = big_mont_sqr(Am, M),
    from_mont(Cm, M).

-spec mont_pow(Meth::mont_meth(),
	       A::unsigned(), N::unsigned(), P::unsigned()) ->
	  unsigned().

%% A^N mod P
mont_pow(Meth, A, N, P) when ?is_meth(Meth) ->
    M = mont(Meth,P),
    Am = to_mont(A, M),
    Rm = big_mont_pow(Am, N, M),
    from_mont(Rm, M).


random({Min,Max}, N) when Min =< Max ->
    (Min + rand:uniform((Max-Min)+1) - 1) rem N;
random(Max, N) ->
    (rand:uniform(Max) rem N).

