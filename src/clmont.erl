%%% @author Tony Rogvall <tony@rogvall.se>
%%% @copyright (C) 2021, Tony Rogvall
%%% @doc
%%%    Run mongomery operations on openCL
%%% @end
%%% Created : 12 Feb 2021 by Tony Rogvall <tony@rogvall.se>

-module(clmont).

-compile(export_all).

-include_lib("cl/include/cl.hrl").
-include("mont.hrl").


build() ->
    E = clu:setup(all),
    Filename = filename:join(code:priv_dir(mpa),"clmont.cl"),
    Include  = filename:join(code:lib_dir(mpa),"c_src"),
    case clu:build_source_file(E,Filename,"-I"++Include) of
	{error,{[{ok,error}],
		[Message]}} ->
	    io:put_chars(Message),
	    error;
	Res ->
	    Res
    end.

-define(P_1024, 1191703890297837857254846218124820162520314254482239260141586246493315566589245659462156276340012962327654624865776671922725912417154643528357403702766406672783187741039499777500937664819366321506835371609274218842538110523885904400885445461904752292635899168049169243216400297218378136654191604761801220538347).

mont() ->
    %% N = (1 bsl 63)+13,
    N = ?P_1024,
    mpz:mont_w(fips, N, 32).

format() ->
    mpz:format_mont(mont()).

%% Run 1000 pow operations

test() ->
    test(all).

test(DevType) ->
    test(DevType, 1000).

test(DevType, Count) ->
    M = mont(),

    As = [random_a(M) || _ <- lists:seq(1, Count)],
    Ams = [mpz:to_mont(Ai,M) || Ai <- As],
    AsData = list_to_binary([encode_number(Am,32,M#mont.k) || Am <- Ams]),

    X = random_x(M),
    XBits = gmp_nif:big_bits(X),
    XData = encode_number(X, 32, XBits),
    case run(AsData,XData,XBits,DevType,Count) of
	{ok,ResData} ->
	    Rs = decode_numbers(ResData,32,M#mont.k),
	    %% FIXME: does not return correct result!
	    lists:foreach(
	      fun({Rm,Ai}) ->
		      io:format("Ai=~w, X=~w\n", [Ai,X]),
		      io:format("Rm=~w\n", [Rm]),
		      Ri = mpz:from_mont(Rm,M),
		      io:format("Ri=~w\n", [Ri]),
		      R = mpz:powm(Ai, X, M#mont.n),
		      io:format("R=~w\n", [R])
	      end, lists:zip(Rs, As)),
	    Rs;
	Error ->
	    Error
    end.

run(AsData,XData,XBits,DevType,Count) ->
    E = clu:setup(DevType),
    io:format("platform created\n"),

    Filename = filename:join(code:priv_dir(mpa),"clmont.cl"),
    Include  = filename:join(code:lib_dir(mpa),"c_src"),

    io:format("build: ~s\n", [Filename]),
    {ok,Program} = clu:build_source_file(E, Filename, "-I"++Include),
    io:format("program built\n"),
    
    N = byte_size(AsData),  %% match!

    %% Create input/output data memory (implicit copy_host_ptr)
    {ok,Amem} = cl:create_buffer(E#cl.context,[read_write],N),
    io:format("input Amem memory created\n"),

    %% Create input data memory (implicit copy_host_ptr)
    NX = byte_size(XData),
    {ok,Xmem} = cl:create_buffer(E#cl.context,[read_only],NX),
    io:format("input Xmem memory created\n"),

    %% Create the command queue for the first device
    {ok,Queue} = cl:create_queue(E#cl.context,hd(E#cl.devices),[]),
    io:format("queue created\n"),

    %% Create the squre kernel object
    {ok,Kernel} = cl:create_kernel(Program, "pow"),
    io:format("kernel created: ~p\n", [Kernel]),

    %% dump_data(AsData),

    %% Write data into input array 
    {ok,Event1} = cl:enqueue_write_buffer(Queue, Amem, 0, N, AsData, []),
    io:format("write ~w bytes to Amem data enqueued\n", [N]),

    {ok,Event2} = cl:enqueue_write_buffer(Queue, Xmem, 0, NX, XData, []),
    io:format("write ~w bytes X data enqueued\n", [NX]),

    %% Set kernel arguments
    clu:apply_kernel_args(Kernel, [Amem,Xmem,{uint,XBits},{uint,Count}]),
    io:format("kernel args set\n"),

    Device = hd(E#cl.devices),
    {ok,MaxUnits} = cl:get_device_info(Device, max_compute_units),
    io:format("max_compute_units = ~w\n", [MaxUnits]),
    {ok,WGSize} = cl:get_kernel_workgroup_info(Kernel, Device, work_group_size),

    EGlobal = global_size(Count,WGSize),
    ELocal = local_size(EGlobal,WGSize),

    T0 = erlang:monotonic_time(),

    {ok,Event3} = cl:enqueue_nd_range_kernel(Queue, Kernel,
					     [EGlobal], [ELocal], 
					     [Event1,Event2]),
    io:format("nd range [~w, ~w] kernel enqueued\n",
	      [[EGlobal],[ELocal]]),
    
    %% Enqueue the read from device memory (wait for kernel to finish)
    {ok,Event4} = cl:enqueue_read_buffer(Queue,Amem,0,N,[Event3]),
    io:format("read buffer enqueued\n"),

    %% Now flush the queue to make things happend 
    ok = cl:flush(Queue),
    io:format("flushed\n"),

    %% Wait for Result buffer to be written
    io:format("wait\n"),
    io:format("Event1 = ~p\n", [cl:wait(Event1,1000)]),
    io:format("Event2 = ~p\n", [cl:wait(Event2,1000)]),
    Res3 = cl:wait(Event3,10000),
    T1 = erlang:monotonic_time(),
    io:format("Event3 = ~p\n", [Res3]),
    Time = erlang:convert_time_unit(T1-T0,native,microsecond),
    TimeS = Time/1000000,
    PPS = Count/TimeS,
    io:format("~s: time=~f,  #pow/s = ~f\n", [fips,TimeS,PPS]),

    Res = cl:wait(Event4,10000),

    %%
    cl:release_mem_object(Amem),
    cl:release_mem_object(Xmem),
    cl:release_queue(Queue),
    cl:release_kernel(Kernel),
    cl:release_program(Program),

    clu:teardown(E),
    Res.


local_size(N, WorkGroupSize) when N > WorkGroupSize -> WorkGroupSize;
local_size(N, _WorkGroupSize) -> N.

global_size(N, WorkGroupSize) when N > WorkGroupSize -> 
    ((N+WorkGroupSize+1) div WorkGroupSize)*WorkGroupSize;
global_size(N, _WorkGroupSize) -> N.


random_a(M) ->
    %% P = M#mont.n,
    %% uniform(0, P-1).
    %% uniform(1,10).
    1.

random_x(M) ->
    %% Q = (M#mont.n - 1) div 2,
    %% uniform(1, Q).
    60.

uniform(Min, Max) ->
    Min1 = Min - 1,
    N = Max-Min1,
    R = rand:uniform(N),
    R+Min1.

%% encode a K bit number in W bit digits
encode_number(N, W, K) ->
    encode_number_(N, W, K, []).

encode_number_(_N, _W, K, Acc) when K =< 0 ->
    list_to_binary(lists:reverse(Acc));
encode_number_(N, W, K, Acc) ->
    encode_number_(N bsr W, W, K-W, [<<N:W/native>>|Acc]).

%% decode numbers
decode_numbers(Data, W, K) ->
    Size = (((K+W-1) div W)*W) div 8,  %% byte size
    decode_numbers_(Data, Size, W, K, []).

decode_numbers_(Data, Size, W, K, Acc) ->
    case Data of
	<<>> ->
	    lists:reverse(Acc);
	<<Bin:Size/binary,Data1/binary>> ->
	    decode_numbers_(Data1,Size,W,K,[decode_number(Bin,W)|Acc])
    end.

decode_number(Bin, W) ->
    decode_number(Bin, W, 0, 0).

decode_number(<<>>, _W, _Shift, Num) ->
    Num;
decode_number(Bin, W, Shift, Num) ->
    <<N:W/native,Bin1/binary>> = Bin,
    decode_number(Bin1, W, Shift+W, Num bor (N bsl Shift)).
    
