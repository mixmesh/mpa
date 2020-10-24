#@BEGIN-C_SRC-DEFAULT-RULES@
ERL=$(shell readlink -f `which erl`)
ERL_TOP=$(ERL:%/bin/erl=%)
OS=$(shell uname -s)
CFLAGS=-I$(ERL_TOP)/usr/include
#@END-C_SRC-DEFAULT-RULES@

CFLAGS+=-finline-functions -std=c99 -Wall -Werror -O # -DDEBUG
LDFLAGS+=-shared

OBJS = gmp_nif.o dloglib.o erl_nif_big.o

all: ../priv/gmp_nif.so

../priv/gmp_nif.so: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) $(ERL_TOP)/usr/lib/libei.a -lgmp -o ../priv/gmp_nif.so

%.o:	%.c
	$(CC) $(CFLAGS) -c -fPIC $<

# For test purposes
dlog: dlog.c
	$(CC) -std=c99 -Wall -Werror dlog.c -lgmp -o dlog

clean:
	rm -rf *.o ../priv/gmp_nif.so dlog
