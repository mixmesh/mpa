#@BEGIN-C_SRC-DEFAULT-RULES@
ERL=$(shell readlink -f `which erl`)
ERL_TOP=$(ERL:%/bin/erl=%)
OS=$(shell uname -s)
CFLAGS=-I$(ERL_TOP)/usr/include
#@END-C_SRC-DEFAULT-RULES@

CFLAGS+=-finline-functions -std=c99 -Wall -Werror -O # -DDEBUG
LDFLAGS+=-shared

all: ../priv/gmp_nif.so

../priv/gmp_nif.so: gmp_nif.o dloglib.o
	$(CC) $(LDFLAGS) gmp_nif.o dloglib.o $(ERL_TOP)/usr/lib/libei.a -lgmp -o ../priv/gmp_nif.so

gmp_nif.o: gmp_nif.c
	$(CC) $(CFLAGS) -c -fPIC gmp_nif.c -o gmp_nif.o

dloglib.o: dloglib.c
	$(CC) $(CFLAGS) -c -fPIC dloglib.c -o dloglib.o

# For test purposes
dlog: dlog.c
	$(CC) -std=c99 -Wall -Werror dlog.c -lgmp -o dlog

clean:
	rm -rf *.o ../priv/gmp_nif.so dlog
