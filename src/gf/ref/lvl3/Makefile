
CC=gcc
CFLAGS= -O3 -std=gnu11 -Wall -march=native -Wno-missing-braces -Wno-logical-not-parentheses 
LDFLAGS=-lm
AR=ar rcs
RANLIB=ranlib

OBJECTS=objs/fp_p47441.o objs/fp.o objs/fp2.o objs/random.o

all: lib tests
	
objs/fp_p47441.o: fp_p47441.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) fp_p47441.c -o objs/fp_p47441.o
	
objs/fp.o: fp.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) fp.c -o objs/fp.o
	
objs/fp2.o: fp2.c
	@mkdir -p $(@D)
	$(CC) -c $(CFLAGS) fp2.c -o objs/fp2.o

objs/random.o: ../../../common/generic/randombytes_system.c
	$(CC) -c $(CFLAGS) ../../../common/generic/randombytes_system.c -o objs/random.o

lib: $(OBJECTS)
	rm -rf lib
	mkdir lib
	$(AR) lib/libtest.a $^
	$(RANLIB) lib/libtest.a

tests: lib
	$(CC) $(CFLAGS) -L./lib test/test_fp.c test/test_extras.c -ltest $(LDFLAGS) -o test_fp -lgmp
	$(CC) $(CFLAGS) -L./lib test/test_fp2.c test/test_extras.c -ltest $(LDFLAGS) -o test_fp2 -lgmp

check: tests

.PHONY: clean

clean:
	rm -rf *.req objs lib test_fp*

