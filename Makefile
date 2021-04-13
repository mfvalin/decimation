
FC      = gfortran
FFLAGS  = 
CC      = gcc
CFLAGS  = -O3 -Wall
DEFINES =

all:	decimate.c.o decimate.Abs

decimate.c.o: decimate.c
	$(CC) -I. -c $(CFLAGS) $< -o $@

decimate.Abs: decimate.c
	$(CC) -I. -DSELF_TEST $(DEFINES) $(CFLAGS) $< -o $@

test:	decimate.Abs
	./$<

test_ftn: decimate_test.F90 decimate.c.o
	$(FC) -I. $(FFLAGS)  $(DEFINES) decimate_test.F90 decimate.c.o -o decimate_ftn.Abs
	./decimate_ftn.Abs

clean:
	rm -f *.o *.Abs a.out
