
FC      = gfortran
FFLAGS  = 
CC      = gcc
CFLAGS  = -O3 -Wall
DEFINES =

all:	decimate.c.o decimate.Abs decimate_ftn.Abs

decimate.c.o: decimate.c
	$(CC) -I. $(DEFINES) -c $(CFLAGS) $< -o $@

decimate.Abs: decimate.c
	$(CC) -I. -DSELF_TEST $(DEFINES) $(CFLAGS) $< -o $@

test:	decimate.Abs
	./$<
	rm -f ./$<

test_ftn: decimate_ftn.Abs
	./$<
	rm -f ./$<

decimate_ftn.Abs:	decimate_test.F90 decimate.c.o decimate_array.F90
	$(FC) -I. $(FFLAGS)  $(DEFINES) decimate_array.F90 decimate_test.F90 decimate.c.o -o decimate_ftn.Abs
	$(CC) -I. $(DEFINES) -c $(CFLAGS) -DDEBUG=1 decimate.c -o decimated.c.o
	$(FC) -I. $(FFLAGS)  $(DEFINES) -DDEBUG=1 decimate_array.F90 decimate_test.F90 decimated.c.o -o decimate_ftn.AbsD

clean:
	rm -f *.o *.Abs a.out *.mod *.AbsD
