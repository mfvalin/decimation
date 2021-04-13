
CC      = cc
CFLAGS  = -O3 -Wall
DEFINES =

all:	decimate.c.o decimate.Abs

decimate.c.o: decimate.c
	$(CC) -c $(CFLAGS) $< -o $@

decimate.Abs: decimate.c
	$(CC) -DSELF_TEST $(DEFINES) $(CFLAGS) $< -o $@

test:	decimate.Abs
	./$<

clean:
	rm -f *.o *.Abs
