CC = gcc
#CC = icc

#OPTIMIZE  = -O3
OPTIMIZE  = -Ofast -march=native
#OPTIMIZE  = -g -Wall
#  -g    adds debugging information to the executable file
#  -Wall turns on compiler warnings

### office version
#INCLUDES = -I/usr/include/lapacke
#LIBS = -llapacke -llapack -lblas -lgfortran

### home version
#INCLUDES = -I/usr/local/include/lapacke
#LIBS = -llapacke -lblas -lgfortran

CFLAGS = $(OPTIMIZE) -fopenmp
LIBS = -lm -fopenmp

.PHONY: all
all : libsky.so

### Choose a source for libsky.so
SOURCE = sky.c

OBJS = sky.o LLGsolver.o
INCL = LLGsolver.h Makefile
   
### -fPIC: position-independent code
LLGsolver.o : LLGsolver.c $(INCL)
	$(CC) $(CFLAGS) -c -fPIC LLGsolver.c -o LLGsolver.o

### -fPIC: position-independent code
sky.o : $(SOURCE) $(INCL)
	$(CC) $(CFLAGS) -c -fPIC $(SOURCE) -o sky.o

### -Wl,xx,yy: pass option as an options xx,yy to the linker
libsky.so : LLGsolver.o sky.o
	$(CC) $(CFLAGS) -shared -Wl,-soname,libsky.so -o libsky.so $(OBJS) $(LIBS)

### C interface for drift
cint-drift.o : cint-drift.c drift.h
	$(CC) $(CFLAGS) -c cint-drift.c -o cint-drift.o

drift.o : drift.c LLGsolver.h
	$(CC) $(CFLAGS) -c drift.c -o drift.o

cint-drift : cint-drift.o drift.o LLGsolver.o
	$(CC) $(CFLAGS) cint-drift.o drift.o LLGsolver.o -o cint-drift

### C interface for switch
cint-switch.o : cint-switch.c switch.h
	$(CC) $(CFLAGS) -c cint-switch.c -o cint-switch.o

switch.o : switch.c LLGsolver.h
	$(CC) $(CFLAGS) -c switch.c -o switch.o

cint-switch : cint-switch.o switch.o LLGsolver.o
	$(CC) $(CFLAGS) cint-switch.o switch.o LLGsolver.o -o cint-switch -lm

### C interface for switch range
cint-switch-range.o : cint-switch-range.c switch.h
	$(CC) $(CFLAGS) -c cint-switch-range.c -o cint-switch-range.o

cint-switch-range : cint-switch-range.o switch.o LLGsolver.o
	$(CC) $(CFLAGS) cint-switch-range.o switch.o LLGsolver.o -o cint-switch-range -lm

### C interface for switch range
cint-stripe-range.o : cint-stripe-range.c switch.h
	$(CC) $(CFLAGS) -c cint-stripe-range.c -o cint-stripe-range.o

cint-stripe-range : cint-stripe-range.o switch.o LLGsolver.o
	$(CC) $(CFLAGS) cint-stripe-range.o switch.o LLGsolver.o -o cint-stripe-range -lm
###

.PHONY : clean
clean :
	-rm -vf libsky.so sky.o LLGsolver.o cint-drift.o drift.o switch.o cint-switch.o cint-switch-range.o

