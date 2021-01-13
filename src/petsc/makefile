PETSC_DIR=/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/petsc
PETSC_ARCH=petsc-3.012.4-openmpi
TETGEN_DIR=/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/fem-codes/fem-working/ext/tetgen1.6.0

#CC=mpicc
CC=mpic++
#CC=g++
PI=-I${PETSC_DIR}/${PETSC_ARCH}/include
PLD=-L${PETSC_DIR}/${PETSC_ARCH}/lib
TI=-I${TETGEN_DIR}
TLD=-L${TETGEN_DIR}
INCLUDE=$(PI) $(TI)
LDFLAGS=$(PLD) $(TLD)
#LDLIBS=-lm -w -ggdb
LDLIBS=-lm -w -ggdb\
	-lpetsc -ltet
OBJC=fem.c\
	../readtetgen.c ../output.c ../calcerror.c ../gaussfunc.c ../listtool.c ../d3n4.c \
	psolve.c meshvolume.c ../readgrummp.c \
	./proto/calcmass.c ./proto/movemesh.c
OBJCX= ./runtetgen.cxx

all: main

#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/rules

objc = $(patsubst %.c,%.o,$(OBJC))
objcx = $(patsubst %.cxx,%.o,$(OBJCX))

main:	$(objc) $(objcx)
#main:	$(objc)
	$(CC) -o fem $(objc) $(objcx) $(INCLUDE) $(LDFLAGS) $(LDLIBS)
#	$(CC) -o fem $(objc) $(INCLUDE) $(LDFLAGS) $(LDLIBS)
$(objc):	%.o:%.c
		$(CC) -c $(INCLUDE) $(LDFLAGS) $(LDLIBS) $< -o $@
$(objcx):	%.o:%.cxx
		$(CC) -c $(INCLUDE) $(LDFLAGS) $(LDLIBS) $< -o $@

clean::
	rm fem *.o ../*.o ./proto/*.o
