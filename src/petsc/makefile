#PETSC_DIR=/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/petsc/3.012.4
#PETSC_ARCH=arch-linux2-c-debug
PETSC_DIR=/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/petsc
PETSC_ARCH=petsc-3.012.4-openmpi

CC=mpicc
#CC=gcc
INCLUDE=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LDFLAGS=-L${PETSC_DIR}/${PETSC_ARCH}/lib
#LDLIBS=-lm -ggdb
LDLIBS=-lm -w\
	-lpetsc
OBJS=fem.c\
	../readtetgen.c ../output.c ../calcerror.c ../addfunc.c ../listtool.c ../d3n4.c \
	psolve.c \
	./proto/meshvolume.c
	#./proto/meshlengthscale.c 

all: main

#include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/rules

obj = $(patsubst %.c,%.o,$(OBJS))

main: $(obj)
	$(CC) $(INCLUDE) $(LDFLAGS) $(LDLIBS) -o fem $(obj)
$(obj): %.o:%.c
	$(CC) -c $(INCLUDE) $(LDFLAGS) $(LDLIBS) $< -o $@

clean::
	rm fem *.o ../*.o ./proto/*.o
