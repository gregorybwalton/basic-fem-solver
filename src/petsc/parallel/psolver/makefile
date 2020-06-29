PETSC_DIR=/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/petsc/3.012.4
PETSC_ARCH=arch-linux2-c-debug

CC=mpicc
#CC=gcc
INCLUDE=-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include
LDFLAGS=-L${PETSC_DIR}/${PETSC_ARCH}/lib
LDLIBS=-lm -ggdb\
	-lpetsc
OBJS=mpisolve.c

all: main

#include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

obj = $(patsubst %.c,%.o,$(OBJS))

main: $(obj)
	$(CC) $(INCLUDE) $(LDFLAGS) $(LDLIBS) -o solve $(obj)
$(obj): %.o:%.c
	$(CC) -c $(INCLUDE) $(LDFLAGS) $(LDLIBS) $< -o $@

clean::
	rm solve *.o ../*.o
