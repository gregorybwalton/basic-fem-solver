CC=gcc
INCLUDE=-I/home/csunix/linux/apps/install/gsl/2.6/include
LDFLAGS=-L/home/csunix/linux/apps/install/gsl/2.6/lib
LDLIBS=-lm -ggdb -lgsl -lgslcblas
OBJS=fem.c \
../readmesh.c ../output.c ../calcerror.c ../addfunc.c ../listtool.c ../sparsemat.c \
gslsolve.c printspmatrices.c pc.c

obj = $(patsubst %.c,%.o,$(OBJS))

all: $(obj) fem

$(obj):%.o:%.c
	$(CC) -c $(INCLUDE) $(LDFLAGS) -o $@ $<
fem:
	$(CC) $(obj) $(INCLUDE) $(LDFLAGS) $(LDLIBS) -o fem 

clean:
	rm fem *.o ../*.o
