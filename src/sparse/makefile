CC=gcc
LIBS= -L/home/csunix/linux/apps/install/gsl/2.6/lib -I/home/csunix/linux/apps/install/gsl/2.6/include
CFLAGS=-lm -ggdb -lgsl -lgslcblas #-pg
OBJS=fem.c \
../readmesh.c ../output.c ../calcerror.c ../addfuncs.c ../listtool.c\
spmatrices.c spsolve.c printspmatrices.c pc.c

obj = $(patsubst %.c,%.o,$(OBJS))

fem: $(obj)
	$(CC) $(obj) $(LIBS) -o fem $(CFLAGS)
$(obj):%.o:%.c
	$(CC) -c $(LIBS) -o $@ $< $(CFLAGS)

clean:
	rm fem *.o ../*.o
