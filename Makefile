# Names of source files included:
include Makefile.srcs

LIBS += -L/usr/local/lib
CPPFLAGS = -I. -I./test

CC=gcc -g -Wall

CCFLAGS = -DTEST
CCFLAGS += -I/usr/include/cfitsio
LDFLAGS=-lm -lgsl -lgslcblas -lcfitsio -lcpgplot

SRCS = ${PROJ_SOURCES} ${TEST_SOURCES}
INCS = ${PROJ_INCLUDES} ${TEST_INCLUDES}
OBJS = $(SRCS:.c=.o)

%.o: %.c
	${CC} ${CCFLAGS} ${CPPFLAGS} -c -o $@ $<

all: mytest

mytest: ${OBJS}
	${CC} -o $@ $^ ${LIBS} ${LDFLAGS}  

.PHONY: clean

clean:
	rm -f src/*.o test/*.o *~ src/*~ test/*~ core mytest ./junk*

