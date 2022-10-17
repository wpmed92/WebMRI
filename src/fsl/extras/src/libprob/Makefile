include ${FSLCONFDIR}/default.mk

PROJNAME = libprob

OBJS = bdtr.o btdtr.o chdtr.o drand.o fdtr.o gamma.o gdtr.o igam.o \
igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o

XFILES = 
AFILES = libprob.a
HFILES = *.h

all: libprob.a

libprob.a: ${OBJS}
	${AR} -r libprob.a ${OBJS}

