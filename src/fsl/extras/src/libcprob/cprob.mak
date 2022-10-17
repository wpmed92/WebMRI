# Makefile for probability integrals.
# Be sure to set the type of computer and endianness in mconf.h.

CC = gcc
CFLAGS = -O2 -Wall
INCS = mconf.h

OBJS = incbetl.o incbil.o incbi.o qincb.o qincbi.o igami.o fdtrl.o
OBJS = bdtr.o btdtr.o chdtr.o drand.o fdtr.o gamma.o gdtr.o igam.o \
igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o

libprob.a: $(OBJS) $(INCS)
	ar rv libprob.a $(OBJS)
