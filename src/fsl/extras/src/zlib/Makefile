include ${FSLCONFDIR}/default.mk

PROJNAME = zlib

OBJS = adler32.o compress.o crc32.o gzio.o uncompr.o deflate.o trees.o \
       zutil.o inflate.o infblock.o inftrees.o infcodes.o infutil.o inffast.o

DBGFLAGS = 

XFILES = 
AFILES = libz.a

all: zlib.a

zlib.a: ${OBJS}
	${AR} -r libz.a ${OBJS}

