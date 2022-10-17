include ${FSLCONFDIR}/default.mk

PROJNAME = libpng

USRINCFLAGS = -I${INC_ZLIB}

DBGFLAGS = 


OBJS = png.o pngset.o pngget.o pngrutil.o pngtrans.o pngwutil.o \
	pngread.o pngrio.o pngwio.o pngwrite.o pngrtran.o \
	pngwtran.o pngmem.o pngerror.o pngpread.o

LIBS = -lpng -lz -lm

XFILES = 

all: libpng.a 

libpng.a: ${OBJS}
	${AR} -r libpng.a ${OBJS}

pngtest: pngtest.o libpng.a
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ pngtest.o ${LIBS}

