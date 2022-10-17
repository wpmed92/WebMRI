# The Makefile invokes emcc to generate the final wasm and JS wrappers for bet2
# FORCE_FILESYSTEM is needed to be able to use FS.mkdir() and other FS calls in the browsers
# ALLOW_MEMORY_GROWTH helps prevent OOM
include ${FSLCONFDIR}/default.mk

PROJNAME = bet2

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS=-lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

BETOBJS=bet2.o 

all: bet2.html

bet2.html: ${BETOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS}  -sFORCE_FILESYSTEM=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['FS']" -s ALLOW_MEMORY_GROWTH=1 -o $@ ${BETOBJS} ${LIBS}


