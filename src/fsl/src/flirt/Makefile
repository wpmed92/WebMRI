# The Makefile invokes emcc to generate the final wasm and JS wrappers for flirt
# FORCE_FILESYSTEM is needed to be able to use FS.mkdir() and other FS calls in the browsers
# ALLOW_MEMORY_GROWTH helps prevent OOM

include ${FSLCONFDIR}/default.mk

PROJNAME = flirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lm -lz

FL_OBJS = globaloptions.o flirt.o 

all: flirt.html

flirt.html: ${FL_OBJS}
	        $(CXX)  ${CXXFLAGS} ${LDFLAGS} -sFORCE_FILESYSTEM=1 -s "EXTRA_EXPORTED_RUNTIME_METHODS=['FS']" -s ALLOW_MEMORY_GROWTH=1  -o $@ ${FL_OBJS} ${LIBS}
