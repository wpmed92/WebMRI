# $Id: Makefile,v 1.7 2005/11/03 15:17:45 duncan Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = meshclass

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB}

LIBS=-lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lz

MESHOBJS=point.o mpoint.o triangle.o mesh.o pt_special.o profile.o 

DBGFLAGS=-g

all: meshclass.a ${XFILES}

meshclass.a: ${MESHOBJS}
	${AR} -r libmeshclass.a ${MESHOBJS}



