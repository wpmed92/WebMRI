# $Id: vars.mk,v 1.8 2005/12/01 15:33:30 duncan Exp $

INCDIR = ${FSLDIR}/include
LIBDIR = ${FSLDIR}/lib

DEVINCDIR = ${FSLDEVDIR}/include
DEVLIBDIR = ${FSLDEVDIR}/lib

DESTDIR = ${FSLDEVDIR}

dest_INCDIR = ${DESTDIR}/include
dest_LIBDIR = ${DESTDIR}/lib
dest_BINDIR = ${DESTDIR}/bin
dest_MANDIR = ${DESTDIR}/man
dest_TCLDIR = ${DESTDIR}/tcl
dest_DOCDIR = ${DESTDIR}/doc
dest_REFDOCDIR = ${DESTDIR}/refdoc


PROJNAME =

USRLDFLAGS =
USRINCFLAGS = 
USRCFLAGS = 
USRCXXFLAGS =

LDFLAGS = ${ARCHLDFLAGS} ${USRLDFLAGS} -L. -L${DEVLIBDIR} -L${LIBDIR}

AccumulatedIncFlags = ${USRINCFLAGS} -I. -I${DEVINCDIR} -I${INCDIR}

CFLAGS = ${ANSI_FLAGS} ${DBGFLAGS} ${USEDCSTATICFLAGS} ${USRCFLAGS} ${ARCHFLAGS} ${OPTFLAGS}  \
	${AccumulatedIncFlags}

CXXFLAGS = ${ANSI_FLAGS} ${DBGFLAGS} ${USEDCXXSTATICFLAGS} ${USRCXXFLAGS} ${ARCHFLAGS} ${OPTFLAGS}  \
	${AccumulatedIncFlags}

HFILES = *.h
AFILES = *.a
XFILES = 
SCRIPTS =
TCLFILES = *.tcl
MANFILES = man/*

TESTXFILES =

DATATYPES =     if [ $$FDT = "8UI" ] ;  then STR="unsigned char" ; fi ; \
                if [ $$FDT = "8SI" ] ;  then STR="signed char" ; fi ; \
                if [ $$FDT = "16UI" ] ; then STR="unsigned short" ; fi ; \
                if [ $$FDT = "16SI" ] ; then STR="signed short" ; fi ; \
                if [ $$FDT = "32UI" ] ; then STR="unsigned int" ; fi ; \
                if [ $$FDT = "32SI" ] ; then STR="signed int" ; fi ; \
                if [ $$FDT = "32R" ] ;  then STR="float" ; fi

