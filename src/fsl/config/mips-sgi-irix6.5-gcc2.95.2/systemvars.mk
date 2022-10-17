# $Id: systemvars.mk,v 1.4 2003/06/04 08:00:28 mark Exp $

# System dependent paths

RM = /bin/rm
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
CP = /bin/cp
INSTALL = ginstall -p
TCLSH = ${FSLDIR}/bin/tclsh
RANLIB = echo

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++
CSTATICFLAGS =
CXXSTATICFLAGS =


ARCHFLAGS = 

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}
