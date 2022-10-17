# $Id: systemvars.mk,v 1.1 2004/08/28 12:47:25 steve Exp $

# System dependent paths

RM = /bin/rm
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
CP = /bin/cp
INSTALL = install -p
TCLSH = ${FSLDIR}/bin/tclsh
RANLIB = echo

FSLML = ${FSLDIR}/bin/fslml

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = -march=pentium -mcpu=pentium

DEPENDFLAGS = -MM

OPTFLAGS = -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -ansi -pedantic -Wno-deprecated
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

