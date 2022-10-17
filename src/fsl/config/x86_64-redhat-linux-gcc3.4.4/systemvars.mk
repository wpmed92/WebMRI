# $Id: systemvars.mk,v 1.2 2006/01/13 10:48:46 flitney Exp $

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

ARCHFLAGS = -m64

DEPENDFLAGS = -MM

OPTFLAGS = -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -ansi -pedantic -Wno-deprecated -Wno-long-long
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

