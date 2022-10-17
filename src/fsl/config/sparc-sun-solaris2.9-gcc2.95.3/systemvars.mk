# $Id: systemvars.mk,v 1.3 2004/09/09 09:54:42 steve Exp $

# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# Compiler dependent variables

CC = gcc
CXX = c++

ARCHFLAGS = -mv8 -ffast-math -fomit-frame-pointer
ARCHLDFLAGS = -static
DEPENDFLAGS = -MM

OPTFLAGS = -O6 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}

INSTALL = ginstall
RM = /bin/rm
CP = /bin/cp
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
RANLIB = echo
TCLSH = ${FSLDIR}/bin/tclsh

