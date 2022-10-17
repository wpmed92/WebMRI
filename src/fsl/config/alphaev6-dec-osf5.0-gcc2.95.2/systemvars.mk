 


# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh


# Compiler dependent variables

CC = gcc
CXX = c++
CSTATICFLAGS = -static
CXXSTATICFLAGS = -static

ARCHFLAGS = -mieee -mfp-trap-mode=sui

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}
 
 
# Variables determined with AUTOCONFIG: 
 
INSTALL = ginstall -p
RM = /bin/rm
CP = /bin/cp
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
RANLIB = ranlib
TCLSH = ${FSLDIR}/bin/tclsh
