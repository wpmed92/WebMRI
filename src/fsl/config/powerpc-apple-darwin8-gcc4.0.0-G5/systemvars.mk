 


# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh


# Compiler dependent variables

CC = cc
CXX = c++
CSTATICFLAGS = 
CXXSTATICFLAGS = 


ARCHFLAGS = -arch ppc64
ARCHLDFLAGS = -Wl,-search_paths_first -arch ppc64

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations -mcpu=G5 -mtune=G5 -mpowerpc64 -mpowerpc-gpopt
MACHDBGFLAGS = -g
GNU_ANSI_FLAGS = -Wall -Wno-long-long -Wno-long-double -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS}
 
 
# Variables determined with AUTOCONFIG: 
 
INSTALL = install -p -c
RM = /bin/rm
CP = /bin/cp
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
RANLIB = ranlib
TCLSH = ${FSLDIR}/bin/tclsh
