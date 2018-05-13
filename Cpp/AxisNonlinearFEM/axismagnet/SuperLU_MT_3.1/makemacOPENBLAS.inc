############################################################################
#
#  Program:         SuperLU_MT
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   August 15, 1997
#
#  Modified:	    September 1, 1999 version 1.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _PTHREAD

#
#  The name of the libraries to be created/linked to
#
TMGLIB       	= libtmglib$(PLAT).a
SUPERLULIB   	= libsuperlu_mt$(PLAT).a
BLASDEF	        = -DUSE_VENDOR_BLAS
BLASLIB		= -L/usr/lib -lf77blas -latlas
BLASLIB =  -I ../../OpenBLAS/lib/include ../../OpenBLAS/lib/lib/libopenblas.a
MATHLIB		= -lm
MPLIB		= -lpthread

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH          	= /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ar
ARCHFLAGS     	= cr
RANLIB       	= /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib

############################################################################
#
# Definitions used by CPP:
PREDEFS      	= -D__PTHREAD  ## -DDEBUGlevel=0 -DPRNTlevel=0
#
CC           	= /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc
CFLAGS       	= -O3 -DNDEBUG -DPRNTlevel=0 -DDEBUGlevel=0 -mmacosx-version-min=10.7
# Uncomment the following to use 64-bit integer
#It seems that mac doesn't support 64-bit integer, so comment it.
#CFLAGS          += -D_LONGINT

#
# There are a few files that should NOT be optimized, which appear in the
# end of SRC/Makefile, taking the flag $(NOOPTS)
# NOOPTS should be defined to explicitly turn off any optimization.
NOOPTS	     	= -O0

BLASOPTS	=

############################################################################

FORTRAN		= /usr/local/bin/gfortran
FFLAGS		= -O -fdefault-integer-8

############################################################################

LOADER       	= $(CC)
LOADOPTS     	= -O #-pthread

############################################################################
#
#  C preprocessor defs for compilation (-DNoChange, -DAdd_, or -DUpCase)
#
CDEFS        	= -DAdd_
