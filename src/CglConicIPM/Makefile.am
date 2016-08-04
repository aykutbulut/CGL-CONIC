# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University, 2015-06-04

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libCglConicIPM                              #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglConicIPM.la

# List all source files for this library, including headers
libCglConicIPM_la_SOURCES = CglConicIPM.cpp\
		CglConicIPM.hpp\
		CglConicIPMParam.cpp\
		CglConicIPMParam.hpp

# This is for libtool (on Windows)
libConicIPM_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(COINUTILS_CFLAGS)\
	$(OSI_CFLAGS)\
	$(OSICONIC_CFLAGS)\
	$(CGL_CFLAGS)\
	$(OSIMOSEK_CFLAGS)

if IPOPT_IPM_SOLVER
  AM_CPPFLAGS += -D__IPOPT__
endif

if MOSEK_IPM_SOLVER
  AM_CPPFLAGS += -D__OSI_MOSEK__
endif

if CPLEX_IPM_SOLVER
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = CglConicIPM.hpp\
		CglConicIPMParam.hpp
