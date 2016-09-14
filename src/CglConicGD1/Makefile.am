# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University, 2015-06-04

AUTOMAKE_OPTIONS = foreign

########################################################################
#                          libCglConicGD1                              #
########################################################################

# Name of the library compiled in this directory.  We don't want it to be
# installed since it will be collected into the libCgl library
noinst_LTLIBRARIES = libCglConicGD1.la

# List all source files for this library, including headers
libCglConicGD1_la_SOURCES = CglConicGD1.cpp\
		CglConicGD1.hpp\
		CglConicGD1Help.hpp\
		CglConicGD1Help.cpp\
		CglConicGD1Cut.cpp\
		CglConicGD1Cut.hpp\
		CglConicGD1Param.cpp\
		CglConicGD1Param.hpp

# This is for libtool (on Windows)
libConicGD1_la_LDFLAGS = $(LT_LDFLAGS)

# Here list all include flags, relative to this "srcdir" directory.  This
# "cygpath" stuff is necessary to compile with native compilers on Windows.
AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` $(COINUTILS_CFLAGS)\
	$(OSI_CFLAGS)\
	$(OSICONIC_CFLAGS)\
	$(CGL_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I..

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = CglConicGD1.hpp\
		CglConicGD1Help.hpp\
		CglConicGD1Cut.hpp\
		CglConicGD1Param.hpp
