# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University, 2015-06-04

AUTOMAKE_OPTIONS = foreign

# We need to recurse into all Cgl subproject subdirectories
SUBDIRS = $(CGLCONIC_SUBDIRS)

########################################################################
#                           Common Flags                               #
########################################################################

LDADD = $(CGLCONICLIB_LIBS)
DEPENDENCIES = $(CGLCONICLIB_DEPENDENCIES)
AM_CPPFLAGS = $(CGLCONICLIB_CFLAGS)

########################################################################
#                    Define Macros for IPM Solver                      #
########################################################################

if IPOPT_IPM_SOLVER
  AM_CPPFLAGS += -D__IPOPT__
endif

if MOSEK_IPM_SOLVER
  AM_CPPFLAGS += -D__OSI_MOSEK__
endif

if CPLEX_IPM_SOLVER
  AM_CPPFLAGS += -D__OSI_CPLEX__
endif

########################################################################
#                               libCglConic                            #
########################################################################

# Name of the library compiled in this directory.  We want it to be installed
# in the $libdir directory
lib_LTLIBRARIES = libCglConic.la

# List all source files for this library, including headers
libCglConic_la_SOURCES = \
	CglConicConfig.h \
	CglConicCutGenerator.cpp CglConicCutGenerator.hpp\
	CglConicParam.cpp CglConicParam.hpp

# We want to have all the sublibraries from the CglConic subprojects
# collected into this library
libCglConic_la_LIBADD = $(CGLCONIC_SUBLIBS)
if DEPENDENCY_LINKING
libCglConic_la_LIBADD += $(CGLCONICLIB_LIBS)
endif

# This is for libtool (on Windows)
libCglConic_la_LDFLAGS = $(LT_LDFLAGS)

# Since automake is not doing this on its own, we need to declare the
# dependencies to the subdirectory libraries here
libCglConic_la_DEPENDENCIES = $(CGLCONIC_SUBLIBS)


# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)`

########################################################################
#                Headers that need to be installed                     #
########################################################################

# Here list all the header files that are required by a user of the library,
# and that therefore should be installed in 'include/coin'
includecoindir = $(includedir)/coin
includecoin_HEADERS = \
	CglConicCutGenerator.hpp \
	CglConicParam.hpp

install-exec-local:
	$(install_sh_DATA) config_cglconic.h $(DESTDIR)$(includecoindir)/CglConicConfig.h

uninstall-local:
	rm -f $(DESTDIR)$(includecoindir)/CglConicConfig.h
