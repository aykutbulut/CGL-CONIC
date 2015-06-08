# Copyright (C) 2015 Lehigh University.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Author:  Aykut Bulut, Lehigh University, 2015-06-04

#############################################################################
#                       Names and other basic things                        #
#############################################################################

AC_PREREQ(2.59)

AC_INIT([CglConic],[0.1],[https://github.com/aykutbulut/CGL-CONIC/issues/new])

AC_COPYRIGHT([
/*===========================================================================*/
/*                                                                           */
/* CglConic is developed by Aykut Bulut (aykutblt@gmail.com)                 */
/*                                                                           */
/* (c) Copyright 2015 Lehigh University. All Rights Reserved.                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/
])

# List one file in the package so that the configure script can test
# whether the package is actually there
AC_CONFIG_SRCDIR(src/CglConicCutGenerator.hpp)

# Where should everything be installed by default?  Here, we want it
# to be installed directly in 'bin', 'lib', 'include' subdirectories
# of the directory where configure is run.  The default would be
# /usr/local.
AC_PREFIX_DEFAULT([`pwd`])

#############################################################################
#                         Standard build tool stuff                         #
#############################################################################

# Get the system type
AC_CANONICAL_BUILD

# If this project depends on external projects, the Externals file in
# the source root directory contains definition of where to find those
# externals.  The following macro ensures that those externals are
# retrieved by svn if they are not there yet.
AC_COIN_PROJECTDIR_INIT(CglConic)

# Check if user wants to produce debugging code
AC_COIN_DEBUG_COMPILE(CglConic)

# Get the name of the C++ compiler and appropriate compiler options
AC_COIN_PROG_CXX

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

#############################################################################
#                            COIN-OR components                             #
#############################################################################
AC_COIN_CHECK_PACKAGE(CoinUtils, [coinutils], [CglConicLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils not available.])
fi
AC_COIN_CHECK_PACKAGE(Osi, [osi], [CglConicLib])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi not available.])
fi
AC_COIN_CHECK_PACKAGE(OsiClp,    [osi-clp],   [CglConicLib])
if test $coin_has_osiclp != yes ; then
  AC_MSG_ERROR([Required package OsiClp not available.])
fi
AC_COIN_CHECK_PACKAGE(OsiConic,  [osiconic],  [CglConicLib])
if test $coin_has_osiconic != yes ; then
  AC_MSG_ERROR([Required package OsiConic not available.])
fi
AC_COIN_CHECK_PACKAGE(Cgl,  [cgl],  [CglConicLib])
if test $coin_has_cgl != yes ; then
  AC_MSG_ERROR([Required package Cgl not available.])
fi
dnl AC_COIN_CHECK_PACKAGE(Sample, [coindatasample])

AC_COIN_CHECK_PACKAGE_BLAS(CglConicLib)
AC_COIN_CHECK_PACKAGE_LAPACK(CglConicLib)

#############################################################################
#                       Solvers for unittesting                             #
#############################################################################

# These are not really required for the compilation of Cgl, but we
# need to know which are available in order to compile the unitTest
# correctly.

AC_COIN_CHECK_PACKAGE(OsiCpx,  [osi-cplex])
AC_COIN_CHECK_PACKAGE(OsiGlpk, [osi-glpk])
#AC_COIN_CHECK_PACKAGE(OsiGrb,  [osi-gurobi])
AC_COIN_CHECK_PACKAGE(OsiMsk,  [osi-mosek])
AC_COIN_CHECK_PACKAGE(OsiXpr,  [osi-xpress])
AC_COIN_CHECK_PACKAGE(OsiVol,  [osi-vol])
AC_COIN_CHECK_PACKAGE(OsiDyLP, [osi-dylp])
#AC_COIN_CHECK_PACKAGE(OsiSpx,  [osi-soplex])

#############################################################################
#                  System header                                            #
#############################################################################

AC_COIN_CHECK_CXX_CHEADER(math)
AC_COIN_CHECK_CXX_CHEADER(float)
AC_COIN_CHECK_CXX_CHEADER(ieeefp)

#############################################################################
#                   Determine list of all CglConic subprojects              #
#############################################################################

# In order to avoid that one has to update this configure.ac and the
# Makefile.am whenever a new cut generator is added, we automatically determine
# the list of all subdirectories of src with a Makefile.in (CGL_SUBDIRS),
# and also a list of all sublibraries relative to src (CGL_SUBLIBS) that
# can be used in the src/Makefile.am.  (It is still necessary to manually
# add an entry to write the Makefile in the AC_CONFIG_FILES list below.)

AC_SUBST(CGLCONIC_SUBDIRS)
AC_SUBST(CGLCONIC_SUBLIBS)

curr_dir=`pwd`
cd $srcdir/src
for file in `ls */Makefile.in`; do
  newdir=`echo $file | sed -e s%/Makefile.in%% | sed -e s%^./%% `
  CGLCONIC_SUBDIRS="$CGLCONIC_SUBDIRS $newdir"
  CGLCONIC_SUBLIBS="$CGLCONIC_SUBLIBS $newdir/lib${newdir}.la"
done
cd $curr_dir

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi OsiConic Cgl)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Install files for unitTest
dnl AC_COIN_VPATH_LINK([test/CglTestData/capPlan1.mps test/CglTestData/l152lav.mps test/CglTestData/egout.mps test/CglTestData/lseu.mps])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/CglConicGD1/Makefile
                 src/CglConicGD2/Makefile
                 src/CglConicMIR/Makefile
                 cglconic.pc
                 cglconic-uninstalled.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cglconic.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
