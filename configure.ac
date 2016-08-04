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
AC_LANG_PUSH(C++)

# Initialize automake and libtool
AC_COIN_INIT_AUTO_TOOLS

# Get the EXEXT variable for CYGWIN
AC_EXEEXT

#############################################################################
#                            COIN-OR components                             #
#############################################################################
AC_COIN_CHECK_PACKAGE(CoinUtils, [coinutils], [CglConicLib])
if test $coin_has_coinutils != yes ; then
  AC_MSG_ERROR([Required package CoinUtils is not available.])
fi

AC_COIN_CHECK_PACKAGE(Osi, [osi], [CglConicLib])
if test $coin_has_osi != yes ; then
  AC_MSG_ERROR([Required package Osi is not available.])
fi

AC_COIN_CHECK_PACKAGE(OsiClp, [osi-clp], [CglConicLib])
if test $coin_has_osiclp != yes ; then
  AC_MSG_ERROR([Required package OsiClp is not available.])
fi

AC_COIN_CHECK_PACKAGE(Cgl, [cgl], [CglConicLib])
if test $coin_has_cgl != yes ; then
  AC_MSG_ERROR([Required package Cgl is not available.])
fi

AC_COIN_CHECK_PACKAGE(OsiConic, [osiconic], [CglConicLib])
if test $coin_has_osiconic != yes ; then
  AC_MSG_ERROR([Required package OsiConic is not available.])
fi

AC_COIN_CHECK_PACKAGE(Cola, [cola], [CglConicLib])
if test $coin_has_cola != yes ; then
  AC_MSG_ERROR([Required package Cola is not available.])
fi

AC_COIN_CHECK_PACKAGE(OsiIpopt, [osiipopt], [CglConicLib])
if test $coin_has_osiipopt != yes ; then
  AC_MSG_ERROR([Required package OsiIpopt is not available.])
fi

#############################################################################
#                            Third party solvers                            #
#############################################################################

# Check which third party solvers are available
AC_COIN_CHECK_PACKAGE(Cplex, [osicplex], [CglConicLib])
AC_COIN_CHECK_PACKAGE(Mosek, [osimosek], [CglConicLib])

#############################################################################
#                      Check for IPM solver                                 #
#############################################################################
AC_MSG_CHECKING([for CglConic default IPM solver])
AC_ARG_WITH([ipm-solver],
[AC_HELP_STRING([--with-ipm-solver@<:@=ipmsolver@:>@],
                [specify the IPM solver in small
                letters(default ipmsolver=ipopt)])],
                [cglconic_ipm_solver=$withval],[cglconic_ipm_solver=ipopt])

found_package=true
case $cglconic_ipm_solver in
   ipopt)
     AC_MSG_RESULT(Ipopt)
     if [test $coin_has_osiipopt = unavailable || test $coin_has_osiipopt = skipping]; \
     then
         found_package=false
     fi
   ;;
   mosek)
     AC_MSG_RESULT(Mosek)
     if test $coin_has_mosek = false; \
     then
         found_package=false
     fi
   ;;
   cplex)
     AC_MSG_RESULT(Cplex)
     if test $coin_has_cplex = false; \
     then
         found_package=false
     fi
   ;;
   *)
   AC_MSG_ERROR(Unknown IPM solver!)
   ;;
esac

if test $found_package = false; then
  AC_MSG_ERROR(Couldn't find the $cglconic_ipm_solver package!)
fi

AM_CONDITIONAL(IPOPT_IPM_SOLVER,[test $cglconic_ipm_solver = ipopt &&
test $coin_has_coindepend = yes])
AM_CONDITIONAL(MOSEK_IPM_SOLVER,[test $cglconic_ipm_solver = mosek &&
test $coin_has_mosek = yes])
AM_CONDITIONAL(CPLEX_IPM_SOLVER,[test $cglconic_ipm_solver = cplex &&
test $coin_has_cplex = yes])

# Check thirdparty libraries
AC_COIN_CHECK_PACKAGE(Sample, [coindatasample])
AC_COIN_CHECK_PACKAGE_BLAS(CglConicLib)
AC_COIN_CHECK_PACKAGE_LAPACK(CglConicLib)

#############################################################################
#                       Solvers for unittesting                             #
#############################################################################

# These are not really required for the compilation of Cgl, but we
# need to know which are available in order to compile the unitTest
# correctly.

dnl AC_COIN_CHECK_PACKAGE(OsiCpx,  [osi-cplex])
dnl AC_COIN_CHECK_PACKAGE(OsiGlpk, [osi-glpk])
dnl #AC_COIN_CHECK_PACKAGE(OsiGrb,  [osi-gurobi])
dnl AC_COIN_CHECK_PACKAGE(OsiMsk,  [osi-mosek])
dnl AC_COIN_CHECK_PACKAGE(OsiXpr,  [osi-xpress])
dnl AC_COIN_CHECK_PACKAGE(OsiVol,  [osi-vol])
dnl AC_COIN_CHECK_PACKAGE(OsiDyLP, [osi-dylp])
dnl #AC_COIN_CHECK_PACKAGE(OsiSpx,  [osi-soplex])

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

##############################################################################
#                   VPATH links for example input files                      #
##############################################################################
# In case this is a VPATH configuration we need to make sure that the
# input files for the examples are available in the VPATH directory.
AC_COIN_VPATH_LINK(examples/ex1.mps)

#############################################################################
#                  Check for doxygen                                        #
#############################################################################

AC_COIN_DOXYGEN(CoinUtils Osi OsiConic Cgl OsiIpopt)

##############################################################################
#                   Finishing up by writing all the output                   #
##############################################################################

# Install files for unitTest
dnl AC_COIN_VPATH_LINK([test/CglTestData/capPlan1.mps test/CglTestData/l152lav.mps test/CglTestData/egout.mps test/CglTestData/lseu.mps])

# Here list all the files that configure should create (except for the
# configuration header file)
AC_CONFIG_FILES([Makefile
                 examples/Makefile
                 src/Makefile
                 src/CglConicGD1/Makefile
                 src/CglConicGD2/Makefile
                 src/CglConicMIR/Makefile
                 src/CglConicOA/Makefile
                 src/CglConicIPM/Makefile
                 src/CglConicIPMint/Makefile
                 cglconic.pc
                 cglconic-uninstalled.pc])
AC_CONFIG_FILES([doxydoc/doxygen.conf])

# Here put the location and name of the configuration header file
AC_CONFIG_HEADER([src/config.h src/config_cglconic.h])

# Finally, we let configure write all the output...
AC_COIN_FINALIZE
