AC_PREREQ([2.69])
AC_INIT([cola], [0.1], [aykut@lehigh.edu])
AC_CONFIG_SRCDIR([config.h.in])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_CONFIG_HEADERS([config.h])

# for creating library
AC_PROG_LIBTOOL

# Checks for programs.
AC_PROG_CXX

# Installs to ./bin, changes default behavior
AC_PREFIX_DEFAULT([`pwd`])

# Checks for header files.
AC_CHECK_HEADERS([math.h])

# Checks for typedefs, str
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Check dependent libraries
AC_CHECK_LIB([m], [sqrt])
# Check dependent libraries through pkg-config
# cgl osiconic osi-clp osi coinutils
PKG_CHECK_MODULES([CGL], [cgl >= 0.0])
PKG_CHECK_MODULES([OSICONIC], [osiconic >= 0.0])
PKG_CHECK_MODULES([OSICLP], [osi-clp >= 0.105])
PKG_CHECK_MODULES([OSI], [osi >= 0.105])
PKG_CHECK_MODULES([COINUTILS], [coinutils >= 2.8])




#############################################################################
#                   Determine list of all Cgl subprojects                   #
#############################################################################

# In order to avoid that one has to update this configure.ac and the
# Makefile.am whenever a new cut generator is added, we automatically determine
# the list of all subdirectories of src with a Makefile.in (CGL_SUBDIRS),
# and also a list of all sublibraries relative to src (CGL_SUBLIBS) that
# can be used in the src/Makefile.am.  (It is still necessary to manually
# add an entry to write the Makefile in the AC_CONFIG_FILES list below.)
#
# modified from COIN-OR CGL's configure.ac

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


AC_CONFIG_FILES([Makefile
	src/Makefile
	src/CglConicMIR/Makefile
	src/CglConicGD1/Makefile
	src/CglConicGD2/Makefile
	cglconic.pc
	cglconic-uninstalled.pc])

AC_OUTPUT
