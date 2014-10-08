noinst_LTLIBRARIES = libCglConicGD2.la

libCglConicGD2_la_SOURCES=CglConicGD2.hpp\
	CglConicGD2.cpp\
	CglConicGD2Param.hpp\
	CglConicGD2Param.cpp

libCglConicGD2_la_LIBADD=\
	$(CGL_LIBS)\
	$(OSICONIC_LIBS)\
	$(OSICLP_LIBS)\
	$(OSI_LIBS)\
	$(COINUTILS_LIBS)\
	$(LIBS)

AM_CPPFLAGS=\
	-I../../include\
	$(CGL_CFLAGS)\
	$(OSICONIC_CFLAGS)\
	$(OSICLP_CFLAGS)\
	$(OSI_CFLAGS)\
	$(CLP_CFLAGS)\
	$(COINUTILS_CFLAGS)


# header files to be installed to $(includedir)
include_HEADERS=CglConicGD2.hpp\
	CglConicGD2Param.hpp

# this are in the disctribution but not in include
#noinst_HEADERS=
