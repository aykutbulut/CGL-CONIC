noinst_LTLIBRARIES = libCglConicGD1.la

libCglConicGD1_la_SOURCES=CglConicGD1.hpp\
	CglConicGD1.cpp\
	CglConicGD1Param.hpp\
	CglConicGD1Param.cpp

libCglConicGD1_la_LIBADD=\
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
include_HEADERS=CglConicGD1.hpp\
	CglConicGD1Param.hpp

# this are in the disctribution but not in include
#noinst_HEADERS=
