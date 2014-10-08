noinst_LTLIBRARIES = libCglConicMIR.la

libCglConicMIR_la_SOURCES=CglConicMIR.hpp\
	CglConicMIR.cpp\
	CglConicMIRParam.hpp\
	CglConicMIRParam.cpp

libCglConicMIR_la_LIBADD=\
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
include_HEADERS=CglConicMIR.hpp\
	CglConicMIRParam.hpp

# this are in the disctribution but not in include
#noinst_HEADERS=
