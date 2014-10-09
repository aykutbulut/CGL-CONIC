# We need to recurse into all Cgl subproject subdirectories
SUBDIRS=$(CGLCONIC_SUBDIRS)

# compile all usual Cgl cuts and new conic cuts to CglConic
lib_LTLIBRARIES=libCglConic.la

# empty resources, we do not have files common to all CglConic
libCglConic_la_SOURCES=

libCglConic_la_LIBADD=\
	$(CGLCONIC_SUBLIBS)

# no cpp flags
AM_CPPFLAGS=
