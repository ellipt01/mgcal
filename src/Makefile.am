#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libmgcal.la

#######################################
# Build information for each library

# Sources for mgcal
libmgcal_la_SOURCES = \
	private/util.c\
	mgcal.c\
	vector3d.c\
	grid.c\
	scattered.c\
	data_array.c\
	source.c\
	calc.c\
	calc_bhattacharyya.c\
	kernel.c\
	mapped.c\
	io.c

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libmgcal_la_LIBADD = -L$(top_srcdir)/src
libmgcal_la_CPPFLAGS = -I$(top_srcdir)/include
libmgcal_la_LDFLAGS = -version-info $(MGCAL_LT_VERSION)
