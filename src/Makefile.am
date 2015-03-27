#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
lib_LTLIBRARIES = libmgcal.la

#######################################
# Build information for each library

# Sources for mgcal
libmgcal_la_SOURCES = \
	util.c\
	vector.c\
	grid.c\
	data_array.c\
	source.c\
	calc.c\
	kernel.c\
	io.c

# Linker options libTestProgram
libmgcal_la_LDFLAGS = 

libmgcal_la_LIBADD = -lm

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libmgcal_la_CPPFLAGS = -I$(top_srcdir)/include

