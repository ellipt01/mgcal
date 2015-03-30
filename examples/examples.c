/*
 ============================================================================
 Name        : exampleProgram.c
 Author      : Mitsuru Utsugi
 Version     :
 Copyright   : none
 Description : Uses shared library to print greeting
               To run the resulting executable the LD_LIBRARY_PATH must be
               set to ${project_loc}/grid/.libs
               Alternatively, libtool creates a wrapper shell script in the
               build directory of this program which can be used to run it.
               Here the script will be called exampleProgram.
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "mgcal.h"
#include "examples.h"

const int		nx = 50;
const int		ny = 50;
const int		nz = 25;
const double	x[] = {-10., 10.};
const double	y[] = {-10., 10.};
const double	z[] = {0., -10.};
const double	zobs = 0.1;

int
main (void)
{
	FILE	*fp;

	if ((fp = fopen ("dipole.data", "w"))) {
		example_dipole (fp, nx, ny, 1, x, y, &zobs);
		fclose (fp);
	}

	if ((fp = fopen ("prism.data", "w"))) {
		example_prism (fp, nx, ny, 1, x, y, &zobs);
		fclose (fp);
	}

	if ((fp = fopen ("dipole_irregular.data", "w"))) {
		example_dipole_irregular_surface (fp, nx, ny, 1, x, y, &zobs);
		fclose (fp);
	}

	if ((fp = fopen ("dipole_multi_sources.data", "w"))) {
		example_dipole_multi_sources (fp, nx, ny, 1, x, y, &zobs);
		fclose (fp);
	}

	if ((fp = fopen ("kernel.data", "w"))) {
		double	*a = example_kernel (nx, ny, nz, x, y, z, zobs);
		int		i;
		for (i = 0; i < nx * ny * nz; i++) fprintf (fp, "%f\n", a[i]);
		fclose (fp);
	}

	return EXIT_SUCCESS;
}
