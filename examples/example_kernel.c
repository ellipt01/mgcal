/*
 * example_kernel.c
 *
 *  Created on: 2015/03/26
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <omp.h>

#include "mgcal.h"

data_array *
create_data_array (const grid *g)
{
	int			i, j, k;
	data_array	*array = data_array_new (g->nh);

	k = 0;
	for (j = 0; j < g->ny; j++) {
		for (i = 0; i < g->nx; i++) {
			array->x[k] = g->x[i];
			array->y[k] = g->y[j];
			array->z[k] = g->z[0];
			array->data[k] = 0.;
			k++;
		}
	}
	return array;
}

double *
example_kernel (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs)
{
	data_array	*array;
	grid		*gobs;
	grid		*gsrc;
	mgcal_func	*f;
	cvector		*exf;
	cvector		*mgz;
	double		*a;

	double		t;

	gobs = grid_new (nx, ny, 1, x, y, &zobs);
	array = create_data_array (gobs);

	gsrc = grid_new (nx, ny, nz, x, y, z);
	exf = cvector_new_with_geodesic_poler (1., 45., 0.);
	mgz = cvector_new_with_geodesic_poler (5., 45., 0.);
	f = mgcal_func_new (total_force_dipole, NULL);

	t = omp_get_wtime ();
	a = kernel_matrix (array, gsrc, mgz, exf, f);
	fprintf (stderr, "time = %.4e\n", omp_get_wtime () - t);

	data_array_free (array);
	grid_free (gsrc);
	grid_free (gobs);
	mgcal_func_free (f);
	cvector_free (exf);
	cvector_free (mgz);

	return a;
}
