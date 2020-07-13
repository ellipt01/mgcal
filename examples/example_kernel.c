/*
 * example_kernel.c
 *
 *  Created on: 2015/03/26
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include "mgcal.h"

extern const double	dec;
extern const double	inc;

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
	vector3d	*exf;
	vector3d	*mgz;
	double		*a;

	double		t;

	gobs = grid_new (nx, ny, 1, x, y, &zobs);
	array = create_data_array (gobs);

	gsrc = grid_new (nx, ny, nz, x, y, z);
	exf = vector3d_new_with_geodesic_poler (1., inc, dec);
	mgz = vector3d_new_with_geodesic_poler (1., inc, dec);
	f = mgcal_func_new (total_force_dipole, NULL);
//	f = mgcal_func_new (total_force_prism, NULL);

	t = omp_get_wtime ();
	a = kernel_matrix (array, gsrc, mgz, exf, f);
	fprintf (stderr, "time1 = %.4e\n", omp_get_wtime () - t);

	data_array_free (array);
	grid_free (gsrc);
	grid_free (gobs);
	mgcal_func_free (f);
	vector3d_free (exf);
	vector3d_free (mgz);

	return a;
}

double *
example_kernel_nth_grid (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs)
{
	data_array	*array;
	grid		*gobs;
	grid		*gsrc;
	mgcal_func	*f;
	vector3d	*exf;
	vector3d	*mgz;
	double		*a;

	double		t;

	gobs = grid_new (nx, ny, 1, x, y, &zobs);
	array = create_data_array (gobs);

	gsrc = grid_new (nx, ny, nz, x, y, z);
	exf = vector3d_new_with_geodesic_poler (1., inc, dec);
	mgz = vector3d_new_with_geodesic_poler (1., inc, dec);
	f = mgcal_func_new (total_force_dipole, NULL);

	a = (double *) malloc (array->n * gsrc->n * sizeof (double));
	t = omp_get_wtime ();
	{
		int		i, j;
#pragma omp parallel for
		for (j = 0; j < gsrc->n; j++) {
			double	*aj = kernel_matrix_nth_grid (j, array, gsrc, mgz, exf, f);
			for (i = 0; i < array->n; i++) a[i + j * array->n] = aj[i];
			free (aj);
		}
	}
	fprintf (stderr, "time2 = %.4e\n", omp_get_wtime () - t);

	data_array_free (array);
	grid_free (gsrc);
	grid_free (gobs);
	mgcal_func_free (f);
	vector3d_free (exf);
	vector3d_free (mgz);

	return a;
}

double *
example_kernel_mth_site (const int nx, const int ny, const int nz, const double x[], const double y[], const double z[], const double zobs)
{
	data_array	*array;
	grid		*gobs;
	grid		*gsrc;
	mgcal_func	*f;
	vector3d	*exf;
	vector3d	*mgz;
	double		*a;

	double		t;

	gobs = grid_new (nx, ny, 1, x, y, &zobs);
	array = create_data_array (gobs);

	gsrc = grid_new (nx, ny, nz, x, y, z);
	exf = vector3d_new_with_geodesic_poler (1., inc, dec);
	mgz = vector3d_new_with_geodesic_poler (1., inc, dec);
	f = mgcal_func_new (total_force_dipole, NULL);

	a = (double *) malloc (array->n * gsrc->n * sizeof (double));
	t = omp_get_wtime ();
	{
		int		i, j;
#pragma omp parallel for
		for (i = 0; i < array->n; i++) {
			double	*ai = kernel_matrix_mth_site (i, array, gsrc, mgz, exf, f);
			for (j = 0; j < gsrc->n; j++) a[i + j * array->n] = ai[j];
			free (ai);
		}
	}
	fprintf (stderr, "time3 = %.4e\n", omp_get_wtime () - t);

	data_array_free (array);
	grid_free (gsrc);
	grid_free (gobs);
	mgcal_func_free (f);
	vector3d_free (exf);
	vector3d_free (mgz);

	return a;
}
